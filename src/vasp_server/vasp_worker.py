import os
import shutil
import asyncio
import traceback
from typing import Dict, Any, Optional
import requests
from pathlib import Path
import subprocess
import time

from .mp import download_with_criteria
from .base import cif_to_poscar
from .Config import get_path_config, get_kpoints_config

class VaspWorker:
    """VASP计算工作器"""
    
    def __init__(self, user_id: str, base_work_dir: str = "/data/home/ysl9527/vasp_calculations"):
        self.user_id = user_id
        self.base_work_dir = Path(base_work_dir) / user_id  # 为每个用户创建独立目录
        self.base_work_dir.mkdir(parents=True, exist_ok=True)
    
    async def run_structure_optimization(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """
        运行结构优化计算
        
        Args:
            task_id: 任务ID
            params: 任务参数
            progress_callback: 进度回调函数
            
        Returns:
            Dict: 计算结果
        """
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # 更新进度: 开始处理
            if progress_callback:
                await progress_callback(5, "开始处理输入参数...")
            
            # 1. 获取CIF文件
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            
            # 2. 转换为POSCAR
            if progress_callback:
                await progress_callback(10, "转换CIF为POSCAR...")
            poscar_path = await self._convert_cif_to_poscar(cif_path, work_dir, params)
            
            # 3. 生成VASP输入文件
            if progress_callback:
                await progress_callback(20, "生成VASP输入文件...")
            await self._generate_vasp_inputs(work_dir, params)
            
            # 4. 运行VASP计算
            if progress_callback:
                await progress_callback(30, "开始VASP计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)
            
            # 5. 分析结果
            if progress_callback:
                await progress_callback(90, "分析计算结果...")
            final_result = await self._analyze_results(work_dir, result)
            
            if progress_callback:
                await progress_callback(100, "计算完成！")
                
            return final_result
            
        except Exception as e:
            error_msg = f"结构优化计算失败: {str(e)}"
            print(f"[ERROR] {error_msg}")
            print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
            raise Exception(error_msg)
    
    async def run_scf_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """
        运行自洽场计算
        
        Args:
            task_id: 任务ID
            params: 任务参数
            progress_callback: 进度回调函数
            
        Returns:
            Dict: 计算结果
        """
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # 更新进度: 开始处理
            if progress_callback:
                await progress_callback(5, "开始自洽场计算...")
            
            # 1. 获取结构文件
            poscar_path = await self._get_structure_for_scf(work_dir, params, progress_callback)
            if not poscar_path:
                raise Exception("无法获取结构文件")
            
            # 2. 生成自洽场VASP输入文件
            if progress_callback:
                await progress_callback(30, "生成自洽场VASP输入文件...")
            await self._generate_scf_inputs(work_dir, params)
            
            # 3. 运行VASP自洽场计算
            if progress_callback:
                await progress_callback(40, "开始VASP自洽场计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)
            
            # 4. 分析自洽场结果
            if progress_callback:
                await progress_callback(90, "分析自洽场计算结果...")
            final_result = await self._analyze_scf_results(work_dir, result)
            
            if progress_callback:
                await progress_callback(100, "自洽场计算完成！")
                
            return final_result
            
        except Exception as e:
            error_msg = f"自洽场计算失败: {str(e)}"
            print(f"[ERROR] {error_msg}")
            print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
            raise Exception(error_msg)
    
    async def run_dos_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """
        运行态密度计算
        
        Args:
            task_id: 任务ID
            params: 任务参数
            progress_callback: 进度回调函数
            
        Returns:
            Dict: 计算结果
        """
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # 更新进度: 开始处理
            if progress_callback:
                await progress_callback(5, "开始态密度计算...")
            
            # 1. 获取结构文件和准备文件
            dos_files = await self._prepare_dos_files(work_dir, params, progress_callback)
            if not dos_files:
                raise Exception("无法准备DOS计算文件")
            
            # 2. 生成态密度VASP输入文件
            if progress_callback:
                await progress_callback(30, "生成态密度VASP输入文件...")
            await self._generate_dos_inputs(work_dir, params, dos_files)
            
            # 3. 运行VASP态密度计算
            if progress_callback:
                await progress_callback(40, "开始VASP态密度计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)
            
            # 4. 分析态密度结果
            if progress_callback:
                await progress_callback(90, "分析态密度计算结果...")
            final_result = await self._analyze_dos_results(work_dir, result)
            
            if progress_callback:
                await progress_callback(100, "态密度计算完成！")
                
            return final_result
            
        except Exception as e:
            error_msg = f"态密度计算失败: {str(e)}"
            print(f"[ERROR] {error_msg}")
            print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
            raise Exception(error_msg)
    
    async def _get_cif_file(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[str]:
        """获取CIF文件"""
        if params.get('formula'):
            # 从Materials Project下载
            formula = params['formula']
            if progress_callback:
                await progress_callback(10, f"从Materials Project下载 {formula}...")
            
            # 构建搜索条件
            search_kwargs = {}
            for key in ['spacegroup', 'max_energy_above_hull', 'min_band_gap', 
                       'max_band_gap', 'max_nsites', 'min_nsites', 'stable_only', 'selection_mode']:
                if key in params and params[key] is not None:
                    search_kwargs[key] = params[key]
            
            cif_path = await download_with_criteria(
                formula=formula,
                save_path=str(work_dir),
                task_id=params.get('task_id', 'unknown'),
                **search_kwargs
            )
            return str(cif_path) if cif_path else None
            
        elif params.get('cif_url'):
            # 从URL下载
            cif_url = params['cif_url']
            if progress_callback:
                await progress_callback(15, f"从URL下载CIF: {cif_url}")
            
            cif_path = work_dir / "structure.cif"
            response = requests.get(str(cif_url))
            response.raise_for_status()
            
            with open(cif_path, 'wb') as f:
                f.write(response.content)
            
            return str(cif_path)
        
        return None
    
    async def _get_structure_for_scf(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[str]:
        """为自洽场计算获取结构文件"""
        
        if params.get('optimized_task_id'):
            # 从已完成的结构优化任务获取CONTCAR
            if progress_callback:
                await progress_callback(15, "从结构优化任务获取优化后结构...")
            
            optimized_task_id = params['optimized_task_id']
            # 构建优化任务的工作目录路径
            opt_work_dir = self.base_work_dir / optimized_task_id
            contcar_path = opt_work_dir / "CONTCAR"
            
            if not contcar_path.exists():
                raise Exception(f"优化任务 {optimized_task_id} 的CONTCAR文件不存在")
            
            # 复制CONTCAR作为POSCAR
            poscar_path = work_dir / "POSCAR"
            import shutil
            shutil.copy(str(contcar_path), str(poscar_path))
            
            # 修改第一行为计算类型
            calc_type = params.get('calc_type', 'OXC')
            with open(poscar_path, 'r') as f:
                lines = f.readlines()
            lines[0] = f"{calc_type}\n"
            with open(poscar_path, 'w') as f:
                f.writelines(lines)
            
            return str(poscar_path)
            
        elif params.get('formula'):
            # 从化学式下载CIF然后转换
            if progress_callback:
                await progress_callback(10, f"从Materials Project下载 {params['formula']}...")
            
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            
            return await self._convert_cif_to_poscar(cif_path, work_dir, params)
            
        elif params.get('cif_url'):
            # 从CIF URL下载
            if progress_callback:
                await progress_callback(15, f"从URL下载CIF: {params['cif_url']}")
            
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            
            return await self._convert_cif_to_poscar(cif_path, work_dir, params)
        
        else:
            raise Exception("必须提供 formula、cif_url 或 optimized_task_id 中的一个")
    
    async def _prepare_dos_files(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[Dict[str, str]]:
        """为态密度计算准备文件"""
        
        if params.get('scf_task_id'):
            # 从已完成的自洽场计算任务获取文件
            if progress_callback:
                await progress_callback(15, "从自洽场计算任务获取结果文件...")
            
            scf_task_id = params['scf_task_id']
            # 构建自洽场任务的工作目录路径
            scf_work_dir = self.base_work_dir / scf_task_id
            
            # 需要复制的文件列表 (按照vasp(1).py中的逻辑)
            required_files = ["POSCAR", "POTCAR", "CHG", "CHGCAR", "WAVECAR"]
            copied_files = {}
            
            import shutil
            for filename in required_files:
                src_path = scf_work_dir / filename
                dst_path = work_dir / filename
                
                if src_path.exists():
                    shutil.copy(str(src_path), str(dst_path))
                    copied_files[filename] = str(dst_path)
                    print(f"复制文件: {filename}")
                else:
                    print(f"⚠️ 文件不存在: {src_path}")
                    if filename in ["POSCAR", "POTCAR"]:  # 关键文件
                        raise Exception(f"关键文件 {filename} 不存在于SCF任务 {scf_task_id}")
            
            return copied_files
            
        elif params.get('formula'):
            # 从化学式进行完整的计算流程 (先SCF再DOS)
            if progress_callback:
                await progress_callback(10, f"从Materials Project下载 {params['formula']}...")
            
            # 这种情况需要先进行SCF计算，然后基于结果进行DOS
            # 暂时抛出异常，建议用户先完成SCF计算
            raise Exception("建议先完成自洽场计算，然后使用scf_task_id参数进行DOS计算")
            
        elif params.get('cif_url'):
            # 从CIF URL进行完整的计算流程
            raise Exception("建议先完成自洽场计算，然后使用scf_task_id参数进行DOS计算")
        
        else:
            raise Exception("必须提供 formula、cif_url 或 scf_task_id 中的一个")
    
    async def _convert_cif_to_poscar(self, cif_path: str, work_dir: Path, params: Dict[str, Any]) -> str:
        """转换CIF为POSCAR"""
        from .base import cif_to_poscar
        print(f"🔍 转换CIF为POSCAR: {cif_path}")
        # 修改POSCAR第一行为计算类型
        poscar_path = cif_to_poscar(cif_path, str(work_dir))
        
        # 读取原POSCAR内容
        with open(poscar_path, 'r') as f:
            lines = f.readlines()
        
        # 修改第一行为计算类型
        calc_type = self._get_calc_type_from_params(params)
        lines[0] = f"{calc_type}\n"
        
        # 写回文件
        with open(poscar_path, 'w') as f:
            f.writelines(lines)
        
        return poscar_path
    
    def _get_calc_type_from_params(self, params: Dict[str, Any]) -> str:
        """从参数中获取计算类型"""
        calc_type = params.get('calc_type', 'OXC')
        return calc_type
    
    async def _generate_vasp_inputs(self, work_dir: Path, params: Dict[str, Any]):
        """生成VASP输入文件"""
        from .base import generate_kpoints, generate_potcar, generate_incar
        
        calc_type = self._get_calc_type_from_params(params)
        kpoint_density = params.get('kpoint_density', 30.0)
        
        # 生成KPOINTS
        generate_kpoints(str(work_dir))
        
        # 生成POTCAR
        generate_potcar(str(work_dir))
        
        # 生成INCAR
        generate_incar(str(work_dir), calc_type)
    
    async def _generate_scf_inputs(self, work_dir: Path, params: Dict[str, Any]):
        """生成自洽场VASP输入文件"""
        from .base import generate_kpoints, generate_potcar, generate_incar
        
        calc_type = self._get_calc_type_from_params(params)
        precision = params.get('precision', 'Accurate')
        
        # 生成KPOINTS (自洽场计算通常使用更密的K点网格)
        generate_kpoints(str(work_dir))
        
        # 生成POTCAR
        generate_potcar(str(work_dir))
        
        # 生成基础INCAR
        generate_incar(str(work_dir), calc_type)
        
        # 修改INCAR为自洽场计算设置
        await self._modify_incar_for_scf(work_dir, precision)
    
    async def _modify_incar_for_scf(self, work_dir: Path, precision: str):
        """修改INCAR文件用于自洽场计算"""
        incar_path = work_dir / "INCAR"
        
        # 读取原INCAR
        with open(incar_path, 'r') as f:
            lines = f.readlines()
        
        new_lines = []
        scf_settings_added = False
        
        for line in lines:
            stripped = line.strip().upper()
            
            # 修改基础设置
            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = SCF\n")
            elif stripped.startswith("PREC"):
                new_lines.append(f"PREC = {precision}\n")
            elif stripped.startswith("NSW"):
                new_lines.append("NSW = 0\n")  # 自洽场不优化结构
            elif stripped.startswith("ISIF"):
                new_lines.append("ISIF = 2\n")
            elif stripped.startswith("IBRION"):
                new_lines.append("IBRION = -1\n")  # 不做离子步
            elif stripped.startswith("POTIM"):
                new_lines.append("# POTIM = 0\n")  # 注释掉
            elif stripped.startswith("EDIFFG"):
                new_lines.append("# EDIFFG = -0.01\n")  # 注释掉
            elif stripped.startswith("LWAVE"):
                new_lines.append("LWAVE = .TRUE.\n")  # 保存波函数
            elif stripped.startswith("LCHARG"):
                new_lines.append("LCHARG = .TRUE.\n")  # 保存电荷密度
            else:
                new_lines.append(line)
        
        # 添加自洽场专用设置
        if not scf_settings_added:
            new_lines.append("\n# 自洽场计算专用设置\n")
            new_lines.append("EDIFF = 1E-6\n")  # 更严格的电子收敛
            new_lines.append("NELMIN = 4\n")   # 最小电子步数
            new_lines.append("NELM = 200\n")   # 更多电子步数
            new_lines.append("ISMEAR = 0\n")   # Gaussian展宽
            new_lines.append("SIGMA = 0.05\n") # 展宽参数
        
        # 写回文件
        with open(incar_path, 'w') as f:
            f.writelines(new_lines)
        
        print(f"自洽场INCAR已生成于 {incar_path}")
    
    async def _generate_dos_inputs(self, work_dir: Path, params: Dict[str, Any], dos_files: Dict[str, str]):
        """生成态密度VASP输入文件"""
        
        # DOS计算不需要重新生成POSCAR和POTCAR，直接使用从SCF复制的文件
        
        # 1. 修改INCAR文件用于DOS计算
        await self._modify_incar_for_dos(work_dir, params)
        
        # 2. 生成DOS专用KPOINTS（基于优化计算倍增）
        await self._generate_dos_kpoints(work_dir, params)
    
    async def _modify_incar_for_dos(self, work_dir: Path, params: Dict[str, Any]):
        """修改INCAR文件用于态密度计算"""
        
        # 查找源INCAR文件
        scf_task_id = params.get('scf_task_id')
        if scf_task_id:
            # 从SCF任务复制INCAR
            scf_work_dir = self.base_work_dir / scf_task_id
            scf_incar_path = scf_work_dir / "INCAR"
            
            if not scf_incar_path.exists():
                raise Exception(f"SCF任务 {scf_task_id} 的INCAR文件不存在")
            
            # 读取SCF的INCAR
            with open(scf_incar_path, 'r') as f:
                lines = f.readlines()
        else:
            # 生成基础INCAR
            from .base import generate_incar
            calc_type = self._get_calc_type_from_params(params)
            generate_incar(str(work_dir), calc_type)
            
            # 读取生成的INCAR
            with open(work_dir / "INCAR", 'r') as f:
                lines = f.readlines()
        
        new_lines = []
        
        for line in lines:
            stripped = line.strip().upper()
            
            # 修改DOS计算专用设置
            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = DOS\n")
            elif stripped.startswith("NSW"):
                new_lines.append("NSW = 0\n")  # DOS不进行离子步
            elif stripped.startswith("IBRION"):
                new_lines.append("IBRION = -1\n")  # 不做离子运动
            elif stripped.startswith("ICHARG"):
                new_lines.append("ICHARG = 11\n")  # 从CHGCAR读取电荷密度
            elif stripped.startswith("ISMEAR"):
                new_lines.append("ISMEAR = -5\n")  # 四面体方法
            elif stripped.startswith("SIGMA"):
                new_lines.append("# SIGMA = 0.05\n")  # 注释掉，四面体方法不需要
            elif stripped.startswith("LWAVE"):
                new_lines.append("LWAVE = .FALSE.\n")  # DOS计算不需要保存波函数
            elif stripped.startswith("LCHARG"):
                new_lines.append("LCHARG = .FALSE.\n")  # DOS计算不需要保存电荷密度
            else:
                new_lines.append(line)
        
        # 添加DOS专用设置
        new_lines.append("\n# 态密度计算专用设置\n")
        new_lines.append("LORBIT = 11\n")  # 计算轨道分辨态密度
        new_lines.append("NEDOS = 2000\n")  # 能量网格点数
        new_lines.append("EMIN = -20\n")   # 能量范围最小值
        new_lines.append("EMAX = 10\n")    # 能量范围最大值
        
        # 写入DOS INCAR文件
        dos_incar_path = work_dir / "INCAR"
        with open(dos_incar_path, 'w') as f:
            f.writelines(new_lines)
        
        print(f"DOS INCAR已生成于 {dos_incar_path}")
    
    async def _generate_dos_kpoints(self, work_dir: Path, params: Dict[str, Any]):
        """生成DOS计算的KPOINTS文件"""
        
        # 获取K点倍增因子
        kpoint_multiplier = params.get('kpoint_multiplier', 2.0)
        
        # 查找优化任务的KPOINTS作为基础
        scf_task_id = params.get('scf_task_id')
        if scf_task_id:
            scf_work_dir = self.base_work_dir / scf_task_id
            
            # 尝试找到对应的优化任务KPOINTS
            # 按照vasp(1).py的逻辑，需要从"1-opt"目录获取KPOINTS
            # 这里简化处理，直接从SCF任务目录获取KPOINTS，然后倍增
            scf_kpoints_path = scf_work_dir / "KPOINTS"
            
            if scf_kpoints_path.exists():
                # 读取原KPOINTS
                with open(scf_kpoints_path, 'r') as f:
                    lines = f.readlines()
                
                # 修改网格密度（第4行）
                if len(lines) >= 4:
                    # 解析网格
                    grid_line = lines[3].strip().split()
                    if len(grid_line) >= 3:
                        try:
                            nx, ny, nz = map(int, grid_line[:3])
                            # 应用倍增因子
                            new_nx = max(1, int(nx * kpoint_multiplier))
                            new_ny = max(1, int(ny * kpoint_multiplier))
                            new_nz = max(1, int(nz * kpoint_multiplier))
                            
                            lines[3] = f"{new_nx} {new_ny} {new_nz}\n"
                            
                            # 写入DOS KPOINTS
                            dos_kpoints_path = work_dir / "KPOINTS"
                            with open(dos_kpoints_path, 'w') as f:
                                f.writelines(lines)
                            
                            print(f"DOS KPOINTS已生成: {new_nx}x{new_ny}x{new_nz} (倍增因子: {kpoint_multiplier})")
                            return
                        except ValueError:
                            pass
        
        # 如果无法从原有KPOINTS倍增，则生成新的
        from .base import generate_kpoints
        generate_kpoints(str(work_dir))
        print("已生成默认DOS KPOINTS")
    
    async def _run_vasp_calculation(self, work_dir: Path, progress_callback=None) -> Dict[str, Any]:
        """运行VASP计算"""
        start_time = time.time()
        
        # 提交作业
        if progress_callback:
            await progress_callback(35, "提交VASP作业...")
        
        # 目前使用直接运行的方式
        vasp_path = get_path_config()["VASP_PATH"]
        shell_command = f"""
        source /etc/profile.d/modules.sh
        module load vasp/6.3.2-intel
        srun {vasp_path}
        """
        
        try:
            # 运行VASP
            process = await asyncio.create_subprocess_shell(
                shell_command,
                cwd=str(work_dir),
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            pid = process.pid
            print(f"🔍 VASP进程ID: {pid}")
            
            # 通过特殊回调传递PID
            if progress_callback:
                await progress_callback(35, f"VASP进程已启动，PID: {pid}", pid=pid)
            
            # 等待计算完成，同时更新进度
            progress = 36
            while process.returncode is None:
                if progress_callback:
                    await progress_callback(min(progress, 90), "VASP计算进行中...")
                
                await asyncio.sleep(10)  # 每10秒检查一次
                progress = min(progress + 2, 90)
                
                # 检查进程状态
                try:
                    await asyncio.wait_for(process.wait(), timeout=1.0)
                except asyncio.TimeoutError:
                    continue
            
            stdout, stderr = await process.communicate()
            
            end_time = time.time()
            computation_time = end_time - start_time
            
            if process.returncode != 0:
                error_msg = f"VASP执行失败，返回码: {process.returncode}\n"
                error_msg += f"标准错误输出: {stderr.decode()}"
                raise Exception(error_msg)
            
            return {
                'success': True,
                'computation_time': computation_time,
                'stdout': stdout.decode(),
                'stderr': stderr.decode(),
                'process_id': pid
            }
            
        except Exception as e:
            raise Exception(f"VASP计算执行失败: {str(e)}")
    
    def _create_slurm_job(self,num_nodes=2, total_tasks=56, tasks_per_node=28, partition="normal3", cmd="srun /path/to/vasp_std"):
        script = f"""#!/bin/bash
    #SBATCH -N {num_nodes}
    #SBATCH -n {total_tasks}
    #SBATCH --ntasks-per-node={tasks_per_node}
    #SBATCH --partition={partition}
    #SBATCH --output=%j.out
    #SBATCH --error=%j.err

    {cmd}
    """
        return script


    
    async def _analyze_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析计算结果"""
        try:
            # 检查收敛性
            outcar_path = work_dir / "OUTCAR"
            convergence = self._check_convergence(outcar_path)
            
            # 提取能量
            energy = self._extract_energy(outcar_path)
            
            # 提取力
            forces = self._extract_forces(outcar_path)
            
            # 复制优化后的结构
            contcar_path = work_dir / "CONTCAR"
            optimized_structure = None
            if contcar_path.exists():
                optimized_structure = str(contcar_path)
            
            return {
                'success': True,
                'convergence': convergence,
                'energy': energy,
                'final_forces': forces,
                'optimized_structure': optimized_structure,
                'computation_time': vasp_result.get('computation_time'),
                'process_id': vasp_result.get('process_id'),
                'work_directory': str(work_dir)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f"结果分析失败: {str(e)}",
                'work_directory': str(work_dir)
            }
    
    def _check_convergence(self, outcar_path: Path) -> bool:
        """检查计算是否收敛"""
        try:
            with open(outcar_path, 'rb') as f:
                f.seek(-1024, os.SEEK_END)
                last_lines = f.readlines()[-10:]
                last_content = b''.join(last_lines).decode('utf-8', errors='ignore')
                return 'reached required accuracy' in last_content or 'Voluntary' in last_content
        except Exception:
            return False
    
    def _extract_energy(self, outcar_path: Path) -> Optional[float]:
        """从OUTCAR提取最终能量"""
        try:
            with open(outcar_path, 'r') as f:
                lines = f.readlines()
            
            for line in reversed(lines):
                if 'free energy    TOTEN' in line:
                    parts = line.split()
                    return float(parts[4])
            return None
        except Exception:
            return None
    
    def _extract_forces(self, outcar_path: Path) -> Optional[list]:
        """从OUTCAR提取最终力矩阵"""
        try:
            with open(outcar_path, 'r') as f:
                lines = f.readlines()
            
            forces = []
            reading_forces = False
            
            for line in reversed(lines):
                if 'TOTAL-FORCE' in line:
                    reading_forces = True
                    continue
                
                if reading_forces:
                    if line.strip() and not line.startswith('-'):
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                force = [float(parts[3]), float(parts[4]), float(parts[5])]
                                forces.insert(0, force)
                            except (ValueError, IndexError):
                                break
                    elif line.startswith('-'):
                        break
                        
            return forces if forces else None
        except Exception:
            return None
    
    async def _analyze_scf_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析自洽场计算结果"""
        try:
            # 检查收敛性
            outcar_path = work_dir / "OUTCAR"
            convergence = self._check_convergence(outcar_path)
            
            # 提取总能量
            total_energy = self._extract_energy(outcar_path)
            
            # 提取费米能级
            fermi_energy = self._extract_fermi_energy(outcar_path)
            
            # 提取带隙
            band_gap = self._extract_band_gap(outcar_path)
            
            # 提取电子步数
            electronic_steps = self._extract_electronic_steps(outcar_path)
            
            # SCF结构文件路径
            poscar_path = work_dir / "POSCAR"
            scf_structure = str(poscar_path) if poscar_path.exists() else None
            
            return {
                'success': True,
                'convergence': convergence,
                'total_energy': total_energy,
                'fermi_energy': fermi_energy,
                'band_gap': band_gap,
                'electronic_steps': electronic_steps,
                'scf_structure': scf_structure,
                'computation_time': vasp_result.get('computation_time'),
                'process_id': vasp_result.get('process_id'),
                'work_directory': str(work_dir)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f"自洽场结果分析失败: {str(e)}",
                'work_directory': str(work_dir)
            }
    
    def _extract_fermi_energy(self, outcar_path: Path) -> Optional[float]:
        """从OUTCAR提取费米能级"""
        try:
            with open(outcar_path, 'r') as f:
                lines = f.readlines()
            
            for line in reversed(lines):
                if 'E-fermi' in line:
                    # 格式: E-fermi :   -1.2345     XC(G=0):  -10.2345     alpha+bet : -11.2345
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'E-fermi' and i + 2 < len(parts):
                            return float(parts[i + 2])
            return None
        except Exception:
            return None
    
    def _extract_band_gap(self, outcar_path: Path) -> Optional[float]:
        """从OUTCAR提取带隙"""
        try:
            with open(outcar_path, 'r') as f:
                content = f.read()
            
            # 在OUTCAR中查找带隙信息
            if 'band gap' in content.lower():
                lines = content.split('\n')
                for line in lines:
                    if 'band gap' in line.lower():
                        # 尝试提取数值
                        import re
                        numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                        if numbers:
                            return float(numbers[0])
            
            # 如果没有直接的带隙信息，返回None（可能是金属）
            return None
        except Exception:
            return None
    
    def _extract_electronic_steps(self, outcar_path: Path) -> Optional[int]:
        """从OUTCAR提取电子步数"""
        try:
            with open(outcar_path, 'r') as f:
                lines = f.readlines()
            
            electronic_steps = 0
            for line in lines:
                if 'RMM:' in line or 'DAV:' in line:
                    electronic_steps += 1
            
            return electronic_steps if electronic_steps > 0 else None
        except Exception:
            return None
    
    async def _analyze_dos_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析态密度计算结果"""
        try:
            # 检查收敛性
            outcar_path = work_dir / "OUTCAR"
            convergence = self._check_convergence(outcar_path)
            
            # 提取总能量
            total_energy = self._extract_energy(outcar_path)
            
            # 提取费米能级
            fermi_energy = self._extract_fermi_energy(outcar_path)
            
            # 提取带隙
            band_gap = self._extract_band_gap(outcar_path)
            
            # 检查DOSCAR文件
            doscar_path = work_dir / "DOSCAR"
            doscar_exists = doscar_path.exists()
            
            # 提取DOS数据
            dos_data = None
            kpoints_used = None
            if doscar_exists:
                dos_data = self._extract_dos_data(doscar_path)
                kpoints_used = self._extract_kpoints_info(work_dir / "KPOINTS")
            
            # DOS结构文件路径
            poscar_path = work_dir / "POSCAR"
            dos_structure = str(poscar_path) if poscar_path.exists() else None
            
            return {
                'success': True,
                'convergence': convergence,
                'total_energy': total_energy,
                'fermi_energy': fermi_energy,
                'band_gap': band_gap,
                'dos_structure': dos_structure,
                'doscar_path': str(doscar_path) if doscar_exists else None,
                'dos_data': dos_data,
                'kpoints_used': kpoints_used,
                'computation_time': vasp_result.get('computation_time'),
                'process_id': vasp_result.get('process_id'),
                'work_directory': str(work_dir)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f"态密度结果分析失败: {str(e)}",
                'work_directory': str(work_dir)
            }
    
    def _extract_dos_data(self, doscar_path: Path) -> Optional[dict]:
        """从DOSCAR文件提取态密度数据"""
        try:
            with open(doscar_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 6:
                return None
            
            # 读取DOSCAR头部信息
            natoms = int(lines[0].split()[0])
            fermi_line = lines[5].split()
            fermi_energy = float(fermi_line[3])
            
            # 提取总DOS数据 (从第7行开始)
            total_dos = {
                'energy': [],
                'dos_total': [],
                'dos_integrated': []
            }
            
            # 查找总DOS数据结束位置
            start_line = 6
            end_line = start_line
            for i in range(start_line, len(lines)):
                line = lines[i].strip()
                if not line or line.startswith('#'):
                    end_line = i
                    break
                try:
                    data = line.split()
                    if len(data) >= 3:
                        total_dos['energy'].append(float(data[0]))
                        total_dos['dos_total'].append(float(data[1]))
                        total_dos['dos_integrated'].append(float(data[2]))
                    end_line = i + 1
                except (ValueError, IndexError):
                    end_line = i
                    break
            
            # 返回处理后的数据
            result = {
                'natoms': natoms,
                'fermi_energy': fermi_energy,
                'total_dos': total_dos,
                'data_points': len(total_dos['energy'])
            }
            
            return result
            
        except Exception as e:
            print(f"提取DOS数据失败: {str(e)}")
            return None
    
    def _extract_kpoints_info(self, kpoints_path: Path) -> Optional[list]:
        """从KPOINTS文件提取K点网格信息"""
        try:
            with open(kpoints_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) >= 4:
                grid_line = lines[3].strip().split()
                if len(grid_line) >= 3:
                    return [int(x) for x in grid_line[:3]]
            
            return None
        except Exception:
            return None 