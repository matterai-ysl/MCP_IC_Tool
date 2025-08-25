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
            
            # 2. 生成DOS计算输入文件
            if progress_callback:
                if params.get('scf_task_id'):
                    await progress_callback(30, "生成态密度VASP输入文件...")
                    await self._generate_dos_inputs(work_dir, params, dos_files)
                else:
                    await progress_callback(25, "单点自洽+DOS输入文件已准备完成")
            else:
                if params.get('scf_task_id'):
                    await self._generate_dos_inputs(work_dir, params, dos_files)
            
            # 3. 运行VASP计算
            if progress_callback:
                if params.get('scf_task_id'):
                    await progress_callback(40, "开始VASP态密度计算...")
                else:
                    await progress_callback(30, "开始单点自洽+DOS计算...")
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
    
    async def run_md_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """
        运行分子动力学计算
        
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
                await progress_callback(5, "开始分子动力学计算...")
            
            # 1. 获取结构文件和准备文件
            md_files = await self._prepare_md_files(work_dir, params, progress_callback)
            if not md_files:
                raise Exception("无法准备MD计算文件")
            
            # 2. 生成MD计算输入文件
            if progress_callback:
                if params.get('scf_task_id'):
                    await progress_callback(30, "生成分子动力学VASP输入文件...")
                    await self._generate_md_inputs(work_dir, params, md_files)
                else:
                    await progress_callback(25, "纯MD输入文件已准备完成")
            else:
                if params.get('scf_task_id'):
                    await self._generate_md_inputs(work_dir, params, md_files)
            
            # 3. 运行VASP计算
            if progress_callback:
                if params.get('scf_task_id'):
                    await progress_callback(40, "开始VASP分子动力学计算...")
                else:
                    await progress_callback(30, "开始纯MD计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)
            
            # 4. 分析MD结果
            if progress_callback:
                await progress_callback(90, "分析分子动力学计算结果...")
            final_result = await self._analyze_md_results(work_dir, result)
            
            if progress_callback:
                await progress_callback(100, "分子动力学计算完成！")
                
            return final_result
            
        except Exception as e:
            error_msg = f"分子动力学计算失败: {str(e)}"
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
            # 从化学式进行单点自洽+DOS计算（一步完成）
            if progress_callback:
                await progress_callback(10, f"从Materials Project下载 {params['formula']}...")
            
            # 获取CIF并转换为POSCAR
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            poscar_path = await self._convert_cif_to_poscar(cif_path, work_dir, params)
            
            # 生成单点自洽+DOS的输入文件
            if progress_callback:
                await progress_callback(20, "准备单点自洽+DOS计算文件...")
            await self._prepare_single_point_dos_files(work_dir, params)
            
            return {"POSCAR": str(poscar_path)}
            
        elif params.get('cif_url'):
            # 从CIF URL进行单点自洽+DOS计算（一步完成）
            if progress_callback:
                await progress_callback(10, f"从URL下载CIF: {params['cif_url']}")
            
            # 获取CIF并转换为POSCAR
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            poscar_path = await self._convert_cif_to_poscar(cif_path, work_dir, params)
            
            # 生成单点自洽+DOS的输入文件
            if progress_callback:
                await progress_callback(20, "准备单点自洽+DOS计算文件...")
            await self._prepare_single_point_dos_files(work_dir, params)
            
            return {"POSCAR": str(poscar_path)}
        
        else:
            raise Exception("必须提供 formula、cif_url 或 scf_task_id 中的一个")
    
    async def _prepare_md_files(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[Dict[str, str]]:
        """为分子动力学计算准备文件"""
        
        if params.get('scf_task_id'):
            # 从已完成的自洽场计算任务获取文件
            if progress_callback:
                await progress_callback(15, "从自洽场计算任务获取结果文件...")
            
            scf_task_id = params['scf_task_id']
            scf_work_dir = self.base_work_dir / scf_task_id
            
            # MD计算只需要POSCAR和POTCAR (按照vasp(1).py的逻辑)
            required_files = ["POSCAR", "POTCAR"]
            copied_files = {}
            
            import shutil
            for filename in required_files:
                src_path = scf_work_dir / filename
                dst_path = work_dir / filename
                
                if src_path.exists():
                    shutil.copy(str(src_path), str(dst_path))
                    copied_files[filename] = str(dst_path)
                    print(f"复制MD文件: {filename}")
                else:
                    print(f"⚠️ 文件不存在: {src_path}")
                    raise Exception(f"关键文件 {filename} 不存在于SCF任务 {scf_task_id}")
            
            return copied_files
            
        elif params.get('formula'):
            # 从化学式进行纯MD计算（一步完成）
            if progress_callback:
                await progress_callback(10, f"从Materials Project下载 {params['formula']}...")
            
            # 获取CIF并转换为POSCAR
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            poscar_path = await self._convert_cif_to_poscar(cif_path, work_dir, params)
            
            # 生成纯MD的输入文件
            if progress_callback:
                await progress_callback(20, "准备纯MD计算文件...")
            await self._prepare_single_point_md_files(work_dir, params)
            
            return {"POSCAR": str(poscar_path)}
            
        elif params.get('cif_url'):
            # 从CIF URL进行纯MD计算（一步完成）
            if progress_callback:
                await progress_callback(10, f"从URL下载CIF: {params['cif_url']}")
            
            # 获取CIF并转换为POSCAR
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            poscar_path = await self._convert_cif_to_poscar(cif_path, work_dir, params)
            
            # 生成纯MD的输入文件
            if progress_callback:
                await progress_callback(20, "准备纯MD计算文件...")
            await self._prepare_single_point_md_files(work_dir, params)
            
            return {"POSCAR": str(poscar_path)}
        
        else:
            raise Exception("必须提供 formula、cif_url 或 scf_task_id 中的一个")
    
    async def _prepare_single_point_md_files(self, work_dir: Path, params: Dict[str, Any]):
        """准备纯MD计算的输入文件"""
        from .base import generate_potcar
        
        # 1. 生成POTCAR
        generate_potcar(str(work_dir))
        
        # 2. 生成固定的MD KPOINTS (1 1 1)
        await self._generate_md_kpoints(work_dir)
        
        # 3. 生成纯MD的INCAR
        await self._generate_single_point_md_incar(work_dir, params)
        
        print("纯MD输入文件已准备完成")
    
    async def _generate_md_inputs(self, work_dir: Path, params: Dict[str, Any], md_files: Dict[str, str]):
        """生成分子动力学VASP输入文件"""
        
        # 1. 生成固定的MD KPOINTS (1 1 1)
        await self._generate_md_kpoints(work_dir)
        
        # 2. 生成MD专用INCAR
        await self._generate_md_incar(work_dir, params)
    
    async def _generate_md_kpoints(self, work_dir: Path):
        """生成MD计算的固定KPOINTS (1 1 1)"""
        kpoints_path = work_dir / "KPOINTS"
        
        # MD计算使用固定的1x1x1 K点网格
        kpoints_content = """Automatic mesh
        0
        Gamma
        1 1 1
        0.0 0.0 0.0
        """
        
        with open(kpoints_path, 'w') as f:
            f.write(kpoints_content)
        
        print("MD KPOINTS已生成: 1x1x1 (固定)")
    
    async def _generate_md_incar(self, work_dir: Path, params: Dict[str, Any]):
        """生成分子动力学的INCAR文件"""
        
        # MD计算的基础INCAR内容（基于vasp(1).py的MD_INCAR_CONTENT）
        md_steps = params.get('md_steps', 1000)
        temperature = params.get('temperature', 300.0)
        time_step = params.get('time_step', 1.0)
        ensemble = params.get('ensemble', 'NVT')
        precision = params.get('precision', 'Normal')
        calc_type = self._get_calc_type_from_params(params)
        
        incar_content = f"""SYSTEM = MD-{calc_type}
        PREC = {precision}
        ISMEAR = 0
        SIGMA = 0.1
        IBRION = 0
        NSW = {md_steps}
        POTIM = {time_step}
        TEBEG = {temperature}
        TEEND = {temperature}
        SMASS = 0
        NBLOCK = 1
        ISYM = 0
        LCHARG = .FALSE.
        LWAVE = .FALSE.
        """

        # 根据系综类型添加特定设置
        if ensemble.upper() == 'NVT':
            incar_content += """
        # NVT系综设置
        MDALGO = 2
        ANDERSEN_PROB = 0.1
        """
        elif ensemble.upper() == 'NVE':
            incar_content += """
        # NVE系综设置  
        MDALGO = 1
        """
        elif ensemble.upper() == 'NPT':
            incar_content += """
        # NPT系综设置
        MDALGO = 3
        PSTRESS = 0.0
        LANGEVIN_GAMMA = 10.0
        """
        
        # 写入INCAR文件
        incar_path = work_dir / "INCAR"
        with open(incar_path, 'w') as f:
            f.write(incar_content.strip())
        
        print(f"MD INCAR已生成于 {incar_path} ({ensemble}系综, {md_steps}步, {temperature}K)")
    
    async def _generate_single_point_md_incar(self, work_dir: Path, params: Dict[str, Any]):
        """生成纯MD的INCAR文件（直接进行分子动力学计算，无需自洽场）"""
        
        # 直接调用纯MD的INCAR生成方法
        await self._generate_md_incar(work_dir, params)
        print("纯MD INCAR已生成（无需自洽场计算）")
    
    async def _analyze_md_results(self, work_dir: Path, run_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析分子动力学计算结果"""
        
        result = {
            "md_structure": None,
            "xdatcar_path": None,
            "oszicar_path": None,
            "final_energy": None,
            "average_temperature": None,
            "total_md_steps": None,
            "convergence": False,
            "computation_time": run_result.get("computation_time"),
            "trajectory_data": None,
            "error_message": run_result.get("error_message")
        }
        
        try:
            # 检查VASP计算是否成功
            if not run_result.get("success", False):
                result["error_message"] = run_result.get("error_message", "VASP计算失败")
                return result
            
            # 1. 检查POSCAR文件（初始结构）
            poscar_path = work_dir / "POSCAR"
            if poscar_path.exists():
                result["md_structure"] = str(poscar_path)
                print("✅ 找到初始结构文件: POSCAR")
            
            # 2. 检查XDATCAR文件（轨迹文件）
            xdatcar_path = work_dir / "XDATCAR"
            if xdatcar_path.exists():
                result["xdatcar_path"] = str(xdatcar_path)
                print("✅ 找到轨迹文件: XDATCAR")
                
                # 分析轨迹数据
                try:
                    trajectory_data = await self._extract_trajectory_data(xdatcar_path)
                    result["trajectory_data"] = trajectory_data
                    result["total_md_steps"] = trajectory_data.get("total_steps", 0)
                except Exception as e:
                    print(f"⚠️ 分析轨迹数据失败: {e}")
            
            # 3. 检查OSZICAR文件（能量和温度信息）
            oszicar_path = work_dir / "OSZICAR"
            if oszicar_path.exists():
                result["oszicar_path"] = str(oszicar_path)
                print("✅ 找到能量文件: OSZICAR")
                
                # 分析能量和温度数据
                try:
                    energy_temp_data = await self._extract_energy_temperature_data(oszicar_path)
                    result["final_energy"] = energy_temp_data.get("final_energy")
                    result["average_temperature"] = energy_temp_data.get("average_temperature")
                except Exception as e:
                    print(f"⚠️ 分析能量温度数据失败: {e}")
            
            # 4. 检查OUTCAR文件获取更多信息
            outcar_path = work_dir / "OUTCAR"
            if outcar_path.exists():
                print("✅ 找到输出文件: OUTCAR")
                try:
                    # 检查计算是否正常完成
                    with open(outcar_path, 'r') as f:
                        outcar_content = f.read()
                        if "General timing and accounting informations for this job:" in outcar_content:
                            result["convergence"] = True
                            print("✅ MD计算正常完成")
                        else:
                            print("⚠️ MD计算可能未正常完成")
                except Exception as e:
                    print(f"⚠️ 分析OUTCAR失败: {e}")
            
            # 总结结果
            if result["convergence"]:
                print(f"🎉 MD计算成功完成!")
                if result["total_md_steps"]:
                    print(f"   完成步数: {result['total_md_steps']}")
                if result["average_temperature"]:
                    print(f"   平均温度: {result['average_temperature']:.2f} K")
                if result["final_energy"]:
                    print(f"   最终能量: {result['final_energy']:.6f} eV")
            else:
                print("❌ MD计算未能正常完成")
            
            return result
            
        except Exception as e:
            error_msg = f"分析MD结果失败: {str(e)}"
            print(f"[ERROR] {error_msg}")
            result["error_message"] = error_msg
            return result
    
    async def _extract_trajectory_data(self, xdatcar_path: Path) -> Dict[str, Any]:
        """提取轨迹数据统计信息"""
        
        trajectory_data = {
            "total_steps": 0,
            "lattice_parameters": [],
            "volume_data": [],
            "step_intervals": []
        }
        
        try:
            with open(xdatcar_path, 'r') as f:
                lines = f.readlines()
            
            step_count = 0
            current_step = 0
            
            for i, line in enumerate(lines):
                stripped = line.strip()
                
                # 检查是否是新的MD步
                if stripped.startswith("Direct configuration="):
                    step_count += 1
                    # 提取步数信息
                    parts = stripped.split()
                    if len(parts) >= 2:
                        try:
                            current_step = int(parts[2])
                            trajectory_data["step_intervals"].append(current_step)
                        except (ValueError, IndexError):
                            pass
            
            trajectory_data["total_steps"] = step_count
            
            print(f"轨迹分析: 共 {step_count} 个MD步")
            
        except Exception as e:
            print(f"提取轨迹数据失败: {e}")
            raise
        
        return trajectory_data
    
    async def _extract_energy_temperature_data(self, oszicar_path: Path) -> Dict[str, Any]:
        """提取能量和温度数据"""
        
        energy_temp_data = {
            "final_energy": None,
            "average_temperature": None,
            "energy_series": [],
            "temperature_series": []
        }
        
        try:
            with open(oszicar_path, 'r') as f:
                lines = f.readlines()
            
            energies = []
            temperatures = []
            
            for line in lines:
                stripped = line.strip()
                
                # 解析MD步的能量和温度信息
                # OSZICAR格式: DAV:   1    -0.123456E+02    -0.12345E-02   -0.123E-03  1234   0.123E-01    0.123E+02
                if 'DAV:' in stripped or 'RMM:' in stripped:
                    parts = stripped.split()
                    if len(parts) >= 3:
                        try:
                            # 第三列通常是总能量
                            energy = float(parts[2])
                            energies.append(energy)
                        except (ValueError, IndexError):
                            pass
                
                # 查找温度信息 (T= 或 Temperature=)
                if 'T=' in stripped:
                    parts = stripped.split('T=')
                    if len(parts) > 1:
                        temp_part = parts[1].split()[0]
                        try:
                            temperature = float(temp_part)
                            temperatures.append(temperature)
                        except ValueError:
                            pass
            
            # 计算统计数据
            if energies:
                energy_temp_data["final_energy"] = energies[-1]
                energy_temp_data["energy_series"] = energies[-min(100, len(energies)):]  # 保存最后100个数据点
                print(f"能量分析: 最终能量 = {energies[-1]:.6f} eV")
            
            if temperatures:
                energy_temp_data["average_temperature"] = sum(temperatures) / len(temperatures)
                energy_temp_data["temperature_series"] = temperatures[-min(100, len(temperatures)):]  # 保存最后100个数据点
                print(f"温度分析: 平均温度 = {energy_temp_data['average_temperature']:.2f} K")
            
        except Exception as e:
            print(f"提取能量温度数据失败: {e}")
            raise
        
        return energy_temp_data
    
    async def _prepare_single_point_dos_files(self, work_dir: Path, params: Dict[str, Any]):
        """准备单点自洽+DOS计算的输入文件"""
        from .base import generate_kpoints, generate_potcar
        
        # 1. 生成KPOINTS (DOS计算使用更密的网格)
        generate_kpoints(str(work_dir))
        kpoint_multiplier = params.get('kpoint_multiplier', 2.0)
        await self._apply_kpoint_multiplier(work_dir, kpoint_multiplier)
        
        # 2. 生成POTCAR
        generate_potcar(str(work_dir))
        
        # 3. 生成单点自洽+DOS的INCAR
        await self._generate_single_point_dos_incar(work_dir, params)
        
        print("单点自洽+DOS输入文件已准备完成")
    
    async def _apply_kpoint_multiplier(self, work_dir: Path, multiplier: float):
        """应用K点倍增因子"""
        kpoints_path = work_dir / "KPOINTS"
        
        if kpoints_path.exists():
            with open(kpoints_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) >= 4:
                grid_line = lines[3].strip().split()
                if len(grid_line) >= 3:
                    try:
                        nx, ny, nz = map(int, grid_line[:3])
                        new_nx = max(1, int(nx * multiplier))
                        new_ny = max(1, int(ny * multiplier))
                        new_nz = max(1, int(nz * multiplier))
                        
                        lines[3] = f"{new_nx} {new_ny} {new_nz}\n"
                        
                        with open(kpoints_path, 'w') as f:
                            f.writelines(lines)
                        
                        print(f"K点网格已调整: {new_nx}x{new_ny}x{new_nz} (倍增: {multiplier})")
                    except ValueError:
                        pass
    
    async def _generate_single_point_dos_incar(self, work_dir: Path, params: Dict[str, Any]):
        """生成单点自洽+DOS的INCAR文件"""
        from .base import generate_incar
        
        # 获取计算类型和精度
        calc_type = self._get_calc_type_from_params(params)
        precision = params.get('precision', 'Accurate')
        
        # 先生成基础INCAR
        generate_incar(str(work_dir), calc_type)
        
        # 读取生成的INCAR
        incar_path = work_dir / "INCAR"
        with open(incar_path, 'r') as f:
            lines = f.readlines()
        
        new_lines = []
        
        for line in lines:
            stripped = line.strip().upper()
            
            # 修改基础设置
            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = Single-point SCF+DOS\n")
            elif stripped.startswith("PREC"):
                new_lines.append(f"PREC = {precision}\n")
            elif stripped.startswith("NSW"):
                new_lines.append("NSW = 0\n")  # 单点计算
            elif stripped.startswith("IBRION"):
                new_lines.append("IBRION = -1\n")  # 不做离子运动
            elif stripped.startswith("LWAVE"):
                new_lines.append("LWAVE = .TRUE.\n")  # 保存波函数
            elif stripped.startswith("LCHARG"):
                new_lines.append("LCHARG = .TRUE.\n")  # 保存电荷密度
            elif stripped.startswith("ISMEAR"):
                new_lines.append("ISMEAR = -5\n")  # DOS计算推荐四面体方法
            elif stripped.startswith("SIGMA"):
                new_lines.append("# SIGMA = 0.05\n")  # 四面体方法不需要
            else:
                new_lines.append(line)
        
        # 添加DOS专用设置
        new_lines.append("\n# 自洽场设置\n")
        new_lines.append("EDIFF = 1E-6\n")    # 严格的电子收敛
        new_lines.append("NELMIN = 4\n")      # 最小电子步数
        new_lines.append("NELM = 200\n")      # 更多电子步数
        
        new_lines.append("\n# 态密度计算设置\n")
        new_lines.append("LORBIT = 11\n")     # 轨道分辨态密度
        new_lines.append("NEDOS = 2000\n")    # 能量网格点数
        new_lines.append("EMIN = -20\n")      # 能量范围最小值
        new_lines.append("EMAX = 10\n")       # 能量范围最大值
        
        # 写入INCAR文件
        with open(incar_path, 'w') as f:
            f.writelines(new_lines)
        
        print(f"单点自洽+DOS INCAR已生成于 {incar_path}")
    
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
            new_lines.append("LAECHG = .TRUE.\n") # 保存电荷密度
            new_lines.append("LELF = .TRUE.\n") 
            new_lines.append("LORBIT = 11\n") 
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
        import re
        start_time = time.time()
        
        # 提交作业
        if progress_callback:
            await progress_callback(35, "提交VASP作业...")
        
        # SLURM作业调度参数
        nodes = 2                    # 节点数
        total_tasks = 56             # 总任务数
        tasks_per_node = 28          # 每节点任务数
        
        script = f"""#!/bin/bash
#SBATCH --job-name={work_dir.name}
#SBATCH --partition=p1
#SBATCH -N {nodes}
#SBATCH -n {total_tasks}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load vasp/6.3.2-intel
source /data/app/intel/oneapi-2023.2/setvars.sh >/dev/null 2>&1
ulimit -s unlimited
ulimit -l unlimited

echo "=== 作业信息 ==="
echo "作业ID: $SLURM_JOB_ID"
echo "分区: $SLURM_JOB_PARTITION"
echo "节点数: $SLURM_JOB_NUM_NODES"
echo "总任务数: $SLURM_NPROCS"
echo "每节点任务数: $SLURM_NTASKS_PER_NODE"
echo "节点列表: $SLURM_JOB_NODELIST"

echo "=== 开始VASP计算 ==="
mpirun -np $SLURM_NPROCS vasp_std>result.log 2>&1
echo "VASP计算完成
        """

        # 使用.sh扩展名
        script_file = work_dir / "vasp_job.sh"
        with open(script_file, "w") as f:
            f.write(script)

        try:
            # 提交SLURM作业
            submit_process = await asyncio.create_subprocess_shell(
                f"sbatch {script_file.name}",
                cwd=str(work_dir),
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            # 等待提交结果
            submit_stdout, submit_stderr = await submit_process.communicate()
            
            if submit_process.returncode != 0:
                error_msg = f"SLURM作业提交失败，返回码: {submit_process.returncode}\n"
                error_msg += f"错误信息: {submit_stderr.decode()}"
                raise Exception(error_msg)
            
            # 解析SLURM作业ID
            output = submit_stdout.decode().strip()
            print(f"✅ 作业提交成功: {output}")
            
            job_match = re.search(r'(\d+)', output)
            if not job_match:
                raise Exception(f"无法解析SLURM作业ID: {output}")
            
            slurm_job_id = job_match.group(1)
            print(f"🆔 SLURM作业ID: {slurm_job_id}")
            
            # 通过回调传递作业ID
            if progress_callback:
                await progress_callback(40, f"VASP作业已提交，作业ID: {slurm_job_id}", pid=slurm_job_id)
            
            # 监控作业状态
            progress = 45
            job_completed = False
            
            while not job_completed:
                # 检查作业状态
                status_process = await asyncio.create_subprocess_shell(
                    f"squeue -j {slurm_job_id} --noheader --format='%T'",
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE
                )
                
                status_stdout, status_stderr = await status_process.communicate()
                
                if status_process.returncode == 0:
                    status = status_stdout.decode().strip()
                    
                    if status == "":
                        # 作业不在队列中，可能已完成
                        job_completed = True
                        print("✅ 作业已完成（不在队列中）")
                    elif status in ["COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"]:
                        job_completed = True
                        print(f"✅ 作业状态: {status}")
                        
                        if status != "COMPLETED":
                            # 检查错误日志
                            error_files = list(work_dir.glob("*.err"))
                            error_msg = f"作业以状态 {status} 结束"
                            if error_files:
                                try:
                                    with open(error_files[0], 'r') as f:
                                        error_content = f.read()
                                    if error_content.strip():
                                        error_msg += f"\n错误日志:\n{error_content}"
                                except:
                                    pass
                            raise Exception(error_msg)
                    else:
                        # 作业仍在运行
                        status_msg = {
                            "PENDING": "排队中",
                            "RUNNING": "计算中", 
                            "CONFIGURING": "配置中"
                        }.get(status, f"状态: {status}")
                        
                        if progress_callback:
                            await progress_callback(min(progress, 90), f"VASP{status_msg}...")
                        
                        print(f"🔄 作业状态: {status}")
                else:
                    # 查询失败，可能作业已完成
                    print("⚠️  无法查询作业状态，检查是否已完成")
                    job_completed = True
                
                if not job_completed:
                    await asyncio.sleep(30)  # 每30秒检查一次
                    progress = min(progress + 3, 90)
            
            # 检查输出文件
            outcar_file = work_dir / "OUTCAR"
            if not outcar_file.exists():
                # 查找输出文件
                output_files = list(work_dir.glob("*.out"))
                error_msg = "VASP计算可能失败，未找到OUTCAR文件"
                if output_files:
                    try:
                        with open(output_files[0], 'r') as f:
                            output_content = f.read()
                        error_msg += f"\n作业输出:\n{output_content}"
                    except:
                        pass
                raise Exception(error_msg)
            
            # 读取结果
            result_log = ""
            result_log_file = work_dir / "result.log"
            if result_log_file.exists():
                with open(result_log_file, 'r') as f:
                    result_log = f.read()
            
            end_time = time.time()
            computation_time = end_time - start_time
            
            return {
                'success': True,
                'computation_time': computation_time,
                'stdout': result_log,
                'stderr': '',
                'slurm_job_id': slurm_job_id,
                'output_files': [str(f) for f in work_dir.glob("*") if f.is_file()]
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
            
            # 生成可视化分析报告（仅对结构优化任务）
            html_report_path = None
            analysis_data = None
            try:
                from .optimization_analyzer import generate_optimization_report, OUTCARAnalyzer
                if outcar_path.exists():
                    # 生成分析数据
                    analyzer = OUTCARAnalyzer(str(work_dir), task_id="optimization")
                    analysis_data = analyzer.analyze()
                    
                    # 生成HTML报告
                    html_report_path = generate_optimization_report(str(work_dir), "optimization")
                    print(f"📊 结构优化分析报告已生成: {html_report_path}")
            except Exception as e:
                print(f"⚠️ 生成可视化分析报告失败: {e}")
            
            result = {
                'success': True,
                'convergence': convergence,
                'energy': energy,
                'final_forces': forces,
                'optimized_structure': optimized_structure,
                'computation_time': vasp_result.get('computation_time'),
                'process_id': vasp_result.get('process_id'),
                'work_directory': str(work_dir)
            }
            
            # 如果生成了HTML报告，添加到结果中
            if html_report_path:
                result['html_analysis_report'] = html_report_path
            
            # 如果生成了分析数据，添加到结果中
            if analysis_data:
                result['analysis_data'] = analysis_data
            
            return result
            
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

    def _run_bader_analysis(self, work_dir: Path):
        """运行Bader电荷分析"""
        CHGSUM_PL_PATH = get_path_config()["vasp_path"]["chgsum_pl"]
        BADER_PATH = get_path_config()["vasp_path"]["bader"]
        for f in ["AECCAR0", "AECCAR2", "CHGCAR"]:
            if not os.path.exists(os.path.join(work_dir, f)): 
                raise Exception("  - 错误: 未找到Bader分析所需文件 {}。".format(f))
        chgsum_cmd = ["perl", CHGSUM_PL_PATH, "AECCAR0", "AECCAR2"]
        result = subprocess.run(
            chgsum_cmd, cwd=work_dir, stdout=subprocess.STDOUT, stderr=subprocess.STDOUT, text=True, check=True
        )
        if result.returncode != 0:
            raise Exception("  - 错误: 生成 CHGCAR_sum 文件失败。")
        if not os.path.exists(os.path.join(work_dir, "CHGCAR_sum")): 
            raise Exception("  - 错误: 未生成 CHGCAR_sum 文件。")
        bader_cmd = [BADER_PATH, "CHGCAR", "-ref", "CHGCAR_sum"]
        result = subprocess.run(
            bader_cmd, cwd=work_dir, stdout=subprocess.STDOUT, stderr=subprocess.STDOUT, text=True, check=True
        )
        if result.returncode != 0:
            raise Exception("  - 错误: Bader分析失败。")
        return True
    
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
            # 运行Bader电荷分析
            self._run_bader_analysis(work_dir)
            # 生成可视化分析报告（使用新的SCF分析器）
            html_report_path = None
            analysis_data = None
            try:
                from .scf_analyzer import generate_scf_report, SCFAnalyzer
                if outcar_path.exists():
                    # 生成分析数据
                    analyzer = SCFAnalyzer(str(work_dir), task_id="scf")
                    analysis_data = analyzer.analyze()
                    
                    # 生成HTML报告
                    html_report_path = generate_scf_report(str(work_dir), "scf")
                    print(f"📊 SCF计算分析报告已生成: {html_report_path}")
            except Exception as e:
                print(f"⚠️ 生成SCF可视化分析报告失败: {e}")
            
            # SCF结构文件路径
            poscar_path = work_dir / "POSCAR"
            scf_structure = str(poscar_path) if poscar_path.exists() else None
            
            result = {
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
            
            # 如果生成了HTML报告，添加到结果中
            if html_report_path:
                result['html_analysis_report'] = html_report_path
            
            # 如果生成了分析数据，添加到结果中
            if analysis_data:
                result['analysis_data'] = analysis_data
            
            return result
            
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
            
            # 生成可视化分析报告（使用SCF分析器，因为DOS计算也包含SCF过程）
            html_report_path = None
            analysis_data = None
            try:
                from .scf_analyzer import generate_scf_report, SCFAnalyzer
                if outcar_path.exists():
                    # 生成分析数据
                    analyzer = SCFAnalyzer(str(work_dir), task_id="dos")
                    analysis_data = analyzer.analyze()
                    
                    # 生成HTML报告
                    html_report_path = generate_scf_report(str(work_dir), "dos")
                    print(f"📊 DOS计算分析报告已生成: {html_report_path}")
            except Exception as e:
                print(f"⚠️ 生成DOS可视化分析报告失败: {e}")
            
            # DOS结构文件路径
            poscar_path = work_dir / "POSCAR"
            dos_structure = str(poscar_path) if poscar_path.exists() else None
            
            result = {
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
            
            # 如果生成了HTML报告，添加到结果中
            if html_report_path:
                result['html_analysis_report'] = html_report_path
            
            # 如果生成了分析数据，添加到结果中
            if analysis_data:
                result['analysis_data'] = analysis_data
            
            return result
            
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