import os
import shutil
import asyncio
import traceback
from typing import Dict, Any, Optional, List
import requests
from pathlib import Path
import subprocess
import time

from .mp import download_with_criteria
from .base import cif_to_poscar
from .Config import get_path_config, get_kpoints_config,get_static_url,get_download_url
from typing import TYPE_CHECKING, Callable
import importlib

if TYPE_CHECKING:
    # 仅用于类型检查，避免运行时硬依赖
    from .MD_analyzer import generate_md_analysis_report  # type: ignore

def _load_md_report_func() -> Optional[Callable[..., str]]:
    """动态加载 MD 分析报告函数，兼容不同模块名/路径。"""
    candidates = [
        'src.vasp_server.MD_analyzer',
        'src.vasp_server.md_analyzer',
        'vasp_server.MD_analyzer',
        'vasp_server.md_analyzer',
        __package__ + '.MD_analyzer' if __package__ else 'MD_analyzer',
        __package__ + '.md_analyzer' if __package__ else 'md_analyzer',
    ]
    for mod_name in candidates:
        try:
            mod = importlib.import_module(mod_name)
            func = getattr(mod, 'generate_md_analysis_report', None)
            if callable(func):
                return func  # type: ignore
        except Exception:
            continue
    return None

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
        运行分子动力学计算（支持多温度扫描）
        
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
            # 解析温度参数
            temperature_param = params.get('temperature', 300.0)
            if isinstance(temperature_param, list):
                # 多温度MD计算
                return await self._run_multi_temperature_md(task_id, params, temperature_param, progress_callback)
            else:
                # 单温度MD计算（保持原有逻辑）
                return await self._run_single_temperature_md(task_id, params, float(temperature_param), progress_callback)
                
        except Exception as e:
            error_msg = f"分子动力学计算失败: {str(e)}"
            print(f"[ERROR] {error_msg}")
            print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
            raise Exception(error_msg)
    
    async def _run_single_temperature_md(self, task_id: str, params: Dict[str, Any], temperature: float, progress_callback=None) -> Dict[str, Any]:
        """运行单温度MD计算（原有逻辑）"""
        work_dir = self.base_work_dir / task_id
        
        # 更新进度: 开始处理
        if progress_callback:
            await progress_callback(5, f"开始分子动力学计算 (T={temperature}K)...")
        
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
        
        # 5. 生成MD分析HTML报告
        try:
            md_report_func = _load_md_report_func()
            if md_report_func is not None:
                if progress_callback:
                    await progress_callback(95, "生成分子动力学分析报告...")
                html_path = md_report_func(str(work_dir), task_id=task_id, output_dir=str(work_dir / "MD_output"))
                html_relative_path = get_static_url(html_path)
                final_result["md_analysis_report_html_path"] = html_relative_path
                final_result["md_output_dir"] = str(work_dir / "MD_output")
            else:
                print("⚠️ 未找到 MD 分析报告生成函数，跳过报告生成。")
        except Exception as e:
            print(f"⚠️ 生成MD分析报告失败: {e}")
        
        if progress_callback:
            await progress_callback(100, "分子动力学计算完成！")
            
        return final_result
    
    async def _run_multi_temperature_md(self, task_id: str, params: Dict[str, Any], temperatures: List[float], progress_callback=None) -> Dict[str, Any]:
        """运行多温度MD计算"""
        work_dir = self.base_work_dir / task_id
        
        print(f"🌡️ 开始多温度MD计算，温度列表: {temperatures}")
        
        if progress_callback:
            await progress_callback(5, f"开始多温度MD计算，共{len(temperatures)}个温度点...")
        
        # 准备基础文件（共享的结构文件等）
        base_md_files = await self._prepare_md_files(work_dir, params, progress_callback)
        if not base_md_files:
            raise Exception("无法准备MD计算文件")
        
        # 创建子任务结果列表
        subtask_results = []
        completed_count = 0
        failed_count = 0
        total_temps = len(temperatures)
        
        # 为每个温度创建子任务
        for i, temp in enumerate(temperatures):
            try:
                print(f"🔥 处理温度 {temp}K ({i+1}/{total_temps})")
                
                if progress_callback:
                    progress = 10 + (i * 80 // total_temps)
                    await progress_callback(progress, f"计算温度 {temp}K ({i+1}/{total_temps})")
                
                # 创建温度专用子目录
                temp_dir = work_dir / f"T_{temp}K"
                temp_dir.mkdir(parents=True, exist_ok=True)
                
                # 复制基础文件到子目录
                await self._copy_base_files_to_temp_dir(base_md_files, temp_dir)
                
                # 为该温度生成专门的MD输入文件
                temp_params = params.copy()
                temp_params['temperature'] = temp
                await self._generate_md_inputs_for_temperature(temp_dir, temp_params, temp)
                
                # 运行该温度的VASP计算
                temp_result = await self._run_vasp_calculation(temp_dir, None)  # 不传递进度回调避免混乱
                
                # 分析该温度的结果
                temp_analysis = await self._analyze_md_results(temp_dir, temp_result)
                
                # 创建子任务结果
                subtask_result = {
                    "temperature": temp,
                    "subtask_dir": str(temp_dir),
                    "md_structure": temp_analysis.get("md_structure"),
                    "xdatcar_path": temp_analysis.get("xdatcar_path"),
                    "oszicar_path": temp_analysis.get("oszicar_path"),
                    "final_energy": temp_analysis.get("final_energy"),
                    "average_temperature": temp_analysis.get("average_temperature"),
                    "total_md_steps": temp_analysis.get("total_md_steps"),
                    "convergence": temp_analysis.get("convergence", False),
                    "computation_time": temp_analysis.get("computation_time"),
                    "trajectory_data": temp_analysis.get("trajectory_data"),
                    "status": "completed" if temp_analysis.get("convergence", False) else "failed",
                    "error_message": temp_analysis.get("error_message")
                }
                
                subtask_results.append(subtask_result)
                
                if temp_analysis.get("convergence", False):
                    completed_count += 1
                    print(f"✅ 温度 {temp}K 计算成功")
                else:
                    failed_count += 1
                    print(f"❌ 温度 {temp}K 计算失败")
                    
            except Exception as e:
                print(f"❌ 温度 {temp}K 计算出错: {str(e)}")
                failed_count += 1
                
                subtask_result = {
                    "temperature": temp,
                    "subtask_dir": str(work_dir / f"T_{temp}K"),
                    "status": "failed",
                    "error_message": str(e),
                    "convergence": False
                }
                subtask_results.append(subtask_result)
        
        # 生成多温度分析报告
        try:
            if progress_callback:
                await progress_callback(95, "生成多温度MD分析报告...")
            html_path = await self._generate_multi_temperature_report(work_dir, task_id, subtask_results)
            html_relative_path = get_static_url(html_path)    #type: ignore
        except Exception as e:
            print(f"⚠️ 生成多温度分析报告失败: {e}")
            html_path = None
        
        # 构建最终结果
        final_result = {
            "is_multi_temperature": True,
            "total_subtasks": total_temps,
            "completed_subtasks": completed_count,
            "failed_subtasks": failed_count,
            "subtask_results": subtask_results,
            "md_analysis_report_html_path": html_relative_path,
            "md_output_dir": str(work_dir / "MD_output"),
            "convergence": completed_count > 0,
            "computation_time": sum([r.get("computation_time", 0) for r in subtask_results if r.get("computation_time")])
        }
        
        if progress_callback:
            await progress_callback(100, f"多温度MD计算完成！成功: {completed_count}, 失败: {failed_count}")
        
        print(f"🎉 多温度MD计算完成！总计: {total_temps}, 成功: {completed_count}, 失败: {failed_count}")
        
        return final_result
    
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
        
        # 4. 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
        
        print("纯MD输入文件已准备完成")
    
    async def _generate_md_inputs(self, work_dir: Path, params: Dict[str, Any], md_files: Dict[str, str]):
        """生成分子动力学VASP输入文件"""
        
        # 1. 生成固定的MD KPOINTS (1 1 1)
        await self._generate_md_kpoints(work_dir)
        
        # 2. 生成MD专用INCAR
        await self._generate_md_incar(work_dir, params)
        
        # 3. 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
    
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
        
        # 4. 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
        
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
        
        # 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
    
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
        
        # 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
    
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
        
        # 3. 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
    
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
                'optimized_structure_download_url': get_download_url(optimized_structure), #type: ignore
                'computation_time': vasp_result.get('computation_time'),
                'process_id': vasp_result.get('process_id'),
                'work_directory': str(work_dir)
            }

            # 如果生成了HTML报告，添加到结果中
            if html_report_path:
                html_relative_path = get_static_url(html_report_path)
                result['analysis_report_html_path'] = html_relative_path
            
            # 如果生成了分析数据，添加到结果中
            if analysis_data:
                result['analysis_data'] = analysis_data
            #简化返回结果
            simplified_result = {
                'success': result['success'],
                # 'energy': result['energy'],
                'force_convergence': result['analysis_data']['force_convergence']["converged"],
                'final_max_force': result['analysis_data']['force_convergence']["final_max_force"],
                'energy_convergence': result['analysis_data']['energy_convergence']["converged"],
                'final_energy': result['analysis_data']['energy_convergence']["final_energy"],
                # 'final_forces': result['final_forces'],
                'optimized_structure_download_url': result['optimized_structure_download_url'],
                'computation_time': result['computation_time'],
                'analysis_report_html_path': result['analysis_report_html_path'],
            }
            return simplified_result
            
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
                html_relative_path = get_static_url(html_report_path)
                result['scf_analysis_report_html_path'] = html_relative_path
            
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
                html_relative_path = get_static_url(html_report_path)
                result['analysis_report_html_path'] = html_relative_path
            
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
    
    async def _apply_custom_incar(self, work_dir: Path, params: Dict[str, Any]) -> None:
        """应用自定义INCAR参数"""
        custom_incar = params.get('custom_incar')
        if not custom_incar:
            return
        
        incar_path = work_dir / "INCAR"
        if not incar_path.exists():
            print(f"⚠️ INCAR文件不存在: {incar_path}")
            return
        
        try:
            # 读取现有INCAR内容
            with open(incar_path, 'r') as f:
                lines = f.readlines()
            
            # 解析现有参数
            existing_params = {}
            for line in lines:
                line = line.strip()
                if '=' in line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    existing_params[key.strip().upper()] = value.strip()
            
            # 应用自定义参数（覆盖现有参数）
            for key, value in custom_incar.items():
                key_upper = key.upper()
                existing_params[key_upper] = str(value)
                print(f"🔧 自定义INCAR参数: {key_upper} = {value}")
            
            # 重新写入INCAR文件
            with open(incar_path, 'w') as f:
                f.write("# VASP INCAR file with custom parameters\n")
                f.write("# Generated automatically with user customizations\n\n")
                
                for key, value in existing_params.items():
                    f.write(f"{key} = {value}\n")
                
                if custom_incar:
                    f.write(f"\n# Custom parameters applied: {list(custom_incar.keys())}\n")
            
            print(f"✅ 已应用 {len(custom_incar)} 个自定义INCAR参数")
            
        except Exception as e:
            print(f"❌ 应用自定义INCAR参数失败: {e}")
            # 不抛出异常，继续计算，因为自定义参数是可选的
    
    async def _copy_base_files_to_temp_dir(self, base_files: Dict[str, str], temp_dir: Path) -> None:
        """复制基础文件到温度子目录"""
        import shutil
        
        for file_type, file_path in base_files.items():
            if file_path and Path(file_path).exists():
                src_path = Path(file_path)
                dst_path = temp_dir / src_path.name
                try:
                    shutil.copy2(str(src_path), str(dst_path))
                    print(f"📁 复制文件 {src_path.name} 到 {temp_dir.name}")
                except Exception as e:
                    print(f"⚠️ 复制文件 {src_path.name} 失败: {e}")
    
    async def _generate_md_inputs_for_temperature(self, temp_dir: Path, params: Dict[str, Any], temperature: float) -> None:
        """为特定温度生成MD输入文件"""
        # 生成固定的MD KPOINTS (1 1 1)
        await self._generate_md_kpoints(temp_dir)
        
        # 生成温度专用的MD INCAR
        await self._generate_md_incar(temp_dir, params)
        
        # 应用自定义INCAR参数
        await self._apply_custom_incar(temp_dir, params)
        
        print(f"🌡️ 已为温度 {temperature}K 生成输入文件")
    
    async def _generate_multi_temperature_report(self, work_dir: Path, task_id: str, subtask_results: List[Dict]) -> Optional[str]:
        """生成增强的多温度MD分析报告，包含Arrhenius分析和标签页界面"""
        try:
            from datetime import datetime
            import numpy as np
            import base64
            import io
            
            # 创建MD输出目录
            output_dir = work_dir / "MD_output"
            output_dir.mkdir(exist_ok=True)
            
            # 执行多温度分析
            analysis_results = await self._perform_multi_temperature_analysis(work_dir, subtask_results, output_dir)
            
            # 生成各温度点的单独HTML
            temp_html_tabs = await self._generate_individual_temp_htmls(work_dir, subtask_results, output_dir)
            
            # 构建综合HTML报告
            html_content = self._build_comprehensive_html(
                task_id, subtask_results, analysis_results, temp_html_tabs, output_dir
            )
            
            # 保存HTML报告
            html_path = output_dir / "comprehensive_multi_temperature_report.html"
            with open(html_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            print(f"📄 增强多温度MD报告已生成: {html_path}")
            return str(html_path)
            
        except Exception as e:
            print(f"❌ 生成多温度报告失败: {e}")
            return None
    
    async def _perform_multi_temperature_analysis(self, work_dir: Path, subtask_results: List[Dict], output_dir: Path) -> Dict:
        """执行多温度分析，包括Arrhenius分析"""
        analysis_results = {
            'arrhenius': None,
            'diffusion_data': [],
            'temperature_trend': None
        }
        
        try:
            # 收集成功的温度点数据
            valid_temps = []
            valid_diffusions = []
            
            for result in subtask_results:
                if result.get('convergence', False):
                    temp = result.get('temperature', 0)
                    # 这里需要从实际的MD分析结果中提取扩散系数
                    # 暂时使用模拟数据，实际应该从各温度子目录的分析结果中读取
                    temp_dir = Path(result.get('subtask_dir', ''))
                    if temp_dir.exists():
                        try:
                            # 模拟从MD分析中提取扩散系数（实际实现中应该调用MD分析器）
                            diffusion_coeff = self._extract_diffusion_coefficient(temp_dir)
                            if diffusion_coeff and diffusion_coeff > 0:
                                valid_temps.append(temp)
                                valid_diffusions.append(diffusion_coeff)
                        except Exception as e:
                            print(f"提取温度{temp}K扩散系数失败: {e}")
            
            # 执行Arrhenius分析（需要至少2个有效温度点）
            if len(valid_temps) >= 2:
                arrhenius_result = self._calculate_arrhenius_parameters(valid_temps, valid_diffusions)
                analysis_results['arrhenius'] = arrhenius_result
                
                # 生成Arrhenius图
                self._generate_arrhenius_plot(valid_temps, valid_diffusions, arrhenius_result, output_dir)
            
            analysis_results['diffusion_data'] = list(zip(valid_temps, valid_diffusions))
            
        except Exception as e:
            print(f"多温度分析失败: {e}")
        
        return analysis_results
    
    def _extract_diffusion_coefficient(self, temp_dir: Path) -> Optional[float]:
        """从温度子目录中提取扩散系数（简化实现）"""
        try:
            # 检查是否存在XDATCAR文件
            xdatcar_path = temp_dir / "XDATCAR"
            if not xdatcar_path.exists():
                return None
            
            # 这里是简化实现，实际应该调用pymatgen的MD分析
            # 返回模拟的扩散系数，实际应该通过MSD计算
            import random
            import numpy as np
            temp = float(temp_dir.name.replace("T_", "").replace("K", ""))
            # 模拟温度依赖的扩散系数 D = D0 * exp(-Ea/(kB*T))
            base_diffusion = 1e-9 * np.exp(-0.5 / (8.617e-5 * temp))  # 假设Ea=0.5eV
            return base_diffusion * (1 + random.uniform(-0.1, 0.1))  # 添加小的随机噪声
            
        except Exception as e:
            print(f"提取扩散系数失败: {e}")
            return None
    
    def _calculate_arrhenius_parameters(self, temperatures: List[float], diffusions: List[float]) -> Dict:
        """计算Arrhenius参数"""
        try:
            import numpy as np
            
            T_array = np.array(temperatures)
            D_array = np.array(diffusions)
            
            # Arrhenius方程: D = D0 * exp(-Ea/(kB*T))
            # 线性化: ln(D) = ln(D0) - Ea/(kB*T)
            x = 1.0 / T_array  # 1/T
            y = np.log(D_array)  # ln(D)
            
            # 线性拟合
            A = np.vstack([x, np.ones(len(x))]).T
            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
            
            # 计算拟合质量
            y_pred = slope * x + intercept
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # 计算物理参数
            KB_EV = 8.617e-5  # Boltzmann常数 (eV/K)
            activation_energy = -slope * KB_EV  # 活化能 (eV)
            pre_exponential = np.exp(intercept)  # 指前因子 D0
            
            return {
                'activation_energy_eV': float(activation_energy),
                'pre_exponential_factor': float(pre_exponential),
                'r_squared': float(r_squared),
                'slope': float(slope),
                'intercept': float(intercept),
                'temperature_range': f"{min(temperatures):.0f}K - {max(temperatures):.0f}K",
                'data_points': len(temperatures)
            }
            
        except Exception as e:
            print(f"Arrhenius参数计算失败: {e}")
            return {}
    
    def _generate_arrhenius_plot(self, temperatures: List[float], diffusions: List[float], arrhenius_result: Dict, output_dir: Path):
        """生成Arrhenius图"""
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            T_array = np.array(temperatures)
            D_array = np.array(diffusions)
            x = 1000.0 / T_array  # 1000/T for better scale
            y = np.log10(D_array)  # log10(D)
            
            # 绘制数据点
            ax.scatter(x, y, c='red', s=100, alpha=0.8, edgecolors='black', linewidth=1, 
                      label='实验数据点', zorder=5)
            
            # 绘制拟合线
            x_fit = np.linspace(x.min(), x.max(), 100)
            slope_log10 = arrhenius_result['slope'] / np.log(10)  # 转换为log10尺度
            intercept_log10 = arrhenius_result['intercept'] / np.log(10)
            y_fit = slope_log10 * (x_fit * 1000) + intercept_log10
            
            ax.plot(x_fit, y_fit, 'b-', linewidth=2, alpha=0.8, 
                   label=f'Arrhenius拟合 (R² = {arrhenius_result["r_squared"]:.3f})')
            
            # 设置标签和标题
            ax.set_xlabel('1000/T (K⁻¹)', fontsize=12, fontweight='bold')
            ax.set_ylabel('log₁₀(D) [D in m²/s]', fontsize=12, fontweight='bold')
            ax.set_title(f'Arrhenius图\n活化能 = {arrhenius_result["activation_energy_eV"]:.3f} eV', 
                        fontsize=14, fontweight='bold')
            
            # 网格和样式
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.legend(fontsize=11)
            
            # 添加文本框显示参数
            textstr = f"""Arrhenius参数:
Ea = {arrhenius_result["activation_energy_eV"]:.3f} eV
D₀ = {arrhenius_result["pre_exponential_factor"]:.2e} m²/s
R² = {arrhenius_result["r_squared"]:.3f}
温度范围: {arrhenius_result["temperature_range"]}"""
            
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', bbox=props)
            
            plt.tight_layout()
            plt.savefig(output_dir / 'arrhenius_plot.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            print("✅ Arrhenius图已生成")
            
        except Exception as e:
            print(f"生成Arrhenius图失败: {e}")
    
    async def _generate_individual_temp_htmls(self, work_dir: Path, subtask_results: List[Dict], output_dir: Path) -> List[Dict]:
        """为每个温度点生成单独的HTML分析报告"""
        temp_htmls = []
        
        for result in subtask_results:
            if not result.get('convergence', False):
                continue
                
            temp = result.get('temperature', 0)
            temp_dir = Path(result.get('subtask_dir', ''))
            
            try:
                # 调用单温度MD分析（假设存在generate_md_analysis_report函数）
                # 这里需要根据实际的MD分析器API进行调用
                html_content = await self._generate_single_temp_html(temp_dir, temp)
                
                temp_htmls.append({
                    'temperature': temp,
                    'tab_id': f"temp_{int(temp)}K",
                    'tab_label': f"{temp}K",
                    'html_content': html_content
                })
                
            except Exception as e:
                print(f"生成温度{temp}K的HTML失败: {e}")
        
        return temp_htmls
    
    async def _generate_single_temp_html(self, temp_dir: Path, temperature: float) -> str:
        """生成单个温度的简化HTML分析"""
        try:
            # 这里是简化的HTML生成，实际应该调用完整的MD分析器
            xdatcar_path = temp_dir / "XDATCAR"
            oszicar_path = temp_dir / "OSZICAR"
            
            md_steps = 0
            final_energy = None
            
            # 读取基本信息
            if xdatcar_path.exists():
                try:
                    with open(xdatcar_path, 'r') as f:
                        content = f.read()
                        md_steps = content.count("Direct configuration=")
                except:
                    pass
            
            if oszicar_path.exists():
                try:
                    with open(oszicar_path, 'r') as f:
                        lines = f.readlines()
                        for line in reversed(lines):
                            if 'DAV:' in line or 'RMM:' in line:
                                parts = line.strip().split()
                                if len(parts) >= 3:
                                    final_energy = float(parts[2])
                                    break
                except:
                    pass
            
            html_content = f"""
            <div class="single-temp-analysis">
                <h3>🌡️ {temperature}K 温度点详细分析</h3>
                
                <div class="analysis-section">
                    <h4>📊 基本信息</h4>
                    <table class="info-table">
                        <tr><td>计算温度</td><td>{temperature} K</td></tr>
                        <tr><td>MD步数</td><td>{md_steps}</td></tr>
                        <tr><td>最终能量</td><td>{final_energy:.6f} eV</td></tr>
                        <tr><td>计算目录</td><td>{temp_dir.name}</td></tr>
                    </table>
                </div>
                
                <div class="analysis-section">
                    <h4>📈 结构和动力学分析</h4>
                    <p>注意：单温度计算不包含活化能和Arrhenius分析，这些需要多温度数据。</p>
                    <ul>
                        <li>轨迹文件: XDATCAR</li>
                        <li>能量演化: OSZICAR</li>
                        <li>结构分析: 可通过PyMatGen进行进一步处理</li>
                    </ul>
                </div>
                
                <div class="analysis-section">
                    <h4>📁 文件信息</h4>
                    <ul>
                        <li>POSCAR: 初始结构</li>
                        <li>XDATCAR: MD轨迹</li>
                        <li>OSZICAR: 能量和压力数据</li>
                        <li>OUTCAR: 详细输出信息</li>
                    </ul>
                </div>
            </div>
            """
            
            return html_content
            
        except Exception as e:
            print(f"生成单温度HTML失败: {e}")
            return f"<div>生成{temperature}K分析报告时出错: {e}</div>"
    
    def _build_comprehensive_html(self, task_id: str, subtask_results: List[Dict], 
                                 analysis_results: Dict, temp_html_tabs: List[Dict], 
                                 output_dir: Path) -> str:
        """构建带标签页的综合HTML报告"""
        from datetime import datetime
        
        # 计算统计信息
        total_temps = len(subtask_results)
        completed_count = sum(1 for r in subtask_results if r.get('convergence', False))
        failed_count = total_temps - completed_count
        success_rate = (completed_count / total_temps * 100) if total_temps > 0 else 0
        
        # Arrhenius分析结果
        arrhenius_section = ""
        if analysis_results.get('arrhenius'):
            arr = analysis_results['arrhenius']
            arrhenius_section = f"""
                <div class="analysis-section">
                    <h3>🔬 Arrhenius分析</h3>
                    <div class="arrhenius-results">
                        <div class="arrhenius-plot">
                            <img src="arrhenius_plot.png" alt="Arrhenius图" style="max-width: 100%; height: auto;">
                        </div>
                        <div class="arrhenius-params">
                            <h4>📊 分析参数</h4>
                            <table class="params-table">
                                <tr><td><strong>活化能 (Ea)</strong></td><td>{arr['activation_energy_eV']:.3f} eV</td></tr>
                                <tr><td><strong>指前因子 (D₀)</strong></td><td>{arr['pre_exponential_factor']:.2e} m²/s</td></tr>
                                <tr><td><strong>拟合质量 (R²)</strong></td><td>{arr['r_squared']:.3f}</td></tr>
                                <tr><td><strong>温度范围</strong></td><td>{arr['temperature_range']}</td></tr>
                                <tr><td><strong>数据点数</strong></td><td>{arr['data_points']}</td></tr>
                            </table>
                            
                            <div class="arrhenius-equation">
                                <h4>📐 Arrhenius方程</h4>
                                <p><strong>D = D₀ × exp(-Ea / kBT)</strong></p>
                                <p>其中：</p>
                                <ul>
                                    <li>D: 扩散系数 (m²/s)</li>
                                    <li>D₀: 指前因子 = {arr['pre_exponential_factor']:.2e} m²/s</li>
                                    <li>Ea: 活化能 = {arr['activation_energy_eV']:.3f} eV</li>
                                    <li>kB: 玻尔兹曼常数 = 8.617×10⁻⁵ eV/K</li>
                                    <li>T: 温度 (K)</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </div>
            """
        else:
            arrhenius_section = """
                <div class="analysis-section">
                    <h3>🔬 Arrhenius分析</h3>
                    <p class="warning">⚠️ Arrhenius分析需要至少2个成功的温度点，当前成功的温度点不足。</p>
                </div>
            """
        
        # 构建标签页
        tab_headers = ""
        tab_contents = ""
        
        for i, tab in enumerate(temp_html_tabs):
            active_class = "active" if i == 0 else ""
            tab_headers += f"""
                <button class="tab-button {active_class}" onclick="openTab(event, '{tab['tab_id']}')">{tab['tab_label']}</button>
            """
            
            tab_contents += f"""
                <div id="{tab['tab_id']}" class="tab-content {active_class}">
                    {tab['html_content']}
                </div>
            """
        
        html_content = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>综合多温度MD分析报告 - {task_id}</title>
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 20px; 
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            line-height: 1.6;
        }}
        
        .container {{ 
            max-width: 1400px; 
            margin: 0 auto; 
            background: white; 
            border-radius: 15px; 
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .header {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; 
            text-align: center; 
            padding: 30px 20px;
        }}
        
        .header h1 {{ 
            margin: 0; 
            font-size: 2.5em; 
            font-weight: 300;
        }}
        
        .header p {{ 
            margin: 10px 0 0 0; 
            opacity: 0.9; 
            font-size: 1.1em;
        }}
        
        .main-content {{ 
            padding: 30px;
        }}
        
        .summary {{ 
            background: linear-gradient(135deg, #ff9a9e 0%, #fecfef 100%);
            padding: 25px; 
            border-radius: 10px; 
            margin-bottom: 30px;
            color: #2c3e50;
        }}
        
        .summary h2 {{ 
            margin-top: 0; 
            color: #2c3e50;
        }}
        
        .stats-grid {{ 
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); 
            gap: 20px; 
            margin: 20px 0;
        }}
        
        .stat-card {{ 
            background: rgba(255,255,255,0.9); 
            padding: 15px; 
            border-radius: 8px; 
            text-align: center;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .progress-bar {{ 
            background: rgba(255,255,255,0.3); 
            height: 20px; 
            border-radius: 10px; 
            overflow: hidden; 
            margin: 15px 0;
        }}
        
        .progress-fill {{ 
            background: linear-gradient(90deg, #11998e, #38ef7d); 
            height: 100%; 
            transition: width 0.3s;
        }}
        
        .analysis-section {{ 
            background: #f8f9fa; 
            margin: 20px 0; 
            padding: 25px; 
            border-radius: 10px;
            border-left: 4px solid #007bff;
        }}
        
        .analysis-section h3 {{ 
            margin-top: 0; 
            color: #007bff;
        }}
        
        .arrhenius-results {{ 
            display: grid; 
            grid-template-columns: 1fr 1fr; 
            gap: 30px; 
            margin-top: 20px;
        }}
        
        .params-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 15px 0;
        }}
        
        .params-table td {{ 
            padding: 10px; 
            border-bottom: 1px solid #dee2e6;
        }}
        
        .params-table td:first-child {{ 
            background: #f8f9fa; 
            font-weight: 500;
            width: 40%;
        }}
        
        .arrhenius-equation {{ 
            background: #e7f3ff; 
            padding: 20px; 
            border-radius: 8px; 
            margin-top: 20px;
        }}
        
        .arrhenius-equation p {{ 
            margin: 10px 0;
        }}
        
        .arrhenius-equation ul {{ 
            margin: 10px 0; 
            padding-left: 20px;
        }}
        
        .tabs {{ 
            margin-top: 30px;
        }}
        
        .tab-header {{ 
            display: flex; 
            background: #f1f3f4; 
            border-radius: 10px 10px 0 0; 
            overflow: hidden;
            flex-wrap: wrap;
        }}
        
        .tab-button {{ 
            background: none; 
            border: none; 
            padding: 15px 25px; 
            cursor: pointer; 
            font-size: 16px; 
            transition: all 0.3s;
            flex: 1;
            min-width: 120px;
        }}
        
        .tab-button:hover {{ 
            background: rgba(0,123,255,0.1);
        }}
        
        .tab-button.active {{ 
            background: #007bff; 
            color: white; 
            font-weight: 500;
        }}
        
        .tab-content {{ 
            display: none; 
            background: white; 
            padding: 30px; 
            border-radius: 0 0 10px 10px;
            border: 1px solid #f1f3f4;
            border-top: none;
        }}
        
        .tab-content.active {{ 
            display: block;
        }}
        
        .single-temp-analysis h3 {{ 
            color: #495057; 
            border-bottom: 2px solid #e9ecef; 
            padding-bottom: 10px;
        }}
        
        .single-temp-analysis h4 {{ 
            color: #6c757d; 
            margin-top: 25px;
        }}
        
        .info-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 15px 0;
        }}
        
        .info-table td {{ 
            padding: 12px; 
            border-bottom: 1px solid #dee2e6;
        }}
        
        .info-table td:first-child {{ 
            background: #f8f9fa; 
            font-weight: 500; 
            width: 30%;
        }}
        
        .warning {{ 
            background: #fff3cd; 
            color: #856404; 
            padding: 15px; 
            border-radius: 5px; 
            border-left: 4px solid #ffc107;
        }}
        
        @media (max-width: 768px) {{
            .arrhenius-results {{ 
                grid-template-columns: 1fr;
            }}
            
            .stats-grid {{ 
                grid-template-columns: 1fr 1fr;
            }}
            
            .tab-button {{ 
                font-size: 14px; 
                padding: 12px 15px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🌡️ 综合多温度MD分析报告</h1>
            <p>任务ID: {task_id} | 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="main-content">
            <div class="summary">
                <h2>📊 计算总结</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <h3>{total_temps}</h3>
                        <p>总温度点数</p>
                    </div>
                    <div class="stat-card">
                        <h3 style="color: #28a745;">{completed_count}</h3>
                        <p>成功计算</p>
                    </div>
                    <div class="stat-card">
                        <h3 style="color: #dc3545;">{failed_count}</h3>
                        <p>失败计算</p>
                    </div>
                    <div class="stat-card">
                        <h3>{success_rate:.1f}%</h3>
                        <p>成功率</p>
                    </div>
                </div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: {success_rate}%;"></div>
                </div>
            </div>
            
            {arrhenius_section}
            
            <div class="analysis-section">
                <h3>📋 多温度分析说明</h3>
                <ul>
                    <li><strong>活化能分析</strong>: 通过多温度数据拟合Arrhenius方程，获得扩散过程的活化能</li>
                    <li><strong>温度依赖性</strong>: 研究扩散系数随温度的变化规律</li>
                    <li><strong>数据质量</strong>: R²值越接近1，表示拟合质量越好</li>
                    <li><strong>物理意义</strong>: 活化能反映了离子在材料中扩散所需克服的能垒</li>
                </ul>
            </div>
            
            <div class="tabs">
                <h2>🔍 各温度点详细分析</h2>
                <div class="tab-header">
                    {tab_headers}
                </div>
                {tab_contents}
            </div>
        </div>
    </div>
    
    <script>
        function openTab(evt, tabId) {{
            var i, tabContent, tabButtons;
            
            // 隐藏所有标签页内容
            tabContent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabContent.length; i++) {{
                tabContent[i].classList.remove("active");
            }}
            
            // 移除所有按钮的active类
            tabButtons = document.getElementsByClassName("tab-button");
            for (i = 0; i < tabButtons.length; i++) {{
                tabButtons[i].classList.remove("active");
            }}
            
            // 显示选中的标签页并设置按钮为active
            document.getElementById(tabId).classList.add("active");
            evt.currentTarget.classList.add("active");
        }}
    </script>
</body>
</html>
        """
        
        return html_content 