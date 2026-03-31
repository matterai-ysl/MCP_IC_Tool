import os
import shutil
import asyncio
import traceback
import logging
from typing import Dict, Any, Optional, List
import requests
from pathlib import Path
import subprocess
import time

from .base import cif_to_poscar
from .Config import get_path_config, get_kpoints_config,get_static_url,get_download_url, DOWNLOAD_URL
from .input_resolver import UpstreamInputResolver

# 配置日志
logger = logging.getLogger(__name__)

from .analyzers.md import generate_md_analysis_report

class VaspWorker:
    """VASP计算工作器"""
    
    def __init__(self, user_id: str, base_work_dir: Optional[str] = None):
        from .Config import DOWNLOAD_URL
        self.user_id = user_id
        base = base_work_dir or DOWNLOAD_URL
        self.base_work_dir = Path(base) / user_id  # 为每个用户创建独立目录
        self.base_work_dir.mkdir(parents=True, exist_ok=True)
        self.input_resolver = UpstreamInputResolver()
    
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
                
            logger.info("Final result: %s", final_result)
            return final_result
            
        except Exception as e:
            error_msg = f"结构优化计算失败: {str(e)}"
            logger.error(error_msg)
            logger.error("详细错误信息: %s", traceback.format_exc())
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
            logger.info("自洽场VASP输入文件生成完成")
            # 3. 运行VASP自洽场计算
            if progress_callback:
                await progress_callback(40, "开始VASP自洽场计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)
            logger.info("VASP自洽场计算完成")
            logger.debug("result: %s", result)
            # 4. 分析自洽场结果
            if progress_callback:
                await progress_callback(90, "分析自洽场计算结果...")
            final_result = await self._analyze_scf_results(work_dir, result)
            logger.info("自洽场计算结果分析完成")
            logger.debug("final_result: %s", final_result)
            if progress_callback:
                await progress_callback(100, "自洽场计算完成！")
                
            return final_result
            
        except Exception as e:
            error_msg = f"自洽场计算失败: {str(e)}"
            logger.error(error_msg)
            logger.error("详细错误信息: %s", traceback.format_exc())
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
            logger.error(error_msg)
            logger.error("详细错误信息: %s", traceback.format_exc())
            raise Exception(error_msg)

    async def run_band_structure_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """
        运行能带结构计算

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
            if progress_callback:
                await progress_callback(5, "开始能带结构计算...")

            # 1. 准备文件（与DOS计算类似，需要SCF的CHGCAR）
            bs_files = await self._prepare_band_structure_files(work_dir, params, progress_callback)
            if not bs_files:
                raise Exception("无法准备能带结构计算文件")

            # 2. 生成能带结构计算输入文件
            if progress_callback:
                await progress_callback(30, "生成能带结构VASP输入文件...")
            await self._generate_band_structure_inputs(work_dir, params, bs_files)

            # 3. 运行VASP计算
            if progress_callback:
                await progress_callback(40, "开始VASP能带结构计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)

            # 4. 分析能带结构结果
            if progress_callback:
                await progress_callback(90, "分析能带结构计算结果...")
            final_result = await self._analyze_band_structure_results(work_dir, result)

            if progress_callback:
                await progress_callback(100, "能带结构计算完成!")

            return final_result

        except Exception as e:
            error_msg = f"能带结构计算失败: {str(e)}"
            logger.error(error_msg)
            logger.error("详细错误信息: %s", traceback.format_exc())
            raise Exception(error_msg)

    async def run_md_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """
        运行单温度分子动力学计算

        Args:
            task_id: 任务ID
            params: 任务参数（temperature 为单个 float）
            progress_callback: 进度回调函数

        Returns:
            Dict: 计算结果
        """
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)
        temperature = float(params.get('temperature', 300.0))
        
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

        # 确保 success 字段存在 — 修复 MD 成功判定 bug
        if 'success' not in final_result:
            final_result['success'] = final_result.get('convergence', False)

        # 5. 生成MD分析HTML报告
        try:
            if progress_callback:
                await progress_callback(95, "生成分子动力学分析报告...")
            html_path = generate_md_analysis_report(str(work_dir), task_id=task_id, output_dir=str(work_dir / "MD_output"))
            html_relative_path = get_static_url(html_path)
            final_result["md_analysis_report_html_path"] = html_relative_path
            final_result["md_output_dir"] = str(work_dir / "MD_output")
        except Exception as e:
            logger.warning(f"生成MD分析报告失败: {e}")
        
        # 确保 work_directory 存在
        final_result['work_directory'] = str(work_dir)

        if progress_callback:
            await progress_callback(100, "分子动力学计算完成！")

        return final_result

    async def _get_cif_file(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[str]:
        """获取CIF文件"""
        if params.get('cif_url'):
            # 从URL下载
            cif_url = params['cif_url']
            if progress_callback:
                await progress_callback(15, f"从URL下载结构文件: {cif_url}")
            
            cif_path = work_dir / "structure.cif"
            response = requests.get(str(cif_url), timeout=60)
            response.raise_for_status()
            
            with open(cif_path, 'wb') as f:
                f.write(response.content)
            
            return str(cif_path)
        
        return None
    
    async def _get_structure_for_scf(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[str]:
        """为自洽场计算获取结构文件"""
        upstream_artifacts = params.get("upstream_artifact_manifest") or []
        
        if params.get('optimized_task_id') and upstream_artifacts:
            if progress_callback:
                await progress_callback(15, "从上游产物清单获取优化后结构...")

            resolved = self.input_resolver.resolve_for_scf(upstream_artifacts, work_dir)
            return resolved["POSCAR"]

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
            shutil.copy(str(contcar_path), str(poscar_path))
            
            return str(poscar_path)
            
        elif params.get('cif_url'):
            # 从CIF URL下载
            if progress_callback:
                await progress_callback(15, f"从URL下载CIF: {params['cif_url']}")
            
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            
            return await self._convert_cif_to_poscar(cif_path, work_dir, params)
        
        else:
            raise Exception("必须提供 cif_url 或 optimized_task_id 中的一个")
    
    async def _prepare_dos_files(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[Dict[str, str]]:
        """为态密度计算准备文件"""
        upstream_artifacts = params.get("upstream_artifact_manifest") or []
        
        if params.get('scf_task_id') and upstream_artifacts:
            if progress_callback:
                await progress_callback(15, "从上游产物清单获取自洽场结果文件...")

            return self.input_resolver.resolve_for_dos(upstream_artifacts, work_dir)

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
            
            for filename in required_files:
                src_path = scf_work_dir / filename
                dst_path = work_dir / filename
                
                if src_path.exists():
                    shutil.copy(str(src_path), str(dst_path))
                    copied_files[filename] = str(dst_path)
                    logger.info("复制文件: %s", filename)
                else:
                    logger.warning("文件不存在: %s", src_path)
                    if filename in ["POSCAR", "POTCAR"]:  # 关键文件
                        raise Exception(f"关键文件 {filename} 不存在于SCF任务 {scf_task_id}")
            
            return copied_files
            
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
            raise Exception("必须提供 cif_url 或 scf_task_id 中的一个")
    
    async def _prepare_band_structure_files(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[Dict[str, str]]:
        """为能带结构计算准备文件"""
        upstream_artifacts = params.get("upstream_artifact_manifest") or []

        if params.get('scf_task_id') and upstream_artifacts:
            if progress_callback:
                await progress_callback(15, "从上游产物清单获取自洽场结果文件...")

            return self.input_resolver.resolve_for_band_structure(upstream_artifacts, work_dir)

        if params.get('scf_task_id'):
            # 从已完成的自洽场计算任务获取文件
            if progress_callback:
                await progress_callback(15, "从自洽场计算任务获取结果文件...")

            scf_task_id = params['scf_task_id']
            scf_work_dir = self.base_work_dir / scf_task_id

            required_files = ["POSCAR", "POTCAR", "CHG", "CHGCAR"]
            copied_files = {}

            for filename in required_files:
                src_path = scf_work_dir / filename
                dst_path = work_dir / filename

                if src_path.exists():
                    shutil.copy(str(src_path), str(dst_path))
                    copied_files[filename] = str(dst_path)
                    logger.info("复制文件: %s", filename)
                else:
                    logger.warning("文件不存在: %s", src_path)
                    if filename in ["POSCAR", "POTCAR", "CHGCAR"]:
                        raise Exception(f"关键文件 {filename} 不存在于SCF任务 {scf_task_id}")

            return copied_files

        elif params.get('cif_url'):
            if progress_callback:
                await progress_callback(10, f"从URL下载CIF: {params['cif_url']}")

            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取CIF文件")
            poscar_path = await self._convert_cif_to_poscar(cif_path, work_dir, params)

            if progress_callback:
                await progress_callback(20, "准备SCF+能带结构计算文件...")
            await self._prepare_scf_then_band_structure_files(work_dir, params)

            return {"POSCAR": str(poscar_path)}

        else:
            raise Exception("必须提供 cif_url 或 scf_task_id 中的一个")

    async def _prepare_md_files(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> Optional[Dict[str, str]]:
        """为分子动力学计算准备文件"""
        upstream_artifacts = params.get("upstream_artifact_manifest") or []

        if params.get('scf_task_id') and upstream_artifacts:
            if progress_callback:
                await progress_callback(15, "从上游产物清单获取自洽场结果文件...")

            return self.input_resolver.resolve_for_md(upstream_artifacts, work_dir)
        
        if params.get('scf_task_id'):
            # 从已完成的自洽场计算任务获取文件
            if progress_callback:
                await progress_callback(15, "从自洽场计算任务获取结果文件...")
            
            scf_task_id = params['scf_task_id']
            scf_work_dir = self.base_work_dir / scf_task_id
            
            # MD计算只需要POSCAR和POTCAR (按照vasp(1).py的逻辑)
            required_files = ["POSCAR", "POTCAR"]
            copied_files = {}
            
            for filename in required_files:
                src_path = scf_work_dir / filename
                dst_path = work_dir / filename
                
                if src_path.exists():
                    shutil.copy(str(src_path), str(dst_path))
                    copied_files[filename] = str(dst_path)
                    logger.info("复制MD文件: %s", filename)
                else:
                    logger.warning("文件不存在: %s", src_path)
                    raise Exception(f"关键文件 {filename} 不存在于SCF任务 {scf_task_id}")
            
            return copied_files
            
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
            raise Exception("必须提供 cif_url 或 scf_task_id 中的一个")
    
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
        
        logger.info("纯MD输入文件已准备完成")
    
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
        
        logger.info("MD KPOINTS已生成: 1x1x1 (固定)")
    
    async def _generate_md_incar(self, work_dir: Path, params: Dict[str, Any]):
        """生成分子动力学的INCAR文件"""
        
        # MD计算的基础INCAR内容（基于vasp(1).py的MD_INCAR_CONTENT）
        md_steps = params.get('md_steps', 1000)
        temperature = params.get('temperature', 300.0)
        time_step = params.get('time_step', 1.0)
        ensemble = params.get('ensemble', 'NVT')
        precision = params.get('precision', 'Normal')
        incar_content = f"""SYSTEM = AIMD
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
        
        logger.info("MD INCAR已生成于 %s (%s系综, %s步, %sK)", incar_path, ensemble, md_steps, temperature)
    
    async def _generate_single_point_md_incar(self, work_dir: Path, params: Dict[str, Any]):
        """生成纯MD的INCAR文件（直接进行分子动力学计算，无需自洽场）"""
        await self._generate_md_incar(work_dir, params)
    
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
                logger.info("找到初始结构文件: POSCAR")
            
            # 2. 检查XDATCAR文件（轨迹文件）
            xdatcar_path = work_dir / "XDATCAR"
            if xdatcar_path.exists():
                result["xdatcar_path"] = str(xdatcar_path)
                logger.info("找到轨迹文件: XDATCAR")
                
                # 分析轨迹数据
                try:
                    trajectory_data = await self._extract_trajectory_data(xdatcar_path)
                    result["trajectory_data"] = trajectory_data
                    result["total_md_steps"] = trajectory_data.get("total_steps", 0)
                except Exception as e:
                    logger.warning("分析轨迹数据失败: %s", e)
            
            # 3. 检查OSZICAR文件（能量和温度信息）
            oszicar_path = work_dir / "OSZICAR"
            if oszicar_path.exists():
                result["oszicar_path"] = str(oszicar_path)
                logger.info("找到能量文件: OSZICAR")
                
                # 分析能量和温度数据
                try:
                    energy_temp_data = await self._extract_energy_temperature_data(oszicar_path)
                    result["final_energy"] = energy_temp_data.get("final_energy")
                    result["average_temperature"] = energy_temp_data.get("average_temperature")
                except Exception as e:
                    logger.warning("分析能量温度数据失败: %s", e)
            
            # 4. 检查OUTCAR文件获取更多信息
            outcar_path = work_dir / "OUTCAR"
            if outcar_path.exists():
                logger.info("找到输出文件: OUTCAR")
                try:
                    # 检查计算是否正常完成
                    with open(outcar_path, 'r') as f:
                        outcar_content = f.read()
                        if "General timing and accounting informations for this job:" in outcar_content:
                            result["convergence"] = True
                            logger.info("MD计算正常完成")
                        else:
                            logger.warning("MD计算可能未正常完成")
                except Exception as e:
                    logger.warning("分析OUTCAR失败: %s", e)
            
            # 总结结果
            if result["convergence"]:
                logger.info("MD计算成功完成!")
                if result["total_md_steps"]:
                    logger.info("   完成步数: %s", result['total_md_steps'])
                if result["average_temperature"]:
                    logger.info("   平均温度: %.2f K", result['average_temperature'])
                if result["final_energy"]:
                    logger.info("   最终能量: %.6f eV", result['final_energy'])
            else:
                logger.error("MD计算未能正常完成")
            
            return result
            
        except Exception as e:
            error_msg = f"分析MD结果失败: {str(e)}"
            logger.error(error_msg)
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
            
            logger.info("轨迹分析: 共 %s 个MD步", step_count)

        except Exception as e:
            logger.error("提取轨迹数据失败: %s", e)
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
                logger.info("能量分析: 最终能量 = %.6f eV", energies[-1])
            
            if temperatures:
                energy_temp_data["average_temperature"] = sum(temperatures) / len(temperatures)
                energy_temp_data["temperature_series"] = temperatures[-min(100, len(temperatures)):]  # 保存最后100个数据点
                logger.info("温度分析: 平均温度 = %.2f K", energy_temp_data['average_temperature'])
            
        except Exception as e:
            logger.error("提取能量温度数据失败: %s", e)
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
        
        logger.info("单点自洽+DOS输入文件已准备完成")
    
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
                        
                        logger.info("K点网格已调整: %sx%sx%s (倍增: %s)", new_nx, new_ny, new_nz, multiplier)
                    except ValueError:
                        pass
    
    async def _generate_single_point_dos_incar(self, work_dir: Path, params: Dict[str, Any]):
        """生成单点自洽+DOS的INCAR文件"""
        from .base import generate_incar
        
        precision = params.get('precision', 'Accurate')

        # 先生成基础INCAR
        generate_incar(str(work_dir))
        
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
        
        logger.info("单点自洽+DOS INCAR已生成于 %s", incar_path)
    
    async def _convert_cif_to_poscar(self, cif_path: str, work_dir: Path, params: Dict[str, Any]) -> str:
        """转换CIF为POSCAR"""
        from .base import cif_to_poscar
        logger.info("转换CIF为POSCAR: %s", cif_path)
        # 修改POSCAR第一行为计算类型
        poscar_path = cif_to_poscar(cif_path, str(work_dir))
        
        # 读取原POSCAR内容
        with open(poscar_path, 'r') as f:
            lines = f.readlines()
        
        return poscar_path

    async def _generate_vasp_inputs(self, work_dir: Path, params: Dict[str, Any]):
        """生成VASP输入文件"""
        from .base import generate_kpoints, generate_potcar, generate_incar

        # 生成KPOINTS
        generate_kpoints(str(work_dir))

        # 生成POTCAR
        generate_potcar(str(work_dir))

        # 生成INCAR
        generate_incar(str(work_dir))
        
        # 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)
    
    async def _generate_scf_inputs(self, work_dir: Path, params: Dict[str, Any]):
        """生成自洽场VASP输入文件"""
        from .base import generate_kpoints, generate_potcar, generate_incar

        precision = params.get('precision', 'Accurate')

        # 生成KPOINTS (自洽场计算通常使用更密的K点网格)
        generate_kpoints(str(work_dir))

        # 生成POTCAR
        generate_potcar(str(work_dir))

        # 生成基础INCAR
        generate_incar(str(work_dir))
        
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
        
        logger.info("自洽场INCAR已生成于 %s", incar_path)
    
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
            generate_incar(str(work_dir))
            
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
        
        logger.info("DOS INCAR已生成于 %s", dos_incar_path)
    
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
                            
                            logger.info("DOS KPOINTS已生成: %sx%sx%s (倍增因子: %s)", new_nx, new_ny, new_nz, kpoint_multiplier)
                            return
                        except ValueError:
                            pass
        
        # 如果无法从原有KPOINTS倍增，则生成新的
        from .base import generate_kpoints
        generate_kpoints(str(work_dir))
        logger.info("已生成默认DOS KPOINTS")

    async def _generate_band_structure_inputs(self, work_dir: Path, params: Dict[str, Any], bs_files: Dict[str, str]):
        """生成能带结构VASP输入文件"""

        # 1. 修改INCAR文件用于能带结构计算
        await self._modify_incar_for_band_structure(work_dir, params)

        # 2. 生成高对称k路径KPOINTS
        await self._generate_band_structure_kpoints(work_dir, params)

        # 3. 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)

    async def _modify_incar_for_band_structure(self, work_dir: Path, params: Dict[str, Any]):
        """修改INCAR文件用于能带结构计算"""

        scf_task_id = params.get('scf_task_id')
        if scf_task_id:
            scf_work_dir = self.base_work_dir / scf_task_id
            scf_incar_path = scf_work_dir / "INCAR"

            if not scf_incar_path.exists():
                raise Exception(f"SCF任务 {scf_task_id} 的INCAR文件不存在")

            with open(scf_incar_path, 'r') as f:
                lines = f.readlines()
        else:
            from .base import generate_incar
            generate_incar(str(work_dir))

            with open(work_dir / "INCAR", 'r') as f:
                lines = f.readlines()

        new_lines = []

        for line in lines:
            stripped = line.strip().upper()

            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = BAND\n")
            elif stripped.startswith("NSW"):
                new_lines.append("NSW = 0\n")
            elif stripped.startswith("IBRION"):
                new_lines.append("IBRION = -1\n")
            elif stripped.startswith("ICHARG"):
                new_lines.append("ICHARG = 11\n")  # 从CHGCAR读取电荷密度
            elif stripped.startswith("ISMEAR"):
                new_lines.append("ISMEAR = 0\n")  # Gaussian展宽（能带结构不用四面体）
            elif stripped.startswith("SIGMA"):
                new_lines.append("SIGMA = 0.05\n")
            elif stripped.startswith("LWAVE"):
                new_lines.append("LWAVE = .FALSE.\n")
            elif stripped.startswith("LCHARG"):
                new_lines.append("LCHARG = .FALSE.\n")
            else:
                new_lines.append(line)

        # 添加能带结构专用设置
        new_lines.append("\n# 能带结构计算专用设置\n")
        new_lines.append("LORBIT = 11\n")  # 轨道投影

        incar_path = work_dir / "INCAR"
        with open(incar_path, 'w') as f:
            f.writelines(new_lines)

        logger.info("能带结构 INCAR 已生成于 %s", incar_path)

    async def _generate_band_structure_kpoints(self, work_dir: Path, params: Dict[str, Any]):
        """使用 pymatgen HighSymmKpath 生成高对称k路径的KPOINTS文件"""

        line_density = params.get('line_density', 20)

        poscar_path = work_dir / "POSCAR"
        if not poscar_path.exists():
            raise Exception("POSCAR文件不存在，无法生成高对称k路径")

        try:
            from pymatgen.core import Structure
            from pymatgen.symmetry.bandstructure import HighSymmKpath

            structure = Structure.from_file(str(poscar_path))
            kpath = HighSymmKpath(structure)
            kpts = kpath.kpath

            # 生成线模式KPOINTS
            kpoints_lines = ["K-Path generated by pymatgen\n"]
            kpoints_lines.append(f"{line_density}\n")
            kpoints_lines.append("Line-mode\n")
            kpoints_lines.append("Reciprocal\n")

            path_segments = kpts['path']
            kpts_coords = kpts['kpoints']

            for segment in path_segments:
                for i in range(len(segment) - 1):
                    start_label = segment[i]
                    end_label = segment[i + 1]
                    start_coord = kpts_coords[start_label]
                    end_coord = kpts_coords[end_label]

                    kpoints_lines.append(
                        f"  {start_coord[0]:.10f}  {start_coord[1]:.10f}  {start_coord[2]:.10f}  ! {start_label}\n"
                    )
                    kpoints_lines.append(
                        f"  {end_coord[0]:.10f}  {end_coord[1]:.10f}  {end_coord[2]:.10f}  ! {end_label}\n"
                    )
                    kpoints_lines.append("\n")

            kpoints_path = work_dir / "KPOINTS"
            with open(kpoints_path, 'w') as f:
                f.writelines(kpoints_lines)

            # 保存k路径信息供后续分析使用
            import json
            kpath_info = {
                "path": path_segments,
                "kpoints": {k: list(v) for k, v in kpts_coords.items()},
                "line_density": line_density,
            }
            with open(work_dir / "kpath_info.json", 'w') as f:
                json.dump(kpath_info, f, indent=2)

            logger.info("能带结构 KPOINTS 已生成 (line_density=%d, 路径段数=%d)",
                        line_density, sum(len(s) - 1 for s in path_segments))

        except ImportError:
            raise Exception("需要安装 pymatgen 来生成高对称k路径: pip install pymatgen")

    async def _prepare_scf_then_band_structure_files(self, work_dir: Path, params: Dict[str, Any]):
        """准备从头开始的SCF+能带结构计算输入文件

        当没有scf_task_id时，先做一次SCF计算得到CHGCAR，
        然后再做能带结构计算。这里生成的是能带结构的输入文件，
        ICHARG=2表示从头计算电荷密度（非自洽读取模式）。
        """
        from .base import generate_kpoints, generate_potcar, generate_incar

        # 1. 生成POTCAR
        generate_potcar(str(work_dir))

        # 2. 生成能带结构k路径
        await self._generate_band_structure_kpoints(work_dir, params)

        # 3. 生成INCAR（直接用能带结构模式，ICHARG=2从头算）
        generate_incar(str(work_dir))

        incar_path = work_dir / "INCAR"
        with open(incar_path, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            stripped = line.strip().upper()
            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = BAND\n")
            elif stripped.startswith("NSW"):
                new_lines.append("NSW = 0\n")
            elif stripped.startswith("IBRION"):
                new_lines.append("IBRION = -1\n")
            elif stripped.startswith("ISMEAR"):
                new_lines.append("ISMEAR = 0\n")
            elif stripped.startswith("LWAVE"):
                new_lines.append("LWAVE = .FALSE.\n")
            elif stripped.startswith("LCHARG"):
                new_lines.append("LCHARG = .TRUE.\n")  # 保存电荷密度
            else:
                new_lines.append(line)

        new_lines.append("\n# 能带结构计算\n")
        new_lines.append("LORBIT = 11\n")

        with open(incar_path, 'w') as f:
            f.writelines(new_lines)

        # 4. 应用自定义INCAR参数
        await self._apply_custom_incar(work_dir, params)

        logger.info("SCF+能带结构输入文件已准备完成")

    async def _run_vasp_calculation(self, work_dir: Path, progress_callback=None) -> Dict[str, Any]:
        """运行VASP计算"""
        import re
        start_time = time.time()
        
        # 提交作业
        if progress_callback:
            await progress_callback(35, "提交VASP作业...")
        
        # SLURM作业调度参数（可通过环境变量覆盖）
        nodes = int(os.getenv("SLURM_NODES", "2"))  # 节点数
        total_tasks = int(os.getenv("SLURM_NTASKS", "80"))  # 总任务数
        tasks_per_node = int(os.getenv("SLURM_TASKS_PER_NODE", "40"))  # 每节点任务数
        partition = os.getenv("SLURM_PARTITION", "p1")
        time_limit = os.getenv("SLURM_TIME_LIMIT", "24:00:00")
        module_load = os.getenv("SLURM_MODULE_LOAD", "vasp/6.3.2-intel")
        oneapi_env_script = os.getenv("ONEAPI_ENV_SCRIPT", "/data/app/intel/oneapi-2023.2/setvars.sh")
        run_line = os.getenv("VASP_SLURM_RUN_LINE", "mpirun -np $SLURM_NPROCS vasp_std>result.log 2>&1")
        
        script = f"""#!/bin/bash
#SBATCH -N {nodes}
#SBATCH -n {total_tasks}
#SBATCH --job-name={work_dir.name}
#SBATCH --partition={partition}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --time={time_limit}
#SBATCH --output=%j.out
#SBATCH --error=%j.err

{"module load " + module_load if module_load else ""}
{"source " + oneapi_env_script + " >/dev/null 2>&1" if oneapi_env_script else ""}
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
{run_line}
echo "VASP计算完成"
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
            logger.info("作业提交成功: %s", output)
            
            job_match = re.search(r'(\d+)', output)
            if not job_match:
                raise Exception(f"无法解析SLURM作业ID: {output}")
            
            slurm_job_id = job_match.group(1)
            logger.info("SLURM作业ID: %s", slurm_job_id)
            
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
                        logger.info("作业已完成（不在队列中）")
                    elif status in ["COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"]:
                        job_completed = True
                        logger.info("作业状态: %s", status)
                        
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
                        
                        logger.info("作业状态: %s", status)
                else:
                    # 查询失败，可能作业已完成
                    logger.warning("无法查询作业状态，检查是否已完成")
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
    
    def _create_slurm_job(self,num_nodes=2, total_tasks=80, tasks_per_node=40, partition="normal3", cmd="srun /path/to/vasp_std"):
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
            logger.info(f"🔍 分析计算结果: {work_dir}")
            outcar_path = work_dir / "OUTCAR"
            convergence = self._check_convergence(outcar_path)
            logger.info(f"📊 收敛性: {convergence}")

            # 提取能量
            logger.debug(f"提取能量: {outcar_path}")
            energy = self._extract_energy(outcar_path)
            logger.info(f"⚡ 能量: {energy} eV")

            # 提取力
            logger.debug(f"提取力: {outcar_path}")
            forces = self._extract_forces(outcar_path)
            logger.debug(f"💪 力: {forces}")

            # 复制优化后的结构
            logger.debug(f"复制优化后的结构: {work_dir}")
            contcar_path = work_dir / "CONTCAR"
            optimized_structure = None
            if contcar_path.exists():
                optimized_structure = str(contcar_path)
                logger.info(f"📄 优化后的结构: {optimized_structure}")

            # 生成可视化分析报告（仅对结构优化任务）
            logger.info(f"📊 生成可视化分析报告: {work_dir}")
            html_report_path = None
            analysis_data = None
            try:
                from .optimization_analyzer import generate_optimization_report, OUTCARAnalyzer
                if outcar_path.exists():
                    # 生成分析数据
                    logger.debug(f"生成分析数据: {work_dir}")
                    analyzer = OUTCARAnalyzer(str(work_dir), task_id="optimization")
                    analysis_data = analyzer.analyze()
                    logger.debug(f"分析数据长度: {len(str(analysis_data))} 字符")
                    # 生成HTML报告
                    logger.debug(f"生成HTML报告: {work_dir}")
                    html_report_path = generate_optimization_report(str(work_dir), "optimization")
                    logger.info(f"📊 结构优化分析报告已生成: {html_report_path}")
            except Exception as e:
                logger.warning(f"⚠️ 生成可视化分析报告失败: {e}")

            # 安全地生成下载URL（处理路径不在DOWNLOAD_URL下的情况）
            optimized_structure_url = None
            if optimized_structure:
                try:
                    # 检查路径是否在DOWNLOAD_URL下
                    Path(optimized_structure).relative_to(DOWNLOAD_URL)
                    optimized_structure_url = get_download_url(optimized_structure)
                except ValueError:
                    # 路径不在DOWNLOAD_URL下，使用绝对路径
                    optimized_structure_url = optimized_structure

            result = {
                'success': True,
                'convergence': convergence,
                'energy': energy,
                'final_forces': forces,
                'optimized_structure_download_url': optimized_structure_url,
                'computation_time': vasp_result.get('computation_time',None),
                'process_id': vasp_result.get('process_id',None),
                'work_directory': str(work_dir)
            }

            # 如果生成了HTML报告，添加到结果中
            if html_report_path:
                try:
                    # 检查路径是否在DOWNLOAD_URL下
                    Path(html_report_path).relative_to(DOWNLOAD_URL)
                    html_relative_path = get_static_url(html_report_path)
                    result['analysis_report_html_path'] = html_relative_path
                except ValueError:
                    # 路径不在DOWNLOAD_URL下，使用绝对路径
                    result['analysis_report_html_path'] = html_report_path
            
            # 如果生成了分析数据，添加到结果中
            if analysis_data:
                result['analysis_data'] = analysis_data

            # 简化返回结果
            if analysis_data and 'convergence_analysis' in analysis_data:
                # 使用详细分析数据
                conv_analysis = analysis_data['convergence_analysis']
                logger.info(f"📦 简化结果收敛分析: 力收敛={conv_analysis.get('force_convergence', {}).get('converged', False)}, 能量收敛={conv_analysis.get('energy_convergence', {}).get('converged', False)}")
                simplified_result = {
                    'success': result['success'],
                    'force_convergence': conv_analysis.get('force_convergence', {}).get("converged", False),
                    'final_max_force': conv_analysis.get('force_convergence', {}).get("final_max_force", None),
                    'energy_convergence': conv_analysis.get('energy_convergence', {}).get("converged", False),
                    'final_energy': conv_analysis.get('energy_convergence', {}).get("final_energy", None),
                    'optimized_structure_download_url': result['optimized_structure_download_url'],
                    'computation_time': result['computation_time'],
                    'analysis_report_html_path': result.get('analysis_report_html_path', None),
                }
                return simplified_result
            else:
                # 如果没有分析数据，返回基础结果
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
        logger.info("运行Bader电荷分析")
        CHGSUM_PL_PATH = get_path_config()["CHGSUM_PL_PATH"]
        BADER_PATH = get_path_config()["BADER_PATH"]
        logger.debug("CHGSUM_PL_PATH: %s", CHGSUM_PL_PATH)
        logger.debug("BADER_PATH: %s", BADER_PATH)
        for f in ["AECCAR0", "AECCAR2", "CHGCAR"]:
            if not os.path.exists(os.path.join(work_dir, f)):
                logger.error("未找到Bader分析所需文件: %s", f)
                raise Exception("  - 错误: 未找到Bader分析所需文件 {}。".format(f))
        chgsum_cmd = ["perl", CHGSUM_PL_PATH, "AECCAR0", "AECCAR2"]
        logger.debug("chgsum_cmd: %s", chgsum_cmd)
        result = subprocess.run(
            chgsum_cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True
        )
        logger.debug("result: %s", result)
        if result.returncode != 0:
            raise Exception("  - 错误: 生成 CHGCAR_sum 文件失败。")
        if not os.path.exists(os.path.join(work_dir, "CHGCAR_sum")): 
            raise Exception("  - 错误: 未生成 CHGCAR_sum 文件。")
        bader_cmd = [BADER_PATH, "CHGCAR", "-ref", "CHGCAR_sum"]
        result = subprocess.run(
            bader_cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True
        )
        if result.returncode != 0:
            raise Exception("  - 错误: Bader分析失败。")
        return True
    
    async def _analyze_scf_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析自洽场计算结果"""
        try:
            # 检查收敛性
            logger.debug("检查收敛性")
            outcar_path = work_dir / "OUTCAR"
            convergence = self._check_convergence(outcar_path)
            logger.debug("收敛性: %s", convergence)
            # 提取总能量
            total_energy = self._extract_energy(outcar_path)
            logger.debug("总能量: %s", total_energy)
            # 提取费米能级
            fermi_energy = self._extract_fermi_energy(outcar_path)
            logger.debug("费米能级: %s", fermi_energy)
            # 提取带隙
            band_gap = self._extract_band_gap(outcar_path)
            logger.debug("带隙: %s", band_gap)
            # 提取电子步数
            electronic_steps = self._extract_electronic_steps(outcar_path)
            logger.debug("电子步数: %s", electronic_steps)
            # 运行Bader电荷分析
            logger.info("运行Bader电荷分析")
            self._run_bader_analysis(work_dir)
            logger.info("Bader电荷分析完成")
            # 生成可视化分析报告（使用新的SCF分析器）
            logger.info("生成可视化分析报告")
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
                    logger.info("SCF计算分析报告已生成: %s", html_report_path)
            except Exception as e:
                logger.warning("生成SCF可视化分析报告失败: %s", e)
            
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
            
            # 简化返回结果
            simplified_result = {
                'success': result['success'],
                'convergence': result['convergence'],
                'total_energy': result['total_energy'],
                'fermi_energy': result['fermi_energy'],
                'band_gap': result['band_gap'],
                'electronic_steps': result['electronic_steps'],
                'scf_structure': get_download_url(result['scf_structure']),
                'scf_analysis_report_html_path': result.get('scf_analysis_report_html_path', None),
                'calculation_settings': result.get('analysis_data', {}).get('calculation_settings', None),
                'Note': "More analysis results such as Bader charge analysis, etc. are available in the analysis html report.",
            }

            return simplified_result
            
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

            # 生成可视化分析报告（使用DOS分析器）
            html_report_path = None
            analysis_data = None
            try:
                from .dos_analyzer import generate_pymatgen_dos_report, PyMatGenDOSAnalyzer
                if outcar_path.exists():
                    # 生成分析数据
                    analyzer = PyMatGenDOSAnalyzer(str(work_dir), task_id="dos")
                    analysis_data = analyzer.analyze()
                    # 生成HTML报告
                    html_report_path = generate_pymatgen_dos_report(str(work_dir), task_id="dos")
                    logger.info("DOS计算分析报告已生成: %s", html_report_path)
            except Exception as e:
                logger.warning("生成DOS可视化分析报告失败: %s", e)
            
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
            if analysis_data is not None:
                analysis_data.pop('visualizations', None)
                
                analysis_data.pop('task_info', None)

            simplified_result = {
                'success': result['success'],
                'convergence': result['convergence'],
                'total_energy': result['total_energy'],
                'fermi_energy': result['fermi_energy'],
                'band_gap': result['band_gap'],
                'dos_structure': get_download_url(result['dos_structure']),
                'dos_analysis_report_html_path': result.get('analysis_report_html_path', None),
                "analysis_data": analysis_data,
                'Note': "More analysis results and visualization are available in the analysis html report.",
            }
            return simplified_result
            
        except Exception as e:
            return {
                'success': False,
                'error': f"态密度结果分析失败: {str(e)}",
                'work_directory': str(work_dir)
            }
    
    async def _analyze_band_structure_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析能带结构计算结果"""
        try:
            outcar_path = work_dir / "OUTCAR"
            convergence = self._check_convergence(outcar_path)
            total_energy = self._extract_energy(outcar_path)
            fermi_energy = self._extract_fermi_energy(outcar_path)

            # 使用pymatgen解析能带结构
            band_gap = None
            is_direct = None
            vbm = None
            cbm = None
            kpath_info = None

            try:
                from pymatgen.io.vasp import Vasprun
                from pymatgen.electronic_structure.plotter import BSPlotter

                vasprun_path = work_dir / "vasprun.xml"
                if vasprun_path.exists():
                    vasprun = Vasprun(str(vasprun_path), parse_projected_eigen=True)
                    bs = vasprun.get_band_structure(line_mode=True)

                    band_gap = bs.get_band_gap()['energy']
                    is_direct = bs.get_band_gap()['direct']

                    if band_gap > 0:
                        vbm = bs.get_vbm()['energy']
                        cbm = bs.get_cbm()['energy']

                    logger.info("能带结构解析完成: band_gap=%.4f eV, direct=%s", band_gap, is_direct)
            except Exception as e:
                logger.warning("pymatgen 能带结构解析失败: %s", e)
                band_gap = self._extract_band_gap(outcar_path)

            # 读取k路径信息
            kpath_json = work_dir / "kpath_info.json"
            if kpath_json.exists():
                import json
                with open(kpath_json, 'r') as f:
                    kpath_info = json.load(f)

            # 生成可视化分析报告
            html_report_path = None
            analysis_data = None
            try:
                from .analyzers.band_structure import BandStructureAnalyzer, generate_band_structure_report
                if outcar_path.exists():
                    analyzer = BandStructureAnalyzer(str(work_dir))
                    analysis_data = analyzer.analyze()
                    html_report_path = generate_band_structure_report(str(work_dir))
                    logger.info("能带结构分析报告已生成: %s", html_report_path)
            except Exception as e:
                logger.warning("生成能带结构可视化分析报告失败: %s", e)

            poscar_path = work_dir / "POSCAR"
            band_structure_path = str(poscar_path) if poscar_path.exists() else None

            result = {
                'success': True,
                'convergence': convergence,
                'total_energy': total_energy,
                'fermi_energy': fermi_energy,
                'band_gap': band_gap,
                'is_direct': is_direct,
                'vbm': vbm,
                'cbm': cbm,
                'band_structure_path': band_structure_path,
                'kpath_info': kpath_info,
                'computation_time': vasp_result.get('computation_time'),
                'process_id': vasp_result.get('process_id'),
                'work_directory': str(work_dir),
            }

            if html_report_path:
                result['band_structure_report_html_path'] = get_static_url(html_report_path)

            simplified_result = {
                'success': result['success'],
                'convergence': result['convergence'],
                'total_energy': result['total_energy'],
                'fermi_energy': result['fermi_energy'],
                'band_gap': result['band_gap'],
                'is_direct': result['is_direct'],
                'vbm': result['vbm'],
                'cbm': result['cbm'],
                'band_structure_path': get_download_url(result['band_structure_path']) if result['band_structure_path'] else None,
                'band_structure_report_html_path': result.get('band_structure_report_html_path'),
                'kpath_info': result['kpath_info'],
                'analysis_data': analysis_data,
                'Note': "More analysis results and visualization are available in the analysis html report.",
            }
            return simplified_result

        except Exception as e:
            return {
                'success': False,
                'error': f"能带结构结果分析失败: {str(e)}",
                'work_directory': str(work_dir),
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
            logger.error("提取DOS数据失败: %s", str(e))
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
            logger.warning("INCAR文件不存在: %s", incar_path)
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
                logger.debug("自定义INCAR参数: %s = %s", key_upper, value)
            
            # 重新写入INCAR文件
            with open(incar_path, 'w') as f:
                f.write("# VASP INCAR file with custom parameters\n")
                f.write("# Generated automatically with user customizations\n\n")
                
                for key, value in existing_params.items():
                    f.write(f"{key} = {value}\n")
                
                if custom_incar:
                    f.write(f"\n# Custom parameters applied: {list(custom_incar.keys())}\n")
            
            logger.info("已应用 %d 个自定义INCAR参数", len(custom_incar))
            
        except Exception as e:
            logger.error("应用自定义INCAR参数失败: %s", e)
            # 不抛出异常，继续计算，因为自定义参数是可选的

    # ================================================================== #
    #  NEB 计算
    # ================================================================== #

    async def run_neb_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """运行 NEB（过渡态）计算"""
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)

        try:
            if progress_callback:
                await progress_callback(5, "准备 NEB 初始/终态结构...")

            # 1. 获取初始和终态 POSCAR
            initial_poscar, final_poscar = await self._prepare_neb_endpoints(work_dir, params, progress_callback)

            # 2. 生成中间图像
            if progress_callback:
                await progress_callback(20, "生成 NEB 中间图像...")
            await self._generate_neb_images(work_dir, initial_poscar, final_poscar, params)

            # 3. 生成 INCAR / KPOINTS / POTCAR
            if progress_callback:
                await progress_callback(30, "生成 NEB VASP 输入文件...")
            await self._generate_neb_inputs(work_dir, params)

            # 4. 运行 VASP
            if progress_callback:
                await progress_callback(40, "开始 VASP NEB 计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)

            # 5. 分析结果
            if progress_callback:
                await progress_callback(90, "分析 NEB 计算结果...")
            final_result = await self._analyze_neb_results(work_dir, result)

            if progress_callback:
                await progress_callback(100, "NEB 计算完成!")
            return final_result

        except Exception as e:
            error_msg = f"NEB 计算失败: {str(e)}"
            logger.error(error_msg)
            logger.error("详细错误: %s", traceback.format_exc())
            raise Exception(error_msg)

    async def _prepare_neb_endpoints(self, work_dir: Path, params: Dict[str, Any], progress_callback=None):
        """获取 NEB 初始结构和终态结构，返回 (initial_poscar_path, final_poscar_path)。"""
        from .base import generate_potcar

        async def _get_poscar(task_id_key: str, cif_url_key: str, manifest_key: str, dest_name: str) -> str:
            """从 task_id / cif_url 中获取 POSCAR，保存为 dest_name。"""
            task_id_val = params.get(task_id_key)
            cif_url_val = params.get(cif_url_key)
            dest = work_dir / dest_name

            if task_id_val:
                upstream_artifacts = params.get(manifest_key) or []
                if upstream_artifacts:
                    resolved = self.input_resolver.resolve_single_structure(
                        upstream_artifacts,
                        work_dir,
                        dest_name=dest_name,
                    )
                    logger.info("NEB 端点结构来自上游产物清单: %s", task_id_val)
                    return resolved[dest_name]

                src_dir = self.base_work_dir / task_id_val
                # 优先使用 CONTCAR（优化后结构），回退 POSCAR
                for fname in ("CONTCAR", "POSCAR"):
                    src = src_dir / fname
                    if src.exists():
                        shutil.copy(str(src), str(dest))
                        logger.info("NEB 端点结构来自 %s/%s", task_id_val, fname)
                        return str(dest)
                raise Exception(f"任务 {task_id_val} 中未找到 POSCAR 或 CONTCAR")

            tmp_dir = work_dir / f"_tmp_{dest_name}"
            tmp_dir.mkdir(exist_ok=True)
            tmp_params = dict(params)
            if cif_url_val:
                tmp_params['cif_url'] = cif_url_val
            else:
                raise Exception(f"必须提供 {task_id_key} 或 {cif_url_key}")

            cif_path = await self._get_cif_file(tmp_dir, tmp_params, None)
            if not cif_path:
                raise Exception(f"无法获取 {dest_name} 的 CIF 文件")
            poscar_path = await self._convert_cif_to_poscar(cif_path, tmp_dir, tmp_params)
            shutil.copy(poscar_path, str(dest))
            shutil.rmtree(str(tmp_dir), ignore_errors=True)
            return str(dest)

        initial_poscar = await _get_poscar(
            "initial_task_id",
            "initial_cif_url",
            "initial_upstream_artifact_manifest",
            "POSCAR_initial",
        )
        final_poscar = await _get_poscar(
            "final_task_id",
            "final_cif_url",
            "final_upstream_artifact_manifest",
            "POSCAR_final",
        )
        return initial_poscar, final_poscar

    async def _generate_neb_images(self, work_dir: Path, initial_poscar: str, final_poscar: str, params: Dict[str, Any]):
        """使用 pymatgen 线性插值生成中间图像，写入数字命名子目录。"""
        from pymatgen.core import Structure

        n_images = int(params.get('n_images', 5))
        initial = Structure.from_file(initial_poscar)
        final = Structure.from_file(final_poscar)

        # interpolate 返回 n_images+2 个结构（含端点）
        images = initial.interpolate(final, nimages=n_images + 1, autosort_tol=0.5)

        for i, img in enumerate(images):
            img_dir = work_dir / f"{i:02d}"
            img_dir.mkdir(exist_ok=True)
            poscar_out = img_dir / "POSCAR"
            img.to(fmt="poscar", filename=str(poscar_out))
            logger.info("NEB 图像 %02d 已写入: %s", i, poscar_out)

        logger.info("共生成 %d 个图像（含端点）", len(images))

    async def _generate_neb_inputs(self, work_dir: Path, params: Dict[str, Any]):
        """生成 NEB 计算所需的 INCAR、KPOINTS、POTCAR（POTCAR 复制到各子目录）。"""
        from .base import generate_potcar, generate_kpoints

        n_images = int(params.get('n_images', 5))

        # POSCAR 放 00/ 中，用于生成 POTCAR
        poscar_00 = work_dir / "00" / "POSCAR"
        if not poscar_00.exists():
            raise Exception("NEB 图像目录 00/POSCAR 不存在")
        shutil.copy(str(poscar_00), str(work_dir / "POSCAR"))

        # 生成 POTCAR 和 KPOINTS 到 work_dir
        generate_potcar(str(work_dir))
        generate_kpoints(str(work_dir))

        # 将 POTCAR 复制到每个图像子目录
        potcar_src = work_dir / "POTCAR"
        for i in range(n_images + 2):
            img_dir = work_dir / f"{i:02d}"
            if img_dir.exists():
                shutil.copy(str(potcar_src), str(img_dir / "POTCAR"))

        # 生成 NEB INCAR
        incar_path = work_dir / "INCAR"
        with open(incar_path, 'w') as f:
            f.write(f"""SYSTEM = NEB
ICHAIN = 0
IMAGES = {n_images}
SPRING = -5
LCLIMB = .TRUE.
IBRION = 3
NSW = 500
POTIM = 0.5
EDIFF = 1E-4
EDIFFG = -0.05
ISIF = 2
ISYM = 0
ENCUT = 520
PREC = Normal
ISMEAR = 0
SIGMA = 0.05
LWAVE = .FALSE.
LCHARG = .FALSE.
ISPIN = 2
GGA = PE
LREAL = Auto
""")
        await self._apply_custom_incar(work_dir, params)
        logger.info("NEB 输入文件已生成（IMAGES=%d）", n_images)

    async def _analyze_neb_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析 NEB 计算结果，生成 HTML 报告。"""
        try:
            html_report_path = None
            neb_data = {}
            try:
                from .analyzers.neb import NEBAnalyzer, generate_neb_report
                analyzer = NEBAnalyzer(str(work_dir), task_id="neb")
                analysis_data = analyzer.analyze()
                neb_data = analysis_data.get('neb', {})
                html_report_path = generate_neb_report(str(work_dir), task_id="neb")
                logger.info("NEB 分析报告已生成: %s", html_report_path)
            except Exception as e:
                logger.warning("生成 NEB 报告失败: %s", e)

            result = {
                'success': vasp_result.get('success', False),
                'forward_barrier_eV': neb_data.get('forward_barrier_eV'),
                'backward_barrier_eV': neb_data.get('backward_barrier_eV'),
                'reaction_energy_eV': neb_data.get('reaction_energy_eV'),
                'n_images': neb_data.get('n_images'),
                'ts_image_index': neb_data.get('ts_image_index'),
                'computation_time': vasp_result.get('computation_time'),
                'work_directory': str(work_dir),
            }
            if html_report_path:
                result['neb_report_html_path'] = get_static_url(html_report_path)
            return result

        except Exception as e:
            return {
                'success': False,
                'error': f"NEB 结果分析失败: {str(e)}",
                'work_directory': str(work_dir),
            }

    # ================================================================== #
    #  Phonon 计算
    # ================================================================== #

    async def run_phonon_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """运行声子计算（IBRION=6 有限位移+对称性，Gamma 点）"""
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)

        try:
            if progress_callback:
                await progress_callback(5, "准备声子计算结构文件...")

            # 1. 获取结构文件
            await self._prepare_phonon_files(work_dir, params, progress_callback)

            # 2. 生成 INCAR / KPOINTS
            if progress_callback:
                await progress_callback(25, "生成声子计算 VASP 输入文件...")
            await self._generate_phonon_inputs(work_dir, params)

            # 3. 运行 VASP
            if progress_callback:
                await progress_callback(35, "开始 VASP 声子计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)

            # 4. 分析结果
            if progress_callback:
                await progress_callback(90, "分析声子计算结果...")
            final_result = await self._analyze_phonon_results(work_dir, result)

            if progress_callback:
                await progress_callback(100, "声子计算完成!")
            return final_result

        except Exception as e:
            error_msg = f"声子计算失败: {str(e)}"
            logger.error(error_msg)
            logger.error("详细错误: %s", traceback.format_exc())
            raise Exception(error_msg)

    async def _prepare_phonon_files(self, work_dir: Path, params: Dict[str, Any], progress_callback=None):
        """为声子计算准备 POSCAR 和 POTCAR。"""
        from .base import generate_potcar
        upstream_artifacts = params.get("upstream_artifact_manifest") or []

        if params.get('scf_task_id') and upstream_artifacts:
            if progress_callback:
                await progress_callback(15, "从上游产物清单获取声子输入文件...")
            self.input_resolver.resolve_for_phonon(upstream_artifacts, work_dir)
            if not (work_dir / "POTCAR").exists():
                generate_potcar(str(work_dir))
            return

        if params.get('scf_task_id'):
            if progress_callback:
                await progress_callback(15, "从 SCF 任务获取结构文件...")
            scf_dir = self.base_work_dir / params['scf_task_id']
            for fname in ("POSCAR", "POTCAR"):
                src = scf_dir / fname
                if fname == "POTCAR" and not src.exists():
                    # POTCAR 由 CONTCAR/POSCAR 重新生成
                    continue
                if not src.exists():
                    raise Exception(f"SCF 任务 {params['scf_task_id']} 中缺少 {fname}")
                shutil.copy(str(src), str(work_dir / fname))
            # 优先使用优化后的 CONTCAR 作为 POSCAR
            contcar = scf_dir / "CONTCAR"
            if contcar.exists():
                shutil.copy(str(contcar), str(work_dir / "POSCAR"))
            if not (work_dir / "POTCAR").exists():
                generate_potcar(str(work_dir))

        elif params.get('cif_url'):
            if progress_callback:
                await progress_callback(10, "获取 CIF 结构文件...")
            cif_path = await self._get_cif_file(work_dir, params, progress_callback)
            if not cif_path:
                raise Exception("无法获取 CIF 文件")
            await self._convert_cif_to_poscar(cif_path, work_dir, params)
            generate_potcar(str(work_dir))
        else:
            raise Exception("必须提供 cif_url 或 scf_task_id 中的一个")

    async def _generate_phonon_inputs(self, work_dir: Path, params: Dict[str, Any]):
        """生成声子计算 INCAR 和 KPOINTS。"""
        from .base import generate_kpoints

        displacement = float(params.get('displacement', 0.015))
        kpoint_density = float(params.get('kpoint_density', 30.0))

        # KPOINTS（普通 Monkhorst-Pack，声子不需要高密度）
        generate_kpoints(str(work_dir))

        incar_path = work_dir / "INCAR"
        with open(incar_path, 'w') as f:
            f.write(f"""SYSTEM = Phonon
IBRION = 6
NSW = 1
NFREE = 2
POTIM = {displacement}
EDIFF = 1E-8
ENCUT = 520
PREC = Accurate
ISMEAR = 0
SIGMA = 0.01
LREAL = .FALSE.
ADDGRID = .TRUE.
LWAVE = .FALSE.
LCHARG = .FALSE.
ISPIN = 2
GGA = PE
""")
        await self._apply_custom_incar(work_dir, params)
        logger.info("声子计算 INCAR 已生成（displacement=%.4f Å）", displacement)

    async def _analyze_phonon_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """分析声子计算结果，生成 HTML 报告。"""
        try:
            html_report_path = None
            phonon_data = {}
            try:
                from .analyzers.phonon import PhononAnalyzer, generate_phonon_report
                analyzer = PhononAnalyzer(str(work_dir), task_id="phonon")
                analysis_data = analyzer.analyze()
                phonon_data = analysis_data.get('phonon', {})
                html_report_path = generate_phonon_report(str(work_dir), task_id="phonon")
                logger.info("声子分析报告已生成: %s", html_report_path)
            except Exception as e:
                logger.warning("生成声子报告失败: %s", e)

            result = {
                'success': vasp_result.get('success', False),
                'n_modes': phonon_data.get('n_modes'),
                'n_imaginary': phonon_data.get('n_imaginary'),
                'dynamically_stable': phonon_data.get('dynamically_stable'),
                'max_imaginary_freq_cm1': phonon_data.get('max_imaginary_freq_cm1'),
                'max_real_freq_cm1': phonon_data.get('max_real_freq_cm1'),
                'min_real_freq_cm1': phonon_data.get('min_real_freq_cm1'),
                'computation_time': vasp_result.get('computation_time'),
                'work_directory': str(work_dir),
            }
            if html_report_path:
                result['phonon_report_html_path'] = get_static_url(html_report_path)
            return result

        except Exception as e:
            return {
                'success': False,
                'error': f"声子结果分析失败: {str(e)}",
                'work_directory': str(work_dir),
            }

    # ================================================================== #
    #  通用自定义计算
    # ================================================================== #

    async def run_custom_calculation(self, task_id: str, params: Dict[str, Any], progress_callback=None) -> Dict[str, Any]:
        """运行通用自定义 VASP 计算 — 用户完全控制 INCAR"""
        work_dir = self.base_work_dir / task_id
        work_dir.mkdir(parents=True, exist_ok=True)

        try:
            if progress_callback:
                await progress_callback(5, "准备结构文件...")
            await self._prepare_custom_structure(work_dir, params, progress_callback)

            if progress_callback:
                await progress_callback(15, "生成POTCAR...")
            from .base import generate_potcar
            generate_potcar(str(work_dir))

            if progress_callback:
                await progress_callback(18, "生成KPOINTS...")
            await self._generate_custom_kpoints(work_dir, params)

            if progress_callback:
                await progress_callback(20, "写入INCAR...")
            await self._write_custom_incar(work_dir, params)

            if progress_callback:
                await progress_callback(30, "提交VASP计算...")
            result = await self._run_vasp_calculation(work_dir, progress_callback)

            if progress_callback:
                await progress_callback(90, "分析计算结果...")
            final_result = await self._analyze_custom_results(work_dir, result)
            final_result['work_directory'] = str(work_dir)

            if progress_callback:
                await progress_callback(100, "计算完成！")
            return final_result

        except Exception as e:
            error_msg = f"自定义计算失败: {str(e)}"
            logger.error(error_msg)
            logger.error("详细错误信息: %s", traceback.format_exc())
            raise Exception(error_msg)

    async def _prepare_custom_structure(self, work_dir: Path, params: Dict[str, Any], progress_callback=None) -> None:
        """为自定义计算准备 POSCAR（优先 CONTCAR）"""
        upstream_artifacts = params.get("upstream_artifact_manifest") or []

        if params.get('from_task_id') and upstream_artifacts:
            resolved = self.input_resolver.resolve_for_custom(upstream_artifacts, work_dir)
            logger.info("从上游产物清单复制 POSCAR（来自任务 %s）", params['from_task_id'])
            if progress_callback:
                await progress_callback(10, "从上游产物清单获取结构文件...")
            if resolved.get("POSCAR"):
                return

        if params.get('from_task_id'):
            src_dir = self.base_work_dir / params['from_task_id']
            for fname in ("CONTCAR", "POSCAR"):
                src = src_dir / fname
                if src.exists() and src.stat().st_size > 0:
                    shutil.copy(str(src), str(work_dir / "POSCAR"))
                    logger.info("复制 %s → POSCAR（来自任务 %s）", fname, params['from_task_id'])
                    return
            raise Exception(f"任务 {params['from_task_id']} 中未找到有效的 CONTCAR 或 POSCAR")

        cif_path = await self._get_cif_file(work_dir, params, progress_callback)
        if not cif_path:
            raise Exception("无法获取 CIF 文件")
        await self._convert_cif_to_poscar(cif_path, work_dir, params)

    async def _generate_custom_kpoints(self, work_dir: Path, params: Dict[str, Any]) -> None:
        """生成 K 点文件（mesh 或 gamma）"""
        mode = params.get('kpoint_mode', 'mesh').lower()
        if mode == 'gamma':
            kpoints_content = "Automatic mesh\n0\nGamma\n1 1 1\n0.0 0.0 0.0\n"
            with open(work_dir / "KPOINTS", 'w') as f:
                f.write(kpoints_content)
            logger.info("自定义KPOINTS: Gamma-only (1×1×1)")
        else:
            from .base import generate_kpoints
            generate_kpoints(str(work_dir))
            logger.info("自定义KPOINTS: MP网格（密度 %.1f）", params.get('kpoint_density', 30.0))

    async def _write_custom_incar(self, work_dir: Path, params: Dict[str, Any]) -> None:
        """将用户提供的 incar 字典直接写入 INCAR 文件"""
        incar_dict = dict(params.get('incar', {}))
        incar_dict.setdefault('SYSTEM', 'CUSTOM')

        lines = []
        for key, value in incar_dict.items():
            k = key.upper()
            if isinstance(value, bool):
                v = '.TRUE.' if value else '.FALSE.'
            elif isinstance(value, float) and abs(value) != 0 and abs(value) < 1e-3:
                v = f'{value:.2E}'
            else:
                v = str(value)
            lines.append(f"{k} = {v}\n")

        with open(work_dir / "INCAR", 'w') as f:
            f.writelines(lines)
        logger.info("自定义INCAR已写入 %s（%d 个参数）", work_dir / "INCAR", len(incar_dict))

    async def _analyze_custom_results(self, work_dir: Path, vasp_result: Dict[str, Any]) -> Dict[str, Any]:
        """从 OUTCAR 提取通用结果摘要"""
        import re
        result: Dict[str, Any] = {
            'success': False,
            'convergence': False,
            'final_energy': None,
            'fermi_energy': None,
            'n_ionic_steps': None,
            'n_electronic_steps': None,
            'output_files': [],
            'computation_time': vasp_result.get('computation_time'),
            'process_id': vasp_result.get('process_id'),
        }

        if not vasp_result.get('success', False):
            result['error'] = vasp_result.get('error_message', 'VASP计算失败')
            return result

        outcar_path = work_dir / "OUTCAR"
        if not outcar_path.exists():
            result['error'] = "未找到OUTCAR文件"
            return result

        try:
            with open(outcar_path, 'r', errors='ignore') as f:
                content = f.read()

            result['convergence'] = 'General timing and accounting informations for this job:' in content
            result['success'] = result['convergence']
            result['final_energy'] = self._extract_energy(outcar_path)

            fermi_match = re.search(r'E-fermi\s*:\s*([-\d.]+)', content)
            if fermi_match:
                result['fermi_energy'] = float(fermi_match.group(1))

            ionic_steps = len(re.findall(r'Iteration\s+\d+\(\s*1\)', content))
            result['n_ionic_steps'] = ionic_steps if ionic_steps > 0 else None

            elec_matches = re.findall(r'Iteration\s+\d+\(\s*(\d+)\)', content)
            if elec_matches:
                result['n_electronic_steps'] = int(elec_matches[-1])

        except Exception as e:
            logger.warning("分析自定义计算OUTCAR时出错: %s", e)

        for fname in ("OUTCAR", "CONTCAR", "POSCAR", "EIGENVAL", "DOSCAR",
                      "CHGCAR", "WAVECAR", "PROCAR", "LOCPOT", "vasprun.xml"):
            if (work_dir / fname).exists():
                result['output_files'].append(fname)

        return result
