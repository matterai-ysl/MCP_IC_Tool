from fastapi import FastAPI, BackgroundTasks, HTTPException, Depends
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import HTTPBearer
from typing import List
import uvicorn

from .schemas import (
    StructOptRequest, StructOptResponse, TaskStatusResponse, 
    TaskStatus, StructOptResult, SCFRequest, SCFResponse, SCFResult,
    DOSRequest, DOSResponse, DOSResult, MDRequest, MDResponse, MDResult
)
from .task_manager.manager import TaskManager
from .task_manager.models import Task
from .task_manager.database import check_and_init_db

# 初始化数据库
print("🚀 启动VASP计算服务...")
print("🔧 检查并初始化数据库...")
check_and_init_db()

app = FastAPI(title="VASP计算服务API", version="1.0.0")
security = HTTPBearer()

# 配置CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 允许所有来源，生产环境中应该设置具体的域名
    allow_credentials=True,
    allow_methods=["*"],  # 允许所有方法
    allow_headers=["*"],  # 允许所有头部
)

# 创建全局任务管理器实例
task_manager = TaskManager()

@app.get("/")
async def root():
    return {"message": "VASP计算服务API", "version": "1.0.0"}


@app.post("/vasp/structure-optimization", response_model=StructOptResponse)
async def submit_structure_optimization(request: StructOptRequest):
    """
    提交结构优化任务
    
    支持两种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件
    2. CIF URL：直接从指定URL下载CIF文件
    
    Returns:
        StructOptResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        if not request.formula and not request.cif_url:
            raise HTTPException(status_code=400, detail="必须提供 formula 或 cif_url 中的一个")
        
        if request.formula and request.cif_url:
            raise HTTPException(status_code=400, detail="不能同时提供 formula 和 cif_url")
        
        # 准备任务参数
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "calc_type": request.calc_type.value,
            "kpoint_density": request.kpoint_density,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        if request.user_id is None:
            request.user_id = "123"
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="structure_optimization",
            params=task_params
        )
        
        return StructOptResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"结构优化任务已提交，任务ID: {task_id}"
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交任务失败: {str(e)}")


@app.post("/vasp/scf-calculation", response_model=SCFResponse)
async def submit_scf_calculation(request: SCFRequest):
    """
    提交自洽场计算任务
    
    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件
    2. CIF URL：直接从指定URL下载CIF文件
    3. 优化任务ID：基于已完成的结构优化任务的CONTCAR文件
    
    Returns:
        SCFResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url), 
            bool(request.optimized_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400, 
                detail="必须提供 formula、cif_url 或 optimized_task_id 中的一个"
            )
        
        # 如果基于优化任务，验证任务存在性
        if request.optimized_task_id:
            opt_task = task_manager.get_task(request.optimized_task_id, request.user_id)
            if not opt_task:
                raise HTTPException(
                    status_code=404, 
                    detail=f"结构优化任务 {request.optimized_task_id} 未找到或无权限访问"
                )
            if str(opt_task.status) != "completed":  # type: ignore
                raise HTTPException(
                    status_code=400, 
                    detail=f"结构优化任务 {request.optimized_task_id} 尚未完成"
                )
        
        # 准备任务参数
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "optimized_task_id": request.optimized_task_id,
            "calc_type": request.calc_type.value,
            "kpoint_density": request.kpoint_density,
            "precision": request.precision,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="scf_calculation",
            params=task_params
        )
        
        input_source = "化学式" if request.formula else "CIF URL" if request.cif_url else "结构优化任务"
        
        return SCFResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"自洽场计算任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交自洽场计算任务失败: {str(e)}")


@app.post("/vasp/dos-calculation", response_model=DOSResponse)
async def submit_dos_calculation(request: DOSRequest):
    """
    提交态密度计算任务
    
    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件（需要先完成自洽场计算）
    2. CIF URL：直接从指定URL下载CIF文件（需要先完成自洽场计算）
    3. 自洽场任务ID：基于已完成的自洽场计算任务结果
    
    Returns:
        DOSResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url), 
            bool(request.scf_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400, 
                detail="必须提供 formula、cif_url 或 scf_task_id 中的一个"
            )
        
        # 如果基于自洽场任务，验证任务存在性
        if request.scf_task_id:
            scf_task = task_manager.get_task(request.scf_task_id, request.user_id)
            if not scf_task:
                raise HTTPException(
                    status_code=404, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 未找到或无权限访问"
                )
            if str(scf_task.status) != "completed":  # type: ignore
                raise HTTPException(
                    status_code=400, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 尚未完成"
                )
        
        # 准备任务参数
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,
            "calc_type": request.calc_type.value,
            "kpoint_density": request.kpoint_density,
            "kpoint_multiplier": request.kpoint_multiplier,
            "precision": request.precision,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="dos_calculation",
            params=task_params
        )
        
        if request.scf_task_id:
            input_source = "自洽场计算任务"
            calc_mode = "态密度计算"
        elif request.formula:
            input_source = "化学式"
            calc_mode = "单点自洽+DOS计算"
        else:
            input_source = "CIF URL"
            calc_mode = "单点自洽+DOS计算"
        
        return DOSResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"{calc_mode}任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交态密度计算任务失败: {str(e)}")


@app.post("/vasp/md-calculation", response_model=MDResponse)
async def submit_md_calculation(request: MDRequest):
    """
    提交分子动力学计算任务
    
    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件（纯MD计算）
    2. CIF URL：直接从指定URL下载CIF文件（纯MD计算）
    3. 自洽场任务ID：基于已完成的自洽场计算任务结果
    
    Returns:
        MDResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url), 
            bool(request.scf_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400, 
                detail="必须提供 formula、cif_url 或 scf_task_id 中的一个"
            )
        
        # 如果基于自洽场任务，验证任务存在性
        if request.scf_task_id:
            scf_task = task_manager.get_task(request.scf_task_id, request.user_id)
            if not scf_task:
                raise HTTPException(
                    status_code=404, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 未找到或无权限访问"
                )
            if str(scf_task.status) != "completed":  # type: ignore
                raise HTTPException(
                    status_code=400, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 尚未完成"
                )
        
        # 准备任务参数
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,
            "calc_type": request.calc_type.value,
            "md_steps": request.md_steps,
            "temperature": request.temperature,
            "time_step": request.time_step,
            "ensemble": request.ensemble,
            "precision": request.precision,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="md_calculation",
            params=task_params
        )
        
        if request.scf_task_id:
            input_source = "自洽场计算任务"
            calc_mode = "分子动力学计算"
        elif request.formula:
            input_source = "化学式"
            calc_mode = "纯MD计算"
        else:
            input_source = "CIF URL"
            calc_mode = "纯MD计算"
        
        return MDResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"{calc_mode}任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交分子动力学计算任务失败: {str(e)}")

 
@app.get("/vasp/task/{task_id}", response_model=TaskStatusResponse)
async def get_task_status(task_id: str, user_id: str):
    """
    查询任务状态
    
    Args:
        task_id: 任务ID
        user_id: 用户ID
        
    Returns:
        TaskStatusResponse: 任务状态信息
    """
    try:
        task = task_manager.get_task(task_id, user_id)
        
        if not task:
            raise HTTPException(status_code=404, detail="任务未找到或无权限访问")
        
        return TaskStatusResponse(
            task_id=task.id,  # type: ignore
            user_id=task.user_id,  # type: ignore
            task_type=task.task_type,  # type: ignore
            status=TaskStatus(task.status),  # type: ignore
            progress=task.progress,  # type: ignore
            params=task.params,  # type: ignore
            result_path=task.result_path,  # type: ignore
            external_job_id=task.external_job_id,  # type: ignore
            process_id=task.process_id,  # type: ignore
            error_message=task.error_message,  # type: ignore
            result_data=task.result_data,  # type: ignore
            created_at=task.created_at.isoformat() if task.created_at else "",  # type: ignore
            updated_at=task.updated_at.isoformat() if task.updated_at else ""  # type: ignore
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"查询任务状态失败: {str(e)}")


@app.post("/vasp/task/{task_id}/cancel")
async def cancel_task(task_id: str, user_id: str):
    """
    取消任务
    
    Args:
        task_id: 任务ID
        user_id: 用户ID
        
    Returns:
        dict: 取消结果
    """
    try:
        success = task_manager.cancel_task(task_id, user_id)
        
        if not success:
            raise HTTPException(status_code=404, detail="任务未找到或无法取消")
        
        return {"message": f"任务 {task_id} 已请求取消"}
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"取消任务失败: {str(e)}")


@app.get("/vasp/tasks", response_model=List[TaskStatusResponse])
async def list_user_tasks(user_id: str):
    """
    列出用户的所有任务
    
    Args:
        user_id: 用户ID
        
    Returns:
        List[TaskStatusResponse]: 任务列表
    """
    try:
        tasks = task_manager.list_tasks(user_id)
        
        return [
            TaskStatusResponse(
                task_id=task.id,  # type: ignore
                user_id=task.user_id,  # type: ignore
                task_type=task.task_type,  # type: ignore
                status=TaskStatus(task.status),  # type: ignore
                progress=task.progress,  # type: ignore
                params=task.params,  # type: ignore
                result_path=task.result_path,  # type: ignore
                external_job_id=task.external_job_id,  # type: ignore
                process_id=task.process_id,  # type: ignore
                error_message=task.error_message,  # type: ignore
                result_data=task.result_data,  # type: ignore
                created_at=task.created_at.isoformat() if task.created_at else "",  # type: ignore
                updated_at=task.updated_at.isoformat() if task.updated_at else ""  # type: ignore
            )
            for task in tasks
        ]
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取任务列表失败: {str(e)}")


@app.get("/vasp/task/{task_id}/result")
async def get_task_result(task_id: str, user_id: str):
    """
    获取任务结果文件
    
    Args:
        task_id: 任务ID
        user_id: 用户ID
        
    Returns:
        FileResponse: 结果文件
    """
    try:
        task = task_manager.get_task(task_id, user_id)
        
        if not task:
            raise HTTPException(status_code=404, detail="任务未找到或无权限访问")
        
        if str(task.status) != "completed":  # type: ignore
            raise HTTPException(status_code=400, detail="任务尚未完成")
        
        if not str(task.result_path or ""):  # type: ignore
            raise HTTPException(status_code=404, detail="结果文件不存在")
        
        # 这里可以返回压缩包或者特定的结果文件
        # 暂时返回工作目录路径信息
        return {"result_path": task.result_path, "message": "结果文件路径"}
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取结果失败: {str(e)}")


@app.get("/vasp/task/{task_id}/md-result", response_model=MDResult)
async def get_md_result(task_id: str, user_id: str):
    """
    获取分子动力学计算的详细结果
    
    Args:
        task_id: 任务ID
        user_id: 用户ID
        
    Returns:
        MDResult: 分子动力学计算结果详情
    """
    try:
        task = task_manager.get_task(task_id, user_id)
        
        if not task:
            raise HTTPException(status_code=404, detail="任务未找到或无权限访问")
        
        if str(task.task_type) != "md_calculation":  # type: ignore
            raise HTTPException(status_code=400, detail="该任务不是分子动力学计算任务")
        
        if str(task.status) != "completed":  # type: ignore
            raise HTTPException(status_code=400, detail="任务尚未完成")
        
        if not str(task.result_path or ""):  # type: ignore
            raise HTTPException(status_code=404, detail="结果文件不存在")
        
        # 解析MD计算结果
        from pathlib import Path
        
        work_dir = Path(task.result_path)  # type: ignore
        
        # 读取结果文件
        md_result = {
            "md_structure": None,
            "xdatcar_path": None,
            "oszicar_path": None,
            "final_energy": None,
            "average_temperature": None,
            "total_md_steps": None,
            "convergence": False,
            "computation_time": None,
            "trajectory_data": None
        }
        
        # 检查POSCAR文件
        poscar_path = work_dir / "POSCAR"
        if poscar_path.exists():
            md_result["md_structure"] = str(poscar_path)
        
        # 检查XDATCAR文件
        xdatcar_path = work_dir / "XDATCAR"
        if xdatcar_path.exists():
            md_result["xdatcar_path"] = str(xdatcar_path)
            # 快速统计MD步数
            try:
                with open(xdatcar_path, 'r') as f:
                    content = f.read()
                    step_count = content.count("Direct configuration=")
                    md_result["total_md_steps"] = step_count
            except Exception as e:
                print(f"读取XDATCAR失败: {e}")
        
        # 检查OSZICAR文件
        oszicar_path = work_dir / "OSZICAR"
        if oszicar_path.exists():
            md_result["oszicar_path"] = str(oszicar_path)
            # 快速提取最终能量
            try:
                with open(oszicar_path, 'r') as f:
                    lines = f.readlines()
                    for line in reversed(lines):
                        if 'DAV:' in line or 'RMM:' in line:
                            parts = line.strip().split()
                            if len(parts) >= 3:
                                md_result["final_energy"] = float(parts[2])
                                break
            except Exception as e:
                print(f"读取OSZICAR失败: {e}")
        
        # 检查OUTCAR文件
        outcar_path = work_dir / "OUTCAR"
        if outcar_path.exists():
            try:
                with open(outcar_path, 'r') as f:
                    content = f.read()
                    if "General timing and accounting informations for this job:" in content:
                        md_result["convergence"] = True
            except Exception as e:
                print(f"读取OUTCAR失败: {e}")
        
        return MDResult(**md_result)
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取MD结果失败: {str(e)}")


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)


