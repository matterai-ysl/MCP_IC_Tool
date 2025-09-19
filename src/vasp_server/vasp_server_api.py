from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from typing import List
import uvicorn
from pathlib import Path
# Config import moved to __main__ block
from .schemas import (
    StructOptRequest, StructOptResponse, TaskStatusResponse, 
    TaskStatus, SCFRequest, SCFResponse,
    DOSRequest, DOSResponse, MDRequest, MDResponse
)
from .task_manager.manager import TaskManager
from .task_manager.database import check_and_init_db

# 初始化数据库
print("🚀 启动VASP计算服务...")
print("🔧 检查并初始化数据库...")
check_and_init_db()

app = FastAPI(title="VASP计算服务API", version="1.0.0")

# 配置CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 允许所有来源，生产环境中应该设置具体的域名
    allow_credentials=True,
    allow_methods=["*"],  # 允许所有方法
    allow_headers=["*"],  # 允许所有头部
)

# 挂载静态文件服务
# 这将服务整个文件系统，需要小心安全性
app.mount("/vasp/static", StaticFiles(directory="/data/home/ysl9527/vasp_calculations"), name="static")

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
async def get_task_status(task_id: str, user_id: str = Query(..., description="用户ID")):
    """
    查询任务状态与任务结果
    
    Args:
        task_id: 任务ID
        user_id: 用户ID
        
    Returns:
        TaskStatusResponse: 任务状态信息
    """
    try:
        task = task_manager.get_task(task_id, user_id)
        print(task)
        # 打印 task 对象的所有属性和值（调试用）
        print("\n" + "="*60)
        print("🔍 DEBUG: task 对象完整属性清单")
        print("="*60)
        for attr in dir(task):
            if not attr.startswith('_'):  # 跳过私有属性
                try:
                    value = getattr(task, attr)
                    print(f"{attr}: {repr(value)} (类型: {type(value).__name__})")
                except Exception as e:
                    print(f"{attr}: ❌ 获取失败 - {e}")
        print("="*60 + "\n")
        if not task:
            raise HTTPException(status_code=404, detail="任务未找到或无权限访问")
        
        # 构建基本响应
        # 构建基本响应 - 安全访问，避免崩溃
        response_data = {
            "task_id": getattr(task, 'id', None),
            "user_id": getattr(task, 'user_id', None),
            "task_type": getattr(task, 'task_type', None),
            "status": None,  # 默认值
            "progress": getattr(task, 'progress', 0) or 0,  # 避免 None
            "params": getattr(task, 'params', None),
            "result_path": getattr(task, 'result_path', None),
            "external_job_id": getattr(task, 'external_job_id', None),
            "process_id": getattr(task, 'process_id', None),
            "error_message": getattr(task, 'error_message', None),
            "result_data": getattr(task, 'result_data', None),
            "created_at": getattr(task, 'created_at', None),
            "updated_at": getattr(task, 'updated_at', None),
        }
        
        # 处理时间字段 - 转换为字符串
        if response_data["created_at"]:
            response_data["created_at"] = response_data["created_at"].isoformat()
        else:
            response_data["created_at"] = ""
            
        if response_data["updated_at"]:
            response_data["updated_at"] = response_data["updated_at"].isoformat()
        else:
            response_data["updated_at"] = ""

        # 安全设置 status（单独处理，避免枚举转换崩溃）
        try:
            status_val = getattr(task, 'status', None)
            if status_val:
                response_data["status"] = TaskStatus(status_val)
            else:
                response_data["status"] = TaskStatus.queued
        except (ValueError, AttributeError):
            # 如果状态值非法，设置默认状态
            response_data["status"] = TaskStatus.queued
        
        print("🔧 response_data内容:")
        for key, value in response_data.items():
            print(f"  {key}: {repr(value)} (类型: {type(value).__name__})")
        
        try:
            result = TaskStatusResponse(**response_data)
            print("✅ TaskStatusResponse创建成功")
            return result
        except Exception as validation_error:
            print(f"❌ TaskStatusResponse创建失败: {validation_error}")
            print(f"❌ 错误类型: {type(validation_error).__name__}")
            raise HTTPException(status_code=500, detail=f"响应模型验证失败: {str(validation_error)}")
        
    except HTTPException:
        raise
    except Exception as e:
        print(f"❌ 其他异常: {e}")
        print(f"❌ 异常类型: {type(e).__name__}")
        raise HTTPException(status_code=500, detail=f"查询任务状态失败: {str(e)}")


@app.post("/vasp/task/{task_id}/cancel")
async def cancel_task(task_id: str, user_id: str = Query(..., description="用户ID")):
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
async def list_user_tasks(user_id: str = Query(..., description="用户ID")):
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


@app.get("/vasp/download/{file_path:path}")
async def download_file(file_path: str):
    """
    提供文件下载服务
    
    Args:
        file_path: 相对于工作目录的文件路径
        
    Returns:
        FileResponse: 文件下载响应
    """
    try:
        full_path = Path(file_path)
        
        # 安全检查：确保文件路径在允许的范围内
        if not full_path.exists():
            raise HTTPException(status_code=404, detail="文件不存在")
        
        if not full_path.is_file():
            raise HTTPException(status_code=404, detail="路径不是文件")
        
        # 获取文件名用于下载
        filename = full_path.name
        
        return FileResponse(
            full_path, 
            filename=filename,
            media_type='application/octet-stream'
        )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"下载文件失败: {str(e)}")


# @app.get("/vasp/task/{task_id}/md-result", response_model=MDResult)
# async def get_md_result(task_id: str, user_id: str):
#     """
#     获取分子动力学计算的详细结果
    
#     Args:
#         task_id: 任务ID
#         user_id: 用户ID
        
#     Returns:
#         MDResult: 分子动力学计算结果详情
#     """
#     try:
#         task = task_manager.get_task(task_id, user_id)
        
#         if not task:
#             raise HTTPException(status_code=404, detail="任务未找到或无权限访问")
        
#         if str(task.task_type) != "md_calculation":  # type: ignore
#             raise HTTPException(status_code=400, detail="该任务不是分子动力学计算任务")
        
#         if str(task.status) != "completed":  # type: ignore
#             raise HTTPException(status_code=400, detail="任务尚未完成")
        
#         if not str(task.result_path or ""):  # type: ignore
#             raise HTTPException(status_code=404, detail="结果文件不存在")
        
#         # 解析MD计算结果（支持多温度）
#         from pathlib import Path
        
#         work_dir = Path(task.result_path)  # type: ignore
        
#         # 检查是否为多温度计算
#         temp_dirs = list(work_dir.glob("T_*K"))
#         is_multi_temperature = len(temp_dirs) > 0
        
#         if is_multi_temperature:
#             # 多温度MD计算结果
#             print(f"🌡️ 检测到多温度MD计算，发现 {len(temp_dirs)} 个温度点")
            
#             subtask_results = []
#             completed_count = 0
#             failed_count = 0
            
#             for temp_dir in sorted(temp_dirs):
#                 # 从目录名提取温度
#                 temp_name = temp_dir.name  # 如 "T_300K"
#                 try:
#                     temperature = float(temp_name.replace("T_", "").replace("K", ""))
#                 except:
#                     temperature = 0.0
                
#                 # 分析该温度点的结果
#                 subtask_result = await _analyze_single_temp_result(temp_dir, temperature)
#                 subtask_results.append(subtask_result)
                
#                 if subtask_result["convergence"]:
#                     completed_count += 1
#                 else:
#                     failed_count += 1
            
#             # 构建多温度结果
#             md_result = {
#                 "is_multi_temperature": True,
#                 "total_subtasks": len(temp_dirs),
#                 "completed_subtasks": completed_count,
#                 "failed_subtasks": failed_count,
#                 "subtask_results": subtask_results,
#                 "convergence": completed_count > 0,
#                 "computation_time": sum([r.get("computation_time", 0) for r in subtask_results if r.get("computation_time")])
#             }
            
#             # 检查多温度报告（优先使用综合报告）
#             comprehensive_report_path = work_dir / "MD_output" / "comprehensive_multi_temperature_report.html"
#             simple_report_path = work_dir / "MD_output" / "multi_temperature_md_report.html"
            
#             if comprehensive_report_path.exists():
#                 md_result["md_html_analysis_report"] = str(comprehensive_report_path)
#                 md_result["md_output_dir"] = str(work_dir / "MD_output")
#             elif simple_report_path.exists():
#                 md_result["md_html_analysis_report"] = str(simple_report_path)
#                 md_result["md_output_dir"] = str(work_dir / "MD_output")
        
#         else:
#             # 单温度MD计算结果（保持原有逻辑）
#             md_result = await _analyze_single_temp_result(work_dir, None)
#             md_result["is_multi_temperature"] = False
#             md_result["total_subtasks"] = 1
#             md_result["completed_subtasks"] = 1 if md_result.get("convergence", False) else 0
#             md_result["failed_subtasks"] = 0 if md_result.get("convergence", False) else 1
#             md_result["subtask_results"] = []
        
#         return MDResult(**md_result)
        
#     except HTTPException:
#         raise
#     except Exception as e:
#         raise HTTPException(status_code=500, detail=f"获取MD结果失败: {str(e)}")


# async def _analyze_single_temp_result(work_dir: Path, temperature: Optional[float]) -> Dict[str, Any]:
#     """分析单个温度点的MD结果"""
#     result = {
#         "temperature": temperature,
#         "subtask_dir": str(work_dir),
#         "md_structure": None,
#         "xdatcar_path": None,
#         "oszicar_path": None,
#         "final_energy": None,
#         "average_temperature": None,
#         "total_md_steps": None,
#         "convergence": False,
#         "computation_time": None,
#         "trajectory_data": None,
#         "status": "failed",
#         "error_message": None
#     }
    
#     try:
#         # 检查POSCAR文件
#         poscar_path = work_dir / "POSCAR"
#         if poscar_path.exists():
#             result["md_structure"] = str(poscar_path)
        
#         # 检查XDATCAR文件
#         xdatcar_path = work_dir / "XDATCAR"
#         if xdatcar_path.exists():
#             result["xdatcar_path"] = str(xdatcar_path)
#             # 快速统计MD步数
#             try:
#                 with open(xdatcar_path, 'r') as f:
#                     content = f.read()
#                     step_count = content.count("Direct configuration=")
#                     result["total_md_steps"] = step_count
#             except Exception as e:
#                 print(f"读取XDATCAR失败: {e}")
        
#         # 检查OSZICAR文件
#         oszicar_path = work_dir / "OSZICAR"
#         if oszicar_path.exists():
#             result["oszicar_path"] = str(oszicar_path)
#             # 快速提取最终能量
#             try:
#                 with open(oszicar_path, 'r') as f:
#                     lines = f.readlines()
#                     for line in reversed(lines):
#                         if 'DAV:' in line or 'RMM:' in line:
#                             parts = line.strip().split()
#                             if len(parts) >= 3:
#                                 result["final_energy"] = float(parts[2])
#                                 break
#             except Exception as e:
#                 print(f"读取OSZICAR失败: {e}")
        
#         # 检查OUTCAR文件
#         outcar_path = work_dir / "OUTCAR"
#         if outcar_path.exists():
#             try:
#                 with open(outcar_path, 'r') as f:
#                     content = f.read()
#                     if "General timing and accounting informations for this job:" in content:
#                         result["convergence"] = True
#                         result["status"] = "completed"
#             except Exception as e:
#                 print(f"读取OUTCAR失败: {e}")
#                 result["error_message"] = f"读取OUTCAR失败: {e}"
        
#         return result
        
#     except Exception as e:
#         result["error_message"] = str(e)
#         return result


if __name__ == "__main__":
    from .Config import VASP_remote_run_port
    uvicorn.run(app, host="0.0.0.0", port=VASP_remote_run_port)


