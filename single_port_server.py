#!/usr/bin/env python3
"""
单端口MCP服务器 - streamable-http + 文件服务
    
这个服务器使用FastMCP的streamable-http transport，并在同一个端口上集成文件服务功能：
1. MCP协议（streamable-http）
2. 文件下载服务
3. 静态文件服务
4. API端点

所有服务都运行在单个端口（默认8100）上。
"""

import logging
from pathlib import Path
import sys

# FastAPI相关
from starlette.routing import Route, Mount
from starlette.responses import JSONResponse, FileResponse
from starlette.requests import Request
from starlette.staticfiles import StaticFiles
from starlette.exceptions import HTTPException
import uvicorn

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "src"))

# 导入MCP服务器和FastMCP内部模块
from src.mcp_ic_tool.mcp_server import mcp
from fastmcp.server.http import create_streamable_http_app
from src.mcp_ic_tool.config import MCP_PORT
# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SinglePortMCPServer:
    """单端口MCP服务器类，使用streamable-http + 文件服务"""
    
    def __init__(self, port: int = MCP_PORT, host: str = "0.0.0.0"):
        self.port = port
        self.host = host
        self.trained_models_dir = Path("./trained_models")
        self.reports_dir = Path("./reports")
        
        # 确保目录存在
        self.trained_models_dir.mkdir(exist_ok=True)
        self.reports_dir.mkdir(exist_ok=True)
    
    def _create_file_service_routes(self):
        """创建文件服务路由"""
        routes = []
        
        async def get_server_info(_: Request):
            """获取服务器信息"""
            return JSONResponse({
                "service": "Single Port MCP BO Tool Server",
                "version": "3.0.0",
                "transport": "streamable-http",
                "port": self.port,
                "endpoints": {
                    "mcp_endpoint": f"http://{self.host}:{self.port}/mcp",  # 更新MCP端点路径
                    "server_info": f"http://{self.host}:{self.port}/api/info",
                    "health_check": f"http://{self.host}:{self.port}/api/health",
                    "model_list": f"http://{self.host}:{self.port}/api/models/list",
                    "file_download": f"http://{self.host}:{self.port}/api/download/file/{{file_path}}",
                    "static_models": f"http://{self.host}:{self.port}/static/models/{{path}}",
                    "static_reports": f"http://{self.host}:{self.port}/static/reports/{{path}}",
                },
                "directories": {
                    "trained_models": str(self.trained_models_dir.absolute()),
                    "reports": str(self.reports_dir.absolute())
                }
            })
        
        async def health_check(_: Request):
            """健康检查端点"""
            return JSONResponse({
                "status": "healthy",
                "mcp_status": "active",
                "transport": "streamable-http",
                "directories": {
                    "trained_models_exists": self.trained_models_dir.exists(),
                    "reports_exists": self.reports_dir.exists()
                },
                "port": self.port,
                "host": self.host
            })
        
        async def list_model_files(_: Request):
            """列出所有训练模型文件"""
            if not self.trained_models_dir.exists():
                return JSONResponse({"files": [], "message": "训练模型目录不存在"})
            
            files = []
            for item in self.trained_models_dir.rglob("*"):
                if item.is_file():
                    relative_path = item.relative_to(self.trained_models_dir)
                    files.append({
                        "path": str(relative_path),
                        "name": item.name,
                        "size": item.stat().st_size,
                        "download_url": f"http://{self.host}:{self.port}/api/download/file/{relative_path}",
                        "static_url": f"http://{self.host}:{self.port}/static/models/{relative_path}"
                    })
            
            return JSONResponse({"files": files, "count": len(files)})
        
        async def download_model_file(request: Request):
            """下载训练模型文件"""
            file_path = request.path_params["file_path"]
            
            # 构建完整文件路径
            full_path = self.trained_models_dir / file_path
            print(full_path)
            # 安全检查：确保路径在允许的目录内
            try:
                resolved_path = full_path.resolve()
                models_resolved = self.trained_models_dir.resolve()
                
                if not str(resolved_path).startswith(str(models_resolved)):
                    raise HTTPException(status_code=403, detail="访问被拒绝：路径超出允许范围")
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"无效的文件路径: {str(e)}")
            
            # 检查文件是否存在
            if not resolved_path.exists():
                raise HTTPException(status_code=404, detail=f"文件未找到: {file_path}")
            
            # 检查是否为文件
            if not resolved_path.is_file():
                raise HTTPException(status_code=400, detail=f"路径不是文件: {file_path}")
            
            # 返回文件
            return FileResponse(
                path=str(resolved_path),
                filename=resolved_path.name,
                media_type='application/octet-stream'
            )
        
        async def welcome_page(_: Request):
            """欢迎页面"""
            return JSONResponse({
                "message": "欢迎使用 Neural Network MCP 服务器",
                "service": "Single Port MCP BO Tool Server",
                "version": "3.0.0", 
                "mcp_endpoint": f"http://{self.host}:{self.port}/mcp",
                "server_info": f"http://{self.host}:{self.port}/api/info",
                "health_check": f"http://{self.host}:{self.port}/api/health"
            })
        
        # 添加API路由
        routes.extend([
            Route("/", endpoint=welcome_page, methods=["GET"]),  # 根路径欢迎页面
            Route("/api/info", endpoint=get_server_info, methods=["GET"]),
            Route("/api/health", endpoint=health_check, methods=["GET"]), 
            Route("/api/models/list", endpoint=list_model_files, methods=["GET"]),
            Route("/api/download/file/{file_path:path}", endpoint=download_model_file, methods=["GET"]),
            # 添加不带/api/前缀的下载路由，方便使用
            Route("/download/file/{file_path:path}", endpoint=download_model_file, methods=["GET"]),
        ])
        
        # 添加静态文件路由
        if self.trained_models_dir.exists():
            routes.append(
                Mount("/static", StaticFiles(directory=str(self.trained_models_dir)), name="static_models")
            )
            # 添加单数形式的路径以支持 /static/model/ 
            routes.append(
                Mount("/static/model", StaticFiles(directory=str(self.trained_models_dir)), name="static_model")
            )
        
        if self.reports_dir.exists():
            routes.append(
                Mount("/static/reports", StaticFiles(directory=str(self.reports_dir)), name="static_reports")
            )
        
        return routes
    
    def create_app(self):
        """创建集成了文件服务的streamable-http应用"""
        # 创建文件服务路由
        file_routes = self._create_file_service_routes()
        
        # 使用FastMCP的create_streamable_http_app函数，传入自定义路由
        app = create_streamable_http_app(
            server=mcp,  # 直接传入FastMCP对象
            streamable_http_path="/mcp",  # 修改：MCP协议使用专用路径，避免拦截所有请求
            routes=file_routes,  # 添加文件服务路由
            debug=False,
            stateless_http=True  # 关键：启用无状态HTTP模式，每个请求独立处理
        )
        
        return app
    
    async def start(self):
        """启动单端口服务器（异步版本）"""
        logger.info("启动单端口MCP服务器（streamable-http + 文件服务）...")
        logger.info(f"服务地址: http://{self.host}:{self.port}")
        logger.info(f"MCP端点: http://{self.host}:{self.port}/mcp")
        logger.info(f"服务器信息: http://{self.host}:{self.port}/api/info")
        logger.info(f"健康检查: http://{self.host}:{self.port}/api/health") 
        logger.info(f"模型文件列表: http://{self.host}:{self.port}/api/models/list")
        logger.info(f"训练模型目录: {self.trained_models_dir.absolute()}")
        logger.info(f"报告目录: {self.reports_dir.absolute()}")
        logger.info("="*60)
        
        
        # 创建应用
        app = self.create_app()
        
        # 创建uvicorn配置
        config = uvicorn.Config(
            app=app,
            host=self.host,
            port=self.port,
            log_level="info",
            access_log=True
        )
        
        # 创建服务器实例并异步运行
        server = uvicorn.Server(config)
        await server.serve()
    
    def start_sync(self):
        """同步启动方法（向后兼容）"""
        import asyncio
        asyncio.run(self.start())


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="单端口MCP服务器 - streamable-http + 文件服务",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  %(prog)s                           # 默认配置（localhost:8090）
  %(prog)s --host 0.0.0.0            # 绑定所有接口
  %(prog)s --port 9000               # 自定义端口
  %(prog)s --host 0.0.0.0 --port 9000 # 自定义主机和端口

服务端点（单端口）:
  MCP端点:               http://host:port/mcp
  欢迎页面:              http://host:port/
  服务器信息:            http://host:port/api/info
  健康检查:              http://host:port/api/health
  模型文件列表:          http://host:port/api/models/list
  文件下载:              http://host:port/api/download/file/{path}
  静态模型文件:          http://host:port/static/models/{path}
  静态报告文件:          http://host:port/static/reports/{path}
        """
    )
    
    parser.add_argument(
        "--host",
        default="localhost",
        help="绑定主机地址（默认: localhost）"
    )
    
    parser.add_argument(
        "--port",
        type=int,
        default=MCP_PORT,
        help="端口号（默认: 8130）"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="启用调试模式"
    )
    
    args = parser.parse_args()
    
    # 设置日志级别
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info("调试模式已启用")
    
    try:
        # 创建并启动服务器
        server = SinglePortMCPServer(port=args.port, host=args.host)
        server.start_sync()  # 使用同步启动方法
        
    except KeyboardInterrupt:
        logger.info("\n服务器已停止")
    except Exception as e:
        logger.error(f"服务器启动失败: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()