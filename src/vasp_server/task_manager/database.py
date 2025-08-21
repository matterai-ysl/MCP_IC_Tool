import os
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker, declarative_base

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(BASE_DIR, "tasks.db")
SQLALCHEMY_DATABASE_URL = f"sqlite:///{DB_PATH}"

engine = create_engine(
    SQLALCHEMY_DATABASE_URL,
    connect_args={"check_same_thread": False},
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()


def init_db() -> None:
    """初始化数据库表"""
    # 延迟导入以避免循环依赖
    from . import models  # noqa: F401

    Base.metadata.create_all(bind=engine)
    print("✅ 数据库表已初始化")


def check_and_init_db() -> None:
    """检查数据库表是否存在，如果不存在则自动创建"""
    try:
        # 检查数据库文件是否存在
        if not os.path.exists(DB_PATH):
            print(f"📁 数据库文件不存在，将创建: {DB_PATH}")
            init_db()
            return

        # 检查表是否存在
        inspector = inspect(engine)
        tables = inspector.get_table_names()
        
        if 'tasks' not in tables:
            print("⚠️  数据库表不存在，开始自动初始化...")
            init_db()
        else:
            print("✅ 数据库表已存在")
            
    except Exception as e:
        print(f"⚠️  数据库检查失败，尝试初始化: {e}")
        try:
            init_db()
        except Exception as init_error:
            print(f"❌ 数据库初始化失败: {init_error}")
            raise


# 注意：数据库初始化现在在API服务启动时手动调用 