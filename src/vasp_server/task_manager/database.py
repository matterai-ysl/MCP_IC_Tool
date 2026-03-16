import os
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker, declarative_base

# 云端数据库配置
DB_HOST = os.getenv("DB_HOST", "pgm-uf69uij17z9vh123jo.pg.rds.aliyuncs.com")
DB_PORT = os.getenv("DB_PORT", "5432")
DB_NAME = os.getenv("DB_NAME", "ADK")
DB_USER = os.getenv("DB_USER", "a2252222223")
DB_PASSWORD = os.getenv("DB_PASSWORD", "Jixiaobei123")

SQLALCHEMY_DATABASE_URL = f"postgresql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

engine = create_engine(
    SQLALCHEMY_DATABASE_URL,
    pool_pre_ping=True,
    pool_recycle=300,
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
        print(f"🔗 连接云端数据库: {DB_HOST}:{DB_PORT}/{DB_NAME}")

        # 检查表是否存在
        inspector = inspect(engine)
        tables = inspector.get_table_names()

        required_tables = {'tasks', 'execution_attempts', 'analysis_runs', 'artifacts'}
        missing = required_tables - set(tables)
        if missing:
            print(f"⚠️  缺少数据库表 {missing}，开始自动初始化...")
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