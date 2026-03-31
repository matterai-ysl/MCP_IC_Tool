import logging
from sqlalchemy import create_engine, event, inspect
from sqlalchemy.orm import declarative_base, scoped_session, sessionmaker

from ..settings import settings

logger = logging.getLogger(__name__)

SQLALCHEMY_DATABASE_URL = settings.effective_database_url

_engine_kwargs: dict = dict(pool_pre_ping=True, pool_recycle=300)

if SQLALCHEMY_DATABASE_URL.startswith("sqlite"):
    _engine_kwargs["connect_args"] = {"check_same_thread": False}

engine = create_engine(SQLALCHEMY_DATABASE_URL, **_engine_kwargs)

# Enable WAL mode for SQLite to allow concurrent readers/writers.
if SQLALCHEMY_DATABASE_URL.startswith("sqlite"):
    @event.listens_for(engine, "connect")
    def _set_sqlite_wal(dbapi_connection, connection_record):
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA journal_mode=WAL")
        cursor.close()

# Use scoped_session for thread-safe session management.
SessionLocal = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))
Base = declarative_base()


def database_dialect_name() -> str:
    return engine.dialect.name


def init_db() -> None:
    """初始化数据库表"""
    # 延迟导入以避免循环依赖
    from . import models  # noqa: F401

    Base.metadata.create_all(bind=engine)
    logger.info("数据库表已初始化")


def check_and_init_db() -> None:
    """检查数据库表是否存在，如果不存在则自动创建"""
    try:
        logger.info("Checking database schema...")

        # 检查表是否存在
        inspector = inspect(engine)
        tables = inspector.get_table_names()

        required_tables = {'tasks', 'execution_attempts', 'analysis_runs', 'artifacts'}
        missing = required_tables - set(tables)
        if missing:
            logger.warning("缺少数据库表 %s，开始自动初始化...", missing)
            init_db()
        else:
            logger.info("数据库表已存在")

    except Exception as e:
        logger.warning("数据库检查失败，尝试初始化: %s", e)
        try:
            init_db()
        except Exception as init_error:
            logger.error("数据库初始化失败: %s", init_error)
            raise


# 注意：数据库初始化现在在API服务启动时手动调用 
