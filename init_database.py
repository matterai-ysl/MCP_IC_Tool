#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据库初始化脚本

用于在超算上初始化VASP API所需的SQLite数据库
"""

import os
import sys
from pathlib import Path

# 添加项目路径到sys.path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

try:
    from src.vasp_server.task_manager.database import engine, Base, init_db
    from src.vasp_server.task_manager.models import Task
    print("✅ 成功导入数据库模块")
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保在项目根目录运行此脚本")
    sys.exit(1)


def create_database():
    """创建数据库表"""
    print("🔧 开始创建数据库表...")
    
    try:
        # 创建所有表
        Base.metadata.create_all(bind=engine)
        print("✅ 数据库表创建成功")
        
        # 验证表是否创建成功
        from sqlalchemy import inspect
        inspector = inspect(engine)
        tables = inspector.get_table_names()
        
        print(f"📋 已创建的表: {tables}")
        
        if 'tasks' in tables:
            print("✅ tasks 表创建成功")
            
            # 显示表结构
            columns = inspector.get_columns('tasks')
            print("📊 tasks 表结构:")
            for col in columns:
                print(f"   - {col['name']}: {col['type']}")
        else:
            print("❌ tasks 表创建失败")
            return False
            
        return True
        
    except Exception as e:
        print(f"❌ 创建数据库失败: {e}")
        return False


def check_database():
    """检查数据库状态"""
    print("🔍 检查数据库状态...")
    
    try:
        from sqlalchemy import text
        from src.vasp_server.task_manager.database import SessionLocal
        
        db = SessionLocal()
        
        # 检查表是否存在
        result = db.execute(text("SELECT name FROM sqlite_master WHERE type='table';"))
        tables = [row[0] for row in result]
        
        print(f"📋 数据库中的表: {tables}")
        
        if 'tasks' in tables:
            # 检查tasks表的记录数
            result = db.execute(text("SELECT COUNT(*) FROM tasks;"))
            row = result.fetchone()
            count = row[0] if row else 0
            print(f"📊 tasks 表中有 {count} 条记录")
        
        db.close()
        return True
        
    except Exception as e:
        print(f"❌ 检查数据库失败: {e}")
        return False


def reset_database():
    """重置数据库（删除所有表并重新创建）"""
    print("🗑️  重置数据库...")
    
    try:
        # 删除所有表
        Base.metadata.drop_all(bind=engine)
        print("✅ 已删除所有表")
        
        # 重新创建表
        return create_database()
        
    except Exception as e:
        print(f"❌ 重置数据库失败: {e}")
        return False


def show_database_info():
    """显示数据库信息"""
    try:
        # 从engine获取数据库URL
        db_url = str(engine.url)
        print(f"📍 数据库URL: {db_url}")
        
        # 如果是SQLite，显示文件路径
        if db_url.startswith("sqlite:///"):
            db_path = db_url.replace("sqlite:///", "")
            if os.path.exists(db_path):
                size = os.path.getsize(db_path)
                print(f"📁 数据库文件: {db_path}")
                print(f"📊 文件大小: {size} 字节")
            else:
                print(f"⚠️  数据库文件不存在: {db_path}")
        
    except Exception as e:
        print(f"❌ 获取数据库信息失败: {e}")


def main():
    """主函数"""
    print("🚀 VASP API 数据库初始化工具")
    print("=" * 50)
    
    show_database_info()
    
    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
        
        if command == "create":
            success = create_database()
        elif command == "check":
            success = check_database()
        elif command == "reset":
            print("⚠️  警告: 这将删除所有现有数据!")
            confirm = input("确认重置数据库? (yes/NO): ").strip().lower()
            if confirm == "yes":
                success = reset_database()
            else:
                print("❌ 已取消重置操作")
                return
        else:
            print(f"❌ 未知命令: {command}")
            print("可用命令: create, check, reset")
            return
    else:
        # 默认：检查数据库，如果表不存在则创建
        print("🔍 检查数据库状态...")
        
        try:
            from sqlalchemy import inspect
            inspector = inspect(engine)
            tables = inspector.get_table_names()
            
            if 'tasks' in tables:
                print("✅ 数据库已存在且正常")
                check_database()
            else:
                print("⚠️  数据库表不存在，开始创建...")
                success = create_database()
                if success:
                    print("🎉 数据库初始化完成！")
                else:
                    print("❌ 数据库初始化失败")
                    sys.exit(1)
                    
        except Exception as e:
            print(f"❌ 数据库操作失败: {e}")
            print("\n尝试创建数据库...")
            success = create_database()
            if not success:
                sys.exit(1)


if __name__ == "__main__":
    main()


# =============================================================================
# 使用说明
# =============================================================================
"""
使用方法:

1. 初始化数据库（默认）:
   python init_database.py

2. 仅创建表:
   python init_database.py create

3. 检查数据库状态:
   python init_database.py check

4. 重置数据库（危险操作）:
   python init_database.py reset

在超算上的使用步骤:
1. 上传此脚本到项目根目录
2. 运行: python init_database.py
3. 启动API服务: python -m src.vasp_server.vasp_server_api
""" 