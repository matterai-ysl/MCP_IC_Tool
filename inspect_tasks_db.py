import sqlite3

# 1. 设置数据库路径（请修改为你的实际路径）
db_path = "/Users/ysl/Desktop/Code/MCP_IC_Tool/tasks.db"  # 例如："/Users/ysl/project/tasks.db"

# 2. 连接数据库
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# 3. 查询 id 和 task_type
cursor.execute("SELECT id, task_type FROM tasks;")
rows = cursor.fetchall()

# 4. 打印结果
print(f"{'ID':<8} {'任务类型'}")
print("-" * 30)
for row in rows:
    task_id = row[0]
    task_type = row[1] if row[1] else "NULL"  # 防止为空
    print(f"{task_id:<8} {task_type}")

# 5. 关闭连接
conn.close()