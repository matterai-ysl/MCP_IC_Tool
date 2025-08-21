#!/usr/bin/env python3
"""
自洽场计算使用示例

演示如何使用VASP API进行自洽场计算的三种方式：
1. 从化学式进行自洽场计算
2. 从CIF URL进行自洽场计算  
3. 基于已完成的结构优化结果进行自洽场计算
"""

import requests
import json
import time

def submit_scf_from_formula():
    """从化学式提交自洽场计算"""
    print("🔬 示例1: 从化学式进行自洽场计算")
    
    # API请求数据
    data = {
        "user_id": "example_user",
        "formula": "Li2O",
        "calc_type": "SSE",  # 自旋轨道耦合
        "precision": "Accurate",
        "kpoint_density": 30.0,
        "stable_only": True,
        "selection_mode": "most_stable"
    }
    
    try:
        response = requests.post(
            "http://localhost:9000/vasp/scf-calculation",
            json=data,
            timeout=30
        )
        response.raise_for_status()
        result = response.json()
        
        print(f"✅ 任务提交成功:")
        print(f"   任务ID: {result['task_id']}")
        print(f"   状态: {result['status']}")
        print(f"   消息: {result['message']}")
        
        return result['task_id']
        
    except Exception as e:
        print(f"❌ 提交失败: {e}")
        return None

def submit_scf_from_optimization(opt_task_id):
    """基于结构优化结果进行自洽场计算"""
    print(f"🔗 示例2: 基于优化任务 {opt_task_id[:8]}... 进行自洽场计算")
    
    data = {
        "user_id": "example_user",
        "optimized_task_id": opt_task_id,
        "calc_type": "OXC",  # 氧化计算
        "precision": "Accurate",
        "kpoint_density": 35.0
    }
    
    try:
        response = requests.post(
            "http://localhost:9000/vasp/scf-calculation",
            json=data,
            timeout=30
        )
        response.raise_for_status()
        result = response.json()
        
        print(f"✅ 自洽场计算任务提交成功:")
        print(f"   任务ID: {result['task_id']}")
        print(f"   状态: {result['status']}")
        
        return result['task_id']
        
    except Exception as e:
        print(f"❌ 提交失败: {e}")
        return None

def monitor_task(task_id, user_id="example_user"):
    """监控任务进度"""
    print(f"\n📊 监控任务 {task_id[:8]}...")
    
    for i in range(20):  # 监控5分钟
        try:
            response = requests.get(
                f"http://localhost:9000/vasp/task/{task_id}",
                params={"user_id": user_id},
                timeout=10
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"  状态: {result['status']} | 进度: {result['progress']}%", end="")
            if result.get('process_id'):
                print(f" | PID: {result['process_id']}", end="")
            print()
            
            if result['status'] in ['completed', 'failed', 'canceled']:
                if result['status'] == 'completed':
                    print("🎉 计算完成!")
                    print(f"   结果路径: {result.get('result_path', 'N/A')}")
                else:
                    print(f"❌ 任务状态: {result['status']}")
                    if result.get('error_message'):
                        print(f"   错误信息: {result['error_message']}")
                break
            
            time.sleep(15)  # 每15秒检查一次
            
        except Exception as e:
            print(f"❌ 监控出错: {e}")
            break
    else:
        print("⏰ 监控超时，任务可能仍在进行中")

def main():
    """主函数"""
    print("="*60)
    print("🧪 VASP自洽场计算示例")
    print("="*60)
    
    # 测试连接
    try:
        response = requests.get("http://localhost:9000/", timeout=5)
        print("✅ API服务连接正常")
    except Exception as e:
        print(f"❌ 无法连接到API服务: {e}")
        print("请确保:")
        print("1. VASP服务正在运行")
        print("2. 端口转发已设置: ssh -L 9000:localhost:8000 username@supercomputer")
        return
    
    # 示例1: 从化学式进行自洽场计算
    task_id = submit_scf_from_formula()
    if task_id:
        monitor_task(task_id)
    
    print("\n" + "-"*60)
    
    # 示例2: 查找已完成的优化任务，基于其结果进行自洽场计算
    print("🔍 查找已完成的结构优化任务...")
    try:
        response = requests.get(
            "http://localhost:9000/vasp/tasks",
            params={"user_id": "example_user"},
            timeout=10
        )
        response.raise_for_status()
        tasks = response.json()
        
        completed_opt_tasks = [
            task for task in tasks 
            if task['task_type'] == 'structure_optimization' and task['status'] == 'completed'
        ]
        
        if completed_opt_tasks:
            opt_task = completed_opt_tasks[0]  # 选择第一个
            print(f"📋 找到优化任务: {opt_task['task_id'][:8]}...")
            
            scf_task_id = submit_scf_from_optimization(opt_task['task_id'])
            if scf_task_id:
                monitor_task(scf_task_id)
        else:
            print("⚠️ 没有找到已完成的结构优化任务")
            print("建议先运行结构优化任务")
            
    except Exception as e:
        print(f"❌ 查询任务失败: {e}")
    
    print("\n🎯 示例完成!")
    print("提示:")
    print("- 自洽场计算通常在结构优化后进行")
    print("- 可以用于精确计算电子结构、能带、态密度等")
    print("- 支持不同精度设置 (Normal, High, Accurate)")

if __name__ == "__main__":
    main() 