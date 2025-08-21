#!/usr/bin/env python3
"""
态密度(DOS)计算使用示例

演示如何使用VASP API进行态密度计算，包括：
1. 基于已完成的自洽场计算进行态密度计算
2. 完整的工作流程：结构优化 → 自洽场计算 → 态密度计算
3. 态密度结果分析
"""

import requests
import json
import time

def submit_dos_from_scf(scf_task_id: str):
    """基于自洽场计算结果提交态密度计算"""
    print(f"📊 示例: 基于自洽场任务 {scf_task_id[:8]}... 进行态密度计算")
    
    data = {
        "user_id": "dos_example_user",
        "scf_task_id": scf_task_id,
        "calc_type": "OXC",  # 氧化计算
        "kpoint_multiplier": 2.0,  # K点倍增因子
        "precision": "Accurate"
    }
    
    try:
        response = requests.post(
            "http://localhost:9000/vasp/dos-calculation",
            json=data,
            timeout=30
        )
        response.raise_for_status()
        result = response.json()
        
        print(f"✅ 态密度计算任务提交成功:")
        print(f"   任务ID: {result['task_id']}")
        print(f"   状态: {result['status']}")
        print(f"   消息: {result['message']}")
        
        return result['task_id']
        
    except Exception as e:
        print(f"❌ 提交失败: {e}")
        return None

def find_completed_scf_tasks(user_id: str = "dos_example_user"):
    """查找已完成的自洽场计算任务"""
    print("🔍 查找已完成的自洽场计算任务...")
    
    try:
        response = requests.get(
            "http://localhost:9000/vasp/tasks",
            params={"user_id": user_id},
            timeout=10
        )
        response.raise_for_status()
        tasks = response.json()
        
        completed_scf_tasks = [
            task for task in tasks 
            if task['task_type'] == 'scf_calculation' and task['status'] == 'completed'
        ]
        
        if completed_scf_tasks:
            print(f"✅ 找到 {len(completed_scf_tasks)} 个已完成的自洽场任务:")
            for i, task in enumerate(completed_scf_tasks, 1):
                formula = task.get('params', {}).get('formula', 'unknown')
                print(f"   {i}. {task['task_id'][:8]}... ({formula})")
            
            return completed_scf_tasks
        else:
            print("⚠️ 没有找到已完成的自洽场计算任务")
            return []
            
    except Exception as e:
        print(f"❌ 查询任务失败: {e}")
        return []

def monitor_dos_task(task_id: str, user_id: str = "dos_example_user"):
    """监控态密度计算任务"""
    print(f"\n📊 监控态密度计算任务 {task_id[:8]}...")
    
    for i in range(30):  # 监控10分钟
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
                    print("🎉 态密度计算完成!")
                    print(f"   结果路径: {result.get('result_path', 'N/A')}")
                    
                    # 显示DOS计算的特殊结果
                    print("📈 态密度计算结果:")
                    print("   - DOSCAR文件已生成")
                    print("   - 可进行后续的能带分析、电子结构分析")
                    print("   - 支持态密度可视化和轨道分辨态密度分析")
                else:
                    print(f"❌ 任务状态: {result['status']}")
                    if result.get('error_message'):
                        print(f"   错误信息: {result['error_message']}")
                break
            
            time.sleep(20)  # 每20秒检查一次
            
        except Exception as e:
            print(f"❌ 监控出错: {e}")
            break
    else:
        print("⏰ 监控超时，任务可能仍在进行中")

def demonstrate_full_workflow():
    """演示完整的工作流程"""
    print("="*70)
    print("🧪 VASP 完整工作流程演示: 结构优化 → 自洽场计算 → 态密度计算")
    print("="*70)
    
    user_id = "workflow_demo_user"
    
    print("\n第一阶段: 结构优化")
    print("-" * 30)
    
    # 这里假设结构优化已经完成，实际使用中需要等待
    print("💡 通常需要先完成结构优化，获得优化后的结构")
    print("   示例: Li2O 结构优化 → 获得稳定的几何结构")
    
    print("\n第二阶段: 自洽场计算")
    print("-" * 30)
    
    print("💡 基于优化后的结构进行自洽场计算")
    print("   目的: 获得收敛的电子密度和波函数")
    print("   生成: CHG, CHGCAR, WAVECAR 等文件")
    
    print("\n第三阶段: 态密度计算")
    print("-" * 30)
    
    # 查找已完成的自洽场任务
    scf_tasks = find_completed_scf_tasks(user_id)
    
    if scf_tasks:
        # 选择第一个自洽场任务进行DOS计算
        selected_task = scf_tasks[0]
        print(f"📋 选择自洽场任务: {selected_task['task_id'][:8]}...")
        
        dos_task_id = submit_dos_from_scf(selected_task['task_id'])
        if dos_task_id:
            monitor_dos_task(dos_task_id, user_id)
    else:
        print("⚠️ 需要先完成自洽场计算才能进行态密度计算")
        print("\n建议步骤:")
        print("1. 提交结构优化任务并等待完成")
        print("2. 基于优化结果提交自洽场计算并等待完成")
        print("3. 基于自洽场结果提交态密度计算")

def main():
    """主函数"""
    print("="*60)
    print("🧪 VASP态密度(DOS)计算示例")
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
    
    # 演示完整工作流程
    demonstrate_full_workflow()
    
    print("\n🎯 示例完成!")
    print("\n📚 关于态密度计算:")
    print("- 态密度(DOS)描述了材料中电子态在能量空间的分布")
    print("- 可以分析材料的电子结构、导电性、磁性等性质")
    print("- LORBIT=11 可以计算轨道分辨态密度(PDOS)")
    print("- K点密度通常是结构优化的2-3倍以获得更精确的结果")
    print("- DOSCAR文件包含完整的态密度数据，可用于可视化分析")

if __name__ == "__main__":
    main() 