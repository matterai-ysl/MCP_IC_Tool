#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VASP API 调试程序

用于测试部署在超算上的VASP结构优化API服务
端口转发地址: http://localhost:9000
"""

import requests
import json
import time
import sys
from typing import Dict, Any, Optional


class VASPAPIDebugger:
    """VASP API调试器"""
    
    def __init__(self, base_url: str = "http://localhost:9000"):
        self.base_url = base_url
        self.session = requests.Session()
        self.timeout = 30  # 30秒超时
        
    def test_connection(self) -> bool:
        """测试API连接"""
        print("🔗 测试API连接...")
        try:
            response = self.session.get(f"{self.base_url}/", timeout=self.timeout)
            response.raise_for_status()
            result = response.json()
            print(f"✅ 连接成功: {result}")
            return True
        except requests.exceptions.ConnectionError:
            print("❌ 连接失败: 无法连接到API服务")
            print("请检查:")
            print("1. 超算上的服务是否正在运行")
            print("2. 端口转发是否正确设置")
            print("3. 防火墙设置")
            return False
        except Exception as e:
            print(f"❌ 连接测试失败: {e}")
            return False
    
    def submit_structure_optimization(self, user_id: str, **kwargs) -> Optional[str]:
        """提交结构优化任务"""
        print(f"\n📤 提交结构优化任务 (用户: {user_id})...")
        
        # 构建请求数据
        data = {
            "user_id": user_id,
            **kwargs
        }
        
        print(f"请求数据: {json.dumps(data, indent=2, ensure_ascii=False)}")
        
        try:
            response = self.session.post(
                f"{self.base_url}/vasp/structure-optimization",
                json=data,
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 任务提交成功:")
            print(f"   任务ID: {result['task_id']}")
            print(f"   状态: {result['status']}")
            print(f"   消息: {result['message']}")
            
            return result['task_id']
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 任务提交失败: {e}")
            return None
    
    def submit_scf_calculation(self, user_id: str, **kwargs) -> Optional[str]:
        """提交自洽场计算任务"""
        print(f"\n🔬 提交自洽场计算任务 (用户: {user_id})...")
        
        # 构建请求数据
        data = {
            "user_id": user_id,
            **kwargs
        }
        
        print(f"请求数据: {json.dumps(data, indent=2, ensure_ascii=False)}")
        
        try:
            response = self.session.post(
                f"{self.base_url}/vasp/scf-calculation",
                json=data,
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 自洽场计算任务提交成功:")
            print(f"   任务ID: {result['task_id']}")
            print(f"   状态: {result['status']}")
            print(f"   消息: {result['message']}")
            
            return result['task_id']
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 自洽场计算任务提交失败: {e}")
            return None
    
    def submit_dos_calculation(self, user_id: str, **kwargs) -> Optional[str]:
        """提交态密度计算任务"""
        print(f"\n📊 提交态密度计算任务 (用户: {user_id})...")
        
        # 构建请求数据
        data = {
            "user_id": user_id,
            **kwargs
        }
        
        print(f"请求数据: {json.dumps(data, indent=2, ensure_ascii=False)}")
        
        try:
            response = self.session.post(
                f"{self.base_url}/vasp/dos-calculation",
                json=data,
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 态密度计算任务提交成功:")
            print(f"   任务ID: {result['task_id']}")
            print(f"   状态: {result['status']}")
            print(f"   消息: {result['message']}")
            
            return result['task_id']
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 态密度计算任务提交失败: {e}")
            return None
    
    def submit_md_calculation(self, user_id: str, **kwargs) -> Optional[str]:
        """提交分子动力学计算任务"""
        print(f"\n🧬 提交分子动力学计算任务 (用户: {user_id})...")
        
        # 构建请求数据
        data = {
            "user_id": user_id,
            **kwargs
        }
        
        print(f"请求数据: {json.dumps(data, indent=2, ensure_ascii=False)}")
        
        try:
            response = self.session.post(
                f"{self.base_url}/vasp/md-calculation",
                json=data,
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 分子动力学计算任务提交成功:")
            print(f"   任务ID: {result['task_id']}")
            print(f"   状态: {result['status']}")
            print(f"   消息: {result['message']}")
            
            return result['task_id']
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 分子动力学计算任务提交失败: {e}")
            return None
    
    def get_task_status(self, task_id: str, user_id: str) -> Optional[Dict[str, Any]]:
        """获取任务状态"""
        print(f"\n📊 查询任务状态: {task_id[:8]}...")
        
        try:
            response = self.session.get(
                f"{self.base_url}/vasp/task/{task_id}",
                params={"user_id": user_id},
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 任务状态:")
            print(f"   ID: {result['task_id'][:8]}...")
            print(f"   状态: {result['status']}")
            print(f"   进度: {result['progress']}%")
            print(f"   类型: {result['task_type']}")
            if result.get('process_id'):
                print(f"   进程ID: {result['process_id']}")
            if result.get('error_message'):
                print(f"   错误信息: {result['error_message']}")
            if result.get('result_path'):
                print(f"   结果路径: {result['result_path']}")
            print(f"   创建时间: {result['created_at']}")
            
            return result
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 查询状态失败: {e}")
            return None
    
    def list_user_tasks(self, user_id: str) -> Optional[list]:
        """列出用户任务"""
        print(f"\n📋 列出用户任务 (用户: {user_id})...")
        
        try:
            response = self.session.get(
                f"{self.base_url}/vasp/tasks",
                params={"user_id": user_id},
                timeout=self.timeout
            )
            response.raise_for_status()
            tasks = response.json()
            
            print(f"✅ 找到 {len(tasks)} 个任务:")
            for i, task in enumerate(tasks, 1):
                pid_info = f" | PID:{task['process_id']}" if task.get('process_id') else ""
                print(f"   {i}. {task['task_id'][:8]}... | {task['status']} | {task['progress']}% | {task['task_type']}{pid_info}")
            
            return tasks
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 列出任务失败: {e}")
            return None
    
    def cancel_task(self, task_id: str, user_id: str) -> bool:
        """取消任务"""
        print(f"\n🚫 取消任务: {task_id[:8]}...")
        
        try:
            response = self.session.post(
                f"{self.base_url}/vasp/task/{task_id}/cancel",
                params={"user_id": user_id},
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ {result['message']}")
            return True
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return False
        except Exception as e:
            print(f"❌ 取消任务失败: {e}")
            return False
    
    def get_task_result(self, task_id: str, user_id: str) -> Optional[Dict[str, Any]]:
        """获取任务结果"""
        print(f"\n📁 获取任务结果: {task_id[:8]}...")
        
        try:
            response = self.session.get(
                f"{self.base_url}/vasp/task/{task_id}/result",
                params={"user_id": user_id},
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 任务结果:")
            print(f"   结果路径: {result.get('result_path')}")
            print(f"   消息: {result.get('message')}")
            
            return result
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 获取结果失败: {e}")
            return None
    
    def get_md_result(self, task_id: str, user_id: str) -> Optional[Dict[str, Any]]:
        """获取分子动力学计算的详细结果"""
        print(f"\n🧬 获取分子动力学计算结果: {task_id[:8]}...")
        
        try:
            response = self.session.get(
                f"{self.base_url}/vasp/task/{task_id}/md-result",
                params={"user_id": user_id},
                timeout=self.timeout
            )
            response.raise_for_status()
            result = response.json()
            
            print(f"✅ 分子动力学计算结果:")
            print(f"   初始结构: {result.get('md_structure')}")
            print(f"   轨迹文件: {result.get('xdatcar_path')}")
            print(f"   能量文件: {result.get('oszicar_path')}")
            print(f"   最终能量: {result.get('final_energy')} eV")
            print(f"   平均温度: {result.get('average_temperature')} K")
            print(f"   MD步数: {result.get('total_md_steps')}")
            print(f"   是否收敛: {result.get('convergence')}")
            print(f"   计算时间: {result.get('computation_time')} 秒")
            
            return result
            
        except requests.exceptions.HTTPError as e:
            print(f"❌ HTTP错误 {e.response.status_code}: {e.response.text}")
            return None
        except Exception as e:
            print(f"❌ 获取MD结果失败: {e}")
            return None
    
    def monitor_task(self, task_id: str, user_id: str, max_time: int = 600, interval: int = 10):
        """监控任务执行"""
        print(f"\n👀 开始监控任务 {task_id[:8]}... (最大等待时间: {max_time}秒)")
        
        start_time = time.time()
        last_progress = -1
        
        while time.time() - start_time < max_time:
            status = self.get_task_status(task_id, user_id)
            
            if not status:
                print("⚠️ 无法获取任务状态，停止监控")
                break
            
            current_progress = status['progress']
            if current_progress != last_progress:
                print(f"📈 进度更新: {current_progress}% | 状态: {status['status']}")
                last_progress = current_progress
            
            if status['status'] in ['completed', 'failed', 'canceled']:
                print(f"🏁 任务结束: {status['status']}")
                if status['status'] == 'completed':
                    self.get_task_result(task_id, user_id)
                elif status['status'] == 'failed':
                    print(f"错误信息: {status.get('error_message', '未知错误')}")
                break
            
            time.sleep(interval)
        else:
            print(f"⏰ 监控超时 ({max_time}秒)")


def test_formula_submission():
    """测试化学式提交"""
    print("\n" + "="*60)
    print("🧪 测试1: 化学式提交 (LiFePO4)")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    
    if not debugger.test_connection():
        return None
    
    # 提交LiFePO4优化任务
    task_id = debugger.submit_structure_optimization(
        user_id="test_user_001",
        formula="LiFePO4",
        calc_type="OXC",
        stable_only=True,
        max_energy_above_hull=0.1,
        selection_mode="auto",
        kpoint_density=30.0
    )
    
    return task_id


def test_cif_url_submission():
    """测试CIF URL提交"""
    print("\n" + "="*60)
    print("🌐 测试2: CIF URL提交")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    
    # 使用Materials Project的CIF URL (示例)
    cif_url = "https://contribs-api.materialsproject.org/projects/dtu/structures/mp-19017.cif"
    
    task_id = debugger.submit_structure_optimization(
        user_id="test_user_002", 
        cif_url=cif_url,
        calc_type="ORC",
        kpoint_density=25.0
    )
    
    return task_id


def test_task_management():
    """测试任务管理功能"""
    print("\n" + "="*60)
    print("📊 测试3: 任务管理功能")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "test_user_001"
    
    # 列出所有任务
    tasks = debugger.list_user_tasks(user_id)
    
    if tasks:
        # 查询第一个任务的详细状态
        first_task_id = tasks[0]['task_id']
        debugger.get_task_status(first_task_id, user_id)
        
        # 如果任务还在运行，可以选择取消
        if tasks[0]['status'] in ['queued', 'running']:
            print(f"\n是否取消任务 {first_task_id[:8]}...? (y/N): ", end="")
            if input().lower() == 'y':
                debugger.cancel_task(first_task_id, user_id)


def test_quick_submission():
    """快速测试提交"""
    print("\n" + "="*60)
    print("⚡ 快速测试: Li2O 提交")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    
    task_id = debugger.submit_structure_optimization(
        user_id="debug_user",
        formula="Li2O",
        calc_type="SSE",
        stable_only=True,
        kpoint_density=20.0
    )
    
    if task_id:
        # 监控5分钟
        debugger.monitor_task(task_id, "debug_user", max_time=300, interval=15)
    
    return task_id


def test_scf_from_formula():
    """测试从化学式进行自洽场计算"""
    print("\n" + "="*60)
    print("🔬 测试: 从化学式进行自洽场计算 (Li2O)")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    
    if not debugger.test_connection():
        return None
    
    # 提交自洽场计算任务
    task_id = debugger.submit_scf_calculation(
        user_id="scf_test_user",
        formula="Li2O",
        calc_type="SSE",
        stable_only=True,
        precision="Accurate",
        kpoint_density=25.0
    )
    
    return task_id


def test_scf_from_optimization():
    """测试基于结构优化结果进行自洽场计算"""
    print("\n" + "="*60)
    print("🔗 测试: 基于结构优化结果进行自洽场计算")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "workflow_test_user"
    
    if not debugger.test_connection():
        return None, None
    
    # 第一步：提交结构优化
    print("第一步：提交结构优化任务...")
    # opt_task_id = debugger.submit_structure_optimization(
    #     user_id=user_id,
    #     formula="LiFePO4",
    #     calc_type="OXC",
    #     stable_only=True,
    #     kpoint_density=20.0
    # )
    opt_task_id = "7daaf464596d4adc9d4a82c9b5a1ba9b"
    if not opt_task_id:
        print("❌ 结构优化任务提交失败")
        return None, None
    
    print(f"✅ 结构优化任务已提交: {opt_task_id}")
    print("注意：实际使用中需要等待结构优化完成后再提交自洽场计算")
    
    # 第二步：提交自洽场计算（基于优化结果）
    print("\n第二步：提交基于优化结果的自洽场计算...")
    scf_task_id = debugger.submit_scf_calculation(
        user_id=user_id,
        optimized_task_id=opt_task_id,
        calc_type="OXC",
        precision="Accurate",
        kpoint_density=30.0
    )
    
    if scf_task_id:
        print(f"✅ 自洽场计算任务已提交: {scf_task_id}")
    else:
        print("❌ 自洽场计算任务提交失败")
    
    return opt_task_id, scf_task_id


def test_scf_workflow():
    """测试完整的工作流程：优化→自洽场"""
    print("\n" + "="*60)
    print("🔄 测试: 完整工作流程 (优化→自洽场)")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "workflow_user"
    
    # 列出现有任务，找到已完成的优化任务
    tasks = debugger.list_user_tasks(user_id)
    
    completed_opt_task = None
    if tasks:
        for task in tasks:
            if (task['task_type'] == 'structure_optimization' and 
                task['status'] == 'completed'):
                completed_opt_task = task['task_id']
                break
    
    if completed_opt_task:
        print(f"📋 找到已完成的优化任务: {completed_opt_task[:8]}...")
        
        # 基于已完成的优化任务进行自洽场计算
        scf_task_id = debugger.submit_scf_calculation(
            user_id=user_id,
            optimized_task_id=completed_opt_task,
            calc_type="OXC",
            precision="Accurate"
        )
        
        if scf_task_id:
            print("🎉 自洽场计算工作流程启动成功!")
            # 监控自洽场计算
            debugger.monitor_task(scf_task_id, user_id, max_time=600)
        
        return scf_task_id
    else:
        print("⚠️ 没有找到已完成的结构优化任务")
        print("建议先运行结构优化任务并等待完成")
        return None


def test_dos_from_scf():
    """测试基于自洽场计算结果进行态密度计算"""
    print("\n" + "="*60)
    print("📊 测试: 基于自洽场计算结果进行态密度计算")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "dos_test_user"
    
    if not debugger.test_connection():
        return None
    
    # 查找已完成的自洽场计算任务
    tasks = debugger.list_user_tasks(user_id)
    
    completed_scf_task = None
    if tasks:
        for task in tasks:
            if (task['task_type'] == 'scf_calculation' and 
                task['status'] == 'completed'):
                completed_scf_task = task['task_id']
                break
    
    if completed_scf_task:
        print(f"📋 找到已完成的自洽场任务: {completed_scf_task[:8]}...")
        
        # 基于自洽场计算进行态密度计算
        dos_task_id = debugger.submit_dos_calculation(
            user_id=user_id,
            scf_task_id=completed_scf_task,
            calc_type="OXC",
            kpoint_multiplier=2.0,
            precision="Accurate"
        )
        
        if dos_task_id:
            print("🎉 态密度计算任务提交成功!")
            # 监控态密度计算
            debugger.monitor_task(dos_task_id, user_id, max_time=600)
        
        return dos_task_id
    else:
        print("⚠️ 没有找到已完成的自洽场计算任务")
        print("建议先运行自洽场计算任务并等待完成")
        return None


def test_single_point_dos():
    """测试单点自洽+DOS计算（一步完成）"""
    print("\n" + "="*60)
    print("⚡ 测试: 单点自洽+DOS计算 (Li2O) - 一步搞定")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "single_dos_user"
    
    if not debugger.test_connection():
        return None
    
    # 直接从化学式进行单点自洽+DOS计算（一次VASP运行完成）
    dos_task_id = debugger.submit_dos_calculation(
        user_id=user_id,
        formula="Li2O",
        calc_type="SSE",
        kpoint_multiplier=2.5,
        precision="Accurate",
        stable_only=True,
        selection_mode="most_stable"
    )
    
    if dos_task_id:
        print("🎉 单点自洽+DOS计算任务提交成功!")
        print("⚡ 该任务特点：")
        print("   • 一次VASP运行完成自洽场+DOS计算")
        print("   • INCAR同时包含自洽场和DOS设置")
        print("   • 无需分步操作，一步搞定")
        print("📊 执行流程：")
        print("   1. 下载Li2O的CIF文件")
        print("   2. 转换为POSCAR格式") 
        print("   3. 生成包含DOS设置的INCAR")
        print("   4. 一次运行VASP得到自洽场+DOS结果")
        
        # 监控任务
        debugger.monitor_task(dos_task_id, user_id, max_time=1200)
        
        return dos_task_id
    else:
        print("❌ 单点自洽+DOS计算任务提交失败")
        return None


def test_full_workflow():
    """测试完整工作流程：优化→自洽场→态密度"""
    print("\n" + "="*60)
    print("🔄 测试: 完整工作流程 (优化→自洽场→态密度)")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "workflow_dos_user"
    
    if not debugger.test_connection():
        return None, None, None
    
    print("第一步：查找或提交结构优化任务...")
    
    # 查找已完成的优化任务
    tasks = debugger.list_user_tasks(user_id)
    completed_opt_task = None
    
    if tasks:
        for task in tasks:
            if (task['task_type'] == 'structure_optimization' and 
                task['status'] == 'completed'):
                completed_opt_task = task['task_id']
                break
    
    if not completed_opt_task:
        print("📤 提交新的结构优化任务...")
        completed_opt_task = debugger.submit_structure_optimization(
            user_id=user_id,
            formula="Li2O",
            calc_type="SSE",
            stable_only=True,
            kpoint_density=20.0
        )
        if not completed_opt_task:
            print("❌ 结构优化任务提交失败")
            return None, None, None
        print("⚠️ 需要等待结构优化完成后再进行后续步骤")
        return completed_opt_task, None, None
    
    print(f"✅ 找到已完成的优化任务: {completed_opt_task[:8]}...")
    
    print("\n第二步：查找或提交自洽场计算任务...")
    
    # 查找基于该优化任务的自洽场计算
    completed_scf_task = None
    if tasks:  # 确保tasks不为None
        for task in tasks:
            if (task['task_type'] == 'scf_calculation' and 
                task['status'] == 'completed' and
                task.get('params', {}).get('optimized_task_id') == completed_opt_task):
                completed_scf_task = task['task_id']
                break
    
    if not completed_scf_task:
        print("📤 基于优化结果提交自洽场计算...")
        completed_scf_task = debugger.submit_scf_calculation(
            user_id=user_id,
            optimized_task_id=completed_opt_task,
            calc_type="SSE",
            precision="Accurate"
        )
        if not completed_scf_task:
            print("❌ 自洽场计算任务提交失败")
            return completed_opt_task, None, None
        print("⚠️ 需要等待自洽场计算完成后再进行态密度计算")
        return completed_opt_task, completed_scf_task, None
    
    print(f"✅ 找到已完成的自洽场任务: {completed_scf_task[:8]}...")
    
    print("\n第三步：提交态密度计算...")
    
    dos_task_id = debugger.submit_dos_calculation(
        user_id=user_id,
        scf_task_id=completed_scf_task,
        calc_type="SSE",
        kpoint_multiplier=2.5,
        precision="Accurate"
    )
    
    if dos_task_id:
        print("🎉 完整工作流程启动成功!")
        print("🔍 监控态密度计算...")
        debugger.monitor_task(dos_task_id, user_id, max_time=900)
    
    return completed_opt_task, completed_scf_task, dos_task_id


def test_md_from_scf():
    """测试从自洽场计算结果运行分子动力学计算"""
    print("\n" + "="*60)
    print("🧬 测试: 从自洽场任务运行分子动力学计算")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "md_scf_user"
    
    if not debugger.test_connection():
        return None
    
    # 查找已完成的自洽场任务
    tasks = debugger.list_user_tasks(user_id)
    completed_scf_task = None
    
    if tasks:
        for task in tasks:
            if (task['task_type'] == 'scf_calculation' and 
                task['status'] == 'completed'):
                completed_scf_task = task['task_id']
                break
    
    if not completed_scf_task:
        print("⚠️ 未找到已完成的自洽场任务，先提交一个:")
        scf_task_id = debugger.submit_scf_calculation(
            user_id=user_id,
            formula="Li2O",
            calc_type="SSE",
            precision="Normal",
            stable_only=True,
            selection_mode="most_stable"
        )
        
        if scf_task_id:
            print("⏳ 等待自洽场计算完成...")
            debugger.monitor_task(scf_task_id, user_id, max_time=600)
            completed_scf_task = scf_task_id
        else:
            print("❌ 自洽场任务提交失败")
            return None
    
    print(f"✅ 找到已完成的自洽场任务: {completed_scf_task[:8]}...")
    
    # 提交分子动力学计算
    md_task_id = debugger.submit_md_calculation(
        user_id=user_id,
        scf_task_id=completed_scf_task,
        calc_type="SSE",
        md_steps=500,
        temperature=300.0,
        time_step=1.0,
        ensemble="NVT",
        precision="Normal"
    )
    
    if md_task_id:
        print("🎉 分子动力学计算任务提交成功!")
        print("🔍 监控分子动力学计算...")
        debugger.monitor_task(md_task_id, user_id, max_time=1800)
        
        # 获取MD结果
        print("\n🔍 获取分子动力学计算结果...")
        md_result = debugger.get_md_result(md_task_id, user_id)
        
        return md_task_id
    else:
        print("❌ 分子动力学计算任务提交失败")
        return None


def test_single_point_md():
    """测试单点自洽+MD计算（一步完成）"""
    print("\n" + "="*60)
    print("⚡ 测试: 单点自洽+MD计算 (Li2O) - 一步搞定")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "single_md_user"
    
    if not debugger.test_connection():
        return None
    
    # 直接从化学式进行单点自洽+MD计算（一次VASP运行完成）
    md_task_id = debugger.submit_md_calculation(
        user_id=user_id,
        formula="Li2O",
        calc_type="SSE",
        md_steps=1000,
        temperature=300.0,
        time_step=1.0,
        ensemble="NVT",
        precision="Normal",
        stable_only=True,
        selection_mode="most_stable"
    )
    
    if md_task_id:
        print("🎉 单点自洽+MD计算任务提交成功!")
        print("⚡ 该任务特点：")
        print("   • 一次VASP运行完成自洽场+MD计算")
        print("   • INCAR同时包含自洽场和MD设置")
        print("   • 无需分步操作，一步搞定")
        print("📊 执行流程：")
        print("   1. 下载Li2O的CIF文件")
        print("   2. 转换为POSCAR格式") 
        print("   3. 生成包含MD设置的INCAR")
        print("   4. 一次运行VASP得到自洽场+MD结果")
        
        # 监控任务
        debugger.monitor_task(md_task_id, user_id, max_time=2400)
        
        # 获取MD结果
        print("\n🔍 获取分子动力学计算结果...")
        md_result = debugger.get_md_result(md_task_id, user_id)
        
        return md_task_id
    else:
        print("❌ 单点自洽+MD计算任务提交失败")
        return None


def test_full_md_workflow():
    """测试完整工作流程：优化→自洽场→分子动力学"""
    print("\n" + "="*60)
    print("🔄 测试: 完整工作流程 (优化→自洽场→分子动力学)")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    user_id = "workflow_md_user"
    
    if not debugger.test_connection():
        return None, None, None
    
    print("第一步：查找或提交结构优化任务...")
    
    # 查找已完成的优化任务
    tasks = debugger.list_user_tasks(user_id)
    completed_opt_task = None
    
    if tasks:
        for task in tasks:
            if (task['task_type'] == 'structure_optimization' and 
                task['status'] == 'completed'):
                completed_opt_task = task['task_id']
                break
    
    if not completed_opt_task:
        print("⚠️ 未找到已完成的优化任务，先提交一个:")
        opt_task_id = debugger.submit_structure_optimization(
            user_id=user_id,
            formula="Li2O",
            calc_type="OXC",
            precision="Normal",
            stable_only=True,
            selection_mode="most_stable"
        )
        
        if opt_task_id:
            print("⏳ 等待结构优化完成...")
            debugger.monitor_task(opt_task_id, user_id, max_time=900)
            completed_opt_task = opt_task_id
        else:
            print("❌ 结构优化任务提交失败")
            return None, None, None
    
    print(f"✅ 找到已完成的优化任务: {completed_opt_task[:8]}...")
    
    print("\n第二步：查找或提交自洽场计算...")
    
    # 查找已完成的自洽场任务
    completed_scf_task = None
    
    if tasks:
        for task in tasks:
            if (task['task_type'] == 'scf_calculation' and 
                task['status'] == 'completed'):
                completed_scf_task = task['task_id']
                break
    
    if not completed_scf_task:
        print("⚠️ 未找到已完成的自洽场任务，先提交一个:")
        scf_task_id = debugger.submit_scf_calculation(
            user_id=user_id,
            optimized_task_id=completed_opt_task,
            calc_type="SSE",
            precision="Normal"
        )
        
        if scf_task_id:
            print("⏳ 等待自洽场计算完成...")
            debugger.monitor_task(scf_task_id, user_id, max_time=600)
            completed_scf_task = scf_task_id
        else:
            print("❌ 自洽场任务提交失败")
            return completed_opt_task, None, None
    
    print(f"✅ 找到已完成的自洽场任务: {completed_scf_task[:8]}...")
    
    print("\n第三步：提交分子动力学计算...")
    
    md_task_id = debugger.submit_md_calculation(
        user_id=user_id,
        scf_task_id=completed_scf_task,
        calc_type="SSE",
        md_steps=1000,
        temperature=300.0,
        time_step=1.0,
        ensemble="NVT",
        precision="Normal"
    )
    
    if md_task_id:
        print("🎉 完整工作流程启动成功!")
        print("🔍 监控分子动力学计算...")
        debugger.monitor_task(md_task_id, user_id, max_time=1800)
        
        # 获取MD结果
        print("\n🔍 获取分子动力学计算结果...")
        md_result = debugger.get_md_result(md_task_id, user_id)
    
    return completed_opt_task, completed_scf_task, md_task_id


def test_error_cases():
    """测试错误情况"""
    print("\n" + "="*60)
    print("❌ 测试4: 错误情况处理")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    
    # 测试1: 既不提供formula也不提供cif_url
    print("\n测试错误1: 缺少必要参数")
    debugger.submit_structure_optimization(
        user_id="test_user",
        calc_type="OXC"
    )
    
    # 测试2: 同时提供formula和cif_url
    print("\n测试错误2: 同时提供formula和cif_url")
    debugger.submit_structure_optimization(
        user_id="test_user",
        formula="Li2O",
        cif_url="https://example.com/test.cif",
        calc_type="OXC"
    )
    
    # 测试3: 查询不存在的任务
    print("\n测试错误3: 查询不存在的任务")
    debugger.get_task_status("nonexistent_task_id", "test_user")


def interactive_mode():
    """交互式调试模式"""
    print("\n" + "="*60)
    print("🎮 交互式调试模式")
    print("="*60)
    
    debugger = VASPAPIDebugger()
    
    if not debugger.test_connection():
        return
    
    while True:
        print("\n请选择操作:")
        print("1. 提交LiFePO4优化任务")
        print("2. 提交Li2O优化任务")
        print("3. 提交Li2O自洽场计算")
        print("4. 基于优化结果进行自洽场计算")
        print("5. 基于自洽场结果进行态密度计算")
        print("6. 提交Li2O单点自洽+DOS计算 (一步搞定)")
        print("7. 基于自洽场结果进行分子动力学计算")
        print("8. 提交Li2O单点自洽+MD计算 (一步搞定)")
        print("9. 列出所有任务")
        print("a. 查询任务状态")
        print("b. 取消任务")
        print("c. 获取任务结果")
        print("d. 获取MD计算结果")
        print("e. 监控任务")
        print("0. 退出")
        
        choice = input("\n请输入选择 (0-9, a-e): ").strip()
        
        if choice == "0":
            print("👋 退出调试程序")
            break
        elif choice == "1":
            task_id = debugger.submit_structure_optimization(
                user_id="interactive_user",
                formula="LiFePO4",
                calc_type="OXC",
                stable_only=True
            )
            if task_id:
                print(f"任务ID: {task_id}")
        elif choice == "2":
            task_id = debugger.submit_structure_optimization(
                user_id="interactive_user",
                formula="Li2O",
                calc_type="SSE",
                stable_only=True
            )
            if task_id:
                print(f"任务ID: {task_id}")
        elif choice == "3":
            task_id = debugger.submit_scf_calculation(
                user_id="interactive_user",
                formula="Li2O",
                calc_type="SSE",
                precision="Accurate",
                stable_only=True
            )
            if task_id:
                print(f"自洽场计算任务ID: {task_id}")
        elif choice == "4":
            # 列出任务找到已完成的优化任务
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            tasks = debugger.list_user_tasks(user_id)
            
            completed_opt_tasks = []
            if tasks:
                for task in tasks:
                    if (task['task_type'] == 'structure_optimization' and 
                        task['status'] == 'completed'):
                        completed_opt_tasks.append(task)
            
            if completed_opt_tasks:
                print("\n已完成的结构优化任务:")
                for i, task in enumerate(completed_opt_tasks, 1):
                    print(f"  {i}. {task['task_id'][:8]}... ({task.get('params', {}).get('formula', 'unknown')})")
                
                try:
                    idx = int(input("选择优化任务序号: ").strip()) - 1
                    if 0 <= idx < len(completed_opt_tasks):
                        opt_task_id = completed_opt_tasks[idx]['task_id']
                        
                        scf_task_id = debugger.submit_scf_calculation(
                            user_id=user_id,
                            optimized_task_id=opt_task_id,
                            calc_type="OXC",
                            precision="Accurate"
                        )
                        if scf_task_id:
                            print(f"自洽场计算任务ID: {scf_task_id}")
                    else:
                        print("❌ 无效选择")
                except ValueError:
                    print("❌ 请输入有效数字")
            else:
                print("❌ 没有找到已完成的结构优化任务")
        elif choice == "5":
            # 基于自洽场结果进行态密度计算
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            tasks = debugger.list_user_tasks(user_id)
            
            completed_scf_tasks = []
            if tasks:
                for task in tasks:
                    if (task['task_type'] == 'scf_calculation' and 
                        task['status'] == 'completed'):
                        completed_scf_tasks.append(task)
            
            if completed_scf_tasks:
                print("\n已完成的自洽场计算任务:")
                for i, task in enumerate(completed_scf_tasks, 1):
                    formula = task.get('params', {}).get('formula', 'unknown')
                    print(f"  {i}. {task['task_id'][:8]}... ({formula})")
                
                try:
                    idx = int(input("选择自洽场任务序号: ").strip()) - 1
                    if 0 <= idx < len(completed_scf_tasks):
                        scf_task_id = completed_scf_tasks[idx]['task_id']
                        
                        dos_task_id = debugger.submit_dos_calculation(
                            user_id=user_id,
                            scf_task_id=scf_task_id,
                            calc_type="OXC",
                            kpoint_multiplier=2.0,
                            precision="Accurate"
                        )
                        if dos_task_id:
                            print(f"态密度计算任务ID: {dos_task_id}")
                    else:
                        print("❌ 无效选择")
                except ValueError:
                    print("❌ 请输入有效数字")
            else:
                print("❌ 没有找到已完成的自洽场计算任务")
        elif choice == "6":
            # 单点自洽+DOS计算
            task_id = debugger.submit_dos_calculation(
                user_id="interactive_user",
                formula="Li2O",
                calc_type="SSE",
                kpoint_multiplier=2.5,
                precision="Accurate",
                stable_only=True
            )
            if task_id:
                print(f"单点自洽+DOS计算任务ID: {task_id}")
                print("⚡ 该任务会在一次VASP运行中完成自洽场和DOS计算")
        elif choice == "7":
            # 基于自洽场结果进行分子动力学计算
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            tasks = debugger.list_user_tasks(user_id)
            
            completed_scf_tasks = []
            if tasks:
                for task in tasks:
                    if (task['task_type'] == 'scf_calculation' and 
                        task['status'] == 'completed'):
                        completed_scf_tasks.append(task)
            
            if completed_scf_tasks:
                print("\n已完成的自洽场计算任务:")
                for i, task in enumerate(completed_scf_tasks, 1):
                    print(f"  {i}. {task['task_id'][:8]}... ({task.get('params', {}).get('formula', 'unknown')})")
                
                try:
                    idx = int(input("选择自洽场任务序号: ").strip()) - 1
                    if 0 <= idx < len(completed_scf_tasks):
                        scf_task_id = completed_scf_tasks[idx]['task_id']
                        
                        md_task_id = debugger.submit_md_calculation(
                            user_id=user_id,
                            scf_task_id=scf_task_id,
                            calc_type="SSE",
                            md_steps=500,
                            temperature=300.0,
                            time_step=1.0,
                            ensemble="NVT",
                            precision="Normal"
                        )
                        if md_task_id:
                            print(f"分子动力学计算任务ID: {md_task_id}")
                    else:
                        print("❌ 无效选择")
                except ValueError:
                    print("❌ 请输入有效数字")
            else:
                print("❌ 没有找到已完成的自洽场计算任务")
        elif choice == "8":
            # 单点自洽+MD计算
            task_id = debugger.submit_md_calculation(
                user_id="interactive_user",
                formula="Li2O",
                calc_type="SSE",
                md_steps=1000,
                temperature=300.0,
                time_step=1.0,
                ensemble="NVT",
                precision="Normal",
                stable_only=True
            )
            if task_id:
                print(f"单点自洽+MD计算任务ID: {task_id}")
                print("⚡ 该任务会在一次VASP运行中完成自洽场和MD计算")
        elif choice == "9":
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            debugger.list_user_tasks(user_id)
        elif choice == "a":
            task_id = input("输入任务ID: ").strip()
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            debugger.get_task_status(task_id, user_id)
        elif choice == "b":
            task_id = input("输入要取消的任务ID: ").strip()
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            debugger.cancel_task(task_id, user_id)
        elif choice == "c":
            task_id = input("输入任务ID: ").strip()
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            debugger.get_task_result(task_id, user_id)
        elif choice == "d":
            task_id = input("输入MD任务ID: ").strip()
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            debugger.get_md_result(task_id, user_id)
        elif choice == "e":
            task_id = input("输入要监控的任务ID: ").strip()
            user_id = input("输入用户ID (默认: interactive_user): ").strip() or "interactive_user"
            max_time = int(input("最大监控时间(秒, 默认300): ").strip() or "300")
            debugger.monitor_task(task_id, user_id, max_time)
        else:
            print("❌ 无效选择，请重试")


def main():
    """主函数"""
    print("🚀 VASP API 调试程序")
    print("=" * 60)
    print("目标服务: http://localhost:9000")
    print("=" * 60)
    
    if len(sys.argv) > 1:
        test_type = sys.argv[1].lower()
        
        if test_type == "formula":
            test_formula_submission()
        elif test_type == "url":
            test_cif_url_submission()
        elif test_type == "manage":
            test_task_management()
        elif test_type == "quick":
            test_quick_submission()
        elif test_type == "scf":
            test_scf_from_formula()
        elif test_type == "scf_opt":
            test_scf_from_optimization()
        elif test_type == "workflow":
            test_scf_workflow()
        elif test_type == "dos":
            test_dos_from_scf()
        elif test_type == "single_dos":
            test_single_point_dos()
        elif test_type == "full_workflow":
            test_full_workflow()
        elif test_type == "md":
            test_md_from_scf()
        elif test_type == "single_md":
            test_single_point_md()
        elif test_type == "full_md_workflow":
            test_full_md_workflow()
        elif test_type == "error":
            test_error_cases()
        elif test_type == "interactive":
            interactive_mode()
        else:
            print(f"❌ 未知测试类型: {test_type}")
            print("可用类型: formula, url, manage, quick, scf, scf_opt, workflow, dos, single_dos, full_workflow, md, single_md, full_md_workflow, error, interactive")
    else:
        # 默认运行所有测试
        print("运行所有基础测试...")
        
        # 基础连接测试
        debugger = VASPAPIDebugger()
        if not debugger.test_connection():
            return
        
        # 运行基础测试
        test_formula_submission()
        test_task_management()
        test_error_cases()
        
        # 询问是否进入交互模式
        print(f"\n是否进入交互式调试模式? (y/N): ", end="")
        if input().lower() == 'y':
            interactive_mode()


if __name__ == "__main__":
    main()


# =============================================================================
# 使用说明
# =============================================================================
"""
使用方法:

1. 运行所有测试:
   python debug_vasp_api.py

2. 运行特定测试:
   python debug_vasp_api.py formula     # 测试化学式提交（结构优化）
   python debug_vasp_api.py url         # 测试CIF URL提交
   python debug_vasp_api.py manage      # 测试任务管理
   python debug_vasp_api.py quick       # 快速测试
   python debug_vasp_api.py scf         # 测试自洽场计算（从化学式）
   python debug_vasp_api.py scf_opt     # 测试自洽场计算（基于优化结果）
   python debug_vasp_api.py workflow    # 测试部分工作流程（优化→自洽场）
   python debug_vasp_api.py dos         # 测试态密度计算（基于自洽场结果）
   python debug_vasp_api.py single_dos  # 测试单点自洽+DOS计算（一步搞定）
   python debug_vasp_api.py full_workflow  # 测试完整工作流程（优化→自洽场→态密度）
   python debug_vasp_api.py md          # 测试分子动力学计算（基于自洽场结果）
   python debug_vasp_api.py single_md   # 测试单点自洽+MD计算（一步搞定）
   python debug_vasp_api.py full_md_workflow  # 测试完整MD工作流程（优化→自洽场→分子动力学）
   python debug_vasp_api.py error       # 测试错误处理
   python debug_vasp_api.py interactive # 交互式模式

3. 交互式模式功能:
   - 选项1-6: 各种基础计算任务（结构优化、自洽场、态密度）
   - 选项7-8: 分子动力学计算（基于已有结果或一步完成）
   - 选项9: 列出所有任务
   - 选项a-e: 任务管理（状态查询、取消、结果获取、监控等）

4. 分子动力学计算特点:
   - 支持NVT、NVE、NPT三种系综
   - 可配置温度、时间步长、MD步数
   - 输出XDATCAR轨迹文件和OSZICAR能量文件
   - 支持基于自洽场结果或单点自洽+MD一步完成

5. 常见问题排查:
   - 连接失败: 检查端口转发和服务状态
   - HTTP 500错误: 检查超算上的服务日志
   - 任务一直pending: 检查VASP环境和依赖

6. 端口转发设置 (示例):
   ssh -L 9000:localhost:8000 username@supercomputer
""" 