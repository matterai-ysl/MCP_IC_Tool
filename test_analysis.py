"""
测试结构弛豫结果分析 - 直接调用分析方法
"""
import sys
import os
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

def test_analyze():
    """测试分析已完成的结构弛豫结果"""

    # 测试文件夹路径
    work_dir = Path("/Users/ysl/Desktop/Code/MCP_IC_Tool/02f98add7eb74967b436ce512725b184")

    # 模拟 vasp_result（从运行阶段传过来的结果）
    vasp_result = {
        'computation_time': 100.0,
        'process_id': 'test_pid'
    }

    print("=" * 60)
    print("开始测试结构弛豫结果分析...")
    print(f"工作目录: {work_dir}")
    print("=" * 60)

    try:
        # 直接执行分析逻辑（模拟 _analyze_results 方法）
        from src.vasp_server.Config import get_download_url, get_static_url

        # 检查收敛性
        print("\n[1/7] 检查收敛性...")
        outcar_path = work_dir / "OUTCAR"
        with open(outcar_path, 'rb') as f:
            f.seek(-1024, os.SEEK_END)
            last_lines = f.readlines()[-10:]
            last_content = b''.join(last_lines).decode('utf-8', errors='ignore')
            convergence = 'reached required accuracy' in last_content or 'Voluntary' in last_content
        print(f"  收敛性: {convergence}")

        # 提取能量
        print("\n[2/7] 提取能量...")
        energy = None
        with open(outcar_path, 'r') as f:
            lines = f.readlines()
        for line in reversed(lines):
            if 'free energy    TOTEN' in line:
                parts = line.split()
                energy = float(parts[4])
                break
        print(f"  能量: {energy}")

        # 提取力
        print("\n[3/7] 提取力...")
        forces = None  # 简化测试，不提取
        print(f"  力: {forces}")

        # 检查优化后的结构
        print("\n[4/7] 检查优化后的结构...")
        contcar_path = work_dir / "CONTCAR"
        optimized_structure = None
        if contcar_path.exists():
            optimized_structure = str(contcar_path)
        print(f"  CONTCAR存在: {optimized_structure is not None}")

        # 生成可视化分析报告
        print("\n[5/7] 生成可视化分析报告...")
        html_report_path = None
        analysis_data = None
        try:
            from src.vasp_server.optimization_analyzer import generate_optimization_report, OUTCARAnalyzer
            if outcar_path.exists():
                analyzer = OUTCARAnalyzer(str(work_dir), task_id="optimization")
                analysis_data = analyzer.analyze()
                print(f"  分析数据生成成功: {analysis_data is not None}")

                html_report_path = generate_optimization_report(str(work_dir), "optimization")
                print(f"  HTML报告生成成功: {html_report_path}")
        except Exception as e:
            print(f"  ⚠️ 生成可视化分析报告失败: {e}")
            import traceback
            traceback.print_exc()

        # 构建结果字典
        print("\n[6/7] 构建结果字典...")

        # 安全地生成下载URL（处理路径不在DOWNLOAD_URL下的情况）
        from src.vasp_server.Config import DOWNLOAD_URL
        optimized_structure_url = None
        if optimized_structure:
            try:
                # 检查路径是否在DOWNLOAD_URL下
                Path(optimized_structure).relative_to(DOWNLOAD_URL)
                optimized_structure_url = get_download_url(optimized_structure)
                print(f"  使用下载URL: {optimized_structure_url}")
            except ValueError:
                # 路径不在DOWNLOAD_URL下，使用绝对路径
                optimized_structure_url = optimized_structure
                print(f"  使用绝对路径: {optimized_structure_url}")

        result = {
            'success': True,
            'convergence': convergence,
            'energy': energy,
            'final_forces': forces,
            'optimized_structure_download_url': optimized_structure_url,
            'computation_time': vasp_result.get('computation_time',None),
            'process_id': vasp_result.get('process_id',None),
            'work_directory': str(work_dir)
        }

        if html_report_path:
            try:
                # 检查路径是否在DOWNLOAD_URL下
                Path(html_report_path).relative_to(DOWNLOAD_URL)
                html_relative_path = get_static_url(html_report_path)
                result['analysis_report_html_path'] = html_relative_path
                print(f"  HTML报告URL: {html_relative_path}")
            except ValueError:
                # 路径不在DOWNLOAD_URL下，使用绝对路径
                result['analysis_report_html_path'] = html_report_path
                print(f"  HTML报告绝对路径: {html_report_path}")

        if analysis_data:
            result['analysis_data'] = analysis_data

        print(f"  结果字典构建成功")
        print(f"  'analysis_data' 是否在结果中: {'analysis_data' in result}")

        if 'analysis_data' in result:
            print(f"\n  analysis_data 的键: {list(result['analysis_data'].keys())}")

        # 简化返回结果
        print("\n[7/7] 简化返回结果...")

        # 使用修复后的逻辑
        if analysis_data and 'convergence_analysis' in analysis_data:
            conv_analysis = analysis_data['convergence_analysis']
            simplified_result = {
                'success': result['success'],
                'force_convergence': conv_analysis.get('force_convergence', {}).get("converged", False),
                'final_max_force': conv_analysis.get('force_convergence', {}).get("final_max_force", None),
                'energy_convergence': conv_analysis.get('energy_convergence', {}).get("converged", False),
                'final_energy': conv_analysis.get('energy_convergence', {}).get("final_energy", None),
                'optimized_structure_download_url': result['optimized_structure_download_url'],
                'computation_time': result['computation_time'],
                'analysis_report_html_path': result.get('analysis_report_html_path', None),
            }
            result = simplified_result
        else:
            print("  ⚠️ 没有convergence_analysis数据，使用基础结果")

        print("\n✅ 分析成功!")
        print("\n结果:")
        for key, value in result.items():
            print(f"  {key}: {value}")

    except Exception as e:
        print(f"\n❌ 分析失败!")
        print(f"错误类型: {type(e).__name__}")
        print(f"错误信息: {str(e)}")

        import traceback
        print("\n完整错误堆栈:")
        print(traceback.format_exc())

if __name__ == "__main__":
    test_analyze()
