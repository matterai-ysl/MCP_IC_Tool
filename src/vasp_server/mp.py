import os
import shutil
import asyncio
import traceback
from mp_api.client import MPRester
API_KEY = 'g8j3tc9BUugnPSzgJ2ppCaxPgEo5W8H7'

def _sanitize_filename(text: str) -> str:
    """
    清洗文件名中的非法字符，确保可作为文件名/路径的一部分。
    - 替换: / \ : * ? " < > | 为空或为 '-'
    - 去除首尾空格
    """
    if not isinstance(text, str):
        text = str(text)
    invalid = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']
    for ch in invalid:
        text = text.replace(ch, '-')
    return text.strip()

async def download_cif_by_formula(formula, save_path, task_id, criteria=None, selection_mode="auto"):
    """
    异步下载指定化学式的CIF文件，支持多种筛选条件
    
    Args:
        formula: 化学式，如 "LiFePO4"
        save_path: 保存路径
        task_id: 任务ID
        criteria: 额外筛选条件字典，例如:
                 {
                     'spacegroup_symbol': 'Pnma',  # 空间群
                     'energy_above_hull': (None, 0.1),  # 稳定性 (min, max)
                     'band_gap': (0.5, None),  # 带隙 (min, max)
                     'num_sites': (None, 20)  # 原子数 (min, max)
                 }
        selection_mode: 选择模式
                       - "auto": 自动选择最稳定的材料
                       - "interactive": 交互式选择
                       - "all": 下载所有匹配的材料
                       - "first": 选择第一个
    
    Returns:
        str|list|None: 成功时返回保存路径或路径列表，失败时返回None
    """
    # 保存当前工作目录
    original_dir = os.getcwd()
    
    print(f"[INFO] 开始查询化学式 {formula} 的材料...")
    print(f"[INFO] 保存路径: {save_path}")
    print(f"[INFO] 任务ID: {task_id}")
    print(f"[INFO] 选择模式: {selection_mode}")
    print(f"[INFO] 使用新版 mp-api")
    
    try:
        # 创建保存目录
        try:
            os.makedirs(save_path, exist_ok=True)
            print(f"[INFO] 保存目录已创建/确认: {save_path}")
        except OSError as e:
            print(f"[ERROR] 创建保存目录失败: {e}")
            return None
        
        # 构建查询参数
        search_params = {'formula': formula}
        
        # 添加用户自定义筛选条件
        if criteria:
            search_params.update(convert_criteria_to_search_params(criteria))
            print(f"[INFO] 应用筛选条件: {criteria}")
        
        try:
            with MPRester(API_KEY) as mpr:
                print(f"[INFO] 已连接到Materials Project API")
                
                # 执行查询
                try:
                    docs = mpr.materials.summary.search(
                        **search_params,
                        fields=['material_id', 'formula_pretty', 'symmetry', 'nsites', 
                               'energy_above_hull', 'formation_energy_per_atom', 'band_gap', 
                               'density', 'is_stable']
                    )
                
                    print(f"[INFO] 查询到 {len(docs)} 个匹配的材料")
                
                except Exception as e:
                    print(f"[ERROR] 查询失败: {e}")
                    print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
                    return None
                
                if not docs:
                    print(f"[WARNING] 未找到与条件匹配的材料")
                    return None
                
                # 显示找到的材料信息
                print(f"\n[INFO] 找到的材料:")
                for i, doc in enumerate(docs):
                    stability = "稳定" if doc.is_stable else f"E_hull={doc.energy_above_hull:.3f}eV"  # type: ignore
                    print(f"  {i+1}. {doc.material_id} - {doc.formula_pretty} "  # type: ignore
                          f"(空间群: {doc.symmetry.symbol}, "  # type: ignore
                          f"原子数: {doc.nsites}, {stability})")  # type: ignore
                
                # 根据选择模式确定要下载的材料
                selected_docs = await select_materials(docs, selection_mode)
                
                if not selected_docs:
                    print(f"[WARNING] 没有选择任何材料下载")
                    return None
                
                # 下载选中的材料
                download_results = []
                for doc in selected_docs:
                    result = await download_single_material(mpr, doc, formula, save_path)
                    if result:
                        download_results.append(result)
                
                if download_results:
                    if len(download_results) == 1:
                        return download_results[0]
                    else:
                        return download_results
                else:
                    return None
                    
        except Exception as e:
            print(f"[ERROR] API连接失败: {e}")
            print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
            return None
            
    except Exception as e:
        print(f"[ERROR] 下载过程中发生未预期的错误: {e}")
        print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
        return None
        
    finally:
        # 确保恢复原始目录
        try:
            os.chdir(original_dir)
            print(f"[INFO] 已恢复原始工作目录: {original_dir}")
        except Exception as e:
            print(f"[WARNING] 恢复工作目录时出错: {e}")

def convert_criteria_to_search_params(criteria):
    """
    将用户的筛选条件转换为新版API的搜索参数
    
    Args:
        criteria: 用户筛选条件字典
    
    Returns:
        dict: 转换后的搜索参数
    """
    search_params = {}
    
    for key, value in criteria.items():
        if key == 'spacegroup_symbol':
            search_params['spacegroup_symbol'] = value
        elif key == 'energy_above_hull':
            if isinstance(value, tuple):
                search_params['energy_above_hull'] = value
            else:
                search_params['energy_above_hull'] = (None, value)
        elif key == 'band_gap':
            if isinstance(value, tuple):
                search_params['band_gap'] = value
            else:
                search_params['band_gap'] = (value, None)
        elif key == 'nsites' or key == 'num_sites':
            # 统一转换为num_sites（作为搜索参数），但返回字段使用nsites
            rng = value if isinstance(value, tuple) else (None, value)
            search_params['num_sites'] = rng
        elif key == 'density':
            if isinstance(value, tuple):
                search_params['density'] = value
            else:
                search_params['density'] = (value, None)
        elif key == 'is_stable':
            search_params['is_stable'] = value
        elif key == 'formation_energy_per_atom':
            if isinstance(value, tuple):
                search_params['formation_energy_per_atom'] = value
            else:
                search_params['formation_energy_per_atom'] = (None, value)
        else:
            # 直接传递其他条件
            search_params[key] = value
    
    return search_params

async def select_materials(docs, selection_mode):
    """
    根据选择模式确定要下载的材料
    
    Args:
        docs: 查询到的材料文档列表
        selection_mode: 选择模式
    
    Returns:
        list: 选中的材料文档列表
    """
    if selection_mode == "auto":
        # 自动选择：优先稳定材料，其次是能量最低的
        stable_docs = [doc for doc in docs if doc.is_stable]  # type: ignore
        if stable_docs:
            # 如果有稳定材料，选择形成能最低的
            selected = min(stable_docs, key=lambda x: x.formation_energy_per_atom)  # type: ignore
            print(f"[AUTO] 自动选择稳定材料: {selected.material_id}")  # type: ignore
            return [selected]
        else:
            # 没有稳定材料，选择能量最接近凸包的
            selected = min(docs, key=lambda x: x.energy_above_hull)  # type: ignore
            print(f"[AUTO] 自动选择最接近稳定的材料: {selected.material_id} (E_hull={selected.energy_above_hull:.3f}eV)")  # type: ignore
            return [selected]
    
    elif selection_mode == "interactive":
        # 交互式选择（简化版，实际可以更复杂）
        print(f"\n[INTERACTIVE] 请选择要下载的材料:")
        print(f"输入材料编号 (1-{len(docs)})，多个用逗号分隔，或输入 'all' 下载全部:")
        
        # 在实际应用中，这里可以接收用户输入
        # 为了演示，这里使用自动选择逻辑
        print(f"[INTERACTIVE] 演示模式：自动选择最稳定的材料")
        return await select_materials(docs, "auto")
    
    elif selection_mode == "all":
        print(f"[ALL] 选择所有 {len(docs)} 个材料")
        return docs
    
    elif selection_mode == "first":
        print(f"[FIRST] 选择第一个材料: {docs[0].material_id}")  # type: ignore
        return [docs[0]]
    
    else:
        print(f"[ERROR] 未知的选择模式: {selection_mode}")
        return []

async def download_single_material(mpr, doc, formula, save_path):
    """
    下载单个材料的CIF文件
    
    Args:
        mpr: MPRester对象
        doc: 材料文档对象
        formula: 化学式
        save_path: 保存路径
    
    Returns:
        str|None: 成功时返回保存路径，失败时返回None
    """
    material_id = doc.material_id  # type: ignore
    spacegroup = doc.symmetry.symbol  # type: ignore
    
    try:
        print(f"[INFO] 开始下载材料 {material_id}...")
        
        # 获取材料结构
        structure = mpr.get_structure_by_material_id(material_id)
        print(f"[INFO] 成功获取材料 {material_id} 的结构数据")
        
        # 转换为CIF格式
        cif_data = structure.to(fmt="cif")
        print(f"[INFO] 成功转换为CIF格式")
        
        # 构建保存路径，包含材料ID和空间群信息（清洗非法字符）
        safe_formula = _sanitize_filename(formula)
        safe_spacegroup = _sanitize_filename(spacegroup)
        filename = f"{safe_formula}_{material_id}_{safe_spacegroup}.cif"
        full_save_path = os.path.join(save_path, filename)
        print(f"[INFO] 保存路径: {full_save_path}")
        
        # 保存CIF文件
        with open(full_save_path, 'w', encoding='utf-8') as f:
            f.write(cif_data)
        
        # 验证文件
        if os.path.exists(full_save_path) and os.path.getsize(full_save_path) > 0:
            file_size = os.path.getsize(full_save_path)
            print(f"[SUCCESS] 材料 {material_id} CIF文件保存成功，大小: {file_size} 字节")
            return full_save_path
        else:
            print(f"[ERROR] 材料 {material_id} 文件保存失败或文件为空")
            return None
            
    except Exception as e:
        print(f"[ERROR] 下载材料 {material_id} 时失败: {e}")
        print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
        return None

async def download_with_criteria(formula, save_path, task_id, **kwargs):
    """
    便捷函数：使用常见筛选条件下载材料
    
    Args:
        formula: 化学式
        save_path: 保存路径  
        task_id: 任务ID
        **kwargs: 筛选条件，支持：
                 - spacegroup: 空间群符号
                 - max_energy_above_hull: 最大能量上凸包距离
                 - min_band_gap: 最小带隙
                 - max_band_gap: 最大带隙
                 - max_nsites: 最大原子数
                 - min_nsites: 最小原子数
                 - stable_only: 只选择稳定材料
                 - selection_mode: 选择模式
    
    Returns:
        下载结果
    """
    criteria = {}
    save_path = os.path.join(save_path, "mp_cif")
    # 构建筛选条件
    if 'spacegroup' in kwargs:
        criteria['spacegroup_symbol'] = kwargs['spacegroup']
    
    if 'max_energy_above_hull' in kwargs:
        criteria['energy_above_hull'] = (None, kwargs['max_energy_above_hull'])
    
    if 'min_band_gap' in kwargs or 'max_band_gap' in kwargs:
        min_bg = kwargs.get('min_band_gap', None)
        max_bg = kwargs.get('max_band_gap', None)
        criteria['band_gap'] = (min_bg, max_bg)
    
    if 'min_nsites' in kwargs or 'max_nsites' in kwargs:
        min_ns = kwargs.get('min_nsites', None)
        max_ns = kwargs.get('max_nsites', None)
        criteria['num_sites'] = (min_ns, max_ns)
    
    if kwargs.get('stable_only', False):
        criteria['is_stable'] = True
    
    selection_mode = kwargs.get('selection_mode', 'auto')
    
    return await download_cif_by_formula(formula, save_path, task_id, criteria, selection_mode)

async def batch_download_cifs(formulas, base_save_path, task_id, criteria=None, selection_mode="auto"):
    """
    批量下载多个化学式的CIF文件
    
    Args:
        formulas: 化学式列表
        base_save_path: 基础保存路径
        task_id: 任务ID
        criteria: 筛选条件
        selection_mode: 选择模式
    
    Returns:
        dict: 下载结果字典 {formula: save_path_or_None}
    """
    print(f"[INFO] 开始批量下载 {len(formulas)} 个化学式的CIF文件")
    
    results = {}
    
    for i, formula in enumerate(formulas, 1):
        print(f"\n[INFO] 处理第 {i}/{len(formulas)} 个化学式: {formula}")
        
        try:
            # 为每个化学式创建单独的保存目录
            formula_save_path = os.path.join(base_save_path, task_id, formula)
            
            result = await download_cif_by_formula(formula, formula_save_path, task_id, criteria, selection_mode)
            results[formula] = result
            
            if result:
                print(f"[SUCCESS] {formula} 下载成功")
            else:
                print(f"[FAILED] {formula} 下载失败")
                
            # 添加短暂延迟，避免API请求过于频繁
            await asyncio.sleep(0.5)
            
        except Exception as e:
            print(f"[ERROR] 处理化学式 {formula} 时发生错误: {e}")
            print(f"[ERROR] 详细错误信息: {traceback.format_exc()}")
            results[formula] = None
    
    # 统计结果
    successful = sum(1 for result in results.values() if result is not None)
    failed = len(formulas) - successful
    
    print(f"\n[SUMMARY] 批量下载完成:")
    print(f"[SUMMARY] 成功: {successful}, 失败: {failed}, 总计: {len(formulas)}")
    
    return results

async def download_with_retry(formula, save_path, task_id, max_retries=3, retry_delay=2, **kwargs):
    """
    带重试机制的下载函数
    
    Args:
        formula: 化学式
        save_path: 保存路径
        task_id: 任务ID
        max_retries: 最大重试次数
        retry_delay: 重试延迟（秒）
        **kwargs: 其他参数
    
    Returns:
        str|None: 成功时返回保存路径，失败时返回None
    """
    for attempt in range(max_retries + 1):
        if attempt > 0:
            print(f"[RETRY] 第 {attempt} 次重试下载 {formula}")
            await asyncio.sleep(retry_delay)
        
        try:
            result = await download_cif_by_formula(formula, save_path, task_id, **kwargs)
            if result:
                if attempt > 0:
                    print(f"[SUCCESS] {formula} 在第 {attempt} 次重试后下载成功")
                return result
        except Exception as e:
            print(f"[ERROR] 第 {attempt + 1} 次尝试失败: {e}")
            if attempt == max_retries:
                print(f"[FAILED] {formula} 在 {max_retries + 1} 次尝试后仍然失败")
    
    return None

# 示例使用和测试函数
async def main():
    """主函数示例"""
    print("=== Materials Project CIF 智能下载工具 (新版 mp-api) ===\n")
    
    # # 示例1: 基本下载
    # print("1. 基本下载测试 (自动选择最稳定):")
    # result = await download_cif_by_formula("LiFePO4", "downloads", "test_001")
    # if result:
    #     print(f"基本下载成功: {result}")
    # else:
    #     print("基本下载失败")
    
    # print("\n" + "="*60 + "\n")
    
    # 示例2: 指定空间群下载
    # print("2. 指定空间群下载:")
    # result2 = await download_with_criteria(
    #     "LiFePO4", 
    #     "downloads", 
    #     "test_002",
    #     spacegroup="Pnma",  # 指定空间群
    #     selection_mode="auto"
    # )
    # if result2:
    #     print(f"指定空间群下载成功: {result2}")
    
    # print("\n" + "="*60 + "\n")
    
    # 示例3: 只下载稳定材料
    # print("3. 只下载稳定材料:")
    # result3 = await download_with_criteria(
    #     "Li2O", 
    #     "downloads", 
    #     "test_003",
    #     stable_only=True,
    #     selection_mode="all"
    # )
    # if result3:
    #     print(f"稳定材料下载成功: {result3}")
    
    # print("\n" + "="*60 + "\n")
    
    # 示例4: 复合条件筛选
    print("4. 复合条件筛选:")
    result4 = await download_with_criteria(
        "TiO2",
        "downloads", 
        "test_004",
        max_energy_above_hull=0.1,  # 接近稳定
        min_band_gap=2.0,           # 宽带隙
        max_nsites=12,              # 小晶胞
        selection_mode="all"
    )
    if result4:
        print(f"复合条件筛选成功: {result4}")

if __name__ == "__main__":
    # 运行异步主函数
    asyncio.run(main())