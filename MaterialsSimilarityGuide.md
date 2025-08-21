# Materials Project 相似度查询指南

## 概述

是的！Materials Project确实提供了强大的相似度查询功能，非常适合新材料设计中寻找参考结构。本指南将介绍几种主要的相似度查询方法。

## 主要相似度查询方法

### 1. **结构相似度查询**

Materials Project使用基于局部配位环境的算法来计算结构相似度：

- **原理**: 基于61维结构指纹，包含配位数、几何形状等特征
- **相似度计算**: 使用欧几里德距离，距离越小相似度越高
- **阈值**: 距离 < 0.9 通常认为是相似结构

### 2. **通过API进行相似度查询的几种方式**

#### 方式1: 新版API的similarity端点
```python
from mp_api.client import MPRester

with MPRester("your_api_key") as mpr:
    # 查询相似材料（如果可用）
    similarity_docs = mpr.materials.similarity.search(
        material_ids=["mp-66"],  # 目标材料ID
        fields=["material_id", "similar_structures"]
    )
```

#### 方式2: 使用matminer进行结构指纹比较
```python
import numpy as np
from mp_api.client import MPRester
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint

with MPRester() as mpr:
    # 获取结构
    diamond = mpr.get_structure_by_material_id("mp-66")
    gaas = mpr.get_structure_by_material_id("mp-2534")

# 计算结构指纹
ssf = SiteStatsFingerprint(
    CrystalNNFingerprint.from_preset('ops', distance_cutoffs=None, x_diff_weight=0),
    stats=('mean', 'std_dev', 'minimum', 'maximum')
)

v_diamond = np.array(ssf.featurize(diamond))
v_gaas = np.array(ssf.featurize(gaas))

# 计算相似度
distance = np.linalg.norm(v_diamond - v_gaas)
similarity = np.exp(-distance) * 100  # 转换为百分比
```

#### 方式3: 基于空间群和化学组成的简单相似度
```python
from pymatgen.ext.matproj import MPRester

def find_similar_by_composition_and_spacegroup(target_mp_id, api_key):
    """
    基于化学组成和空间群查找相似材料
    """
    with MPRester(api_key) as mpr:
        # 获取目标材料信息
        target_data = mpr.get_data(target_mp_id, prop="elements,spacegroup,pretty_formula")
        target_info = target_data[0]
        
        # 构建查询条件
        chemsys = "-".join(target_info['elements'])
        target_spacegroup = target_info['spacegroup']['symbol']
        
        # 查询相同化学系统的材料
        candidates = mpr.get_data(chemsys, prop="material_id,pretty_formula,spacegroup,e_above_hull")
        
        # 筛选相同空间群且稳定的材料
        similar_materials = []
        for material in candidates:
            if (material['spacegroup']['symbol'] == target_spacegroup and 
                material.get('e_above_hull', 1.0) < 0.1):  # 稳定性过滤
                similar_materials.append(material)
        
        return similar_materials
```

### 3. **实际应用示例**

#### 示例1: 寻找钙钛矿类型结构
```python
# 以CaTiO3为原型
prototype_id = "mp-5827"
similar_perovskites = find_similar_by_composition_and_spacegroup(prototype_id, api_key)
```

#### 示例2: 查找特定化学系统中的稳定材料
```python
def find_stable_materials_in_system(elements, api_key):
    """
    在指定化学系统中查找稳定材料
    """
    chemsys = "-".join(elements)
    with MPRester(api_key) as mpr:
        materials = mpr.get_data(
            chemsys, 
            prop="material_id,pretty_formula,e_above_hull,formation_energy_per_atom"
        )
        
        # 筛选稳定材料 (e_above_hull < 0.05 eV)
        stable_materials = [m for m in materials if m.get('e_above_hull', 1.0) < 0.05]
        
        # 按稳定性排序
        stable_materials.sort(key=lambda x: x['e_above_hull'])
        
        return stable_materials

# 示例使用
li_fe_o_stable = find_stable_materials_in_system(['Li', 'Fe', 'O'], api_key)
```

### 4. **新材料设计的实用策略**

#### 策略1: 原型替换法
1. 找到一个已知稳定的原型结构（如钙钛矿、尖晶石等）
2. 查询该原型的空间群信息
3. 搜索具有相同空间群但不同元素组成的材料
4. 评估这些材料的稳定性和性质

#### 策略2: 化学系统扩展法
1. 从已知的稳定化合物开始
2. 逐步替换或添加元素
3. 在新的化学系统中寻找稳定的化合物
4. 比较结构特征和性质

#### 策略3: 性质导向搜索
```python
def find_materials_by_properties(band_gap_range, density_range, api_key):
    """
    基于性质范围查找材料
    """
    with MPRester(api_key) as mpr:
        # 使用query方法进行复杂查询
        materials = mpr.query(
            criteria={
                'band_gap': {'$gte': band_gap_range[0], '$lte': band_gap_range[1]},
                'density': {'$gte': density_range[0], '$lte': density_range[1]},
                'e_above_hull': {'$lte': 0.1}
            },
            properties=['material_id', 'pretty_formula', 'band_gap', 'density', 'spacegroup']
        )
        
        return materials
```

## 5. **相似度指标解释**

### 结构距离指标
- **距离 = 0**: 完全相同的结构（如金刚石 vs GaAs）
- **距离 < 0.9**: 高度相似的结构
- **距离 > 3.0**: 结构差异很大

### 实际案例
- 金刚石 vs GaAs: 距离 = 0 (高度相似)
- 金刚石 vs 岩盐结构: 距离 = 3.57 (差异很大)
- 钙钛矿 vs 岩盐: 距离 = 2.74 (中等差异)

## 6. **推荐的工作流程**

### 新材料设计流程
1. **确定设计目标**: 明确需要的性质（带隙、密度等）
2. **寻找原型**: 在已知材料中找到具有相似性质的化合物
3. **结构分析**: 分析原型的结构特征和空间群
4. **相似度搜索**: 使用上述方法查找结构相似的材料
5. **稳定性评估**: 检查候选材料的热力学稳定性
6. **性质预测**: 基于相似材料预测新材料的性质

### 具体实施步骤
```python
# 1. 获取原型信息
prototype_id = "mp-5827"  # CaTiO3钙钛矿
with MPRester(api_key) as mpr:
    prototype_structure = mpr.get_structure_by_material_id(prototype_id)
    prototype_info = mpr.get_data(prototype_id, prop="spacegroup,elements")

# 2. 设计新的元素组合
new_elements = ['Ba', 'Ti', 'O']  # 尝试Ba替换Ca

# 3. 查找新化学系统中的材料
new_chemsys = "-".join(new_elements)
candidates = mpr.get_data(new_chemsys, prop="material_id,pretty_formula,spacegroup,e_above_hull")

# 4. 筛选具有相同结构的材料
similar_structure_candidates = [
    mat for mat in candidates 
    if mat['spacegroup']['symbol'] == prototype_info[0]['spacegroup']['symbol']
]
```

## 7. **注意事项**

1. **API限制**: 避免过于频繁的查询，注意API使用限制
2. **数据质量**: 优先选择计算收敛且稳定的材料
3. **实验验证**: 理论预测需要实验验证
4. **版本更新**: Materials Project数据库会定期更新，注意版本差异

## 8. **有用的材料ID参考**

### 常见原型结构
- **金刚石**: mp-66
- **GaAs (闪锌矿)**: mp-2534  
- **CaTiO3 (钙钛矿)**: mp-5827
- **NaCl (岩盐)**: mp-22862
- **尖晶石**: mp-1408976

这些方法可以帮助您在新材料设计中快速找到结构相似的参考材料，为进一步的理论计算和实验验证提供指导。 