import os

def get_vasp_config():
    """
    获取VASP计算的完整配置参数
    
    Returns:
        dict: 包含所有VASP配置的字典
    """
    
    # ========================
    # 基础路径和文件配置
    # ========================
    paths_config = {
        "POSCAR_FILE": "POSCAR",
        "VASP_PATH": "/data/app/vasp/6.3.2-intel/bin/vasp_std",
        "PSEUDO_PATH": '/data/home/ysl9527/software/psudopotential',
        "U_VALUES_JSON": "/data/home/ysl9527/software/u_values.json",
        "BADER_PATH": "/data/home/ysl9527/software/bader",
        "CHGSUM_PL_PATH": "/data/home/ysl9527/software/chgsum.pl",
    }
    
    # ========================
    # KPOINTS 参数
    # ========================
    kpoints_config = {
        "TARGET_KPOINT_PRODUCT": 30.0  # OPT 使用的目标乘积
    }
    
    # ========================
    # POTCAR 参数
    # ========================
    potcar_config = {
        "PAW_PICK_LIST": os.path.join(paths_config["PSEUDO_PATH"], 'PAWPickList'),
        "PAW_PATH": 'paw_pbe'
    }
    
    # ========================
    # INCAR 模板（集中管理）
    # ========================
    base_incars = {
"OXC": """
# 基础控制
SYSTEM = OPT
PREC = High
ENCUT = 520
ISTART = 0
ICHARG = 2
# 电子迭代
EDIFF = 1E-5
EDIFFG = -0.01
NELM = 100
NELMIN = 2
ALGO = Normal
ISMEAR = 0
SIGMA = 0.05
# 结构优化
IBRION = 2
NSW = 500
ISIF = 3
POTIM = 0.2
# 自旋极化与磁矩
ISPIN = 2
# 交换关联泛函
GGA = PE
# 并行与加速 
LREAL = Auto
# 输出控制
LWAVE  = .FALSE.
LCHARG = .FALSE.
""",
"ORC": """
# 基础控制
SYSTEM = OPT
PREC = High
ENCUT = 520
ISTART = 0
ICHARG = 2
# 电子迭代
EDIFF = 1E-5
EDIFFG = -0.01
NELM = 100
NELMIN = 2
ALGO = Normal
ISMEAR = 0
SIGMA = 0.05
# 结构优化
IBRION = 2
NSW = 500
ISIF = 3
POTIM = 0.2
# 自旋极化与磁矩
ISPIN = 2
# 交换关联泛函
GGA = PE
# 范德华力
IVDW = 12
VDW_RADIUS = 50.0
VDW_S8 = 0.7875
VDW_A1 = 0.4289
VDW_A2 = 4.4407
# 并行与加速 
LREAL = Auto
# 输出控制
LWAVE  = .FALSE.
LCHARG = .FALSE.
"""
    }
    
    # ========================
    # 分子动力学 INCAR 模板
    # ========================
    md_incar_content = """
# 基础控制
SYSTEM = AIMD             # 体系标识
PREC = Medium             # 中等精度模式
ENCUT = 400               # 截断能（最大的ENMIN）
# 电子迭代
EDIFF = 1E-4              # 电子步收敛阈值
LREAL = Auto              # 实空间投影
ISMEAR = 0                # Gaussian 展宽
SIGMA = 0.0407            # 展宽宽度，K*0.086*10-3
ALGO = Fast               # 快速算法
ISYM = 0                  # 关闭对称性（MD 必需）
# 分子动力学控制
IBRION = 0                # MD 模式（自由动力学）
NSW = 30000               # MD 总步数
POTIM = 0.5               # 时间步长（单位 fs）
TEBEG = 300               # 初始温度（单位 K）
TEEND = 300               # 最终温度（单位 K）
SMASS = 2                 # Nose-Hoover 热浴（SMASS=2 对应固定温度）
# 输出控制
LCHARG = .FALSE.          # 不输出 CHGCAR（节省存储）
LWAVE = .FALSE.           # 不输出 WAVECAR（若需续算则置为.TRUE.）
NWRITE = 1                # 输出详细程度（1=默认，2=详细）
# 收敛与性能
NELMIN = 4                # 最小电子步数（防止过早终止）
    """
    
    # ========================
    # 返回完整配置字典
    # ========================
    config = {
        "paths": paths_config,
        "kpoints": kpoints_config,
        "potcar": potcar_config,
        "base_incars": base_incars,
        "md_incar": md_incar_content.strip()
    }
    
    return config

def get_path_config():
    """
    获取路径相关配置
    
    Returns:
        dict: 路径配置字典
    """
    return get_vasp_config()["paths"]

def get_kpoints_config():
    """
    获取KPOINTS相关配置
    
    Returns:
        dict: KPOINTS配置字典
    """
    return get_vasp_config()["kpoints"]

def get_potcar_config():
    """
    获取POTCAR相关配置
    
    Returns:
        dict: POTCAR配置字典
    """
    return get_vasp_config()["potcar"]

def get_incar_templates():
    """
    获取INCAR模板
    
    Returns:
        dict: INCAR模板字典，包含不同类型的模板
    """
    return get_vasp_config()["base_incars"]

def get_md_incar_template():
    """
    获取分子动力学INCAR模板
    
    Returns:
        str: MD INCAR模板字符串
    """
    return get_vasp_config()["md_incar"]

def get_incar_template(template_type="OXC"):
    """
    获取指定类型的INCAR模板
    
    Args:
        template_type (str): 模板类型，支持 "OXC" 或 "ORC"
    
    Returns:
        str: INCAR模板字符串
    
    Raises:
        ValueError: 当模板类型不支持时抛出异常
    """
    templates = get_incar_templates()
    
    if template_type not in templates:
        available_types = list(templates.keys())
        raise ValueError(f"不支持的INCAR模板类型: {template_type}. 可用类型: {available_types}")
    
    return templates[template_type].strip()

# ========================
# 向后兼容的全局变量
# ========================
def _init_global_vars():
    """初始化全局变量以保持向后兼容性"""
    config = get_vasp_config()
    
    globals().update({
        # 路径配置
        "POSCAR_FILE": config["paths"]["POSCAR_FILE"],
        "VASP_PATH": config["paths"]["VASP_PATH"],
        "PSEUDO_PATH": config["paths"]["PSEUDO_PATH"],
        "U_VALUES_JSON": config["paths"]["U_VALUES_JSON"],
        
        # KPOINTS配置
        "TARGET_KPOINT_PRODUCT": config["kpoints"]["TARGET_KPOINT_PRODUCT"],
        
        # POTCAR配置
        "PAW_PICK_LIST": config["potcar"]["PAW_PICK_LIST"],
        "PAW_PATH": config["potcar"]["PAW_PATH"],
        
        # INCAR模板
        "BASE_INCARS": config["base_incars"],
        "MD_INCAR_CONTENT": config["md_incar"]
    })

# 初始化全局变量（保持向后兼容）
_init_global_vars()

# ========================
# 使用示例和测试函数
# ========================
def print_config_summary():
    """打印配置摘要"""
    config = get_vasp_config()
    
    print("=== VASP 配置摘要 ===")
    print(f"VASP路径: {config['paths']['VASP_PATH']}")
    print(f"POSCAR文件: {config['paths']['POSCAR_FILE']}")
    print(f"赝势路径: {config['paths']['PSEUDO_PATH']}")
    print(f"U值配置: {config['paths']['U_VALUES_JSON']}")
    print(f"K点乘积目标: {config['kpoints']['TARGET_KPOINT_PRODUCT']}")
    print(f"可用INCAR模板: {list(config['base_incars'].keys())}")
    print("=" * 50)

if __name__ == "__main__":
    # 测试配置函数
    print_config_summary()
    
    # 测试获取特定模板
    try:
        oxc_template = get_incar_template("OXC")
        print("\n=== OXC 模板预览 ===")
        print(oxc_template[:200] + "...")
        
        orc_template = get_incar_template("ORC")
        print("\n=== ORC 模板预览 ===")
        print(orc_template[:200] + "...")
        
        md_template = get_md_incar_template()
        print("\n=== MD 模板预览 ===")
        print(md_template[:200] + "...")
        
    except Exception as e:
        print(f"错误: {e}")