import os
import logging
from pathlib import Path

from .settings import settings

logger = logging.getLogger(__name__)

BASE_URL = settings.vasp_server_base_url
DOWNLOAD_URL = settings.vasp_work_root
MCP_PORT = settings.mcp_port
VASP_remote_run_url = settings.vasp_remote_run_url
VASP_remote_run_port = settings.vasp_remote_run_port
static_url = settings.vasp_public_base_url
OSS_ENDPOINT = settings.s3_endpoint or settings.oss_endpoint
OSS_BUCKET = settings.oss_bucket
OSS_REGION = settings.oss_region
OSS_ACCESS_KEY_ID = settings.oss_access_key_id
OSS_ACCESS_KEY_SECRET = settings.oss_access_key_secret
ARTIFACT_URL_EXPIRE_SECONDS = settings.artifact_url_expire_seconds


def get_download_url(path: str) -> str:
    p = Path(path)
    try:
        rel = p.relative_to(DOWNLOAD_URL).as_posix()
    except ValueError:
        rel = p.as_posix().lstrip("/")
    return f"{static_url}/download/file/{rel}"


def get_static_url(path: str) -> str:
    p = Path(path)
    try:
        rel = p.relative_to(DOWNLOAD_URL).as_posix()
    except ValueError:
        rel = p.as_posix().lstrip("/")
    return f"{static_url}/static/{rel}"






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
        "VASP_PATH": settings.vasp_path,
        "PSEUDO_PATH": settings.pseudo_path,
        "U_VALUES_JSON": settings.u_values_json,
        "BADER_PATH": settings.bader_path,
        "CHGSUM_PL_PATH": settings.chgsum_pl_path,
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
"DEFAULT": """
# 基础控制
SYSTEM = OPT
PREC = Normal
ENCUT = 450
ISTART = 0
ICHARG = 2
# 电子迭代
EDIFF = 1E-4
EDIFFG = -0.02
NELM = 100
NELMIN = 2
ALGO = Normal
ISMEAR = 0
SIGMA = 0.05
# 结构优化
IBRION = 2
NSW = 500
ISIF = 2
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

def get_incar_template():
    """
    获取默认INCAR模板

    Returns:
        str: INCAR模板字符串
    """
    return get_incar_templates()["DEFAULT"].strip()

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

    logger.info("=== VASP 配置摘要 ===")
    logger.info("VASP路径: %s", config['paths']['VASP_PATH'])
    logger.info("POSCAR文件: %s", config['paths']['POSCAR_FILE'])
    logger.info("赝势路径: %s", config['paths']['PSEUDO_PATH'])
    logger.info("U值配置: %s", config['paths']['U_VALUES_JSON'])
    logger.info("K点乘积目标: %s", config['kpoints']['TARGET_KPOINT_PRODUCT'])
    logger.info("可用INCAR模板: %s", list(config['base_incars'].keys()))
    logger.info("=" * 50)

if __name__ == "__main__":
    print_config_summary()

    try:
        default_template = get_incar_template()
        print("\n=== DEFAULT 模板预览 ===")
        print(default_template[:200] + "...")

        md_template = get_md_incar_template()
        print("\n=== MD 模板预览 ===")
        print(md_template[:200] + "...")
    except Exception as e:
        print(f"错误: {e}")
