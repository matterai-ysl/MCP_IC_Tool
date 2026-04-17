import os
import io
import subprocess
import sys
import math
import json
import logging
from .Config import get_kpoints_config,get_vasp_config,get_incar_template

logger = logging.getLogger(__name__)

_DEFAULT_MAGMOM_BY_ELEMENT = {
    "TI": 1.0,
    "V": 3.0,
    "CR": 5.0,
    "MN": 3.0,
    "FE": 5.0,
    "CO": 3.0,
    "NI": 2.0,
    "CU": 1.0,
    "MO": 1.0,
    "TC": 2.0,
    "RU": 2.0,
    "RH": 1.0,
    "W": 1.0,
    "RE": 2.0,
    "OS": 2.0,
    "IR": 1.0,
}

# ========================
# 工具函数
# ========================


def _normalize_cif_text_for_pymatgen(cif_text: str) -> str:
    """补齐常见但非严格规范的 CIF 字段，避免 pymatgen 直接报错。"""
    lines = cif_text.splitlines()
    normalized: list[str] = []
    index = 0

    while index < len(lines):
        line = lines[index]
        if line.strip() != "loop_":
            normalized.append(line)
            index += 1
            continue

        header_start = index
        header_end = index + 1
        while header_end < len(lines) and lines[header_end].lstrip().startswith("_"):
            header_end += 1

        headers = lines[index + 1:header_end]
        is_atom_site_loop = any("_atom_site_type_symbol" in header for header in headers) and any(
            "_atom_site_fract_" in header for header in headers
        )
        has_atom_site_label = any("_atom_site_label" in header for header in headers)

        if not is_atom_site_loop or has_atom_site_label:
            normalized.extend(lines[header_start:header_end])
            index = header_end
            continue

        normalized.append("loop_")
        inserted_label = False
        for header in headers:
            if "_atom_site_type_symbol" in header and not inserted_label:
                indent = header[: len(header) - len(header.lstrip())]
                normalized.append(f"{indent}_atom_site_label")
                inserted_label = True
            normalized.append(header)

        row_index = 1
        data_index = header_end
        while data_index < len(lines):
            row = lines[data_index]
            stripped = row.strip()
            if not stripped:
                normalized.append(row)
                data_index += 1
                break
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            if stripped.startswith("#"):
                normalized.append(row)
                data_index += 1
                continue

            tokens = stripped.split()
            species = tokens[0].strip("'\"") if tokens else f"site{row_index}"
            leading = row[: len(row) - len(row.lstrip())]
            normalized.append(f"{leading}{species}{row_index} {stripped}")
            row_index += 1
            data_index += 1

        index = data_index

    return "\n".join(normalized) + ("\n" if cif_text.endswith("\n") else "")


def _sanitize_structure_for_vasp(structure):
    """清理数值噪声，避免 POSCAR 中出现接近 0 的伪分量。"""
    try:
        from pymatgen.core import Structure
    except ImportError as exc:
        raise RuntimeError("需要安装 pymatgen 以处理结构文件") from exc

    zero_tol = 1e-12
    lattice_matrix = []
    for vector in structure.lattice.matrix:
        cleaned_vector = []
        for value in vector:
            cleaned_vector.append(0.0 if abs(float(value)) < zero_tol else float(value))
        lattice_matrix.append(cleaned_vector)

    frac_coords = []
    for coords in structure.frac_coords:
        cleaned_coords = []
        for value in coords:
            numeric = float(value)
            if abs(numeric) < zero_tol or abs(numeric - 1.0) < zero_tol:
                numeric = 0.0
            cleaned_coords.append(numeric % 1.0)
        frac_coords.append(cleaned_coords)

    return Structure(
        lattice_matrix,
        [site.species for site in structure],
        frac_coords,
        coords_are_cartesian=False,
        site_properties=structure.site_properties,
        to_unit_cell=True,
    )


def run_command_background(command, work_dir=None, log_file="result.log"):
    """后台运行命令，立即返回"""
    full_log_path = os.path.join(work_dir, log_file) if work_dir else log_file
    
    try:
        logger.debug("后台执行命令: %s", ' '.join(command))
        
        with open(full_log_path, 'w') as f:
            # 使用 Popen 而不是 run，不等待完成
            process = subprocess.Popen(
                command,
                cwd=work_dir,
                stdout=f,
                stderr=subprocess.STDOUT
            )
            
            return {
                "pid": process.pid,
                "status": "started",
                "log_file": full_log_path
            }
            
    except Exception as e:
        return {"error": str(e)}
def load_u_values(config_path):
    """从JSON配置文件加载U值"""
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        raise SystemExit(f"错误：配置文件 {config_path} 未找到")
    except json.JSONDecodeError:
        raise SystemExit(f"错误：配置文件 {config_path} 格式不正确")
def cif_to_poscar(cif_path, work_dir='./'):
    """将CIF转换为POSCAR（vasp5格式）"""
    logger.info("正在转换CIF到POSCAR...")
    if not os.path.exists(cif_path):
        sys.exit("错误：未找到文件 {}".format(cif_path))

    poscar_path = os.path.join(work_dir, "POSCAR")
    try:
        from pymatgen.io.cif import CifParser
        from pymatgen.io.vasp import Poscar
    except ImportError as exc:
        raise RuntimeError("需要安装 pymatgen 以将 CIF 转换为 POSCAR") from exc

    with open(cif_path, "r", encoding="utf-8", errors="ignore") as fh:
        cif_text = fh.read()

    parser = CifParser(io.StringIO(_normalize_cif_text_for_pymatgen(cif_text)))
    structures = parser.parse_structures(primitive=False)
    if not structures:
        raise RuntimeError(f"无法从 CIF 解析结构: {cif_path}")

    sanitized_structure = _sanitize_structure_for_vasp(structures[0])
    Poscar(sanitized_structure).write_file(poscar_path)
    logger.info("POSCAR生成成功")
    return poscar_path

def read_poscar_elements(poscar_path="POSCAR"):
    """从POSCAR读取元素顺序"""
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    elements = lines[5].strip().split()
    unique_elements = []
    for elem in elements:
        if elem not in unique_elements:
            unique_elements.append(elem)
    return unique_elements


def read_poscar_species_counts(poscar_path="POSCAR"):
    """从 POSCAR 读取元素及其个数"""
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    elements = lines[5].strip().split()
    counts = [int(value) for value in lines[6].strip().split()]
    return list(zip(elements, counts))


def generate_ldau_params(elements, u_values):
    """生成LDAU参数"""
    ldaul = []
    ldauu = []
    ldauj = []
    for elem in elements:
        if elem not in u_values:
            logger.warning("元素 %s 未定义 U 值，使用 U=0.0", elem)
            u_values[elem] = 0.0
        if u_values[elem] > 0:
            ldaul.append(2)
        else:
            ldaul.append(-1)
        ldauu.append(u_values[elem])
        ldauj.append(0.0)
    return ldaul, ldauu, ldauj


def generate_magmom_params(species_counts):
    """根据元素种类生成一个保守的 MAGMOM 初始猜测"""
    return " ".join(
        f"{count}*{_DEFAULT_MAGMOM_BY_ELEMENT.get(element.upper(), 0.0):.1f}"
        for element, count in species_counts
    )


def generate_kpoints(work_dir='./', target_product=None):

    """生成 KPOINTS 文件"""
    poscar_path = os.path.join(work_dir, "POSCAR")
    target_product = float(target_product if target_product is not None else get_kpoints_config()["TARGET_KPOINT_PRODUCT"])
    def get_real_lattice_lengths(path):
        """从 POSCAR 获取实空间晶格常数"""
        try:
            with open(path, 'r') as f:
                lines = [f.readline() for _ in range(5)][2:5]
            return [
                math.sqrt(sum(float(x) ** 2 for x in line.strip().split()))
                for line in lines if line.strip()
            ]
        except IOError:
            raise RuntimeError(f"无法读取文件: {path}")
        except ValueError:
            raise RuntimeError("POSCAR 格式不正确")

    lattice_lengths = get_real_lattice_lengths(poscar_path)
    kmesh = [max(int(round(target_product / l)), 1) for l in lattice_lengths]
    kpoints_path = os.path.join(work_dir, 'KPOINTS')
    try:
        with open(kpoints_path, 'w') as f:
            f.write("Automatically generated K-points\n")
            f.write("0\nGamma\n")
            f.write("%3d %3d %3d\n" % tuple(kmesh))
            f.write("0 0 0\n")
        logger.info("KPOINTS 文件已生成于 %s", kpoints_path)
        return True
    except IOError:
        raise RuntimeError("无法写入 KPOINTS 文件")

def generate_potcar(work_dir='./'):
    """生成 POTCAR 文件"""
    poscar_path = os.path.join(work_dir, "POSCAR")
    try:
        with open(poscar_path, 'r') as f:
            lines = f.readlines()
        elements = lines[5].strip().split()
    except Exception as e:
        raise RuntimeError(f"读取 POSCAR 失败: {str(e)}")
    potcar_paths = []
    vasp_config = get_vasp_config()["paths"]
    paw_pick_list = os.path.join(vasp_config["PSEUDO_PATH"], "PAWPickList")
    paw_path = 'paw_pbe'
    pseudo_path = vasp_config["PSEUDO_PATH"]
    for element in elements:
        try:
            cmd = ["awk", f'$1 == "{element}" {{print $2}}', paw_pick_list]
            output = subprocess.check_output(cmd, universal_newlines=True).strip()
            if not output:
                raise RuntimeError(f"未在 {paw_pick_list} 中找到元素 {element}")
            potcar = os.path.join(pseudo_path, paw_path, output, 'POTCAR')
            if not os.path.exists(potcar):
                raise RuntimeError(f"POTCAR 文件不存在: {potcar}")
            potcar_paths.append(potcar)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"查找元素 {element} 时发生错误: {e}")
    potcar_path = os.path.join(work_dir, 'POTCAR')
    try:
        with open(potcar_path, 'wb') as outfile:
            for potcar in potcar_paths:
                with open(potcar, 'rb') as infile:
                    outfile.write(infile.read())
        logger.info("POTCAR 文件已生成于 %s", potcar_path)
        return True
    except IOError as e:
        raise RuntimeError(f"写入 POTCAR 文件失败: {str(e)}")

def generate_incar(work_dir):
    """生成默认 INCAR 文件"""
    incar_path = os.path.join(work_dir, "INCAR")
    with open(incar_path, 'w') as f:
        f.write(get_incar_template())
    u_values = load_u_values(get_vasp_config()["paths"]["U_VALUES_JSON"])
    poscar_path = os.path.join(work_dir, "POSCAR")
    species_counts = read_poscar_species_counts(poscar_path)
    elements = [element for element, _ in species_counts]
    ldaul, ldauu, ldauj = generate_ldau_params(elements, u_values)
    magmom = generate_magmom_params(species_counts)
    with open(incar_path, 'a') as f:
        f.write("\n# LDA+U参数（自动添加）\n")
        f.write("LDAU = .TRUE.\n")
        f.write("LDAUTYPE = 2\n")
        f.write("LMAXMIX = 4\n")
        f.write("LDAUL = {}\n".format(' '.join(map(str, ldaul))))
        f.write("LDAUU = {}\n".format(' '.join(map(str, ldauu))))
        f.write("LDAUJ = {}\n".format(' '.join(map(str, ldauj))))
        f.write("\n# 初始磁矩（自动添加）\n")
        f.write(f"MAGMOM = {magmom}\n")
    logger.info("INCAR 文件已生成于 %s", incar_path)
    return True

if __name__ == "__main__":
    cif_to_poscar("/Users/ysl/Desktop/Code/MCP_IC_Tool/downloads/Li2O_mp-1960_Fm-3m.cif")
