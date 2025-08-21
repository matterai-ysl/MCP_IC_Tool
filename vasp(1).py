#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import math
import json

# ========================
# 配置参数（可根据需要修改）
# ========================
POSCAR_FILE = "POSCAR"
VASP_PATH = "/gpfs/home/lixingyu/software/vasp.6.3.2/vasp.6.3.2/bin/vasp_std"
# KPOINTS 参数
TARGET_KPOINT_PRODUCT = 30.0  # OPT 使用的目标乘积
# POTCAR 参数
PSEUDO_PATH = '/gpfs/home/lixingyu/bin/psudopotential'
PAW_PICK_LIST = os.path.join(PSEUDO_PATH, 'PAWPickList')
PAW_PATH = 'paw_pbe'
# U值配置文件路径
U_VALUES_JSON = "/gpfs/home/lixingyu/bin/u_values.json"

# ========================
# INCAR 模板（集中管理）
# ========================
BASE_INCARS = {
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

# 固定 MD 的 INCAR 内容
MD_INCAR_CONTENT = """
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
# 工具函数
# ========================
def run_command(command, cwd=None, log_file="result.log"):
    """
    运行系统命令并记录输出到日志文件

    参数:
        command: 要执行的命令列表
        cwd: 工作目录
        log_file: 日志文件名（相对于cwd）

    返回:
        bool: 成功返回 True，失败返回 False
    """
    full_log_path = os.path.join(cwd, log_file) if cwd else log_file

    try:
        print(f"执行命令: {' '.join(command)}")
        with open(full_log_path, 'w') as f:
            # 将标准输出和标准错误都重定向到日志文件
            result = subprocess.run(
                command,
                cwd=cwd,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )

        if result.returncode != 0:
            print(f"命令执行失败，错误日志已保存到: {full_log_path}")
            return False
        return True
    except Exception as e:
        print(f"运行命令时发生错误: {str(e)}")
        return False


def check_convergence(outcar_path):
    """检查 OUTCAR 文件最后一行是否包含 'Voluntary'"""
    try:
        with open(outcar_path, 'rb') as f:
            f.seek(-1024, os.SEEK_END)
            last_line = f.readlines()[-1].decode()
            return 'Voluntary' in last_line
    except Exception as e:
        print(f"无法读取 OUTCAR: {str(e)}")
        return False


# ========================
# KPOINTS 生成函数
# ========================
def generate_kpoints(caldir='./', target_product=30.0):
    """生成 KPOINTS 文件"""
    poscar_path = os.path.join(caldir, POSCAR_FILE)

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
    kpoints_path = os.path.join(caldir, 'KPOINTS')
    try:
        with open(kpoints_path, 'w') as f:
            f.write("Automatically generated K-points\n")
            f.write("0\nGamma\n")
            f.write("%3d %3d %3d\n" % tuple(kmesh))
            f.write("0 0 0\n")
        print(f"KPOINTS 文件已生成于 {kpoints_path}")
        return True
    except IOError:
        raise RuntimeError("无法写入 KPOINTS 文件")


def generate_kpoints_from_opt(opt_dir, target_dir):
    """
    从 opt/KPOINTS 读取并生成 DOS 用的 KPOINTS（每个方向翻倍）
    """
    opt_kpoints_path = os.path.join(opt_dir, "KPOINTS")
    target_kpoints_path = os.path.join(target_dir, "KPOINTS")
    try:
        with open(opt_kpoints_path, 'r') as f:
            lines = f.readlines()
        kmesh_line = None
        for line in lines:
            parts = line.strip().split()
            if len(parts) == 3 and all(p.isdigit() for p in parts):
                kmesh_line = line.strip()
                break
        if not kmesh_line:
            raise RuntimeError("无法从 opt/KPOINTS 读取有效的 k 点数")
        kmesh = list(map(int, kmesh_line.split()))
        kmesh_dos = [k * 2 for k in kmesh]
        with open(target_kpoints_path, 'w') as f:
            f.write("Automatically generated K-points (DOS)\n")
            f.write("0\nGamma\n")
            f.write("%3d %3d %3d\n" % tuple(kmesh_dos))
            f.write("0 0 0\n")
        print(f"DOS KPOINTS 已生成于 {target_kpoints_path}")
        return True
    except Exception as e:
        print(f"生成 DOS KPOINTS 失败: {str(e)}")
        return False


def generate_kpoints_fixed_md(target_dir):
    """生成固定为 1 1 1 的 KPOINTS（用于 MD）"""
    kpoints_path = os.path.join(target_dir, "KPOINTS")
    try:
        with open(kpoints_path, 'w') as f:
            f.write("Fixed KPOINTS for MD\n")
            f.write("0\nGamma\n")
            f.write(" 1  1  1\n")
            f.write("0 0 0\n")
        print(f"MD KPOINTS 已生成于 {kpoints_path}")
        return True
    except IOError:
        raise RuntimeError("无法写入 MD KPOINTS 文件")


# ========================
# POTCAR 生成函数
# ========================
def generate_potcar(caldir='./'):
    """生成 POTCAR 文件"""
    poscar_path = os.path.join(caldir, POSCAR_FILE)
    try:
        with open(poscar_path, 'r') as f:
            lines = f.readlines()
        elements = lines[5].strip().split()
    except Exception as e:
        raise RuntimeError(f"读取 POSCAR 失败: {str(e)}")
    potcar_paths = []
    for element in elements:
        try:
            cmd = ["awk", f'$1 == "{element}" {{print $2}}', PAW_PICK_LIST]
            output = subprocess.check_output(cmd, universal_newlines=True).strip()
            if not output:
                raise RuntimeError(f"未在 {PAW_PICK_LIST} 中找到元素 {element}")
            potcar = os.path.join(PSEUDO_PATH, PAW_PATH, output, 'POTCAR')
            if not os.path.exists(potcar):
                raise RuntimeError(f"POTCAR 文件不存在: {potcar}")
            potcar_paths.append(potcar)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"查找元素 {element} 时发生错误: {e}")
    potcar_path = os.path.join(caldir, 'POTCAR')
    try:
        with open(potcar_path, 'wb') as outfile:
            for potcar in potcar_paths:
                with open(potcar, 'rb') as infile:
                    outfile.write(infile.read())
        print(f"POTCAR 文件已生成于 {potcar_path}")
        return True
    except IOError as e:
        raise RuntimeError(f"写入 POTCAR 文件失败: {str(e)}")


# ========================
# 体系类型识别函数
# ========================
def read_calc_type(poscar_path="POSCAR"):
    """
    读取 POSCAR 第一行，确定计算类型（OXC/SSE 或 ORC/ECAT_OER/ECAT_HER）
    支持的体系类型：
        - OXC         → 使用 OXC 模板
        - SSE         → 使用 OXC 模板
        - ORC         → 使用 ORC 模板
        - ECAT_OER    → 使用 ORC 模板
        - ECAT_HER    → 使用 ORC 模板
    返回:
        str: "OXC" 或 "ORC"，表示使用的 INCAR 模板类型
    """
    try:
        with open(poscar_path, 'r') as f:
            first_line = f.readline().strip().upper()
        if not first_line:
            raise ValueError("POSCAR 第一行为空")
        first_word = first_line.split()[0]
        oxc_types = ["OXC", "SSE"]
        orc_types = ["ORC", "ECAT_OER", "ECAT_HER"]
        if first_word in oxc_types:
            return "OXC"
        elif first_word in orc_types:
            return "ORC"
        else:
            raise ValueError(f"""
未识别的体系类型: {first_word}
请确保 POSCAR 第一行为以下关键词之一：
- OXC 或 SSE → 使用 OXC 模板
- ORC、ECAT_OER 或 ECAT_HER → 使用 ORC 模板
""")
    except FileNotFoundError:
        raise SystemExit("错误：未找到 POSCAR 文件，请确认当前目录下存在 POSCAR")
    except Exception as e:
        raise SystemExit(f"错误读取 POSCAR: {str(e)}")


# ========================
# U值自动设置函数
# ========================
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


def load_u_values(config_path):
    """从JSON配置文件加载U值"""
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        raise SystemExit(f"错误：配置文件 {config_path} 未找到")
    except json.JSONDecodeError:
        raise SystemExit(f"错误：配置文件 {config_path} 格式不正确")


def generate_ldau_params(elements, u_values):
    """生成LDAU参数"""
    ldaul = []
    ldauu = []
    ldauj = []
    for elem in elements:
        if elem not in u_values:
            print(f"警告：元素 {elem} 未定义 U 值，使用 U=0.0")
            u_values[elem] = 0.0
        if u_values[elem] > 0:
            ldaul.append(2)
        else:
            ldaul.append(-1)
        ldauu.append(u_values[elem])
        ldauj.append(0.0)
    return ldaul, ldauu, ldauj


def generate_incar(caldir, calc_type):
    """根据 calc_type 生成 INCAR 文件"""
    incar_path = os.path.join(caldir, "INCAR")
    with open(incar_path, 'w') as f:
        f.write(BASE_INCARS[calc_type].strip())
    u_values = load_u_values(U_VALUES_JSON)
    elements = read_poscar_elements(os.path.join(caldir, POSCAR_FILE))
    ldaul, ldauu, ldauj = generate_ldau_params(elements, u_values)
    with open(incar_path, 'a') as f:
        f.write("\n# LDA+U参数（自动添加）\n")
        f.write("LDAU = .TRUE.\n")
        f.write("LDAUTYPE = 2\n")
        f.write("LMAXMIX = 4\n")
        f.write("LDAUL = {}\n".format(' '.join(map(str, ldaul))))
        f.write("LDAUU = {}\n".format(' '.join(map(str, ldauu))))
        f.write("LDAUJ = {}\n".format(' '.join(map(str, ldauj))))
    print(f"INCAR 文件已生成于 {incar_path}")
    return True


# ========================
# INCAR 修改函数
# ========================
def get_total_atoms(poscar_path):
    """从 POSCAR 获取总原子数"""
    try:
        with open(poscar_path, 'r') as f:
            lines = f.readlines()
        atom_counts = list(map(int, lines[6].strip().split()))
        return sum(atom_counts)
    except Exception as e:
        raise RuntimeError(f"读取原子数失败: {str(e)}")


def modify_incar_for_scf(opt_incar_path, scf_incar_path, poscar_path):
    """基于 opt/INCAR 修改生成 scf/INCAR（只改不加）"""
    try:
        with open(opt_incar_path, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            stripped = line.strip().upper()
            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = SCF\n")
            elif stripped.startswith("PREC"):
                new_lines.append("PREC = Accurate\n")
            elif stripped.startswith("NSW"):
                new_lines.append("NSW = 0\n")
            elif stripped.startswith("ISIF"):
                new_lines.append("ISIF = 2\n")
            elif stripped.startswith("POTIM"):
                new_lines.append("POTIM = 0\n")
            elif stripped.startswith("EDIFFG"):
                new_lines.append("#EDIFFG = -0.01\n")
            elif stripped.startswith("LWAVE"):
                new_lines.append("#LWAVE = .FALSE.\n")
            elif stripped.startswith("LCHARG"):
                new_lines.append("#LCHARG = .FALSE.\n")
            else:
                new_lines.append(line)
        with open(scf_incar_path, 'w') as f:
            f.writelines(new_lines)
        print(f"SCF INCAR 已生成于 {scf_incar_path}")
        return True
    except Exception as e:
        print(f"修改 INCAR 失败: {str(e)}")
        return False


def modify_incar_for_dos(scf_incar_path, dos_incar_path, poscar_path):
    """基于 scf/INCAR 修改生成 dos/INCAR"""
    try:
        with open(scf_incar_path, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            stripped = line.strip().upper()
            if stripped.startswith("SYSTEM"):
                new_lines.append("SYSTEM = DOS\n")
            elif stripped.startswith("ISTART"):
                new_lines.append("ISTART = 1\n")
            elif stripped.startswith("ICHARG"):
                new_lines.append("ICHARG = 11\n")
            else:
                new_lines.append(line)
        # 添加 DOS 专用参数
        new_lines.append("NEDOS = 5000\n")
        new_lines.append("LORBIT = 11\n")
        with open(dos_incar_path, 'w') as f:
            f.writelines(new_lines)
        print(f"DOS INCAR 已生成于 {dos_incar_path}")
        return True
    except Exception as e:
        print(f"修改 INCAR 失败: {str(e)}")
        return False


# ========================
# 步骤函数
# ========================
def opt_step(calc_type):
    """结构优化步骤（带收敛检测和重试机制）"""
    step_dir = "1-opt"
    os.makedirs(step_dir, exist_ok=True)

    # 初始文件准备
    shutil.copy(POSCAR_FILE, step_dir)

    max_retries = 5
    retry_count = 0
    converged = False

    while retry_count < max_retries:
        # 清理旧数据（除POSCAR外）
        for f in ['CHGCAR', 'WAVECAR', 'CONTCAR', 'OSZICAR', 'OUTCAR', 'INCAR', 'KPOINTS', 'POTCAR']:
            path = os.path.join(step_dir, f)
            if os.path.exists(path):
                if f == 'POSCAR':
                    continue
                os.remove(path)

        # 生成 KPOINTS 和 POTCAR
        if not generate_kpoints(step_dir, TARGET_KPOINT_PRODUCT):
            return False
        if not generate_potcar(step_dir):
            return False

        # 生成 INCAR
        if not generate_incar(step_dir, calc_type):
            return False

        # 运行 VASP
        print(f"运行 VASP（尝试次数: {retry_count + 1}/{max_retries}）...")
        if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir):
            print("VASP 运行失败")
            return False

        # 检查收敛
        outcar_path = os.path.join(step_dir, "OUTCAR")
        if check_convergence(outcar_path):
            print("结构优化已收敛！")
            converged = True
            break

        print("结构优化未收敛，准备重试...")
        contcar_path = os.path.join(step_dir, "CONTCAR")
        poscar_path = os.path.join(step_dir, "POSCAR")
        if not os.path.exists(contcar_path):
            print("错误：未找到 CONTCAR 文件")
            return False

        os.remove(poscar_path)
        shutil.move(contcar_path, poscar_path)
        retry_count += 1

    if not converged:
        print("达到最大重试次数，结构优化仍未收敛！")
        return False

    return True


def scf_step(calc_type):
    """自洽场计算"""
    step_dir = "2-scf"
    os.makedirs(step_dir, exist_ok=True)

    # 复制文件
    for f in ["CONTCAR", "POTCAR", "KPOINTS"]:
        src_path = os.path.join("1-opt", f)
        dst_path = os.path.join(step_dir, f)
        if os.path.exists(src_path):
            shutil.copy(src_path, dst_path)

    os.rename(os.path.join(step_dir, "CONTCAR"), os.path.join(step_dir, "POSCAR"))

    # 修改 INCAR
    opt_incar_path = os.path.join("1-opt", "INCAR")
    scf_incar_path = os.path.join(step_dir, "INCAR")
    if not modify_incar_for_scf(opt_incar_path, scf_incar_path, os.path.join(step_dir, "POSCAR")):
        return False

    # 运行 VASP
    print("运行 VASP（自洽场）...")
    if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir):
        return False

    return True


def dos_step(calc_type):
    """态密度计算"""
    step_dir = "3-dos"
    os.makedirs(step_dir, exist_ok=True)

    # 复制文件
    for f in ["POSCAR", "POTCAR", "CHG", "CHGCAR", "WAVECAR"]:
        src_path = os.path.join("2-scf", f)
        dst_path = os.path.join(step_dir, f)
        if os.path.exists(src_path):
            shutil.copy(src_path, dst_path)

    # 修改 INCAR
    scf_incar_path = os.path.join("2-scf", "INCAR")
    dos_incar_path = os.path.join(step_dir, "INCAR")
    if not modify_incar_for_dos(scf_incar_path, dos_incar_path, os.path.join(step_dir, "POSCAR")):
        return False

    # 生成 KPOINTS（基于 opt 倍增）
    print("生成 DOS KPOINTS（基于 OPT 倍增）...")
    if not generate_kpoints_from_opt("1-opt", step_dir):
        return False

    # 运行 VASP
    print("运行 VASP（态密度）...")
    if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir):
        return False

    return True


def md_step(calc_type):
    """分子动力学计算"""
    step_dir = "4-md"
    os.makedirs(step_dir, exist_ok=True)

    # 复制文件
    for f in ["POSCAR", "POTCAR"]:
        src_path = os.path.join("2-scf", f)
        dst_path = os.path.join(step_dir, f)
        if os.path.exists(src_path):
            shutil.copy(src_path, dst_path)

    # 生成 KPOINTS（固定为 1 1 1）
    print("生成 MD KPOINTS（固定为 1 1 1）...")
    if not generate_kpoints_fixed_md(step_dir):
        return False

    # 写入固定内容的 INCAR
    incar_path = os.path.join(step_dir, "INCAR")
    with open(incar_path, 'w') as f:
        f.write(MD_INCAR_CONTENT.strip())

    # 运行 VASP
    print("运行 VASP（分子动力学）...")
    if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir):
        return False

    return True


# ========================
# 主函数
# ========================
def main():
    print("开始执行 VASP 计算流程...")

    # 自动识别计算类型
    calc_type = read_calc_type()
    print(f"识别到计算类型: {calc_type}")

    steps = [
        ("结构优化", lambda: opt_step(calc_type)),
        ("自洽场计算", lambda: scf_step(calc_type)),
        ("态密度计算", lambda: dos_step(calc_type)),
        ("分子动力学", lambda: md_step(calc_type)),
    ]

    for step_name, step_func in steps:
        print(f"\n{'=' * 50}")
        print(f"开始执行: {step_name}")
        print(f"{'=' * 50}")

        if not step_func():
            print(f"错误: {step_name} 执行失败，停止流程")
            return False

        print(f"{step_name} 成功完成")

    print("\n所有计算步骤已完成！")
    return True


if __name__ == "__main__":
    main()