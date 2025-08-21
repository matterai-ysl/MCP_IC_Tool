import os
import subprocess
import sys
import math
import json
from .Config import get_kpoints_config,get_vasp_config,get_incar_templates
# ========================
# 工具函数
# ========================


def run_command_background(command, work_dir=None, log_file="result.log"):
    """后台运行命令，立即返回"""
    full_log_path = os.path.join(work_dir, log_file) if work_dir else log_file
    
    try:
        print(f"后台执行命令: {' '.join(command)}")
        
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
def cif_to_poscar(cif_path,work_dir='./'):
    """将CIF转换为POSCAR（vasp5格式）"""
    print("正在转换CIF到POSCAR...")
    if not os.path.exists(cif_path):
        sys.exit("错误：未找到文件 {}".format(cif_path))
    
    # 获取CIF文件所在目录
    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    poscar_path = os.path.join(work_dir, "POSCAR")
    
    # 保存当前目录
    original_dir = os.getcwd()
    try:    
        # 严格按照用户提供的固定命令格式
        command = [
            "cif2cell",
            "-f", cif_path,
            "-p", "vasp",
            "-o", poscar_path,
            "--no-reduce",
            "--vasp-format=5",
            "--vasp-cartesian-lattice-vectors"
        ]
        subprocess.run(
                command,
                cwd=work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
        print("POSCAR生成成功")
        return poscar_path
    finally:
        # 确保无论发生什么都恢复原始目录
        os.chdir(original_dir)

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

def generate_kpoints(work_dir='./'):

    """生成 KPOINTS 文件"""
    poscar_path = os.path.join(work_dir, "POSCAR")
    target_product = get_kpoints_config()["TARGET_KPOINT_PRODUCT"]
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
        print(f"KPOINTS 文件已生成于 {kpoints_path}")
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
        print(f"POTCAR 文件已生成于 {potcar_path}")
        return True
    except IOError as e:
        raise RuntimeError(f"写入 POTCAR 文件失败: {str(e)}")

def generate_incar(work_dir, calc_type):
    """根据 calc_type 生成 INCAR 文件"""
    incar_path = os.path.join(work_dir, "INCAR")
    with open(incar_path, 'w') as f:
        f.write(get_incar_templates()[calc_type].strip())
    u_values = load_u_values(get_vasp_config()["paths"]["U_VALUES_JSON"])
    elements = read_poscar_elements(os.path.join(work_dir, "POSCAR"))
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

if __name__ == "__main__":
    cif_to_poscar("/Users/ysl/Desktop/Code/MCP_IC_Tool/downloads/Li2O_mp-1960_Fm-3m.cif")