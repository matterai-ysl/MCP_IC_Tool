#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import sys
import csv
import numpy as np
import datetime
import logging

# ========================
# 配置参数（请根据您的环境修改）
# ========================
# VASP可执行文件路径
VASP_PATH = "/gpfs/home/lixingyu/software/vasp.6.3.2/vasp.6.3.2/bin/vasp_std"

# Bader分析相关程序路径
CHGSUM_PL_PATH = "/gpfs/home/lixingyu/bin/chgsum.pl"
BADER_PATH = "/gpfs/home/lixingyu/bin/bader"

# 结构优化最大重试次数
MAX_RETRIES = 3


# ========================
# 通用工具函数
# ========================
def run_command(command, cwd=None, log_file="vasp_run.log"):
    """运行系统命令并将详细输出重定向到指定文件。"""
    full_log_path = os.path.join(cwd, log_file) if cwd else log_file
    try:
        logging.info("  - 执行命令: {} (详细日志: {})".format(' '.join(command), full_log_path))
        with open(full_log_path, 'w') as f:
            result = subprocess.run(
                command, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, text=True, check=True
            )
        return True
    except FileNotFoundError:
        logging.error("  - 错误: 命令未找到 -> {}。请检查路径配置。".format(command[0]));
        return False
    except subprocess.CalledProcessError:
        logging.error("  - 命令执行失败。请检查详细日志文件: {}".format(full_log_path));
        return False
    except Exception as e:
        logging.error("  - 运行命令时发生未知错误: {}".format(str(e)));
        return False


def check_convergence(outcar_path):
    """检查 OUTCAR 文件是否收敛。"""
    if not os.path.exists(outcar_path): return False
    try:
        with open(outcar_path, 'rb') as f:
            f.seek(-2048, os.SEEK_END)
            for line in reversed(f.readlines()):
                if b"Voluntary context switch" in line:
                    return True
        return False
    except Exception as e:
        logging.error("  - 读取 OUTCAR 时出错: {}".format(str(e)));
        return False


def get_kmesh(kpoints_path):
    try:
        with open(kpoints_path, 'r') as f:
            lines = f.readlines()
        if len(lines) >= 4:
            parts = lines[3].strip().split()
            if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
                return [int(p) for p in parts[:3]]
        return None
    except Exception as e:
        logging.error("  - 读取KPOINTS失败: {}".format(str(e)));
        return None


def generate_doubled_kpoints(src_kpoints, dest_kpoints, comment="Automatically generated"):
    kmesh = get_kmesh(src_kpoints)
    if not kmesh:
        logging.error("  - 错误: 无法从源KPOINTS {} 提取k点网格。".format(src_kpoints));
        return False
    kmesh_doubled = [k * 2 for k in kmesh]
    try:
        with open(dest_kpoints, 'w') as f:
            f.write("{}\n0\nGamma\n {}\n 0 0 0\n".format(comment, ' '.join(map(str, kmesh_doubled))))
        logging.info("  - 新的 KPOINTS 已生成于 {}".format(dest_kpoints));
        return True
    except Exception as e:
        logging.error("  - 生成 KPOINTS 文件失败: {}".format(str(e)));
        return False


def parse_generic_poscar(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        elements = lines[5].strip().split()
        num_atoms = [int(x) for x in lines[6].strip().split()]
        if len(elements) != len(num_atoms) or any(char.isdigit() for char in "".join(elements)):
            raise ValueError("POSCAR格式可能为V4或第6行非元素行")
        atom_list = [symbol for symbol, count in zip(elements, num_atoms) for _ in range(count)]
        logging.info("  - 从 {} 成功解析V5格式原子信息。".format(filename))
        return elements, atom_list
    except (ValueError, IndexError):
        try:
            num_atoms = [int(x) for x in lines[5].strip().split()]
            elements = ["Elem{}".format(i + 1) for i in range(len(num_atoms))]
            atom_list = [symbol for symbol, count in zip(elements, num_atoms) for _ in range(count)]
            logging.info("  - 从 {} 按V4格式解析原子信息 (无元素符号)。".format(filename))
            return elements, atom_list
        except Exception as e:
            logging.error("  - 错误: 解析 {} 失败: {}".format(filename, e));
            return None, None


def parse_potcar_zval(elements, filename):
    zval_map = {}
    try:
        with open(filename, 'r') as f:
            content = f.read()
        potentials = content.split("End of Dataset")
        logging.info("  - 从 {} 解析 ZVAL 信息:".format(filename))
        for i, symbol in enumerate(elements):
            for line in potentials[i].splitlines():
                if "ZVAL" in line:
                    zval_part = line.strip().split(';')[1]
                    zval_string = zval_part.split()[2]
                    zval = float(zval_string)
                    zval_map[symbol] = zval
                    logging.info("    - 元素: {}, ZVAL = {}".format(symbol, zval))
                    break
        return zval_map
    except Exception as e:
        logging.error("  - 错误: 解析 {} 失败: {}".format(filename, e));
        return None


# ========================
# 第 I 部分: VASP 计算步骤
# ========================
def opt_step():
    step_dir = "1-opt"
    os.makedirs(step_dir, exist_ok=True)
    logging.info("  - 准备文件: 复制POSCAR，移动INCAR, KPOINTS, POTCAR...")
    shutil.copy("POSCAR", os.path.join(step_dir, "POSCAR"))
    for f in ["INCAR", "KPOINTS", "POTCAR"]: shutil.move(f, os.path.join(step_dir, f))

    for i in range(MAX_RETRIES):
        logging.info("  - --- 结构优化尝试: {}/{} ---".format(i + 1, MAX_RETRIES))
        if i > 0:
            contcar, poscar = os.path.join(step_dir, "CONTCAR"), os.path.join(step_dir, "POSCAR")
            if os.path.exists(contcar) and os.path.getsize(contcar) > 0:
                shutil.copy(contcar, poscar);
                logging.info("    - CONTCAR 已复制为 POSCAR, 准备重试...")
            else:
                logging.error("    - 错误: CONTCAR不存在或为空, 无法重试。");
                return False

        if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir, log_file="opt_vasp.log"):
            logging.warning("    - VASP 运行失败或提前退出")

        if check_convergence(os.path.join(step_dir, "OUTCAR")):
            logging.info("  - 结构优化已收敛！")
            logging.info("  - 清理不必要的文件...")
            files_to_keep = ["INCAR", "KPOINTS", "POTCAR", "POSCAR", "CONTCAR", "OUTCAR", "opt_vasp.log", "OSZICAR"]
            for f in os.listdir(step_dir):
                if f not in files_to_keep: os.remove(os.path.join(step_dir, f))
            return True
        logging.warning("  - 结构优化未收敛...")

    logging.error("  - 达到最大重试次数，结构优化失败！");
    return False


def scf_step():
    step_dir, opt_dir = "2-scf", "1-opt"
    os.makedirs(step_dir, exist_ok=True)
    opt_contcar = os.path.join(opt_dir, "CONTCAR")
    if not (os.path.exists(opt_contcar) and os.path.getsize(opt_contcar) > 0):
        logging.error("  - 错误: {} 不存在或为空！".format(opt_contcar));
        return False
    shutil.copy(opt_contcar, os.path.join(step_dir, "POSCAR"))
    for f in ["INCAR", "POTCAR"]: shutil.copy(os.path.join(opt_dir, f), step_dir)
    if not generate_doubled_kpoints(os.path.join(opt_dir, "KPOINTS"), os.path.join(step_dir, "KPOINTS")): return False
    logging.info("  - 修改 INCAR 用于静态自洽计算...")
    scf_incar = os.path.join(step_dir, "INCAR")
    try:
        with open(scf_incar, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            s = line.strip().upper()
            if s.startswith("SYSTEM"):
                new_lines.append("SYSTEM = SCF\n")
            elif s.startswith("NSW"):
                new_lines.append("NSW = 0\n")
            elif s.startswith("IBRION"):
                new_lines.append("IBRION = -1\n")
            elif s.startswith("ISIF"):
                new_lines.append("ISIF = 2\n")
            elif s.startswith("POTIM"):
                new_lines.append("POTIM = 0\n")
            elif s.startswith("EDIFFG"):
                new_lines.append("#EDIFFG = -0.01\n")
            elif s.startswith("LWAVE"):
                new_lines.append("LWAVE = .TRUE.\n")
            elif s.startswith("LCHARG"):
                new_lines.append("LCHARG = .TRUE.\n")
            else:
                new_lines.append(line)
        if not any("LAECHG" in l.upper() for l in new_lines): new_lines.append("LAECHG = .TRUE.\n")
        if not any("LELF" in l.upper() for l in new_lines): new_lines.append("LELF = .TRUE.\n")
        with open(scf_incar, 'w') as f:
            f.writelines(new_lines)
    except Exception as e:
        logging.error("  - 修改 INCAR 失败: {}".format(str(e)));
        return False
    if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir, log_file="scf_vasp.log"): return False
    return True


def dos_step():
    step_dir, scf_dir = "3-dos", "2-scf"
    os.makedirs(step_dir, exist_ok=True)
    logging.info("  - 为DOS计算准备输入文件...")
    for f in ["POSCAR", "POTCAR", "CHGCAR", "WAVECAR", "INCAR", "KPOINTS"]:
        src, dest = os.path.join(scf_dir, f), os.path.join(step_dir, f)
        if os.path.exists(src):
            shutil.copy(src, dest);
            logging.info("    - 已复制: {}".format(f))
        elif f in ["POSCAR", "POTCAR", "INCAR", "KPOINTS"]:
            logging.error("    - 错误: 关键文件 {} 不存在，无法进行DOS计算！".format(src));
            return False
    logging.info("  - 修改 INCAR 用于态密度计算...")
    dos_incar = os.path.join(step_dir, "INCAR")
    try:
        with open(dos_incar, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            s = line.strip().upper()
            if s.startswith("SYSTEM"):
                new_lines.append("SYSTEM = DOS\n")
            elif s.startswith("ISTART"):
                new_lines.append("ISTART = 1\n")
            elif s.startswith("ICHARG"):
                new_lines.append("ICHARG = 11\n")
            elif "LAECHG" in s or "LELF" in s:
                continue
            else:
                new_lines.append(line)
        if not any("NEDOS" in l.upper() for l in new_lines): new_lines.append("NEDOS = 5000\n")
        if not any("LORBIT" in l.upper() for l in new_lines): new_lines.append("LORBIT = 11\n")
        with open(dos_incar, 'w') as f:
            f.writelines(new_lines)
    except Exception as e:
        logging.error("  - 修改 INCAR 失败: {}".format(str(e)));
        return False
    if not run_command(["mpirun", "-bootstrap", "lsf", VASP_PATH], cwd=step_dir, log_file="dos_vasp.log"): return False
    return True


# ========================
# 第 II 部分: 数据处理步骤
# ========================
def run_bader_analysis_step():
    scf_dir = "2-scf"
    logging.info("  - 在 {} 目录内运行Bader电荷分析程序...".format(scf_dir))
    for f in ["AECCAR0", "AECCAR2", "CHGCAR"]:
        if not os.path.exists(os.path.join(scf_dir, f)): logging.error(
            "  - 错误: 未找到Bader分析所需文件 {}。".format(f)); return False
    chgsum_cmd = ["perl", CHGSUM_PL_PATH, "AECCAR0", "AECCAR2"]
    if not run_command(chgsum_cmd, cwd=scf_dir, log_file="chgsum.log"): return False
    if not os.path.exists(os.path.join(scf_dir, "CHGCAR_sum")): logging.error(
        "  - 错误: 未生成 CHGCAR_sum 文件。"); return False
    bader_cmd = [BADER_PATH, "CHGCAR", "-ref", "CHGCAR_sum"]
    if not run_command(bader_cmd, cwd=scf_dir, log_file="bader.log"): return False
    return True


def generate_bader_csv_step():
    scf_dir = "2-scf"
    poscar_path, potcar_path, acf_path = os.path.join(scf_dir, "POSCAR"), os.path.join(scf_dir, "POTCAR"), os.path.join(
        scf_dir, "ACF.dat")
    output_csv = "bader_charge_analysis.csv"
    for f in [poscar_path, potcar_path, acf_path]:
        if not os.path.exists(f): logging.error("  - 错误: 分析所需文件 {} 不存在！".format(f)); return False
    elements, atom_list = parse_generic_poscar(poscar_path)
    if not elements: return False
    zval_map = parse_potcar_zval(elements, potcar_path)
    if not zval_map: return False
    try:
        with open(acf_path, 'r') as f_in, open(output_csv, 'w', newline='') as f_out:
            reader = csv.reader(f_in, delimiter=' ', skipinitialspace=True)
            writer = csv.writer(f_out)
            writer.writerow(["ATOM_INDEX", "ELEMENT", "X", "Y", "Z", "ZVAL", "CHARGE", "CHARGE_TRANSFER"])
            for row in reader:
                row_data = [item for item in row if item]
                if not row_data: continue
                try:
                    atom_index = int(row_data[0])
                except (ValueError, IndexError):
                    continue
                symbol = atom_list[atom_index - 1]
                zval = zval_map[symbol]
                charge = float(row_data[4])
                transfer = zval - charge
                writer.writerow([atom_index, symbol, row_data[1], row_data[2], row_data[3], zval, charge,
                                 "{:.6f}".format(transfer)])
        logging.info("  - Bader分析完成！结果已保存到根目录的 {} 文件中。".format(output_csv));
        return True
    except Exception as e:
        logging.error("  - 生成Bader CSV报告时出错: {}".format(e));
        return False


def generate_dos_data_step():
    dos_dir, output_dir = "3-dos", "dos-data"
    os.makedirs(output_dir, exist_ok=True)
    poscar_path, doscar_path, eigenval_path = os.path.join(dos_dir, "POSCAR"), os.path.join(dos_dir,
                                                                                            "DOSCAR"), os.path.join(
        dos_dir, "EIGENVAL")
    if not os.path.exists(poscar_path) or not os.path.exists(doscar_path):
        logging.error("  - 错误: DOS分析所需文件 {} 或 {} 未找到".format(poscar_path, doscar_path));
        return False
    _, atom_types = parse_generic_poscar(poscar_path)
    if not atom_types: return False
    with open(doscar_path, 'r') as f:
        lines = f.readlines()
    natom, nl, efermi = int(lines[0].split()[0]), int(lines[5].split()[2]), float(lines[5].split()[3])
    dos_raw, start_line = [], 6
    for _ in range(natom + 1):
        dos_raw.append(np.loadtxt(lines[start_line:start_line + nl]));
        start_line += nl + 1
    nspin = 2 if dos_raw[0].shape[1] > 2 else 1
    energy = dos_raw[0][:, 0] - efermi
    element_pdos = {'TOTAL': {'energy': energy, 'total_up': dos_raw[0][:, 1],
                              'total_down': dos_raw[0][:, 2] if nspin == 2 else np.zeros_like(dos_raw[0][:, 1])}}
    for i in range(1, natom + 1):
        atom_dos, elem = dos_raw[i], atom_types[i - 1]
        s_up = atom_dos[:, 1];
        p_up = np.sum(atom_dos[:, 2:5], axis=1);
        d_up = np.sum(atom_dos[:, 5:10], axis=1)
        s_down, p_down, d_down = (atom_dos[:, 2], np.sum(atom_dos[:, 3:9:2], axis=1),
                                  np.sum(atom_dos[:, 9::2], axis=1)) if nspin == 2 and atom_dos.shape[1] > 9 else (
            np.zeros_like(s_up), np.zeros_like(p_up), np.zeros_like(d_up))
        if elem not in element_pdos:
            element_pdos[elem] = {k: np.zeros_like(energy) for k in
                                  ['s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']}
        for k, v in zip(['s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down'],
                        [s_up, s_down, p_up, p_down, d_up, d_down]):
            element_pdos[elem][k] += v
    for elem, pdos in element_pdos.items():
        if elem == 'TOTAL': continue
        pdos['energy'] = energy
        pdos['total_up'] = pdos['s_up'] + pdos['p_up'] + pdos['d_up']
        pdos['total_down'] = pdos['s_down'] + pdos['p_down'] + pdos['d_down']
    for elem, pdos in element_pdos.items():
        with open(os.path.join(output_dir, "{}-DOS.csv".format(elem)), 'w', newline='') as f:
            writer = csv.writer(f)
            if elem == 'TOTAL':
                writer.writerow(['Energy (eV)', 'Total DOS (up)', 'Total DOS (down)'])
                writer.writerows(np.column_stack([pdos['energy'], pdos['total_up'], -pdos['total_down']]))
            else:
                writer.writerow(
                    ['Energy (eV)', 'Total_up', 'Total_down', 's_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down'])
                writer.writerows(np.column_stack(
                    [pdos['energy'], pdos['total_up'], -pdos['total_down'], pdos['s_up'], -pdos['s_down'], pdos['p_up'],
                     -pdos['p_down'], pdos['d_up'], -pdos['d_down']]))
    logging.info("  - DOS数据已成功保存到 {} 目录".format(output_dir))
    if os.path.exists(eigenval_path):
        try:
            with open(eigenval_path, 'r') as f:
                lines = f.readlines()
            nkpoints, nbands = int(lines[5].split()[1]), int(lines[5].split()[2])
            vbm, cbm = -np.inf, np.inf
            line_idx = 7
            for _ in range(nkpoints):
                for _ in range(nbands):
                    energy = float(lines[line_idx].split()[1]);
                    line_idx += 1
                    if energy <= efermi and energy > vbm: vbm = energy
                    if energy > efermi and energy < cbm: cbm = energy
                line_idx += 2
            with open(os.path.join(output_dir, "gap.txt"), 'w') as f:
                f.write(
                    "Fermi Level (from DOSCAR): {:.6f} eV\nVBM: {:.6f} eV\nCBM: {:.6f} eV\nBand Gap: {:.6f} eV\n".format(
                        efermi, vbm, cbm, cbm - vbm))
            logging.info("  - 能带间隙信息已保存到 {}/gap.txt".format(output_dir))
        except Exception as e:
            logging.error("  - 解析 EIGENVAL 文件时出错: {}".format(e))
    return True


# ========================
# 主函数
# ========================
def main():
    log_filename = 'workflow.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',  # 简化格式，不再每行都加时间戳
        handlers=[logging.FileHandler(log_filename, mode='w'), logging.StreamHandler(sys.stdout)]
    )

    if not all(os.path.exists(f) for f in ["POSCAR", "INCAR", "KPOINTS", "POTCAR"]):
        logging.error("错误：启动前请确保 POSCAR, INCAR, KPOINTS, POTCAR 文件都在当前目录。");
        return

    logging.info("====== VASP自动化流程启动 ======")
    total_start_time = datetime.datetime.now()

    calculation_steps = [("1. 结构优化", opt_step), ("2. 自洽计算", scf_step), ("3. 态密度计算", dos_step)]
    analysis_steps = [("A. 运行Bader分析程序", run_bader_analysis_step),
                      ("B. 生成Bader电荷转移CSV", generate_bader_csv_step),
                      ("C. 生成DOS分析数据", generate_dos_data_step)]
    durations = {}

    logging.info("\n>>>>>>>>> 第 I 部分: VASP 计算 <<<<<<<<<")
    logging.info(total_start_time.strftime("%Y-%m-%d %H:%M:%S"))  # 只在开头打印一次时间
    for name, func in calculation_steps:
        logging.info("\n--- 开始计算步骤: {} ---".format(name))
        start_time = datetime.datetime.now()
        if not func():
            logging.error("\n错误: 计算步骤 '{}' 失败，流程中止。".format(name));
            return
        end_time = datetime.datetime.now()
        durations[name] = end_time - start_time
        logging.info("--- 计算步骤 '{}' 成功完成。 ---".format(name))

    logging.info("\n>>>>>>>>> 第 II 部分: 数据处理与分析 <<<<<<<<<")
    analysis_start_time = datetime.datetime.now()
    for name, func in analysis_steps:
        logging.info("\n--- 开始分析步骤: {} ---".format(name))
        if not func():
            logging.error("\n错误: 分析步骤 '{}' 失败，流程中止。".format(name));
            return
        logging.info("--- 分析步骤 '{}' 成功完成。 ---".format(name))
    analysis_end_time = datetime.datetime.now()
    durations["数据分析"] = analysis_end_time - analysis_start_time
    total_end_time = datetime.datetime.now()
    total_duration = total_end_time - total_start_time

    logging.info("\n" + "=" * 50)
    logging.info("祝贺您！所有计算和分析步骤已成功完成！")
    logging.info("=" * 50)
    logging.info("               流程耗时总结")
    logging.info("-" * 50)
    logging.info("总耗时: {}".format(total_duration))
    logging.info("-" * 50)
    for name, func in calculation_steps:
        logging.info("计算步骤 '{}' 耗时: {}".format(name, durations.get(name, "N/A")))
    logging.info("所有数据分析步骤 总耗时: {}".format(durations.get("数据分析", "N/A")))
    logging.info("=" * 50)


if __name__ == "__main__":
    main()