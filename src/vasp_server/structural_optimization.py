import os
import shutil
from .base import run_command_background,generate_kpoints, generate_potcar,generate_incar
from .Config import get_vasp_config




def run_opt_step(work_dir,calc_type):
    """结构优化步骤（带收敛检测和重试机制）"""
    step_dir = os.path.join(work_dir, "1-opt")
    os.makedirs(step_dir, exist_ok=True)

    # 初始文件准备
    shutil.copy(os.path.join(work_dir, "POSCAR"), step_dir)

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
        if not generate_kpoints(step_dir):
            return False
        if not generate_potcar(step_dir):
            return False

        # 生成 INCAR
        if not generate_incar(step_dir, calc_type):
            return False

        # 运行 VASP
        print(f"运行 VASP（尝试次数: {retry_count + 1}/{max_retries}）...")
        result = run_command_background(["mpirun", "-bootstrap", "lsf", get_vasp_config()["VASP_PATH"]], work_dir=step_dir)
        if "error" in result.keys():
            print(f"VASP 运行失败: {result['error']}")
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