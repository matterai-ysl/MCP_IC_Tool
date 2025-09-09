import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import json
from numba import njit, prange

# ==================================================================================
# ✅ 用户配置部分：请在此修改文件路径和参数
# ==================================================================================
BASE_PATH = r"C:\Users\Administrator\Desktop\data\DFT\MSD"
XDATCAR_FILE = os.path.join(BASE_PATH, "XDATCAR")
INCAR_FILE = os.path.join(BASE_PATH, "INCAR")
OUTPUT_DIR = os.path.join(BASE_PATH, "MSD_output")

USE_LOG_SCALE = True
TIME_UNIT = "ps"

# ==================================================================================
# 🧠 内部函数定义（无需修改）
# ==================================================================================

def read_potim_from_incar(incar_file):
    """
    从 INCAR 文件中读取 POTIM 参数值（时间步长）。
    """
    potim = None
    try:
        with open(incar_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if 'POTIM' in line and '=' in line:
                    value = line.split('=')[1].strip().split()[0]
                    potim = float(value)
                    break
        if potim is None:
            raise ValueError("INCAR 文件中未找到 POTIM 参数，请确认是否设置了 POTIM。")
        print(f"成功读取 POTIM = {potim} {TIME_UNIT}")
        return potim
    except Exception as e:
        raise RuntimeError(f"读取 POTIM 时发生错误: {e}")

def read_xdatcar(xdatcar_file):
    """
    从 XDATCAR 文件中读取晶格信息、元素种类、原子数量和所有帧的坐标。
    """
    with open(xdatcar_file, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    idx = 0
    title = lines[idx].strip()
    idx += 1

    scale = float(lines[idx])
    idx += 1

    lattice = []
    for _ in range(3):
        lattice.append([float(x) for x in lines[idx].split()])
        idx += 1
    lattice = np.array(lattice) * scale

    elements = lines[idx].split()
    idx += 1
    atom_counts = [int(x) for x in lines[idx].split()]
    idx += 1

    atom_types = []
    for element, count in zip(elements, atom_counts):
        atom_types.extend([element] * count)

    frames_fractional = []
    frame_indices = []

    # 找出所有帧头的位置
    while idx < len(lines):
        if lines[idx].startswith("Direct configuration="):
            frame_indices.append(idx)
        idx += 1

    print("Reading XDATCAR frames...")
    for start_idx in tqdm(frame_indices, desc="Processing Frames"):
        idx = start_idx + 1
        coords = []
        while idx < len(lines):
            line = lines[idx].strip()
            if not line:
                idx += 1
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    x, y, z = map(float, parts[:3])
                    coords.append([x, y, z])
                    idx += 1
                except ValueError:
                    break
            else:
                break
        if coords:
            if len(coords) != len(atom_types):
                raise ValueError("帧原子数不一致")
            frames_fractional.append(np.array(coords))

    frames_cartesian = [np.dot(frame, lattice) for frame in frames_fractional]

    return elements, atom_counts, atom_types, frames_cartesian

@njit(parallel=True)
def compute_msd_for_element(traj, n_frames):
    """
    使用 Numba 加速计算 MSD（向量化+并行化）
    traj: (T, M, 3) 的轨迹数组
    """
    msd = np.zeros(n_frames)
    for dt in range(1, n_frames):
        count = 0
        sum_sq = 0.0
        for t0 in range(n_frames - dt):
            for i in range(traj.shape[1]):  # 遍历所有原子
                dx = traj[t0 + dt, i, 0] - traj[t0, i, 0]
                dy = traj[t0 + dt, i, 1] - traj[t0, i, 1]
                dz = traj[t0 + dt, i, 2] - traj[t0, i, 2]
                sum_sq += dx*dx + dy*dy + dz*dz
                count += 1
        if count > 0:
            msd[dt] = sum_sq / count
    return msd

def compute_msd_by_element(frames, atom_types):
    """
    按元素分类计算 MSD（使用 Numba 加速）
    """
    n_frames = len(frames)
    elements = set(atom_types)
    msd_dict = {}

    all_positions = np.array(frames)  # (T, N, 3)

    # 删除了重复的 print 语句
    for element in tqdm(elements, desc="Processing Elements"):  # 修改了 tqdm 的 desc
        indices = np.array([i for i, t in enumerate(atom_types) if t == element])
        if len(indices) == 0:
            continue

        element_traj = all_positions[:, indices, :]  # (T, M, 3)
        msd = compute_msd_for_element(element_traj, n_frames)
        msd_dict[element] = msd

    return msd_dict

def plot_msd_by_element(msd_dict, dt=1.0, log_scale=True, output_dir="MSD_output"):
    os.makedirs(output_dir, exist_ok=True)
    times = np.arange(len(next(iter(msd_dict.values())))) * dt
    plt.figure(figsize=(10, 6))

    for element, msd in msd_dict.items():
        plt.plot(times, msd, 'o-', label=element, markersize=4)

    plt.xlabel(f"Time ({TIME_UNIT})")
    plt.ylabel("MSD (Å²)")
    plt.title("MSD by Element")
    if log_scale:
        plt.xscale('log')
        plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # 保存图像
    plot_path = os.path.join(output_dir, "msd_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.show()

def save_msd_to_files(msd_dict, dt=1.0, output_dir="MSD_output"):
    os.makedirs(output_dir, exist_ok=True)
    times = np.arange(len(next(iter(msd_dict.values())))) * dt
    print(f"Saving MSD data to {output_dir}...")

    for element, msd in tqdm(msd_dict.items(), desc="Saving MSD Files"):
        data = {
            "Time (ps)": times.tolist(),
            "MSD (Å²)": msd.tolist()
        }
        filename = os.path.join(output_dir, f"msd_{element}.json")
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4)
    print("MSD 数据保存完成。")

# ==================================================================================
# 🚀 主程序入口
# ==================================================================================
if __name__ == "__main__":
    try:
        # 读取 POTIM（时间步长）
        dt = read_potim_from_incar(INCAR_FILE)

        # 读取 XDATCAR 数据
        print("Reading XDATCAR file...")
        elements, atom_counts, atom_types, frames = read_xdatcar(XDATCAR_FILE)

        # 计算 MSD（按元素分类）
        print("Calculating MSD by element...")
        msd_dict = compute_msd_by_element(frames, atom_types)

        # 保存 MSD 数据到文件
        save_msd_to_files(msd_dict, dt=dt, output_dir=OUTPUT_DIR)

        # 绘制 MSD 曲线
        print("Plotting MSD by element...")
        plot_msd_by_element(msd_dict, dt=dt, log_scale=USE_LOG_SCALE)

    except Exception as e:
        print(f"ERROR: {e}")