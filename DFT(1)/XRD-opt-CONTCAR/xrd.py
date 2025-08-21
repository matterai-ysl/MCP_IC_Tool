from pymatgen.core import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import json
import os

# 1. 读取 CONTCAR 文件
file_path = r"C:\Users\Administrator\Desktop\data\DFT\XRD-opt-CONTCAR\CONTCAR"
structure = Structure.from_file(file_path)

# 2. 初始化 XRD 计算器
xrd_calculator = XRDCalculator()

# 3. 计算衍射图谱
pattern = xrd_calculator.get_pattern(structure, two_theta_range=(10, 80))

# 4. 提取数据并扩展为完整角度范围
two_theta_full = np.linspace(10, 80, 1000)
intensity_full = np.zeros_like(two_theta_full)

# 将原始数据映射到完整数组中（最近邻插值）
for theta, intensity in zip(pattern.x, pattern.y):
    idx = np.abs(two_theta_full - theta).argmin()
    intensity_full[idx] = intensity

two_thetas = two_theta_full
intensities = intensity_full

# 5. 创建输出目录
output_dir = os.path.join(os.path.dirname(file_path), "XRD_output")
os.makedirs(output_dir, exist_ok=True)

# 6. 创建平滑曲线（样条插值）
spline = UnivariateSpline(two_thetas, intensities, s=0.5)
x_smooth = np.linspace(min(two_thetas), max(two_thetas), 1000)
y_smooth = spline(x_smooth)

# 7. 绘图配置（严格遵循需求：灰线+红曲线）
plt.figure(figsize=(12, 5))

# ✅ XRD Peaks：灰色垂直线（无标记，仅线段）
for x, y in zip(two_thetas, intensities):
    if y > 0:  # 仅绘制有信号的峰
        plt.plot([x, x], [0, y],
                 color='gray',       # 灰色垂线
                 linewidth=1.2,
                 linestyle='-',
                 zorder=2)           # 确保在曲线之下

# ✅ 手动添加灰色填充+黑色边框的标记
plt.scatter(two_thetas[intensities > 0], intensities[intensities > 0],
            color='lightgray',      # 填充色（浅灰）
            edgecolor='black',      # 边框色（黑）
            linewidth=0.8,
            marker='o',             # 圆形标记
            s=40,                   # 标记大小
            zorder=3,               # 确保在垂线上层
            label='XRD Peaks')      # 图例标签

# ✅ Fitted Curve：红色平滑曲线
plt.plot(x_smooth, y_smooth,
         color='#FF3B30',          # 珊瑚红色
         linewidth=3,              # 加粗线条
         linestyle='-',            # 实线
         alpha=0.95,               # 微调透明度
         zorder=1,                 # 确保在垂线下层
         label='Fitted Curve')     # 图例标签

# 设置坐标轴标签和标题
plt.xlabel("2θ (degrees)", fontsize=12, weight='bold')
plt.ylabel("Intensity (a.u.)", fontsize=12, weight='bold')
plt.title("XRD Peaks + Fitted Curve", fontsize=14, pad=20, weight='bold')

# 设置网格和范围
plt.grid(True, linestyle='--', alpha=0.6, color='gray')
plt.ylim(0, max(intensities) * 1.15)
plt.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, shadow=True)
plt.tight_layout()

# 保存图像
image_path = os.path.join(output_dir, "xrd_peaks.png")
plt.savefig(image_path, dpi=600, format='png', bbox_inches='tight')
print(f"图像已保存至: {image_path}")

# 保存数据为 JSON
data = {
    "two_theta": two_thetas.tolist(),
    "intensity": intensities.tolist(),
    "smooth_two_theta": x_smooth.tolist(),
    "smooth_intensity": y_smooth.tolist()
}

json_path = os.path.join(output_dir, "xrd_data.json")
with open(json_path, 'w') as f:
    json.dump(data, f, indent=2)
print(f"数据已保存至: {json_path}")

# 显示图像
plt.show()