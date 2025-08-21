# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from matplotlib.font_manager import FontProperties
import os

# ================================
# 文件路径配置
# ================================
DATA_FILE_PATH = r"C:\Users\Administrator\Desktop\python\Raman_tool\Raman.xlsx"

# ================================
# 配置区
# ================================
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False


class Config:
    FIXED_PEAKS = [1200, 1500, 1620]
    VARIABLE_RANGES = [(1250, 1450), (1550, 1600)]
    MAX_SIGMA = {'fixed': 100, 'variable': 80}
    MIN_SIGMA = {'fixed': 2, 'variable': 5}
    PEAK_SPACING = 30
    REG_WEIGHT = 1e-6
    FIG_SIZE = (16, 14)
    COLOR_THEME = {
        'plot_bg': '#f5f5f5',
        'grid_color': '#dddddd',
        'table_header': '#606060',
        'total_row': '#f0f0f0'
    }
    SAVE_DPI = 300
    EXPORT_DATA = True


# ================================
# 数据加载模块
# ================================
def load_raman_data():
    try:
        file_path = DATA_FILE_PATH
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"找不到文件: {file_path}")

        data = pd.read_excel(file_path, header=None)
        x = data[0].values.astype(float)
        y = data[1].values.astype(float)

        if len(x) != len(y):
            raise ValueError("X和Y数据长度不一致")

        data_scale = np.max(np.abs(y)) + 1e-9
        return x, y / data_scale, data_scale, file_path
    except Exception as e:
        print(f"数据加载失败: {str(e)}")
        exit(1)


# ================================
# 优化算法模块
# ================================
class PeakOptimizer:
    def __init__(self, x, y_normalized):
        self.x = x
        self.y = y_normalized
        self.bounds = self._init_bounds()

    def _init_bounds(self):
        bounds = (
                [(0.05, 15)] * 3 +
                [(Config.MIN_SIGMA['fixed'] ** 2, Config.MAX_SIGMA['fixed'] ** 2)] * 3 +
                [(r[0], r[1]) for r in Config.VARIABLE_RANGES] +
                [(0.05, 15)] * 2 +
                [(Config.MIN_SIGMA['variable'] ** 2, Config.MAX_SIGMA['variable'] ** 2)] * 2
        )
        return bounds

    def optimize(self):
        print("正在进行全局优化...")
        result_de = differential_evolution(
            self.objective, self.bounds, args=(self.x, self.y),
            strategy='best1bin', maxiter=1500, popsize=75,
            mutation=(0.7, 1.2), recombination=0.9, seed=42, polish=True
        )
        print("正在进行局部优化...")
        result = minimize(
            self.objective, result_de.x, args=(self.x, self.y),
            method='L-BFGS-B', bounds=self.bounds,
            options={'maxiter': 5000, 'ftol': 1e-8}
        )
        return result

    @staticmethod
    def objective(params, x, y_normalized):
        n_fixed = len(Config.FIXED_PEAKS)
        fixed_weights = params[:n_fixed]
        fixed_vars = params[n_fixed:2 * n_fixed]
        var_pos = params[2 * n_fixed:2 * n_fixed + 2]
        var_weights = params[2 * n_fixed + 2:2 * n_fixed + 4]
        var_vars = params[2 * n_fixed + 4:2 * n_fixed + 6]

        gaussian = lambda x, H, mu, s: H * np.exp(-(x - mu) ** 2 / (2 * s))
        total = sum(gaussian(x, H, Config.FIXED_PEAKS[i], s)
                    for i, (H, s) in enumerate(zip(fixed_weights, fixed_vars)))
        total += sum(gaussian(x, H, mu, s)
                     for H, mu, s in zip(var_weights, var_pos, var_vars))

        weights = np.ones_like(x)
        for mu in [*Config.FIXED_PEAKS, *var_pos]:
            peak_region = np.abs(x - mu) < 3 * np.sqrt(Config.MAX_SIGMA['variable'])
            weights[peak_region] = 5.0

        mse = np.average((total - y_normalized) ** 2, weights=weights)
        reg = Config.REG_WEIGHT * (np.sum(fixed_weights ** 2) + np.sum(var_weights ** 2))

        penalty = 0
        all_peaks = sorted([*Config.FIXED_PEAKS, *var_pos])
        for i in range(1, len(all_peaks)):
            if (spacing := all_peaks[i] - all_peaks[i - 1]) < Config.PEAK_SPACING:
                penalty += (Config.PEAK_SPACING - spacing) ** 2 * 10

        for s in [*fixed_vars, *var_vars]:
            if s < (Config.MIN_SIGMA['fixed'] ** 2):
                penalty += 100 * ((Config.MIN_SIGMA['fixed'] ** 2) - s)
            if s > (Config.MAX_SIGMA['variable'] ** 2):
                penalty += 100 * (s - (Config.MAX_SIGMA['variable'] ** 2))

        return mse + reg + penalty


# ================================
# 统计量计算模块
# ================================
def calculate_statistics(x, y, y_fit):
    n = len(y)
    p = 12
    dof = n - p

    residuals = y - y_fit
    ss_res = np.sum(residuals ** 2)
    y_mean = np.mean(y)
    ss_tot = np.sum((y - y_mean) ** 2)
    ss_reg = np.sum((y_fit - y_mean) ** 2)

    r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else np.inf
    adj_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - p) if n > p else np.inf
    chi2 = np.sum((residuals ** 2) / y_fit) if np.all(y_fit > 0) else np.inf

    return {
        'n': n, 'p': p, 'dof': dof, 'ss_res': ss_res,
        'ss_tot': ss_tot, 'ss_reg': ss_reg, 'r2': r_squared,
        'adj_r2': adj_r_squared, 'chi2': chi2
    }


# ================================
# 可视化与数据导出模块
# ================================
class RamanVisualizer:
    def __init__(self, peaks_info, x, y, data_scale, file_path, y_fit_original_x):
        self.peaks_info = peaks_info
        self.x = x
        self.y = y
        self.data_scale = data_scale
        self.file_path = file_path
        self.save_dir = os.path.dirname(file_path)
        self.base_name = os.path.splitext(os.path.basename(file_path))[0]
        self.x_fine = np.linspace(x.min(), x.max(), 1000)
        self.stats = None
        self.residuals = y - y_fit_original_x
        self.y_fit_original_x = y_fit_original_x

    def create_figure(self, stats=None):
        self.stats = stats
        plt.figure(figsize=Config.FIG_SIZE)
        gs = plt.GridSpec(2, 1, height_ratios=[4, 1], hspace=0.12)

        ax1 = plt.subplot(gs[0])
        self._setup_spectrum_plot(ax1)

        ax2 = plt.subplot(gs[1])
        self._setup_result_table(ax2)

        plt.tight_layout(pad=4)
        plt.subplots_adjust(top=0.95)
        self.save_plot()
        plt.show()

    def _setup_spectrum_plot(self, ax):
        ax.set_facecolor(Config.COLOR_THEME['plot_bg'])
        ax.set_title('拉曼光谱多峰拟合（含残差曲线）', fontsize=18, pad=15)
        ax.set_xlabel('拉曼位移', fontsize=16, labelpad=12)
        ax.set_ylabel('强度', fontsize=16, labelpad=12)
        ax.grid(True, color=Config.COLOR_THEME['grid_color'], linestyle=':')
        ax.tick_params(axis='both', which='major', labelsize=14)

        total_fine = np.zeros_like(self.x_fine)
        line_styles = ['--', '-.', ':']

        # 绘制各峰分量
        for idx, peak in enumerate(self.peaks_info):
            y_peak = peak['height'] * np.exp(
                -(self.x_fine - peak['position']) ** 2 / (2 * peak['sigma'] ** 2))
            ax.plot(self.x_fine, y_peak, linestyle=line_styles[idx % 3],
                    linewidth=1.8, alpha=0.85, label=f'峰{idx + 1}')
            total_fine += y_peak

        # 绘制总拟合曲线和原始数据
        ax.plot(self.x_fine, total_fine, 'r-', lw=2.5, label='总拟合')
        ax.plot(self.x, self.y, 'k-', alpha=0.7, lw=1.2, label='原始数据')

        # 绘制残差曲线（橙色虚线）
        ax.plot(self.x, self.residuals, '--', color='#FFA500',
                alpha=0.8, lw=1.5, label='残差')

        if self.stats is not None:
            self._add_stat_text(ax)

        ax.legend(loc='upper right', ncol=1, fontsize=16,
                  bbox_to_anchor=(0.995, 0.995), frameon=False,
                  handlelength=2.5)

    def _add_stat_text(self, ax):
        stat_text = (
            f"数据点(n): {self.stats['n']}\n"
            f"自由度: {self.stats['dof']}\n"
            f"调整后R²: {self.stats['adj_r2']:.4f}\n"
            f"SS: {self.stats['ss_res']:.4e}"
        )
        ax.text(0.02, 0.98, stat_text, transform=ax.transAxes,
                fontsize=14, verticalalignment='top', horizontalalignment='left')

    def _setup_result_table(self, ax):
        ax.axis('off')
        table_data = self._prepare_table_data()
        table = ax.table(cellText=table_data,
                         colLabels=['编号', '峰面积', 'FWHM', '高度', '位置', '面积占比'],
                         colColours=['#f8f8f8'] * 6,
                         cellLoc='center',
                         loc='center',
                         bbox=[0.05, 0.15, 0.9, 0.75])
        table.auto_set_font_size(False)
        table.set_fontsize(16)
        table.scale(1.1, 1.6)
        for (row, col), cell in table.get_celld().items():
            cell.set_edgecolor('#cccccc')
            if row == 0:
                cell.set_facecolor(Config.COLOR_THEME['table_header'])
                cell.set_text_props(color='white', fontproperties=FontProperties(weight='bold', size=12))
            elif row >= len(self.peaks_info) + 1:
                cell.set_facecolor(Config.COLOR_THEME['total_row'])
                if col == 4:
                    cell.set_text_props(fontproperties=FontProperties(weight='bold', size=12))
                if col == 5:
                    cell.set_text_props(color='#d62728', fontproperties=FontProperties(weight='bold', size=12))

    def _prepare_table_data(self):
        total_area = sum(peak['area'] for peak in self.peaks_info)
        table_data = []
        for idx, peak in enumerate(self.peaks_info, 1):
            table_data.append([
                f'峰{idx}',
                f"{peak['area']:.4f}",
                f"{2.355 * peak['sigma']:.2f}",
                f"{peak['height']:.4f}",
                f"{peak['position']:.2f}",
                f"{peak['area_percent']:.2f}%"
            ])
        total_percent = sum(p['area_percent'] for p in self.peaks_info)
        table_data.append(['', '', '', '', '合 计', f"{total_percent:.2f}%"])
        return table_data

    def save_plot(self, dpi=Config.SAVE_DPI):
        os.makedirs(self.save_dir, exist_ok=True)
        output_path = os.path.join(self.save_dir, f"{self.base_name}_fit.png")
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"图像已保存至: {output_path}")

    def export_data(self):
        """导出数据到Excel（分工作表存储不同分辨率数据）"""
        # 创建Excel写入器
        output_path = os.path.join(self.save_dir, f"{self.base_name}_fit_data.xlsx")
        with pd.ExcelWriter(output_path) as writer:
            # 工作表1: 高分辨率拟合曲线
            self._export_highres_curves(writer)

            # 工作表2: 原始数据及残差
            self._export_original_data(writer)

    def _export_highres_curves(self, writer):
        """导出高分辨率拟合曲线"""
        fit_data = {
            'x_fine': self.x_fine,
            '总拟合': np.zeros_like(self.x_fine)
        }

        # 收集各峰数据
        for idx, peak in enumerate(self.peaks_info, 1):
            y_peak = peak['height'] * np.exp(
                -(self.x_fine - peak['position']) ** 2 / (2 * peak['sigma'] ** 2))
            fit_data[f'峰{idx}'] = y_peak
            fit_data['总拟合'] += y_peak

        pd.DataFrame(fit_data).to_excel(
            writer, sheet_name='拟合曲线', index=False)

    def _export_original_data(self, writer):
        """导出原始数据及残差"""
        original_data = {
            '原始_x': self.x,
            '原始_y': self.y,
            '拟合值': self.y_fit_original_x,
            '残差': self.residuals
        }
        pd.DataFrame(original_data).to_excel(
            writer, sheet_name='原始数据', index=False)


# ================================
# 主程序
# ================================
if __name__ == "__main__":
    x, y_normalized, data_scale, file_path = load_raman_data()

    optimizer = PeakOptimizer(x, y_normalized)
    result = optimizer.optimize()
    params = result.x

    peaks_info = []
    sqrt_2pi = np.sqrt(2 * np.pi)
    n_fixed = len(Config.FIXED_PEAKS)

    # 处理固定峰
    for i in range(n_fixed):
        H = params[i] * data_scale
        mu = Config.FIXED_PEAKS[i]
        var = params[n_fixed + i]
        sigma = np.sqrt(var)
        area = H * sigma * sqrt_2pi
        peaks_info.append({
            'height': H, 'position': mu,
            'sigma': sigma, 'area': area,
            'area_percent': (area / (sum(p['area'] for p in peaks_info) + area)) * 100
        })

    # 处理可变峰
    var_pos = params[2 * n_fixed:2 * n_fixed + 2]
    for i in range(2):
        H = params[2 * n_fixed + 2 + i] * data_scale
        mu = var_pos[i]
        var = params[2 * n_fixed + 4 + i]
        sigma = np.sqrt(var)
        area = H * sigma * sqrt_2pi
        total_area = sum(p['area'] for p in peaks_info) + area
        for p in peaks_info:
            p['area_percent'] = p['area'] / total_area * 100
        peaks_info.append({
            'height': H, 'position': mu,
            'sigma': sigma, 'area': area,
            'area_percent': area / total_area * 100
        })

    # 计算拟合值和统计量
    y_fit_original_x = np.zeros_like(x)
    for peak in peaks_info:
        y_fit_original_x += peak['height'] * np.exp(
            -(x - peak['position']) ** 2 / (2 * peak['sigma'] ** 2))
    stats = calculate_statistics(x, y_normalized * data_scale, y_fit_original_x)

    # 可视化与导出
    visualizer = RamanVisualizer(
        peaks_info, x, y_normalized * data_scale, data_scale, file_path, y_fit_original_x
    )
    visualizer.create_figure(stats)

    if Config.EXPORT_DATA:
        visualizer.export_data()

    print("\n=== 拟合统计量 ===")
    print(f"数据点(n): {stats['n']}")
    print(f"自由度: {stats['dof']}")
    print(f"调整R²: {stats['adj_r2']:.4f}")
    print(f"SS残差: {stats['ss_res']:.4e}")