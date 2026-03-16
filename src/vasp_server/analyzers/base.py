"""
BaseAnalyzer — 所有 VASP 分析器的公共基类。

提取了 optimization_analyzer 和 scf_analyzer 中重复的：
- 工作目录/OUTCAR 路径解析
- POSCAR/CONTCAR 结构文件读取
- OUTCAR 文件信息解析
- OUTCAR 计算设置解析
"""

import re
import os
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime

from .file_utils import (
    resolve_work_dir_and_outcar,
    find_vasp_file,
    extract_composition_from_poscar,
    read_structure_files,
)


class BaseAnalyzer:
    """
    VASP 分析器基类。

    子类只需覆写:
        - _init_data(): 返回初始 data dict
        - _do_analyze(): 执行具体的分析步骤
    """

    def __init__(self, input_path: str, task_id: Optional[str] = None):
        self.task_id = task_id or "unknown"
        self.work_dir, self.outcar_path = resolve_work_dir_and_outcar(input_path)
        self.poscar_path = find_vasp_file(self.work_dir, "POSCAR")
        self.contcar_path = find_vasp_file(self.work_dir, "CONTCAR")
        self.data: Dict[str, Any] = self._init_data()

    def _init_data(self) -> Dict[str, Any]:
        """子类覆写以提供自己的 data 结构。"""
        return {
            'file_info': {},
            'calculation_settings': {},
            'final_results': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {},
        }

    def analyze(self) -> Dict[str, Any]:
        """执行完整分析，返回 data dict。"""
        self._parse_structure_files()
        self._parse_outcar_common()
        self._do_analyze()
        return self.data

    def _do_analyze(self):
        """子类覆写以执行具体分析。"""
        pass

    # ------------------------------------------------------------------ #
    #  结构文件解析 (共享)
    # ------------------------------------------------------------------ #
    def _parse_structure_files(self):
        """解析 POSCAR 和 CONTCAR。"""
        files = read_structure_files(self.work_dir)
        self.data['structure_files'] = files

        poscar_content = files.get('poscar', '')
        if poscar_content:
            composition = extract_composition_from_poscar(poscar_content)
            self.data['task_info']['composition'] = composition
        else:
            self.data['task_info']['composition'] = "Unknown"

    # ------------------------------------------------------------------ #
    #  OUTCAR 公共解析
    # ------------------------------------------------------------------ #
    def _parse_outcar_common(self):
        """解析 OUTCAR 的通用部分: file_info 和 calculation_settings。"""
        if not self.outcar_path or not self.outcar_path.exists():
            return

        try:
            with open(self.outcar_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
        except Exception:
            return

        self._parse_file_info(content)
        self._parse_calculation_settings(content)
        # 子类在 _do_analyze 中调用自己的专用解析
        self._outcar_content = content

    def _parse_file_info(self, content: str):
        """解析 OUTCAR 的文件信息部分。"""
        info: Dict[str, Any] = {}

        # 文件大小
        try:
            info['file_size'] = os.path.getsize(self.outcar_path)
            info['file_size_mb'] = round(info['file_size'] / (1024 * 1024), 2)
        except Exception:
            pass

        # VASP 版本
        m = re.search(r'vasp\.\s*(\S+)', content[:2000])
        if m:
            info['vasp_version'] = m.group(1)

        # 运行日期
        m = re.search(r'executed on.*date\s+(\S+)', content[:2000])
        if m:
            info['run_date'] = m.group(1)

        # 总行数 (近似)
        info['total_lines'] = content.count('\n')

        self.data['file_info'] = info

    def _parse_calculation_settings(self, content: str):
        """解析 OUTCAR 中的 INCAR 计算设置。"""
        settings: Dict[str, Any] = {}

        # 常见 INCAR 参数的正则
        param_patterns = {
            'ENCUT': r'ENCUT\s*=\s*([\d.]+)',
            'EDIFF': r'EDIFF\s*=\s*([\d.E+-]+)',
            'EDIFFG': r'EDIFFG\s*=\s*([\d.E+-]+)',
            'ISMEAR': r'ISMEAR\s*=\s*([-\d]+)',
            'SIGMA': r'SIGMA\s*=\s*([\d.]+)',
            'ISPIN': r'ISPIN\s*=\s*(\d+)',
            'NSW': r'NSW\s*=\s*(\d+)',
            'IBRION': r'IBRION\s*=\s*([-\d]+)',
            'ISIF': r'ISIF\s*=\s*(\d+)',
            'PREC': r'PREC\s*=\s*(\w+)',
            'ALGO': r'ALGO\s*=\s*(\w+)',
            'LREAL': r'LREAL\s*=\s*(\S+)',
            'GGA': r'GGA\s*=\s*(\w+)',
            'NELM': r'NELM\s*=\s*(\d+)',
            'POTIM': r'POTIM\s*=\s*([\d.]+)',
            'IVDW': r'IVDW\s*=\s*(\d+)',
        }

        # 只搜索前 10000 字符（设置通常在文件头部）
        header = content[:10000]
        for param, pattern in param_patterns.items():
            m = re.search(pattern, header, re.IGNORECASE)
            if m:
                settings[param] = m.group(1)

        # K 点网格
        m = re.search(r'k-points\s+NKPTS\s*=\s*(\d+)', header)
        if m:
            settings['NKPTS'] = int(m.group(1))

        # 原子数
        m = re.search(r'NIONS\s*=\s*(\d+)', header)
        if m:
            settings['NIONS'] = int(m.group(1))

        self.data['calculation_settings'] = settings

    # ------------------------------------------------------------------ #
    #  常用 OUTCAR 提取方法 (子类可复用)
    # ------------------------------------------------------------------ #
    def _get_outcar_content(self) -> str:
        """获取已读取的 OUTCAR 内容。"""
        return getattr(self, '_outcar_content', '')
