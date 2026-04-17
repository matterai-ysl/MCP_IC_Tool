"""
BaseHTMLGenerator — 所有 VASP HTML 报告生成器的公共基类。

提取了各 HTMLGenerator 中重复的：
- CSS 基础样式
- HTML shell（doctype / head / body）
- 页面头部、页脚生成
- generate_html_report 文件写入逻辑
"""

import base64
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional


class BaseHTMLGenerator:
    """
    VASP HTML 报告生成器基类。

    子类需覆写:
        - report_title: str          — 页面 <title>
        - report_emoji: str          — 标题前的 emoji
        - report_heading: str        — <h1> 标题文本
        - report_attribution: str    — 页脚来源说明
        - header_gradient: str       — header 背景渐变 CSS
        - _generate_body_sections()  — 返回各专属段落 HTML
        - _generate_javascript()     — 返回 Chart.js 或其它 JS 代码
    """

    # ---------- 子类必须/可覆写的属性 ----------
    report_title: str = "VASP分析报告"
    report_emoji: str = "📊"
    report_heading: str = "VASP分析报告"
    report_attribution: str = "由VASP API可视化分析模块生成"
    header_gradient: str = "linear-gradient(135deg, #667eea 0%, #764ba2 100%)"
    use_chartjs: bool = True

    def __init__(self, analysis_data: Dict[str, Any]):
        self.data = analysis_data

    # ------------------------------------------------------------------ #
    #  公共入口
    # ------------------------------------------------------------------ #
    def generate_html_report(self, output_path: str) -> str:
        """生成 HTML 报告并写入文件。"""
        print(f"📄 正在生成HTML报告: {output_path}")

        html_content = self._generate_html_content()

        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"✅ HTML报告已生成: {output_file}")
        return str(output_file)

    # ------------------------------------------------------------------ #
    #  HTML 整体框架
    # ------------------------------------------------------------------ #
    def _generate_html_content(self) -> str:
        chartjs_tag = (
            '<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>'
            if self.use_chartjs else ''
        )
        return f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.report_title}</title>
    {chartjs_tag}
    <style>
        {self._get_css_styles()}
    </style>
</head>
<body>
    <div class="container">
        {self._generate_header()}
        {self._generate_body_sections()}
        {self._generate_footer()}
    </div>

    <script>
        {self._generate_javascript()}
    </script>
</body>
</html>
"""

    # ------------------------------------------------------------------ #
    #  共享 CSS
    # ------------------------------------------------------------------ #
    def _get_css_styles(self) -> str:
        """基础 CSS 样式，子类可通过 _get_extra_css_styles() 追加。"""
        base = f"""
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f5f5f5;
        }}

        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}

        .header {{
            background: {self.header_gradient};
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .section {{
            background: white;
            padding: 25px;
            margin-bottom: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .section h2 {{
            color: #4a5568;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #e2e8f0;
            font-size: 1.5em;
        }}

        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}

        .summary-card {{
            background: #f7fafc;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #4299e1;
        }}

        .summary-card h3 {{
            color: #2d3748;
            margin-bottom: 10px;
        }}

        .summary-card .value {{
            font-size: 1.5em;
            font-weight: bold;
            color: #4299e1;
        }}

        .convergence-status {{
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
            color: white;
        }}

        .converged {{
            background-color: #48bb78;
        }}

        .not-converged {{
            background-color: #f56565;
        }}

        .chart-container {{
            position: relative;
            height: 400px;
            margin: 20px 0;
        }}

        .data-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}

        .data-table th,
        .data-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e2e8f0;
        }}

        .data-table th {{
            background-color: #edf2f7;
            font-weight: bold;
            color: #4a5568;
        }}

        .data-table tr:hover {{
            background-color: #f7fafc;
        }}

        .footer {{
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #718096;
            font-size: 0.9em;
        }}

        .grid-2 {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }}

        .grid-3 {{
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 20px;
        }}

        .download-bar {{
            text-align: center;
            margin: 6px 0 12px 0;
        }}

        .download-link {{
            display: inline-block;
            padding: 4px 14px;
            font-size: 0.82em;
            color: #4299e1;
            border: 1px solid #bee3f8;
            border-radius: 20px;
            background: #ebf8ff;
            text-decoration: none;
            transition: background 0.2s, color 0.2s;
        }}

        .download-link:hover {{
            background: #4299e1;
            color: white;
        }}

        @media (max-width: 768px) {{
            .grid-2, .grid-3 {{
                grid-template-columns: 1fr;
            }}

            .summary-grid {{
                grid-template-columns: 1fr;
            }}
        }}
        """
        extra = self._get_extra_css_styles()
        if extra:
            base += "\n" + extra
        return base

    def _get_extra_css_styles(self) -> str:
        """子类可覆写以追加额外 CSS。"""
        return ""

    # ------------------------------------------------------------------ #
    #  页面头部 / 页脚
    # ------------------------------------------------------------------ #
    def _generate_header(self) -> str:
        file_info = self.data.get('file_info', {})
        task_info = self.data.get('task_info', {})
        return f"""
        <div class="header">
            <h1>{self.report_emoji} {self.report_heading}</h1>
            <p>任务ID: <strong>{task_info.get('task_id', 'Unknown')}</strong> |
               材料组成: <strong>{task_info.get('composition', 'Unknown')}</strong></p>
            <p>生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>VASP版本: {file_info.get('vasp_version', 'Unknown')} |
               计算日期: {file_info.get('calculation_date', 'Unknown')} |
               核数: {file_info.get('total_cores', 'Unknown')}</p>
        </div>
        """

    def _generate_footer(self) -> str:
        return f"""
        <div class="footer">
            <p>{self.report_emoji} {self.report_heading} | 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>{self.report_attribution}</p>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  工具方法：生成 CSV 下载链接
    # ------------------------------------------------------------------ #
    def _csv_download_link(
        self,
        filename: str,
        headers: List[str],
        rows: List[List],
        label: str = "📥 下载原始数据 (CSV)",
    ) -> str:
        """将给定的表头+数据行编码为 base64 CSV，返回一个可直接内嵌的下载链接 HTML。"""
        lines = [','.join(headers)]
        for row in rows:
            lines.append(','.join('' if v is None else str(v) for v in row))
        csv_bytes = '\n'.join(lines).encode('utf-8')
        csv_b64 = base64.b64encode(csv_bytes).decode('ascii')
        return (
            f'<div class="download-bar">'
            f'<a href="data:text/csv;base64,{csv_b64}" download="{filename}" class="download-link">{label}</a>'
            f'</div>'
        )

    def _binary_download_link(
        self,
        filename: str,
        media_type: str,
        base64_payload: str,
        label: str,
    ) -> str:
        if not base64_payload:
            return ""
        return (
            f'<div class="download-bar">'
            f'<a href="data:{media_type};base64,{base64_payload}" download="{filename}" class="download-link">{label}</a>'
            f'</div>'
        )

    def _png_download_link(
        self,
        filename: str,
        base64_payload: str,
        label: str = "🖼️ 下载图片 (PNG)",
    ) -> str:
        return self._binary_download_link(filename, "image/png", base64_payload, label)

    # ------------------------------------------------------------------ #
    #  子类必须覆写
    # ------------------------------------------------------------------ #
    def _generate_body_sections(self) -> str:
        """子类覆写，返回所有内容段落的 HTML 拼接。"""
        return ""

    def _generate_javascript(self) -> str:
        """子类覆写，返回 Chart.js 或其它 JavaScript 代码。"""
        return ""
