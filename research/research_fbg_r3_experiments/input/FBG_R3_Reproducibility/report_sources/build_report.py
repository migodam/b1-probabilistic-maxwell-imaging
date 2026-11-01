"""Build the Chinese mathematical report and its portable reproduction package.

The mathematical body uses XeLaTeX. The designed cover is generated with
ReportLab and open Noto glyph outlines, then combined using pypdf.
No literature full texts are distributed in the reproduction archive.
"""
from pathlib import Path
import csv
import json
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'output'
OUT.mkdir(exist_ok=True)


def escape_cell(s):
    s = s.replace('|', '\\|').replace('*', '\\*').replace('_', '\\_')
    # Very long DOI tokens otherwise exceed narrow table columns.
    s = re.sub(r'10\.\d{4,9}/[A-Za-z0-9.()/_-]+',
               lambda m: '\\url{' + m.group(0) + '}', s)
    return s


def build_matrix():
    with (HERE/'ALGORITHM_COMPARISON.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f))
    assert len(rows) == 19 and all(len(r) == 14 for r in rows)
    panels = [
        ('A：变量、电流、非线性与目标', [0, 1, 2, 3, 4]),
        ('B：子空间、照明、搜索与数学成熟度', [0, 5, 6, 7, 8]),
        ('C：前人、增量、风险、成本与优先级', [0, 9, 10, 11, 12, 13]),
    ]
    text = '\n```{=latex}\n\\clearpage\n\\begin{landscape}\n\\renewcommand{\\ReportTableSize}{\\small}\n```\n\n# 附录：18 条路线的完整 14 列比较矩阵\n\n'
    text += '三个面板使用相同的算法名称逐行对应；同名行合并后即为一张完整 14 列矩阵。另附 CSV 保留全部字段，便于筛选与逐格比较。\n\n'
    for k, (title, cols) in enumerate(panels):
        if k:
            text += '\n\\clearpage\n\n'
        text += '## 面板 ' + title + '\n\n'
        text += '|' + '|'.join(escape_cell(rows[0][i]) for i in cols) + '|\n'
        text += '|' + '|'.join('---' for _ in cols) + '|\n'
        for row in rows[1:]:
            text += '|' + '|'.join(escape_cell(row[i]) for i in cols) + '|\n'
    text += '\n```{=latex}\n\\end{landscape}\n```\n'
    (HERE/'matrix_appendix.md').write_text(text, encoding='utf-8')
    shutil.copy2(HERE/'ALGORITHM_COMPARISON.csv', OUT/'FBG_R3_Algorithm_Comparison.csv')
    return text


def cover():
    from reportlab.pdfgen import canvas
    from reportlab.lib.colors import HexColor, Color
    from fontTools.ttLib import TTFont
    from fontTools.pens.basePen import BasePen

    fonts = {}
    for weight in ('Regular', 'Bold'):
        font = TTFont(HERE/'fonts'/f'NotoSansCJKsc-{weight}.otf')
        fonts[weight] = (font, font.getGlyphSet(), font.getBestCmap())
    target = HERE/'cover.pdf'
    c = canvas.Canvas(str(target), pagesize=(595.276, 841.89))
    c.setTitle('Research FBG-R3-T')

    class CanvasPen(BasePen):
        def __init__(self, gs, px, py, scale):
            super().__init__(gs)
            self.p = c.beginPath()
            self.px, self.py, self.sc = px, py, scale
        def point(self, a):
            return (self.px+a[0]*self.sc, self.py+a[1]*self.sc)
        def _moveTo(self, a): self.p.moveTo(*self.point(a))
        def _lineTo(self, a): self.p.lineTo(*self.point(a))
        def _curveToOne(self, a, b, d):
            self.p.curveTo(*self.point(a), *self.point(b), *self.point(d))
        def _closePath(self): self.p.close()
        def _endPath(self): pass

    def zh(s, x, y, size, color='#183B40', bold=False):
        font, gs, cmap = fonts['Bold' if bold else 'Regular']
        scale = size/font['head'].unitsPerEm
        c.setFillColor(HexColor(color))
        for char in s:
            name = cmap.get(ord(char), '.notdef')
            pen = CanvasPen(gs, x, y, scale)
            gs[name].draw(pen)
            c.drawPath(pen.p, stroke=0, fill=1)
            x += font['hmtx'][name][0]*scale

    c.setFillColor(HexColor('#FAFBF8'))
    c.rect(0, 0, 595.276, 841.89, fill=1, stroke=0)
    c.setFillColor(HexColor('#124E55'))
    c.rect(0, 584, 595.276, 258, fill=1, stroke=0)
    c.setFillColor(HexColor('#96C3BB'))
    c.rect(44, 779, 44, 4, fill=1, stroke=0)
    c.setFillColor(HexColor('#FFFFFF'))
    c.setFont('Helvetica-Bold', 29)
    c.drawString(44, 731, 'RESEARCH FBG-R3-T')
    c.setFont('Helvetica', 11)
    c.drawString(45, 708, 'MAXWELL FEEDBACK GEOMETRY AND NONLINEAR INVERSION')
    zh('Maxwell 反馈几何', 44, 657, 25, '#FFFFFF', True)
    zh('辅助电流与多分支材料反演', 44, 621, 25, '#FFFFFF', True)

    zh('研究判断', 44, 542, 12, '#52716F', True)
    zh('以共享材料的可实现性几何', 44, 507, 19, bold=True)
    zh('指导电流选择、有限释放与物理恢复', 44, 476, 19, bold=True)
    zh('保留成熟 FFT-TSOM，独立检验小支持代数多根提案。', 44, 442, 11.5, '#52716F')
    c.setStrokeColor(HexColor('#CEDDDA'))
    c.line(44, 413, 551, 413)

    for number, title, desc, y in [
        ('01', '有限材料读出与 0 / 2 / 4 阶次', '完整公式、条件、证明及反例', 369),
        ('02', '材料任务相关的辅助电流释放', '两通道材料步、软释放与同空间负控制', 298),
        ('03', '精确多根提案与全材料校正', '有限响应、合法分支、真实成本与决定性对照', 227),
    ]:
        c.setFillColor(HexColor('#124E55'))
        c.setFont('Helvetica-Bold', 17)
        c.drawString(44, y, number)
        zh(title, 87, y, 14, bold=True)
        zh(desc, 87, y-23, 10.5, '#52716F')

    c.setFillColor(HexColor('#EAF2EF'))
    c.roundRect(44, 100, 507, 65, 5, fill=1, stroke=0)
    zh('理论推导与小型代数核验已完成', 59, 140, 11.5, bold=True)
    zh('未运行新的 576 参数成像基准；后续实验协议已附。', 59, 119, 10.5, '#52716F')
    zh('完整中文研究报告 · 文献与来源审计 · 算法比较矩阵', 44, 68, 9.5, '#52716F')
    c.setFillColor(HexColor('#52716F'))
    c.setFont('Helvetica', 9.5)
    c.drawRightString(551, 42, '31 OCTOBER 2026')
    c.save()
    return target


def main():
    cover()
    body = (HERE/'report_body.md').read_text(encoding='utf-8')
    for dash in ('\u2010','\u2011','\u2012','\u2013','\u2014'):
        body = body.replace(dash, '-')
    assert '\r' not in body
    complete = HERE/'report_complete.md'
    complete.write_text(body+build_matrix(), encoding='utf-8')
    subprocess.run([
        'pandoc', complete.name, '--from', 'markdown+tex_math_dollars', '--to', 'latex',
        '--standalone', '--template=report_template.tex', '--no-highlight',
        '-o', 'report_complete.tex'], cwd=HERE, check=True)
    for k in range(2):
        with (HERE/f'build_pass{k+1}.log').open('w') as log:
            subprocess.run(['xelatex','-interaction=nonstopmode','-halt-on-error',
                            'report_complete.tex'], cwd=HERE, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    from pypdf import PdfReader
    final = OUT/'Research_FBG_R3_T_Report_ZH.pdf'
    shutil.copy2(HERE/'report_complete.pdf', final)
    print(json.dumps({'pdf':str(final),'pages':len(PdfReader(final).pages),
                      'matrix_rows':18,'matrix_columns':14},ensure_ascii=False))


if __name__ == '__main__':
    main()
