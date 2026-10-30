"""Isolated current-atlas rendering; never calls a historical builder main."""
from pathlib import Path
import html
import importlib.util
import json
import os
import re

ROOT = Path(__file__).resolve().parents[3]
AT = ROOT / 'research/maxwell_formula_atlas'
P = ROOT / 'research/p3l_r12_bilateral_fixedpoint'
spec = importlib.util.spec_from_file_location(
    'r12_atlas_document_reader',
    ROOT / 'research/p3l_r9_stabilized_subspace/src/build_reading_site.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


class Reader(base.Reader):
    def slug(self, title):
        match = re.search(r'F-(?:BF|SS|BI|BASE|OP|OPT|ERR|MAT|PROB|VOL|POST)-\d+', title)
        return match[0].lower() if match else super().slug(title)

    def shell(self, title, body, page, subtitle='物理 → 原先验与材料几何 → 完整导数 → 有限证据'):
        links = [(AT / 'index.html', '公式目录'),
                 (AT / '12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.html', 'R12新增九条'),
                 (AT / 'FORMULA_ATLAS.html', '完整公式总集'),
                 (AT / 'FORMULA_ATLAS.md', '原始Markdown'),
                 (P / 'index.html', 'R12研究报告')]
        nav = '<nav>' + ''.join('<a href="' + os.path.relpath(path, page.parent) + '">' + label + '</a>'
                                for path, label in links) + '</nav>'
        style = base.old.STYLE + '''
        .math,mjx-container{max-width:100%;overflow-x:auto;overflow-y:hidden}
        .formula-list{columns:2}.formula-list a{display:block;break-inside:avoid;padding:6px}
        .formula-list a[hidden]{display:none}#formula-filter{width:100%;padding:12px;font-size:17px}
        @media(max-width:800px){.formula-list{columns:1}}
        '''
        return ('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1"><title>'
                + html.escape(title) + '</title><style>' + style + '</style>' + base.MATH
                + '</head><body><header>' + nav + '<h1>' + html.escape(title) + '</h1><p>'
                + html.escape(subtitle) + '</p></header><main>' + body
                + '</main><footer style="padding:30px;max-width:1100px;margin:auto">'
                '稳定编号；数学推导、实现、有限实验与新颖性分开。'
                '原先验、完整材料RMS及历史失败记录按各章来源保留。'
                '</footer></body></html>')


reader = Reader(AT)
documents = [AT / '12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.md',
             AT / 'FORMULA_ATLAS.md', P / 'agent_a_max/CONTINUATION_NOTES.md']
reader.render(documents)
registry = json.loads((AT / 'FORMULA_REGISTRY.json').read_text())
entries = ''.join('<a data-text="' + html.escape((entry['id'] + ' ' + entry['title']).lower(), quote=True)
                  + '" href="' + entry['chapter'].replace('.md', '.html') + '#' + entry['anchor']
                  + '"><b>' + entry['id'] + '</b> ' + html.escape(entry['title']) + '</a>'
                  for entry in registry['formulas'])
body = ('<section class="panel"><p>当前 <b>' + str(registry['count'])
        + '</b> 条可追溯公式。R12新增双边材料模态、原先验增广、完整固定点导数与有限经验诊断九条；原有158条编号保持稳定。</p>'
        '<p><a href="12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.html">第12章：双边材料模态与固定点</a> · '
        '<a href="FORMULA_ATLAS.html">连续阅读全文</a> · <a href="FORMULA_ATLAS.md">完整Markdown</a> · '
        '<a href="FORMULA_REGISTRY.json">JSON索引</a> · <a href="SOURCE_FORMULA_INDEX.csv">源码追溯CSV</a></p>'
        '<p>原先验支持、材料可能性、经验覆盖和后验认证分别声明。历史各轮费用与失败状态保留在对应章节。</p>'
        '<input id="formula-filter" type="search" placeholder="搜索公式编号 / 双边 / 先验 / 固定点 / 材料" aria-label="搜索公式">'
        '<div class="formula-list">' + entries + '</div></section>'
        '<script>document.getElementById("formula-filter").addEventListener("input",e=>{'
        'const q=e.target.value.toLowerCase();document.querySelectorAll(".formula-list a").forEach(a=>a.hidden=!a.dataset.text.includes(q));});</script>')
(AT / 'index.html').write_text(reader.shell('Maxwell成像：持续维护的公式集', body, AT / 'index.html'))
print(json.dumps(dict(rendered_documents=len(documents), formulas=registry['count'],
                      historical_builder_main_called=False, physics_calls=0)))
