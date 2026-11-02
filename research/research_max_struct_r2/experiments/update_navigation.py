"""Add R2 navigation without changing earlier reports."""
from pathlib import Path
import re
R=Path(__file__).resolve().parents[3]
base='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r2/'
raw='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r2/'
card='''<!-- MAX_STRUCT_R2_PUBLICATION -->
## MAX-STRUCT-R2：结构审计、新盲测与新增照明

三新对象四方法：成熟CIE/global均2/3，精确profile1/3，未建立超越成熟CIE的新成像机制。满容量事后控制仍失败；三维新增照明只支持公开一参数族中的相容估计。代码、原始结果、独立审查与全部费用保留。

'''+f'[中文网页]({base}) · [给Chat的入口]({raw}CHAT_HANDOFF.md) · [完整代码/证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_max_struct_r2)\n<!-- /MAX_STRUCT_R2_PUBLICATION -->\n\n'
for name in ['README.md','START_HERE.md','research/README.md']:
    p=R/name;text=p.read_text()
    if '<!-- MAX_STRUCT_R2_PUBLICATION -->' not in text:p.write_text(card+text)
html=f'<section style="margin:24px 0;padding:24px;border:1px solid #bdd0dd;border-radius:10px"><h2>MAX-STRUCT-R2：结构审计与新盲测</h2><p>成熟CIE/global均2/3；未建立新增成像优势。三维新增照明与强公开prior单独归因。</p><p><a href="{base}">中文完整报告</a> · <a href="{raw}CHAT_HANDOFF.md">Chat读取入口</a></p></section>'
for name in ['index.html','research/index.html']:
    p=R/name;text=p.read_text()
    if 'MAX-STRUCT-R2：结构审计与新盲测' not in text:
        text=re.sub(r'(<main\b[^>]*>)',lambda m:m[0]+html,text,count=1)
        p.write_text(text)
