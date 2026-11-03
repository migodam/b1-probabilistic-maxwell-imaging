"""Add R3 entry only; frozen R1/R2 directories are not modified."""
from pathlib import Path
import re
R=Path(__file__).resolve().parents[3]
base='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r3/'
raw='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r3/'
card='''<!-- MAX_STRUCT_R3_PUBLICATION -->
## MAX-STRUCT-R3：公式驱动的机制研究循环

噪声硬锁代价、被动单元圆弧与更新秩层级；真实边界停滞和3D适用性控制。11个同一公开对象的开发控制未建立新增成像优势，全部失败、成本和后续方向保留。

'''+f'[中文网页]({base}) · [给Chat的入口]({raw}CHAT_HANDOFF.md) · [数学推导]({raw}THEORY_DERIVATIONS.md)\n<!-- /MAX_STRUCT_R3_PUBLICATION -->\n\n'
for name in ['README.md','START_HERE.md','research/README.md']:
 p=R/name;s=p.read_text()
 if '<!-- MAX_STRUCT_R3_PUBLICATION -->' not in s:p.write_text(card+s)
cardhtml=f'<section style="margin:24px 0;padding:24px;border:1px solid #bdd0dd;border-radius:10px"><h2>MAX-STRUCT-R3：公式驱动的机制研究循环</h2><p>噪声硬锁、被动响应圆弧与更新秩的物理边界；未建立新增成像优势，完整正负证据保留。</p><p><a href="{base}">中文完整报告</a> · <a href="{raw}CHAT_HANDOFF.md">Chat读取入口</a></p></section>'
for name in ['index.html','research/index.html']:
 p=R/name;s=p.read_text()
 if 'MAX-STRUCT-R3：公式驱动的机制研究循环' not in s:p.write_text(re.sub(r'(<main\b[^>]*>)',lambda m:m[0]+cardhtml,s,count=1))
