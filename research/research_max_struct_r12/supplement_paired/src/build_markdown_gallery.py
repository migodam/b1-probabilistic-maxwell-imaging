"""GitHub-rendered gallery independent of Pages/Actions availability."""
from pathlib import Path
root=Path(__file__).resolve().parents[1]
lines=['# R12 同起点配对补充：完整图集','',
'[中文分析](RESEARCH_REPORT_ZH.md) · [成本](COST_LEDGER.md) · [HTML图集](GALLERY.html) · [原始证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-paired-20261112)','',
'44条运行已封存评分。本图集由GitHub直接渲染，无需等待Pages构建。A为3个旧坏终点，B为完整24对象（4触发、20旁路），C为4个冻结NN对象；均为暴露对象上的后续配对，不是新盲确认。GT只用于封存后评分。','',
'材料图统一χ色标；相同编号尾号不代表同一个物理对象。TRIGGERED表示启动救援，不表示恢复成功。direct没有辅助proposal，其raw是共同初始材料。CPU横轴使用新增whole-child费用，B父运行费用另外保留。','']
for title,pattern in [('对象对应关系','cohort*'),('A：三个真实旧坏终点','A_*'),('B：同一批24对象的G19条件救援','B_*'),('C：冻结NN的接续接口','C_*'),('既有二维与三维corrector背景（非新实验）','original_*')]:
 lines+=['## '+title,'']
 for p in sorted((root/'FIGURES').glob(pattern+'.png')):
  label=p.stem.replace('_',' ');lines+=['### '+label,'',f'![{label}](FIGURES/{p.name})','',f'[原尺寸PNG](FIGURES/{p.name}) · [PDF](FIGURES/{p.stem}.pdf)','']
(root/'GALLERY.md').write_text('\n'.join(lines)+'\n');print('21-figure Markdown gallery generated')
