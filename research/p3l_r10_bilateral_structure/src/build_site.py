"""Render delivered Markdown and append the living formula atlas; zero physics."""
from pathlib import Path
import importlib.util, json, re, html, os, shutil
P=Path(__file__).resolve().parents[1];R=P.parents[1];AT=R/'research/maxwell_formula_atlas'
spec=importlib.util.spec_from_file_location('r9_reader',R/'research/p3l_r9_stabilized_subspace/src/build_reading_site.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
class Reader(base.Reader):
 def slug(self,t):
  m=re.search(r'F-(?:BI|BASE|OP|OPT|ERR|MAT|PROB|VOL|POST)-\d+',t)
  return m[0].lower() if m else super().slug(t)
 def url(self,u,source,page,image=False):
  if '/b1-p3l-r10-bilateral-structure/' in u:
   tail=u.split('/b1-p3l-r10-bilateral-structure/',1)[1]
   return 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/'+tail
  return super().url(u,source,page,image)
 def shell(self,title,body,page,subtitle='固定结构 → 可复用响应 → 同任务成本 → 材料覆盖的真实边界'):
  links=[(P/'index.html','R10研究'),(AT/'index.html','持续公式集'),(P/'GO_NO_GO.html','科学决策'),(P/'RESEARCH_REPORT_ZH.md','原始Markdown'),(P/'PAPER_NARRATIVE_V2_EN.html','英文论文草稿')]
  nav='<nav>'+''.join('<a href="'+os.path.relpath(f,page.parent)+'">'+t+'</a>' for f,t in links)+'</nav>'
  return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+base.old.STYLE+' .math,mjx-container{max-width:100%;overflow-x:auto;overflow-y:hidden}table{display:block;overflow-x:auto}#formula-filter{width:100%;padding:12px;font-size:17px;box-sizing:border-box}.formula-list{columns:2}.formula-list a{display:block;break-inside:avoid;padding:5px}.formula-list a[hidden]{display:none}@media(max-width:800px){.formula-list{columns:1}} </style>'+base.MATH+'</head><body><header>'+nav+'<h1>'+html.escape(title)+'</h1><p>'+html.escape(subtitle)+'</p></header><main>'+body+'</main><footer style="padding:30px;max-width:1100px;margin:auto">原P1完整非线性Maxwell。物理核费用、材料可能性与条件概率分别评价。新结果不改写历史结论。</footer></body></html>'

def build():
 chapters=[AT/'09_BILATERAL_SHAPE_AMPLITUDE.md',AT/'10_TSOM_ON_BILATERAL_TRANSFER.md']
 registry=json.loads((AT/'FORMULA_REGISTRY.json').read_text())
 existing={x['id']:x for x in registry['formulas']}
 for chapter in chapters:
  for m in re.finditer(r'^##\s+(F-BI-\d+)\s*·\s*(.+)$',chapter.read_text(),re.M):
   existing[m[1]]=dict(id=m[1],title=m[2],chapter=chapter.name,anchor=m[1].lower(),source='R10 derived formula chapter; precise assumptions, standard theory and executed evidence separated')
 registry.update(version='2026-10-30 R10',count=len(existing),formulas=list(existing.values()))
 (AT/'FORMULA_REGISTRY.json').write_text(json.dumps(registry,ensure_ascii=False,indent=2)+'\n')
 combined=AT/'FORMULA_ATLAS.md';s=combined.read_text();marker='# 09 · 双边传递'
 if marker in s:s=s.split(marker)[0].rstrip().removesuffix('---').rstrip()
 s=s.replace('P1—R9主链','P1—R10主链')+'\n\n---\n\n'+'\n\n---\n\n'.join(c.read_text() for c in chapters)
 combined.write_text(s)
 log=AT/'CHANGELOG.md';s=log.read_text()
 if 'F-BI-001' not in s:log.write_text(s+'\n2026-10-30 R10：追加F-BI-001—018，双边Gram、形状幅度有理传递、随机目标条件协方差、损耗稳定性与residual任务界。原122条编号与历史实验保留；总计140条。明确共同不变空间与对角化条件，数值费用不作为新恒等式证明。\n')
 rd=Reader(P)
 docs=list((P/'followup_shape_reuse').glob('*.md'))+list(P.glob('*.md'))+list((P/'GPT_PRO_HANDOFF').glob('*.md'))+list(AT.glob('[0-9][0-9]_*.md'))+[combined]
 rd.render(docs)
 shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'index.html');shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'report.html')
 entries=''.join('<a data-text="'+html.escape((x['id']+' '+x['title']).lower(),quote=True)+'" href="'+x['chapter'].replace('.md','.html')+'#'+x['anchor']+'"><b>'+x['id']+'</b> '+html.escape(x['title'])+'</a>' for x in registry['formulas'])
 body='<section class="panel"><p>当前'+str(registry['count'])+'条可追溯推导。标准恒等式与实验增量分开；TSOM追问新增6条解析推导，零新增物理计算。</p><p><a href="10_TSOM_ON_BILATERAL_TRANSFER.html">TSOM怎样作用于中间逆 →</a> · <a href="09_BILATERAL_SHAPE_AMPLITUDE.html">R10双边结构18条</a> · <a href="FORMULA_ATLAS.html">连续阅读全文</a> · <a href="FORMULA_ATLAS.md">原始Markdown</a> · <a href="FORMULA_REGISTRY.json">JSON索引</a></p><input id="formula-filter" type="search" placeholder="检索 Born / Gram / Schur / 概率 / 形状 / 公式编号" aria-label="检索公式"><div class="formula-list">'+entries+'</div></section><script>document.getElementById("formula-filter").addEventListener("input",e=>{const q=e.target.value.toLowerCase();document.querySelectorAll(".formula-list a").forEach(a=>a.hidden=!a.dataset.text.includes(q));});</script>'
 (AT/'index.html').write_text(rd.shell('Maxwell成像：持续维护的公式集',body,AT/'index.html'))
 print(json.dumps(dict(formulas=registry['count'],documents=len(docs))))
if __name__=='__main__':build()
