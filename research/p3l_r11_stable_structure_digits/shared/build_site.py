"""R11 documents/formula rendering only; no experiment or physics invocation."""
from pathlib import Path
import importlib.util,html,re,json,shutil,os
P=Path(__file__).resolve().parents[1];R=P.parents[1];AT=R/'research/maxwell_formula_atlas'
spec=importlib.util.spec_from_file_location('r11_base_reader',R/'research/p3l_r9_stabilized_subspace/src/build_reading_site.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
class Reader(base.Reader):
 def slug(self,t):
  m=re.search(r'F-(?:SS|BI|BASE|OP|OPT|ERR|MAT|PROB|VOL|POST)-\d+',t)
  return m[0].lower() if m else super().slug(t)
 def shell(self,title,body,page,subtitle='严格支撑 · 固定总场 · 条件耗散界 · 数字1–9真实逆成像'):
  links=[(P/'index.html','R11报告'),(P/'DIGITS_VISUAL_REVIEW_ZH.html','逐数字视觉审查'),(P/'GO_NO_GO.html','科学裁决'),(P/'RESEARCH_REPORT_ZH.md','原始Markdown'),(AT/'index.html','158条公式')]
  nav='<nav>'+''.join('<a href="'+os.path.relpath(a,page.parent)+'">'+t+'</a>' for a,t in links)+'<a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr11-stable-digits-20261030">完整证据下载</a></nav>'
  warn='<div class="notice"><strong>R11资源协议 FAIL_CAP_EXCEEDED：</strong>全阶LU 28,727 / 12,000。数值与图像证据保留；预算未改写，新增物理已停止。</div>'
  return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+base.old.STYLE+' .math,mjx-container{max-width:100%;overflow-x:auto;overflow-y:hidden}.formula-list{columns:2}.formula-list a{display:block;break-inside:avoid;padding:5px}.formula-list a[hidden]{display:none}#formula-filter{width:100%;padding:12px;font-size:17px}@media(max-width:800px){.formula-list{columns:1}}</style>'+base.MATH+'</head><body><header>'+nav+'<h1>'+html.escape(title)+'</h1><p>'+html.escape(subtitle)+'</p></header><main>'+warn+body+'</main><footer style="padding:30px;max-width:1100px;margin:auto">N/P材料模型、物理核与整成像收益分别评价。标准恒等式、数值条件界与实测结果分开。未认证Full66 posterior。</footer></body></html>'

def main():
 rd=Reader(P);docs=list(P.glob('*.md'))+list((P/'GPT_PRO_HANDOFF').glob('*.md'))
 for d in ['agent_a','agent_b','agent_c']:docs+=list((P/d).glob('*.md'))
 docs+=[AT/'11_ACTIVE_SUPPORT_TOTAL_FIELD_AND_GLOBAL_BOUNDS.md',AT/'FORMULA_ATLAS.md']
 rd.render(docs)
 shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'index.html');shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'report.html')
 reg=json.loads((AT/'FORMULA_REGISTRY.json').read_text())
 entries=''.join('<a data-text="'+html.escape((z['id']+' '+z['title']).lower(),quote=True)+'" href="'+z['chapter'].replace('.md','.html')+'#'+z['anchor']+'"><b>'+z['id']+'</b> '+html.escape(z['title'])+'</a>' for z in reg['formulas'])
 body='<section class="panel"><p>当前 <b>'+str(reg['count'])+'</b> 条可追溯推导。新增R11严格支撑完整导数、总场固定注入与全形状条件耗散界。</p><p><a href="11_ACTIVE_SUPPORT_TOTAL_FIELD_AND_GLOBAL_BOUNDS.html">R11新增12条</a> · <a href="10_TSOM_ON_BILATERAL_TRANSFER.html">TSOM怎样作用于中间逆</a> · <a href="FORMULA_ATLAS.html">连续阅读全文</a> · <a href="FORMULA_REGISTRY.json">JSON索引</a></p><input id="formula-filter" type="search" placeholder="搜索公式编号/支撑/TSOM/概率/Green" aria-label="检索公式"><div class="formula-list">'+entries+'</div></section><script>document.getElementById("formula-filter").addEventListener("input",e=>{const q=e.target.value.toLowerCase();document.querySelectorAll(".formula-list a").forEach(a=>a.hidden=!a.dataset.text.includes(q));});</script>'
 (AT/'index.html').write_text(rd.shell('Maxwell成像：持续维护的公式集',body,AT/'index.html'))
 print(json.dumps(dict(documents=len(docs),formulas=reg['count'],new_physics_RHS=0)))
if __name__=='__main__':main()
