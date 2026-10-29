"""Documentation-only rendering. No Maxwell, model, GT or tests."""
from pathlib import Path
import importlib.util,html,json,re,shutil,datetime,os
P=Path(__file__).resolve().parents[1];R=P.parents[1];AT=R/'research/maxwell_formula_atlas'
spec=importlib.util.spec_from_file_location('oldreader',R/'research/research_z_alg_asymmetric_validation/build_site.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
HEADS={}
manifest=P/'SOURCE_MANIFEST.json'
if manifest.exists():
 for a in json.loads(manifest.read_text()).get('files',[]):HEADS[a.get('repo')]=a.get('git_head','main')
MATH=r'''<script>window.MathJax={tex:{inlineMath:[['$','$'],['\\(','\\)']],displayMath:[['$$','$$'],['\\[','\\]']],processEscapes:true},svg:{fontCache:'global'}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>'''
class Reader(old.Builder):
 def __init__(self,root):
  super().__init__();self.root=root
 def slug(self,t):
  m=re.search(r'F-(?:BASE|OP|OPT|ERR|MAT|PROB|VOL|POST)-\d+',t)
  return m[0].lower() if m else old.Builder.slug(t)
 def url(self,u,source,page,image=False):
  u=u.strip('<>')
  if u.startswith(('https://','http://','mailto:','#')):return u
  if u.startswith('/Volumes/'):
   m=re.search(r'/B1/([^/]+)/(.*)',u)
   if m:
    repo,tail=m.groups();line=re.search(r':(\d+)$',tail);tail=re.sub(r':\d+$','',tail)
    if repo=='b1-p3l-r9-stabilized-subspace':rev='main'
    else:rev=HEADS.get(repo,'main')
    return 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/'+rev+'/'+tail+('#L'+line[1] if line else '')
   return '#source-formula-index'
  target=(source.parent/u.split('#')[0]).resolve();fragment='#'+u.split('#',1)[1] if '#' in u else ''
  if source.parent==AT and not target.exists():
   original=(P/'formula_collection'/source.name).parent/u.split('#')[0]
   if original.exists():target=original.resolve()
  if target in self.pages:target=self.pages[target]
  elif target.suffix=='.md' and target.with_suffix('.html').exists():target=target.with_suffix('.html')
  return os.path.relpath(target,page.parent)+fragment
 def shell(self,title,body,page,subtitle='推导 → 固定与变化算子 → 计算决策 → 实际证据'):
  links=[(AT/'index.html','公式总集'),(P/'index.html','R9研究'),(P/'RESEARCH_REPORT_ZH.md','原始报告'),(AT/'FORMULA_ATLAS.md','完整公式Markdown')]
  nav='<nav>'+''.join('<a href="'+os.path.relpath(f,page.parent)+'">'+t+'</a>' for f,t in links)+'</nav>'
  return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>'+html.escape(title)+'</title><style>'+old.STYLE+''' .math{overflow-x:auto;overflow-wrap:normal} mjx-container{max-width:100%;overflow-x:auto;overflow-y:hidden} #formula-filter{width:100%;padding:12px;border:1px solid #9eb8c7;border-radius:6px;font-size:17px}.formula-list{columns:2}.formula-list a{display:block;break-inside:avoid;padding:5px}.formula-list a[hidden]{display:none}@media(max-width:800px){.formula-list{columns:1}} </style>'''+MATH+'</head><body><header>'+nav+'<h1>'+html.escape(title)+'</h1><p>'+html.escape(subtitle)+'</p></header><main>'+body+'</main><footer style="padding:30px;max-width:1100px;margin:auto">公式编号保持稳定；证明、数值证据和待验证机制分开。近似候选权重不自动成为已认证posterior。每条数据以其来源和范围为准。</footer></body></html>'
 def render(self,docs):
  self.sources={p.resolve():p.read_text() for p in docs};self.pages={p.resolve():p.with_suffix('.html') for p in docs}
  for source,text in self.sources.items():
   page=self.pages[source];body,toc=self.markdown(text,source,page)
   toc='<details><summary>章节目录</summary><div class="toc">'+''.join('<a href="#'+a+'">'+html.escape(l)+'</a>' for _,l,a in toc)+'</div></details>'
   title=self.title(source);page.write_text(self.shell(title,toc+'<article>'+body+'</article>',page))

if __name__=='__main__':
 AT.mkdir(parents=True,exist_ok=True)
 chapters=sorted((P/'formula_collection').glob('[0-9][0-9]_*.md'))
 for f in chapters:shutil.copy2(f,AT/f.name)
 for f in (P/'formula_collection').glob('*.csv'):shutil.copy2(f,AT/f.name)
 intro='''# Maxwell 成像公式总集\n\n这是持续维护的推导主干：固定传播算子、完整非线性材料、投影式SOM/TSOM、合法候选体积，以及按需增强的概率任务积分。公式可被解析证明时不再安排机器精度小实验。先明确哪些量固定、哪些量变化，再决定需要哪一项论文级实验。\n\n当前版本覆盖实际已读、已使用的P1—R9主链和直接相关的F/F3、A17-X、Z-ALG公式；没有声称穷尽未提供的全部历史推导。每章保留源文件、成立条件与资格。完整posterior认证不作为材料可能性体积与风险成像的前置门。\n\n'''
 combined=intro+'\n\n---\n\n'.join(f.read_text() for f in chapters)
 (AT/'FORMULA_ATLAS.md').write_text(combined)
 ids=[]
 for f in chapters:
  for m in re.finditer(r'^#{2,4}\s+(F-[A-Z]+-\d+)\s*([^\n]*)',f.read_text(),re.M):ids.append(dict(id=m[1],title=m[2].lstrip('·｜ |'),chapter=f.name,anchor=m[1].lower(),source='R9 formula-first derivation snapshot; per-item sources inside chapter'))
 (AT/'FORMULA_REGISTRY.json').write_text(json.dumps(dict(version='2026-10-30',count=len(ids),formulas=ids),ensure_ascii=False,indent=2)+'\n')
 if not (AT/'CHANGELOG.md').exists():(AT/'CHANGELOG.md').write_text('# 修订记录\n\n2026-10-30：建立从物理、投影、材料几何、概率任务到候选体积的第一版推导主干。保留稳定编号，区分精确公式、条件结论、数值证据和未运行算法。后续更新追加日期与受影响编号，不覆盖历史实验。\n')
 rr=Reader(AT);rr.render([AT/f.name for f in chapters]+[AT/'FORMULA_ATLAS.md'])
 entries=''.join('<a data-text="'+html.escape((x['id']+' '+x['title']).lower(),quote=True)+'" href="'+x['chapter'].replace('.md','.html')+'#'+x['anchor']+'"><b>'+x['id']+'</b> '+html.escape(x['title'])+'</a>' for x in ids)
 body='<section class="panel"><p>当前收录 <b>'+str(len(ids))+'</b> 个可追溯公式条目。先看推导与计算含义，再决定实验。主目标是可计算的材料可能性体积；概率校准是独立增强目标。</p><p><a href="FORMULA_ATLAS.html">连续阅读全文</a> · <a href="FORMULA_ATLAS.md">原始Markdown</a> · <a href="FORMULA_REGISTRY.json">编号索引JSON</a> · <a href="SOURCE_FORMULA_INDEX.csv">源码追溯CSV</a></p><input id="formula-filter" type="search" placeholder="搜索：Born / 投影 / 体积 / prior / ROI / 误差 / 公式编号" aria-label="搜索公式"><div class="formula-list">'+entries+'</div></section><script>document.getElementById("formula-filter").addEventListener("input",e=>{const q=e.target.value.toLowerCase();document.querySelectorAll(".formula-list a").forEach(a=>a.hidden=!a.dataset.text.includes(q));});</script>'
 (AT/'index.html').write_text(rr.shell('Maxwell 成像：持续维护的公式集',body,AT/'index.html'))
 docs=[p for p in P.glob('*.md') if p.name!='START_HERE.md']+[p for p in (P/'GPT_PRO_HANDOFF').glob('*.md')]
 rd=Reader(P);rd.render(docs)
 if (P/'RESEARCH_REPORT_ZH.html').exists():
  shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'index.html');shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'report.html')
 print(json.dumps(dict(formulas=len(ids),chapters=len(chapters),report_docs=len(docs))))
