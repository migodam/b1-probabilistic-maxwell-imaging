"""R6静态科研阅读页生成器；root在最终收尾时执行，本文件不运行任何实验。

Required: RESEARCH_REPORT_ZH.md, RESULTS_SUMMARY.json, GO_NO_GO.md,
PRIOR_EVIDENCE_AND_FAILURES.md, and every referenced/required figure.
Uses stdlib only. No network, solver, training, model load, raw-array read, rglob,
or D-truth access. Writes only managed HTML and SITE_BUILD_RECEIPT.json.

Optional RESULTS_SUMMARY.json presentation metadata (all values author-supplied):
  "site": {
    "title": "...", "subtitle": "...", "intro_markdown": "...",
    "highlights": [{"label":"...", "value":"...", "note":"...",
                    "identity":"material_risk|lowdim_local_probability|physical_accelerator"}],
    "figures": [{"path":"relative.png", "title":"...", "caption":"...",
                 "identity":"...", "sources":["relative.csv"]}],
    "required_figures": ["relative/D_final.png"],
    "appendices": [{"path":"relative.md", "title":"..."}],
    "downloads": [{"path":"relative.zip", "label":"..."}]
  }
If site.figures is omitted, explicit image references in the main Markdown are
used. Numeric results and GO statuses are never inferred from raw data. Missing
figures/required inputs stop before any HTML output. --check-inputs is read-only.
Existing unmarked HTML is never overwritten; original MD/JSON/images are never
modified. Nested appendix pages preserve relative source/image/download links.
"""
from __future__ import annotations
import argparse,datetime,html,json,math,os,posixpath,re
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit,unquote,quote

P=Path(__file__).resolve().parent
MARKER='<!-- R6_STATIC_SITE_GENERATED: build_site.py -->'
REQUIRED=('RESEARCH_REPORT_ZH.md','RESULTS_SUMMARY.json','GO_NO_GO.md','PRIOR_EVIDENCE_AND_FAILURES.md')
IDENTITIES={
 'material_risk':('材料风险动作','材料图或功能值按独立parent的任务损失评价；不自动称为posterior mean或概率认证。'),
 'lowdim_local_probability':('低维局部概率','指定2D/3D条件子模型的积分与latent功能；Q1/Q2数值收敛不等于Full66或全局尾部证书。'),
 'physical_accelerator':('物理计算加速','比较真实锚点路径与完整成本，并保持相同估计器的风险合同；模块微秒不替代系统wall。')}
FIXED_APPENDICES=(
 'PRIOR_EVIDENCE_AND_FAILURES.md','SHARED_STATISTICAL_CONTRACT.md','THEORY_AND_COUNTEREXAMPLES.md',
 'RISK_COST_PARETO.md','RESOURCE_SUMMARY.md','RESOURCE_AND_FAILURE_LEDGER.md','RESOURCE_ACCOUNTING.md',
 'FAILURE_ANALYSIS.md','EVIDENCE_PACKAGE.md','REPRODUCE.md',
 'track_a/analysis/REPORT_ZH.md','track_a/analysis/POSTHOC_LATENT_Q_REPORT_ZH.md','track_a/COPULA_COUNTEREXAMPLE_REPORT_ZH.md',
 'track_b/REPORT_ZH.md','track_b/b6_moments/REPORT_ZH.md','track_c/kkt_closure/REPORT_ZH.md',
 'audit/r6_evidence/CONFORMAL_THEOREM_REVIEW.md')
STYLE='''
:root{--ink:#173346;--muted:#536a78;--blue:#075e7c;--line:#d8e3e9;--bg:#f2f6f9;--paper:#fff;--warn:#a77b25}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--bg);color:var(--ink);font:17px/1.9 -apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}a{color:var(--blue);text-underline-offset:3px}a:hover{color:#003e58}.skip{position:absolute;left:-9999px}.skip:focus{left:12px;top:10px;background:white;padding:8px;z-index:10}
header{background:#102f43;color:#f6fbff;padding:35px max(20px,calc((100vw - 1140px)/2))}header h1{font-size:clamp(27px,4vw,42px);line-height:1.4;margin:18px 0}header p{color:#cfdee8;max-width:1000px}nav{display:flex;gap:10px;flex-wrap:wrap}header nav a{color:white;border:1px solid #668393;border-radius:6px;padding:3px 11px;text-decoration:none}main{max-width:1200px;margin:25px auto;padding:0 24px 55px}.panel,article,figure{background:var(--paper);border:1px solid var(--line);border-radius:11px;padding:clamp(18px,3vw,35px);margin:22px 0;overflow-wrap:anywhere}.grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}.card{background:white;border:1px solid var(--line);border-radius:9px;padding:20px}.card h3{margin:0 0 10px}.card .value{font-size:27px;line-height:1.45;font-weight:700}.muted,figcaption,.source{font-size:14px;color:var(--muted)}.notice{border-left:4px solid var(--warn);background:#fff8e9;padding:15px 20px}.tag{font-size:13px;background:#e9f1f6;border-radius:5px;padding:3px 7px;display:inline-block}h1,h2,h3,h4{line-height:1.5}h2{font-size:25px;margin-top:34px}h3{font-size:20px}article h2{border-top:1px solid var(--line);padding-top:22px}p{margin:14px 0}li{margin:5px 0}.table-wrap{overflow-x:auto;max-width:100%;margin:18px 0}table{border-collapse:collapse;min-width:100%;font-size:15px}th,td{border:1px solid var(--line);padding:10px 13px;text-align:left;vertical-align:top}th{background:#edf3f7}code{background:#edf2f5;border-radius:3px;padding:2px 4px;font-size:.88em}pre{overflow-x:auto;background:#edf2f5;padding:16px;line-height:1.65;font-size:13px}pre code{padding:0}blockquote{border-left:3px solid #93adbc;margin:18px 0;padding:4px 20px;color:var(--muted)}img{display:block;max-width:100%;height:auto;margin:0 auto}figure{padding:20px}figcaption{margin-top:12px}.toc{columns:2;column-gap:28px}.toc a{display:block;margin:6px 0;break-inside:avoid}.toc .sub{padding-left:14px}.downloads{display:flex;gap:12px;flex-wrap:wrap}.downloads a{display:inline-block}details{border:1px solid var(--line);border-radius:7px;padding:12px 16px;background:white}summary{cursor:pointer}footer{border-top:1px solid var(--line);padding:22px 0;color:var(--muted);font-size:14px}.unavailable{color:var(--muted);text-decoration:underline dotted}.math{white-space:pre-wrap;overflow-wrap:anywhere}a:focus-visible{outline:3px solid #9ecce0;outline-offset:3px}@media(max-width:800px){body{font-size:16px}.grid{grid-template-columns:1fr}.toc{columns:1}main{padding:0 12px 35px}header{padding:25px 18px}th,td{padding:8px 10px}h2{font-size:22px}}@media print{body{background:white;color:black;font-size:11pt}header{background:white;color:black;padding:0}header p{color:black}header nav,.skip{display:none}main{max-width:none;padding:0}.panel,article,figure,.card{border-color:#bbb;break-inside:avoid-page}a{color:black}img{max-height:none}.grid{display:block}}
'''

class Builder:
 def __init__(self):
  self.warnings=[];self.missing_images=[];self.missing_links=[];self.sources={};self.pages={};self.summary={};self.site={};self.generated={}
 def text(self,path):
  path=Path(path)
  if path.stat().st_size>2_000_000:raise ValueError('reading document exceeds2MB; do not point renderer at raw assets: '+str(path))
  return path.read_text(encoding='utf-8')
 def local(self,value,reading=False):
  path=(P/str(value)).resolve()
  if not path.is_relative_to(P.resolve()):raise ValueError('presentation input must remain within R6 package: '+str(value))
  if reading and any(part in ('raw','offline_truth','inputs','contexts','scaffold_v1','scaffold_v3') for part in path.relative_to(P).parts):raise ValueError('raw/truth/source-context input not a reading document: '+str(value))
  return path
 def load(self):
  missing=[name for name in REQUIRED if not (P/name).is_file()]
  if missing:raise FileNotFoundError('final inputs not ready: '+', '.join(missing))
  self.summary=json.loads(self.text(P/'RESULTS_SUMMARY.json'))
  self.site=(self.summary.get('site') or {}) if isinstance(self.summary,dict) else {}
  if not isinstance(self.site,dict):raise ValueError('RESULTS_SUMMARY.site must be an object if supplied')
  docs=[P/'RESEARCH_REPORT_ZH.md',P/'GO_NO_GO.md']
  docs += [P/name for name in FIXED_APPENDICES if (P/name).is_file()]
  # Only direct report filenames; never a recursive directory/raw scan.
  docs += sorted(P.glob('STAGE_*.md'))
  handoff=P/'GPT_PRO_HANDOFF'
  if handoff.is_dir():docs+=sorted(handoff.glob('*.md'))
  for item in self.site.get('appendices',[]):
   name=item['path'] if isinstance(item,dict) else item;path=self.local(name,reading=True)
   if path.suffix.lower()!='.md' or not path.is_file():raise FileNotFoundError('declared appendix missing/notMarkdown: '+str(name))
   docs.append(path)
  for rel in self.site.get('section_sources',{}).values():
   path=self.local(rel,reading=True)
   if path.suffix!='.md' or not path.is_file():raise FileNotFoundError('declared section source missing/notMarkdown: '+str(rel))
   docs.append(path)
  for path in dict.fromkeys(x.resolve() for x in docs):
   self.sources[path]=self.text(path);self.pages[path]=P/'report.html' if path==P/'RESEARCH_REPORT_ZH.md' else P/'read'/path.relative_to(P).with_suffix('.html')
  for value in self.site.get('required_figures',[]):
   path=self.local(value)
   if not path.is_file():self.missing_images.append(str(value))
  self.figures=self.site.get('figures') or self.main_figures()
  if not self.figures:raise ValueError('declare actual main images in RESEARCH_REPORT_ZH.md or RESULTS_SUMMARY.site.figures; no image/number is invented')
  for entry in self.figures:
   if not isinstance(entry,dict) or not entry.get('path'):raise ValueError('site.figures entries need relative path')
   path=self.local(entry['path'])
   if not path.is_file():self.missing_images.append(entry['path'])
   if entry.get('identity') and entry['identity'] not in IDENTITIES:raise ValueError('unknown output identity: '+str(entry['identity']))
  if self.missing_images:raise FileNotFoundError('required/main figure not ready: '+', '.join(sorted(set(self.missing_images))))
 def main_figures(self):
  return [dict(path=url,title=label,caption=label) for label,url in re.findall(r'!\[([^\]]*)\]\(([^\s)]+)\)',re.sub(r'(?ms)^```.*?^```[^\n]*','',self.sources[P/'RESEARCH_REPORT_ZH.md'])) if not re.match(r'^(?:https?:|data:)',url)]
 def title(self,source):
  match=re.search(r'^#\s+(.+)$',self.sources.get(source,''),re.M);return re.sub(r'[`*]','',match[1]).strip() if match else source.stem
 def rel(self,path,page):return quote(os.path.relpath(path,page.parent).replace(os.sep,'/'),safe='/._-')
 def url(self,url,source,page,image=False):
  url=html.unescape(url.strip().strip('<>'));parts=urlsplit(url)
  if parts.scheme or parts.netloc:
   if image:raise ValueError('figures must use relative package paths, not remote/data resources: '+url)
   if parts.scheme not in ('http','https','mailto'):return '#'
   return url
  if not parts.path:return '#'+quote(unquote(parts.fragment),safe='-_') if parts.fragment else '#'
  path=(source.parent/unquote(parts.path)).resolve()
  if image and not path.is_relative_to(P):raise ValueError('image outside package: '+url)
  if image and not path.is_file():self.missing_images.append(str(path.relative_to(P)))
  if not path.exists():self.missing_links.append(dict(source=str(source.relative_to(P)),target=url))
  target=self.pages.get(path,path)
  if path.is_absolute() and not path.is_relative_to(P) and parts.path.startswith('/'):
   self.warnings.append('outside-package absolute source link omitted: '+url);return '#'
  return urlunsplit(('', '', self.rel(target,page),parts.query,quote(unquote(parts.fragment),safe='-_')))
 def inline(self,text,source,page):
  tokens=[]
  def token(value):tokens.append(value);return '\x00'+str(len(tokens)-1)+'\x00'
  text=re.sub(r'`([^`]+)`',lambda m:token('<code>'+html.escape(m[1])+'</code>'),text)
  def linked(m):
   image,label,url=m[1],m[2],m[3];href=self.url(url,source,page,bool(image));label=html.escape(label)
   if image:return token('<img loading="lazy" decoding="async" src="'+html.escape(href,quote=True)+'" alt="'+label+'">')
   return token('<a href="'+html.escape(href,quote=True)+'">'+label+'</a>')
  text=re.sub(r'(!?)\[([^\]]+)\]\((<[^>]+>|[^\s)]+)\)',linked,text);text=html.escape(text)
  text=re.sub(r'\*\*(.+?)\*\*',r'<strong>\1</strong>',text);text=text.replace('&lt;br&gt;','<br>').replace('&lt;br/&gt;','<br>')
  return re.sub(r'\x00(\d+)\x00',lambda m:tokens[int(m[1])],text)
 @staticmethod
 def slug(text):
  text=re.sub(r'[`*]','',text);text=re.sub(r'[^\w\s-]','',text.lower());return re.sub(r'\s+','-',text.strip())
 @staticmethod
 def cells(line):
  result=[];current=[];in_code=False;escaped=False
  for char in line.strip().strip('|'):
   if escaped:current.append(char);escaped=False;continue
   if char=='\\':escaped=True;current.append(char);continue
   if char=='`':in_code=not in_code
   if char=='|' and not in_code:result.append(''.join(current).replace('\\|','|').strip());current=[]
   else:current.append(char)
  result.append(''.join(current).replace('\\|','|').strip());return result
 def markdown(self,text,source,page):
  lines=text.splitlines();out=[];paragraph=[];toc=[];ids={};i=0
  def flush():
   if paragraph:out.append('<p>'+self.inline(' '.join(paragraph),source,page)+'</p>');paragraph.clear()
  while i<len(lines):
   line=lines[i];stripped=line.strip()
   if re.match(r'^\s*```',line):
    flush();language=stripped[3:].strip();code=[];i+=1
    while i<len(lines) and not lines[i].strip().startswith('```'):code.append(lines[i]);i+=1
    out.append('<pre><code data-language="'+html.escape(language,quote=True)+'">'+html.escape('\n'.join(code))+'</code></pre>')
   elif stripped=='$$':
    flush();mathlines=[];i+=1
    while i<len(lines) and lines[i].strip()!='$$':mathlines.append(lines[i]);i+=1
    out.append('<pre class="math" aria-label="数学公式">'+html.escape('\n'.join(mathlines))+'</pre>')
   elif re.match(r'^#{1,6}\s+',line):
    flush();m=re.match(r'^(#{1,6})\s+(.*)',line);level=len(m[1]);label=m[2];base=self.slug(label);n=ids.get(base,0);ids[base]=n+1;anchor=base+('-'+str(n) if n else '');out.append(f'<h{level} id="{html.escape(anchor,quote=True)}">'+self.inline(label,source,page)+f'</h{level}>')
    if 2<=level<=3:toc.append((level,re.sub(r'[`*]','',label),anchor))
   elif stripped.startswith('|') and i+1<len(lines) and re.match(r'^\s*\|?\s*:?-{2,}',lines[i+1]):
    flush();table=[]
    while i<len(lines) and lines[i].strip().startswith('|'):
     cells=self.cells(lines[i])
     if not all(re.match(r'^:?-{2,}:?$',c.replace(' ','')) for c in cells):table.append(cells)
     i+=1
    parts=[]
    for n,row in enumerate(table):tag='th' if n==0 else 'td';parts.append('<tr>'+''.join('<'+tag+'>'+self.inline(c,source,page)+'</'+tag+'>' for c in row)+'</tr>')
    out.append('<div class="table-wrap"><table>'+''.join(parts)+'</table></div>');i-=1
   elif re.match(r'^\s*(?:[-*+] |\d+\. )',line):
    flush();ordered=bool(re.match(r'^\s*\d+\. ',line));tag='ol' if ordered else 'ul';items=[]
    while i<len(lines) and re.match(r'^\s*(?:[-*+] |\d+\. )',lines[i]):items.append(re.sub(r'^\s*(?:[-*+] |\d+\. )','',lines[i]));i+=1
    out.append('<'+tag+'>'+''.join('<li>'+self.inline(v,source,page)+'</li>' for v in items)+'</'+tag+'>');i-=1
   elif stripped.startswith('>'):
    flush();quote_lines=[]
    while i<len(lines) and lines[i].strip().startswith('>'):quote_lines.append(re.sub(r'^\s*>\s?','',lines[i]));i+=1
    out.append('<blockquote>'+self.inline(' '.join(quote_lines),source,page)+'</blockquote>');i-=1
   elif re.fullmatch(r'!\[([^\]]*)\]\(([^\s)]+)\)',stripped):
    flush();m=re.fullmatch(r'!\[([^\]]*)\]\(([^\s)]+)\)',stripped);url=self.url(m[2],source,page,True);out.append('<figure><a href="'+html.escape(url,quote=True)+'"><img loading="lazy" decoding="async" src="'+html.escape(url,quote=True)+'" alt="'+html.escape(m[1],quote=True)+'"></a><figcaption>'+html.escape(m[1])+'</figcaption></figure>')
   elif re.fullmatch(r'[-*_]{3,}',stripped):flush();out.append('<hr>')
   elif not stripped:flush()
   else:paragraph.append(stripped)
   i+=1
  flush();return '\n'.join(out),toc
 def nav(self,page):
  links=[(P/'index.html','阅读入口'),(P/'report.html','完整报告'),(self.pages[P/'GO_NO_GO.md'],'路线判定'),(P/'appendices.html','附录与证据'),(P/'RESEARCH_REPORT_ZH.md','原始Markdown'),(P/'RESULTS_SUMMARY.json','结果JSON')]
  return '<nav aria-label="研究导航">'+''.join('<a href="'+self.rel(t,page)+'">'+label+'</a>' for t,label in links)+'</nav>'
 def shell(self,title,body,page,subtitle=None):
  subtitle=subtitle or '实际数据与范围分别呈现 · 原始失败保留 · 网页生成不增加实验资格'
  return '<!doctype html>\n'+MARKER+'\n<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="R6任务最优Maxwell材料计算：材料风险、低维条件概率、真实物理成本分别呈现"><title>'+html.escape(title)+'</title><style>'+STYLE+'</style></head><body><a class="skip" href="#main">跳到正文</a><header>'+self.nav(page)+'<h1>'+html.escape(title)+'</h1><p>'+html.escape(subtitle)+'</p></header><main id="main">'+body+'<footer>本页只读呈现实际来源；静态网页生成不改变实验、统计、物理或论文资格。JSON/Markdown/图表保留相对下载链接；未记录的资源不能填成零。</footer></main></body></html>'
 def doc_page(self,source,page):
  body,toc=self.markdown(self.sources[source],source,page);toc_html='<div class="toc">'+''.join('<a class="'+('sub' if n==3 else '')+'" href="#'+html.escape(a,quote=True)+'">'+html.escape(label)+'</a>' for n,label,a in toc)+'</div>' if toc else ''
  sourcebar='<p class="source">来源：<a href="'+self.rel(source,page)+'" download>原始Markdown</a>。附录按来源版本呈现，最终状态以完整报告及路线判定为准。</p>'
  return self.shell(self.title(source),sourcebar+toc_html+'<article>'+body+'</article>',page)
 def overview(self):
  md=self.sources[P/'RESEARCH_REPORT_ZH.md'];md=re.sub(r'^#\s+.*\n','',md,count=1);parts=re.split(r'^##\s+',md,maxsplit=1,flags=re.M);lead=parts[0].strip()
  return lead if len(lead)<16000 else '\n\n'.join(lead.split('\n\n')[:5])
 def figure_html(self,entry,page):
  path=self.local(entry['path']);url=self.rel(path,page);title=entry.get('title') or path.stem;caption=entry.get('caption','');tag=IDENTITIES[entry['identity']][0] if entry.get('identity') else '图的适用范围以原文为准';sources=[]
  for item in entry.get('sources',[]):
   name=item['path'] if isinstance(item,dict) else item;target=self.local(name)
   if not target.is_file():raise FileNotFoundError('declared figure source missing: '+str(name))
   label=item.get('label',Path(name).name) if isinstance(item,dict) else Path(name).name;sources.append('<a href="'+self.rel(target,page)+'">'+html.escape(label)+'</a>')
  return '<figure><p class="tag">'+html.escape(tag)+'</p><h3>'+html.escape(title)+'</h3><a href="'+url+'"><img loading="lazy" decoding="async" src="'+url+'" alt="'+html.escape(title,quote=True)+'"></a><figcaption>'+html.escape(caption)+('</figcaption><p class="source">图源：'+' · '.join(sources)+'</p>' if sources else '</figcaption>')+'</figure>'
 def section_source(self,name,defaults,pattern,page):
  override=self.site.get('section_sources',{}).get(name)
  candidates=([override] if override else [])+list(defaults)
  for rel in candidates:
   source=self.local(rel,reading=True)
   if source in self.sources:return self.markdown(self.sources[source],source,page)[0]
  # Extract actual matching main-report sections; never synthesize outcomes.
  text=self.sources[P/'RESEARCH_REPORT_ZH.md'];chunks=re.split(r'(?=^##\s+)',text,flags=re.M);chosen=[v for v in chunks if re.search(pattern,v.split('\n',1)[0])]
  if chosen:return self.markdown('\n\n'.join(chosen),P/'RESEARCH_REPORT_ZH.md',page)[0]
  return '<p>请按<a href="report.html">完整报告</a>与<a href="RESULTS_SUMMARY.json">原始结果JSON</a>阅读此项；此入口不补写缺失结论。</p>'
 def landing(self):
  page=P/'index.html';title=self.site.get('title') or self.title(P/'RESEARCH_REPORT_ZH.md');intro=self.site.get('intro_markdown') or '从复数散射测量恢复材料，需要同时考虑材料图误差、条件概率的适用范围，以及实际求解Maxwell方程的费用。这三个输出各有自己的证据，不能互相替代。'
  intro_html=self.markdown(intro,P/'RESEARCH_REPORT_ZH.md',page)[0];cards='<div class="grid">'+''.join('<section class="card"><h3>'+label+'</h3><p>'+desc+'</p></section>' for label,desc in IDENTITIES.values())+'</div>'
  highlights=[]
  for item in self.site.get('highlights',[]):
   if not isinstance(item,dict) or 'label' not in item or 'value' not in item:raise ValueError('highlights require author-supplied label/value')
   value=item['value'];value='未记录' if value is None or isinstance(value,float) and not math.isfinite(value) else str(value);ident=IDENTITIES.get(item.get('identity'),('',))[0]
   highlights.append('<div class="card"><span class="tag">'+html.escape(ident)+'</span><div class="value">'+html.escape(value)+'</div><p>'+html.escape(str(item['label']))+'</p><p class="muted">'+html.escape(str(item.get('note','')))+'</p></div>')
  history=self.markdown(self.sources[P/'PRIOR_EVIDENCE_AND_FAILURES.md'],P/'PRIOR_EVIDENCE_AND_FAILURES.md',page)[0];overview=self.markdown(self.overview(),P/'RESEARCH_REPORT_ZH.md',page)[0];decision=self.markdown(self.sources[P/'GO_NO_GO.md'],P/'GO_NO_GO.md',page)[0]
  body='<section class="panel"><h2>先分清物理问题与三类输出</h2>'+intro_html+cards+'</section><section class="panel" id="history"><h2>本轮起点与历史事实</h2>'+history+'</section><section class="panel" id="results"><h2>本轮真实结果与适用范围</h2>'+overview+('<div class="grid">'+''.join(highlights)+'</div>' if highlights else '')+'<p><a href="report.html">阅读全文、所有区间与口径 →</a></p></section><section id="figures"><h2>实际主图</h2>'+''.join(self.figure_html(f,page) for f in self.figures)+'</section><section class="panel" id="negative"><h2>负结果与未运行项目</h2>'+self.section_source('negative',('FAILURE_ANALYSIS.md',),'失败|负结果|未成立|限制|NOT_RUN',page)+'</section><section class="panel" id="cost"><h2>完整成本与资源缺口</h2>'+self.section_source('cost',('RESOURCE_SUMMARY.md','RESOURCE_AND_FAILURE_LEDGER.md','RESOURCE_ACCOUNTING.md'),'成本|费用|资源',page)+'</section><section class="panel" id="judgments"><h2>各路线判定：原文呈现</h2>'+decision+'</section><section class="panel"><h2>GPT Pro与证据入口</h2>'+self.appendix_links(page,True)+'</section>'
  return self.shell(title,body,page,self.site.get('subtitle'))
 def appendix_links(self,page,brief=False):
  items=[]
  for source,target in self.pages.items():
   if source.name=='RESEARCH_REPORT_ZH.md':continue
   if brief and source.parent.name!='GPT_PRO_HANDOFF' and source.name not in ('EVIDENCE_PACKAGE.md','REPRODUCE.md','THEORY_AND_COUNTEREXAMPLES.md','SHARED_STATISTICAL_CONTRACT.md'):continue
   items.append('<a href="'+self.rel(target,page)+'">'+html.escape(self.title(source))+'</a>')
  if brief:items+=['<a href="'+self.rel(P/'appendices.html',page)+'">全部阶段报告、附录与原始Markdown</a>']
  for item in self.site.get('downloads',[]):
   value=item['path'];parts=urlsplit(value)
   if parts.scheme in ('http','https'):url=value
   else:
    target=self.local(value)
    if not target.is_file():raise FileNotFoundError('declared download missing: '+value)
    url=self.rel(target,page)
   items.append('<a href="'+html.escape(url,quote=True)+'" download>'+html.escape(item.get('label',value))+'</a>')
  return '<div class="toc">'+''.join(items)+'</div>'
 def assemble(self):
  self.load()
  for source,page in self.pages.items():self.generated[page]=self.doc_page(source,page)
  self.generated[P/'index.html']=self.landing();page=P/'appendices.html';table='<div class="table-wrap"><table><tr><th>阅读页</th><th>原始Markdown</th></tr>'
  for source,target in self.pages.items():table+='<tr><td><a href="'+self.rel(target,page)+'">'+html.escape(self.title(source))+'</a></td><td><a href="'+self.rel(source,page)+'" download>'+html.escape(str(source.relative_to(P)))+'</a></td></tr>'
  table+='</table></div>';self.generated[page]=self.shell('R6：阶段报告、附录与证据',self.appendix_links(page)+table,page)
  if self.missing_images:raise FileNotFoundError('referenced appendix image missing: '+', '.join(sorted(set(self.missing_images))))
  for output in self.generated:
   if output.exists() and MARKER not in output.read_text(encoding='utf-8')[:256]:raise FileExistsError('refuse to overwrite unmarked/frozen HTML: '+str(output))
  return dict(status='INPUTS_CHECKED',required=list(REQUIRED),reading_documents=len(self.sources),HTML_pages=len(self.generated),main_figures=len(self.figures),missing_local_links=self.missing_links,warnings=self.warnings,raw_or_truth_arrays_read=0,physics_RHS=0,training_jobs=0)
 def commit(self,receipt):
  prior=P/'SITE_BUILD_RECEIPT.json'
  if prior.exists() and json.loads(prior.read_text()).get('producer')!='build_site.py':raise FileExistsError('refuse to overwrite unowned build receipt')
  for output,body in self.generated.items():
   output.parent.mkdir(parents=True,exist_ok=True);temp=output.with_suffix(output.suffix+'.tmp');temp.write_text(body,encoding='utf-8');temp.replace(output)
  receipt.update(status='GENERATED_STATIC_READING_SITE',producer='build_site.py',generated_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),outputs=[str(v.relative_to(P)) for v in self.generated],sources=[str(v.relative_to(P)) for v in self.sources],result_JSON='RESULTS_SUMMARY.json',JS_bytes=0,remote_assets_added=0,source_claims_or_frozen_files_modified=False)
  prior.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check-inputs',action='store_true',help='assemble/check in memory only; no output files');args=parser.parse_args();builder=Builder();receipt=builder.assemble()
 if not args.check_inputs:builder.commit(receipt)
 print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
