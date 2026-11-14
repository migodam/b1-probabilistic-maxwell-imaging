"""Build static supplementary pages from existing report, summary rows and figures.
No optimization, scoring, scientific selection, deployment, or invented downloads.
"""
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit
from html.parser import HTMLParser
import argparse, ast, csv, html, json, os, re, time
from datetime import datetime, timezone

LABELS={'cie':'CIE','dual_cie':'Dual-CIE','g19':'G19','dr':'DR','a_plus_mm_r13':'A+-MM/R13','a_plus_arc_r7':'A+-Arc/R7-readout'}
SOURCES={}
def esc(x):return html.escape(str(x),quote=True)
def relative(p,base):return quote(os.path.relpath(p,base).replace(os.sep,'/'),safe='/_.-')
def readcsv(p):
    with p.open(newline='') as f:rows=list(csv.DictReader(f))
    SOURCES[str(p.relative_to(ROOT))]={'rows':len(rows)}
    return rows

def unique(rows,**query):
    found=[r for r in rows if all(r.get(k)==v for k,v in query.items())]
    if len(found)!=1:raise ValueError(f'Expected one exact row: {query}; found {len(found)}')
    return found[0]

def inline(text,base):
    pattern=r'(`[^`]+`|\*\*[^*]+\*\*|\[[^\]]+\]\([^)]+\)|\\\(.+?\\\))'
    out=[]
    for token in re.split(pattern,text):
        if token.startswith('`') and token.endswith('`'):out.append('<code>'+esc(token[1:-1])+'</code>')
        elif token.startswith('**') and token.endswith('**'):out.append('<strong>'+inline(token[2:-2],base)+'</strong>')
        elif token.startswith('\\(') and token.endswith('\\)'):out.append('<code class="math">'+esc(token[2:-2])+'</code>')
        elif re.fullmatch(r'\[[^\]]+\]\([^)]+\)',token):
            m=re.fullmatch(r'\[([^\]]+)\]\(([^)]+)\)',token);label,url=m.groups();u=urlsplit(url)
            if u.scheme in ('https','http') or (not u.scheme and (base/unquote(u.path)).exists()):out.append('<a href="'+esc(url)+'">'+inline(label,base)+'</a>')
            else:out.append(inline(label,base)+' <span class="note">（附件尚未生成）</span>')
        else:out.append(esc(token))
    return ''.join(out)

def markdown(text,base):
    lines=text.splitlines();out=[];i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line.startswith('```'):
            code=[];i+=1
            while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
            out.append('<pre><code>'+esc('\n'.join(code))+'</code></pre>');i+=1;continue
        m=re.match(r'^(#{1,6})\s+(.+)$',line)
        if m:out.append(f'<h{len(m[1])}>'+inline(m[2],base)+f'</h{len(m[1])}>');i+=1;continue
        if line.startswith('|') and i+1<len(lines) and re.match(r'^\s*\|[\s:|\-]+\|?\s*$',lines[i+1]):
            headers=[x.strip() for x in line.strip('|').split('|')];data=[];i+=2
            while i<len(lines) and lines[i].strip().startswith('|'):
                data.append([x.strip() for x in lines[i].strip().strip('|').split('|')]);i+=1
            out.append('<div class="comparison"><table><thead><tr>'+''.join('<th>'+inline(x,base)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+inline(x,base)+'</td>' for x in row)+'</tr>' for row in data)+'</tbody></table></div>');continue
        if re.match(r'^(?:[-*]|\d+\.)\s+',line):
            tag='ol' if re.match(r'^\d+\.',line) else 'ul';items=[]
            while i<len(lines) and re.match(r'^(?:[-*]|\d+\.)\s+',lines[i].strip()):
                items.append(re.sub(r'^(?:[-*]|\d+\.)\s+','',lines[i].strip()));i+=1
            out.append('<'+tag+'>'+''.join('<li>'+inline(x,base)+'</li>' for x in items)+'</'+tag+'>');continue
        para=[line];i+=1
        while i<len(lines) and lines[i].strip() and not re.match(r'^(?:#|\||[-*]\s|\d+\.\s|```)',lines[i].strip()):para.append(lines[i].strip());i+=1
        out.append('<p>'+inline(' '.join(para),base)+'</p>')
    return '\n'.join(out)

def page(title,body,base):
    nav='<nav><a href="'+relative(ROOT/'index.html',base)+'">补充首页</a><a href="'+relative(ROOT/'report.html',base)+'">完整中文报告</a><a href="'+relative(ROOT.parent/'index.html',base)+'">原R13入口</a></nav>'
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+esc(title)+'</title><style>'+CSS+'</style></head><body><main><header><h1>'+esc(title)+'</h1>'+nav+'</header>'+body+'<footer>R13补充 · 既有公开48对象的后续配对研究，非新的盲确认。正例、负例与共同失败完整保留。</footer></main></body></html>\n'

def table(headers,rows):return '<div class="comparison"><table><thead><tr>'+''.join('<th scope="col">'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def fmt(row,key):return f'{float(row[key]):.5f}'
def figure(p,base):
    url=relative(p,base);pdf=p.with_suffix('.pdf');label=p.stem
    return '<figure id="'+esc(label)+'"><a href="'+url+'"><img loading="lazy" decoding="async" src="'+url+'" alt="'+esc(label)+'"></a><figcaption><strong>'+esc(label)+'</strong> · <a href="'+url+'">原尺寸PNG</a>'+(' · <a href="'+relative(pdf,base)+'">PDF</a>' if pdf.is_file() else '')+'<span class="path">'+esc(p.relative_to(ROOT))+'</span></figcaption></figure>'

def gallery(name,title,files,note,grouped=False):
    target=ROOT/'galleries'/f'{name}.html';target.parent.mkdir(exist_ok=True);base=target.parent
    body='<p>'+esc(note)+'</p><p class="note">'+str(len(files))+'张实际图；保留原文件名与图中完整对象ID。点击图片查看原尺寸，手机可横屏阅读。</p>'
    if grouped:
        groups={}
        for p in files:groups.setdefault(p.stem.split('__')[0],[]).append(p)
        for key,items in groups.items():body+='<details><summary>'+esc(key)+' · '+str(len(items))+'张诊断图</summary><div class="gallery">'+''.join(figure(p,base) for p in items)+'</div></details>'
    else:
        body+='<nav class="pages">'+''.join('<a href="#'+esc(p.stem)+'">'+str(i+1)+'</a>' for i,p in enumerate(files))+'</nav><div class="gallery">'+''.join(figure(p,base) for p in files)+'</div>'
    target.write_text(page(title,body,base));GENERATED.append(target)
    return {'name':name,'title':title,'page':target,'figures':[str(p.relative_to(ROOT)) for p in files]}

def download_html():
    p=ROOT/'DOWNLOADS.json'
    if not p.exists():return '<p class="empty">本补充的实际下载清单尚未生成；发布附件链接将在清单落地后显示。</p>',[]
    data=json.loads(p.read_text());found=[]
    def walk(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if isinstance(v,str) and urlsplit(v).scheme in ('https','http') and ('url' in k.lower() or k.lower() in ('href','download')):found.append((str(x.get('name',x.get('filename',x.get('original_name',k)))),v))
                elif isinstance(v,(dict,list)):walk(v)
        elif isinstance(x,list):
            for v in x:walk(v)
    walk(data);found=list(dict.fromkeys(found));links='<div class="doclist"><a href="DOWNLOADS.json">实际DOWNLOADS.json</a>'+''.join('<a href="'+esc(url)+'">'+esc(label)+'</a>' for label,url in found)+'</div>'
    return links,found

class Links(HTMLParser):
    def __init__(self):super().__init__();self.urls=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ('href','src') and v:self.urls.append(v)

def main():
    global ROOT,CSS,GENERATED
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);args=ap.parse_args();ROOT=args.root.resolve();GENERATED=[];start=time.perf_counter();cpu=time.process_time()
    tree=ast.parse((ROOT.parent/'src/build_site.py').read_text());CSS=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='CSS' for t in n.targets))
    CSS+=' .gallery{display:grid;grid-template-columns:1fr;gap:24px}.gallery figure{padding:14px}.pages a{min-width:34px;padding:4px 9px;border:1px solid #dce3eb;background:white;border-radius:6px}th{white-space:nowrap}td{font-variant-numeric:tabular-nums}.report{max-width:1100px;margin:auto}.report p{max-width:none}.report li{margin:8px 0}pre{white-space:pre-wrap;overflow:auto;background:#fff;padding:16px}code{font-size:.9em;background:#e9eff6;padding:1px 4px;border-radius:3px}.math{white-space:normal}.banner{border-left:5px solid #a65c29;background:#fff4e9;padding:14px 20px}.cards h3{margin-top:0}'
    report=(ROOT/'RESEARCH_REPORT_ZH.md').read_text();report_target=ROOT/'report.html';report_target.write_text(page('R13补充｜完整中文研究报告','<article class="report">'+markdown(report,ROOT)+'</article>',ROOT));GENERATED.append(report_target)
    summary_path=ROOT/'analysis/stack_comparison/METHOD_VIEW_SUMMARY.csv';common_path=ROOT/'analysis/stack_comparison/COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY.csv';rows=readcsv(summary_path);common=readcsv(common_path)
    # CPU display strings copied from the parent-authored report, not re-averaged.
    main_block=report.split('|方法|平均NRMSE|中位NRMSE|平均IoU|球内/48|实际均值CPU秒|')[1].split('\n\n')[0];cpu_strings={}
    for line in main_block.splitlines():
        c=[x.strip().replace('**','') for x in line.strip().strip('|').split('|')]
        if len(c)==6 and c[0] in LABELS.values():cpu_strings[c[0]]=c[5]
    primary=[];audit_rows=[]
    for solver,label in LABELS.items():
        r=unique(rows,cohort='uniform48',solver=solver,initialization='uniform',view='primary');f=unique(common,cohort='uniform48',solver=solver,initialization='uniform',view='common_actual_CPU_six',budget_group='uniform_six')
        assert r['available_cases']=='48' and r['NRMSE_n']=='48'
        primary.append([label,fmt(r,'NRMSE_mean'),fmt(r,'NRMSE_median'),fmt(r,'IoU_mean'),r['ball_inside']+'/48',cpu_strings[label],fmt(f,'NRMSE_mean')]);audit_rows.append({'table':'uniform48','query':{k:r[k] for k in ['cohort','solver','initialization','view']},'exact_primary_row':r,'common_CPU_row':f,'CPU_display_from_parent_report':cpu_strings[label]})
    nn=[]
    for solver,label in LABELS.items():
        values=[];selected=[]
        for initial in ['uniform','NN']:
            r=unique(rows,cohort='factorial16',solver=solver,initialization=initial,view='primary');assert r['available_cases']=='16' and r['NRMSE_n']=='16';values.append(fmt(r,'NRMSE_mean'));selected.append(r)
        for initial in ['uniform','NN']:
            r=unique(common,cohort='factorial16',solver=solver,initialization=initial,view='common_actual_CPU_twelve',budget_group='factorial_twelve');values.append(fmt(r,'NRMSE_mean'));selected.append(r)
        nn.append([label]+values);audit_rows.append({'table':'NN16','solver':solver,'exact_selected_rows':selected})
    specs=[('five_absolute','原五方法 · 全48绝对色标','comparison48','absolute_scale_page*.png','12页全部48对象；原五方法与A+-MM/R13，绝对色标按原图保留。'),('five_row_diagnostic','原五方法 · 全48逐行诊断色标','comparison48','diagnostic_row_scale_page*.png','12页全部48对象；逐行色标仅作结构诊断，不能替代绝对幅度比较。'),('five_curves','原五方法 · 误差与已付费CPU曲线','comparison48','*.png','原comparison48汇总图；使用各图原有预算与指标。'),('arc_mm48','A+-Arc/R7-readout 与 A+-MM/R13 · 全48','stack_comparison','arc_vs_MM_all48_page*.png','8页全部48对象；两版本名称分列，不能混称同一实现。'),('nn16','六solver × 两种初始化 · 全16','stack_comparison','NN_interfaces_all16_page*.png','8页全部16对象；均匀与冻结NN初值，保留正例、退化和共同失败。'),('nn_curves','NN16 · 配对与已付费CPU曲线','stack_comparison','NN_paired16*.png','固定16配对曲线；同12条流水线预算与各six预算分别阅读。'),('dynamics','失败动力学 · 汇总','dynamics_summary','*.png','只展示已保存诊断；原Phi变化与材料误差分开，不据转动大小新增预测器。'),('native','Native状态 · 完整对象ID分组','native_states','*.png','记录副本不替换未插桩主结果；每个对象展开全部实际诊断图。')]
    galleries=[]
    for name,title,campaign,pattern,note in specs:
        files=sorted((ROOT/'analysis'/campaign/'FIGURES').glob(pattern))
        if name=='five_curves':files=[p for p in files if not p.name.startswith(('absolute_scale_page','diagnostic_row_scale_page'))]
        assert all(not any(x in str(p.relative_to(ROOT)).lower() for x in ['private_operational','fixture','fake']) for p in files)
        galleries.append(gallery(name,title,files,note,grouped=name=='native'))
    required={'five_absolute':12,'five_row_diagnostic':12,'arc_mm48':8,'nn16':8,'nn_curves':2,'dynamics':2,'native':48}
    for g in galleries:
        if g['name'] in required:assert len(g['figures'])==required[g['name']],(g['name'],len(g['figures']))
    lead=report.split('## 最终判断',1)[1].split('## 1.',1)[0]
    body='<p class="banner"><strong>既有48对象上的后续配对、版本和接口研究；不是新的盲确认。</strong> 二维scalar-TM、同原观测/prior/box与冻结选根，重复运行不增加独立材料数。</p>'+markdown(lead,ROOT)
    body+='<h2 id="main48">同48 · 六方法primary</h2><p>主规则：预算内、合法、已完成原Maxwell审计的最小原Phi=D+R。球内数是噪声资格，不是材料成功数。数值取原CSV精确primary行，展示保留五位小数；CPU两位小数沿用正文。</p>'+table(['方法','平均NRMSE','中位NRMSE','平均IoU','球内','实际均值CPU秒','同实际CPU平均NRMSE'],primary)
    body+='<p class="note">同实际CPU为每对象六方法的共同已付费前缀，不是新测部署wall。Arc与MM终点接近；DR/G19互补正例、A+退化和共同失败均见完整图集。</p><h2 id="nn16">同16 · 冻结NN初始化接口</h2>'+table(['solver','均匀primary','NN primary','共同12流水线CPU：均匀','共同12流水线CPU：NN'],nn)+'<p class="note">表内数值均为平均NRMSE。右两列使用每对象全部12流水线的共同预算；不跨独立six组相减。冻结权重复用，无新训练；训练/推理准备费用另列。</p>'
    body+='<h2 id="galleries">完整图集</h2><div class="cards">'+''.join('<article class="card"><h3><a href="'+relative(g['page'],ROOT)+'">'+esc(g['title'])+'</a></h3><p>'+str(len(g['figures']))+'张原图 · PNG / 实际存在的PDF</p></article>' for g in galleries)+'</div>'
    docs=[('report.html','完整中文报告'),('RESEARCH_REPORT_ZH.md','报告Markdown'),('THEORY2_MATERIALS.md','Theory2素材与边界'),('TRAJECTORY_ANALYSIS.md','完整轨迹报告'),('HISTORICAL_COMPARABILITY.md','历史可比性'),('COST_LEDGER.md','费用总账'),('CHAT_HANDOFF.md','研究交接'),('REPRODUCE.md','复现'),('FAILURES.md','失败与范围'),('ANALYSIS_CONTRACT.md','分析合同'),('STACK_COMPARISON_CONTRACT.md','版本/NN合同'),('analysis/visual_review/CASE_OBSERVATIONS.md','逐例结构观察'),('analysis/stack_audit/GLOBAL_DUAL16_REFERENCE.md','global-Dual16控制'),('analysis/cpu_frequency/CAUSE_REPORT.md','CPU原因审计'),('analysis/cpu_frequency/EVENT_LIMIT_REPORT.md','限频日志审计')]
    body+='<h2 id="documents">报告、理论与来源</h2><div class="doclist">'+''.join('<a href="'+relative(ROOT/path,ROOT)+'">'+esc(label)+'</a>' for path,label in docs if (ROOT/path).is_file())+'</div><p><a href="analysis/stack_comparison/METHOD_VIEW_SUMMARY.csv">精确summary CSV</a> · <a href="analysis/stack_comparison/COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY.csv">共同CPU summary CSV</a> · <a href="SITE_TABLE_SOURCE_AUDIT.json">网站表格来源审计</a></p>'
    downloads,download_links=download_html();body+='<h2 id="downloads">实际下载</h2>'+downloads
    target=ROOT/'index.html';target.write_text(page('R13补充｜A+全48、NN接口与失败动力学',body,ROOT));GENERATED.append(target)
    audit={'scope':'Presentation only; exact summary-row selection, no new scoring or scientific inference.','table_sources':[str(summary_path.relative_to(ROOT)),str(common_path.relative_to(ROOT)),'RESEARCH_REPORT_ZH.md'],'rows':audit_rows,'gallery_inventory':[{k:v for k,v in g.items() if k!='page'} for g in galleries],'exclusions':['private_operational','fixtures','fake'],'Markdown_renderer':'stdlib bounded renderer; report headings/tables/lists/links/emphasis/code rendered; inline TeX preserved as readable formula source, no external JS required.'}
    (ROOT/'SITE_TABLE_SOURCE_AUDIT.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
    broken=[];private=[]
    for p in GENERATED:
        parser=Links();parser.feed(p.read_text())
        for url in parser.urls:
            u=urlsplit(url)
            if not u.scheme and u.path and not (p.parent/unquote(u.path)).exists():broken.append({'page':str(p.relative_to(ROOT)),'url':url})
            if any(x in url.lower() for x in ['private_operational','fixture','fake']):private.append({'page':str(p.relative_to(ROOT)),'url':url})
    assert not broken,broken;assert not private,private
    receipt={'utc':datetime.now(timezone.utc).isoformat(),'status':'COMPLETE_LOCAL_BUILD','generated_pages':[str(p.relative_to(ROOT)) for p in GENERATED],'gallery_figures':sum(len(g['figures']) for g in galleries),'galleries':{g['name']:len(g['figures']) for g in galleries},'primary48_rows':len(primary),'NN16_rows':len(nn),'broken_local_links':broken,'private_or_fixture_links':private,'download_manifest_exists':(ROOT/'DOWNLOADS.json').exists(),'download_URLs_from_actual_manifest':len(download_links),'source_inventory':SOURCES,'controller_CPU_seconds':time.process_time()-cpu,'wall_seconds':time.perf_counter()-start,'new_physics':0,'new_scoring':0,'published':False,'visual_QA':'No browser screenshot reviewed; local link, figure counts and exact-row cardinalities verified.'}
    (ROOT/'SITE_BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps(receipt,ensure_ascii=False))
if __name__=='__main__':main()
