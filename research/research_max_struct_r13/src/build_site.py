"""Build local R13 navigation/galleries from existing analysis figures and docs.

No result values/truth, optimization, deployment or invented release links.
Only analysis/<campaign>/**/FIGURES is eligible; delegated fake fixtures excluded.
"""
import argparse,html,json,re,time
from pathlib import Path
from urllib.parse import quote

CSS=''' .comparison{overflow-x:auto;margin:22px 0}table{border-collapse:collapse;width:100%;background:white}th,td{padding:10px 14px;text-align:left;border:1px solid #dce3eb}th{background:#e9eff6}.lead{font-size:20px;line-height:1.7}.findings li{margin:10px 0} *{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:#f3f5f8;color:#17253a;font:16px/1.7 -apple-system,BlinkMacSystemFont,"PingFang SC","Microsoft YaHei",sans-serif}main{max-width:1520px;margin:auto;padding:32px 24px 70px}header{padding:24px 0 18px;border-bottom:1px solid #d9e0e8;margin-bottom:24px}h1{font-size:clamp(27px,4vw,43px);line-height:1.25;margin:0 0 14px}h2{font-size:25px;margin-top:38px}h3{font-size:19px;margin:0 0 8px}p{max-width:1000px}a{color:#165a89;text-decoration:none}a:hover{text-decoration:underline}.note{color:#526479}.badge{display:inline-block;font-size:13px;background:#e9eff6;border-radius:7px;padding:3px 9px;margin:2px 5px 2px 0}.cards,.figures{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}.card,figure{background:white;border:1px solid #dce3eb;border-radius:13px;padding:20px;min-width:0}.card p{margin:8px 0}.card a,.card small{overflow-wrap:anywhere}nav{display:flex;flex-wrap:wrap;gap:12px 22px;padding:12px 0}figure{margin:0;padding:12px}figure img{display:block;width:100%;height:auto;border-radius:6px}figcaption{padding:10px 4px 0;font-size:14px;overflow-wrap:anywhere}.path{display:block;color:#69798a;font-size:12px}.doclist{display:flex;flex-wrap:wrap;gap:10px}.doclist a{background:#fff;border:1px solid #dce3eb;border-radius:8px;padding:8px 12px}details{background:#fff;border:1px solid #dce3eb;border-radius:10px;padding:14px;margin:12px 0}summary{cursor:pointer;font-weight:600}footer{border-top:1px solid #dce3eb;padding-top:20px;margin-top:38px;color:#69798a;font-size:13px}.empty{padding:20px;border-left:4px solid #9eb2c8;background:#fff}code{overflow-wrap:anywhere} @media(max-width:760px){main{padding:20px 12px 45px}.cards,.figures{grid-template-columns:1fr;gap:14px}.card{padding:16px}figure{padding:7px}h2{font-size:22px}nav{gap:8px 15px}.doclist a{font-size:14px}header{padding-top:10px}}'''
GROUPS={
 'development':('已暴露对象：开发与表示对照','历史/公开对象上的开发附件，与新对象确认分开阅读。'),
 'fresh48':('新48对象：独立确认','仅收录实际存在的封存后分析附件；运行配置本身不代表已评分。'),
 'mechanism16':('预先指定16对象：机制与NN初值','固定子集的机制与冻结NN初值附件；不是另加16个独立材料。'),
 'rescue48':('条件救援：48对象策略','策略分母、触发与旁路以实际分析报告为准。'),
 'physics':('物理合同：矢量3D、实测与经典对象','不同物理域的图集分别列出，不合并材料指标。'),
 'other':('其他实际生成附件','按真实目录列出，不从文件名推断算法优劣。')}
DOCS=[('BASELINE_AUDIT.md','基线审查'),('METHODS_AND_ABLATIONS.md','方法与消融'),('BENCHMARKS.md','基准与合同'),('CHAT_HANDOFF.md','研究交接'),('THEORY_STATUS.md','理论状态'),('PRIOR_ART.md','先例与范围'),('PRIOR_ART_AND_NOVELTY.md','先例与新颖性审查'),('RESEARCH_REPORT_ZH.md','中文研究报告'),('REPORT.md','研究报告'),('PAPER_DRAFT.md','论文初稿'),('THEORY_DERIVATIONS.md','数学推导'),('REPRODUCE.md','复现与来源'),('COST_LEDGER.md','费用账目'),('FAILURES.md','失败与范围'),('ARTIFACTS.md','原始附件与下载'),('TASK_RESUME.md','任务连续性'),('delegated/evidence/PRIOR_ART.md','最近先例审查'),('delegated/evidence/SOURCE_REUSE.md','历史证据与复用')]
DATA_DOCS=['REPORT.md','SCORING_RECEIPT.json','METHOD_VIEW_SUMMARY.csv','PAIRED_STATS.csv','CORRECTOR_DECOMPOSITION.csv','COST_LEDGER.csv','COST_TOTALS.json','FAILURES.json','TRUTH_SOURCES.json','QUEUE_ATTEMPT_HISTORY.json']


def esc(value):return html.escape(str(value),quote=True)


def link(path,base):return quote(Path(path).relative_to(base).as_posix(),safe='/_.-')


def relative(path,base):
    import os
    return quote(os.path.relpath(path,base).replace(os.sep,'/'),safe='/_.-')


def group(name):
    lower=name.lower()
    if 'rescue' in lower or 'conditional' in lower:return 'rescue48'
    if 'mechanism16' in lower or 'fresh16' in lower or 'nn16' in lower:return 'mechanism16'
    if 'core48' in lower or 'fresh48' in lower:return 'fresh48'
    if 'physics' in lower or 'fresnel' in lower or 'vector' in lower or 'classic' in lower:return 'physics'
    if 'development' in lower or 'representation' in lower:return 'development'
    return 'other'


def heading(path,label):
    """Read only an existing document's first heading; no result synthesis."""
    with path.open(errors='replace') as stream:
        for line in stream:
            if line.startswith('# '):return label+' · '+line[2:].strip()
    return label


def page(title,body,nav=''):
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+esc(title)+'</title><style>'+CSS+'</style></head><body><main><header><h1>'+esc(title)+'</h1>'+nav+'</header>'+body+'<footer>MAX-STRUCT R13 · 报告、真实图件、原始结果与复现入口。科学适用条件及完整失败见研究报告。</footer></main></body></html>\n'


def figure_html(path,base):
    url=relative(path,base);pdf=path.with_suffix('.pdf');extra=' · <a href="'+relative(pdf,base)+'">PDF</a>' if pdf.is_file() else''
    return '<figure><a href="'+url+'" target="_blank" rel="noopener"><img loading="lazy" decoding="async" src="'+url+'" alt="'+esc(path.stem)+'"></a><figcaption><a href="'+url+'">'+esc(path.stem)+'</a>'+extra+'<span class="path">'+esc(path.parent.name+'/'+path.name)+'</span></figcaption></figure>'


def web_summary(root):
    """Render only a parent-authored optional summary; never infer conclusions."""
    source=root/'WEB_SUMMARY.json'
    if not source.is_file():return None,'',[]
    d=json.loads(source.read_text());body='';used=[]
    if d.get('lead'):body+='<p class="lead">'+esc(d['lead'])+'</p>'
    if d.get('findings'):
        body+='<ul class="findings">'
        for item in d['findings']:
            body+='<li>'+('<strong>'+esc(item['label'])+'</strong> · ' if item.get('label') else '')+esc(item['text'])+'</li>'
        body+='</ul>'
    table=d.get('table')
    if table:
        columns=table['columns'];rows=table['rows']
        if not all(isinstance(x,str) for x in columns) or any(len(row)!=len(columns) or not all(isinstance(x,str) for x in row) for row in rows):raise ValueError('WEB_SUMMARY table requires equally sized string rows/columns')
        body+='<div class="comparison"><table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in columns)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
    figures=d.get('figures',[])
    if len(figures)>4:raise ValueError('WEB_SUMMARY permits at most four curated figures; full galleries remain available')
    if figures:body+='<div class="figures">'
    for item in figures:
        raw=Path(item['path'])
        path=(root/raw).resolve();analysis=(root/'analysis').resolve()
        if raw.is_absolute() or not path.is_relative_to(analysis) or not path.is_file() or path.suffix.lower() not in ['.png','.jpg','.jpeg','.svg','.webp'] or any('fixture' in part.lower() or 'fake' in part.lower() for part in path.relative_to(analysis).parts):raise ValueError('Curated figure must be an existing nonfixture analysis image')
        url=relative(path,root);body+='<figure><a href="'+url+'"><img loading="lazy" decoding="async" src="'+url+'" alt="'+esc(item.get('caption',path.stem))+'"></a><figcaption>'+esc(item.get('caption',path.stem))+'</figcaption></figure>';used.append(relative(path,root))
    if figures:body+='</div>'
    if d.get('limitations'):body+='<h2>结论边界</h2><ul>'+''.join('<li>'+esc(x)+'</li>' for x in d['limitations'])+'</ul>'
    return d.get('title'),body,used


def domain_name(directory):
    name=directory.name.removeprefix('DOMAIN_')
    if any(x in name.lower() for x in ['fresnel','measured']):return '实测域 · '+name
    if any(x in name.lower() for x in ['3d','vector']):return '矢量三维域 · '+name
    if any(x in name.lower() for x in ['classic','austria']):return '经典对象域 · '+name
    return name


def campaign_gallery(root,campaign,figures):
    base=campaign;title=campaign.name+' · 实际图集';docs=[base/n for n in DATA_DOCS if(base/n).is_file() and (base/n).stat().st_size<2_000_000]
    nav='<nav><a href="../../index.html">R13首页</a><a href="../../GALLERY.html">全部图集</a><a href="GALLERY.md">GitHub Markdown图集</a></nav>'
    body='<p class="note">只列出本目录实际存在的图和分析文件。点击图片打开原尺寸；每张图的解释以已保存报告为准。</p><div class="doclist">'+''.join('<a href="'+relative(p,base)+'">'+esc(p.name)+'</a>' for p in docs)+'</div>'
    review=root/'analysis/review_tables/README.md'
    if review.is_file():body+='<p><a href="'+relative(review,base)+'">紧凑全对象复核表</a></p>'
    md=['# '+title,'','[R13首页](../../index.html) · [HTML图集](GALLERY.html)','','仅列出实际生成附件；不以运行状态或图名推断科学结论。','']
    if review.is_file():md+=['[紧凑全对象复核表]('+relative(review,base)+')','']
    if docs:md+=[' · '.join('['+p.name+']('+relative(p,base)+')' for p in docs),'']
    bydomain={}
    for f in figures:
        rel=f.relative_to(base);domain=rel.parts[0] if len(rel.parts)>2 else'campaign'
        bydomain.setdefault(domain,[]).append(f)
    for domain,items in sorted(bydomain.items()):
        label=domain_name(base/domain) if domain!='campaign' else'本组图'
        body+='<h2>'+esc(label)+'</h2><div class="figures">'+''.join(figure_html(f,base) for f in items)+'</div>'
        md+=['## '+label,'']
        for f in items:
            md+=['### '+f.stem,'','!['+f.stem+']('+relative(f,base)+')','']
            if f.with_suffix('.pdf').is_file():md+=['[PDF]('+relative(f.with_suffix('.pdf'),base)+')','']
    if not figures:body+='<p class="empty">当前目录没有实际图文件；以上仅为已存在的分析文档。</p>'
    (base/'GALLERY.html').write_text(page(title,body,nav));(base/'GALLERY.md').write_text('\n'.join(md)+'\n')
    return bydomain


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();root=a.root.resolve();started=time.perf_counter();campaigns=[];summary_title,summary_body,summary_figures=web_summary(root)
    for directory in sorted((root/'analysis').iterdir()) if(root/'analysis').is_dir() else[]:
        if not directory.is_dir() or any(x in directory.name.lower() for x in ['fixture','fake']):continue
        figures=sorted(f for f in directory.rglob('*') if f.is_file() and f.suffix.lower() in ['.png','.jpg','.jpeg','.svg','.webp'] and 'FIGURES' in f.relative_to(directory).parts and not any('fixture' in x.lower() or 'fake' in x.lower() for x in f.relative_to(directory).parts))
        documents=[directory/x for x in DATA_DOCS if(directory/x).is_file() and (directory/x).stat().st_size<2_000_000]
        if not figures and not documents:continue
        domains=campaign_gallery(root,directory,figures);campaigns.append(dict(name=directory.name,path=directory,group=group(directory.name),figures=figures,documents=documents,domains=list(domains)))
    docs=[(root/path,label) for path,label in DOCS if(root/path).is_file()]
    # Link only existing download manifests/docs. Never guess a future Release URL.
    downloads=[]
    for directory in [root,root/'archives']:
        if directory.is_dir():
            downloads.extend(f for f in directory.iterdir() if f.is_file() and ('download' in f.name.lower() or 'release' in f.name.lower() and 'manifest' in f.name.lower()) and f.suffix.lower() in ['.json','.md','.csv'])
    body=summary_body+'<p>报告、数学推导、完整分析图集与复现入口。新对象确认、既有开发、经典对象、实测及三维分别列出，完整保留正面和负面结果。</p>'
    body+='<h2 id="documents">实际文档</h2><div class="doclist">'+''.join('<a href="'+link(path,root)+'">'+esc(label)+'</a>' for path,label in docs)+'</div>'
    body+='<h2 id="protocol">执行合同</h2><div class="cards"><article class="card"><h3>C1 · 开发与表示对照</h3><p>历史四P核队列；mod50投影/材料步审计及阶段内原生审计按各配置记录。已暴露对象与后续独立对象分开。</p></article><article class="card"><h3>C2 · 独立确认</h3><p>八P核、每BLAS一线程。统一绝对进程CPU每2秒桶，在下一次完整材料更新后付费审计；原生初始/阶段/终点审计保留。corrector与direct每次已付物理目标评价记录。</p><p class="note">桶可能被不可中断操作跨过，不补造审计；全部实际费用保留。M4本地控制另列，不与i9计时直接混合。</p></article></div>'
    for key,(title,note) in GROUPS.items():
        members=[c for c in campaigns if c['group']==key]
        if not members:continue
        body+='<h2 id="'+key+'">'+title+'</h2><p class="note">'+note+'</p><div class="cards">'
        for c in members:
            url=link(c['path']/'GALLERY.html',root);body+='<article class="card"><h3>'+esc(c['name'])+'</h3><span class="badge">'+str(len(c['figures']))+'张实际图</span><span class="badge">'+str(len(c['domains']))+'个目录域</span><p><a href="'+url+'">打开分页完整图集</a> · <a href="'+link(c['path']/'GALLERY.md',root)+'">GitHub图集</a></p>'
            if(c['path']/'REPORT.md').is_file():body+='<p><a href="'+link(c['path']/'REPORT.md',root)+'">原分析报告</a></p>'
            body+='<small>'+esc(', '.join(c['domains']))+'</small></article>'
        body+='</div>'
    if downloads:body+='<h2>实际下载清单</h2><div class="doclist">'+''.join('<a href="'+link(path,root)+'">'+esc(path.name)+'</a>' for path in downloads)+'</div>'
    manifest=root/'DOWNLOADS.json'
    if manifest.is_file():
        data=json.loads(manifest.read_text())
        if data.get('published') is True:
            body+='<p>以下附件已公开并读回核对。大文件分片按编号二进制连接；下载和解压建议使用外置盘。初次失败、恢复来源与规范结果分别保留，勿覆盖混用。</p><div class="comparison"><table><thead><tr><th>原始附件</th><th>字节</th><th>下载分片（依次连接）</th></tr></thead><tbody>'
            for asset in data['logical_assets']:
                parts=' · '.join('<a href="'+esc(part['url'])+'">'+esc(part['name'])+'</a>' for part in asset['parts'] if part.get('url_readback_verified') is True)
                body+='<tr><td>'+esc(asset['name'])+'</td><td>'+format(asset['bytes'],',')+'</td><td>'+parts+'</td></tr>'
            body+='</tbody></table></div>'
    nav='<nav><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging">GitHub项目首页</a>'+('<a href="RESEARCH_REPORT_ZH.md">中文研究报告</a>' if (root/'RESEARCH_REPORT_ZH.md').is_file() else '')+'<a href="#documents">文档</a><a href="#protocol">执行合同</a><a href="GALLERY.html">全部图集</a><a href="GALLERY.md">GitHub图集</a></nav>'
    (root/'index.html').write_text(page(summary_title or 'MAX-STRUCT R13 · 报告与证据导航',body,nav))
    gallery='<p class="note">此页只汇总真实分析目录。图片采用原文件；科学解释请读取各组报告。虚拟验证夹具不在本图集中。</p>'
    md=['# MAX-STRUCT R13 实际图集','','仅包含analysis/下真实生成图；虚拟夹具不进入本页。','']
    for c in campaigns:
        gallery+='<article class="card"><h2>'+esc(c['name'])+'</h2><p><a href="'+link(c['path']/'GALLERY.html',root)+'">打开完整图集与数据附件</a> · '+str(len(c['figures']))+'张实际图</p></article>'
        md+=['## '+c['name'],'','[单组图集]('+link(c['path']/'GALLERY.md',root)+')','']
        for f in c['figures']:md+=['- ['+f.stem+']('+link(f,root)+')']
        md+=['']
    (root/'GALLERY.html').write_text(page('MAX-STRUCT R13 · 全部实际图集',gallery,'<nav><a href="index.html">首页</a><a href="GALLERY.md">GitHub Markdown图集</a></nav>'));(root/'GALLERY.md').write_text('\n'.join(md)+'\n')
    record=dict(status='LOCAL_BUILD_ONLY_NOT_PUBLISHED',campaigns=[dict(name=c['name'],group=c['group'],figures=len(c['figures']),domains=c['domains']) for c in campaigns],figure_count=sum(len(c['figures']) for c in campaigns),documents=[dict(path=link(path,root),title=heading(path,label)) for path,label in docs],download_manifests=[link(path,root) for path in downloads],fake_fixture_images_included=False,truth_or_numeric_results_read=False,parent_authored_summary_loaded=(root/'WEB_SUMMARY.json').is_file(),curated_summary_figures=summary_figures,wall_seconds=time.perf_counter()-started)
    (root/'SITE_BUILD.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:record[k] for k in ['status','figure_count','fake_fixture_images_included','wall_seconds']}))
if __name__=='__main__':main()
