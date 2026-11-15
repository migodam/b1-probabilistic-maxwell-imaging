"""Render existing R14 prose/figures as static Pages; no science or publication.

Only reads public-facing Markdown, figure names, selected summary CSVs and
DOWNLOADS metadata. Never traverses sealed/private trees or follows symlinks.
"""
from __future__ import annotations
import argparse,csv,html,json,os,re,time,subprocess
from datetime import datetime,timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote,unquote,urlsplit

ROOT=Path(__file__).resolve().parents[1]
DOCS=[('RESEARCH_REPORT_ZH.md','完整研究报告'),('THEORY_DERIVATIONS.md','公式、推导与边界'),
 ('TRAJECTORY_ANALYSIS.md','轨迹分析'),('BENCHMARKS.md','比较实验'),('COST_LEDGER.md','费用总账'),
 ('FAILURES.md','失败与范围'),('REPRODUCE.md','复现说明'),('PRO_HANDOFF.md','Pro研究交接'),
 ('CHAT_HANDOFF.md','研究交接'),('PAPER_DRAFT.md','论文草稿'),
 ('EXPERIMENT_CONTRACT.md','实验合同'),('INDEPENDENT_CONTRACT.md','独立确认合同'),
 ('CONTROL_DERIVATION.md','控制推导'),('CORE_IMPLEMENTATION_NOTES.md','实现说明'),
 ('BASELINE_AUDIT.md','基线审计'),('THEORY_STATUS.md','理论状态'),
 ('PRIOR_ART_AND_NOVELTY.md','先例与新颖性'),('METHODS_AND_ABLATIONS.md','方法与消融'),
 ('TECH_STACK_MAP.md','技术栈与证据范围'),('PRO_PACKET_CONTENTS.md','Pro数据包说明'),
 ('DOWNLOAD_GUIDE.md','下载与恢复指南')]
GROUPS=[('exposed48','已暴露48对象'),('independent24','独立24对象'),('controls','固定控制实验'),
 ('vector3d','三维矢量实验'),('fullstate_root','完整状态干预与局部根'),('other','其他实际分析图')]
PRIVATE=('sealed_truth','private_operational','chat_readback','secrets','credentials')
IMAGE_EXT={'.png','.jpg','.jpeg','.svg','.webp','.gif'}
SUMMARY_NAMES={'METHOD_VIEW_SUMMARY.csv','METHOD_SUMMARY.csv','POLICY_SUMMARY.csv','SUMMARY.csv','COST_LEDGER.csv'}
GIT_IGNORED=set()
GIT_IGNORE_STATUS='NOT_CHECKED'
CSS='''*{box-sizing:border-box}body{margin:0;background:#f3f5f8;color:#17253a;font:16px/1.75 -apple-system,BlinkMacSystemFont,"PingFang SC","Microsoft YaHei",sans-serif}main{max-width:1280px;margin:auto;padding:30px 24px 70px}header{border-bottom:1px solid #d9e0e8;padding:20px 0;margin-bottom:25px}h1{font-size:clamp(27px,4vw,42px);line-height:1.3}h2{font-size:25px;margin-top:36px}h3{font-size:19px}a{color:#165a89;text-decoration:none}a:hover{text-decoration:underline}nav,.doclist{display:flex;flex-wrap:wrap;gap:10px 20px}nav{margin:15px 0}.doclist a{background:white;border:1px solid #dce3eb;border-radius:8px;padding:7px 12px}article{max-width:1100px}article img{max-width:100%;height:auto}.note,footer{color:#526479}.empty{border-left:4px solid #9eb2c8;background:#fff;padding:16px}pre{white-space:pre-wrap;overflow:auto;background:#fff;padding:16px}code{font-size:.9em;background:#e9eff6;padding:1px 4px;border-radius:3px;overflow-wrap:anywhere}.math{white-space:pre-wrap}.math-display{overflow-x:auto;padding:10px 0}mjx-container{max-width:100%}.comparison{overflow:auto;margin:18px 0}table{border-collapse:collapse;width:100%;background:white;font-size:14px}th,td{padding:8px 12px;text-align:left;border:1px solid #dce3eb}th{background:#e9eff6;white-space:nowrap}td{font-variant-numeric:tabular-nums}.gallery{display:grid;grid-template-columns:1fr;gap:24px}figure{margin:0;background:white;border:1px solid #dce3eb;border-radius:10px;padding:12px}figure img{display:block;width:100%;height:auto}figcaption{padding:10px 3px 0;font-size:14px;overflow-wrap:anywhere}.path{display:block;font-size:12px;color:#69798a}details{background:#fff;border:1px solid #dce3eb;border-radius:8px;padding:12px;margin:14px 0}summary{cursor:pointer;font-weight:600}footer{margin-top:40px;border-top:1px solid #dce3eb;padding-top:20px;font-size:13px}.cards{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:15px}.card{background:#fff;border:1px solid #dce3eb;border-radius:10px;padding:16px}@media(max-width:760px){main{padding:18px 12px 45px}.cards{grid-template-columns:1fr}figure{padding:6px}nav{gap:8px 12px}}'''


MATHJAX_HEAD=r"""<script>
window.MathJax={tex:{inlineMath:[['\\(','\\)'],['$','$']],displayMath:[['\\[','\\]'],['$$','$$']],processEscapes:true},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};
</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js"></script>"""


def esc(x):return html.escape(str(x),quote=True)

def private(path):return any(x in unquote(str(path)).lower() for x in PRIVATE)

def safe_file(path,root):
    # Prevent symlink aliases from pulling a virtual D: or private target into
    # the report. A permitted file must resolve beneath the actual study root.
    try:
        if path.resolve() in GIT_IGNORED:return False
        rel=path.relative_to(root)
        if any((root/Path(*rel.parts[:i])).is_symlink() for i in range(1,len(rel.parts)+1)):return False
        return not private(rel) and path.is_file() and path.resolve().is_relative_to(root.resolve())
    except (ValueError,OSError):return False


def relurl(path,base):return quote(os.path.relpath(path,base).replace(os.sep,'/'),safe='/_.-')


def walk_public(directory,root):
    if not directory.is_dir() or directory.is_symlink():return
    for base,dirs,names in os.walk(directory,followlinks=False):
        here=Path(base)
        dirs[:]=[d for d in dirs if not private(here/d) and not (here/d).is_symlink()]
        for name in sorted(names):
            path=here/name
            if safe_file(path,root):yield path


class Renderer:
    def __init__(self,root,documents):self.root=root;self.documents=documents;self.suppressed=[]
    def target(self,url,source_base,output_base,image=False):
        u=urlsplit(url)
        if private(url):self.suppressed.append(dict(url=url,reason='PRIVATE_OR_SEALED'));return None
        if u.scheme in ('http','https') and u.netloc:return url
        if u.scheme or u.netloc:return None
        if not u.path:return '#'+u.fragment if u.fragment else None
        p=(source_base/unquote(u.path)).resolve()
        if not safe_file(p,self.root):
            self.suppressed.append(dict(url=url,reason='MISSING_OR_OUTSIDE_PUBLIC_STUDY'));return None
        if not image and p in self.documents:p=self.documents[p]
        result=relurl(p,output_base)
        return result+('#'+u.fragment if u.fragment else '')
    def inline(self,text,source_base,output_base):
        # Raw HTML is escaped; the only link/image elements are generated here.
        pattern=r'(!?\[[^\]]*\]\([^)]+\)|`[^`]+`|\*\*[^*]+\*\*|\\\(.+?\\\))'
        parts=[]
        for token in re.split(pattern,text):
            m=re.fullmatch(r'(!?)\[([^\]]*)\]\(([^)]+)\)',token)
            if m:
                marker,label,url=m.groups();link=self.target(url,source_base,output_base,bool(marker))
                if link:
                    parts.append(('<img loading="lazy" src="'+esc(link)+'" alt="'+esc(label)+'">') if marker else '<a href="'+esc(link)+'">'+esc(label)+'</a>')
                else:parts.append(esc(label)+' <span class="note">（链接未公开或文件尚未生成）</span>')
            elif token.startswith('`') and token.endswith('`'):parts.append('<code>'+esc(token[1:-1])+'</code>')
            elif token.startswith('**') and token.endswith('**'):parts.append('<strong>'+esc(token[2:-2])+'</strong>')
            elif token.startswith('\\('):parts.append('<span class="math">'+esc(token)+'</span>')
            else:parts.append(esc(token))
        return ''.join(parts)
    def markdown(self,text,source_base,output_base):
        lines=text.splitlines();out=[];i=0
        render=lambda x:self.inline(x,source_base,output_base)
        while i<len(lines):
            line=lines[i].strip()
            if not line:i+=1;continue
            if line.startswith('```') or line in ('\\[','$$'):
                ending='```' if line.startswith('```') else ('\\]' if line=='\\[' else '$$');block=[];i+=1
                while i<len(lines) and not lines[i].strip().startswith(ending):block.append(lines[i]);i+=1
                if ending=='```':out.append('<pre><code>'+esc('\n'.join(block))+'</code></pre>')
                else:
                    opening='\\[' if ending=='\\]' else '$$'
                    out.append('<div class="math math-display">'+esc(opening+'\n'+'\n'.join(block)+'\n'+ending)+'</div>')
                i+=1;continue
            m=re.match(r'^(#{1,6})\s+(.+)$',line)
            if m:out.append(f'<h{len(m[1])}>'+render(m[2])+f'</h{len(m[1])}>');i+=1;continue
            if line.startswith('|') and i+1<len(lines) and re.fullmatch(r'[\s:|\-]+',lines[i+1]):
                headers=[s.strip() for s in line.strip('|').split('|')];data=[];i+=2
                while i<len(lines) and lines[i].strip().startswith('|'):
                    data.append([s.strip() for s in lines[i].strip().strip('|').split('|')]);i+=1
                out.append(table(headers,data,render));continue
            if re.match(r'^(?:[-*]|\d+\.)\s+',line):
                tag='ol' if re.match(r'^\d+\.',line) else 'ul';items=[]
                while i<len(lines) and re.match(r'^(?:[-*]|\d+\.)\s+',lines[i].strip()):
                    items.append(re.sub(r'^(?:[-*]|\d+\.)\s+','',lines[i].strip()));i+=1
                out.append('<'+tag+'>'+''.join('<li>'+render(x)+'</li>' for x in items)+'</'+tag+'>');continue
            if line in ('---','***'):out.append('<hr>');i+=1;continue
            para=[line];i+=1
            while i<len(lines) and lines[i].strip() and not re.match(r'^(?:#|\||[-*]\s|\d+\.\s|```|\\\[|\$\$)',lines[i].strip()):para.append(lines[i].strip());i+=1
            out.append('<p>'+render(' '.join(para))+'</p>')
        return '\n'.join(out)


def table(headers,rows,render=esc):
    return '<div class="comparison"><table><thead><tr>'+''.join('<th scope="col">'+render(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+render(x)+'</td>' for x in r)+'</tr>' for r in rows)+'</tbody></table></div>'


def page(title,body,path,root,docs):
    nav=[(root/'index.html','研究报告'),(root/'galleries.html','全部图集'),(root/'evidence.html','数据与代码')]
    nav += [(target,label) for source,label,target in docs if source.name!='RESEARCH_REPORT_ZH.md']
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+esc(title)+'</title><style>'+CSS+'</style>'+MATHJAX_HEAD+'</head><body><main><header><h1>'+esc(title)+'</h1><nav>'+''.join('<a href="'+relurl(p,path.parent)+'">'+esc(label)+'</a>' for p,label in nav)+'</nav></header>'+body+'<footer>MAX-STRUCT R14 · 此网页呈现已存在的报告、图和文件；实验状态、结果与边界以原始报告及回执为准。网页构建不新增科学结果。</footer></main></body></html>\n'


def group_for(path):
    s=str(path).lower()
    if any(t in s for t in ('vector','3d')):return 'vector3d'
    if any(t in s for t in ('fullstate','interven','local_closed','root','checkpoint')):return 'fullstate_root'
    if 'independent' in s:return 'independent24'
    if any(t in s for t in ('self_only','control','reset','schedule')):return 'controls'
    if any(t in s for t in ('dualmetric48','exposed','pilot12','remaining36','comparison48')):return 'exposed48'
    return 'other'


def figure_label(path):
    ids=list(dict.fromkeys(re.findall(r'(?:r13_fresh_\d{4}|r14_new_\d{4}|[A-Za-z][A-Za-z0-9]*_\d{4})',path.stem)))
    return (' / '.join(ids)+' · ' if ids else '')+path.stem


def gallery_figure(path,target,root,metadata=None):
    label=figure_label(path);href=relurl(path,target.parent);pdf=path.with_suffix('.pdf')
    metadata=metadata or {};case_label=metadata.get('case_ids',metadata.get('case_id',''));scope=metadata.get('scope','')
    if case_label:label=str(case_label)+' · '+label
    manifest=path.parent.parent/'FIGURE_MANIFEST.csv'
    source_link=' · <a href="'+relurl(manifest,target.parent)+'">原图清单</a>' if safe_file(manifest,root) else ''
    extra=' · <a href="'+relurl(pdf,target.parent)+'">PDF</a>' if safe_file(pdf,root) else ''
    return '<figure>'+('<a href="'+href+'"><img loading="lazy" decoding="async" src="'+href+'" alt="'+esc(label)+'"></a>' if path.suffix.lower() in IMAGE_EXT else '<p><a href="'+href+'">打开PDF图</a></p>')+'<figcaption><strong>'+esc(label)+'</strong> · <a href="'+href+'">原图</a>'+extra+source_link+('<span class="path">'+esc(scope)+'</span>' if scope else '')+'<span class="path">'+esc(path.relative_to(root))+'</span></figcaption></figure>'


def downloads(root):
    path=root/'DOWNLOADS.json';links=[];skipped=[]
    if not safe_file(path,root):return '<p class="empty">发布下载清单尚未生成。</p>',links,skipped
    data=json.loads(path.read_text())
    global_verified=(data.get('status')=='PUBLISHED_AND_ANONYMOUSLY_READ_BACK' or data.get('verified') is True)
    for item in data.get('assets',data.get('downloads',[])):
        if not isinstance(item,dict):continue
        url=item.get('url',item.get('download_url',''));u=urlsplit(url)
        verified=(global_verified or item.get('verified') is True or item.get('verification_status') in ('PASS','VERIFIED','ANONYMOUSLY_READ_BACK'))
        if u.scheme!='https' or '/releases/download/' not in u.path or not verified or private(url):
            skipped.append(dict(name=item.get('name',item.get('filename','')),reason='Not a declared verified public release asset'));continue
        label=item.get('name',item.get('filename',u.path.rsplit('/',1)[-1]));links.append(dict(name=label,url=url))
    body='<p class="note">以下仅列DOWNLOADS.json中已标明公开核验的实际release附件；网页不联网验证其内容。</p>'
    body+='<div class="doclist">'+''.join('<a href="'+esc(x['url'])+'">'+esc(x['name'])+'</a>' for x in links)+'</div>'
    if not links:body+='<p class="empty">尚无已核验的公开release附件，未生成占位URL。</p>'
    return body,links,skipped


class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        self.links += [v for k,v in attrs if k in ('href','src') and v]


def audit_links(pages,root):
    broken=[];forbidden=[];checked=0;external=set()
    for path in pages:
        parser=Links();parser.feed(path.read_text())
        for url in parser.links:
            checked+=1;u=urlsplit(url)
            if private(url):forbidden.append(dict(page=str(path.relative_to(root)),url=url));continue
            if u.scheme in ('http','https'):external.add(url);continue
            if u.scheme or u.netloc:forbidden.append(dict(page=str(path.relative_to(root)),url=url));continue
            if u.path and not safe_file((path.parent/unquote(u.path)).resolve(),root):broken.append(dict(page=str(path.relative_to(root)),url=url))
    return dict(checked_href_src=checked,broken_local_links=broken,forbidden_private_or_outside_links=forbidden,
        external_URLs=sorted(external),external_network_verification='NOT_RUN_BY_SITE_BUILDER')


def build(root):
    global GIT_IGNORED,GIT_IGNORE_STATUS
    start=time.perf_counter();cpu=time.process_time();generated=[];sources=[];csv_reads=[]
    # Local Release-only tables must not become broken Pages download links.
    # This only asks Git about names; it never fetches historical blobs.
    candidates=[p for name in ('src','configs','analysis','audit','FIGURES','sources') for p in walk_public(root/name,root)]
    candidates.extend(p for p in root.iterdir() if p.is_file() and not private(p))
    checked=subprocess.run(['git','-C',str(root),'check-ignore','--no-index','-z','--stdin'],
        input=b''.join(os.fsencode(p)+b'\0' for p in candidates),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    GIT_IGNORED={Path(os.fsdecode(p)).resolve() for p in checked.stdout.split(b'\0') if p} if checked.returncode in (0,1) else set()
    GIT_IGNORE_STATUS='APPLIED_TO_PUBLIC_NAMES' if checked.returncode in (0,1) else 'STANDALONE_NOT_A_GIT_TREE'
    docs=[]
    for name,label in DOCS:
        source=root/name
        if safe_file(source,root):docs.append((source,label,root/('index.html' if name=='RESEARCH_REPORT_ZH.md' else name.removesuffix('.md').lower()+'.html')))
    renderer=Renderer(root,{source.resolve():target for source,_,target in docs})
    def write(target,title,body):
        target.parent.mkdir(parents=True,exist_ok=True);target.write_text(page(title,body,target,root,docs));generated.append(target)
    for source,label,target in docs:
        body='<p><a href="'+relurl(source,target.parent)+'">下载原始Markdown</a></p><article>'+renderer.markdown(source.read_text(),source.parent,target.parent)+'</article>'
        if source.name=='RESEARCH_REPORT_ZH.md':
            dl,dl_urls,dl_skip=downloads(root);body+='<h2>完整图集与可下载材料</h2><p><a href="galleries.html">查看全部已生成图</a> · <a href="evidence.html">原始表格、报告与代码</a></p><h2>已核验发布附件</h2>'+dl
        write(target,'MAX-STRUCT R14｜'+label,body);sources.append(str(source.relative_to(root)))
    if not any(s.name=='RESEARCH_REPORT_ZH.md' for s,_,_ in docs):
        write(root/'index.html','MAX-STRUCT R14','<p class="empty">完整研究报告尚未生成。下列入口只展示实际存在的材料，不表示实验已完成。</p><p><a href="galleries.html">已有图集</a> · <a href="evidence.html">已有文档与代码</a></p>')
    discovered=[];archived_materials=[];excluded_v1_figures=[]
    for base in (root/'analysis',root/'audit',root/'FIGURES'):
        for path in walk_public(base,root):
            if path.suffix.lower() in IMAGE_EXT|{'.pdf'} and (base.name in ('audit','FIGURES') or 'FIGURES' in path.relative_to(base).parts):
                if path.suffix.lower()=='.pdf' and any(safe_file(path.with_suffix(ext),root) for ext in IMAGE_EXT):continue
                relative=path.relative_to(root).as_posix()
                if relative.startswith('analysis/exposed48_v1/'):
                    excluded_v1_figures.append(relative);continue
                if relative.startswith('analysis/exposed48_v2/FIGURES/') and (re.fullmatch(r'materials_\d+',path.stem) or path.stem=='all48_materials'):
                    archived_materials.append(path);continue
                discovered.append(path)
    discovered=sorted(set(discovered),key=lambda p:(not p.relative_to(root).as_posix().startswith('FIGURES/exposed48/'),p.suffix.lower() not in IMAGE_EXT,str(p)));groups={key:[] for key,_ in GROUPS}
    figure_meta={}
    for manifest in sorted({p.parent.parent/'FIGURE_MANIFEST.csv' for p in discovered}):
        if not safe_file(manifest,root):continue
        with manifest.open(newline='') as stream:
            for row in csv.DictReader(stream):
                name=Path(row.get('png',row.get('file',row.get('path','')))).name
                if name:figure_meta[(manifest.parent/'FIGURES'/name).resolve()]=row
    for p in discovered:groups[group_for(p.relative_to(root))].append(p)
    cards=[];inventory=[]
    for key,label in GROUPS:
        files=groups[key]
        if not files:continue
        target=root/'galleries'/f'{key}.html';sections={}
        for f in files:sections.setdefault(str(f.parent.relative_to(root)),[]).append(f)
        body='<p>'+str(len(files))+'项实际图像／PDF。完整文件名和对象ID保留；图中方法、视图与指标以原图为准。</p>'
        for folder,items in sections.items():
            body+='<h2>'+esc(folder)+'</h2><div class="gallery">'+''.join(gallery_figure(f,target,root,figure_meta.get(f.resolve())) for f in items)+'</div>'
        write(target,label,body);cards.append('<article class="card"><h3><a href="'+relurl(target,root)+'">'+esc(label)+'</a></h3><p>'+str(len(files))+'张实际图</p></article>')
        inventory.append(dict(group=key,title=label,page=str(target.relative_to(root)),figures=[str(f.relative_to(root)) for f in files]))
    write(root/'galleries.html','R14｜全部图集','<p>从实际analysis/FIGURES、audit图片及根FIGURES发现；不挑选正例，不根据材料评分筛图。</p><div class="cards">'+''.join(cards)+'</div>' if cards else '<p class="empty">尚无实际图文件。</p>')
    public_files=[]
    for folder in ('src','configs','analysis','audit'):
        for f in walk_public(root/folder,root):
            # No result/state/input arrays, no operational machine/session logs.
            # Release bundles can expose the full approved closure later.
            if f.suffix.lower() not in ('.md','.csv','.py'):continue
            if f.relative_to(root).as_posix().startswith('analysis/exposed48_v1/'):continue
            if any(x in f.name.lower() for x in ('truth_source','credential','private','chat_readback')):continue
            public_files.append(f)
    body='<h2>报告与推导原文</h2><div class="doclist">'+''.join('<a href="'+relurl(s,root)+'">'+esc(s.name)+'</a>' for s,_,_ in docs)+'</div>'
    body+='<h2>实际表格与代码</h2><p>原文件下载；未生成科学结论。完整原始状态及数据包通过审核后的release附件提供。</p>'
    for folder in ('src','configs','analysis','audit'):
        fs=[f for f in public_files if f.relative_to(root).parts[0]==folder]
        if fs:body+='<details><summary>'+folder+' · '+str(len(fs))+'个文件</summary><ul>'+''.join('<li><a href="'+relurl(f,root)+'">'+esc(f.relative_to(root))+'</a></li>' for f in fs)+'</ul></details>'
    for f in public_files:
        if f.name not in SUMMARY_NAMES or f.relative_to(root).parts[0]!='analysis':continue
        with f.open(newline='') as stream:
            reader=csv.reader(stream);rows=[]
            for i,row in enumerate(reader):
                if i>=61:break
                rows.append(row)
        csv_reads.append(dict(file=str(f.relative_to(root)),displayed_rows=max(0,len(rows)-1),preview_cap=60))
        if rows:body+='<details><summary>'+esc(f.relative_to(root))+'（最多60行预览）</summary>'+table(rows[0],rows[1:])+'</details>'
    dl,dl_urls,dl_skip=downloads(root);body+='<h2>已核验发布附件</h2>'+dl
    if archived_materials:
        archive_page=root/'galleries/exposed48_materials_archive.html'
        write(archive_page,'R14｜旧版暴露48材料图归档','<p class="note">原analysis/exposed48_v2材料图仅作归档；主画廊使用FIGURES/exposed48的可读色条版本。同一材料结果不在主画廊重复展示。</p><div class="gallery">'+''.join(gallery_figure(f,archive_page,root) for f in sorted(archived_materials))+'</div>')
        body+='<h2>历史图归档</h2><p><a href="'+relurl(archive_page,root)+'">旧版暴露48材料图（归档）</a></p>'
    write(root/'evidence.html','R14｜数据、代码与下载',body)
    linkcheck=audit_links(generated,root)
    receipt=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='COMPLETE_LOCAL_BUILD' if not linkcheck['broken_local_links'] and not linkcheck['forbidden_private_or_outside_links'] else 'LOCAL_LINK_AUDIT_FAILED',
        generated_pages=[str(p.relative_to(root)) for p in generated],report_exists=(root/'RESEARCH_REPORT_ZH.md').is_file(),
        document_sources=sources,missing_optional_documents=[name for name,_ in DOCS if not safe_file(root/name,root)],
        gallery_inventory=inventory,gallery_figures=len(discovered),archived_material_figures=[str(p.relative_to(root)) for p in archived_materials],excluded_bugged_v1_figures=excluded_v1_figures,summary_CSV_previews=csv_reads,math_renderer='MathJax3.2.2CDN;originaldelimitersretainedwithoutJS',
        approved_download_URLs=dl_urls,skipped_download_entries=dl_skip,suppressed_markdown_links=renderer.suppressed,
        link_audit=linkcheck,git_ignore_status=GIT_IGNORE_STATUS,release_only_names_omitted=len(GIT_IGNORED),private_trees_scanned=False,symlinks_followed=False,physics_runs=0,scoring_runs=0,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,published=False,
        visual_QA='NOT_RUN:parentmustinspectactualpagesandimagesbeforepublication',
        source_scope='Presentation-only existingprose/files;noinferredcompletionormaterialverdict')
    (root/'SITE_BUILD_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    a=p.parse_args();receipt=build(a.root.resolve())
    print(json.dumps({k:receipt[k] for k in ('status','generated_pages','gallery_figures','published')},ensure_ascii=False))
    return 0 if receipt['status']=='COMPLETE_LOCAL_BUILD' else 1


if __name__=='__main__':raise SystemExit(main())
