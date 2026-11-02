"""Local report/portable ZIP export. No publishing, no hashes, no physics."""
from pathlib import Path
import html,re,json,time,zipfile
ROOT=Path(__file__).resolve().parents[1]

def inline(s):
    s=html.escape(s)
    s=re.sub(r'!\[([^\]]*)\]\(([^)]+)\)',r'<figure><img alt="\1" src="\2"><figcaption>\1</figcaption></figure>',s)
    s=re.sub(r'\[([^\]]*)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
    s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s)
    return s

def markdown_subset(text):
    out=[];lines=text.splitlines();i=0;paragraph=[]
    def flush():
        if paragraph:out.append('<p>'+inline(' '.join(paragraph))+'</p>');paragraph.clear()
    while i<len(lines):
        s=lines[i]
        if not s.strip():flush();i+=1;continue
        if s.startswith('#'):
            flush();n=len(s)-len(s.lstrip('#'));out.append(f'<h{n}>'+inline(s[n:].strip())+f'</h{n}>')
        elif s.startswith('|'):
            flush();table=[]
            while i<len(lines) and lines[i].startswith('|'):table.append(lines[i]);i+=1
            out.append('<div class="table-scroll"><table>')
            for j,row in enumerate(table):
                cells=[x.strip() for x in row.strip('|').split('|')]
                if all(re.fullmatch(r'[-: ]+',c) for c in cells):continue
                tag='th' if j==0 else 'td'
                out.append('<tr>'+''.join('<'+tag+'>'+inline(c)+'</'+tag+'>' for c in cells)+'</tr>')
            out.append('</table></div>');continue
        elif s.startswith('- '):
            flush();out.append('<ul>')
            while i<len(lines) and lines[i].startswith('- '):out.append('<li>'+inline(lines[i][2:])+'</li>');i+=1
            out.append('</ul>');continue
        else:paragraph.append(s)
        i+=1
    flush();return '\n'.join(out)

def main():
    tic=time.perf_counter();cpu=time.process_time()
    body=markdown_subset((ROOT/'RESEARCH_REPORT_ZH.md').read_text())
    css='''body{margin:0;background:#f3f5f7;color:#172230;font:17px/1.8 system-ui,-apple-system,sans-serif}main{max-width:1080px;margin:32px auto;background:white;padding:48px 56px;border-radius:12px}h1{font-size:32px;line-height:1.4}h2{margin-top:48px;border-top:1px solid #d9e1e8;padding-top:20px;font-size:25px}h3{font-size:20px;margin-top:32px}a{color:#075b93}code{font:14px/1.6 ui-monospace,monospace;background:#eef2f4;padding:2px 5px;overflow-wrap:anywhere}table{width:100%;border-collapse:collapse;font-size:15px}th,td{border:1px solid #d9e1e8;padding:10px;text-align:left;vertical-align:top}th{background:#eaf1f6}.table-scroll{overflow-x:auto}img{max-width:100%;height:auto}figure{margin:24px 0}figcaption{font-size:14px;color:#516476}nav{font-size:14px;color:#516476}.badge{display:inline-block;padding:4px 12px;border-radius:6px;background:#f3e8cb}@media(max-width:700px){main{margin:0;padding:24px 18px}body{font-size:16px}h1{font-size:26px}}'''
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MAX-STRUCT-R1 中文研究报告（本地）</title><style>'+css+'</style><main><nav><span class="badge">本地封存 · 未发布</span>　<a href="README.md">文件索引</a> · <a href="RESEARCH_REPORT_ZH.md">原始 Markdown</a> · <a href="THEORY_DERIVATIONS.md">数学推导</a> · <a href="REPRODUCE.md">复现</a></nav>'+body+'</main></html>'
    (ROOT/'REPORT.html').write_text(page)
    # Validate authored local links/assets; no network or solver.
    missing=[]
    for target in re.findall(r'(?:src|href)="([^"]+)"',page):
        if '://' not in target and not (ROOT/target.split('#')[0]).exists():missing.append(target)
    if missing:raise RuntimeError('Missing local report links: '+str(missing))
    def included(p):
        rel=p.relative_to(ROOT)
        return p.is_file() and '__pycache__' not in rel.parts and p.suffix not in ('.zip','.pyc') and not ('audit/literature/' in str(rel) and p.suffix in ('.pdf','.txt'))
    files=sorted(p for p in ROOT.rglob('*') if included(p))
    manifest=dict(scope='local delivery; no public release',base_git_sha='7317601b43b914525c21f09e5a60020552ea677b',
                  exclusions=['third-party fulltext PDF/txt in audit/literature','Python caches','ZIP itself'],hashes='not computed',
                  files=[str(p.relative_to(ROOT)) for p in files])
    (ROOT/'audit/DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
    receipt=dict(operation='local HTML + ZIP export',new_LU=0,new_RHS=0,missing_local_links=[],publication='NOT_PUBLISHED',
                 wall_s_before_zip=time.perf_counter()-tic,cpu_s_before_zip=time.process_time()-cpu)
    (ROOT/'audit/EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    archive=ROOT/'MAX_STRUCT_R1_delivery.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(ROOT.rglob('*')):
            if included(p):z.write(p,'research_max_struct_r1/'+str(p.relative_to(ROOT)))
    receipt.update(wall_s_including_zip=time.perf_counter()-tic,cpu_s_including_zip=time.process_time()-cpu,archive_bytes=archive.stat().st_size,
                   note='post-ZIP completion timing is in local receipt; archived receipt records pre-ZIP stage to avoid recursive self-packaging')
    (ROOT/'audit/EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt))

if __name__=='__main__':main()
