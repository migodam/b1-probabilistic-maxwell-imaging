"""Reuse R1 report renderer; export own artifacts, excluding third-party fulltext."""
from pathlib import Path
import sys,time,json,zipfile,re
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'research_max_struct_r1/experiments'))
from export_delivery import markdown_subset
def main():
    tick=time.perf_counter();cpu=time.process_time()
    body=markdown_subset((ROOT/'RESEARCH_REPORT_ZH.md').read_text())
    # GitHub renders .md; Pages does not generally render raw .md links in this theme.
    raw='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r2/'
    body=re.sub(r'href="([^":]+\.md)"',lambda m:'href="'+raw+m[1]+'"',body)
    css='body{margin:0;background:#f3f5f7;color:#172230;font:17px/1.8 system-ui,sans-serif}main{max-width:1100px;margin:32px auto;background:white;padding:40px;border-radius:12px}h1{font-size:30px;line-height:1.4}h2{margin-top:40px;border-top:1px solid #dde3e8;padding-top:18px}a{color:#075b93}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border:1px solid #d9e1e8;padding:9px;text-align:left}th{background:#eaf1f6}.table-scroll{overflow:auto}img{max-width:100%;height:auto}figure{margin:24px 0}code{background:#eef2f4;font-size:14px;overflow-wrap:anywhere}@media(max-width:700px){main{margin:0;padding:20px}h1{font-size:25px}}'
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>MAX-STRUCT-R2 中文研究报告</title><style>'+css+'</style><main><nav><a href="'+raw+'CHAT_HANDOFF.md">给Chat的入口</a> · <a href="'+raw+'RESEARCH_REPORT_ZH.md">原始报告</a></nav>'+body+'</main></html>'
    (ROOT/'index.html').write_text(page)
    files=[]
    for p in ROOT.rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts or p.suffix in ('.zip','.pyc'):continue
        if 'literature/receipts/' in str(p.relative_to(ROOT)) and p.suffix in ('.pdf','.txt'):continue
        files.append(p)
    missing=[x for x in re.findall(r'(?:src|href)="([^"#]+)"',page) if '://' not in x and not (ROOT/x).exists()]
    if missing:raise ValueError(missing)
    receipt=dict(scope='HTML/own artifact ZIP; no new physics',LU=0,RHS=0,files=len(files),
         excludes='third-party fulltextPDF/txt; caches',wall_s_before_zip=time.perf_counter()-tick)
    (ROOT/'audit/EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    dest=ROOT.parents[2]/'r2_release';dest.mkdir(exist_ok=True)
    archive=dest/'MAX_STRUCT_R2_reproducibility.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,'research_max_struct_r2/'+str(p.relative_to(ROOT)))
        if ROOT/'audit/EXPORT_RECEIPT.json' not in files:z.write(ROOT/'audit/EXPORT_RECEIPT.json','research_max_struct_r2/audit/EXPORT_RECEIPT.json')
    receipt.update(wall_s=time.perf_counter()-tick,cpu_s=time.process_time()-cpu,archive_bytes=archive.stat().st_size,
        recursive_timing_note='ZIP contains pre-completion exportreceipt; outer final receipt authoritative')
    (ROOT/'audit/EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
if __name__=='__main__':main()
