"""Own R3 evidence export using the existing R1 minimal Markdown renderer."""
from pathlib import Path
import sys,json,time,re,zipfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'research_max_struct_r1/experiments'))
from export_delivery import markdown_subset

def main():
 tick,cpu=time.perf_counter(),time.process_time()
 previous=ROOT/'audit/EXPORT_RECEIPT.json'
 if previous.exists():
  history=ROOT/'audit/export_history';history.mkdir(exist_ok=True)
  (history/f'export_{len(list(history.glob("export_*.json")))+1:02d}.json').write_text(previous.read_text())
 raw='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r3/'
 body=markdown_subset((ROOT/'RESEARCH_REPORT_ZH.md').read_text())
 body=re.sub(r'href="([^":]+\.md)"',lambda m:'href="'+raw+m[1]+'"',body)
 css='body{margin:0;background:#f3f5f7;color:#172230;font:17px/1.8 system-ui,sans-serif}main{max-width:1100px;margin:32px auto;background:white;padding:40px;border-radius:12px}h1{font-size:30px;line-height:1.4}h2{margin-top:40px;border-top:1px solid #dde3e8;padding-top:18px}a{color:#075b93}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border:1px solid #d9e1e8;padding:9px;text-align:left}th{background:#eaf1f6}.table-scroll{overflow:auto}img{max-width:100%;height:auto}figure{margin:24px 0}code{background:#eef2f4;font-size:14px;overflow-wrap:anywhere}@media(max-width:700px){main{margin:0;padding:20px}h1{font-size:25px}}'
 page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>MAX-STRUCT-R3 中文研究报告</title><style>'+css+'</style><main><nav><a href="'+raw+'CHAT_HANDOFF.md">给Chat的入口</a> · <a href="'+raw+'RESEARCH_REPORT_ZH.md">原始报告</a> · <a href="'+raw+'THEORY_DERIVATIONS.md">数学推导</a></nav>'+body+'</main></html>'
 (ROOT/'index.html').write_text(page)
 missing=[u for u in re.findall(r'(?:src|href)="([^"#]+)"',page) if '://' not in u and not(ROOT/u).exists()]
 if missing:raise ValueError(missing)
 files=[p for p in ROOT.rglob('*')if p.is_file()and '__pycache__'not in p.parts and 'local_only'not in p.parts and p.suffix not in ('.pyc','.zip')and p.name!='.DS_Store']
 dest=ROOT.parents[2]/'r3_release';dest.mkdir(exist_ok=True);archive=dest/'MAX_STRUCT_R3_reproducibility.zip'
 receipt=dict(scope='Own HTML/ZIP packaging only, no new physics',physical_LU=0,physical_RHS=0,files=len(files),excludes=['third-party local_only fulltext','Python caches'],complete=False)
 (ROOT/'audit/EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
 if ROOT/'audit/EXPORT_RECEIPT.json'not in files:files.append(ROOT/'audit/EXPORT_RECEIPT.json')
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED)as z:
  for p in files:z.write(p,'research_max_struct_r3/'+str(p.relative_to(ROOT)))
 receipt.update(complete=True,archive_bytes=archive.stat().st_size,wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu,archive_receipt_note='ZIP holds pre-completion receipt; outer final receipt authoritative')
 (ROOT/'audit/EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
if __name__=='__main__':main()
