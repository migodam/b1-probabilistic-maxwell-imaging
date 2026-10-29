"""Anonymous streaming readback of already authorized public research assets."""
from pathlib import Path
import concurrent.futures,datetime,hashlib,json,time,urllib.request
P=Path(__file__).resolve().parents[1]
manifest=json.loads((P/'RELEASE_MANIFEST.json').read_text())
def asset(a):
 t=time.perf_counter();h=hashlib.blake2b(digest_size=32);n=0
 try:
  req=urllib.request.Request(a['url'],headers={'User-Agent':'R10-public-evidence-readback'})
  with urllib.request.urlopen(req,timeout=90) as r:
   code=r.status
   for data in iter(lambda:r.read(4*1024*1024),b''):h.update(data);n+=len(data)
  return dict(name=a['name'],http=code,bytes=n,expected_bytes=a['bytes'],blake2b_256=h.hexdigest(),match=n==a['bytes'] and h.hexdigest()==a['blake2b_256'],seconds=time.perf_counter()-t)
 except Exception as e:return dict(name=a['name'],bytes=n,match=False,error=str(e),seconds=time.perf_counter()-t)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 rows=[]
 for r in pool.map(asset,manifest['assets']):
  rows.append(r);print(json.dumps(r),flush=True)
  (P/'audit/PUBLIC_ASSETS_PARTIAL.json').write_text(json.dumps(rows,indent=2)+'\n')
base='https://migodam.github.io/b1-probabilistic-maxwell-imaging/'
pages=[('research/p3l_r10_bilateral_structure/','R10-BI'),('research/maxwell_formula_atlas/','140'),('research/p3l_r10_bilateral_structure/RESEARCH_REPORT_ZH.md','51.82%')]
checks=[]
for tail,token in pages:
 try:
  with urllib.request.urlopen(base+tail,timeout=30) as r:s=r.read().decode();checks.append(dict(url=base+tail,status=r.status,token=token,match=token in s))
 except Exception as e:checks.append(dict(url=base+tail,match=False,error=str(e)))
result=dict(status='PASS' if all(r['match'] for r in rows+checks) else 'INCOMPLETE',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authenticated=False,remote_bytes_read=sum(r['bytes'] for r in rows),local_extra_archive_copies=0,assets=rows,pages=checks,release_tag=manifest['release_tag'],scientific_commit=manifest['scientific_commit'],scope='Public deliverability and lossless evidence identity; not a new scientific or user-acceptance claim')
(P/'PUBLIC_READBACK.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status=result['status'],bytes=result['remote_bytes_read'])),flush=True)
