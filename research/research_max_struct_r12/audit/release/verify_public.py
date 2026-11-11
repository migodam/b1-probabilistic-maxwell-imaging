import json,time,urllib.request,concurrent.futures
from pathlib import Path
m=json.loads(Path('PUBLICATION_PARTS_MANIFEST.json').read_text())
def check(p):
 req=urllib.request.Request(p['url'],headers={'Range':'bytes=0-15','User-Agent':'R12-public-readback'})
 with urllib.request.urlopen(req,timeout=45) as r:
  b=r.read(16);h=r.headers.get('Content-Range','');ok=r.status==206 and h.endswith('/'+str(p['bytes'])) and len(b)==16
  return {'name':p['name'],'url':p['url'],'status':r.status,'content_range':h,'sample_bytes':len(b),'ok':ok}
parts=[p for a in m['logical_assets'] for p in a['parts']]
t=time.monotonic()
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:r=list(ex.map(check,parts))
out={'status':'PUBLIC_PARTS_VERIFIED' if all(x['ok'] for x in r) else 'FAILED','parts':r,'logical_assets':12,'logical_bytes':sum(a['bytes'] for a in m['logical_assets']),'wall_seconds':time.monotonic()-t,'verification':'Anonymous ranged GET, each complete remote size checked; only 16 bytes read per part. No full repeat download.'}
Path('audit/PUBLICATION_READBACK.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],len(r),flush=True)
