"""Finite publication transport for already verified evidence; no scientific execution."""
import json,subprocess,time,resource,concurrent.futures
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'release_assets/transport_parts';D.mkdir(exist_ok=True)
repo='migodam/b1-probabilistic-maxwell-imaging';tag='research-max-struct-r12-20261111';chunk=128*1024*1024
start=time.monotonic();r0=resource.getrusage(resource.RUSAGE_CHILDREN);c0=time.process_time()
assets=json.loads((ROOT/'UPLOAD_MANIFEST_READY.json').read_text())['assets']
existing=json.loads(subprocess.check_output(['gh','release','view',tag,'--repo',repo,'--json','assets'],text=True))['assets'];uploaded={x['name']:x['size'] for x in existing if x['state']=='uploaded'}
manifest={'status':'TRANSPORT_PREPARED_NOT_PUBLICATION_ACCEPTED','part_size_bytes':chunk,'order':'Parts concatenated in listed numeric order reproduce the exact verified logical file; then zstd decode if suffix .zst, then unzip.','logical_assets':[],'hashes':'NOT_COMPUTED'}
tasks=[]
for a in assets:
 if not a.get('ready'):continue
 src=Path(a['upload_local_path']);assert src.stat().st_size==a['bytes'];parts=[]
 with src.open('rb') as f:
  i=0
  remaining=a['bytes']
  while remaining:
   size=min(chunk,remaining);i+=1;name=f'{a["name"]}.part{i:03d}';target=D/name
   if target.exists() and target.stat().st_size==size:
    f.seek(size,1)  # Already prepared from the unchanged, verified logical asset.
   else:
    b=f.read(size);assert len(b)==size;target.write_bytes(b)
   parts.append({'name':name,'bytes':size,'url':f'https://github.com/{repo}/releases/download/{tag}/{name}'})
   if uploaded.get(name)!=size:tasks.append(target)
   remaining-=size
 manifest['logical_assets'].append({'name':a['name'],'bytes':a['bytes'],'content_index':a['content_index'],'parts':parts})
(ROOT/'PUBLICATION_PARTS_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(f'Prepared {len(manifest["logical_assets"])} logical assets; {len(tasks)} missing parts to upload',flush=True)
def upload(p):
 t=time.monotonic();r=subprocess.run(['gh','release','upload',tag,str(p),'--repo',repo],capture_output=True,text=True)
 return {'name':p.name,'bytes':p.stat().st_size,'returncode':r.returncode,'wall_seconds':time.monotonic()-t,'stderr':r.stderr[-1500:]}
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for r in pool.map(upload,tasks):
  rows.append(r)
  with (ROOT/'audit/release/PART_UPLOAD_ATTEMPTS.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
  print(json.dumps(r),flush=True)
r1=resource.getrusage(resource.RUSAGE_CHILDREN)
with (ROOT/'audit/release/PART_UPLOAD_COSTS.jsonl').open('a') as f:f.write(json.dumps({'parts_attempted':len(tasks),'successful':sum(x['returncode']==0 for x in rows),'wall_seconds':time.monotonic()-start,'controller_CPU_seconds':time.process_time()-c0,'children_CPU_seconds':r1.ru_utime+r1.ru_stime-r0.ru_utime-r0.ru_stime,'source_deleted':False})+'\n')
raise SystemExit(any(x['returncode'] for x in rows))
