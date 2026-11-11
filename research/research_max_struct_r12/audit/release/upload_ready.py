import concurrent.futures,json,subprocess,time,resource
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
manifest=json.loads((ROOT/'UPLOAD_MANIFEST_READY.json').read_text())
repo='migodam/b1-probabilistic-maxwell-imaging';tag='research-max-struct-r12-20261111'
existing=json.loads(subprocess.check_output(['gh','release','view',tag,'--repo',repo,'--json','assets'],text=True))['assets']
sizes={x['name']:x['size'] for x in existing}
tasks=[a for a in manifest['assets'] if a.get('ready') and sizes.get(a['name'])!=a['bytes']]
start=time.monotonic();cpu=resource.getrusage(resource.RUSAGE_CHILDREN);rows=[]
def run(a):
 p=Path(a['upload_local_path']);assert p.stat().st_size==a['bytes']
 t=time.monotonic();r=subprocess.run(['gh','release','upload',tag,str(p),'--repo',repo],capture_output=True,text=True)
 return {'name':a['name'],'bytes':a['bytes'],'returncode':r.returncode,'wall_seconds':time.monotonic()-t,'stderr':r.stderr[-1800:],'stdout':r.stdout[-500:]}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for r in pool.map(run,tasks):
  rows.append(r)
  with (ROOT/'audit/release/UPLOAD_ATTEMPTS.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
  print(json.dumps(r),flush=True)
c=resource.getrusage(resource.RUSAGE_CHILDREN)
with (ROOT/'audit/release/UPLOAD_BATCH_COSTS.jsonl').open('a') as f:f.write(json.dumps({'assets_attempted':len(tasks),'wall_seconds':time.monotonic()-start,'children_CPU_seconds':c.ru_utime+c.ru_stime-cpu.ru_utime-cpu.ru_stime,'successful':sum(r['returncode']==0 for r in rows)})+'\n')
raise SystemExit(any(r['returncode'] for r in rows))
