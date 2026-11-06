"""Transfer and run the finite frozen CPU queue on the authorized existing i9.

No service, drivers, GPU allocation or external job control. SSH credentials
are CLI arguments, never serialized. Result archives contain public predictions.
"""
from pathlib import Path
import argparse,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['launch','run','status','fetch']);ap.add_argument('--host',required=True);ap.add_argument('--identity',required=True);ap.add_argument('--remote-root',default='D:/AI/MAX_STRUCT_R8_20261106');ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--queue-name',default='warm_xinan');ap.add_argument('--manifest',default='research/research_max_struct_r8/configs/warm_xinan_manifest.json');ap.add_argument('--wait-for',default='');a=ap.parse_args()
    base=['ssh','-i',a.identity,'-o','BatchMode=yes','-o','ConnectTimeout=10',a.host,'D:/python/python.exe','-B','-'];start=time.perf_counter()
    def remote(script):
        p=subprocess.run(base,input=script,text=True,capture_output=True,timeout=900)
        if p.returncode:raise RuntimeError(p.stderr+'\n'+p.stdout)
        return json.loads(p.stdout.lstrip('\ufeff'))
    root=repr(a.remote_root)
    if a.action=='run':
        # Windows SSH owns a job object; keep this finite SSH session connected.
        # The previous detached launcher was observed to exit before any task.
        code="""import os,sys,runpy,psutil,json,time
from pathlib import Path
r=Path(REMOTE_ROOT);os.chdir(r);psutil.Process().cpu_affinity(list(range(0,16,2)))
os.environ['PYTHONPATH']=str(r/'.r8_dependencies')+os.pathsep+str(r);os.environ['PYTHONUTF8']='1'
if WAIT_FOR:
 print(json.dumps(dict(status='WAITING_FOR_OWN_PRIMARY_QUEUE',path=WAIT_FOR)),flush=True)
 while True:
  try: previous=json.loads(Path(WAIT_FOR).read_text())
  except json.JSONDecodeError: time.sleep(2);continue
  if previous['status']=='COMPLETED_UNSCORED':break
  if previous['status']=='EXECUTION_EXCEPTION_REQUIRES_REPAIR':raise RuntimeError('Primary queue needs repair')
  time.sleep(10)
sys.argv=['research/research_max_struct_r8/experiments/run_queue.py','--manifest',MANIFEST,'--name',QUEUE_NAME,'--workers',WORKERS]
runpy.run_path(sys.argv[0],run_name='__main__')
""".replace('REMOTE_ROOT',root).replace('WORKERS',repr(str(a.workers))).replace('WAIT_FOR',repr(a.wait_for)).replace('MANIFEST',repr(a.manifest)).replace('QUEUE_NAME',repr(a.queue_name))
        with (ROOT/'audit'/(a.queue_name+'_connected.log')).open('w') as log:
            p=subprocess.run(base,input=code,text=True,stdout=log,stderr=subprocess.STDOUT)
        result=dict(returncode=p.returncode,connected_SSH=True,wall_seconds=time.perf_counter()-start)
        (ROOT/'audit'/('XINAN_CONNECTED_RUN.json' if a.queue_name=='warm_xinan' else a.queue_name+'_CONNECTED_RUN.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result))
        return
    if a.action=='launch':
        prep=remote("from pathlib import Path\nimport json\nr=Path("+root+");r.mkdir(parents=True,exist_ok=False)\nprint(json.dumps(dict(created=str(r))))")
        archive=ROOT/'audit/transfer_warm_xinan.zip';remote_zip=a.remote_root+'/bundle.zip'
        p=subprocess.run(['scp','-i',a.identity,'-o','BatchMode=yes',str(archive),a.host+':'+remote_zip],capture_output=True,text=True,timeout=180)
        if p.returncode:raise RuntimeError(p.stderr)
        code="""from pathlib import Path
import zipfile,json,subprocess,sys,os,psutil
r=Path(REMOTE_ROOT)
with zipfile.ZipFile(r/'bundle.zip') as z:z.extractall(r)
audit=r/'research/research_max_struct_r8/audit';audit.mkdir(parents=True,exist_ok=True)
os.chdir(r);psutil.Process().cpu_affinity(list(range(0,16,2)))
with (audit/'warm_xinan.log').open('w') as log:
 p=subprocess.Popen([sys.executable,'research/research_max_struct_r8/experiments/run_queue.py','--manifest','research/research_max_struct_r8/configs/warm_xinan_manifest.json','--name','warm_xinan','--workers',WORKERS],cwd=r,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NEW_PROCESS_GROUP|subprocess.DETACHED_PROCESS)
rec=dict(status='FINITE_QUEUE_LAUNCHED',pid=p.pid,affinity=list(range(0,16,2)),workers=int(WORKERS),GPU=False,root=str(r),truth_access='NONE')
(audit/'XINAN_LAUNCH.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec))
""".replace('REMOTE_ROOT',root).replace('WORKERS',repr(str(a.workers)))
        result=remote(code);result.update(transfer_and_launch_wall_seconds=time.perf_counter()-start,archive_bytes=archive.stat().st_size)
        (ROOT/'audit/XINAN_LAUNCH.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
    else:
        fetched_path=ROOT/'audit'/(a.queue_name+'_FETCHED_TASKS.json')
        fetched=json.loads(fetched_path.read_text())['task_ids'] if fetched_path.exists() else []
        code="""from pathlib import Path
import json,psutil,zipfile,time,shutil
r=Path(REMOTE_ROOT);root=r/'research/research_max_struct_r8';p=root/'audit'/('queue_'+QUEUE_NAME)/'RECEIPT.json'
d=json.loads(p.read_text());out=dict(status=d['status'],completed=len(d['tasks']),failures=d['failed_attempts'],last=d['tasks'][-3:],elapsed_wall_seconds=d.get('elapsed_wall_seconds'))
if FETCH:
 cpu,wall=time.process_time(),time.perf_counter();copied=[];skip=set(SKIP_IDS)
 path=r/'public_results.zip'
 with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
  for t in d['tasks']:
   if t['task_id'] not in skip and t['status'] in ('SEALED_PUBLIC_PREDICTION','FAILED'):
    for f in (r/t['output_dir']).rglob('*'):
     if f.is_file():z.write(f,f.relative_to(r))
    copied.append(t['task_id'])
  for f in (root/'audit').rglob('*'):
   if f.is_file():z.write(f,f.relative_to(r))
 out['archive_bytes']=path.stat().st_size
 out.update(transferred_tasks=copied,remote_pack_cpu_seconds=time.process_time()-cpu,remote_pack_wall_seconds=time.perf_counter()-wall,remote_free_disk_bytes=shutil.disk_usage(r).free)
print(json.dumps(out))
""".replace('REMOTE_ROOT',root).replace('FETCH',repr(a.action=='fetch')).replace('QUEUE_NAME',repr(a.queue_name)).replace('SKIP_IDS',repr(fetched))
        result=remote(code)
        if a.action=='fetch':
            dest=ROOT/'audit'/(a.queue_name+'_public_results.zip')
            p=subprocess.run(['scp','-i',a.identity,'-o','BatchMode=yes',a.host+':'+a.remote_root+'/public_results.zip',str(dest)],capture_output=True,text=True,timeout=1200)
            if p.returncode:raise RuntimeError(p.stderr)
            import zipfile,shutil
            from datetime import datetime,timezone
            # Remote bundles carry an earlier copy of the controller audit.
            # Import it as provenance, never overwrite a local public seal or
            # its already-frozen evidence merely to refresh a queue snapshot.
            audit_prefix=Path('research/research_max_struct_r8/audit')
            mirror=ROOT/'audit/remote_imports'/a.queue_name/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
            with zipfile.ZipFile(dest) as z:
                for member in z.infolist():
                    if member.is_dir():continue
                    rel=Path(member.filename)
                    target=(mirror/rel.relative_to(audit_prefix)) if rel.is_relative_to(audit_prefix) else REPO/rel
                    target.parent.mkdir(parents=True,exist_ok=True)
                    with z.open(member) as source,target.open('wb') as output:shutil.copyfileobj(source,output)
            result['remote_audit_mirror']=str(mirror.relative_to(REPO))
            fetched_path.write_text(json.dumps(dict(task_ids=sorted(set(fetched)|set(result['transferred_tasks']))),indent=2))
        result['control_wall_seconds']=time.perf_counter()-start
        if a.action=='fetch':
            with (ROOT/'audit'/(a.queue_name+'_FETCH_HISTORY.jsonl')).open('a') as f:f.write(json.dumps(result)+'\n')
        (ROOT/'audit'/((('XINAN_' if a.queue_name=='warm_xinan' else a.queue_name+'_')+a.action.upper())+'.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result))

if __name__=='__main__':main()
