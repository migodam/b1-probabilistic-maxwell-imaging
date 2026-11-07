"""One finite early DEV harvest+old-exposed-GT score; fresh50 unaffected."""
import time,json,subprocess,resource,zipfile,sys
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/dev_early_postseal';OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def childcpu():
 r=resource.getrusage(resource.RUSAGE_CHILDREN);return r.ru_utime+r.ru_stime
steps=[];status='FAILED';failure=None
try:
 log=ROOT/'delegated/vector/remaining_transport/run.log';wait=time.perf_counter()
 while True:
  terminal=None
  if log.exists():
   for line in log.read_text().splitlines():
    try:r=json.loads(line)
    except json.JSONDecodeError:continue
    if r.get('event')=='PHASE_TERMINAL' and r.get('phase')=='dev_denominator':terminal=r
  if terminal:
   if terminal['exit_code']!=0:raise RuntimeError('New DEV phase failed; no GT score')
   break
  time.sleep(10)
 steps.append(dict(operation='wait_for_all_new_DEV6_terminal',wall_seconds=time.perf_counter()-wait,CPU_scope='Own controller only; no CUDA/GT'))
 tasks=json.loads((ROOT/'configs/vector_dev_denominator_queue.json').read_text())['tasks']
 bases=['research/research_max_struct_r11/delegated/vector/dev_denominator_process','research/research_max_struct_r11/delegated/vector/remaining_controller/PHASE_COSTS.json']+['research/research_max_struct_r11/delegated/vector/inversions/'+t['task_id'] for t in tasks]
 code="import zipfile,time;from pathlib import Path;r=Path('"+REMOTE+"');c=time.process_time();t=time.perf_counter();z=zipfile.ZipFile(r/'R11_DEV_EARLY_RESULT.zip','w');files=[p for base in "+repr(bases)+" for p in (([r/base] if (r/base).is_file() else (r/base).rglob('*'))) if p.is_file()];[z.write(p,str(p.relative_to(r)),compress_type=(zipfile.ZIP_STORED if p.suffix=='.npz' else zipfile.ZIP_DEFLATED)) for p in files];z.close();import json;print(json.dumps(dict(CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-t,files=len(files),scope='Whole remote Python earlyDEV archive,NPZ stored without recompression,no physics/GT')))"
 start=time.perf_counter();cpu0=childcpu();p=subprocess.run(SSH+['D:/python/python.exe -c "'+code+'"'],capture_output=True,text=True)
 (OUT/'REMOTE_PACK.log').write_text(p.stdout+p.stderr)
 if p.returncode:raise RuntimeError('Remote early DEV pack failed')
 for line in p.stdout.splitlines():
  try:x=json.loads(line)
  except json.JSONDecodeError:continue
  if 'CPU_seconds' in x:dump(OUT/'REMOTE_PACK_COST.json',x)
 steps.append(dict(operation='remote_pack',wall_seconds=time.perf_counter()-start,local_child_CPU_seconds=childcpu()-cpu0))
 start=time.perf_counter();cpu0=childcpu();p=subprocess.run(['scp','-i',KEY,HOST+':'+REMOTE+'/R11_DEV_EARLY_RESULT.zip',str(OUT/'RESULT.zip')],capture_output=True,text=True)
 (OUT/'DOWNLOAD.log').write_text(p.stdout+p.stderr)
 if p.returncode:raise RuntimeError('Early DEV download failed')
 with zipfile.ZipFile(OUT/'RESULT.zip') as z:z.extractall(REPO)
 steps.append(dict(operation='download_extract',wall_seconds=time.perf_counter()-start,local_child_CPU_seconds=childcpu()-cpu0))
 start=time.perf_counter();cpu0=childcpu()
 command=[sys.executable,str(ROOT/'experiments/vector_postseal_score.py'),'--cohort','dev','--authorization','OLD_EXPOSED_GT_ONLY_AFTER12SEALS','--authorization-note','Root explicitly authorized automatic old-exposed GT score after all12DEV public prediction seals; fresh50 GT remains forbidden']
 with (OUT/'SCORER.log').open('w') as f:p=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,cwd=REPO)
 steps.append(dict(operation='postseal_scorer',exit=p.returncode,whole_child_CPU_seconds=childcpu()-cpu0,wall_seconds=time.perf_counter()-start,scope='Replaces nested SCORING_COST CPU,not additive to it'))
 if p.returncode:raise RuntimeError('OldDEV scorer failed; evidence retained')
 status='COMPLETE'
except BaseException as exc:
 import traceback
 failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
finally:
 dump(OUT/'PROCESS_COST.json',dict(status=status,failure=failure,CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,steps=steps,
  scope='Own finite local controller CPU only; transport children/scorer separately listed; wall overlaps GPU science; oldDEV GT score only',fresh_GT_access=False,new_physics_RHS=0))
print(json.dumps(dict(status=status,failure=failure)),flush=True)
raise SystemExit(0 if status=='COMPLETE' else 1)
