"""Finite sequential DEV/readout then50-confirmation controller; no science/GT."""
import os,time,json,sys,subprocess,ctypes,importlib.util
from pathlib import Path
from ctypes import wintypes as w
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.vector_resources import configure_own_windows_process
OUT=ROOT/'delegated/vector/remaining_controller';OUT.mkdir(parents=True,exist_ok=False)
process_configuration=configure_own_windows_process()
spec=importlib.util.spec_from_file_location('own_job',REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);job=module.job_object();os.environ['R10_PARENT_JOB_ACTIVE']='1'
def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
def child_cpu(child):
 a,b,c,d=(w.FILETIME() for _ in range(4));k=ctypes.WinDLL('kernel32',use_last_error=True)
 k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
 if not k.GetProcessTimes(w.HANDLE(int(child._handle)),ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
 return ((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7
phases=[];status='FAILED';failure=None
try:
 wait=time.perf_counter();prior=ROOT/'delegated/vector/dev_ordinary_process/PROCESS_COST.json'
 print(json.dumps(dict(event='WAITING_FOR_EXISTING_ORDINARY_DEV_TERMINAL',no_CUDA_context=True)),flush=True)
 while True:
  if prior.exists():
   p=json.loads(prior.read_text())
   if p['status']!='COMPLETE':raise RuntimeError('Existing ordinary DEV failed; preserve outputs, no later launch')
   break
  time.sleep(10)
 dump(OUT/'PREVIOUS_PHASE_WAIT.json',dict(wall_seconds=time.perf_counter()-wait,CPU_seconds=time.process_time(),no_CUDA_context=True,scope='Own controller wait for already-running authorized ordinary DEV; not scientific GPU time'))
 for name,config,queue,count in [('dev_denominator','configs/vector_dev_contract.json','configs/vector_dev_denominator_queue.json',6),('confirmation','configs/vector_confirmation_execution.json','configs/vector_confirmation_queue.json',50)]:
  out='delegated/vector/'+name+'_process';command=[sys.executable,str(ROOT/'experiments/vector_run.py'),'--config',config,'--queue',queue,'--max-tasks',str(count),'--output',out]
  start=time.perf_counter();known=set();log=OUT/(name+'.log')
  with log.open('x') as handle:
   child=subprocess.Popen(command,cwd=REPO,stdout=handle,stderr=subprocess.STDOUT)
   dump(OUT/'ACTIVE.json',dict(phase=name,pid=child.pid,command=command,started_unix=time.time(),requested_tasks=count))
   print(json.dumps(dict(event='PHASE_STARTED',phase=name,pid=child.pid,tasks=count)),flush=True)
   while child.poll() is None:
    terminal=ROOT/out/'TERMINALS.json'
    if terminal.exists():
     info=json.loads(terminal.read_text())
     for t in info['terminals']:
      if t['task_id'] not in known:
       known.add(t['task_id']);receipt=json.loads((ROOT/'delegated/vector/inversions'/t['task_id']/'RECEIPT.json').read_text())
       print(json.dumps(dict(event='PHASE_TASK_SEALED',phase=name,completed=len(known),requested=count,task_id=t['task_id'],Phi=receipt['physical_Phi'],CPU_seconds=receipt['CPU_seconds'],physical_LU=receipt['counts'].get('physical_LU'),truth_scoring='NOT_RUN')),flush=True)
    time.sleep(5)
   code=child.returncode;cpu=child_cpu(child)
  row=dict(phase=name,exit_code=code,CPU_seconds=cpu,wall_seconds=time.perf_counter()-start,requested_tasks=count,CPU_source='Windows GetProcessTimes whole scientific child, replace nested CPU',additive=True)
  phases.append(row);dump(OUT/'PHASE_COSTS.json',dict(phases=phases));print(json.dumps(dict(event='PHASE_TERMINAL',**row)),flush=True)
  if code:raise RuntimeError(name+' failed, no retry/extra budget')
 status='COMPLETE'
except BaseException as exc:
 import traceback
 failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
finally:
 dump(OUT/'PROCESS_COST.json',dict(status=status,failure=failure,CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,phases=phases,
   process=process_configuration,job=job,scope='Own finite controller CPU only, scientific child OS CPU in PHASE_COSTS and helper CPU in child receipts; wall overlaps child work and is not additive',truth_access=False,additive=True))
raise SystemExit(0 if status=='COMPLETE' else 1)
