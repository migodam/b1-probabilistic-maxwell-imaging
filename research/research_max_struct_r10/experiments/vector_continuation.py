"""Posthoc public-endpoint physical-budget diagnostic, never reads GT."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,traceback
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.vector_run import dump,stamp,initialize,configure,save_physical

def main():
 out=ROOT/'RESULTS/vector_continuation_v04';out.mkdir(exist_ok=False);model=None;failure=None;status='FAILED'
 try:
  config,job=configure();dump(out/'PROCESS_CONFIGURATION.json',config)
  import torch,numpy as np,psutil
  process=psutil.Process();torch.set_num_threads(1);torch.set_num_interop_threads(1)
  from research.research_max_struct_r10.src.vector_physical import PhysicalObjective,physical_corrector
  contract=json.loads((ROOT/'configs/vector_continuation_v04.json').read_text());physics=json.loads((ROOT/'configs/vector_actual_contract.json').read_text());physics['execution_backend']='cuda-dense'
  dump(out/'FROZEN_CONFIG.json',contract);dump(out/'PHYSICS_CONFIG.json',physics)
  dump(out/'SOURCE_EVIDENCE.json',[stamp(Path(__file__)),stamp(ROOT/'src/vector_physical.py'),stamp(ROOT/'src/vector_cuda_dense.py'),stamp(ROOT/'configs/vector_continuation_v04.json')])
  model=initialize(physics,out,12);public=ROOT/'data_public/vector/v04_matched12'
  with np.load(public/'observations.npz',allow_pickle=False) as z:
   assert z.files==['Y'];Y=z['Y'].copy()
  terminals=[]
  for arm in contract['arms']:
   folder=out/arm;folder.mkdir();model.dense_state=None;objective=PhysicalObjective(model,Y,upper=8.)
   inp=ROOT/contract['initial_materials'][arm]
   with np.load(inp,allow_pickle=False) as z:m=z['material'].copy()
   dump(folder/'INPUT_EVIDENCE.json',[stamp(inp),stamp(public/'observations.npz'),stamp(public/'PUBLIC.json')])
   start=time.perf_counter();c=time.process_time();before=dict(model.counts);event=torch.cuda.Event(enable_timing=True);end=torch.cuda.Event(enable_timing=True);event.record()
   def emit(row):
    with (folder/'TRACE.jsonl').open('a') as h:h.write(json.dumps(row)+'\n')
    print(json.dumps(dict(arm=arm,**row)),flush=True)
   best,corrector=physical_corrector(objective,m,contract['corrector'],arm,lambda label,record:save_physical(folder,label,record),emit)
   save_physical(folder,'prediction',best);end.record();end.synchronize()
   receipt=dict(status='SEALED_POSTHOC_PUBLIC_BUDGET_DIAGNOSTIC',arm=arm,physical_Phi=best['Phi'],PG_inf=best['projected_gradient_inf'],physical_corrector=corrector,counts={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)},CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-start,CUDA_event_stream_interval_ms=event.elapsed_time(end),CUDA_event_scope='Stream timestamp interval including host/gate/IO gaps, not active device time',additive=False,included_in='../PROCESS_COST.json',GT_read=False,prefix_cost='Original18 public prefix reference only; no original prefix CPU charged again',initial_accepted_index=80,reset_LBFGS_history=True)
   dump(folder/'RECEIPT.json',receipt);dump(folder/'SEALED.json',dict(status='SEALED_BEFORE_ROOT_OFFLINE_DIAGNOSTIC_SCORE',files=[stamp(p) for p in sorted(folder.iterdir()) if p.is_file()],GT_read=False))
   terminals.append(dict(arm=arm,status=receipt['status']));dump(out/'TERMINALS.json',dict(terminals=terminals,all_terminal=len(terminals)==3,scoring='ROOT_ONLY_NOT_RUN'))
  status='COMPLETE'
 except BaseException as exc:failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
 finally:
  import psutil
  p=psutil.Process();t=p.cpu_times()
  dump(out/'PROCESS_COST.json',dict(status=status,platform='XINAN',CPU_seconds=t.user+t.system,wall_seconds=time.perf_counter()-ENTRY,peak_RSS_bytes=p.memory_info().peak_wset,additive=True,scope='Whole posthoc3-arm process including coldsetup/alltrials/save, prefix excluded',counts=dict(model.counts) if model else None,GPU_gate_helper_CPU_seconds=model.gate.helper_CPU_seconds if model else None,GPU_gate_helper_scope='Separate query children not included in main CPU',failure=failure,GT_read=False))
 return 0 if status=='COMPLETE' else 1
if __name__=='__main__':raise SystemExit(main())
