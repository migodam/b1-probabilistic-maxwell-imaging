"""Global50 unique-host public seal barrier; no GT or new forward."""
import os,time,json,argparse,shutil
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
p=argparse.ArgumentParser();p.add_argument('--wait',action='store_true');args=p.parse_args()
OUT=ROOT/'delegated/vector/confirmation_host_union';OUT.mkdir(parents=True,exist_ok=True)
def dump(path,x):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(path)
queue=json.loads((ROOT/'configs/vector_confirmation_queue.json').read_text());dispatch=json.loads((ROOT/'configs/vector_confirmation_host_dispatch.json').read_text());mapping={t['task_id']:t['host'] for t in dispatch['tasks']}
if len(mapping)!=50 or set(mapping)!={t['task_id'] for t in queue['tasks']}:raise ValueError('Host partition is not exactly the frozen50 unique tasks')
while True:
 dispatch=json.loads((ROOT/'configs/vector_confirmation_host_dispatch.json').read_text());mapping={t['task_id']:t['host'] for t in dispatch['tasks']}
 ready=[];pending=[]
 for task in queue['tasks']:
  dest=ROOT/'delegated/vector/inversions'/task['task_id']
  if all((dest/s).exists() for s in ['RECEIPT.json','SEALED.json','prediction.npz','COST_PREFIX_SEAL.json','COST_PREFIX_MATERIALS.npz']):ready.append(task)
  else:pending.append(task['task_id'])
 counts={h:sum(mapping[t['task_id']]==h for t in ready) for h in ['MAC_CPU','WINDOWS_CUDA']}
 dump(OUT/'COLLECTION_STATUS.json',dict(status='ALL50_READY_FOR_PUBLIC_VALIDATION' if not pending else 'RUNNING',global_expected=50,expected_by_host=dispatch['host_counts'],completed_local_artifacts=len(ready),completed_by_host=counts,pending=pending,GT_access=False,original_requested_remote50_actual_initial_dispatch25=True))
 if not pending:break
 if not args.wait:raise SystemExit(0)
 time.sleep(10)
import numpy as np
seals=[];prefixes=[];cost_rows=[];groups={};seal_dir=OUT/'public_seals';seal_dir.mkdir(exist_ok=False)
for task in ready:
 dest=ROOT/'delegated/vector/inversions'/task['task_id'];r=json.loads((dest/'RECEIPT.json').read_text());seal=json.loads((dest/'SEALED.json').read_text());prefix=json.loads((dest/'COST_PREFIX_SEAL.json').read_text())
 if r['status']!='SEALED_PUBLIC_PREDICTION_NO_SCORE' or seal['status']!='PREDICTION_SEALED_BEFORE_ROOT_OFFLINE_SCORING' or prefix['status']!='PUBLIC_COST_PHI_SELECTION_SEALED_WITHOUT_GT':raise ValueError('Required scientific public seal status missing')
 if r.get('host_dispatch_role')!=mapping[task['task_id']]:raise ValueError('Host result does not match frozen dispatch')
 with np.load(dest/'prediction.npz',allow_pickle=False) as z:
  shapes={k:list(z[k].shape) for k in ['material','prediction','J','Q']}
  if shapes!={'material':[1728],'prediction':[180,60],'J':[5184,60],'Q':[5184,60]}:raise ValueError('Not the declared realN/true3N all60 model')
  if not np.isfinite(z['material']).all() or z['material'].min()<0 or z['material'].max()>8:raise ValueError('Illegal physical material')
 for name in ['SEALED.json','COST_PREFIX_SEAL.json']:shutil.copy2(dest/name,seal_dir/(task['task_id']+'_'+name))
 row=dict(task=task,host=mapping[task['task_id']],public_seal=seal,prefix_seal=prefix,array_shapes=shapes,truth_access=False)
 seals.append(row);prefixes.extend([dict(task_id=task['task_id'],public_id=task['public_id'],host=mapping[task['task_id']],cost_scope=prefix.get('cost_scope','NORMAL_FIRST_ATTEMPT_END_TO_END_ARM'),end_to_end_equal_cost_qualified=not bool(r.get('infrastructure_attempt_metadata')),**x) for x in prefix['prefixes']])
 cost_rows.append(dict(task_id=task['task_id'],public_id=task['public_id'],host=mapping[task['task_id']],nested_CPU_seconds=r['CPU_seconds'],nested_wall_seconds=r['wall_seconds'],infrastructure_attempt_metadata=r.get('infrastructure_attempt_metadata'),cost_scope='RETRY_ATTEMPT_ONLY_NOT_END_TO_END' if r.get('infrastructure_attempt_metadata') else 'NORMAL_FIRST_ATTEMPT_ARM',physical_LU=r['counts']['physical_LU'],original_current=r['current_metric'],material_readout=r['material_metric'],linalg_library=r['torch_linalg_library']))
 groups.setdefault(task['public_id'],[]).append(task['task_id'])
if len(seals)!=50 or any(len(v)!=5 for v in groups.values()) or len(groups)!=10:raise ValueError('Every frozen observation group must have its5 methods')
dump(OUT/'ALL50_PUBLIC_PREDICTION_AND_PREFIX_SEALS.json',dict(status='ALL50_UNIQUE_HOST_PREDICTIONS_AND_PUBLIC_COST_PREFIXES_SEALED_BEFORE_FRESH_GT',global_tasks=50,hosts=dispatch['host_counts'],matched_tasks=40,crossgrid_tasks=10,independent_materials=8,cross_new_independent_materials=0,seals=seals,fresh_GT_access=False,fresh_GT_scoring_authorized=False,permission_scope='Root must separately authorize fresh GT scoring after this barrier'))
dump(OUT/'ALL_PUBLIC_COST_PREFIX_CHOICES.json',prefixes);dump(OUT/'NESTED_ARM_COST_BY_HOST.json',dict(rows=cost_rows,scope='Nested arm costs only; platform seconds kept separate, whole process wait4/GetProcessTimes replace nested scopes, concurrent wall not additive'))
dump(OUT/'MERGE_COST.json',dict(status='COMPLETE',CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,new_Maxwell_RHS=0,GT_access=False,scope='Own finite local collection monitor/public array shape+box validation/copies, wall overlaps scientific work'))
print(json.dumps(dict(status='ALL50_PUBLICLY_SEALED_NO_FRESH_GT_READ',hosts=dispatch['host_counts'],output=str(OUT))),flush=True)
