"""Seal8×48 public actual-end-CPU-matched trajectory pairs, never read GT."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from twod_secondary_seal import all_terminal,read,dump,resolve,method
from twod_recovery import attempt_cost_info,affected_objects
ROOT=Path(__file__).resolve().parents[1]
def utc():return datetime.now(timezone.utc).isoformat()
def trace(task,receipt):
 d=resolve(task['output_dir']);p=d/'R11_INCUMBENT_TRAJECTORY.npz';records=read(d/'R11_PAID_PHYSICAL_STATES.json') if (d/'R11_PAID_PHYSICAL_STATES.json').exists() else []
 if not p.exists():return None,records
 with np.load(p,allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
 if len(a['m'])!=len(a['CPU_seconds']) or len(a['m'])!=len(a['original_Phi']) or (len(a['m']) and a['m'].shape[1:]!=(576,)):raise ValueError('Malformed public incumbent trajectory')
 if not np.isfinite(a['CPU_seconds']).all() or not np.isfinite(a['original_Phi']).all() or not np.isfinite(a['m']).all():raise ValueError('Nonfinite public trajectory')
 return a,records

def member(task,job,receipt,data,records,B,materials):
 T=job.get('child_process_CPU_seconds');cap=float(receipt.get('whole_child_CPU_budget_seconds',480));cost=receipt.get('physical_model_cost') or job.get('physical_model_cost') or {};physical=[r for r in records if r.get('CPU_seconds') is not None];eligible=[r for r in physical if r.get('eligible_at_cap')]
 before=[r for r in eligible if B is not None and r['CPU_seconds']<=B];last_before=max(before,key=lambda r:r['CPU_seconds']) if before else None;actual_last=max(physical,key=lambda r:r['CPU_seconds']) if physical else None
 row=dict(task_id=task['task_id'],method=method(task),source_queue_receipt=job.get('source_queue_receipt'),cpu_logical=job.get('cpu'),task_status=receipt.get('status',job['status']),actual_full_child_end_CPU_seconds=T,nominal_CPU_allowance=cap,nominal_unconsumed_CPU_allowance=max(0.,cap-T) if T is not None else None,actual_CPU_overrun=max(0.,T-cap) if T is not None else None,actual_end_physical_LU=cost.get('LU'),actual_end_complex_RHS=cost.get('actual_complex_RHS'),actual_last_recorded_paid_physical_CPU_seconds=actual_last['CPU_seconds'] if actual_last else None,last_cap_eligible_paid_physical_CPU_seconds_at_B=last_before['CPU_seconds'] if last_before else None,B_minus_last_cap_eligible_physical_record= B-last_before['CPU_seconds'] if B is not None and last_before else None,selected=None)
 row.update(attempt_cost_info(task['task_id'],T))
 indices=np.flatnonzero((data['CPU_seconds']<=B)&(data['CPU_seconds']<=cap)) if data is not None and B is not None else np.array([],int)
 if indices.size:
  index=int(indices[np.argmin(data['original_Phi'][indices])]);material_index=len(materials);materials.append(data['m'][index].copy());paid=float(data['CPU_seconds'][index])
  row['selected']=dict(material_index=material_index,trajectory_index=index,original_Phi=float(data['original_Phi'][index]),paid_CPU_seconds=paid,physical_LU=int(data['LU'][index]),physical_complex_RHS=int(data['complex_RHS'][index]),B_minus_selected_paid_CPU_seconds=B-paid)
 return row

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--rule',type=Path,default=ROOT/'configs/twod_independent48_actual_cpu_fairness.json');ap.add_argument('--output',type=Path,default=ROOT/'delegated/twod/independent48/actual_cpu_fairness/public_seal');a=ap.parse_args();wall=time.perf_counter();rule=read(a.rule)
 if rule['status']!='FROZEN_SECONDARY_PUBLIC_COST_MATCH_BEFORE48_GT' or rule['baselines']!=['cie_single','cie_dual'] or rule['comparators']!=['fullG19','freezeA','contractedT','global_dual']:raise ValueError('Frozen8pair rule required')
 by_case,jobs,receipts=all_terminal(a.manifest,a.queue);a.output.mkdir(parents=True,exist_ok=False);started=utc();rows=[];materials=[]
 dump(a.output/'PUBLIC_GATE_BEFORE_SELECTION.json',dict(status='PASS_ALL288_TERMINAL',created_utc=started,truth_access='NONE',new_physics_calls=0))
 for case,ts in sorted(by_case.items()):
  # Bound working memory to one object's six trajectories; no private data.
  loaded={m:trace(t,receipts[t['task_id']]) for m,t in ts.items()}
  for baseline in rule['baselines']:
   for comparator in rule['comparators']:
    lt,rt=ts[baseline],ts[comparator];lj,rj=jobs[lt['task_id']],jobs[rt['task_id']];lT,rT=lj.get('child_process_CPU_seconds'),rj.get('child_process_CPU_seconds');B=min(lT,rT) if lT is not None and rT is not None else None
    left=member(lt,lj,receipts[lt['task_id']],*loaded[baseline],B,materials);right=member(rt,rj,receipts[rt['task_id']],*loaded[comparator],B,materials)
    rows.append(dict(engineering_affected_object=case in affected_objects(),exact_end_to_end_cost_qualified=case not in affected_objects(),end_to_end_same_CPU_status='INDETERMINATE_ENGINEERING_ABORT_COST' if case in affected_objects() else 'QUALIFIED_RECORDED_ATTEMPT_COST_VIEW',case_id=case,pair_id=baseline+'__vs__'+comparator,baseline=baseline,comparator=comparator,common_actual_end_CPU_B=B,left=left,right=right,both_have_paid_incumbent=bool(left['selected'] is not None and right['selected'] is not None),selection='Minimum originalPhi withinB, no GT; no late-state backfill',actual_end_CPU_equal=bool(lT==rT) if lT is not None and rT is not None else None,exactly_equal_consumption_claimed=False,view_scope=rule['actual_spend_scope'],gap_scope=rule['state_discreteness']))
 if len(rows)!=384:raise ValueError('Must seal exactly384 public object/pair choices')
 np.savez_compressed(a.output/'PAIRED_MATERIALS.npz',m=np.stack(materials) if materials else np.empty((0,576),dtype=np.float64))
 dump(a.output/'PAIRED_SELECTIONS.json',dict(status='ALL384_PUBLIC_ACTUAL_CPU_PAIRS_SEALED',rows=rows,truth_access='NONE'))
 dump(a.output/'FROZEN_RULE_USED.json',rule)
 receipt=dict(status='ALL384_PUBLIC_ACTUAL_CPU_PAIRS_SEALED_BEFORE48_GT',registered_utc=rule['registered_utc'],registration_context=rule['registration_context'],seal_started_utc=started,seal_finished_utc=utc(),objects=48,pair_types=8,pairs=384,selected_materials=len(materials),pairs_with_both_paid_incumbents=sum(r['both_have_paid_incumbent'] for r in rows),new_inverse_runs=0,new_forward_calls=0,new_adjoint_calls=0,truth_access='NONE',GT_selected_budget=False,GT_selected_state=False,primary_unchanged=True,exactly_equal_consumption_claimed=False,CPU_scope='Absolute whole public selection process_time includes imports and per-object public trace reads through receipt preparation',CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall)
 dump(a.output/'RECEIPT.json',receipt);print(json.dumps({k:receipt[k] for k in ('status','pairs','pairs_with_both_paid_incumbents','CPU_seconds')}))
if __name__=='__main__':main()
