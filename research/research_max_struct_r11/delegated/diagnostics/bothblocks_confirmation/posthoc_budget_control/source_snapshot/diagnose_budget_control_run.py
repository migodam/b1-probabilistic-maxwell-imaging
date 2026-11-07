"""One frozen POSTHOC_BUDGET_CONTROL corrector restart, originalR6 numerical body."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,argparse,traceback,shutil
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control';sys.path.insert(0,str(REPO))
from research.research_max_struct_r11.src.twod_budget import WholeChildBudget,CPUBudgetReached

def dump(p,x):
 def conv(v):
  if hasattr(v,'tolist'):return v.tolist()
  raise TypeError(type(v).__name__)
 p.write_text(json.dumps(x,indent=2,allow_nan=False,default=conv)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pair-id',required=True);a=ap.parse_args();manifest=json.loads((OUT/'FROZEN_MANIFEST.json').read_text());row=next(x for x in manifest['pairs'] if x['pair_id']==a.pair_id)
 if not row['trigger']:raise ValueError('Untriggered pairs must not execute or idle-fill')
 d=Path(row['output_dir']);d.mkdir(parents=True,exist_ok=False);budget=WholeChildBudget(row['remaining_new_child_CPU_allowance_seconds']);model=None;observer=None;corrector=None;failure=None;status='FAILED';started=time.perf_counter()
 try:
  budget.check()
  import numpy as np
  from threadpoolctl import threadpool_limits,threadpool_info
  from research.research_fbg_r3_experiments.src.physics import MaxwellModel
  from research.research_max_struct_r6.src.benchmark_solver import PaidStateObserver,refine_observed
  budget.check()
  with threadpool_limits(limits=1):
   public={}
   for p in row['public_files']:
    with np.load(REPO/p) as z:public.update(dict(z))
   model=MaxwellModel(public);Y=public['Y']
   with np.load(row['original_prediction']) as z:m=z['m'].copy()
   def snapshot(name,state):
    arrays={k:state[k] for k in ('m','j','Q','Y_prediction','gradient') if k in state};np.savez_compressed(d/(name+'.npz'),**arrays);dump(d/(name+'.json'),{k:v for k,v in state.items() if k not in arrays})
   observer=budget.observer_class(PaidStateObserver)(model,Y,snapshot,gradient_checkpoints=(100,500,1000));budget.check()
   _,corrector=refine_observed(observer,Y,m,maximum=1_000_000,accepted_checkpoints=(100,500,1000))
   status='CORRECTOR_RETURNED_UNDER_REMAINING_CPU'
 except CPUBudgetReached:status='REMAINING_CPU_BUDGET_REACHED'
 except BaseException as e:status='FAILED_CONTINUATION';failure=dict(type=type(e).__name__,error=str(e),traceback=traceback.format_exc())
 finally:
  eligible=observer.eligible_best if observer is not None else None;newbetter=bool(eligible is not None and eligible['objective']<row['original_Phi']);selected_phi=float(eligible['objective']) if newbetter else row['original_Phi'];selected_kind='NEW_PAID_ELIGIBLE_LOWER_PHI' if newbetter else 'ORIGINAL_RETAINED'
  selected=d/'SELECTED_PREDICTION.npz'
  if newbetter:
   arrays={k:eligible[k] for k in ('m','j','Q','Y_prediction','gradient') if k in eligible};np.savez_compressed(selected,**arrays)
   chosen_cpu=next((v['CPU_seconds'] for v in observer.paid_states if v['eligible_at_cap'] and v['original_Phi']==selected_phi),None)
  else:shutil.copy2(row['original_prediction'],selected);chosen_cpu=0.
  if eligible is not None:
   np.savez_compressed(d/'NEW_ELIGIBLE_BEST.npz',**{k:eligible[k] for k in ('m','j','Q','Y_prediction','gradient') if k in eligible})
  if observer is not None:
   dump(d/'PAID_STATES.json',observer.paid_states)
   if observer.incumbents:np.savez_compressed(d/'NEW_INCUMBENTS.npz',CPU_seconds=np.asarray([p[0] for p in observer.incumbents]),m=np.asarray([p[1] for p in observer.incumbents]),original_Phi=np.asarray([p[2] for p in observer.incumbents]))
  if corrector is not None:dump(d/'CORRECTOR.json',corrector)
  receipt=dict(status=status,classification='POSTHOC_BUDGET_CONTROL',pair_id=row['pair_id'],case_id=row['case_id'],method=row['method'],truth_access='NONE',scoring='NOT_RUN',target_G19_actual_CPU_seconds=row['target_G19_actual_CPU_seconds'],original_actual_child_CPU_seconds=row['original_actual_child_CPU_seconds'],remaining_new_child_allowance_seconds=budget.seconds,
   selected_kind=selected_kind,selected_original_Phi=selected_phi,original_Phi=row['original_Phi'],selected_new_paid_CPU_seconds=chosen_cpu,selected_cumulative_eligible_CPU_seconds=row['original_actual_child_CPU_seconds']+(chosen_cpu or 0),
   original_prediction=row['original_prediction'],selected_prediction=str(selected),new_physical_cost=dict(model.cost) if model is not None else {},new_paid_evaluations=len(observer.paid_states) if observer else 0,
   near_exit_new_child_CPU_seconds=time.process_time(),new_child_CPU_overrun_near_exit=max(0.,time.process_time()-budget.seconds),wall_seconds=time.perf_counter()-started,cap_stop_CPU_seconds=budget.stop_cpu,failure=failure,
   eligibility='Only paid states completed at absolute new-child CPU<=remaining; originals retained on exact Phi tie; late atomic states excluded even if lowerPhi',cost='Full external actual child CPU measured by queue includes imports/setup/serialization/exit and is charged even if original selected',history='LBFGS restart, no old history or LU reused; not seamless extension',numerical='Unchanged r6.refine_observed: originalR/.02..20/scale1e6/ftol1e-12/gtol.01/maxls30; maxiter1000000 nonbinding',threads=threadpool_info() if 'threadpool_info' in locals() else None)
  dump(d/'RECEIPT.json',receipt)
 print(json.dumps(dict(pair_id=row['pair_id'],status=status,selected=selected_kind)),flush=True);return 1 if status=='FAILED_CONTINUATION' else 0
if __name__=='__main__':raise SystemExit(main())
