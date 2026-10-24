#!/usr/bin/env python3
"""One post-method-freeze/pre-GT exploratory L1 FullJTR baseline; fixed64×2."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,argparse,traceback,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,Target,write,save,from_probit,BudgetStop
from p3r4.evaluation import checkpoint_identity
from p3r4.ordinary_trust_region import solve_control

JOB='ordinary_fullJ_TR_L1_64_addendum_v1';RESULT=OUT/'results/fresh/locked'/JOB;PROTOCOL=OUT/'configs/ORDINARY_FULLJ_TR_L1_64_ADDENDUM_PROTOCOL.json'
PARENTS=tuple(f'R4L1_{i:05d}' for i in range(64));ALG=ROOT/'src/p3r4/ordinary_trust_region.py';DEV_ASRUN=OUT/'II/jobs/ordinary_fullJ_TR_DEV16_v1/source_snapshot/ordinary_trust_region.py'
class ControlBudgetStop(RuntimeError):pass
def check_algorithm():
 if ALG.read_bytes()!=DEV_ASRUN.read_bytes():raise ValueError('mustreuseexactverifiedDEVasrunalgorithm,notnewtuning')
def protocol():
 check_algorithm();dev=json.loads((OUT/'configs/ORDINARY_FULLJ_TR_DEV16_PROTOCOL.json').read_text())
 return dict(status='FROZEN_PRE_GT_EXPLORATORY_STRONG_BASELINE_ADDENDUM',stage='IV',job=JOB,parents=PARENTS,lambda_arms=[0,1],
  DEV_algorithm_identity=checkpoint_identity(DEV_ASRUN),current_identical_algorithm_identity=checkpoint_identity(ALG),exact_algorithm_byte_comparison=True,
  unchanged_DEV_parameters={k:dev[k] for k in ('objectives','coordinates','subproblem','initial_radius','max_radius','accept_ratio_gt','shrink_if_ratio_lt','shrink_factor','expand_if_ratio_gt','expand_factor','expand_only_boundary_fraction','preflight_each_trial_RHS')},
  states_per_parent_arm=16,RHS_per_parent_arm=256,total_RHS_cap=32768,cached_initial_J='originalpaidlocked96+64context',cold_FullJ_RHS_ifnotcached=40,
  prior_and_noise_and_decoder_unchanged=True,target_identity=checkpoint_identity(OUT/'configs/TARGET.json'),max_wall=1800,
  post_method_freeze=True,pre_GT=True,exploratory_not_confirmatory_primary=True,old45risk26densityoutputs_unchanged=True,
  no_GT=True,no_mu=True,no_hyperparameter_search=True,no_extra_aftercaps=True,parent_authorization='explicitfirst64L1lambda0+1sameverifiedDEValgorithm32768RHS')

def run():
 if RESULT.exists():raise FileExistsError('noL1rerun')
 check_algorithm()
 if json.loads(PROTOCOL.read_text())!=json.loads(json.dumps(protocol())):raise ValueError('frozenaddendummismatch')
 b=Book('IV',JOB,max_wall=1800);original=b._update;cost=dict(RHS=0,states=0,LU=0);active=[None];status='COMPLETE';error=None;results=[]
 def tracked(event,**kw):
  if cost['RHS']+kw.get('rhs',0)>32768:raise BudgetStop('32768L1addendumcap')
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiBcap')
  if active[0] is not None:active[0]['require'](states=kw.get('states',0),RHS=kw.get('rhs',0))
  answer=original(event,**kw)
  for key,arg in (('RHS','rhs'),('states','states'),('LU','lu')):
   v=kw.get(arg,0);cost[key]+=v
   if active[0] is not None:active[0]['cost'][key]+=v
  return answer
 b._update=tracked
 try:
  RESULT.mkdir();write(RESULT/'FROZEN_PROTOCOL.json',protocol());t=Target(b,keep_fields=True)
  for parent in PARENTS:
   context=OUT/'results/contexts/locked'/f'{parent}.npz';cache=OUT/'results/contexts/cache'/f'{parent}.npz'
   with np.load(context,allow_pickle=False) as z:x0=z['x'].copy();Jq=z['J'].copy();y=z['observation'].copy()
   for lam in (0,1):
    folder=RESULT/parent/f'lambda{lam}';folder.mkdir(parents=True);guard=dict(cost=dict(RHS=0,states=0,LU=0))
    def require(*,states,RHS):
     if guard['cost']['states']+states>16 or guard['cost']['RHS']+RHS>256:raise ControlBudgetStop('sameDEV16states/256RHSpreflight')
    guard['require']=require;active[0]=guard;b.context=dict(parent=parent,lambda_=lam,role='preGTexploratoryFullJTRaddendum')
    anchor=t.load_cache(cache,owner='paidlocked96+64context,cachedinitialJcold40RHSnotfree')
    if np.max(abs(from_probit(x0)-anchor.q))>1e-12 or Jq.shape!=(512,66) or anchor.T is None:raise ValueError('exactcachednativeJanchorrequired')
    C=np.einsum('qp,ps->qsp',anchor.T,anchor.state.field).reshape(256,576);a=C@t.prior.derivative(anchor.q);check=np.column_stack([t.white(a[:,k].reshape(32,8)) for k in range(66)])
    if not np.allclose(check,Jq,rtol=1e-10,atol=1e-9):raise ValueError('sameDEVcachedJrecompositionfailure')
    anchor.J=Jq.copy();phi0=t.phi(anchor,y);obj0=phi0+.5*lam*float(x0@x0)
    def event(row):
     with (folder/'TRIAL_RATIOS.jsonl').open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
    def derivative(p,x,J,g,regularization,role):save(folder/'derivatives'/f'{p.key}.npz',x=x,q=p.q,Jcoordinate=J,gradient_total=g,lambda_=np.array(regularization),role=np.array(role))
    result,point=solve_control(t,anchor,x0,y,regularization=lam,require=require,used=lambda:guard['cost'],event=event,save_derivative=derivative)
    save(folder/'FINAL.npz',x=result['x'],q=point.q,chi=point.state.chi,current=point.state.current,field=point.state.field,prediction=point.state.data,Phi=np.array(result['Phi']))
    result.update(parent=parent,anchor_Phi=phi0,anchor_objective=obj0,cached_initial_J_source=checkpoint_identity(context),cached_anchor_source=checkpoint_identity(cache),cold_FullJ_setup_RHS=40,
     pre_GT=True,exploratory_addendum_not_primary=True,old45risk26densityoutputs_unchanged=True,GT_read=False,mu_used=False)
    write(folder/'RESULTS.json',result);results.append(result);active[0]=None;t.backend()._transfers.clear()
    write(RESULT/'PROGRESS.json',dict(completed_cases=len(results),actual_cost=cost));print(json.dumps(dict(parent=parent,lambda_=lam,cost=result['actual_cost'],accepted=result['accepted_updates'],stop=result['stop_reason'])),flush=True)
  write(RESULT/'RESULTS.json',dict(status='COMPLETE',cases=results,actual_cost=cost,GT_read=False,mu_used=False,pre_GT=True,exploratory_addendum_not_primary=True,oldoutputs_unchanged=True))
 except BaseException as exc:status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc);write(RESULT/'ERROR.json',dict(error=error,traceback=traceback.format_exc(),actual_cost=cost));raise
 finally:
  b.close(status,error)
  if status=='COMPLETE':
   artifacts=[checkpoint_identity(RESULT/'RESULTS.json'),checkpoint_identity(PROTOCOL)]+[checkpoint_identity(RESULT/p/f'lambda{lam}/FINAL.npz') for p in PARENTS for lam in (0,1)]
   write(RESULT/'CHOICES_SEALED_NO_GT.json',dict(status='CLOSED_CHOICES_SEALED_BEFORE_GT',pid=os.getpid(),GT_read=False,artifacts=artifacts,receipt=checkpoint_identity(b.jobdir/'CLOSED.json'),role='EXPLORATORY_STRONG_BASELINE_NOT_PRIMARY'))

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));a=p.parse_args()
 if a.action=='prepare':
  if PROTOCOL.exists():raise FileExistsError('preservepreGTfreeze')
  write(PROTOCOL,protocol())
 else:run()
if __name__=='__main__':main()
