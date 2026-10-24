#!/usr/bin/env python3
"""Explicit exposed DEV16 Full-J TR control, new files only; no L1/GT/mu."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,argparse,time,traceback,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,Target,write,save,from_probit,BudgetStop
from p3r4.evaluation import checkpoint_identity
from p3r4.ordinary_trust_region import solve_control

JOB='ordinary_fullJ_TR_DEV16_v1';RESULT=OUT/'results/II'/JOB;PROTOCOL=OUT/'configs/ORDINARY_FULLJ_TR_DEV16_PROTOCOL.json'
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(16))
class ControlBudgetStop(RuntimeError):pass
def protocol():return dict(status='FROZEN_EXPOSED_DEV_ONLY',job=JOB,stage='II',parents=PARENTS,lambda_arms=[0,1],
 objectives=['Phi=.5||sqrt2/sigma*pack(F-y)||²','Phi+.5||x66||²'],coordinates='originalfull66probit,exactnativeDandnormalPDFchain',
 subproblem='exactPSD-eigen/secular min g.s+.5s.Hs withinEuclideanradius',initial_radius=.5,max_radius=2.,
 accept_ratio_gt=.1,shrink_if_ratio_lt=.25,shrink_factor=.25,expand_if_ratio_gt=.75,expand_factor=2.,expand_only_boundary_fraction=.99,
 states_per_parent_arm=16,RHS_per_parent_arm=256,total_Maxwell_RHS_cap=8192,accepted_update_J_RHS=32,trialFull_RHS=8,
 preflight_each_trial_RHS=40,initial_J='savedpaidnativewhitenedJ atsameoriginalanchor; no new query',cold_anchor_FullJ_RHS_ifnotcached=40,
 initial_context_owner='originalpaidDEVnative96/probit64contexts,actualperparentcostsretained',max_wall=1200,
 no_GT=True,no_mu=True,no_neural_updates=True,old_representative_or_models_edited=False,new_L1_control_authorized=False,
 purpose='preGTdiagnosticaddendum,notretroactiveprimaryclaim',parameter_search=False)

def run():
 if RESULT.exists():raise FileExistsError('no rerun')
 if json.loads(PROTOCOL.read_text())!=json.loads(json.dumps(protocol())):raise ValueError('frozenprotocolmismatch')
 b=Book('II',JOB,max_wall=1200);original=b._update;cost=dict(RHS=0,states=0,LU=0);active=[None];status='COMPLETE';error=None;results=[]
 def tracked(event,**kw):
  if cost['RHS']+kw.get('rhs',0)>8192:raise BudgetStop('8192globalDEVcap')
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiBordinarycontrolcap')
  if active[0] is not None:active[0]['require'](states=kw.get('states',0),RHS=kw.get('rhs',0))
  result=original(event,**kw)
  for k,a in (('RHS','rhs'),('states','states'),('LU','lu')):
   v=kw.get(a,0);cost[k]+=v
   if active[0] is not None:active[0]['cost'][k]+=v
  return result
 b._update=tracked
 try:
  RESULT.mkdir();t=Target(b,keep_fields=True);write(RESULT/'FROZEN_PROTOCOL.json',protocol())
  for parent in PARENTS:
   cp=OUT/'results/contexts/dev'/f'{parent}.npz';cache=OUT/'results/contexts/cache'/f'{parent}.npz'
   with np.load(cp,allow_pickle=False) as z:x0=z['x'].copy();Jq=z['J'].copy();y=z['observation'].copy()
   for lam in (0,1):
    folder=RESULT/parent/f'lambda{lam}';folder.mkdir(parents=True);guard=dict(cost=dict(RHS=0,states=0,LU=0))
    def require(*,states,RHS):
     if guard['cost']['states']+states>16 or guard['cost']['RHS']+RHS>256:raise ControlBudgetStop('16states/256RHSpreflight')
    guard['require']=require;active[0]=guard;b.context=dict(parent=parent,lambda_=lam,role='newordinaryFullJTRcontrol,noGTormu')
    anchor=t.load_cache(cache,owner='originalpaidDEV96+64context;cachedFullJcold40notfree')
    if np.max(abs(from_probit(x0)-anchor.q))>1e-12 or Jq.shape!=(512,66):raise ValueError('savednativeJprobitanchor mismatch')
    # Validate supplied paid J by array composition with the paid T0/E0.
    if anchor.T is None:raise ValueError('missingpaidanchorT,noautomaticnewJsetup')
    C=np.einsum('qp,ps->qsp',anchor.T,anchor.state.field).reshape(256,576);a=C@t.prior.derivative(anchor.q)
    recomposed=np.column_stack([t.white(a[:,k].reshape(32,8)) for k in range(66)])
    if not np.allclose(recomposed,Jq,rtol=1e-10,atol=1e-9):raise ValueError('cachedJrecompositionfailure')
    anchor.J=Jq.copy();phi0=t.phi(anchor,y);obj0=phi0+.5*lam*float(x0@x0)
    def event(r):
     with (folder/'TRIAL_RATIOS.jsonl').open('a') as f:f.write(json.dumps(r,allow_nan=False)+'\n')
    def derivative(p,x,J,g,regularization,role):save(folder/'derivatives'/f'{p.key}.npz',x=x,q=p.q,Jcoordinate=J,gradient_total=g,lambda_=np.array(regularization),role=np.array(role))
    result,point=solve_control(t,anchor,x0,y,regularization=lam,require=require,used=lambda:guard['cost'],event=event,save_derivative=derivative)
    save(folder/'FINAL.npz',x=result['x'],q=point.q,chi=point.state.chi,current=point.state.current,field=point.state.field,prediction=point.state.data,Phi=np.array(result['Phi']))
    result.update(parent=parent,anchor_Phi=phi0,anchor_objective=obj0,cached_initial_J_source=checkpoint_identity(cp),cached_anchor_source=checkpoint_identity(cache),
     original_context_cost_source=str(cp.with_suffix('.json')),cold_FullJ_setup_RHS=40,GT_read=False,mu_used=False)
    write(folder/'RESULTS.json',result);results.append(result);active[0]=None;t.backend()._transfers.clear()
    write(RESULT/'PROGRESS.json',dict(completed_cases=len(results),actual_cost=cost));print(json.dumps(dict(parent=parent,lambda_=lam,cost=result['actual_cost'],accepted=result['accepted_updates'],stop=result['stop_reason'])),flush=True)
  write(RESULT/'RESULTS.json',dict(status='COMPLETE',cases=results,actual_cost=cost,GT_read=False,mu_used=False,L1_control_launched=False,old_runs_modified=False))
 except BaseException as exc:status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc);write(RESULT/'ERROR.json',dict(error=error,traceback=traceback.format_exc(),actual_cost=cost));raise
 finally:b.close(status,error)

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));a=p.parse_args()
 if a.action=='prepare':
  if PROTOCOL.exists():raise FileExistsError('preservefrozenprotocol')
  write(PROTOCOL,protocol());print('FullJTRDEV16protocolfrozen,no launch.')
 else:run()
if __name__=='__main__':main()
