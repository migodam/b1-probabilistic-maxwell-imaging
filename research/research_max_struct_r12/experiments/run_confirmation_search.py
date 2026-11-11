"""Frozen DR / direct multistart / NN-direct, public input only, whole-child CPU caps."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,argparse,json,time
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
if os.environ.get('R12_LEGACY_REPO'):sys.path.append(os.environ['R12_LEGACY_REPO'])
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
from scipy.stats import chi2
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection

class PhaseEnd(Exception):pass
def load(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def main():
 p=argparse.ArgumentParser();p.add_argument('--public-dir',type=Path,required=True);p.add_argument('--case',required=True)
 p.add_argument('--method',choices=['dr','multistart','nn_direct'],required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--cpu',type=int);p.add_argument('--weights',type=Path)
 a=p.parse_args()
 if a.cpu is not None:
  import psutil
  psutil.Process().cpu_affinity([a.cpu])
 secondary_enabled=a.public_dir.name=='confirm24'
 secondary_freeze=a.public_dir.parents[1]/'SECONDARY_SELECTION_FREEZE.json'
 if secondary_enabled and not secondary_freeze.is_file():
  raise ValueError('Registered secondary selector freeze required before confirmation inverse dispatch')
 out=a.output;out.mkdir(parents=True,exist_ok=False);wall=time.perf_counter()
 public=load(a.public_dir/'operators.npz');public.update(load(a.public_dir/(a.case+'.npz')))
 model=MaxwellModel(public);Y=public['Y'];m0=np.full(model.n,.1)
 sigma=float(public['sigma']);delta=sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
 config=dict(case=a.case,method=a.method,scale=6.,omega=1.,steps=3000,audit_every=50,grouped_ports=True,
  whole_child_CPU_cap=240.,proposal_cap=144.,corrector_allowance=96.,box=[.02,20.],delta_raw=delta,
  selection='minimum original Maxwell Phi plus inherited prior, only completed paid cap-eligible candidates',
  corrector=dict(maxiter=3000,maxfun=12000,scale=1e6,ftol=1e-12,gtol=.01,maxls=30),truth_access='NONE')
 if secondary_enabled:config.update(secondary_selection='CN99_MIN_SAME_PRIOR_AUDIT_ONLY_V1',secondary_freeze_receipt=str(secondary_freeze))
 (out/'CONFIG.json').write_text(json.dumps(config,indent=2))
 extra={}
 if a.method=='nn_direct':
  t=time.process_time()
  import torch
  torch.set_num_threads(1)
  from research.research_max_struct_r8.src.warm_features import FeatureOperator
  from research.research_max_struct_r8.src.warm_network import load_predictor
  features,_,featurecost=FeatureOperator(public)(Y)
  predictor,checkpoint=load_predictor(a.weights,device='cpu')
  m0=np.clip(predictor(features)[0].ravel().astype(float),.02,20.)
  extra=dict(nn_cpu=time.process_time()-t,featurecost=featurecost,prior_information='frozen R8 supervised model')
 history=[];materials=[];best=[float('inf'),m0.copy()];prev={};eligible_deadline=144.
 secondary_best=None;secondary_incumbents=[];secondary_audits=[]
 def audit(m):
  nonlocal prev,secondary_best
  st=model.state(m);phi=model.physical_score(m,Y);clock=time.process_time()
  counts=dict(model.cost);increments={k:v-prev.get(k,0) for k,v in counts.items()};prev=counts
  metrics=dict(relative_data_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),
   prior=float(model.reg.evaluate(m)[0]),noise_feasible=bool(np.linalg.norm(st['Y']-Y)<=delta))
  eligible=clock<=eligible_deadline
  if eligible and phi<best[0]:best[:]=[phi,m.copy()]
  history.append(dict(candidate=len(materials),CPU_seconds=clock,Phi=phi,eligible=eligible,**metrics,cost=counts))
  materials.append(m.copy())
  if secondary_enabled:
   # Only already computed full physical state and the same existing R metric.
   # This side observer never changes best, m, a step, or an optimizer return.
   legal=bool(np.isfinite(m).all() and np.min(m)>=.02 and np.max(m)<=20.)
   score=(float(metrics['prior']),float(phi))
   proposed=legal and metrics['noise_feasible'] and (secondary_best is None or score<secondary_best['score'])
   snapshot=dict(m=m.copy(),j=st['j'].copy(),Q=st['Q'].copy(),Y_prediction=st['Y'].copy()) if proposed else None
   raw_residual_norm=float(np.linalg.norm(st['Y']-Y))
   completed=time.process_time();secondary_eligible=completed<=eligible_deadline
   rec=dict(candidate=len(materials)-1,CPU_seconds=completed,
    phase='proposal' if eligible_deadline==144. else 'corrector',R=score[0],Phi=score[1],
    raw_residual_norm=raw_residual_norm,delta_CN99=float(delta),
    noise_feasible=metrics['noise_feasible'],legal_material=legal,phase_cap_eligible=secondary_eligible,
    phase_deadline_CPU_seconds=eligible_deadline)
   secondary_audits.append(rec)
   if proposed and secondary_eligible:
    secondary_best=dict(score=score,row=snapshot,metadata=rec.copy())
    secondary_incumbents.append(dict(rec,m=snapshot['m']))
  return PhysicalAudit(phi,st['j'].copy(),st['Q'].copy(),metrics,increments)
 def optimize(m,deadline,maxiter=3000):
  def fun(x):
   if time.process_time()>=deadline:raise PhaseEnd()
   v,g=model.data_value_gradient(x,Y);audit(x)
   if time.process_time()>=deadline:raise PhaseEnd()
   return v*1e6,g*1e6
  try:
   z=minimize(fun,m,jac=True,method='L-BFGS-B',bounds=[(.02,20.)]*model.n,
     options=dict(maxiter=maxiter,maxfun=12000,ftol=1e-12,gtol=.01,maxls=30))
   return str(z.message)
  except PhaseEnd:return 'CPU_PHASE_LIMIT'
 audit(m0);run=None;termination=[];setup={}
 if a.method=='dr':
  pr=DenseFieldCurrentProjectors(model.G,model.E,model.S,Y,delta=delta,scale=6.,lower=.02,upper=20.,grouped_ports=True)
  run=run_projection('douglas_rachford',pr,m0,audit,steps=3000,audit_every=50,cpu_budget_seconds=max(0.,144-time.process_time()))
  setup=dict(projection_total_counts=dict(pr.costs),setup_cpu=pr.setup_cpu_seconds,setup_wall=pr.setup_wall_seconds)
  termination.append(run.stop_reason)
 elif a.method=='multistart':
  termination.append(optimize(m0,144.))
  seed=0
  while time.process_time()<144.:
   v=gaussian_filter(np.random.default_rng(120120+seed).normal(size=model.shape),2.).ravel();v/=np.sqrt(np.mean(v*v))
   z=np.clip(m0+[1.,4.,10.][seed%3]*v,.02,20.);seed+=1
   termination.append(optimize(z,min(144.,time.process_time()+24.),1000))
 else:termination.append('NN_FEATURES_AND_INFERENCE')
 proposal_CPU=time.process_time();eligible_deadline=min(240.,proposal_CPU+96.)
 termination.append(optimize(best[1].copy(),eligible_deadline))
 # Frozen candidates are sealed before offline scoring; finalization adds no forward.
 np.savez_compressed(out/'prediction.npz',m=best[1],initial_material=m0,candidates=np.asarray(materials))
 (out/'history.json').write_text(json.dumps(history,indent=2))
 prefix={};prefixrows=[]
 for cap in [30,60,120,240]:
  candidates=[r for r in history if r['eligible'] and r['CPU_seconds']<=cap]
  if candidates:
   r=min(candidates,key=lambda x:x['Phi']);prefix[f'm_{cap}']=materials[r['candidate']]
   prefixrows.append(dict(cap=cap,candidate=r['candidate'],CPU=r['CPU_seconds'],Phi=r['Phi']))
 np.savez_compressed(out/'paid_prefix.npz',**prefix)
 if run:run.save(out/'projection')
 if secondary_enabled:
  info=dict(status='SELECTED_NOISE_FEASIBLE_MIN_SAME_PRIOR' if secondary_best is not None else 'NO_CAP_ELIGIBLE_NOISE_FEASIBLE_PAID_STATE',
   rule='Legal realbox material; raw residual<=globalCN99 radius; minimize original R(m), exact tie by original Phi',
   primary_selection='UNCHANGED_MIN_ORIGINAL_PHI',feedback_into_iterations='NONE',new_forward=0,new_GT=0,
   audits=len(secondary_audits),freeze_receipt=str(secondary_freeze),
   estimator_scope='Pre-registered secondary candidate selector, not a pure DR mechanism benefit')
  if secondary_best is not None:
   np.savez_compressed(out/'NOISE_FEASIBLE_MIN_PRIOR.npz',**secondary_best['row'])
   info['selection']=secondary_best['metadata']
  (out/'NOISE_FEASIBLE_MIN_PRIOR.json').write_text(json.dumps(info,indent=2))
  (out/'NOISE_FEASIBLE_MIN_PRIOR_HISTORY.json').write_text(json.dumps(secondary_audits,indent=2))
  secondary_prefix={};secondary_prefixrows=[]
  for cap in (30,60,120,240):
   choices=[r for r in secondary_incumbents if r['CPU_seconds']<=cap]
   if choices:
    selected=min(choices,key=lambda r:(r['R'],r['Phi']));secondary_prefix['m_'+str(cap)]=selected['m']
    secondary_prefixrows.append(dict(cap=cap,**{k:v for k,v in selected.items() if k!='m'}))
   else:secondary_prefixrows.append(dict(cap=cap,status='NO_CAP_ELIGIBLE_NOISE_FEASIBLE_PAID_STATE'))
  np.savez_compressed(out/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.npz',**secondary_prefix)
  (out/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.json').write_text(json.dumps(secondary_prefixrows,indent=2))
 receipt=dict(status='COMPLETE',case=a.case,method=a.method,objective=best[0],termination=termination,
  proposal_CPU=proposal_CPU,effective_deadline=eligible_deadline,whole_child_CPU_seconds=time.process_time(),
  wall_seconds=time.perf_counter()-wall,physics_cost=dict(model.cost),projection=setup,extra=extra,
  prefixes=prefixrows,truth_access='NONE',cost_limit='Serialization/imports charged; atomic overrun excluded from candidate selection')
 (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
