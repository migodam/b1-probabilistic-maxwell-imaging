from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,argparse
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from research.research_max_struct_r8.src.fresnel_adapter import prepare
from research.research_max_struct_r8.src.fresnel_solver import run
from research.research_fbg_r3_experiments.src.physics import MaxwellModel


def dump(p,x):
 def cv(v):
  if isinstance(v,np.ndarray):return v.tolist()
  if isinstance(v,np.generic):return v.item()
  raise TypeError(type(v).__name__)
 p.write_text(json.dumps(x,indent=2,default=cv,allow_nan=False)+'\n')


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--arm',required=True);ap.add_argument('--illumination',default='multipole');ap.add_argument('--schedule',default='full');ap.add_argument('--attempt',default='v2');ap.add_argument('--initialization',default='uniform',choices=['uniform','neural']);args=ap.parse_args()
 tag='_nn' if args.initialization=='neural' else ''
 out=ROOT/'RESULTS/fresnel'/f'{args.case}_{args.arm}_{args.illumination}_{args.schedule}_{args.attempt}{tag}';out.mkdir(exist_ok=False)
 cpu=time.process_time();wall=time.perf_counter()
 with threadpool_limits(limits=1):
  d=prepare(ROOT/'data_local/fresnel'/(args.case+'.exp'),illumination=args.illumination);model=d['model'];Y=d['Y_train']
  stages=[(3,1.,500),(6,.2,500),(10,.2,500),(12,.2,1500)] if args.schedule=='full' else [(3,1.,100),(6,.2,100),(10,.2,100),(12,.2,300)]
  cfg=dict(host=.1,material_box=[.02,20.],dpcs_rtol=.001,coefficient_tol=0.,miniter=50,
   material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,
   progress_path=str(out/'PROGRESS.json'),regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,aplus_sweeps=2,aplus_arc_candidates=True,aplus_inherited_gtol=True,
   stages=[dict(mf=mf,beta_material=bm,beta_current=bm if args.arm=='cie_single' else 1.,maxiter=it) for mf,bm,it in stages])
  m0=np.full(model.n,.1);nn_receipt=None
  if args.initialization=='neural':
   from research.research_max_struct_r8.src.fresnel_features import features
   from research.research_max_struct_r8.src.warm_network import load_predictor
   import torch
   torch.set_num_threads(1);nn_cpu=time.process_time();xx,nn_receipt=features(model,Y,d['centers']);predict,checkpoint=load_predictor(ROOT/'RESULTS/training/physics/weights.pt');m0=np.clip(predict(xx)[0].ravel().astype(float),.02,20.)
   nn_receipt.update(total_inference_setup_cpu_seconds=time.process_time()-nn_cpu,weights='RESULTS/training/physics/weights.pt',frozen_validation_epoch=int(checkpoint.get('selected_epoch',-1)),initial_min=float(m0.min()),initial_max=float(m0.max()))
   dump(out/'NN_RECEIPT.json',nn_receipt);np.savez_compressed(out/'nn_derived.npz',features=xx,m0=m0)
  frozen=dict(case=args.case,arm=args.arm,frequency_GHz=2.,side=24,width_m=.16,initialization=args.initialization,initialization_receipt=nn_receipt,
   training_channels=int(d['train_mask'].sum()),heldout_channels=int(d['withheld_mask'].sum()),training_sources_onebased=(d['train_sources']+1).tolist(),
   config={k:v for k,v in cfg.items() if k not in ('regularizer','aplus_regularizer')},
   corrector=dict(maxiter=1000,maxfun=4000,ftol=1e-12,gtol=.01,maxls=30,scale=1e6),
   selection='training original normalized physical data loss plus inherited prior only',heldout_selection=False,pixel_GT='UNAVAILABLE',illumination=args.illumination,
   receiver_contract='Per-source SVD DPCS / projected native IDFT APCS and source-specific receive / adjoint; absent channels never observations')
  dump(out/'FROZEN_CONFIG.json',frozen)
  try:
   aux=run(model,Y,m0,args.arm,cfg);m=aux['m'];history=[];best=[m.copy(),model.physical_score(m,Y)]
   dump(out/'AUXILIARY_RECEIPT.json',aux['execution_receipt']);dump(out/'auxiliary_history.json',aux['history']);dump(out/'auxiliary_candidates.json',aux['candidates'])
   aux_cpu=time.process_time()-cpu
   def obj(x):
    v,g=model.data_value_gradient(x,Y);history.append(dict(evaluation=len(history)+1,objective=float(v),cpu_seconds=time.process_time()-cpu))
    if v<best[1]:best[:]=[x.copy(),float(v)]
    return v*1e6,g*1e6
   rr=minimize(obj,m,method='L-BFGS-B',jac=True,bounds=[(.02,20.)]*model.n,options=dict(maxiter=1000,maxfun=4000,ftol=1e-12,gtol=.01,maxls=30))
   full=MaxwellModel(d['public']);state=full.state(best[0]);pred=state['Y']
   rel=lambda mask:float(np.linalg.norm((pred-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))
   sh=np.zeros(d['Y'].shape,bool);sh[:,d['held_sources']]=True;sh &= np.isfinite(d['Y']);rh=d['withheld_mask']&~sh
   receipt=dict(status='COMPLETED_EXPLORATORY_MEASURED_INVERSION',case=args.case,arm=args.arm,initialization=args.initialization,
    train_relative_complex_residual=rel(d['train_mask']),withheld_all_relative_complex_residual=rel(d['withheld_mask']),
    withheld_sources_relative_complex_residual=rel(sh),withheld_receivers_relative_complex_residual=rel(rh),
    training_channels=int(d['train_mask'].sum()),withheld_channels=int(d['withheld_mask'].sum()),training_objective=float(best[1]),
    state_relative_residual=float(state['state_relative_residual']),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,auxiliary_cpu_seconds=aux_cpu,
    auxiliary=aux['execution_receipt'],corrector=dict(message=str(rr.message),iterations=int(rr.nit),evaluations=int(rr.nfev),success=bool(rr.success)),
    physical_cost_training=dict(model.cost),physical_cost_prediction=dict(full.cost),min_contrast=float(best[0].min()),max_contrast=float(best[0].max()),
    qualification='Source multipole illumination calibrated only to empty measurements; receive horn contamination unresolved. No pixel GT.',no_test_refitting=True)
   np.savez_compressed(out/'derived.npz',m=best[0],prediction=pred,train_mask=d['train_mask'],withheld_mask=d['withheld_mask'],centers=d['centers'],auxiliary_endpoint=aux['last_joint_m'])
   dump(out/'RECEIPT.json',receipt);dump(out/'corrector_history.json',history)
   print(json.dumps({k:receipt[k] for k in ('case','arm','train_relative_complex_residual','withheld_all_relative_complex_residual','cpu_seconds')}),flush=True)
  except Exception as exc:
   dump(out/'FAILED.json',dict(error=repr(exc),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost)));raise
if __name__=='__main__':main()
