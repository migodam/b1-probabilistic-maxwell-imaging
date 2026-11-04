from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,argparse
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from research.research_max_struct_r7.src.fresnel_adapter import prepare
from research.research_max_struct_r6.src.benchmark_solver import arm_config,run_auxiliary
from research.research_fbg_r3_experiments.src.physics import MaxwellModel


def dump(p,x):
 def cv(v):
  if isinstance(v,np.ndarray):return v.tolist()
  if isinstance(v,np.generic):return v.item()
  raise TypeError(type(v).__name__)
 p.write_text(json.dumps(x,indent=2,default=cv,allow_nan=False)+'\n')


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--arm',required=True);args=ap.parse_args()
 out=ROOT/'RESULTS'/f'fresnel_{args.case}_{args.arm}';out.mkdir(exist_ok=False)
 cpu=time.process_time();wall=time.perf_counter()
 with threadpool_limits(limits=1):
  d=prepare(ROOT/'data_public'/'fresnel'/(args.case+'.exp'));model=d['model'];Y=d['Y_train']
  common=dict(host=.1,material_box=[.02,20.],dpcs_rtol=.001)
  cfg=arm_config(args.arm,common,Y,model.E,[dict(mf=12,maxiter=1500)])
  cfg.update(regularizer=lambda m:model.reg.evaluate(m)[:2],physical_sigma=float(np.linalg.norm(Y)))
  dump(out/'FROZEN_CONFIG.json',dict(case=args.case,arm=args.arm,frequency_GHz=2.,side=24,width_m=.16,
    training_sources_onebased=(d['train_sources']+1).tolist(),training_receivers_onebased=(d['train_receivers']+1).tolist(),
    config={k:v for k,v in cfg.items() if k not in ('regularizer','data_weights','state_weights')},
    corrector=dict(maxiter=1000,maxfun=4000,ftol=1e-12,gtol=.01,maxls=30,scale=1e6),
    initialization=.1,selection='training original normalized physical data loss plus inherited prior only',
    heldout_selection=False,pixel_GT='UNAVAILABLE',illumination='target-absent on-axis cylindrical approximation'))
  aux=run_auxiliary(model,Y,np.full(model.n,.1),args.arm,cfg)
  m=aux['m'];history=[];best=[m.copy(),model.physical_score(m,Y)]
  aux_cpu=time.process_time()-cpu
  def obj(x):
   v,g=model.data_value_gradient(x,Y)
   history.append(dict(evaluation=len(history)+1,objective=float(v),cpu_seconds=time.process_time()-cpu))
   if v<best[1]:best[:]=[x.copy(),float(v)]
   return v*1e6,g*1e6
  rr=minimize(obj,m,method='L-BFGS-B',jac=True,bounds=[(.02,20.)]*model.n,options=dict(maxiter=1000,maxfun=4000,ftol=1e-12,gtol=.01,maxls=30))
  full=MaxwellModel(d['public']);state=full.state(best[0]);pred=state['Y']
  rel=lambda mask:float(np.linalg.norm((pred-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))
  sourceheld=np.isfinite(d['Y']);sourceheld[:,:9]=False
  rxheld=d['withheld_mask'].copy();rxheld[:,9:]=False
  receipt=dict(status='COMPLETED_EXPLORATORY_MEASURED_INVERSION',case=args.case,arm=args.arm,
   train_relative_complex_residual=rel(d['train_mask']),withheld_all_relative_complex_residual=rel(d['withheld_mask']),
   withheld_sources_relative_complex_residual=rel(sourceheld),withheld_receivers_relative_complex_residual=rel(rxheld),
   training_channels=int(d['train_mask'].sum()),withheld_channels=int(d['withheld_mask'].sum()),
   training_objective=float(best[1]),state_relative_residual=float(state['state_relative_residual']),
   cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,auxiliary_cpu_seconds=aux_cpu,
   auxiliary=aux['execution_receipt'],corrector=dict(message=str(rr.message),iterations=int(rr.nit),evaluations=int(rr.nfev),success=bool(rr.success)),
   physical_cost_training=dict(model.cost),physical_cost_heldout_prediction=dict(full.cost),
   min_contrast=float(best[0].min()),max_contrast=float(best[0].max()),
   GT_metrics='UNAVAILABLE; published nominal geometry/permittivity is not a pixel GT',
   qualification='PDE solution verified by original discrete state residual; horn/antenna approximation not independently calibrated',
   no_test_refitting=True)
  np.savez_compressed(out/'endpoint.npz',m=best[0],prediction=pred,measured_Y=d['Y'],train_mask=d['train_mask'],withheld_mask=d['withheld_mask'],centers=d['centers'])
  dump(out/'RECEIPT.json',receipt);dump(out/'corrector_history.json',history)
  dump(out/'auxiliary_history.json',aux['history'])
  print(json.dumps({k:receipt[k] for k in ('case','arm','train_relative_complex_residual','withheld_all_relative_complex_residual','cpu_seconds')}),flush=True)
if __name__=='__main__':main()
