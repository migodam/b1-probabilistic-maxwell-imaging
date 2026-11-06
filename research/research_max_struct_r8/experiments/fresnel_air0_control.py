"""Post-measurement known-air feasible-set diagnostic; only lowerbound changes.

No A/A+ readout, no new pixeltruth, no frequency/weight/prior selection.
The8GHz anchor/holdout is already exposed;6GHz heldout only scored at endpoints.
"""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import sys,argparse,time,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from research.research_max_struct_r8.src.fresnel_adapter import prepare,BoundedMatrixFreeModel
from fresnel_invert import dump

def main():
 args=argparse.Namespace(mode='air0_control_8')
 out=ROOT/'RESULTS/fresnel'/'FoamTwinDielTM_air0_control_8GHz_dev';out.mkdir(exist_ok=False);cpu=time.process_time();wall=time.perf_counter()
 anchor_path=ROOT/'RESULTS/fresnel/FoamTwinDielTM_8GHz_n60_cie_dual_v3/derived.npz';anchor=np.load(anchor_path)['m'].copy();best=None;history=[]
 with threadpool_limits(limits=1):
  datasets={f:prepare(ROOT/'data_local/fresnel/FoamTwinDielTM.exp',frequency=f,side=60,width=.2,matrix_free=True) for f in [6,8]}
  if not np.array_equal(datasets[6]['train_mask'],datasets[8]['train_mask']):raise ValueError('Frequency masks differ')
  if len(anchor)!=3600 or not np.isfinite(anchor).all() or anchor.min()<.02 or anchor.max()>20:raise ValueError('Anchor outside originalbox')
  active=[8]
  frozen=dict(mode=args.mode,training_frequencies_GHz=active,validation_frequencies_GHz=[6,8],shared_real_material=True,side=60,width_m=.2,box=[0.,20.],
   objective='mean_f[0.5norm(F_f-Y_f)^2/norm(Y_f)^2]+originalprior; each normalized frequency equal weight; prior mathematically once',
   prior='original Regularizer smooth1e-4 tv1e-4 epsilon.01, unchanged',anchor=str(anchor_path.relative_to(ROOT)),anchor_selection='8GHzDualv3bestpaidtrainingphysicalobjective, noGT',
   maximum_accepted_iterations=300,maximum_function_evaluations=1000,corrector=dict(method='original L-BFGS-B',scale=1e6,ftol=1e-12,gtol=.01,maxls=30),
   mask='Same source/receiver actualmask asfinal8GHzfivearm',training_channels_per_frequency=2496,heldout_channels_per_frequency=1842,heldout_selection=False,
   observation_status='DEV_MODEL_DIAGNOSTIC_POST_MEASUREMENT; same8GHzdata; onlyfeasiblelowerbound0instead.02; notA+advantage',diagnostic_status='DEV_MODEL_DIAGNOSTIC',old_control='RESULTS/fresnel/FoamTwinDielTM_multifrequency_control_8/RECEIPT.json',only_optimizer_change='Material lowerbound0instead.02',background_region='max(abs(x),abs(y))>=.085m outerROIband, posthocdistributiononly, notpixelGT/supportconstraint',distribution_thresholds=dict(bound_atol=1e-10,old_floor=.02,near_zero=1e-8,quantiles=[0,.05,.25,.5,.75,.95,1]),eight_GHz_heldout_status='PREVIOUSLY_EXPOSED_FINAL8GHZRESULT',six_GHz_heldout_status='PREVIOUSLY_EXPOSED_BY6PLUS8FOLLOWUP; no new test seal',frequency_sweep=False,
   illumination='Independent target-absent source multipole calibration at each physicalfrequency',physics='SeparatefullVIEperfrequency; sharematerialonly; noMobiustransplant')
  dump(out/'FROZEN_CONFIG.json',frozen)
  for f in datasets:dump(out/f'INCIDENT_CALIBRATION_{f}GHz.json',datasets[f]['calibration'])
  setup_cpu=time.process_time()-cpu;optimization_start=time.process_time();best=[anchor.copy(),float('inf')]
  def obj(m):
   values=[];grads=[]
   for f in active:
    d=datasets[f];v,g=d['model'].data_value_gradient(m,d['Y_train']);values.append(v);grads.append(g)
   # Each term contains the same original R(m), so the average counts R once.
   value=float(np.mean(values));gradient=np.mean(np.stack(grads),axis=0)
   projected=np.where(((m<=0.)&(gradient>0))|((m>=20)&(gradient<0)),0.,gradient)
   history.append(dict(evaluation=len(history)+1,objective=value,per_frequency_objective_with_prior=dict(zip(map(str,active),map(float,values))),projected_gradient_inf=float(np.max(abs(projected))),cpu_seconds=time.process_time()-cpu))
   if value<best[1]:best[:]=[m.copy(),value];np.savez_compressed(out/'paid_incumbent.npz',m=m,objective=value,evaluation=len(history))
   if len(history)%10==0:dump(out/'PROGRESS.json',history[-1]);dump(out/'optimizer_history.json',history)
   return value*1e6,gradient*1e6
  try:
   rr=minimize(obj,anchor,method='L-BFGS-B',jac=True,bounds=[(0.,20.)]*3600,options=dict(maxiter=300,maxfun=1000,ftol=1e-12,gtol=.01,maxls=30))
   optimization_cpu=time.process_time()-optimization_start;metrics={};predictions={};prediction_cost={}
   for f,d in datasets.items():
    full=BoundedMatrixFreeModel(d['public']);st=full.state(best[0]);p=st['Y'];predictions[str(f)]=p;prediction_cost[str(f)]=dict(full.cost)
    sh=np.zeros(d['Y'].shape,bool);sh[:,d['held_sources']]=True;sh &=np.isfinite(d['Y']);rh=d['withheld_mask']&~sh
    rel=lambda mask:float(np.linalg.norm((p-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))
    metrics[str(f)]=dict(train=rel(d['train_mask']),heldout_all=rel(d['withheld_mask']),heldout_sources=rel(sh),heldout_receivers=rel(rh),state_relative_residual=st['state_relative_residual'])
   m=best[0];h=.2/60;foreground=np.maximum(m-.02,0.);concentration=dict(max_contrast=float(m.max()),min_contrast=float(m.min()),contrast_integral_m2=float(m.sum()*h*h),foreground_effective_area_m2=float(h*h*foreground.sum()**2/max(float(np.sum(foreground**2)),1e-300)),prior=float(datasets[8]['model'].reg.evaluate(m)[0]))
   background=np.max(np.abs(datasets[8]['centers']),axis=1)>=.085
   def distribution(x,lower):
    return dict(count=int(len(x)),mean=float(np.mean(x)),std=float(np.std(x)),quantiles=dict(zip(map(str,[0,.05,.25,.5,.75,.95,1]),map(float,np.quantile(x,[0,.05,.25,.5,.75,.95,1])))),fraction_at_declared_lower_bound=float(np.mean(np.abs(x-lower)<=1e-10)),fraction_below_old_0p02=float(np.mean(x<.02-1e-10)),fraction_near_zero=float(np.mean(x<=1e-8)))
   concentration['whole_roi_distribution']=distribution(m,0.);concentration['outer_roi_band_distribution']=distribution(m[background],0.);concentration['effective_area_raw_m2']=float(h*h*m.sum()**2/max(float(np.sum(m*m)),1e-300));concentration['foreground_area_definition']='h2*(sum max(m-.02,0))2/sum max(m-.02,0)2; common old-floor excess for both controls, notGTsupport'
   r=dict(status='COMPLETED_DEV_MODEL_DIAGNOSTIC',mode=args.mode,metrics=metrics,mean_6_8_heldout_relative=float(np.sqrt(np.mean([v['heldout_all']**2 for v in metrics.values()]))),material_diagnostics=concentration,
    corrector=dict(iterations=int(rr.nit),evaluations=int(rr.nfev),success=bool(rr.success),message=str(rr.message)),objective=float(best[1]),setup_cpu_seconds=setup_cpu,optimization_cpu_seconds=optimization_cpu,
    cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost_optimizer={str(f):dict(datasets[f]['model'].cost) for f in active},physical_cost_validation=prediction_cost,
    qualification='Post-measurementDEV_MODEL_DIAGNOSTIC; onlylowerboundchanged; no newobservations; noGTmetrics; notA+advantage; bothfrequencyholdoutsalreadyexposed',no_test_refitting=True)
   dump(out/'RECEIPT.json',r);dump(out/'optimizer_history.json',history);np.savez_compressed(out/'derived.npz',m=m,anchor_m=anchor,prediction6=predictions['6'],prediction8=predictions['8'],centers=datasets[8]['centers'],train_mask=datasets[8]['train_mask'],withheld_mask=datasets[8]['withheld_mask']);print(json.dumps(r),flush=True)
  except BaseException as exc:
   if best is not None:np.savez_compressed(out/'paid_incumbent.npz',m=best[0],objective=best[1])
   dump(out/'optimizer_history.json',history);dump(out/'FAILED.json',dict(error=repr(exc),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost={str(f):dict(datasets[f]['model'].cost) for f in active}));raise
if __name__=='__main__':main()
