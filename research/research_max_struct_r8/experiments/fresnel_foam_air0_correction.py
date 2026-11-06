"""Two authorized fresh-start DEV_MODEL_CORRECTION runs; only material lowerbound changes.
Original3000aux/1000corrector, .1initial, originalprior/calibration/data preserved.
"""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,argparse
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from research.research_max_struct_r8.src.fresnel_adapter import prepare
from research.research_max_struct_r8.src.fresnel_cie_fft import run_cie_fft
from research.research_max_struct_r8.src.fresnel_adapter import BoundedMatrixFreeModel as MatrixFreeModel
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
from fresnel_invert import dump

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--arm',required=True,choices=['cie_dual','a_plus']);ap.add_argument('--side',default=60,type=int,choices=[60]);ap.add_argument('--attempt',default='dev_air0');args=ap.parse_args()
 out=ROOT/'RESULTS/fresnel'/f'FoamTwinDielTM_8GHz_n{args.side}_{args.arm}_{args.attempt}';out.mkdir(exist_ok=False);cpu=time.process_time();wall=time.perf_counter()
 with threadpool_limits(limits=1):
  d=prepare(ROOT/'data_local/fresnel/FoamTwinDielTM.exp',frequency=8.,side=args.side,width=.2,matrix_free=True);model=d['model'];Y=d['Y_train'];m0=np.full(model.n,.1)
  cfg=dict(arm=args.arm,host=.1,material_box=[0.,20.],dpcs_rank=5,coefficient_tol=1e-3,miniter=2,
   material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,
   regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,aplus_sweeps=2,aplus_arc_candidates=True,aplus_inherited_gtol=True,
   fft_geometry=model.paper_geometry,physical_sigma=float(np.linalg.norm(Y)),cache_material_value_gradient=True,
   stages=[dict(mf=mf,beta_material=b,beta_current=b if args.arm=='cie_single' else 6.,maxiter=1000,coefficient_tol=1e-4 if mf==6 else 1e-6) for mf,b in [(6,6.),(10,2.),(16,.5)]],progress_path=str(out/'PROGRESS.json'))
  if args.arm=='a_plus':cfg['material_readout']=mobius_material_readout
  dump(out/'FROZEN_CONFIG.json',dict(case='FoamTwinDielTM',frequency_GHz=8.,side=args.side,width_m=.2,arm=args.arm,
   config={k:v for k,v in cfg.items() if k not in ('regularizer','aplus_regularizer','fft_geometry','material_readout')},
   stopping_contract='Original paper delta2D sqrt metric thresholds.01/.001/.001; inherited squared metric thresholds1e-4/1e-6/1e-6',corrector=dict(maxiter=1000,maxfun=4000,ftol=1e-12,gtol=.01,maxls=30,scale=1e6),
   model='Full equal-area scalar-TM VIE; FFT+boundedGMRES state/analytic-adjoint with1e-9 RHS qualification; exact retained denseLU fallback when needed',
   initialization=.1,diagnostic_status='DEV_MODEL_CORRECTION_POST_MEASUREMENT',only_change='materiallowerbound.02to0inallauxiliaryreadoutsandphysicalcorrector',original_positive_floor_preserved=True,original_budget='3000auxcap(3x1000),1000correctoraccepted/4000eval',heldout_status='Both6and8GHzheldoutpreviouslyexposed; no new blindtest',paper_initialization_gap='NIEpaper air0; correctedcomparison remainscommoninitial.1 box[0,20]',
   training_channels=int(d['train_mask'].sum()),withheld_channels=int(d['withheld_mask'].sum()),heldout_selection=False))
  dump(out/'INCIDENT_CALIBRATION.json',d['calibration'])
  try:
   aux=run_cie_fft(model,Y,m0,cfg);dump(out/'AUXILIARY_RECEIPT.json',aux['execution_receipt']);dump(out/'auxiliary_candidates.json',aux['candidates']);dump(out/'auxiliary_history.json',aux['history'])
   best=[aux['m'].copy(),model.physical_score(aux['m'],Y)];history=[];aux_cpu=time.process_time()-cpu
   def obj(x):
    v,g=model.data_value_gradient(x,Y);history.append(dict(evaluation=len(history)+1,objective=float(v),cpu_seconds=time.process_time()-cpu))
    if v<best[1]:
     best[:]=[x.copy(),v];np.savez_compressed(out/'paid_incumbent.npz',m=best[0],training_objective=best[1],evaluation=len(history),cpu_seconds=time.process_time()-cpu)
    if len(history)%10==0:
     dump(out/'CORRECTOR_PROGRESS.json',history[-1]);dump(out/'corrector_history.json',history)
    return v*1e6,g*1e6
   rr=minimize(obj,best[0],method='L-BFGS-B',jac=True,bounds=[(0.,20.)]*model.n,options=dict(maxiter=1000,maxfun=4000,ftol=1e-12,gtol=.01,maxls=30))
   full=MatrixFreeModel(d['public']);st=full.state(best[0]);pred=st['Y'];rel=lambda mask:float(np.linalg.norm((pred-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))
   sh=np.zeros(d['Y'].shape,bool);sh[:,d['held_sources']]=True;sh &=np.isfinite(d['Y']);rh=d['withheld_mask']&~sh
   # Posthoc6GHz prediction only, after the8GHz endpoint has been paid and fixed.
   d6=prepare(ROOT/'data_local/fresnel/FoamTwinDielTM.exp',frequency=6.,side=args.side,width=.2,matrix_free=True);full6=MatrixFreeModel(d6['public']);st6=full6.state(best[0]);pred6=st6['Y'];dump(out/'INCIDENT_CALIBRATION_6GHz.json',d6['calibration'])
   sh6=np.zeros(d6['Y'].shape,bool);sh6[:,d6['held_sources']]=True;sh6 &=np.isfinite(d6['Y']);rh6=d6['withheld_mask']&~sh6
   rel6=lambda mask:float(np.linalg.norm((pred6-d6['Y'])[mask])/np.linalg.norm(d6['Y'][mask]))
   posthoc6=dict(train=rel6(d6['train_mask']),heldout_all=rel6(d6['withheld_mask']),heldout_sources=rel6(sh6),heldout_receivers=rel6(rh6),state_relative_residual=float(st6['state_relative_residual']))
   m=best[0];h=.2/args.side;foreground=np.maximum(m-.02,0.);background=np.max(np.abs(d['centers']),axis=1)>=.085
   def distribution(x):
    return dict(count=int(len(x)),mean=float(np.mean(x)),std=float(np.std(x)),quantiles=dict(zip(map(str,[0,.05,.25,.5,.75,.95,1]),map(float,np.quantile(x,[0,.05,.25,.5,.75,.95,1])))),fraction_at_declared_lower_bound=float(np.mean(np.abs(x)<=1e-10)),fraction_below_old_0p02=float(np.mean(x<.02-1e-10)),fraction_near_zero=float(np.mean(x<=1e-8)))
   descriptors=dict(max_contrast=float(m.max()),min_contrast=float(m.min()),contrast_integral_m2=float(m.sum()*h*h),foreground_effective_area_m2=float(h*h*foreground.sum()**2/max(float(np.sum(foreground**2)),1e-300)),effective_area_raw_m2=float(h*h*m.sum()**2/max(float(np.sum(m*m)),1e-300)),whole_roi_distribution=distribution(m),outer_roi_band_distribution=distribution(m[background]),prior=float(model.reg.evaluate(m)[0]),foreground_area_definition='Commonmax(m-.02,0)excess; noGTsupport',background_region='max(abs(x),abs(y))>=.085m; posthocdescriptiononly')
   r=dict(status='COMPLETED_DEV_MODEL_CORRECTION',case='FoamTwinDielTM',frequency_GHz=8.,side=args.side,arm=args.arm,initialization='uniform',
    train_relative_complex_residual=rel(d['train_mask']),withheld_all_relative_complex_residual=rel(d['withheld_mask']),withheld_sources_relative_complex_residual=rel(sh),withheld_receivers_relative_complex_residual=rel(rh),
    training_channels=int(d['train_mask'].sum()),withheld_channels=int(d['withheld_mask'].sum()),training_objective=float(best[1]),state_relative_residual=float(st['state_relative_residual']),
    cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,auxiliary_cpu_seconds=aux_cpu,auxiliary=aux['execution_receipt'],
    corrector=dict(message=str(rr.message),iterations=int(rr.nit),evaluations=int(rr.nfev),success=bool(rr.success)),physical_cost_training=dict(model.cost),physical_cost_prediction=dict(full.cost),min_contrast=float(best[0].min()),max_contrast=float(best[0].max()),
    posthoc_6GHz_metrics=posthoc6,physical_cost_prediction_6GHz=dict(full6.cost),material_diagnostics=descriptors,diagnostic_status='DEV_MODEL_CORRECTION_POST_MEASUREMENT',qualification='Onlymateriallowerbound0; freshinitial.1; bothfrequencyholdoutsalreadyexposed; noA+victoryclaim; Actual2005original measuredTM; frozenNIEpapergeometry/grid/frequency/DPCS/Fourier/beta ingredients butCIE/Acomparison notNIEpaper algorithm; unknownhorn/absoluteunits; no pixelGT',no_test_refitting=True)
   np.savez_compressed(out/'derived.npz',m=best[0],prediction=pred,prediction6=pred6,initial_material=m0,centers=d['centers'],train_mask=d['train_mask'],withheld_mask=d['withheld_mask'],auxiliary_endpoint=aux['last_joint_m'])
   dump(out/'RECEIPT.json',r);dump(out/'corrector_history.json',history);print(json.dumps({k:r[k] for k in ['case','arm','side','train_relative_complex_residual','withheld_all_relative_complex_residual','cpu_seconds']}),flush=True)
  except BaseException as exc:
   if 'best' in locals():np.savez_compressed(out/'paid_incumbent.npz',m=best[0],training_objective=best[1],cpu_seconds=time.process_time()-cpu)
   if 'history' in locals():dump(out/'corrector_history.json',history)
   dump(out/'FAILED.json',dict(error=repr(exc),status='INTERRUPTED' if isinstance(exc,KeyboardInterrupt) else 'FAILED',cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost),paid_material_preserved='best' in locals()));raise
if __name__=='__main__':main()
