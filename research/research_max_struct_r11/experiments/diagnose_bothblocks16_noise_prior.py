"""OFFLINE_DIAGNOSTIC: finite actual-y noise-ball and exact inherited-prior audit, zero FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];BASE=ROOT/'delegated/diagnostics/bothblocks_confirmation';OUT=BASE/'finite_noise_prior_audit';sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Regularizer
METHODS=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks')
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
def prior_parts(reg,m):
 value=reg.evaluate(m)[0];d=reg.D@m
 return dict(R=float(value),R_smooth=float(reg.smooth*np.mean(d*d)),R_smoothed_TV=float(reg.tv*np.mean(np.sqrt(d*d+reg.epsilon**2)-reg.epsilon)))
def main():
 wall,cpu=time.perf_counter(),time.process_time()
 if not (BASE/'SCORING_RECEIPT.json').exists():raise RuntimeError('Only root-authorized postseal scoring artifacts may be audited')
 seal=json.loads((BASE/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text());scores=json.loads((BASE/'SINGLE_ARM_SCORES.json').read_text());scoreix={(x['case_id'],x['method']):x for x in scores};config=json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text());meta={x['case_id']:x for x in config['cases']}
 if len(seal['entries'])!=80 or len(meta)!=16:raise RuntimeError('All16x5 required')
 missing=[]
 for cid in meta:
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:
   if 'Y_clean' not in z.files:missing.append(str(cid)+':Y_clean')
 for e in seal['entries']:
  with np.load(e['prediction_path']) as z:
   if 'Y_prediction' not in z.files:missing.append(e['task']['task_id']+':Y_prediction')
 if missing:raise RuntimeError('Missing cached predictions; do not add FW: '+str(missing))
 OUT.mkdir(parents=True,exist_ok=True);reg=Regularizer((24,24));cached={};truthrows=[];rows=[];regcalls=0
 for cid in meta:
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:gt=z['m'].copy();clean=z['Y_clean'].copy()
  with np.load(ROOT/'data_public/bothblocks16'/f'{cid}.npz') as z:Y=z['Y'].copy();sigma=float(z['sigma'])
  if clean.shape!=Y.shape or Y.shape!=(48,12) or sigma<=0:raise ValueError('Inherited observation schema differs')
  scale=float(np.vdot(Y,Y).real);cn99=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2));noise=float(np.linalg.norm(clean-Y)/sigma);rt=prior_parts(reg,gt);regcalls+=1;dt=.5*float(np.vdot(clean-Y,clean-Y).real)/scale;pt=dt+rt['R']
  truthrow=dict(case_id=cid,family=meta[cid]['family'],regime=meta[cid]['regime'],high_chi=meta[cid]['high_chi'],sigma=sigma,Y_norm_squared=scale,CN99_threshold=cn99,
   truth_data_term=dt,truth_R=rt['R'],truth_R_smooth=rt['R_smooth'],truth_R_smoothed_TV=rt['R_smoothed_TV'],truth_Phi=pt,noise_Yclean_minus_Y_over_sigma=noise,GT_in_actual_y_CN99_ball=bool(noise<=cn99),
   truth_low_cells=int(np.count_nonzero(np.isclose(gt,2.5,atol=1e-12,rtol=0))),truth_high_cells=int(np.count_nonzero(np.isclose(gt,meta[cid]['high_chi'],atol=1e-12,rtol=0))),
   source='Cached sealed Y_clean; no new forward',status='OFFLINE_DIAGNOSTIC')
  truthrows.append(truthrow);cached[cid]=(gt,clean,Y,sigma,scale,cn99,truthrow)
 for e in seal['entries']:
  t=e['task'];cid=t['case_id'];method=t['method'];gt,clean,Y,sigma,scale,cn99,tr=cached[cid];primary=scoreix[cid,method]
  with np.load(e['prediction_path']) as z:m=z['m'].copy();yp=z['Y_prediction'].copy()
  rc=prior_parts(reg,m);regcalls+=1;res=yp-Y;dc=.5*float(np.vdot(res,res).real)/scale;pc=dc+rc['R'];cand_noise=float(np.linalg.norm(res)/sigma);clean_dist=float(np.linalg.norm(yp-clean)/sigma);dm=m-gt
  dd=dc-tr['truth_data_term'];dr=rc['R']-tr['truth_R'];dp=pc-tr['truth_Phi'];both=bool(tr['GT_in_actual_y_CN99_ball'] and cand_noise<=cn99)
  row=dict(case_id=cid,method=method,family=tr['family'],regime=tr['regime'],high_chi=tr['high_chi'],status='OFFLINE_DIAGNOSTIC_ONLY',
    original_material_success_unchanged=primary['material_success'],original_amplitude_NRMSE=primary['amplitude_NRMSE'],original_fixed_support_IoU=primary['fixed_support_IoU'],
    truth_data_term=tr['truth_data_term'],truth_R=tr['truth_R'],truth_R_smooth=tr['truth_R_smooth'],truth_R_smoothed_TV=tr['truth_R_smoothed_TV'],truth_Phi=tr['truth_Phi'],
    candidate_data_term=dc,candidate_R=rc['R'],candidate_R_smooth=rc['R_smooth'],candidate_R_smoothed_TV=rc['R_smoothed_TV'],candidate_Phi_reconstructed=pc,candidate_saved_Phi=e['cap_original_Phi'],Phi_reconstruction_absolute_error=abs(pc-e['cap_original_Phi']),
    delta_data_candidate_minus_truth=dd,delta_R_candidate_minus_truth=dr,delta_Phi_candidate_minus_truth=dp,candidate_Phi_lower_than_truth=bool(dp<0),candidate_prior_lower_than_truth=bool(dr<0),
    candidate_data_term_lower_than_truth=bool(dd<0),lower_Phi_despite_higher_data_term_via_prior_reduction=bool(dp<0 and dd>0 and dr<0),
    noise_Yclean_minus_Y_over_sigma=tr['noise_Yclean_minus_Y_over_sigma'],residual_Ypred_minus_Y_over_sigma=cand_noise,prediction_Ypred_minus_Yclean_over_sigma=clean_dist,CN99_threshold=cn99,
    GT_in_actual_y_CN99_ball=tr['GT_in_actual_y_CN99_ball'],candidate_in_actual_y_CN99_ball=bool(cand_noise<=cn99),actual_y_CN99_ball_contains_both=both,both_in_actual_y_ball_but_material_FAIL=bool(both and not primary['material_success']),
    material_L2_distance=float(np.linalg.norm(dm)),material_relative_L2_distance=float(np.linalg.norm(dm)/np.linalg.norm(gt)),material_RMSE=float(np.sqrt(np.mean(dm*dm))),
    material_amplitude_NRMSE_definition='RMSE/(maxGT-.1); copied primary gate unchanged',observed_prediction_distance_scope='Finite cached predictions at actual realized y; not an identifiability theorem')
  for name,mask in [('low',np.isclose(gt,2.5,atol=1e-12,rtol=0)),('high',np.isclose(gt,meta[cid]['high_chi'],atol=1e-12,rtol=0)),('background',np.isclose(gt,.1,atol=1e-12,rtol=0))]:
   n=int(mask.sum());row[name+'_truth_cells']=n;row[name+'_region_RMSE']=float(np.sqrt(np.mean(dm[mask]**2))) if n else None;row[name+'_region_mean_estimate']=float(m[mask].mean()) if n else None;row[name+'_region_median_estimate']=float(np.median(m[mask])) if n else None
   if name=='low':row['low_region_fraction_estimate_above_original_support1']=float(np.mean(m[mask]>1.)) if n else None
  rows.append(row)
 dump(OUT/'TRUTH_OBJECTIVES_ALL16.json',truthrows);csvwrite(OUT/'TRUTH_OBJECTIVES_ALL16.csv',truthrows);dump(OUT/'CANDIDATE_DATA_PRIOR_ALL80.json',rows);csvwrite(OUT/'CANDIDATE_DATA_PRIOR_ALL80.csv',rows)
 summary=[]
 for method in METHODS:
  rr=[x for x in rows if x['method']==method]
  summary.append(dict(method=method,objects=16,original_material_PASS=sum(x['original_material_success_unchanged'] for x in rr),truth_in_CN99_objects=sum(x['GT_in_actual_y_CN99_ball'] for x in rr),candidate_in_CN99=sum(x['candidate_in_actual_y_CN99_ball'] for x in rr),both_in_actual_y_CN99_ball=sum(x['actual_y_CN99_ball_contains_both'] for x in rr),both_in_ball_but_material_FAIL=sum(x['both_in_actual_y_ball_but_material_FAIL'] for x in rr),candidate_Phi_lower_than_truth=sum(x['candidate_Phi_lower_than_truth'] for x in rr),prior_lower_than_truth=sum(x['candidate_prior_lower_than_truth'] for x in rr),lower_Phi_despite_higher_data_via_prior=sum(x['lower_Phi_despite_higher_data_term_via_prior_reduction'] for x in rr),median_candidate_R_to_truth_R=float(np.median([x['candidate_R']/x['truth_R'] for x in rr])),median_clean_prediction_distance_over_sigma=float(np.median([x['prediction_Ypred_minus_Yclean_over_sigma'] for x in rr]))))
 dump(OUT/'SUMMARY_BY_METHOD.json',summary);csvwrite(OUT/'SUMMARY_BY_METHOD.csv',summary)
 objects=[]
 for tr in truthrows:
  rr=[r for r in rows if r['case_id']==tr['case_id']]
  objects.append(dict(case_id=tr['case_id'],family=tr['family'],regime=tr['regime'],truth_noise_over_sigma=tr['noise_Yclean_minus_Y_over_sigma'],truth_R=tr['truth_R'],truth_Phi=tr['truth_Phi'],GT_in_CN99=tr['GT_in_actual_y_CN99_ball'],both_in_ball_methods=[r['method'] for r in rr if r['actual_y_CN99_ball_contains_both']],both_in_ball_material_FAIL_methods=[r['method'] for r in rr if r['both_in_actual_y_ball_but_material_FAIL']],candidate_Phi_below_GT_methods=[r['method'] for r in rr if r['candidate_Phi_lower_than_truth']],outside_ball_methods=[r['method'] for r in rr if not r['candidate_in_actual_y_CN99_ball']]))
 dump(OUT/'OBJECT_PAIR_CLASSIFICATION_ALL16.json',objects);csvwrite(OUT/'OBJECT_PAIR_CLASSIFICATION_ALL16.csv',objects)
 lower=[r for r in rows if r['candidate_Phi_lower_than_truth']]
 aggregate=dict(objects=16,candidates=80,GT_in_actual_y_CN99_objects=sum(t['GT_in_actual_y_CN99_ball'] for t in truthrows),both_in_ball_pairs=sum(r['actual_y_CN99_ball_contains_both'] for r in rows),both_in_ball_material_FAIL_pairs=sum(r['both_in_actual_y_ball_but_material_FAIL'] for r in rows),candidate_outside_ball_with_truth_inside=sum(r['GT_in_actual_y_CN99_ball'] and not r['candidate_in_actual_y_CN99_ball'] for r in rows),lower_Phi_than_GT_pairs=len(lower),lower_data_and_lower_R_and_lower_Phi_pairs=sum(r['candidate_Phi_lower_than_truth'] and r['candidate_data_term_lower_than_truth'] and r['candidate_prior_lower_than_truth'] for r in rows),lower_Phi_despite_higher_data_term_pairs=sum(r['lower_Phi_despite_higher_data_term_via_prior_reduction'] for r in rows),median_prior_reduction_share_of_Phi_improvement=float(np.median([(-r['delta_R_candidate_minus_truth'])/(-r['delta_Phi_candidate_minus_truth']) for r in lower])),noise_truth_over_sigma_min=min(t['noise_Yclean_minus_Y_over_sigma'] for t in truthrows),noise_truth_over_sigma_max=max(t['noise_Yclean_minus_Y_over_sigma'] for t in truthrows),CN99_threshold=truthrows[0]['CN99_threshold'])
 dump(OUT/'AGGREGATE_FACTS.json',aggregate)
 receipt=dict(status='COMPLETE_OFFLINE_DIAGNOSTIC',objects=16,candidates=80,truth_clean_arrays_reused=16,candidate_prediction_arrays_reused=80,prior_evaluations=regcalls,
   prior=dict(source='research.research_z_alg_r4_diversity.src.physics.Regularizer',grid=[24,24],smooth=reg.smooth,tv=reg.tv,epsilon=reg.epsilon,edges=reg.edges,formula='smooth mean((Dm)^2)+tv mean(sqrt((Dm)^2+epsilon^2)-epsilon)'),
   Phi_data_formula='.5 ||prediction-Y||^2/||Y||^2',CN99_scope='Pointwise actual observed-y raw iid circular-complex noise ball; no universal or familywise identifiability claim',
   all80_selection_unchanged=True,new_forward=0,new_adjoint=0,new_optimization=0,missing_cached_Yclean_or_Ypred=missing,max_Phi_reconstruction_absolute_error=max(x['Phi_reconstruction_absolute_error'] for x in rows),
   CPU_seconds=time.process_time()-cpu,whole_child_CPU_seconds_near_exit=time.process_time(),wall_seconds=time.perf_counter()-wall,truth_access='Already root-authorized postseal GT; no online access')
 dump(OUT/'AUDIT_RECEIPT.json',receipt);print(json.dumps(dict(summary=summary,receipt=receipt)))
if __name__=='__main__':main()
