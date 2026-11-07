"""POSTHOC_OFFLINE_DIAGNOSTIC first retained Phi-incumbent CN99 hit; no new FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,argparse,csv
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];BASE=ROOT/'delegated/diagnostics/bothblocks_confirmation';OUT=BASE/'first_cn99_hit';sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Regularizer
METHODS=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks')
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def select():
 wall,cpu=time.perf_counter(),time.process_time();OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'FIRST_HIT_SELECTION_PUBLIC_ONLY.json').exists():raise FileExistsError('First-hit mapping already sealed; refuse GT-dependent reselection')
 seal=json.loads((BASE/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text())
 if len(seal['entries'])!=80:raise RuntimeError('All80 sealed outcomes required')
 reg=Regularizer((24,24));rows=[];materials=[];negative=0;maxround=0.
 for e in seal['entries']:
  t=e['task'];d=REPO/t['output_dir'];path=d/'R11_INCUMBENT_TRAJECTORY.npz';pub={}
  # Only the observation file is needed; no GT/scored-table read in this phase.
  for rel in t['public_files']:
   if Path(rel).name=='operators.npz':continue
   with np.load(REPO/rel) as z:pub.update({k:z[k].copy() for k in ('Y','sigma') if k in z.files})
  Y=pub['Y'];sigma=float(pub['sigma']);scale=float(np.vdot(Y,Y).real);threshold=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2))
  with np.load(path) as z:m=z['m'].copy();phi=z['original_Phi'].copy();cp=z['CPU_seconds'].copy();lu=z['LU'].copy();rhs=z['complex_RHS'].copy()
  Dm=reg.D@m.T;prior=reg.smooth*np.mean(Dm*Dm,axis=0)+reg.tv*np.mean(np.sqrt(Dm*Dm+reg.epsilon**2)-reg.epsilon,axis=0);data=phi-prior
  negative+=int(np.count_nonzero(data<0));maxround=max(maxround,float(max(0.,-data.min())))
  if np.any(data < -1e-12):raise RuntimeError('Saved Phi-prior has a material negative data term; no new FW attempted')
  distance=np.sqrt(np.maximum(data,0.)*2*scale)/sigma
  hit=np.flatnonzero((cp<=480)&(distance<=threshold));i=int(hit[0]) if len(hit) else None
  with np.load(e['prediction_path']) as z:finalm=z['m'].copy();finalyp=z['Y_prediction'].copy()
  finalres=float(np.linalg.norm(finalyp-Y)/sigma);fm=d/'R11_CAP_PREDICTION.json';finalphi=float(json.loads(fm.read_text())['objective']);firstphase=None
  if i is not None:
   paid=json.loads((d/'R11_PAID_PHYSICAL_STATES.json').read_text());match=next((p for p in paid if p['CPU_seconds']==float(cp[i])),None);firstphase=match['phase'] if match else 'UNKNOWN_SAVED_PHASE'
  chosen=m[i] if i is not None else finalm
  materials.append(chosen.copy())
  rows.append(dict(case_id=t['case_id'],method=t['method'],classification='POSTHOC_OFFLINE_DIAGNOSTIC',selection_GT_access='NONE_THIS_PHASE; GT already exposed historically but never read for first-hit choice',rule='First saved CPU<=480 public-Phi incumbent with sqrt(2*(Phi-R)*||Y||^2)/sigma <= fixed CN99 threshold; no threshold scan',
    status='FIRST_RETAINED_INCUMBENT_CN99_HIT' if i is not None else 'NO_HIT_ORIGINAL_ENDPOINT_FAILURE_RETAINED',source=str(path),source_index=i,material_index=len(materials)-1,CN99_threshold=threshold,first_hit_phase=firstphase,
    first_hit_paid_CPU_seconds=float(cp[i]) if i is not None else None,first_hit_LU=int(lu[i]) if i is not None else None,first_hit_complex_RHS=int(rhs[i]) if i is not None else None,
    first_hit_original_Phi=float(phi[i]) if i is not None else None,first_hit_prior_R=float(prior[i]) if i is not None else None,first_hit_reconstructed_data_term=float(data[i]) if i is not None else None,first_hit_reconstructed_residual_over_sigma=float(distance[i]) if i is not None else None,
    final_original_Phi=finalphi,final_residual_over_sigma=finalres,final_CN99_pass=bool(finalres<=threshold),final_actual_complete_child_CPU_seconds=e['terminal_job']['child_process_CPU_seconds'],
    diagnostic_prefix_CPU_proxy_seconds=float(cp[i]) if i is not None else e['terminal_job']['child_process_CPU_seconds'],CPU_scope='Observed paid-state CPU if hit; exact actual full-child CPU if no hit. Not a measured rerun stopped at the hit; serialization/atomic-stop tail unknown.',
    first_scope='First retained Phi-incumbent, not first of every unretained solver iterate'))
 np.savez_compressed(OUT/'FIRST_HIT_SEALED_MATERIALS.npz',m=np.asarray(materials))
 dump(OUT/'FIRST_HIT_SELECTION_PUBLIC_ONLY.json',dict(status='SEALED_PUBLIC_SIGMA_CN99_MAPPING_BEFORE_THIS_DIAGNOSTIC_GT_READ',classification='POSTHOC_OFFLINE_DIAGNOSTIC_NOT_NEW_CONFIRMATION_ALGORITHM',rows=rows,selection_fields='Only saved m/Phi/CPU/LU/RHS and publicY/sigma; originalR; no GT or material score input',negative_rounded_data_terms=negative,max_negative_roundoff=maxround,new_forward=0,new_optimizer=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
 print(json.dumps(dict(status='SEALED_PUBLIC_ONLY',objects=16,methods=5,hits=sum(r['source_index'] is not None for r in rows))))
def score():
 wall,cpu=time.perf_counter(),time.process_time();mapping=json.loads((OUT/'FIRST_HIT_SELECTION_PUBLIC_ONLY.json').read_text())
 if mapping['status']!='SEALED_PUBLIC_SIGMA_CN99_MAPPING_BEFORE_THIS_DIAGNOSTIC_GT_READ':raise RuntimeError('Mapping must be sealed before GT read')
 from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
 cfg=json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text());meta={r['case_id']:r for r in cfg['cases']};final=json.loads((BASE/'SINGLE_ARM_SCORES.json').read_text());fi={(r['case_id'],r['method']):r for r in final};truth={}
 for cid in meta:
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:truth[cid]=z['m'].copy()
 with np.load(OUT/'FIRST_HIT_SEALED_MATERIALS.npz') as z:materials=z['m'].copy()
 rows=[]
 for r in mapping['rows']:
  cid=r['case_id'];gt=truth[cid];f=fi[cid,r['method']];mt=material_metrics(materials[r['material_index']],gt,float(gt.max()-.1),meta[cid]['high_chi']);x=dict(r,family=meta[cid]['family'],regime=meta[cid]['regime'],original_final_material_success=f['material_success'],original_final_amplitude_NRMSE=f['amplitude_NRMSE'],original_final_fixed_support_IoU=f['fixed_support_IoU'],selected_material_success=mt['material_success'],selected_amplitude_NRMSE=mt['amplitude_NRMSE'],selected_fixed_support_IoU=mt['fixed_support_IoU'],delta_amplitude_NRMSE_selected_minus_final=mt['amplitude_NRMSE']-f['amplitude_NRMSE'],delta_fixed_support_IoU_selected_minus_final=mt['fixed_support_IoU']-f['fixed_support_IoU'],GT_access='After public-only first-hit seal; offline posthoc diagnostic only');rows.append(x)
 dump(OUT/'FIRST_HIT_VS_FINAL_ALL80.json',rows);csvwrite(OUT/'FIRST_HIT_VS_FINAL_ALL80.csv',rows)
 summary=[]
 for method in METHODS:
  rr=[r for r in rows if r['method']==method];hit=[r for r in rr if r['source_index'] is not None]
  summary.append(dict(method=method,objects=16,first_hit_count=len(hit),no_hit_original_failures_retained=16-len(hit),final_material_PASS=sum(r['original_final_material_success'] for r in rr),first_hit_or_original_failure_material_PASS=sum(r['selected_material_success'] for r in rr),mean_selected_amplitude_NRMSE=float(np.mean([r['selected_amplitude_NRMSE'] for r in rr])),mean_final_amplitude_NRMSE=float(np.mean([r['original_final_amplitude_NRMSE'] for r in rr])),mean_selected_minus_final_NRMSE=float(np.mean([r['delta_amplitude_NRMSE_selected_minus_final'] for r in rr])),earlier_material_error_smaller=sum(r['delta_amplitude_NRMSE_selected_minus_final']<0 for r in hit),earlier_material_error_larger=sum(r['delta_amplitude_NRMSE_selected_minus_final']>0 for r in hit),material_gate_improvements=sum(r['selected_material_success'] and not r['original_final_material_success'] for r in rr),material_gate_regressions=sum(r['original_final_material_success'] and not r['selected_material_success'] for r in rr),sum_diagnostic_prefix_CPU_proxy_seconds=sum(r['diagnostic_prefix_CPU_proxy_seconds'] for r in rr),sum_actual_final_child_CPU_seconds=sum(r['final_actual_complete_child_CPU_seconds'] for r in rr),CPU_scope='Paid-state prefix proxy, not exact executed early-stop CPU; no method rerun'))
 dump(OUT/'SUMMARY_BY_METHOD.json',summary);csvwrite(OUT/'SUMMARY_BY_METHOD.csv',summary)
 dump(OUT/'SCORE_RECEIPT.json',dict(status='COMPLETE_POSTHOC_OFFLINE_DIAGNOSTIC',objects=16,candidates=80,new_forward=0,new_adjoint=0,new_optimization=0,selection_repeated_after_GT=False,original_confirmation_results_changed=False,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
 print(json.dumps(summary))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--phase',choices=('select','score'),required=True);a=p.parse_args();select() if a.phase=='select' else score()
