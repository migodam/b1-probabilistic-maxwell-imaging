"""Private postfreeze scoring process; no feedback to solver implementations."""
import json,time,argparse
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[2]
def main():
 a=argparse.ArgumentParser();a.add_argument('--freeze-receipt',required=True);args=a.parse_args();receipt=Path(args.freeze_receipt)
 if not receipt.exists():raise RuntimeError('Need implementation+prediction freeze receipt')
 gt=np.load(P/'private/highfreq48/sealed_targets.npz');rows=[]
 for f in sorted((P/'results').glob('highfreq48_HF48*.json')):
  r=json.loads(f.read_text());m=np.asarray(r.get('m',[]));truth=gt[r['case']];a=r['accounting'];v=r.get('validation',{});rms=float(np.linalg.norm(m-truth)/np.sqrt(48)) if m.shape==(48,) else None
  rows.append({'file':f.name,'case':r['case'],'family':r['family'],'noise_mode':r['noise_mode'],'method':r['method'],'status':r['status'],'material_RMS':rms,'material_MAE':float(np.mean(abs(m-truth))) if rms is not None else None,'material_max_absolute_error':float(np.max(abs(m-truth))) if rms is not None else None,'validation':v,'clean_imaging_pass':bool(rms is not None and rms<=1e-5 and v.get('response_fit_pass')) if r['noise_mode']=='clean' else None,'incremental_wall_seconds':r['incremental_wall_seconds'],'incremental_CPU_seconds':r['incremental_CPU_seconds'],'complex_rhs':a['total_complex_rhs'],'full_forward_calls':a['forward']['calls'],'exact_J_calls':a['exact_J']['calls'],'validation_calls':a['validation']['calls'],'nfev':r.get('nfev')})
 setup=json.loads((P/'results/highfreq48_setup.json').read_text());green=sum(e['wall_seconds'] for e in setup['events'] if e['kind']=='GREEN_SETUP');svd=sum(e['wall_seconds'] for e in setup['events'] if e['kind']=='RECEIVER_SVD');totals={}
 for method in ['FullJ_midpoint','exactSVD','SVD_FullJ']:
  selected=[r for r in rows if r['method']==method];need_svd=method!='FullJ_midpoint';incremental=sum(r['incremental_wall_seconds'] for r in selected);setup_cost=green+(svd if need_svd else 0)
  totals[method]={'cases':len(selected),'shared_setup_wall_seconds':setup_cost,'incremental_wall_seconds':incremental,'total_wall_seconds_setup_once':setup_cost+incremental,'CPU_incremental_seconds':sum(r['incremental_CPU_seconds'] for r in selected),'complex_rhs':sum(r['complex_rhs'] for r in selected),'clean_imaging_passes':sum(r['clean_imaging_pass'] is True for r in selected),'clean_material_RMS':[r['material_RMS'] for r in selected if r['noise_mode']=='clean'],'noisy_material_RMS':[r['material_RMS'] for r in selected if r['noise_mode']!='clean']}
 output={'status':'SCORED_AFTER_FROZEN_18RUNS','scored_unix':time.time(),'freeze_receipt':str(receipt),'GT_values_disclosed':False,'rows':rows,'method_totals':totals,'setup_components':{'Green_wall_seconds':green,'receiver_SVD_wall_seconds':svd},'noise_scope':'fixed1e-3 relative observation perturbation, same data for all methods; no postscore tuning; raw exactSVD may be illegal','model_scope':'NEW48-cell highfreq pulse model, not Full66/oldlowfreq applicability','crossover_policy':'same Green setup; exactSVD/hybrid add receiverSVD; per-case incremental and setup-once totals reported, no stable timing claim from tiny samples'}
 (P/'results/highfreq48_scores.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(totals,indent=2))
if __name__=='__main__':main()
