"""Post-score paired-noise reporting only; no tuning or extra PDE solves."""
import json,shutil,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
METHODS=['fullj_native','fullj_reciprocal','globaldual']
LABELS=['Full-J native','Full-J reciprocal','Gram dual + Full-J']
def main():
 if not (ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_NOISE_EXTENSION.json').exists():raise RuntimeError('All27 noisy predictions must freeze')
 noisy=json.loads((ROOT/'results/noise_extension_scores.json').read_text());clean={s:json.loads((ROOT/f'results/{s}_scores.json').read_text()) for s in ['holdout','budget_filled']};lookup={(s,r['observation_id'],r['method']):r for s,data in clean.items() for r in data['rows']};mapping=json.loads((ROOT/'truth_after_freeze/noise_extension/mapping_seedless.json').read_text());paired={r['observation_id']:r['parameters']['paired_source_observation_id'] for r in mapping['records']}
 rows=[]
 for r in noisy['rows']:
  original=paired[r['observation_id']];primary=lookup[('holdout',original,r['method'])];filled=lookup[('budget_filled',original,r['method'])]
  rows.append(dict(noise_observation_id=r['observation_id'],paired_clean_observation_id=original,digit=r['digit'],method=r['method'],clean_primary_RMSE=primary['RMSE'],clean_budget_filled_RMSE=filled['RMSE'],noisy_RMSE=r['RMSE'],noisy_minus_clean_budget_filled_RMSE=r['RMSE']-filled['RMSE'],clean_budget_filled_IoU=filled['IoU'],noisy_IoU=r['IoU'],clean_budget_filled_SSIM=filled['SSIM'],noisy_SSIM=r['SSIM']))
 out=dict(scope='paired 1e-3 noise extension on the SAME nine target materials; not18independent material targets',method_source_and_configs_frozen_before_main_GT=True,main_GT_already_known_during_extension=True,noise_predictions_frozen_before_extension_scoring=True,no_new_PDE_for_noise_generation=True,config_difference='noisy high final-GN iteration cap implements the90s budget; clean primary1000cap and budget-filled supplement retained separately; trajectories differ',rows=rows,aggregate_noisy_means=noisy['aggregate_means'])
 (ROOT/'results/noise_paired_comparison.json').write_text(json.dumps(out,indent=2)+'\n');shutil.copy2(ROOT/'figures/noise_extension_all_digits_main.png',ROOT/'figures/digits_noisy.png')
 fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
 for method,label in zip(METHODS,LABELS):
  mr=[r for r in rows if r['method']==method];ax.scatter([r['clean_budget_filled_RMSE'] for r in mr],[r['noisy_RMSE'] for r in mr],label=label,s=35,alpha=.85)
  for r in mr:ax.annotate(str(r['digit']),(r['clean_budget_filled_RMSE'],r['noisy_RMSE']),xytext=(3,3),textcoords='offset points',fontsize=7)
 maximum=max(max(r['clean_budget_filled_RMSE'],r['noisy_RMSE']) for r in rows)*1.05;ax.plot([0,maximum],[0,maximum],color='gray',linestyle='--',linewidth=1);ax.set_xlim(0,maximum);ax.set_ylim(0,maximum);ax.set_xlabel('Same target, noiseless budget-filled RMSE');ax.set_ylabel('Same target, relative-noise1e-3 RMSE');ax.set_title('Paired same-nine noise extension; digit labels');ax.legend(fontsize=8);ax.grid(alpha=.2);fig.savefig(ROOT/'figures/noise_paired_error.png',dpi=170);plt.close(fig)
 # Full byte comparison against the ten-file pre-generation/pre-main-GT freeze.
 freeze=json.loads((ROOT/'audit/NOISE_EXTENSION_SOURCE_FREEZE.json').read_text());checks=[]
 for f in freeze['files']:
  saved=f.get('snapshot')
  if saved is None:saved=str(Path(freeze.get('source_snapshot','audit/NOISE_EXTENSION_SOURCE_SNAPSHOT'))/f['path'])
  a=(ROOT/f['path']).read_bytes();b=(ROOT/saved).read_bytes();checks.append(dict(path=f['path'],byte_equal=a==b,current_bytes=len(a),snapshot_bytes=len(b)))
 receipt=dict(checked_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),method='exact full file bytes, no hashes',all_noise_source_files_unchanged=all(c['byte_equal'] for c in checks),file_count=len(checks),files=checks)
 (ROOT/'audit/EXPERIMENT_NOISE_SOURCE_EQUALITY.json').write_text(json.dumps(receipt,indent=2)+'\n')
 if not receipt['all_noise_source_files_unchanged']:raise RuntimeError('Noise source drift')
 print(json.dumps(dict(aggregates=noisy['aggregate_means'],source_files_unchanged=True,paired_targets=9,independent_material_targets=9)))
if __name__=='__main__':main()
