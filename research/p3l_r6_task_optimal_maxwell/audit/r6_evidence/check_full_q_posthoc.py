"""Recompute saved-loss posthoc stronger Full-q comparison; zero physics/GT."""
from pathlib import Path
import json
import numpy as np
ROOT=Path.cwd();reported=json.loads((ROOT/'track_d/FULL_Q_REFERENCE_POSTHOC.json').read_text());rows={};errors=[]
for group in ('DLOW','DHIGH'):
 with np.load(ROOT/f'track_d/scores/{group}_PARENT_LOSSES.npz',allow_pickle=False) as z:parents=z['parent'].astype(str);selected=z['selected_deterministic'].copy();reference=z['full_kkt_deterministic'].copy()
 with np.load(ROOT/f'track_d/scores/{group}_BOOTSTRAP_DRAWS.npz',allow_pickle=False) as z:indices=z['parent_indices'].copy()
 assert selected.shape==reference.shape==(256,3) and indices.shape==(20000,256)
 s=selected.mean(1);f=reference.mean(1);c=s-1.01*f;draw=c[indices].mean(1);upper=float(np.quantile(draw,.99375,method='linear'));interval=np.quantile(draw,[.003125,.996875],method='linear').tolist()
 actual=dict(mean_selected=float(s.mean()),mean_full_q=float(f.mean()),relative_change=float(s.mean()/f.mean()-1),NI_linear_contrast=float(c.mean()),upper_at_original_pointwise_level=upper,descriptive_interval=interval,n=len(parents),qualification='posthoc only; pointwise original level does not preserve global8 coverage after adding comparisons',risk_NI_pointwise_indication=upper<0)
 for key in ('mean_selected','mean_full_q','relative_change','NI_linear_contrast','upper_at_original_pointwise_level','descriptive_interval'):
  if not np.allclose(actual[key],reported[group][key],atol=1e-15,rtol=1e-12):errors.append(group+':'+key)
 rows[group]=actual
qa=dict(status='PASS_POSTHOC_SAVED_LOSS_REPLAY' if not errors else 'AUDIT_ERROR',errors=errors,conditions=rows,new_Maxwell_RHS=0,GT_arrays_read=False,new_training=0,new_primary_comparison=False,global8_unchanged=True,conclusion='Original FullActual-reference cost_GO cannot be generalized to the stronger Fullq reference; DLOW 1pct NI remains HOLD there.')
(ROOT/'audit/r6_evidence/FULL_Q_POSTHOC_INDEPENDENT_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa,indent=2))
