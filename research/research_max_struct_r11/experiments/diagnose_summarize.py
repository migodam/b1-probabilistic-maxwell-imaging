from pathlib import Path
import json,csv
ROOT=Path(__file__).absolute().parents[1];OUT=ROOT/'delegated/diagnostics'
s=json.loads((OUT/'START_DIAGNOSTICS.json').read_text());r=json.loads((OUT/'FIXED_STATE_COUNTERFACTUALS.json').read_text())
summary=[]
for a in s:
 row=dict(state_id=a['state_id'],metric_const_norm_ratio=a['metric_const_norm_ratio'],cos_const_metric=a['cos_const_metric'],denominator_cond2=a['denominator_cond2'],
  raw_defect=a['defect']['raw_defect_fro'])
 for k in ('g_const','g_metric','g_full','g_lagged','g_ordinary'):
  for field in ('norm','box_projected_gradient_norm','descent_cos_negative_physical_gradient','projected_descent_cos_negative_physical_gradient'):row[k+'_'+field]=a['gradients'][k][field]
 for c in [z for z in r if z['state_id']==a['state_id']]:
  for k in ('delta_Phi','linear_prediction','step_norm','delta_NRMSE','moving_material_loss','frozen_material_loss'):row[c['method']+'_'+k]=c[k]
 summary.append(row)
with (OUT/'DIAGNOSTIC_SUMMARY.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
text=['# Fixed-state diagnostics: observed scope and counterexamples','',
'Ten existing G19 PGP endpoints (five cases, mf3 and mf12) were used. These are post-readout endpoints. Full R10 re-readout is therefore expected to be close to its box stationarity, and this is not a neutral comparison of full versus lagged trajectories. No inverse replay, new threshold, search or rescue budget was used. All three candidates were retained. GT scoring is OFFLINE_DIAGNOSTIC after candidate generation; it does not choose candidate execution or settings.','',
'At the four early annulus/shell/bridge/confirmation-annulus endpoints, the metric and constitutive gradients are strongly opposed: cosine ranges -0.978 to -0.945 and norm ratio ranges .771 to .958. This is evidence for cancellation in the moving-metric auxiliary gradient. It does not by itself establish that this cancellation causes every bad physical trajectory. The bridge is a direct limit: removing the metric term points against the physical descent direction at both saved states, and its finite lagged steps increase Phi.','',
'Full box projected-gradient norms range 6.35e-8 to 4.97e-5. Lagged norms range 1.71e-5 to .0433. Raw full gradients need not be small because many material coordinates sit at box bounds. Both raw and projected alignments are retained in JSON/CSV. Denominator condition numbers range 6.15 to 38.87; no singularity is needed for the observed cancellation.','',
'|State|metric/const|cos(const,metric)|lagged ΔPhi|lagged ΔNRMSE|ordinary ΔPhi|',
'|---|---:|---:|---:|---:|---:|']
for a in summary:text.append(f"|{a['state_id']}|{a['metric_const_norm_ratio']:.3f}|{a['cos_const_metric']:.3f}|{a['lagged_delta_Phi']:+.6g}|{a['lagged_delta_NRMSE']:+.6g}|{a['ordinary_delta_Phi']:+.6g}|")
text.extend(['','Lagged lowers Phi in 7/10 endpoints, raises it in both bridge endpoints and at the successful main-annulus mf12 endpoint. The confirmation coated object at mf3 lowers Phi but worsens NRMSE by .01944. All ten lagged steps lower raw defect norm; that alone does not guarantee a physical Phi decrease or material improvement. All lagged candidates lower their frozen-reference block objective while their moving-metric block objective rises; the objectives must not be conflated.','',
'An initial ordinary-control wiring mistake used has_reg=False and produced the scalar seed rather than the same-prior ordinary readout. It was corrected to has_reg=True on the same ten states, without repeating full/lagged or original gradients. The initial ten seed controls and their cost are archived as SEED_ONLY_ATTEMPT_* plus *_ordinary_seed_unregularized_CANDIDATE.npz. Final totals include them: 50 physical forwards / LU, 10 data adjoints, 720 complex RHS. The final 30 comparison candidates use the requested prior. Full/lagged use the same scalar seed and 8 matrix iterations; ordinary uses its inherited at-most-60 scalar prior iterations.','',
'Detailed evidence: STATE_INVENTORY.json, START_DIAGNOSTICS.json, FIXED_STATE_COUNTERFACTUALS.json, DIAGNOSTIC_SUMMARY.csv, DIAGNOSTIC_RECEIPT.json and ORDINARY_REPAIR_RECEIPT.json. Candidate arrays preserve every proposed material step. No same-state observed historical next material array exists.'])
(OUT/'DIAGNOSTIC_REVIEW.md').write_text('\n'.join(text)+'\n')
print('\n'.join(text))
