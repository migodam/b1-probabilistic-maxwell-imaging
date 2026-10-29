"""Post-selection analysis only. Imports no solver, spends zero RHS."""
import json,sys,time
from pathlib import Path
from collections import Counter
import numpy as np
P=Path(__file__).resolve().parents[1];OUT=P/'results/exploration';AUD=P/'audit/exploration'
def load(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,allow_nan=False,indent=2)+'\n')
def arr(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
s=load(OUT/'SUMMARY.json');rows=s['rows'];table=[];stops=Counter();cor=[]
for r in rows:
 if r['status']!='EXECUTED':continue
 path=OUT/'parents'/r['parent']/r['arm']/f"RMS_{r['radius']:.4f}"/'CORRECTION.json'
 if not path.exists():continue
 c=load(path);hs=c['history']
 for h in hs:
  if 'stop' in h:stops[(r['arm'],h['stop'])]+=1
 cor.append(dict(parent=r['parent'],arm=r['arm'],radius=r['radius'],executed_states=len(hs),selected_endpoint=c['selected_endpoint'],stop=next((h['stop'] for h in reversed(hs) if 'stop' in h),None),dynamic_geometry_arithmetic_seconds=sum(h.get('dynamic_geometry_arithmetic_seconds',0) for h in hs),substep_wall_seconds=sum(h.get('substep_wall_seconds',0) for h in hs),finite_weak_coordinate_defect_max=max((h.get('finite_weak_coordinate_defect') or 0 for h in hs),default=0.)))
for radius in [.01,.02]:
 for arm in sorted(set(r['arm'] for r in rows)):
  rs=[r for r in rows if r['radius']==radius and r['arm']==arm];valid=[r for r in rs if r['status']=='EXECUTED'];cs=[c for c in cor if c['radius']==radius and c['arm']==arm]
  table.append(dict(arm=arm,radius=radius,n_contexts=len(rs),available=len(valid),qualified=sum(r.get('qualified',False) for r in valid),operational_RHS=sum(r['operational_cost']['RHS'] for r in rs),diagnostic_RHS=sum(r['diagnostic_cost']['RHS'] for r in rs),research_total_RHS=sum(r['cost']['RHS'] for r in rs),operational_states=sum(r['operational_cost']['states'] for r in rs),diagnostic_states=sum(r['diagnostic_cost']['states'] for r in rs),operational_wall_seconds=sum(r['operational_wall_seconds'] for r in rs),diagnostic_wall_seconds=sum(r['diagnostic_wall_seconds'] for r in rs),dynamic_geometry_arithmetic_seconds=sum(c['dynamic_geometry_arithmetic_seconds'] for c in cs),substep_wall_seconds=sum(c['substep_wall_seconds'] for c in cs),one_state_only=sum(c['executed_states']==1 for c in cs),statuses=dict(Counter(r['status'] for r in rs))))
paired=[]
for radius in [.01,.02]:
 for parent in s['parents']:
  q={r['arm']:r for r in rows if r['radius']==radius and r['parent']==parent}
  b,c=q['B_FIXED_WEAK'],q['C_DYNAMIC_WEAK'];bp=arr(OUT/'parents'/parent/'B_FIXED_WEAK'/f'RMS_{radius:.4f}'/'ENDPOINT.npz');cp=arr(OUT/'parents'/parent/'C_DYNAMIC_WEAK'/f'RMS_{radius:.4f}'/'ENDPOINT.npz')
  paired.append(dict(parent=parent,radius=radius,B_qualified=b['qualified'],C_qualified=c['qualified'],same_endpoint_x=bool(np.array_equal(bp['x'],cp['x'])),endpoint_latent_L2=float(np.linalg.norm(bp['x']-cp['x'])),endpoint_material_RMS=float(np.linalg.norm(bp['material']-cp['material'])/24),B_operational_RHS=b['operational_cost']['RHS'],C_operational_RHS=c['operational_cost']['RHS'],B_operational_wall=b['operational_wall_seconds'],C_operational_wall=c['operational_wall_seconds'],B_minus_C_Phi=b['Phi']-c['Phi'],B_pair=b['whitened_pair_L2'],C_pair=c['whitened_pair_L2']))
# Historical A equivalence is an audit only, not independent evidence.
H7=P.parents[2]/'b1-p3l-r7-material-fibers/research/p3l_r7_material_solution_fibers';baseline=[]
for r in rows:
 if r['arm']!='A_NATIVE_FULLJ':continue
 old=load(H7/'E2/parents'/r['parent']/'B_FULLJ'/f"RMS_{r['radius']:.4f}"/'RESULT.json')
 baseline.append(dict(parent=r['parent'],radius=r['radius'],Phi_difference=r['Phi']-old['Phi'],material_RMS_difference=r['material_RMS']-old['material_RMS'],pair_difference=r['whitened_pair_L2']-old['whitened_pair_L2']))
# Source ledger reservation, callback, job wall and disk-I/O receipts.
led=[json.loads(x) for x in (OUT/'RESOURCE_LEDGER.jsonl').read_text().splitlines()];closed=[r for r in led if r['event']=='JOB_END'];actual_f=sum(r['counters'].get('actual_forward_RHS',0) for r in closed);actual_d=sum(r['counters'].get('actual_derivative_RHS',0) for r in closed)
ledger=dict(authoritative=led[-1]['after'],jobs=len(closed),all_jobs_closed=all(r['unused_forward_credit']==0 for r in closed),actual_forward_RHS=actual_f,actual_derivative_RHS=actual_d,callback_equals_reserved=actual_f+actual_d==led[-1]['after']['RHS'],controlled_wall_seconds=sum(r['wall_seconds'] for r in closed),CPU_seconds=sum(r['CPU_seconds'] for r in closed),physics_wall_seconds=sum(r['wall_seconds'] for r in led if r['event']=='PHYSICS_END'),peak_RSS_bytes=max(r['peak_RSS'] for r in closed)*(1 if sys.platform=='darwin' else 1024),array_IO={k:sum(r.get('array_IO',{}).get(k,0) for r in closed) for k in ['array_save_attempts','array_save_completed','array_save_failures','array_save_bytes','array_save_seconds']},solver_failures=[r for r in led if r['event']=='PHYSICS_END' and r['error'] is not None])
# All actual path nodes retain original likelihood/prior gates; finite-node status only.
pathsum=[dict(parent=r['parent'],arm=r['arm'],sampled_path_qualified=r['sampled_path_qualified'],node_count=len(r['nodes']),max_DeltaPhi=max(n['DeltaPhi'] for n in r['nodes']),max_material_projector_operator2_drift=max(n['weak_material_projector_operator2_drift'] for n in r['nodes'])) for r in s['paths']]
report=dict(status='COMPLETE_EXPOSED_DEV',new_RHS=0,table=table,paired_B_C=paired,corrections=cor,stop_counts=[dict(arm=a,stop=b,count=v) for (a,b),v in sorted(stops.items())],historical_native_equivalence=baseline,resource_audit=ledger,paths=pathsum,scope=['eight exposed historical DEV objects; no fresh confirmation','A uses historical radial retraction; B/C/S2 use fixed weak-coordinate decoder Newton retraction','B vs C same initial weak direction; dynamic refresh only current correction complement','a failed step or radius projection is retained, no extra optimizer rescue','SOM failure scoped to chosen dark-current lift, not all SOM null directions','Fourier/Wavelet mapped basis projections change target directions; projection errors disclosed','all current endpoint J diagnostic costs charged after method-selection clock','five-node actual correction path samples do not prove continuous fibers','no posterior or exact/continuous material-fiber claim'])
save(OUT/'E3_ANALYSIS.json',report)
lines=['# R9-S exploration E0/E1/E2/E3 evidence','', 'All eight R7 locked objects are now exposed development objects. New Maxwell work is complete; fresh confirmation has not been run.','', '| Arm | r | available | qualified | operational RHS | audit RHS | operational wall(s) | initial-state-only |','|---|---:|---:|---:|---:|---:|---:|---:|']
for r in table:lines.append('|'+ '|'.join(str(r[k]) if k!='operational_wall_seconds' else f'{r[k]:.4f}' for k in ['arm','radius','available','qualified','operational_RHS','diagnostic_RHS','operational_wall_seconds','one_state_only'])+'|')
lines+=['','B/C paired endpoint comparisons:']
for radius in [.01,.02]:
 rs=[r for r in paired if r['radius']==radius];bt=sum(r['B_operational_wall'] for r in rs);ct=sum(r['C_operational_wall'] for r in rs)
 lines.append(f'- r={radius}: exact same x in {sum(r["same_endpoint_x"] for r in rs)}/8; B/C operational RHS={sum(r["B_operational_RHS"] for r in rs)}/{sum(r["C_operational_RHS"] for r in rs)}; summed operational wall={bt:.4f}/{ct:.4f} s; B relative wall saving={1-bt/ct:.3%}; B/C qualified={sum(r["B_qualified"] for r in rs)}/{sum(r["C_qualified"] for r in rs)}.')
lines+=['','Corrector stops (all retained):']
for r in report['stop_counts']:lines.append(f'- {r["arm"]}: {r["stop"]} x {r["count"]}.')
lines+=['',f'Total actual new physics: {ledger["authoritative"]["RHS"]} complex RHS, {ledger["authoritative"]["states"]} states and LU factorizations, {ledger["authoritative"]["failures"]} solver failures. Reserved and callback counts match: {ledger["callback_equals_reserved"]}. Controlled wall={ledger["controlled_wall_seconds"]:.3f}s; CPU={ledger["CPU_seconds"]:.3f}s; physics wall={ledger["physics_wall_seconds"]:.3f}s; peak RSS={ledger["peak_RSS_bytes"]} bytes. All method, rejected-trial, audit-restoration and path work is charged.','',f'Native A reproduced R7 endpoint Phi to max absolute difference {max(abs(r["Phi_difference"]) for r in baseline):.3g}. This is a reproducibility check, not independent evidence.','', 'Limits and interpretation boundaries:']
lines += ['- '+t+'.' for t in report['scope']]
(OUT/'EXPLORATION_REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(ledger=ledger,paired=paired,table=table),allow_nan=False))
