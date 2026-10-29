"""Zero-RHS delivery audit and historical-reference comparison for REV1."""
import json,time
from pathlib import Path
from collections import Counter
import numpy as np
P=Path(__file__).resolve().parents[1];O=P/'results/exploration_rev1';A=P/'audit/exploration_rev1'
s=json.loads((O/'SUMMARY.json').read_text());old=json.loads((P/'results/exploration/SUMMARY.json').read_text());errors=[]
for r in s['rows']:
 if r['status']!='EXECUTED':continue
 path=O/'parents'/r['parent']/r['arm']/f"RMS_{r['radius']:.4f}"/'ENDPOINT.npz'
 with np.load(path,allow_pickle=False) as z:
  for k,shape in [('x',(66,)),('material',(576,)),('prediction',(32,8)),('current',(576,8)),('field',(576,8)),('J',(512,66))]:
   if z[k].shape!=shape or not np.isfinite(z[k]).all():errors.append(dict(path=str(path),array=k))
  if z['x'].dtype!=np.float64 or z['material'].dtype!=np.complex128:errors.append(dict(path=str(path),error='dtype'))
led=[json.loads(x) for x in (O/'RESOURCE_LEDGER.jsonl').read_text().splitlines()];beg={r['job'] for r in led if r['event']=='JOB_BEGIN'};end={r['job'] for r in led if r['event']=='JOB_END'}
qa=dict(status='PASS' if not errors and beg==end else 'FAIL',utc=time.time(),new_RHS=0,parents=8,endpoint_result_rows=len(s['rows']),full_endpoints_checked=sum(r['status']=='EXECUTED' for r in s['rows']),bad_arrays=errors,all_jobs_closed=beg==end,ledger_sequences_continuous=[r['after']['sequence'] for r in led]==list(range(1,len(led)+1)),frozen_REV1_source_equal=(P/'src/exploration_rev1.py').read_bytes()==(A/'FROZEN_REV1_SOURCE.py').read_bytes(),original_frozen_source_equal=(P/'src/exploration.py').read_bytes()==(P/'audit/exploration/FROZEN_EXPLORATION_SOURCE.py').read_bytes(),root_original_ledger_total=json.loads((P/'results/exploration/RESOURCE_SNAPSHOT.json').read_text())['RHS'],revision_total=s['resource']['total']['RHS'])
(A/'DELIVERY_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
ref=[]
for radius in [.01,.02]:
 for arm in ['A_NATIVE_FULLJ','F_FROZENJ','U_UNCORRECTED','B_RAY','C_RAY']:
  rows=[r for r in (s['rows'] if arm.endswith('_RAY') else old['rows']) if r['radius']==radius and r['arm']==arm];ok=[r for r in rows if r['status']=='EXECUTED'];ref.append(dict(radius=radius,arm=arm,qualified=sum(r.get('qualified',False) for r in ok),pair_le1=sum(r['whitened_pair_L2']<=1 for r in ok),median_pair=float(np.median([r['whitened_pair_L2'] for r in ok])),median_DeltaPhi=float(np.median([r['DeltaPhi'] for r in ok])),operational_RHS=sum(r['operational_cost']['RHS'] for r in rows),operational_wall=sum(r['operational_wall_seconds'] for r in rows),timing_window='REV1' if arm.endswith('_RAY') else 'ORIGINAL_E3'))
pair=[]
for radius in [.01,.02]:
 rs=[r for r in s['paired'] if r['radius']==radius];b=sum(r['B_operational_wall'] for r in rs);c=sum(r['C_operational_wall'] for r in rs)
 pair.append(dict(radius=radius,B_operational_wall=b,C_operational_wall=c,B_wall_relative_saving=1-b/c,B_RHS_relative_saving=1-sum(r['B_operational_RHS'] for r in rs)/sum(r['C_operational_RHS'] for r in rs),efficiency20pct_gate=1-b/c>=.2,distinct_endpoints=sum(r['latent_endpoint_difference']>1e-12 for r in rs),max_material_endpoint_difference_RMS=max(r['material_endpoint_difference_RMS'] for r in rs),B_minus_C_Phi=[r['B_minus_C_Phi'] for r in rs],C_Phi_better=sum(r['B_minus_C_Phi']>0 for r in rs)))
d=dict(new_RHS=0,reference_table=ref,paired_BC=pair,stop_counts=dict(Counter(r['arm']+':'+r['stop'] for r in s['stops'])),claim_boundary='Original E3 and REV1 timing windows differ; use within-REV1 paired B/C clocks only for speed comparison. All are exposed DEV; no fresh confirmation.',qa=qa)
(O/'REFERENCE_COMPARISON.json').write_text(json.dumps(d,indent=2)+'\n')
with (O/'REPORT.md').open('a') as f:
 f.write('\nReference endpoints (original timing window for A/Frozen/U):\n\n|arm|r|qualified|pair<=1|median pair|operational RHS|timing window|\n|---|---:|---:|---:|---:|---:|---|\n')
 for r in ref:f.write('|'+ '|'.join(str(r[k]) for k in ['arm','radius','qualified','pair_le1','median_pair','operational_RHS','timing_window'])+'|\n')
 f.write('\nWithin REV1 paired timing:\n')
 for r in pair:f.write(f'- r={r["radius"]}: fixed B relative wall saving={r["B_wall_relative_saving"]:.3%}; RHS saving={r["B_RHS_relative_saving"]:.3%}; 20% efficiency gate={r["efficiency20pct_gate"]}; endpoints differ in {r["distinct_endpoints"]}/8.\n')
 f.write('\nThis is a development mechanism repair. Original negative results remain part of the evidence; qualifying endpoints do not establish posterior or continuous-fiber validity.\n')
print(json.dumps(d))
