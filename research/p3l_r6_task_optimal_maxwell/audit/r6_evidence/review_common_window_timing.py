"""Existing timing arrays/contexts only; zero new physics and zero GT reads."""
from pathlib import Path
import json,hashlib,collections
import numpy as np
ROOT=Path.cwd();base=ROOT/'track_c/common_window_timing';r=json.loads((base/'RESULTS.json').read_text());protocol=json.loads((base/'PROTOCOL_PRE_TIMING.json').read_text());errors=[];differences=[];closed_sum=0
seen=collections.Counter((x['group'],x['parent'],x['config']) for x in r['records'])
for group,ids in protocol['parents'].items():
 for parent in ids:
  for config in protocol['configs']:
   if seen[(group,parent,config)]!=1:errors.append('missing_or_duplicate_path')
for row in r['records']:
 original=Path(row['original_context_identity']['path']);new=base/'contexts'/row['group']/row['parent']/row['config']/'CONTEXT.npz'
 with np.load(original,allow_pickle=False) as a,np.load(new,allow_pickle=False) as b:
  diff={k:float(np.max(np.abs(a[k]-b[k]))) for k in ('q0','chi0','prediction')}
 if any(v!=0. for v in diff.values()):errors.append('endpoint_changed')
 if diff!=row['endpoint_difference']:errors.append('reported_difference_mismatch')
 differences.append(dict(group=row['group'],parent=row['parent'],config=row['config'],difference=diff))
 closed=json.loads((ROOT/'raw/jobs'/f'TIMING_{row["parent"]}_{row["config"]}'/'CLOSED.json').read_text());counts=closed['actual_counters'];rhs=counts['actual_forward_RHS']+counts['actual_derivative_RHS'];closed_sum+=rhs
 if rhs!=row['RHS'] or closed['status']!='COMPLETE':errors.append('timing_job_cost_or_status')
recomputed={}
for group in protocol['parents']:
 with np.load(base/f'{group}_TIMING_DRAWS.npz',allow_pickle=False) as z:wall=z['wall_seconds'].copy();ix=z['bootstrap_indices'].copy();draw=z['bootstrap_wall_means'].copy();configs=z['configs'].astype(str).tolist()
 if ix.shape!=(20000,16) or not np.array_equal(wall[ix].mean(1),draw):errors.append('bootstrap_draw_replay')
 result={}
 for a,b in [(1,0),(2,0),(2,1)]:
  key=f'{configs[a]}_vs_{configs[b]}';vals=1-draw[:,a]/draw[:,b];out=dict(saving=1-float(wall[:,a].mean()/wall[:,b].mean()),descriptive95_interval=np.quantile(vals,[.025,.975]).tolist());result[key]=out
  old=r['conditions'][group]['comparisons'][key]
  if out!=old:errors.append('summary_replay')
 recomputed[group]=result
if closed_sum!=268544 or closed_sum!=r['actual_new_RHS'] or closed_sum>protocol['max_RHS']:errors.append('RHS_total')
qa=dict(status='PASS_READ_ONLY_TIMING_REPLAY_QA' if not errors else 'AUDIT_ERROR',errors=errors,paths=len(r['records']),parents=32,actual_paid_RHS=closed_sum,max_endpoint_difference=max(v for x in differences for v in x['difference'].values()),recomputed_comparisons=recomputed,endpoint_checks=differences,new_physics_RHS=0,GT_read=False,scientific_limit='post-score descriptive timing control, not new independent risk/selection; original global8 unchanged')
(ROOT/'audit/r6_evidence/TIMING_DRIFT_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps({k:v for k,v in qa.items() if k!='endpoint_checks'},indent=2))
