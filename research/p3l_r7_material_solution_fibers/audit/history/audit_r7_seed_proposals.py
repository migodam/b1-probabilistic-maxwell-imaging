from pathlib import Path
import json,numpy as np
out=Path(__file__).parent;root=Path.cwd();paths=[Path(p) for p in (out/'HISTORICAL_METADATA_PATHS.txt').read_text().splitlines() if 'p3l_r7_' not in p]
ids=set();seeds=set();sources=[];errors=[]
def visit(v,key='',underseed=False):
 seed=underseed or 'seed' in key.lower()
 if isinstance(v,dict):
  for k,x in v.items():visit(x,k,seed)
 elif isinstance(v,list):
  for x in v:visit(x,key,seed)
 elif isinstance(v,str):
  if 'parent' in key.lower() or key.lower() in ('id','ids'):ids.add(v)
 elif isinstance(v,int) and not isinstance(v,bool) and seed and 0<=v<2**128:seeds.add(v)
for p in paths:
 try:visit(json.loads(p.read_text()));sources.append(str(p))
 except Exception as e:errors.append({'path':str(p),'error':str(e)})
oldstreams={tuple(np.random.default_rng(s).bit_generator.random_raw(4).tolist()):s for s in seeds}
rows=[];seen={}
for i in range(8):
 streams={};collisions=[]
 for role,k in [('prior',1),('noise',2)]:
  ss=np.random.SeedSequence([2026102809,i,k]);finger=tuple(np.random.default_rng(ss).bit_generator.random_raw(4).tolist());streams[role]={'entropy':[2026102809,i,k],'seedsequence_uint64_marker':int(ss.generate_state(1,dtype=np.uint64)[0]),'first_four_raw64':list(finger)}
  if finger in oldstreams:collisions.append({'role':role,'historical_scalar_seed':oldstreams[finger]})
  if finger in seen:collisions.append({'role':role,'proposed_stream':seen[finger]})
  seen[finger]=[i,role]
 rows.append({'parent':f'R7_LOCKED_{i:03d}','id_overlap':f'R7_LOCKED_{i:03d}' in ids,'streams':streams,'stream_collisions':collisions})
# Separate E1 master domain; no truth/observations generated.
e1=np.random.default_rng(np.random.SeedSequence(2026102807)).bit_generator.random_raw(4).tolist()
r={'status':'PASS_METADATA_AND_INITIAL_STREAM_DISJOINT_WITH_SCOPE_LIMITS' if not errors and all(not r['id_overlap'] and not r['stream_collisions'] for r in rows) else 'REVIEW_REQUIRED','historical_metadata_files':len(sources),'historical_parent_ids':len(ids),'historical_scalar_seed_values':len(seeds),'sources':sources,'errors':errors,'rows':rows,'E1_master':2026102807,'E1_master_first4_raw64':e1,'limits':['No materials, observations, GT or Maxwell calls generated/read for R7.','First-four raw words detect reused initialization streams; they do not mathematically prove stochastic independence.','Nested entropy SeedSequence([...]) must be passed directly to default_rng; treating its uint64 marker as a new scalar seed gives another stream.','Historical batch RNG states and IDs are audited only as available in 60-name-filter metadata inventory; omitted manifests remain possible.','DEV4 intentionally reuses exposed historical parents, never counted among eight new locked parents.','E1 has distinct master seed; final code must preserve role/domain separation and record NumPy and bit generator versions.']}
(out/'R7_LOCKED_SEED_INDEPENDENCE_AUDIT.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:r[k] for k in ['status','historical_metadata_files','historical_parent_ids','historical_scalar_seed_values','errors']},indent=2))
