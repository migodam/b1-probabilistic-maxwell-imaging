"""Read-only JSON/manifests and NPY headers, no Target imports, no GT arrays."""
from pathlib import Path
import json,collections,zipfile,ast,re,time
ROOT=Path.cwd();OUT=ROOT/'audit/r6_evidence';B1=ROOT.parents[2]
def walk(o,ids,seeds,key=''):
 if isinstance(o,dict):
  for k,v in o.items():walk(v,ids,seeds,k)
 elif isinstance(o,list):
  for v in o:walk(v,ids,seeds,key)
 elif isinstance(o,str) and (key in ('parent','parents','parent_id','parent_ids') or re.match(r'^(R[2456]|P2)',o)):ids.add(o)
 elif isinstance(o,int) and 'seed' in key.lower():seeds.add(o)
def extract(paths):
 ids=set();seeds=set();files=[];errors=[]
 for p in paths:
  try: d=json.loads(p.read_text());walk(d,ids,seeds);files.append(str(p))
  except Exception as e:errors.append(dict(path=str(p),error=repr(e)))
 return ids,seeds,files,errors
newfiles=list((ROOT/'configs').glob('parents_*.json'));ni,ns,nf,ne=extract(newfiles)
ns={int(r[k]) for p in newfiles for r in json.loads(p.read_text())['parents'] for k in ('prior_seed','noise_seed')}
history=[]
for project,sub in [('b1-p3l-r5-b-confirm','p3l_r5_b_independent_confirmation'),('b1-p3l-r4-conditional-transport','p3l_r4_conditional_transport'),('b1-p2-material-fiber','p2_material_fiber_20261019'),('b1-joint-p2','joint_p2')]:
 base=B1/project/'research'/sub
 candidates=list((base/'configs').glob('*.json'))
 candidates += [p for p in base.rglob('*MANIFEST*.json') if 'offline' not in str(p).lower() and p.stat().st_size<1000000]
 # Original P2 material-fiber stores all parent selection/config metadata under results/configs as well.
 ids,seeds,files,errors=extract(dict.fromkeys(candidates));history.append(dict(project=project,parent_ids_count=len(ids),seed_values_count=len(seeds),parent_overlap=sorted(ni&ids),seed_overlap=sorted(ns&seeds),files=files,parse_errors=errors))
new_group=[]
for p in newfiles:
 d=json.loads(p.read_text());ii,ss,_,_=extract([p]);ss={int(r[k]) for r in d['parents'] for k in ('prior_seed','noise_seed')};new_group.append(dict(group=d.get('group'),declared_n=d.get('n'),rows=len(d['parents']),unique_parents=len(ii),unique_seed_values=len(ss)))
# Read a consistent ledger prefix as at one saved snapshot, while campaign may continue.
snapshot=json.loads((ROOT/'RESOURCE_LEDGER_SNAPSHOT.json').read_text());limit=snapshot['sequence'];sums=collections.Counter();lastseq=0;activejobs=set();activeactions=set();errors=[];rows=0
with (ROOT/'RESOURCE_LEDGER.jsonl').open() as f:
 for line in f:
  try:r=json.loads(line)
  except Exception:break
  if r['sequence']>limit:break
  rows+=1
  if r['sequence']!=lastseq+1:errors.append('sequence_gap')
  lastseq=r['sequence'];sums['RHS']+=r['RHS'];sums['states']+=r['states'];sums['LU']+=r['LU'];sums['forward_RHS' if r['kind']=='forward' else 'derivative_RHS']+=r['RHS']
  if r['total_RHS']!=sums['RHS']:errors.append('row_running_total_mismatch')
  if r['event']=='JOB_BEGIN':activejobs.add(r['job'])
  if r['event']=='JOB_END':activejobs.remove(r['job'])
  if r['event']=='PHYSICS_BEGIN':activeactions.add(r['action'])
  if r['event']=='PHYSICS_END':activeactions.remove(r['action'])
for k in sums:
 if sums[k]!=snapshot[k]:errors.append('snapshot_mismatch_'+k)
if activejobs!=set(snapshot['active_jobs']):errors.append('active_jobs_mismatch')
if activeactions!=set(snapshot['active_actions']):errors.append('active_actions_mismatch')
# Online-only NPY headers: never open offline_truth or score arrays.
contexts=[]
for p in sorted((ROOT/'track_c/contexts/CDEV').glob('*/*/CONTEXT.npz')):
 with zipfile.ZipFile(p) as z:
  shapes={}
  for name in z.namelist():
   with z.open(name) as f:
    version=f.read(8);nbytes=2 if version[6]==1 else 4;length=int.from_bytes(f.read(nbytes),'little');header=ast.literal_eval(f.read(length).decode('latin1').strip());shapes[name[:-4]]=dict(shape=header['shape'],dtype=header['descr'])
 report=json.loads(p.with_name('ANCHOR_REPORT.json').read_text());contexts.append(dict(path=str(p),keys_and_shapes=shapes,report_GT_read=report['GT_read'],report_config=report['config']['name'],report_output_context_matches=report['output_context']==str(p),report_path_resource_cost=report['path_resource_cost']))
result=dict(status='READ_ONLY_PREFIX_AUDIT',new_groups=new_group,new_parent_count=len(ni),new_seed_count=len(ns),history_comparisons=history,new_internal_parent_duplicate_count=sum(x['rows'] for x in new_group)-len(ni),new_internal_seed_duplicate_count=2*sum(x['rows'] for x in new_group)-len(ns),ledger=dict(snapshot_sequence=limit,rows_checked=rows,counters=dict(sums),errors=errors,active_jobs_at_snapshot=list(activejobs),active_actions_at_snapshot=list(activeactions),qualification='Live prefix only; active paid work is not closed or failed by this audit'),online_context_headers=contexts,GT_arrays_opened=False,new_Maxwell_RHS=0)
(OUT/'ONLINE_LEDGER_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(groups=new_group,parent_count=len(ni),seed_count=len(ns),history=[{k:v for k,v in x.items() if k not in ('files','parse_errors')} for x in history],ledger=result['ledger'],contexts=len(contexts)),indent=2))
