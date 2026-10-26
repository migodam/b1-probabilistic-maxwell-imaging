"""Read-only metadata audit. No truth payloads, inference, training or physics."""
from pathlib import Path
import json, zipfile, statistics, collections, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
B1=ROOT.parent
R4=B1/'b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport'
R2=B1/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
R5=B1/'b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab'
def read(p):return json.loads(p.read_text())
def write(name,d):(OUT/name).write_text(json.dumps(d,indent=2)+'\n')
def ident(p):
 s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def header(p):
 result={}
 with zipfile.ZipFile(p) as z:
  for info in z.infolist():
   if not info.filename.endswith('.npy'):continue
   with z.open(info) as f:
    v=np.lib.format.read_magic(f)
    shape,order,dtype=(np.lib.format.read_array_header_1_0(f) if v==(1,0) else np.lib.format.read_array_header_2_0(f))
   result[info.filename[:-4]]=dict(shape=list(shape),dtype=str(dtype),array_bytes=info.file_size,compressed_bytes=info.compress_size)
 return result
cohorts=[];historical_ids=set();historical_seeds=[]
def cohort(name,rows,source,seedrule):
 ids=[r['parent'] for r in rows];historical_ids.update(ids)
 cohorts.append(dict(name=name,count=len(ids),unique_ids=len(set(ids)),first=ids[0] if ids else None,last=ids[-1] if ids else None,source=str(source),seed_rule=seedrule,parents=ids))
p=ROOT/'research/p1_material_target_20261018/results/MATERIAL_MANIFEST.json'
rr=read(p)['rows']
for group in sorted({r['group'] for r in rr}):
 rows=[dict(parent=r['parent_id']) for r in rr if r['group']==group];cohort('P1_'+group,rows,p,'batched default_rng(seed); draw normal(n,64), then uniform(n,2), not per-parent seed')
historical_seeds.extend([2026101801,2026101802,2026101803,2026101804])
for name,p,key in [('P2',ROOT/'research/p2_material_fiber_20261019/configs/PARENTS.json','split'),('P3L',ROOT/'research/p3l_local_posterior_imaging/configs/PARENTS.json','group')]:
 rows=read(p)
 for group in sorted({r[key] for r in rows}):
  rs=[r for r in rows if r[key]==group];cohort(name+'_'+group,rs,p,'SeedSequence(row.seed_derivation+[0/1])' if name=='P2' else 'SeedSequence(row.seed+[0/1])')
  for r in rs:
   seed=r.get('seed_derivation',r.get('seed'));historical_seeds.extend([seed+[0],seed+[1]])
for p in sorted(R2.glob('*/jobs/data_*/DATA_MANIFEST.json')):
 d=read(p);rows=[dict(parent=v) for v in d['parents']];cohort('R2_'+p.parent.name,rows,p,'SeedSequence([2026102202,split_code,index,0/1])')
 for parent in d['parents']:
  k=int(parent.rsplit('_',1)[1]);historical_seeds.extend([d['seed']+[k,0],d['seed']+[k,1]])
for group in ['dev','locked','confirmation']:
 rows=[]
 for p in sorted((R4/'inputs/online').glob('*.json')):
  d=read(p)
  if d.get('split')==group:
   rows.append(d);historical_seeds.extend([d['seed']+[0],d['seed']+[1]])
 cohort('R4_'+group,rows,R4/'inputs/online','SeedSequence([cohort_seed,index,0/1])')
def seedstate(seed):return tuple(int(v) for v in np.random.SeedSequence(seed).generate_state(4))
historical_states={seedstate(s) for s in historical_seeds}
pending={};allids=[];alllatent=[];allnoise=[];allmoment=[]
for arm in ['A','B']:
 p=R5/f'shared/pending_confirmation_manifests/confirmation_{arm}_256.json';d=read(p);rows=d['parents'];ids=[r['parent'] for r in rows];lat=[r['latent_seed'] for r in rows];noise=[r['noise_seed'] for r in rows]
 generated=[]
 for base in [R5,ROOT/'research/p3l_r5_b_independent_confirmation']:
  for f in base.rglob('*.npz'):
   if any(parent in f.name for parent in ids):generated.append(str(f))
 moment=[int(np.random.SeedSequence([20261025121,*v.encode('utf8')]).generate_state(1)[0]) for v in ids]
 pending[arm]=dict(source=ident(p),status=d['status'],master_seed=d['master_seed'],count=len(ids),duplicate_ids=len(ids)-len(set(ids)),duplicate_latent_seeds=len(lat)-len(set(lat)),duplicate_noise_seeds=len(noise)-len(set(noise)),latent_noise_overlap=sorted(set(lat)&set(noise)),historical_parent_overlap=sorted(set(ids)&historical_ids),historical_SeedSequence_state128_overlap=sum(seedstate(v) in historical_states for v in lat+noise),generated_npz_filename_matches=generated,K32_parent_seed_collisions=len(moment)-len(set(moment)),parents=rows)
 allids+=ids;alllatent+=lat;allnoise+=noise;allmoment+=moment
write('PARENT_SEED_AUDIT.json',dict(status='PASS_METADATA_DISJOINT_WITH_SCOPE_LIMITS',cohorts=cohorts,pending=pending,historical_unique_parent_ids=len(historical_ids),historical_seed_states128=len(historical_states),AB_parent_overlap=len(set(pending['A']['parents'][i]['parent'] for i in range(256))&set(pending['B']['parents'][i]['parent'] for i in range(256))),AB_all_latent_noise_duplicate_count=len(alllatent+allnoise)-len(set(alllatent+allnoise)),AB_K32_seed_duplicate_count=len(allmoment)-len(set(allmoment)),new_physics_RHS=0,new_training=0,truth_payloads_read=False,hashes=False,qualification='Named historical P1/P2/P3L/R2/R4 metadata only; no full content duplicate proof; pending seed derivation source missing; SeedSequence state128 is RNG metadata, not file hashing'))
assets=[]
for name,a in read(R4/'AUDIT/risk/ASSET_INVENTORY.json')['assets'].items():
 p=Path(a['path']);r=dict(name=name,path=str(p),exists=p.exists())
 if p.is_file():
  r.update(ident(p),headers=header(p))
  if 'parent' in r['headers']:
   with np.load(p,allow_pickle=False) as z:ids=z['parent'].astype(str).tolist()
   r['parent_metadata']=dict(count=len(ids),unique=len(set(ids)),first=ids[0] if ids else None,last=ids[-1] if ids else None,outside_historical_manifest_ids=sorted(set(ids)-historical_ids),pending_overlap=sorted(set(ids)&set(allids)),read_scope='parent string field only; no observation/chi/q payload')
 assets.append(r)
checkpoints=[]
for base in [R4/'models',R5/'stage_b_material_risk']:
 for p in sorted(base.rglob('*.pt')):checkpoints.append(ident(p))
write('ASSET_METADATA_INVENTORY.json',dict(arrays=assets,checkpoints=checkpoints,checkpoint_count=len(checkpoints),checkpoint_bytes=sum(x['bytes'] for x in checkpoints),read_scope='stat, ZIP directory and NPY headers only; checkpoint weights not loaded',truth_payloads_read=False,hashes=False))
cost=[]
for group in ['train','validation','dev','locked','confirmation']:
 rs=[]
 for p in sorted((R4/'results/contexts'/group).glob('*.json')):
  d=read(p)
  if 'new_RHS' in d:rs.append(d)
 vals=[r['new_RHS'] for r in rs]
 cost.append(dict(group=group,count=len(rs),RHS_min=min(vals) if vals else None,RHS_median=statistics.median(vals) if vals else None,RHS_max=max(vals) if vals else None,RHS_sum=sum(vals),wall_sum=sum(r.get('wall_seconds',0) for r in rs),native_failure_count=sum(r.get('native_stationarity',{}).get('status')!='COMPUTED' for r in rs),probit_success_count=sum(r.get('probit_stationarity',{}).get('optimizer',{}).get('success') is True for r in rs)))
write('HISTORICAL_CONTEXT_COSTS.json',dict(groups=cost,new_RHS=0,source=str(R4/'results/contexts'),statistical_scope='historical context JSON only, not new runtime projection'))
oldmoment={}
for parent in sorted(historical_ids):
 seed=int(np.random.SeedSequence([20261025121,*parent.encode('utf8')]).generate_state(1)[0]);oldmoment.setdefault(seed,[]).append(parent)
write('ID_DERIVED_MOMENT_RNG_AUDIT.json',dict(rule='SeedSequence([20261025121,*UTF8(parent)]).generate_state(1)[0] uint32',historical_parent_count=len(historical_ids),historical_collisions=[v for v in oldmoment.values() if len(v)>1],pending_to_historical=[dict(parent=parent,seed=seed,historical=oldmoment[seed]) for parent,seed in zip(allids,allmoment) if seed in oldmoment],pending_within_collisions=len(allmoment)-len(set(allmoment)),scope='Frozen K32 moment Monte Carlo only. These are not latent/noise dataset RNG seeds. No Python built-in hash is used.'))
print(json.dumps(dict(cohorts=len(cohorts),historical_ids=len(historical_ids),pending_counts={k:v['count'] for k,v in pending.items()},checkpoints=len(checkpoints),asset_headers=len(assets),cost=cost)))
