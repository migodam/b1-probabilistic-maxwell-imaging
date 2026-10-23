"""Parent-invoked metadata closure after every registered physics job settles.

No Book/backend imports, no physics, no gate inferred from a short job receipt.
Requires explicit reviewed REDIRECT decision, complete tails, matching receipts,
and all original stage ceilings; never deletes prior seals or resets a budget.
"""
import argparse,collections,datetime as dt,fcntl,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text())
def exclusive(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--reviewed-decision',required=True,choices=['REDIRECT'])
 args=ap.parse_args()
 if (ROOT/'CLOSED.json').exists():raise RuntimeError('Already closed; never rewrite')
 if load(ROOT/'B2_FINISH_QUEUE_CHECKPOINT.json')['status']!='B2_VALIDATION_TAIL_FINISHED_PENDING_PARENT_REVIEW':raise RuntimeError('B2 tail not settled')
 if load(ROOT/'P1B_QUEUE_CHECKPOINT.json')['status']!='P1B_NUMERICAL_WORK_FINISHED_PENDING_SCIENTIFIC_REVIEW':raise RuntimeError('P1B queue not settled')
 if load(ROOT/'results/P1B/green_reuse_control/SUMMARY.json')['status']!='ENGINEERING_CONTROL_COMPLETE':raise RuntimeError('Engineering branch not complete')
 if load(ROOT/'results/P1B/offline_scores/COMPLETE.json')['observations']!=32:raise RuntimeError('Offline population incomplete')
 cfg=load(ROOT/'configs/AUTHORIZATION.json'); path=ROOT/'PHYSICS_LEDGER.jsonl'
 with path.open() as f:
  fcntl.flock(f,fcntl.LOCK_SH);rows=[json.loads(z) for z in f if z.strip()]
 reserves=[r for r in rows if r.get('event')=='reserve']
 ids=[r['id'] for r in reserves]
 if len(ids)!=len(set(ids)):raise RuntimeError('Duplicate reservation IDs')
 jobs=[r for r in rows if r.get('event')=='job_cost'];labels=[(r['stage'],r['label']) for r in jobs]
 startup=[r for r in rows if r.get('event')=='launcher_startup_failure']
 if len(labels)!=len(set(labels)):raise RuntimeError('Duplicate job charge')
 for folder in (ROOT/'results/jobs').iterdir():
  if not folder.is_dir():continue
  if not (folder/'RECEIPT.json').is_file():raise RuntimeError('Unsettled job: '+folder.name)
  receipt=load(folder/'RECEIPT.json')
  if receipt.get('event')=='launcher_startup_failure':
   if receipt.get('status')!='FAILED_PRE_EXEC' or receipt.get('new_rhs')!=0 or receipt.get('new_LU')!=0:raise RuntimeError('Unreviewed startup-failure semantics: '+folder.name)
   source=startup
  else:source=jobs
  matching=[r for r in source if r['stage']==receipt['stage'] and r['label']==folder.name]
  if len(matching)!=1 or abs(matching[0]['charged_cpu_seconds']-receipt['charged_cpu_seconds'])>1e-8:raise RuntimeError('Receipt not settled exactly once: '+folder.name)
 totals={}
 for stage,cap in cfg['stages'].items():
  t=dict(rhs=sum(r.get('rhs',0) for r in reserves if r['stage']==stage),
   factorizations=sum(r.get('factorizations',0) for r in reserves if r['stage']==stage),
   cpu_seconds=sum(r.get('charged_cpu_seconds',0) for r in rows if r.get('stage')==stage),
   controlled_job_wall_seconds=sum(r.get('wall_seconds',0) for r in jobs if r['stage']==stage))
  if t['rhs']>cap['rhs'] or t['factorizations']>cap['factorizations'] or t['cpu_seconds']>cap['cpu_seconds'] or t['controlled_job_wall_seconds']>cap['cpu_seconds']:
   raise RuntimeError('Stage hard cap exceeded: '+stage)
  totals[stage]=dict(used=t,cap=cap,unused_rhs=cap['rhs']-t['rhs'],no_transfer=True)
 utc=dt.datetime.now(dt.timezone.utc).isoformat()
 for stage,t in totals.items():
  target=ROOT/f'{stage}_SEALED.json'
  if target.exists():
   if stage!='R':
    existing=load(target)
    if existing.get('status')!='STAGE_CLOSED' or existing.get('stage')!=stage or existing.get('totals')!=t['used'] or existing.get('caps')!=t['cap']:raise RuntimeError('Existing stage seal differs; inspect, never overwrite')
  else:exclusive(target,dict(status='STAGE_CLOSED',utc=utc,stage=stage,totals=t['used'],caps=t['cap'],remaining_not_new_authorization=True))
 total={k:sum(t['used'][k] for t in totals.values()) for k in totals['R']['used']}
 resources=dict(status='CLOSED_ACCOUNTING_SNAPSHOT',authorization=cfg['stage_authorization'],utc=utc,
  stages=totals,total=total,charged_reservations_including_failures=len(reserves),
  failed_jobs=[dict(label=r['label'],stage=r['stage'],CPU=r['charged_cpu_seconds'],wall=r['wall_seconds']) for r in jobs if r['returncode']!=0],
  pre_exec_failures=[dict(label=r['label'],stage=r['stage'],CPU=r['charged_cpu_seconds'],wall=r.get('wall_seconds'),CPU_attribution=r.get('conservative_cpu_attribution'),new_rhs=r.get('new_rhs'),new_LU=r.get('new_LU')) for r in startup],
  max_completed_wrapper_peak_RSS_bytes=max(r['peak_rss_bytes'] for r in jobs),
  max_sampled_wrapper_RSS_bytes=max((r['peak_sampled_rss_bytes'] for r in jobs if r.get('peak_sampled_rss_bytes') is not None),default=None),
  receipts_without_sampled_RSS=[r['label'] for r in jobs if r.get('peak_sampled_rss_bytes') is None],
  nested_operation_CPU_not_added_twice=True,failed_uncompleted_reservations_remain_charged=True,
  historical_P0_budgets_unchanged=True,GPU=0,new_3D=0,
  independent_parent_counts=dict(P1A=16,P1B=16,noise_replicas_each=2,two_fixed_wall_fields=True),
  postprocessing='Separately profiled pure audits/plots in AUDIT and tables; manual review/document writing not completely CPU-profiled. Not declared zero. Physics/job figures include wrapper persistence/failed work.',
  all_local_raw_fields_uploaded=False,public_fields='First validated joint cal/joint two-world five-frequency subset; Release inventory authoritative')
 dense=[r for r in rows if r.get('event')=='engineering_operation' and r.get('dense_core_LU')]
 resources['Green_reduced_dense_LU_separate_from_full_FE_LU']=dict(
  total=sum(r['dense_core_LU'] for r in dense),
  dimensions=dict(collections.Counter(str(r['dense_core_dimension']) for r in dense)),
  CPU_already_in_wrapper=True)
 resources['Green_sparse_MVP_count']=sum(r.get('matrix_vector_products',0) for r in rows if r.get('event')=='Krylov_MVP_and_base_action_cost')
 (ROOT/'RESOURCE_LEDGER.json').write_text(json.dumps(resources,ensure_ascii=False,indent=2)+'\n')
 closed=dict(status='CAMPAIGN_CLOSED',utc=utc,scientific_decision=args.reviewed_decision,
  total=total,new_authorization_required_for_any_physics=True,old_budgets_unchanged=True,
  original_geometry_prior='minimum_feature unresolved; no retroactive PASS',
  observed_scope='oracle fixed-support witness and newly frozen controlled probability benchmark, not general online imaging')
 exclusive(ROOT/'CLOSED.json',closed)
 # An informational event cannot change any paid-unit or CPU total.
 with path.open('a') as f:
  fcntl.flock(f,fcntl.LOCK_EX);f.write(json.dumps(dict(event='campaign_closed',utc=utc,decision=args.reviewed_decision))+'\n');f.flush()
 print(json.dumps(closed,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
