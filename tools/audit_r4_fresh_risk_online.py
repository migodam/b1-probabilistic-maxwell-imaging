"""Metadata-only closure audit; never materializes observations or heldout truth."""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
from p3r4.evaluation import checkpoint_identity,write_new_json
from r4_fresh_risk import component_check
S=ROOT/'research/p3l_r4_conditional_transport';F=S/'results/fresh';O=F/'risk_online_queue_v1'
def read(p):return json.loads(Path(p).read_text())
def check(v):
 if checkpoint_identity(v['path'])!=v:raise ValueError('sealed identity changed: '+v['path'])
def main():
 start=time.monotonic();c=read(S/'AUDIT/risk/fresh_freeze_v3/RISK_COMPONENT_FROZEN.json');component_check(c);prereq=read(S/'results/FRESH_PHYSICS_AND_CONTEXTS_COMPLETE.json')
 if prereq['status']!='COMPLETE':raise ValueError('global prerequisites not complete')
 marker=checkpoint_identity(S/'results/FRESH_PHYSICS_AND_CONTEXTS_COMPLETE.json');summary={};receipts=[]
 for cohort,n in [('locked',512),('confirmation',256)]:
  base=F/cohort/'risk45';g=read(base/'online/ONLINE_COMPLETE.json');f=read(F/cohort/'BOUND_RISK.json');top=read(S/'configs'/('FROZEN_'+cohort.upper()+'.json'));q=read(base/'QUEUE_COMPLETE.json')
  if g['freeze']!=f or g['truth_read'] or len(g['parents'])!=n or len(f['methods'])!=45 or g['canonical_batch_size']!=64 or g['array_codec']!='stored' or q['scored']:raise ValueError('online contract mismatch')
  if f['risk_rule']!=c['risk_rule'] or top['risk_rule']!=c['risk_rule']:raise ValueError('risk rule adapted')
  for identity in g['sealed_artifact_identities']+g['source_identities']:check(identity)
  featuregate=read(base/'features/FEATURES_COMPLETE.json')
  if any(v['mtime_ns']<=f['freeze_time_ns'] for v in featuregate['observation_sources']):raise ValueError('observation predates bound risk freeze')
  if (base/'scores').exists():raise ValueError('unexpected risk score directory before global GT release')
  for label in ['features','infer','export']:
   job='fresh_risk_'+cohort+'_'+label+'_v1';d=S/'IV/jobs'/job;r=read(d/'CLOSED.json');st=read(d/'START.json')
   if r['status']!='COMPLETE' or st['stage']!='IV' or st['kind']!='analysis' or int(st['start_unix']*1e9)<marker['mtime_ns'] or r['counters'].get('RHS',0)!=0:raise ValueError('receipt/gate/RHS mismatch')
   receipts.append(dict(job=job,identity=checkpoint_identity(d/'CLOSED.json'),wall_seconds=r['controlled_wall_seconds'],CPU_seconds=r['profile']['CPU_seconds'],RSS_peak_bytes=r['profile']['RSS_peak_bytes'],RHS=r['counters'].get('RHS',0)))
  raw=read(base/'exports/RAW_FORWARD_FROZEN.json')
  if raw['stage']!='IV' or len(raw['rows'])!=6*n or raw['planned_RHS']!=6*n*8:raise ValueError('raw diagnostic interface mismatch')
  summary[cohort]=dict(parents=n,methods=45,failures_by_method={k:len(v) for k,v in g['failures'].items()},total_failures=sum(map(len,g['failures'].values())),online_gate_identity=checkpoint_identity(base/'online/ONLINE_COMPLETE.json'),sealed_outputs=len(g['sealed_artifact_identities']),sealed_sources=len(g['source_identities']),raw_forward_spec_identity=checkpoint_identity(base/'exports/RAW_FORWARD_FROZEN.json'),planned_raw_forward_RHS=raw['planned_RHS'],CPU_seconds=g['cpu_seconds'],wall_seconds=g['wall_seconds'],RSS_peak_bytes=g['peak_rss_bytes'],GT_read=False)
 p=F/'locked/representatives_prepared/FROZEN_REQUEST.json';req=read(p)
 if req['stage']!='IV' or req['GT_read'] or req['manifest']['parents']!=[f'R4L1_{i:05d}' for i in range(64)] or req['maximum_total_RHS']!=81920:raise ValueError('representative contract mismatch')
 for v in req['source_identities']:check(v)
 for item in req['qualified']:check(item['mu_identity'])
 d=S/'IV/jobs/fresh_locked_representatives_prepare_v1';r=read(d/'CLOSED.json');st=read(d/'START.json')
 if r['status']!='COMPLETE' or st['stage']!='IV' or st['kind']!='analysis' or r['counters'].get('RHS',0)!=0:raise ValueError('representative prepare receipt mismatch')
 receipts.append(dict(job=d.name,identity=checkpoint_identity(d/'CLOSED.json'),wall_seconds=r['controlled_wall_seconds'],CPU_seconds=r['profile']['CPU_seconds'],RSS_peak_bytes=r['profile']['RSS_peak_bytes'],RHS=r['counters'].get('RHS',0)))
 write_new_json(O/'FINAL_AUDIT.json',dict(status='BOTH_COHORT_RISK_ONLINE_SEALED_NO_GT',global_prerequisite_identity=marker,cohorts=summary,receipts=receipts,total_controlled_wall_seconds=sum(v['wall_seconds'] for v in receipts),total_CPU_seconds=sum(v['CPU_seconds'] for v in receipts),peak_recorded_RSS_bytes=max(v['RSS_peak_bytes'] for v in receipts),new_Maxwell_RHS=0,new_fits=0,heldout_GT_read=False,scoring_run=False,representative=dict(request_identity=checkpoint_identity(p),declared_parents=64,qualified=len(req['qualified']),skips=req['skips'],physical_released=False,raw_clipped_each=True,estimators=['fusion_spline_s2026102501','frozen_r2'],shared_ordinary_per_parent=1),frozen_model_and_source_identities_unchanged=True,density_status='not claimed by risk worker; parent must check both density cohorts before ANYGTscore',audit_wall_seconds=time.monotonic()-start))
 print(json.dumps(dict(status='BOTH_COHORT_RISK_ONLINE_SEALED_NO_GT',total_wall=sum(v['wall_seconds'] for v in receipts),representative_qualified=len(req['qualified']),skips=len(req['skips']))))
if __name__=='__main__':main()
