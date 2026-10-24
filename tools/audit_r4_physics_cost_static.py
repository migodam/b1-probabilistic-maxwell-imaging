#!/usr/bin/env python3
"""Independent static closed-job/cost audit. No experiment/model launch."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
from pathlib import Path
from collections import Counter,defaultdict
import sys,json,time,traceback,fcntl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write,CAPS,STAGE_WALL,TOTAL_RHS

JOB='independent_physics_cost_review_v1'
DEST=OUT/'AUDIT/independent_physics_cost_review'


def main():
 if DEST.exists():raise FileExistsError('independent audit exists; preserve prior audit')
 DEST.mkdir(parents=True);b=Book('II',JOB,max_wall=900,analysis_only=True);status='COMPLETE';error=None;started=time.perf_counter()
 try:
  closed={}
  for p in OUT.glob('*/jobs/*/CLOSED.json'):
   d=json.loads(p.read_text());closed[p.parent.name]=dict(stage=p.parent.parent.parent.name,path=str(p),receipt=d)
  # Briefly synchronize only the committed byte boundary, then release writer lock.
  with (OUT/'budget.lock').open('r') as lock:
   fcntl.flock(lock,fcntl.LOCK_SH);limit=(OUT/'RESOURCE_LEDGER.jsonl').stat().st_size;fcntl.flock(lock,fcntl.LOCK_UN)
  jobs=defaultdict(lambda:dict(counts=Counter(),events=Counter(),roles=Counter(),first_sequence=None,last_sequence=None,start=None,end=None,full_event_lines=[],extra_forward_lines=[]))
  invariants=[];previous=None;last=None;row_count=0;stage_sums=defaultdict(Counter)
  with (OUT/'RESOURCE_LEDGER.jsonl').open('rb') as f:
   while f.tell()<limit:
    raw=f.readline()
    if not raw:break
    row=json.loads(raw);row_count+=1;seq=row['sequence'];job=row['job'];j=jobs[job]
    if previous is not None and seq!=previous['sequence']+1:invariants.append(dict(kind='SEQUENCE_GAP',sequence=seq,previous=previous['sequence']))
    for key in ('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds'):
     amount=row.get(key,0)
     if amount<0:invariants.append(dict(kind='NEGATIVE_INCREMENT',sequence=seq,key=key,value=amount))
     j['counts'][key]+=amount;stage_sums[row['stage']][key]+=amount
    if previous is not None:
     for stage in row['after']['stages']:
      for key in ('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds'):
       expected=previous['after']['stages'][stage][key]+(row.get(key,0) if stage==row['stage'] else 0)
       got=row['after']['stages'][stage][key]
       if abs(got-expected)>1e-7:invariants.append(dict(kind='AFTER_MISMATCH',sequence=seq,stage=stage,key=key,expected=expected,actual=got))
    j['events'][row['event']]+=1;j['roles'][row.get('role','')]+=row.get('RHS',0)
    if j['first_sequence'] is None:j['first_sequence']=seq
    j['last_sequence']=seq
    if row['event']=='JOB_START':j['start']={k:row.get(k) for k in ('sequence','stage','analysis_only','gpu_job','kind','reserve_wall','reserve_wall_seconds')}
    if row['event']=='JOB_END':j['end']=row
    if row['event']=='FULL_FORWARD_ATTEMPT':j['full_event_lines'].append(seq)
    if row['event']=='SOLVER_RHS' and 'full_forward' in row.get('role',''):j['extra_forward_lines'].append(seq)
    previous=row;last=row
  selected=[];snapshot_inventory=[];issues=[]
  for name,meta in sorted(closed.items()):
   if name==JOB:continue
   j=jobs.get(name);receipt=meta['receipt']
   if j is None or j['end'] is None:
    issues.append(dict(kind='CLOSED_WITHOUT_SNAPSHOTTED_JOB_END',job=name,evidence=meta['path']));continue
   for key in ('controlled_wall_seconds',):
    if abs(receipt[key]-j['counts'][key])>1e-7:issues.append(dict(kind='CLOSED_LEDGER_WALL_MISMATCH',job=name,key=key,receipt=receipt[key],ledger=j['counts'][key]))
   if j['events']['JOB_START']!=1 or j['events']['JOB_END']!=1:issues.append(dict(kind='DUPLICATE_START_END',job=name,events=dict(j['events'])))
   if receipt.get('unused_forward_credit',0):issues.append(dict(kind='UNUSED_FORWARD_CREDIT',job=name,value=receipt['unused_forward_credit']))
   if j['extra_forward_lines']:issues.append(dict(kind='POTENTIAL_FORWARD_DOUBLE_BILLING',job=name,sequences=j['extra_forward_lines']))
   profile=receipt.get('profile',{});full_calls=receipt.get('counters',{}).get('full_forward_calls',0)
   if full_calls!=j['events']['FULL_FORWARD_ATTEMPT']:
    # Iterative states are separate; compare only the actual Full event type.
    issues.append(dict(kind='FULL_CALL_EVENT_MISMATCH',job=name,profile_full_calls=full_calls,ledger_full_events=j['events']['FULL_FORWARD_ATTEMPT']))
   selected.append(dict(job=name,stage=meta['stage'],status=receipt['status'],actual_ledger_cost=dict(j['counts']),
    events=dict(j['events']),RHS_by_role=dict(j['roles']),start=j['start'],first_sequence=j['first_sequence'],last_sequence=j['last_sequence'],
    source_receipt=meta['path'],closed_profile=profile,linear_solve_rhs_counter=receipt.get('counters',{}).get('linear_solve_rhs'),
    physics_wall_over_controlled=j['counts']['physics_wall_seconds']/j['counts']['controlled_wall_seconds'] if j['counts']['controlled_wall_seconds'] else None))
   source=Path(meta['path']).parent/'SOURCE_REVISION.json'
   record=dict(job=name,stage=meta['stage'],has_SOURCE_REVISION=source.exists(),exact_old_source_recovery_claim=False)
   if source.exists():
    s=json.loads(source.read_text());base=Path(s['source_snapshot']);argv=s.get('argv',[]);entry=Path(argv[0]).name if argv else None
    tools=base/'tools';record.update(argv=argv,registered_source_files=len(s.get('source_files',[])),
      entrypoint_snapshot_present=bool(entry and (tools/entry).exists()),entrypoint=entry,
      git_revision=s.get('git_revision'),source_snapshot=str(base),has_tools_directory=tools.exists())
    if entry and entry.endswith('.py') and not record['entrypoint_snapshot_present']:
     issues.append(dict(kind='ENTRYPOINT_NOT_SNAPSHOTTED',job=name,entrypoint=entry,source_revision=str(source),has_tools_directory=tools.exists()))
   else:
    result_source=OUT/'results'/meta['stage']/name/'source'
    record['own_result_source_snapshot_exists']=result_source.exists()
   snapshot_inventory.append(record)
  byname={z['job']:z for z in selected}
  checks={}
  for name,expected in [('finite_feedback_reuse_first',5955),('finite_feedback_pcn_decision_audit',2048),
                        ('conservative_reduced_rhs_reconciliation',6144),('decision_refresh_star_v1',62770),
                        ('interval_chain_v1_full',65536),('interval_chain_v1_interval',53184),
                        ('dynamic_geometry_first',291648),('dynamic_geometry_heun_r1',168064)]:
   actual=byname[name]['actual_ledger_cost']['RHS'] if name in byname else None
   checks[name]=dict(expected_RHS=expected,ledger_RHS=actual,PASS=actual==expected)
   if actual!=expected:issues.append(dict(kind='NAMED_RHS_TOTAL_MISMATCH',job=name,expected=expected,actual=actual))
  checks['reuse_all_complex']=dict(full_grid=5955+2048,reduced_reconciliation=6144,total=14147,
      actual=sum(byname[n]['actual_ledger_cost']['RHS'] for n in ('finite_feedback_reuse_first','finite_feedback_pcn_decision_audit','conservative_reduced_rhs_reconciliation')))
  checks['interval_pair']=dict(total=118720,actual=sum(byname[n]['actual_ledger_cost']['RHS'] for n in ('interval_chain_v1_full','interval_chain_v1_interval')))
  checks['geometry_pair']=dict(total=459712,actual=sum(byname[n]['actual_ledger_cost']['RHS'] for n in ('dynamic_geometry_first','dynamic_geometry_heun_r1')))
  active=[dict(job=name,start=j['start'],last_sequence=j['last_sequence'],exclusion='no CLOSED receipt at audit cutoff; committed budget increments only, no live result read') for name,j in jobs.items() if name not in closed and name!=JOB]
  write(DEST/'LEDGER_RECONCILIATION.json',dict(status='STATIC_CLOSED_JOB_RECONCILIATION',cutoff_bytes=limit,last_sequence=last['sequence'],rows=row_count,
     invariants=invariants,issues=issues,named_checks=checks,closed_jobs=selected,excluded_active_jobs=active,
     stage_committed_totals={s:dict(v) for s,v in stage_sums.items()},ledger_after_at_cutoff=last['after'],
     no_new_complex_RHS=0,no_new_states=0,no_new_LU=0))
  write(DEST/'SOURCE_SNAPSHOT_INVENTORY.json',dict(jobs=snapshot_inventory))
  print(json.dumps(dict(closed_jobs=len(selected),excluded_active=len(active),invariant_errors=len(invariants),issues=issues,named_checks=checks,elapsed=time.perf_counter()-started),indent=2))
 except BaseException as exc:
  status='ERROR';error=repr(exc);write(DEST/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:
  b.close(status,error);write(DEST/'FINAL_AUDIT_COST.json',dict(status=status,error=error,new_complex_RHS=0,new_states=0,new_LU=0,elapsed_seconds=time.perf_counter()-started,
   receipt=str(b.jobdir/'CLOSED.json')))

if __name__=='__main__':main()
