"""Read-only W execution/accounting/persistence audit; stdlib only, no physics.

Inputs: root ledger, AUTHORIZATION/B2_ALLOCATION, wrapper RECEIPT metadata,
indexed endpoint file existence plus a bounded safe metadata sample, and small
source snippets documenting accounting/persistence semantics. Never opens
arrays, observations, OFFLINE_TRUTH, GENERATION_PLAN, or job OUTPUT.log files.
Only AUDIT/execution_integrity/{RESULTS.json,REPORT.md} are written.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict,deque
from datetime import datetime,timezone
import fcntl
import json
import math
from pathlib import Path
import time

DEFAULT_ROOT=Path(__file__).resolve().parents[1]


def _utc(value):
    if not isinstance(value,str):return None
    try:return datetime.fromisoformat(value.replace('Z','+00:00')).timestamp()
    except ValueError:return None


def _read(path):return json.loads(path.read_text())


def _write(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def _snapshot(path):
    with path.open('rb') as stream:
        fcntl.flock(stream,fcntl.LOCK_SH)
        raw=stream.read()
    lines=raw.splitlines(keepends=True)
    rows=[];errors=[];partial=None
    for number,line in enumerate(lines,1):
        if not line.strip():continue
        try:
            value=json.loads(line)
            if not isinstance(value,dict):raise ValueError('event must be an object')
            rows.append((number,value))
        except (ValueError,UnicodeDecodeError) as exc:
            info=dict(line=number,error=str(exc))
            if number==len(lines) and not line.endswith(b'\n'):
                partial=info
            else:errors.append(info)
    return rows,dict(bytes_read=len(raw),line_count=len(lines),parsed_events=len(rows),
                     malformed_complete_lines=errors,trailing_partial_line=partial,
                     last_line_terminated=not raw or raw.endswith(b'\n'))


def _role(label,solver=False):
    if solver:
        return 'forward_sources' if label is None else label
    if label=='fem_material_frequency_LU':return 'factor'
    if isinstance(label,str) and label.startswith('fem_'):return label[4:]
    return label


def _context(row):
    return (row.get('stage'),row.get('pid'),row.get('phase'),row.get('scope'))


def _num(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


def _sensitive(path):
    parts=[p.lower() for p in path.parts]
    return ('data' in parts or any('offline_truth' in p or 'generation_plan' in p
                                 or 'observations' in p for p in parts))


def _source_evidence(root):
    """Current source evidence; not a claim all historical snapshots are equal."""
    requests=[
      (root/'src/wr/accounting.py','reserves_are_append_only',['def reserve(',"event='budget_reject' if reason else 'reserve'","f.seek(0,2)"]),
      (root/'src/wr/accounting.py','ledger_is_authoritative',['def totals(',"r.get('event')=='reserve'","sum(r.get('charged_cpu_seconds'"]),
      (root/'src/wr/accounting.py','scope_cap_is_not_reset',['class ScopedBook','scope_rhs_limit']),
      (root/'tools/run_job.py','all_child_process_CPU_and_single_ledger_charge',['resource.getrusage(resource.RUSAGE_CHILDREN)',"cpu=after.ru_utime",'book.event(**row)']),
      (root/'src/wr/physics.py','frequency_files_before_endpoint_completion',["np.savez_compressed(dest/f'f{idx}.npz'","write(dest/'COMPLETE.json'",'max five frequencies']),
      (root.parent/'research_w_continuation/src/wq/fem.py','reserve_before_action_and_nested_wall_timers',["'fem_forward_sources'","factor = splu",'timings = dict(',"self._event('backend_solve'"])]
    evidence=[]
    for path,purpose,tokens in requests:
        if not path.exists():
            evidence.append(dict(purpose=purpose,path=str(path),status='SOURCE_NOT_AVAILABLE'));continue
        lines=path.read_text().splitlines();hits=[]
        for token in tokens:
            matches=[(i+1,line.strip()) for i,line in enumerate(lines) if token in line]
            if matches:
                hits.append(dict(token=token,line=matches[0][0],code=matches[0][1]))
        evidence.append(dict(purpose=purpose,path=str(path),status='CURRENT_SOURCE_ONLY',snippets=hits))
    return evidence


def audit(root=DEFAULT_ROOT,final=False,metadata_samples=24):
    root=Path(root).resolve()
    if not 0<=metadata_samples<=128:raise ValueError('metadata sample bound must be 0..128')
    start=time.perf_counter()
    ledger_path=root/'PHYSICS_LEDGER.jsonl'
    rows,snapshot=_snapshot(ledger_path)
    authorization=_read(root/'configs/AUTHORIZATION.json')
    allocations=_read(root/'configs/B2_ALLOCATION.json')
    findings=[]
    def finding(level,code,message,**details):
        findings.append(dict(level=level,code=code,message=message,**details))
    if snapshot['malformed_complete_lines']:
        finding('FAIL','MALFORMED_LEDGER_LINES','Complete ledger lines cannot be parsed',lines=snapshot['malformed_complete_lines'])
    if snapshot['trailing_partial_line']:
        finding('UNVERIFIED' if final else 'IN_PROGRESS','TRAILING_PARTIAL_APPEND',
                'Trailing unterminated event retained as unsettled, never counted as a successful charge')
    event_counts=Counter(row.get('event','<missing>') for _,row in rows)
    reserves=[(n,r) for n,r in rows if r.get('event')=='reserve']
    ids=defaultdict(list)
    for n,r in rows:
        if r.get('event') in ('reserve','budget_reject') and r.get('id') is not None:ids[str(r['id'])].append(n)
    duplicates=[dict(id=key,lines=value) for key,value in ids.items() if len(value)>1]
    missing_ids=[n for n,r in reserves if not r.get('id')]
    if duplicates:finding('FAIL','DUPLICATE_APPEND_IDS','Reservation/rejection IDs repeat',duplicates=duplicates)
    if missing_ids:finding('UNVERIFIED','RESERVE_WITHOUT_ID','Reservations cannot be uniquely identified',lines=missing_ids)
    invalid_charges=[]
    for n,r in rows:
        fields=('rhs','factorizations') if r.get('event')=='reserve' else ()
        for key in (*fields,'charged_cpu_seconds'):
            value=r.get(key,0)
            if not _num(value) or value<0 or (key in fields and value!=int(value)):
                invalid_charges.append(dict(line=n,field=key,value=repr(value)))
    if invalid_charges:finding('FAIL','INVALID_CHARGE_FIELDS','Negative/nonfinite/nonintegral charges detected',rows=invalid_charges)

    stages={}
    for stage,cap in authorization['stages'].items():
        rr=[r for _,r in reserves if r.get('stage')==stage]
        events=[r for _,r in rows if r.get('stage')==stage]
        totals=dict(rhs=sum(r.get('rhs',0) for r in rr),
                    factorizations=sum(r.get('factorizations',0) for r in rr),
                    cpu_seconds=sum(r.get('charged_cpu_seconds',0) for r in events))
        breached=[key for key in totals if totals[key]>cap[key]+(1e-9 if key=='cpu_seconds' else 0)]
        stages[stage]=dict(ledger_totals=totals,caps=cap,
                          remaining={k:cap[k]-totals[k] for k in totals},
                          reserve_count=len(rr),sealed=(root/f'{stage}_SEALED.json').exists(),
                          CPU_status='settled ledger only; unreceipted live CPU is not yet included')
        if breached:finding('FAIL','STAGE_CAP_EXCEEDED','Root-ledger charges exceed a stage cap',stage=stage,fields=breached)
    unknown_stages=sorted({str(r.get('stage')) for _,r in reserves if r.get('stage') not in stages})
    if unknown_stages:finding('FAIL','UNKNOWN_CHARGED_STAGE','Reserved work has no authorized stage',stages=unknown_stages)
    scopes={}
    b2=[r for _,r in reserves if r.get('stage')=='B2']
    for scope,cap in allocations['scopes'].items():
        rr=[r for r in b2 if r.get('scope')==scope]
        totals=dict(rhs=sum(r.get('rhs',0) for r in rr),factorizations=sum(r.get('factorizations',0) for r in rr))
        scopes[scope]=dict(ledger_totals=totals,caps=cap,remaining_rhs=cap['rhs']-totals['rhs'],
                           CPU_status='no independent scope CPU cap is declared; CPU charged at B2 process/stage level')
        if totals['rhs']>cap['rhs']:finding('FAIL','SCOPE_RHS_CAP_EXCEEDED','Existing B2 scope allowance exceeded',scope=scope)
    b2_unscoped=[n for n,r in reserves if r.get('stage')=='B2' and r.get('scope') not in allocations['scopes']]
    if b2_unscoped:finding('UNVERIFIED','B2_UNSCOPED_RESERVES','B2 charges cannot be assigned to frozen named scope',lines=b2_unscoped)

    # A receipt is evidence of an existing ledger charge, never a second charge.
    job_costs=defaultdict(list)
    for n,r in rows:
        if r.get('event')=='job_cost':job_costs[(r.get('stage'),r.get('label'))].append((n,r))
    for key,values in job_costs.items():
        if len(values)>1:finding('FAIL','DUPLICATE_JOB_COST','One wrapper label is charged more than once',stage=key[0],label=key[1],lines=[n for n,_ in values])
    receipt_records=[];receipts_by_label={}
    jobs=root/'results/jobs'
    for path in sorted(jobs.glob('*/RECEIPT.json')):
        try:record=_read(path)
        except (ValueError,OSError) as exc:
            finding('UNVERIFIED','UNREADABLE_JOB_RECEIPT',str(exc),path=str(path));continue
        label=record.get('label',path.parent.name);stage=record.get('stage')
        linked=job_costs.get((stage,label),[])
        startup=[(n,r) for n,r in rows if r.get('event')=='launcher_startup_failure'
                 and r.get('stage')==stage and r.get('label')==label]
        evidence=linked or startup
        state='MATCHED_ROOT_LEDGER' if len(evidence)==1 else 'UNSETTLED_RECEIPT_WITHOUT_UNIQUE_LEDGER_COST'
        if len(evidence)==1:
            n,charge=evidence[0]
            if not math.isclose(float(record.get('charged_cpu_seconds',0)),float(charge.get('charged_cpu_seconds',0)),rel_tol=1e-9,abs_tol=1e-8):
                state='CPU_MISMATCH';finding('FAIL','RECEIPT_CPU_MISMATCH','Receipt CPU differs from authoritative ledger charge',label=label,receipt=str(path),ledger_line=n)
            if linked and record.get('returncode')!=charge.get('returncode'):
                finding('FAIL','RECEIPT_RETURNCODE_MISMATCH','Receipt returncode differs from ledger',label=label)
        elif final:
            finding('UNVERIFIED','RECEIPT_PENDING_SETTLEMENT','Receipt exists without a unique corresponding ledger event',label=label)
        else:
            finding('IN_PROGRESS','RECEIPT_PENDING_SETTLEMENT','Receipt publication may precede its final ledger append',label=label)
        item=dict(label=label,stage=stage,status=record.get('status'),returncode=record.get('returncode'),
                  charged_cpu_seconds=record.get('charged_cpu_seconds',0),
                  termination_reason=record.get('termination_reason'),association=state,
                  peak_rss_bytes=record.get('peak_rss_bytes'),
                  peak_sampled_rss_bytes=record.get('peak_sampled_rss_bytes'),
                  receipt=str(path),ledger_lines=[n for n,_ in evidence],
                  source_snapshot_exists=(path.parent/'SOURCE_SNAPSHOT').is_dir(),
                  config_snapshot_exists=(path.parent/'CONFIG_SNAPSHOT').is_dir())
        peak=max((v for v in (item['peak_rss_bytes'],item['peak_sampled_rss_bytes']) if _num(v)),default=0)
        if peak>authorization['max_process_RSS_bytes']:
            finding('WARN','RECEIPT_RSS_LIMIT_EXCEEDED','Closed wrapper reports RSS above the configured process cap; retain termination/charge evidence',
                    label=label,reported_peak_bytes=peak,cap_bytes=authorization['max_process_RSS_bytes'])
        receipt_records.append(item);receipts_by_label[(stage,label)]=item
    for key,values in job_costs.items():
        if key not in receipts_by_label:
            finding('UNVERIFIED','LEDGER_JOB_WITHOUT_RECEIPT','Ledger cost remains authoritative but wrapper receipt missing',stage=key[0],label=key[1])
    unreceipted_jobs=[p.name for p in sorted(jobs.iterdir()) if p.is_dir() and not (p/'RECEIPT.json').exists()] if jobs.exists() else []

    # No direct child PID / reserve-ID foreign key exists in wrapper receipts.
    # Conservatively associate whole physical child-PID spans with stage/time
    # windows. Report ambiguity; do not invent stronger provenance.
    processes=defaultdict(list)
    for n,r in rows:
        if r.get('event') in ('reserve','backend_factor','backend_solve','backend_state','material_frequency_attempt') and r.get('pid') is not None:
            processes[(r.get('stage'),r['pid'])].append((n,r))
    process_info={}
    for key,events in processes.items():
        times=[_utc(r.get('utc')) for _,r in events];times=[t for t in times if t is not None]
        candidates=[]
        if times:
            for (stage,label),costs in job_costs.items():
                if stage!=key[0]:continue
                for n,cost in costs:
                    end=_utc(cost.get('utc'));wall=cost.get('wall_seconds')
                    if end is not None and _num(wall) and min(times)>=end-wall-2 and max(times)<=end+1:
                        candidates.append(dict(label=label,returncode=cost.get('returncode'),job_cost_line=n))
        status=('SETTLED_SUCCESS' if candidates[0]['returncode']==0 else 'SETTLED_FAILED') if len(candidates)==1 else ('AMBIGUOUS' if candidates else 'UNSETTLED')
        reservations=[r for _,r in events if r.get('event')=='reserve']
        process_info[key]=dict(stage=key[0],pid=key[1],event_lines=[events[0][0],events[-1][0]],
                               contextual_reserved_rhs=sum(r.get('rhs',0) for r in reservations),
                               contextual_reserved_LU=sum(r.get('factorizations',0) for r in reservations),
                               indexed_material_frequency_attempts=sum(r.get('event')=='material_frequency_attempt' for _,r in events),
                               association=status,job_candidates=candidates,
                               attribution='whole-PID stage/time containment only; no reserve/child-PID receipt foreign key')

    for receipt in receipt_records:
        related=[p for p in process_info.values() if len(p['job_candidates'])==1
                 and p['stage']==receipt['stage'] and p['job_candidates'][0]['label']==receipt['label']]
        receipt['contextual_physical_charges']=dict(
            rhs=sum(p['contextual_reserved_rhs'] for p in related),
            factorizations=sum(p['contextual_reserved_LU'] for p in related),
            child_PIDs=[p['pid'] for p in related],
            attribution='stage/time-span association, not an explicit receipt foreign key',
            no_additional_charge=True)

    pending=defaultdict(deque);matches=[];orphans=[];custom_unknown=[];phase_reconciliations=[]
    for n,r in rows:
        event=r.get('event')
        if event=='reserve':
            role=_role(r.get('label'))
            if r.get('factorizations',0):pending[(_context(r),'factor')].append([n,r,int(r['factorizations'])])
            if r.get('rhs',0):pending[(_context(r),role)].append([n,r,int(r['rhs'])])
        elif event in ('backend_factor','backend_solve','base_LU_action_residual','backend_residual_failure'):
            role=('factor' if event=='backend_factor' else
                  r.get('label') if event=='base_LU_action_residual' else
                  'forward_sources' if event=='backend_residual_failure' else _role(r.get('label'),True))
            count=1 if role=='factor' else r.get('rhs',r.get('RHS'))
            if not isinstance(count,int) or count<1:continue
            context=_context(r)
            # Green's recorded base-LU completion omits phase, whereas its
            # preceding reservation includes it. Reconcile only this schema,
            # and only when a unique earlier unconsumed context matches all
            # of stage/PID/scope/label. No reservation or charge is invented.
            # Multiple possible phases remain unmatched and fail below.
            if event=='base_LU_action_residual' and context[2] is None:
                candidates=[key for key,q in pending.items()
                            if key[1]==role and key[0][0]==context[0]
                            and key[0][1]==context[1] and key[0][3]==context[3]
                            and any(item[2]>0 for item in q)]
                if len(candidates)==1:
                    context=candidates[0][0]
                    phase_reconciliations.append(dict(solver_line=n,
                        recorded_context=list(_context(r)),matched_context=list(context),
                        role=role,rule='unique earlier unconsumed stage/PID/scope/label context'))
            queue=pending.get((context,role),deque())
            remaining=count;linked=[]
            for item in list(queue):
                reserve_line,reserve,left=item
                f0=reserve.get('freq');f1=r.get('freq')
                if f0 is not None and f1 is not None and f0!=f1:continue
                amount=min(left,remaining)
                if amount:
                    item[2]-=amount;remaining-=amount
                    linked.append(dict(reserve_line=reserve_line,reserve_id=reserve.get('id'),units=amount))
                if remaining==0:break
            while queue and queue[0][2]==0:queue.popleft()
            if remaining:
                orphans.append(dict(solver_line=n,event=event,role=role,unmatched_units=remaining,context=list(_context(r))))
            matches.append(dict(solver_line=n,event=event,role=role,matched_reserves=linked))
    if orphans:finding('FAIL','SOLVER_EVENT_WITHOUT_PRIOR_RESERVE','A logged solver completion lacks matching earlier reserved units',events=orphans)
    unmatched=[]
    for (context,role),queue in pending.items():
        for n,r,left in queue:
            if left<=0:continue
            process=process_info.get((r.get('stage'),r.get('pid')),{'association':'UNSETTLED'})
            if role not in ('factor','forward_sources','general_receiver_transpose','pixel_receiver_transpose','material_jvp','material_vjp'):
                custom_unknown.append(dict(line=n,label=role,units=left));continue
            state=process['association']
            classification=('CHARGED_UNCONFIRMED_FAILED_PROCESS' if state=='SETTLED_FAILED' else
                            'UNVERIFIED_ACTION_IN_SUCCESSFUL_PROCESS' if state=='SETTLED_SUCCESS' else
                            'UNSETTLED_ACTION_REVIEW_REQUIRED' if final else 'IN_PROGRESS_RESERVED_ACTION')
            unmatched.append(dict(line=n,id=r.get('id'),role=role,units=left,stage=r.get('stage'),
                                  pid=r.get('pid'),classification=classification))
    if unmatched:
        finding('UNVERIFIED' if final or any(r['classification']=='UNVERIFIED_ACTION_IN_SUCCESSFUL_PROCESS' for r in unmatched) else 'IN_PROGRESS',
                'RESERVED_ACTIONS_WITHOUT_COMPLETION_RECEIPT','Charges are retained; solver completion is not assumed',actions=unmatched)
    if custom_unknown:finding('UNVERIFIED','CUSTOM_ACTION_SCHEMA_UNMATCHED','Custom reservations lack a recognized completion schema',actions=custom_unknown)

    # Inspect only indexed endpoint filenames. Never open .npz/.npy or private
    # generated-data metadata. Earlier unindexed Forward versions are not swept.
    endpoint_events=defaultdict(list)
    for n,r in rows:
        if r.get('event')=='material_frequency_attempt' and isinstance(r.get('endpoint'),str):
            endpoint_events[r['endpoint']].append((n,r))
    endpoints=[]
    for relative,events in sorted(endpoint_events.items()):
        path=(root/relative).resolve()
        if not path.is_relative_to(root/'results'):
            finding('FAIL','ENDPOINT_PATH_OUTSIDE_RESULTS','Indexed endpoint escapes results root',path=relative);continue
        first=events[0][1]
        state=process_info.get((first.get('stage'),first.get('pid')),{'association':'UNSETTLED'})['association']
        complete=(path/'COMPLETE.json').is_file()
        filenames=['PARAMETERS.json','FREQUENCIES.json','RESPONSE.npz','COMPLETE.json']
        missing=[name for name in filenames if not (path/name).is_file() or (path/name).stat().st_size==0]
        seen=[]
        for _,r in events:
            if r.get('frequency') not in seen:seen.append(r.get('frequency'))
        frequency_files=[f'f{i}.npz' for i in range(len(seen))]
        missing_frequency=[name for name in frequency_files if not (path/name).is_file() or (path/name).stat().st_size==0]
        requested=any(r.get('jacobian_requested') is True for _,r in events)
        missing_j=[f'J_f{i}.npy' for i in range(len(seen)) if not (path/f'J_f{i}.npy').is_file()] if requested and complete else []
        upgrade=(path/'J_UPGRADE.json').exists()
        upgrade_missing=upgrade and not (path/'RESPONSE_WITH_J.npz').is_file()
        if complete and not missing and not missing_frequency and not missing_j and not upgrade_missing:
            classification='DECLARED_COMPLETE_FILES_PRESENT'
        elif state=='SETTLED_FAILED':classification='PRESERVED_PARTIAL_FAILED_JOB'
        elif state=='SETTLED_SUCCESS':classification='INCOMPLETE_FILES_IN_SUCCESSFUL_JOB'
        else:classification='UNSETTLED_ENDPOINT_REVIEW_REQUIRED' if final else 'IN_PROGRESS_ENDPOINT'
        entry=dict(endpoint=relative,stage=first.get('stage'),pid=first.get('pid'),process_state=state,
                   indexed_frequency_count=len(seen),declared_complete=complete,
                   missing_core_files=missing,missing_frequency_files=missing_frequency,
                   missing_requested_J_files=missing_j,upgrade_response_missing=upgrade_missing,
                   private_generated_data_metadata_excluded=_sensitive(Path(relative)),classification=classification)
        endpoints.append(entry)
        if classification=='INCOMPLETE_FILES_IN_SUCCESSFUL_JOB':
            finding('FAIL','SUCCESSFUL_JOB_MISSING_DECLARED_ENDPOINT_FILES','Successful wrapper has an incomplete indexed endpoint',**entry)
    # Bounded metadata sample, prioritizing non-complete and recent safe paths.
    safe=[e for e in endpoints if not e['private_generated_data_metadata_excluded']]
    safe.sort(key=lambda e:(e['classification']=='DECLARED_COMPLETE_FILES_PRESENT',-max(n for n,_ in endpoint_events[e['endpoint']])))
    sampled=[]
    for entry in safe[:metadata_samples]:
        path=root/entry['endpoint'];item=dict(endpoint=entry['endpoint'])
        for name,keys in [('PARAMETERS.json',('frequencies','phase','jacobian_requested','identity')),
                          ('COMPLETE.json',('frequency_count','jacobian','identity'))]:
            if (path/name).is_file():
                try:
                    metadata=_read(path/name)
                    item[name]={key:metadata.get(key) for key in keys}
                except (ValueError,OSError) as exc:item[name]=dict(unsettled_read_error=str(exc))
        parameter=item.get('PARAMETERS.json',{});completed=item.get('COMPLETE.json',{})
        frequencies=parameter.get('frequencies')
        if isinstance(frequencies,list) and completed.get('frequency_count') is not None:
            item['frequency_count_consistent']=len(frequencies)==completed['frequency_count']
            if not item['frequency_count_consistent']:
                finding('FAIL','ENDPOINT_FREQUENCY_METADATA_MISMATCH','Sampled endpoint counts differ',endpoint=entry['endpoint'])
        sampled.append(item)

    residuals=[]
    for n,r in rows:
        if r.get('event') in ('backend_solve','backend_state','backend_residual_failure','base_LU_action_residual'):
            value=r.get('residual')
            if _num(value):residuals.append(dict(line=n,value=float(value),event=r.get('event'),
                stage=r.get('stage'),phase=r.get('phase'),label=r.get('label')))
    maximum=max(residuals,key=lambda r:r['value']) if residuals else None
    residual_failures=[dict(line=n,stage=r.get('stage'),phase=r.get('phase'),residual=r.get('residual'),tolerance=r.get('tolerance'))
                       for n,r in rows if r.get('event')=='backend_residual_failure']
    # Memory evidence is taken only from recorded setup/factor events and
    # closed wrapper receipts. Live RSS/opaque workspace usage remains unknown.
    backend_limits={};factor_memory=[];memory_rejections=[]
    for n,r in rows:
        key=(_context(r),r.get('freq'))
        if r.get('event')=='backend_setup':
            backend_limits[key]=r
            if (_num(r.get('ndofs')) and _num(r.get('max_dofs')) and r['ndofs']>r['max_dofs']):
                finding('FAIL','COMPLETED_SETUP_DOF_CAP_MISMATCH','A completed backend setup reports DOFs above its own limit',line=n)
        elif r.get('event')=='backend_factor' and _num(r.get('actual_factor_storage_bytes')):
            limit=backend_limits.get(key,{})
            factor_memory.append(dict(line=n,stage=r.get('stage'),phase=r.get('phase'),
                actual_factor_storage_bytes=r['actual_factor_storage_bytes'],
                configured_max_factor_bytes=limit.get('max_factor_bytes')))
            if _num(limit.get('max_factor_bytes')) and r['actual_factor_storage_bytes']>limit['max_factor_bytes']:
                finding('FAIL','COMPLETED_FACTOR_CAP_MISMATCH','A completed factor exceeds its associated setup limit',line=n)
        elif r.get('event')=='backend_memory_reject':
            memory_rejections.append(dict(line=n,stage=r.get('stage'),phase=r.get('phase'),
                actual_factor_storage_bytes=r.get('actual_factor_storage_bytes'),
                estimated_factor_bytes=r.get('estimated_factor_bytes'),max_factor_bytes=r.get('max_factor_bytes')))
    nested_CPU=defaultdict(float);nested_wall=defaultdict(float)
    for _,r in rows:
        for key,value in r.items():
            if _num(value) and key in ('operation_cpu_seconds','solver_cpu_seconds','factor_cpu_seconds'):
                nested_CPU[key]+=value
        if r.get('event') in ('backend_solve','backend_factor','backend_setup') and _num(r.get('wall_seconds')):
            nested_wall[r['event']]+=r['wall_seconds']
        if r.get('event')=='backend_state' and isinstance(r.get('timings'),dict):
            if _num(r['timings'].get('state_seconds')):nested_wall['backend_state_overlapping_state_seconds']+=r['timings']['state_seconds']
    ledger_after=ledger_path.stat().st_size
    changing=ledger_after!=snapshot['bytes_read']
    unsettled_processes=[p for p in process_info.values() if p['association'] not in ('SETTLED_SUCCESS','SETTLED_FAILED')]
    live=not final and bool(unreceipted_jobs or unsettled_processes or changing or snapshot['trailing_partial_line'])
    unsealed=[stage for stage,row in stages.items() if not row['sealed']]
    if unreceipted_jobs:finding('UNVERIFIED' if final else 'IN_PROGRESS','JOBS_WITHOUT_FINAL_RECEIPT',
        'Existing wrapper directories are preserved; no final CPU charge/status yet',labels=unreceipted_jobs)
    if unsealed:finding('UNVERIFIED' if final else 'IN_PROGRESS','STAGES_NOT_SEALED',
        'Stage totals are a live/unsealed snapshot, not final audited expenditure',stages=unsealed)
    if changing:finding('IN_PROGRESS','LEDGER_ADVANCED_DURING_AUDIT','Filesystem and ledger snapshots are not globally atomic',bytes_at_end=ledger_after)
    limitations=[
        'Root ledger reserve rows, including failed attempts, are authoritative. Rejections are not successful reservations.',
        'Current parse/ID uniqueness and append-mode source evidence do not prove historical append-only immutability; no hash or trusted-prefix comparison performed.',
        'Wrapper receipts have no child PID/reserve-ID foreign key; child attribution is stage/time-span containment and may remain ambiguous.',
        'Endpoint existence/nonempty checks are not array decompression, numerical content validation or scientific correctness.',
        'Only ledger-indexed Forward endpoints are inventoried. Legacy unindexed endpoint persistence is NOT exhaustively verified.',
        'Private generated-data endpoint metadata, OFFLINE_TRUTH, GENERATION_PLAN, observations and job logs are never opened.',
        'Live/unreceipted process CPU is not yet settled in root ledger; no assertion that charged remaining CPU equals live headroom.',
        'Nested backend timers overlap, usually measure wall time, and must not be added to all-process CPU charges.',
        'Current source snippets document semantics; historical source snapshots are checked for existence, not equated to current code.'
    ]
    has_fail=any(r['level']=='FAIL' for r in findings)
    status='FAIL_INTEGRITY_FINDINGS' if has_fail else 'IN_PROGRESS' if live or (unsealed and not final) else 'FINAL_REPLAY_WITH_DECLARED_LIMITATIONS' if final else 'SNAPSHOT_NO_DETECTED_INTEGRITY_VIOLATION'
    result=dict(status=status,scientific_acceptance='NOT_EVALUATED',snapshot_utc=datetime.now(timezone.utc).isoformat(),
        root=str(root),mode='FINAL_REPLAY' if final else 'LIVE_SAFE_SNAPSHOT',ledger_snapshot=snapshot,
        event_counts=dict(event_counts),stage_totals=stages,B2_scope_totals=scopes,
        append_ID_checks=dict(unique_declared_IDs=len(ids),duplicates=duplicates,reserves_without_ID=missing_ids,
                              events_without_declared_ID=sum(not r.get('id') for _,r in rows)),
        reserve_solver_association=dict(matched_solver_events=len(matches),unmatched_solver_events=orphans,
            omitted_phase_unique_context_reconciliations=phase_reconciliations,
            unmatched_reserved_actions=unmatched,custom_unmatched_actions=custom_unknown,
            matches=matches,attribution='ordered stage/PID/phase/scope/role and compatible frequency; no explicit solver reservation foreign key'),
        job_receipts=receipt_records,process_job_associations=list(process_info.values()),
        memory_resource_evidence=dict(process_RSS_cap_bytes=authorization['max_process_RSS_bytes'],
            max_closed_receipt_reported_RSS_bytes=max((r['peak_rss_bytes'] for r in receipt_records if _num(r['peak_rss_bytes'])),default=None),
            max_closed_receipt_sampled_RSS_bytes=max((r['peak_sampled_rss_bytes'] for r in receipt_records if _num(r['peak_sampled_rss_bytes'])),default=None),
            max_completed_factor_storage=max(factor_memory,key=lambda r:r['actual_factor_storage_bytes']) if factor_memory else None,
            backend_memory_rejections=memory_rejections,
            live_process_RSS='NOT_READ_OR_VERIFIED; no process inspection or extra workload launched'),
        jobs_without_final_receipt=unreceipted_jobs,
        CPU_accounting=dict(root_all_process_charged_CPU=sum(r.get('charged_cpu_seconds',0) for _,r in rows),
            root_job_cost_CPU=sum(r.get('charged_cpu_seconds',0) for _,r in rows if r.get('event')=='job_cost'),
            root_non_job_CPU=sum(r.get('charged_cpu_seconds',0) for _,r in rows if r.get('event')!='job_cost'),
            receipt_CPU_not_added_again=True,nested_solver_CPU_informational=dict(nested_CPU),
            nested_solver_CPU_status='INFORMATIONAL_NONADDITIVE' if nested_CPU else 'NOT_SEPARATELY_MEASURED',
            nested_wall_seconds_informational_nonadditive=dict(nested_wall)),
        backward_residuals=dict(max_logged_RHS_relative_residual=maximum,
            definition='backend RHS-relative linear-system residual, not an independently computed normwise backward error',
            logged_count=len(residuals),explicit_backend_residual_failures=residual_failures),
        endpoint_persistence=dict(indexed_endpoints=len(endpoints),metadata_sample_count=len(sampled),
            metadata_samples=sampled,classification_counts=dict(Counter(e['classification'] for e in endpoints)),
            endpoints=endpoints,array_contents_opened=0),
        failures_and_recovery=dict(failed_wrapper_receipts=[r for r in receipt_records if r.get('returncode') not in (0,None) or str(r.get('status','')).startswith('FAILED')],
            budget_rejections=[dict(line=n,stage=r.get('stage'),phase=r.get('phase'),reason=r.get('reason'),requested_rhs=r.get('rhs'))
                               for n,r in rows if r.get('event') in ('budget_reject','scope_budget_stop')],
            partial_endpoints=[e for e in endpoints if e['classification']!='DECLARED_COMPLETE_FILES_PRESENT'],
            no_reset_or_recovery_performed=True),
        findings=findings,limitations=limitations,code_evidence=_source_evidence(root),
        audit_new_rhs=0,audit_new_LU=0,GT_files_read=0,hashes_computed=0,
        physics_source_or_budget_mutations=0,elapsed_wall_seconds=time.perf_counter()-start)
    return result


def report(result):
    lines=['# W execution integrity audit','',f"Status: **{result['status']}**. Snapshot: {result['snapshot_utc']}.",
           '', 'Accounting totals below come only from the root ledger. Receipt CPU is not charged a second time. This report does not judge scientific acceptance.',
           '', '| Stage | RHS / cap | LU / cap | settled all-process CPU / cap (s) | sealed |',
           '|---|---:|---:|---:|---|']
    for stage,row in result['stage_totals'].items():
        t,c=row['ledger_totals'],row['caps']
        lines.append(f"| {stage} | {t['rhs']} / {c['rhs']} | {t['factorizations']} / {c['factorizations']} | {t['cpu_seconds']:.3f} / {c['cpu_seconds']} | {row['sealed']} |")
    lines+=['','| B2 scope | RHS / cap | LU (no separate scope cap declared) |','|---|---:|---:|']
    for scope,row in result['B2_scope_totals'].items():
        lines.append(f"| {scope} | {row['ledger_totals']['rhs']} / {row['caps']['rhs']} | {row['ledger_totals']['factorizations']} |")
    residual=result['backward_residuals']['max_logged_RHS_relative_residual']
    lines+=['','## Recorded checks','',f"- Reservation ID duplicates: {len(result['append_ID_checks']['duplicates'])}.",
            f"- Solver events matched to prior reservations: {result['reserve_solver_association']['matched_solver_events']}.",
            f"- Solver events without prior reservations: {len(result['reserve_solver_association']['unmatched_solver_events'])}.",
            f"- Maximum logged RHS-relative solver residual: {residual['value']:.6g} at ledger line {residual['line']}." if residual else '- No solver residual recorded.',
            f"- Indexed endpoint states: {result['endpoint_persistence']['classification_counts']}.",
            f"- Safe metadata sample: {result['endpoint_persistence']['metadata_sample_count']}; arrays opened: 0.",
            f"- Root settled CPU: {result['CPU_accounting']['root_all_process_charged_CPU']:.6f}s; nested solver CPU: {result['CPU_accounting']['nested_solver_CPU_status']}.",
            '', '## Findings and pending settlement','']
    if not result['findings']:lines.append('- No detected integrity finding within the declared coverage.')
    for finding in result['findings']:
        lines.append(f"- **{finding['level']} / {finding['code']}**: {finding['message']}")
    lines+=['','## Failed runs remain paid','']
    for job in result['failures_and_recovery']['failed_wrapper_receipts']:
        lines.append(f"- `{job['label']}`: {job['status']}, CPU {job['charged_cpu_seconds']}s, contextual RHS {job['contextual_physical_charges']['rhs']} / LU {job['contextual_physical_charges']['factorizations']}, {job['association']}; ledger lines {job['ledger_lines']}.")
    lines+=['','## Concrete current-code evidence','']
    for item in result['code_evidence']:
        lines.append(f"- {item['purpose']}: `{item['path']}` ({item['status']}).")
        for snippet in item.get('snippets',[]):
            lines.append(f"  - Line {snippet['line']}: `{snippet['code']}`")
    lines+=['','## Coverage limits','']+[f'- {item}' for item in result['limitations']]
    lines+=['', 'Audit actions: zero RHS, zero LU, zero GT files, no hashes; no physics/configuration/ledger/job/array changes.',
            '', 'Replay with `python tools/audit_w_execution.py --final` after settlement. Live stages are not failed merely because they are unfinished.','']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    parser.add_argument('--final',action='store_true')
    parser.add_argument('--metadata-samples',type=int,default=24)
    args=parser.parse_args()
    result=audit(args.root,args.final,args.metadata_samples)
    out=Path(args.root)/'AUDIT/execution_integrity';out.mkdir(parents=True,exist_ok=True)
    _write(out/'RESULTS.json',result)
    (out/'REPORT.md').write_text(report(result))
    print(json.dumps(dict(status=result['status'],stage_totals={k:v['ledger_totals'] for k,v in result['stage_totals'].items()},
                         findings=[dict(level=f['level'],code=f['code']) for f in result['findings']],
                         outputs=[str(out/'RESULTS.json'),str(out/'REPORT.md')]),indent=2))


if __name__=='__main__':main()
