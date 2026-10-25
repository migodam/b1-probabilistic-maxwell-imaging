#!/usr/bin/env python3
"""Rebuild evidence exports only. Never mutates budgets or scientific gates."""
import argparse
import csv
import json
import os
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parents[1]
EXCLUDED={'.git','.aws','.codex','.agents','.venv','.deps','node_modules','__pycache__'}
EXPORTS={'COST_SUMMARY.json','EXPORT_STATUS.json','A0_ANALYSIS_PROVENANCE.json','PHYSICS_LEDGER.jsonl','FAILURE_LEDGER.jsonl','PUBLICATION_ARRAY_MANIFEST.csv'}


def atomic_text(path,text):
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as stream:
            stream.write(text);stream.flush();os.fsync(stream.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def write_json(path,value):
    atomic_text(path,json.dumps(value,indent=2,allow_nan=False,sort_keys=True)+'\n')


def write_jsonl(path,values):
    atomic_text(path,''.join(json.dumps(value,allow_nan=False,separators=(',',':'))+'\n' for value in values))


def workspace_files(root):
    """Enumerate only the configured tree; never follow directory symlinks."""
    for directory,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=sorted(d for d in dirs if d not in EXCLUDED and not (Path(directory)/d).is_symlink())
        for name in sorted(files):
            path=Path(directory)/name
            if path.is_symlink():continue
            yield path


def cost_summary(ledger,receipts):
    caps=ledger['caps'];events=ledger['events'];stages={}
    for stage,cap in caps.items():
        selected=[e for e in events if e.get('stage')==stage]
        reservations=[e for e in selected if e.get('kind')=='reservation']
        jobs=[e for e in selected if e.get('kind')=='job_completion']
        # Reservations are the sole physical charges. Observed counts ignored.
        used=dict(rhs=sum(e.get('rhs',0) for e in reservations),
            factorizations=sum(e.get('factorizations',0) for e in reservations),
            process_cpu_seconds=sum(e.get('process_cpu_seconds',0.) for e in selected))
        estimated=[e for e in selected if e.get('estimated_allocation')]
        correction=[e for e in selected if e.get('kind')=='accounting_correction']
        measured=sum(e.get('process_cpu_seconds',0.) for e in jobs)
        prep=sum(e.get('process_cpu_seconds',0.) for e in estimated)
        corrected=sum(e.get('process_cpu_seconds',0.) for e in correction if e not in estimated)
        other=used['process_cpu_seconds']-measured-prep-corrected
        rss=[];seen=set()
        for source,event in [('ledger',e) for e in jobs]+receipts:
            if event.get('stage')!=stage or 'children_maxrss' not in event:continue
            identity=(event.get('id'),event.get('stage'),event.get('label'))
            if identity in seen:continue
            seen.add(identity)
            rss.append(dict(source=source,event_id=event.get('id'),label=event.get('label'),
                value=event['children_maxrss'],unit=event.get('rss_unit','UNKNOWN'),
                semantics=event.get('rss_semantics','NOT_RECORDED')))
        by_unit={}
        for row in rss:by_unit[row['unit']]=max(by_unit.get(row['unit'],0),row['value'])
        stages[stage]=dict(caps=cap,used=used,
            remaining={key:max(0,cap[key]-used[key]) for key in cap},
            cap_exceeded={key:used[key]>cap[key] for key in cap},
            cpu=dict(measured_job_completion_seconds=measured,
                estimated_preparation_seconds=prep,accounting_correction_seconds=corrected,
                other_charged_seconds=other,
                estimated_events=[dict(id=e.get('id'),label=e.get('label'),cpu=e.get('process_cpu_seconds'),scope=e.get('scope')) for e in estimated],
                correction_events=[dict(id=e.get('id'),label=e.get('label'),cpu=e.get('process_cpu_seconds'),measurement_status=e.get('measurement_status'),reason=e.get('reason')) for e in correction],
                running_jobs_cpu='NOT_YET_MEASURED_OR_CHARGED_UNTIL_COMPLETION'),
            summed_job_wall_seconds=sum(e.get('wall_seconds',0.) for e in jobs),
            job_completion_events=len(jobs),failed_job_count=sum(e.get('exit_code',0)!=0 for e in jobs),
            max_child_rss_by_unit=by_unit,rss_receipts=rss,
            rss_rule='Maximum per recorded unit; never sum wrapper lifetime high-water marks',
            wall_rule='Sum job_completion wall once; no solver timers, receipts or estimated wall added')
    return dict(schema='WSLAM_COST_EXPORT_V1',stages=stages,event_count=len(events),
        physical_charge_rule='Only kind=reservation rhs/factorizations; observed_rhs/observed_factorizations never charged',
        shared_cost_rule='Physical work counted once in its stage; standalone policy attribution is not extra actual work',
        measured_cpu_is_not_estimated_cpu=True,budget_caps_unchanged=True)


def brief(node):
    """Keep exact scalar diagnostics and identities; point to full source arrays."""
    if not isinstance(node,dict):return node
    result={}
    for key,value in node.items():
        if value is None or isinstance(value,(str,int,float,bool)):
            result[key]=value
        elif isinstance(value,list) and len(value)<=12 and all(isinstance(x,(str,int,float,bool,type(None))) for x in value):
            result[key]=value
        elif key in ('summary','counts','stage_cost_diff','stage_cost_diffs','budget_start','budget_end','scope','caps','gate') and isinstance(value,dict):
            result[key]=brief(value)
    return result


def walk_nodes(value,pointer=''):
    if isinstance(value,dict):
        yield pointer,value
        for key,child in value.items():
            yield from walk_nodes(child,pointer+'/'+str(key).replace('~','~0').replace('/','~1'))
    elif isinstance(value,list):
        for index,child in enumerate(value):yield from walk_nodes(child,pointer+'/'+str(index))


def failure_category(node,path):
    scalar={str(k).lower():v for k,v in node.items() if isinstance(v,(str,int,float,bool)) or v is None}
    text=' '.join(str(v) for v in scalar.values()).lower()
    if 'broadcast' in text or 'unexpected source layout' in text or 'shape error' in text:
        return 'INVALID_ARRAY_SHAPE_OR_LAYOUT_NOT_NEGATIVE_EVIDENCE'
    if any('preflight' in str(k).lower() or 'preflight' in str(v).lower() for k,v in scalar.items()):
        return 'PREFLIGHT_EVIDENCE'
    statuses=[str(v).upper() for k,v in scalar.items() if any(t in k for t in ('status','decision','gate'))]
    if any('MAX_NFEV' in value or 'NONCONVERGENCE' in value for value in statuses) or 'maximum number of function evaluations' in text:
        return 'NONCONVERGENCE_RETAINED'
    if scalar.get('success') is False or scalar.get('selected_fit_success') is False:
        return 'NONCONVERGENCE_OR_FAILED_SUCCESS_FLAG'
    if any(any(t in value for t in ('FAIL','STOPPED','ERROR')) for value in statuses):
        return 'FAILED_OR_STOPPED_RESULT'
    if any(k.startswith('error') and bool(v) for k,v in scalar.items()):return 'RECORDED_ERROR'
    if any(any(t in value for t in ('NOT_RUN','NOT_ESTABLISHED','UNKNOWN','INCONCLUSIVE','HOLD','NOT_EVALUABLE')) for value in statuses):
        return 'UNRESOLVED_OR_NOT_RUN_RETAINED'
    if 'smc' in path.lower() or 'smc' in text:
        if scalar.get('pass') is False or scalar.get('accurate') is False or scalar.get('within_tolerance') is False:
            return 'SMC_ACCURACY_FAILURE_RECORDED'
    if any(v is True for k,v in scalar.items() if 'exceed' in k or k.endswith('_failed')):
        return 'QUALIFICATION_THRESHOLD_EXCEEDED'
    return None


def analysis_provenance(root):
    original=root/'results/a0/analysis/RESULTS.json';repaired=root/'results/a0/analysis_repaired/RESULTS.json'
    record=dict(original_path=str(original.relative_to(root)),repaired_path=str(repaired.relative_to(root)),
        original_exists=original.exists(),repaired_exists=repaired.exists(),
        original_zero_witness_count_is_not_valid_negative_evidence_when_rows_invalid=True)
    if original.exists():
        data=json.loads(original.read_text());rows=data.get('rows',[])
        invalid=[]
        for i,row in enumerate(rows):
            text=' '.join(str(row.get(key,'')) for key in ('noiseless_pair_error','gradient_unavailable_reason','noiseless_pair_status')).lower()
            if 'broadcast' in text or 'unexpected source layout' in text or 'not_evaluable' in text:invalid.append(i)
        record.update(original_summary=data.get('summary'),original_row_count=len(rows),
            original_invalid_shape_or_layout_rows=len(invalid),original_invalid_row_indices=invalid,
            original_conclusion_scope='Invalid rows cannot establish zero physical ambiguity witnesses; original analysis preserved, never silently replaced.')
    if repaired.exists():
        data=json.loads(repaired.read_text());record['repaired_summary']=data.get('summary')
        record['repaired_scope']='Separate repaired analysis; its reported counts do not erase original implementation failures or independently validate witnesses.'
    else:record['repaired_scope']='NOT_FOUND_AT_EXPORT; no repaired outcome inferred'
    return record


def run(args):
    root=Path(args.root).resolve();ledger_path=root/'results/budget_ledger.json'
    # Atomic ledger replacement makes this read one complete snapshot.
    ledger=json.loads(ledger_path.read_text())
    if ledger.get('schema')!='wslam-budget-v1':raise ValueError('unexpected authoritative ledger schema')
    events=ledger['events'];files=list(workspace_files(root))
    failures=[];receipts=[];source_status=[];read_errors=[]
    # Raw events retain original order and IDs, including zero-PDE failed launches.
    write_jsonl(root/'PHYSICS_LEDGER.jsonl',events)
    for index,event in enumerate(events):
        if event.get('kind')=='job_completion' and event.get('exit_code',0)!=0:
            starts=[i for i in range(index) if events[i].get('kind')=='reservation' and events[i].get('source')=='run_job' and events[i].get('stage')==event.get('stage') and events[i].get('label')==event.get('label')]
            interval=events[starts[-1]:index+1] if starts else None
            charged=None if interval is None else dict(rhs=sum(e.get('rhs',0) for e in interval if e.get('kind')=='reservation'),factorizations=sum(e.get('factorizations',0) for e in interval if e.get('kind')=='reservation'))
            failures.append(dict(category='FAILED_JOB_COMPLETION',source='results/budget_ledger.json',event=event,physical_reservations_during_job=charged,interval_scope='Serialized wrapper start to completion; unknown if matching start absent'))
        elif 'preflight' in json.dumps(event).lower():
            failures.append(dict(category='PREFLIGHT_EVENT',source='results/budget_ledger.json',event=event))
    for path in files:
        rel=path.relative_to(root).as_posix()
        if path.suffix!='.json' or path.name in EXPORTS or path==ledger_path:continue
        # Scientific results and root preregistration/report documents only.
        if not (rel.startswith('results/') or len(path.relative_to(root).parts)==1):continue
        try:data=json.loads(path.read_text())
        except (ValueError,OSError) as exc:
            read_errors.append(dict(source=rel,error=str(exc)));continue
        if path.name=='RECEIPT.json' and isinstance(data,dict):receipts.append((rel,data))
        if isinstance(data,dict) and (path.name.lower()=='gate.json' or any(k in data for k in ('status','decision','overall_gate','scientific_decision'))):
            source_status.append(dict(source=rel,reported=brief(data),scientific_gate_modified=False))
        for pointer,node in walk_nodes(data):
            category=failure_category(node,rel)
            if category:failures.append(dict(category=category,source=rel,json_pointer=pointer,evidence=brief(node)))
        if rel=='results/e0/MESH_RESOURCE_PLAN.json':
            failures.append(dict(category='MEMORY_REFINEMENT_PLAN_EVIDENCE',source=rel,evidence=data,
                rule='Plan/DOF feasibility is distinct from executed physical qualification'))
    # Preserve the original failing qualification history regardless of later PASS.
    for rel in ('results/e0/qualification/RESULTS.json','results/e0/refinement01_executable/RESULTS.json'):
        path=root/rel
        if path.exists():
            data=json.loads(path.read_text())
            failures.append(dict(category='E0_ORIGINAL_QUALIFICATION_AND_REFINEMENT_HISTORY',source=rel,evidence=data,
                scope='Original h=.200/.185 qualification attempts and subsequent refinement are retained independently of final gate; inspect recorded per-row h and outcomes.'))
    # Preserve compiler provenance using the actual failed job log when present.
    compiler_log=root/'results/jobs/finite-bank-001/OUTPUT.log'
    if compiler_log.exists():
        lines=compiler_log.read_text(errors='replace').splitlines()
        matches=[(i,line) for i,line in enumerate(lines) if 'keyword argument repeated' in line or 'SyntaxError' in line]
        if matches:
            context=[]
            for i,_ in matches:context.extend(lines[max(0,i-2):min(len(lines),i+2)])
            failures.append(dict(category='COMPILATION_FAILURE_BEFORE_PDE',source=compiler_log.relative_to(root).as_posix(),
                evidence=context,qualification='AST parsing alone did not establish compilability; failed launch retained, not silently replaced',
                zero_PDE_evidence_rule='Consult original failed job receipt and reservation sequence; no inferred physical solves'))
    provenance=analysis_provenance(root)
    failures.append(dict(category='A0_ORIGINAL_AND_REPAIRED_ANALYSIS_PROVENANCE',source='results/a0/analysis/RESULTS.json',evidence=provenance))
    write_json(root/'A0_ANALYSIS_PROVENANCE.json',provenance)
    write_jsonl(root/'FAILURE_LEDGER.jsonl',failures)
    write_json(root/'COST_SUMMARY.json',cost_summary(ledger,receipts))
    # Enumerate scientific artifact paths and byte sizes only; no hash or content scan.
    manifest=[]
    for path in files:
        if path.suffix.lower() not in ('.npz','.png','.csv') or path.name=='PUBLICATION_ARRAY_MANIFEST.csv':continue
        stat=path.stat();manifest.append((path.relative_to(root).as_posix(),stat.st_size,path.suffix.lower()[1:]))
    import io
    buffer=io.StringIO();writer=csv.writer(buffer,lineterminator='\n');writer.writerow(('source_path','size_bytes','artifact_type'))
    writer.writerows(sorted(manifest));atomic_text(root/'PUBLICATION_ARRAY_MANIFEST.csv',buffer.getvalue())
    write_json(root/'EXPORT_STATUS.json',dict(status='EVIDENCE_EXPORT_REBUILT',ledger_source='results/budget_ledger.json',
        ledger_events=len(events),last_original_event_id=events[-1].get('id') if events else None,
        failure_evidence_entries=len(failures),artifact_manifest_entries=len(manifest),
        scientific_gates_modified=False,budgets_reset=False,ledger_modified=False,
        original_event_order_and_ids_preserved=True,source_status_reports=source_status,source_read_errors=read_errors,
        scope='Deterministic exports from this scientific workspace snapshot; rerun after reports and publication tail.',
        exclusions=sorted(EXCLUDED),directory_symlinks_followed=False,
        credentials_scanned=False,outside_workspace_scanned=False,
        limitations='Export presence does not establish scientific acceptance. Running jobs may have reserved physical work but no completion CPU receipt yet.'))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',default=str(ROOT))
    run(parser.parse_args(argv))


if __name__=='__main__':main()
