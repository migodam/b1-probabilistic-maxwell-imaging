"""Aggregate recorded cost without truth/scorer reads or CPU double counting.

Authorities: outer run_suite RUSAGE_CHILDREN first, future Core queue child
rusage second, inner CLI measured CPU otherwise. Main and explicitly separate
deployment physics ledgers both count. Missing CPU is never replaced by wall.
"""
from pathlib import Path
from collections import defaultdict
import argparse,csv,datetime,json,os
ROOT=Path(__file__).resolve().parents[1]
def read_jsonl(path):
    if not path.exists():return []
    return [json.loads(x) for x in path.read_text().splitlines() if x]
def first_guard_pid(out):
    p=out/'FILE_ACCESS.jsonl'
    if not p.exists():return None
    with p.open() as f:
        try:return json.loads(f.readline()).get('pid')
        except (ValueError,TypeError):return None
def process_alive(pid):
    try:os.kill(pid,0);return True
    except (ProcessLookupError,PermissionError):return False
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'results/COST_AGGREGATION.json');a=ap.parse_args()
    ledgers=[ROOT/'PHYSICS_LEDGER.jsonl']+sorted(p for p in (ROOT/'results').rglob('PHYSICS_LEDGER.jsonl') if 'deployment_accounting' in p.parts)
    bypid=defaultdict(lambda:dict(rhs=0,full_fe=0,first_unix=None,last_unix=None,first_event_unix=None,last_event_unix=None,ledger_paths=set(),phases=set(),outputs=set(),failures=[]))
    completions={};allrows=[]
    for ledger in ledgers:
        for source_line,row in enumerate(read_jsonl(ledger),1):
            row['_source_line']=source_line;row['_ledger']=str(ledger.relative_to(ROOT));allrows.append(row);pid=row.get('pid')
            if pid is None:continue
            d=bypid[pid];t=row.get('time_unix')
            if t is not None:d['first_event_unix']=t if d['first_event_unix'] is None else min(t,d['first_event_unix']);d['last_event_unix']=t if d['last_event_unix'] is None else max(t,d['last_event_unix'])
            d['ledger_paths'].add(str(ledger.relative_to(ROOT)))
            if row.get('kind')=='reservation':
                d['rhs']+=row.get('rhs',0);d['full_fe']+=row.get('factorizations',0);d['phases'].add(row.get('phase'))
                if t is not None:d['first_unix']=t if d['first_unix'] is None else min(t,d['first_unix']);d['last_unix']=t if d['last_unix'] is None else max(t,d['last_unix'])
            if row.get('label')=='job_start' and row.get('output'):d['outputs'].add(row['output'])
            if row.get('kind')=='failure':d['failures'].append(dict(label=row.get('label'),error=row.get('error'),time_unix=t))
            if row.get('label')=='job_completion' and row.get('process_cpu_seconds') is not None:completions[pid]=row
    # Map durable output paths, including renamed invalid/development attempts.
    receipts={}
    delegated_receipts=sorted((ROOT/'research/delegated').rglob('RECEIPT.json')) if (ROOT/'research/delegated').exists() else []
    exploration_receipts=(sorted((ROOT/'exploration').rglob('runtime_cost.json'))+sorted((ROOT/'exploration').rglob('RECEIPT.json'))) if (ROOT/'exploration').exists() else []
    for p in sorted((ROOT/'results').rglob('runtime_cost.json'))+sorted((ROOT/'results').rglob('RECEIPT.json'))+delegated_receipts+exploration_receipts:
        if '/jobs/' in str(p):continue
        try:r=json.loads(p.read_text())
        except (ValueError,OSError):continue
        pid=r.get('pid') or first_guard_pid(p.parent)
        if pid is None and r.get('process_cpu_seconds') is not None:
            matches=[q for q,x in completions.items() if abs(x['process_cpu_seconds']-r['process_cpu_seconds'])<1e-8]
            if len(matches)==1:pid=matches[0]
        if pid is None:continue
        bypid[pid]['outputs'].add(str(p.parent));receipts[pid]=(p,r)
    candidates=[]
    for pid,row in completions.items():candidates.append(dict(key=f'inner:{pid}',priority=1,coverage={pid},cpu=row['process_cpu_seconds'],authority='inner CLI process CPU over recorded job interval',source=row['_ledger']+':'+str(row['_source_line']),label=str(pid)))
    for pid,(p,r) in receipts.items():
        if pid not in completions and r.get('process_cpu_seconds') is not None:candidates.append(dict(key=f'inner_receipt:{pid}',priority=1,coverage={pid},cpu=r['process_cpu_seconds'],authority='inner CLI receipt',source=str(p.relative_to(ROOT)),label=str(pid)))
    queue_rows=[];active_outputs=set()
    for q in sorted((ROOT/'results').glob('CORE_QUEUE_r*.jsonl')):
        for queue_line,row in enumerate(read_jsonl(q),1):
            out=ROOT/'results/online'/row['case']/row['method']
            if row.get('event')=='job_start':active_outputs.add(str(out));continue
            if row.get('event')!='job_finish':continue
            active_outputs.discard(str(out))
            out=ROOT/'results/online'/row['case']/row['method'];pid=first_guard_pid(out);row['_pid']=pid;row['_source']=str(q.relative_to(ROOT));queue_rows.append(row)
            if pid is not None:bypid[pid]['outputs'].add(str(out))
            if row.get('child_process_cpu_seconds') is not None:candidates.append(dict(key=f'queue:{pid}:{q.name}',priority=2,coverage={pid} if pid else set(),cpu=row['child_process_cpu_seconds'],authority='Core queue per-child RUSAGE_CHILDREN delta',source=str(q.relative_to(ROOT))+':'+str(queue_line),label=row['case']+'/'+row['method']))
    outer=read_jsonl(ROOT/'JOB_LEDGER.jsonl')
    for i,row in enumerate(outer):
        cpu=row.get('process_cpu_seconds')
        if cpu is None:continue
        pid=row.get('pid');coverage={pid} if pid is not None else set();finish=row.get('finished_unix');start=finish-row.get('wall_seconds',0) if finish is not None else None
        if not coverage and start is not None:
            possible={q for q,d in bypid.items() if d['first_event_unix'] is not None and d['first_event_unix']>=start-.05 and d['last_event_unix']<=finish+.05}
            if len(possible)==1:coverage=possible
        if 'run_core_queue.py' in ' '.join(row.get('command',[])) and start is not None:
            coverage|={q for q,r in completions.items() if start-.05<=r.get('time_unix',-1)<=finish+.05}
        candidates.append(dict(key=f'outer:{i}:{row.get("label")}',priority=3,coverage=coverage,cpu=cpu,authority=row.get('cost_authority','run_suite child rusage'),source='JOB_LEDGER.jsonl:'+str(i+1),label=row.get('label'),coverage_inference='explicit child PID' if pid is not None else ('unique recorded process lifetime inside wrapper interval' if coverage else 'no physics PID association'),wall_seconds=row.get('wall_seconds'),exit_code=row.get('exit_code')))
    # Highest-authority covering record wins. Inner values remain diagnostic.
    chosen=[];covered=set();suppressed=[]
    for c in sorted(candidates,key=lambda r:-r['priority']):
        if c['coverage'] and c['coverage']&covered:suppressed.append(c);continue
        chosen.append(c);covered|=c['coverage']
    attempts=[]
    for pid,d in sorted(bypid.items()):
        if not(d['rhs'] or d['full_fe']):continue
        rec=receipts.get(pid,(None,{}))[1];q=next((r for r in queue_rows if r['_pid']==pid),{})
        cp=next((r for r in chosen if pid in r['coverage']),None)
        observed_rss=rec.get('peak_rss_bytes',rec.get('polled_peak_RSS_bytes'));rss_authority='inner process self ru_maxrss' if rec.get('peak_rss_bytes') is not None else None
        if observed_rss is None and cp and cp['source'].startswith('JOB_LEDGER.jsonl:'):
            outerrow=outer[int(cp['source'].split(':')[-1])-1];observed_rss=outerrow.get('polled_peak_RSS_bytes');rss_authority='outer run_suite observed polled peak; not full self ru_maxrss'
        attempts.append(dict(pid=pid,rhs=d['rhs'],full_fe=d['full_fe'],ledger_paths=sorted(d['ledger_paths']),phases=sorted(str(x) for x in d['phases']),outputs=sorted(d['outputs']),first_unix=d['first_unix'],last_unix=d['last_unix'],cpu_seconds=cp['cpu'] if cp else None,cpu_status='RECORDED' if cp else ('RUNNING' if any(p in active_outputs for p in d['outputs']) or (not q and not rec.get('scientific_sealed',False) and process_alive(pid)) else 'NOT_RECORDED'),cpu_authority=cp['authority'] if cp else None,cpu_source=cp['source'] if cp else None,wall_seconds=rec.get('wall_seconds',q.get('wall_seconds')),peak_rss_bytes=observed_rss,peak_rss_authority=rss_authority,exit_code=q.get('returncode',rec.get('exit_code')),scientific_sealed=q.get('sealed',rec.get('scientific_sealed')),failures=d['failures']))
    clean=lambda rows:[{k:sorted(v) if isinstance(v,set) else v for k,v in r.items()} for r in rows]
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),physics_ledgers=[str(p.relative_to(ROOT)) for p in ledgers],physics_totals=dict(rhs=sum(d['rhs'] for d in bypid.values()),full_fe=sum(d['full_fe'] for d in bypid.values())),known_recorded_cpu_subtotal_seconds=sum(r['cpu'] for r in chosen),cpu_total_status='KNOWN_SUBTOTAL_ONLY' if any(x['cpu_status']=='NOT_RECORDED' for x in attempts) else 'RECORDED_JOBS_SUBTOTAL',missing_physics_cpu_attempts=[dict(pid=x['pid'],outputs=x['outputs'],rhs=x['rhs'],full_fe=x['full_fe']) for x in attempts if x['cpu_status']=='NOT_RECORDED'],ongoing_physics_attempts=[dict(pid=x['pid'],outputs=x['outputs']) for x in attempts if x['cpu_status']=='RUNNING'],cpu_counted_once=clean(chosen),overlapping_cpu_records_not_added=clean(suppressed),attempts=attempts,notes=['No CPU value is inferred from wall time.','Main and output-local deployment reservation ledgers are both counted.','Outer run_suite child rusage takes precedence over Core queue child rusage and inner CLI CPU.','Score, plot and I/O jobs appear when recorded in JOB_LEDGER; unwrapped preparation/authoring/I/O CPU is not reconstructed.','The two normal_r0 post-seal receipt failures retain unknown CPU/RSS; exact reservation and queue wall costs remain.','This is a live aggregation for root review, not a scientific acceptance gate.'])
    exact_rss=[dict(pid=pid,bytes=r.get('peak_rss_bytes'),authority='process self ru_maxrss',source=str(p.relative_to(ROOT))) for pid,(p,r) in receipts.items() if r.get('peak_rss_bytes') is not None]
    polled_rss=[dict(pid=r.get('pid'),bytes=r['polled_peak_RSS_bytes'],authority='run_suite polled process RSS',source='JOB_LEDGER.jsonl:'+str(i+1)) for i,r in enumerate(outer) if r.get('polled_peak_RSS_bytes') is not None]
    known_rss=exact_rss+polled_rss;maximum=max(known_rss,key=lambda r:r['bytes']) if known_rss else None
    digit_pids={pid for pid,d in bypid.items() if 'digit' in (' '.join(d['outputs'])+' '+' '.join(str(x) for x in d['phases'])).lower()}
    result['physics_totals_by_scope']={scope:dict(rhs=sum(d['rhs'] for pid,d in bypid.items() if (pid in digit_pids)==is_digit),full_fe=sum(d['full_fe'] for pid,d in bypid.items() if (pid in digit_pids)==is_digit)) for scope,is_digit in [('digits_exploration',True),('other_R3_campaign',False)]}
    result['known_cpu_subtotal_by_scope_seconds']={scope:sum(row['cpu'] for row in chosen if (bool(row['coverage']&digit_pids) or 'digit' in (str(row.get('label',''))+' '+row['source']).lower())==is_digit) for scope,is_digit in [('digits_exploration',True),('other_R3_campaign',False)]}
    result['rss_summary']=dict(maximum_known_recorded_peak=maximum,maximum_exact_self_peak_bytes=max((r['bytes'] for r in exact_rss),default=None),maximum_polled_observed_peak_bytes=max((r['bytes'] for r in polled_rss),default=None),missing_physics_rss_attempts=[dict(pid=r['pid'],outputs=r['outputs']) for r in attempts if r['peak_rss_bytes'] is None and r['cpu_status']!='RUNNING'],scope='Known recorded peaks only; two lost normal_r0 RSS receipts remain unknown. Guards were enforced, not used as measurements.')
    result['ignored_mirror_physics_ledgers']=[str(p.relative_to(ROOT)) for p in (ROOT/'exploration').rglob('PHYSICS_LEDGER.jsonl')] if (ROOT/'exploration').exists() else []
    result['physics_ledger_policy']='Only original R3 main ledger plus output-local deployment accounting ledgers are authoritative. Exploration mirrors are excluded.'
    result['notes'].append('Legacy assembly circle_area is auxiliary pi*rx^2, not ellipse area; actual target_area/polygon_area quadrature uses rx*ry correctly. Frozen PDE source unchanged and no redo.')
    normalized=[dict(kind='normalized_cost_view',created_utc=result['created_utc'],original_logs_unchanged=True,source_ledgers=result['physics_ledgers'],cpu_policy='deduplicated authority: outer child rusage > queue child rusage > inner CLI',known_cpu_subtotal_seconds=result['known_recorded_cpu_subtotal_seconds'],rss_summary=result['rss_summary'],legacy_metadata_note=result['notes'][-1])]
    for row in allrows:
        if row.get('kind')=='reservation':normalized.append(dict(kind='physics_reservation',source_ref=row['_ledger']+':'+str(row['_source_line']),pid=row.get('pid'),time_unix=row.get('time_unix'),phase=row.get('phase'),label=row.get('label'),rhs=row.get('rhs',0),factorizations=row.get('factorizations',0),scope='digits_exploration' if row.get('pid') in digit_pids else 'other_R3_campaign',reserved_before_action=True))
    for row in chosen:normalized.append(dict(kind='cpu_record',source_ref=row['source'],coverage_pids=sorted(row['coverage']),cpu_seconds=row['cpu'],authority=row['authority'],label=row.get('label'),coverage_inference=row.get('coverage_inference','direct PID or queue output audit'),counted_once=True))
    for row in attempts:
        if row['cpu_status']!='RECORDED':normalized.append(dict(kind='pending_cpu' if row['cpu_status']=='RUNNING' else 'missing_cpu',pid=row['pid'],outputs=row['outputs'],status=row['cpu_status'],cpu_seconds=None))
    for row in known_rss:normalized.append(dict(kind='rss_observation',source_ref=row['source'],pid=row.get('pid'),rss_bytes=row['bytes'],authority=row['authority']))
    ledger_path=ROOT/'COST_LEDGER.jsonl';temporary=ledger_path.with_suffix('.jsonl.tmp');temporary.write_text(''.join(json.dumps(row)+'\n' for row in normalized));temporary.replace(ledger_path)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');csvpath=a.output.with_name('COST_ATTEMPTS.csv');fields=['pid','rhs','full_fe','cpu_seconds','cpu_status','wall_seconds','peak_rss_bytes','exit_code','scientific_sealed','outputs','cpu_authority','cpu_source']
    with csvpath.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in attempts:w.writerow({k:json.dumps(r[k]) if isinstance(r[k],(list,dict)) else r[k] for k in fields})
    print(json.dumps(dict(physics_totals=result['physics_totals'],known_cpu_seconds=result['known_recorded_cpu_subtotal_seconds'],missing_cpu=len(result['missing_physics_cpu_attempts']),attempts=len(attempts),cpu_records_counted=len(chosen))))
if __name__=='__main__':main()
