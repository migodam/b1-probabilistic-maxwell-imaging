#!/usr/bin/env python3
"""Finalize measured bank setup cost after every finite-bank wrapper terminates.

Reads snapshots only; never modifies or charges the authoritative ledger.
"""
import argparse
import json
import os
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parents[1]


def write_json(path,value):
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as stream:
            json.dump(value,stream,indent=2,allow_nan=False,sort_keys=True)
            stream.write('\n');stream.flush();os.fsync(stream.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def run(args):
    root=Path(args.root).resolve();bank=Path(args.bank or root/'results/e2e3/finite_bank').resolve()
    checkpoint=json.loads((bank/'CHECKPOINT.json').read_text())
    complete=json.loads((bank/'COMPLETE.json').read_text())
    if checkpoint.get('status')!='COMPLETE' or complete.get('status')!='COMPLETE' or int(checkpoint.get('completed',0))!=2950:
        raise RuntimeError('complete physical bank and terminal COMPLETE receipt required')
    ledger=json.loads((root/'results/budget_ledger.json').read_text());events=ledger['events']
    starts=[(i,e) for i,e in enumerate(events) if e.get('kind')=='reservation' and e.get('source')=='run_job' and str(e.get('label','')).startswith('finite-bank-')]
    if not starts:raise RuntimeError('no finite-bank wrapper starts found')
    jobs=[];covered=set();raw_receipts=[]
    for start_index,start in starts:
        label=start['label'];stage=start.get('stage')
        if stage!='E2E3':raise RuntimeError('finite-bank wrapper charged outside E2E3')
        next_start=next((i for i in range(start_index+1,len(events)) if events[i].get('kind')=='reservation' and events[i].get('source')=='run_job' and events[i].get('label')==label),len(events))
        endings=[(i,e) for i,e in enumerate(events[start_index+1:next_start],start_index+1) if e.get('kind')=='job_completion' and e.get('label')==label and e.get('stage')==stage]
        if len(endings)!=1:raise RuntimeError('finite-bank wrapper is not uniquely terminal: '+label)
        end_index,end=endings[0]
        covered.update(range(start_index,end_index+1))
        receipt_path=root/'results/jobs'/label/'RECEIPT.json'
        if not receipt_path.exists():raise RuntimeError('missing durable wrapper receipt: '+label)
        receipt=json.loads(receipt_path.read_text())
        if receipt.get('id')!=end.get('id') or receipt.get('exit_code')!=end.get('exit_code'):
            raise RuntimeError('wrapper receipt does not match authoritative completion: '+label)
        raw_receipts.append(dict(path=receipt_path.relative_to(root).as_posix(),receipt=receipt))
        jobs.append(dict(label=label,start_event_id=start.get('id'),completion_event_id=end.get('id'),
            exit_code=end.get('exit_code'),measured_process_cpu_seconds=end.get('process_cpu_seconds',0.),
            wall_seconds=end.get('wall_seconds',0.),children_maxrss=end.get('children_maxrss'),rss_unit=end.get('rss_unit','UNKNOWN'),
            rss_semantics=end.get('rss_semantics','NOT_RECORDED')))
    # Discover receipt-only finite-bank jobs, rather than silently omitting them.
    declared={row['label'] for row in jobs}
    for directory in (root/'results/jobs').glob('finite-bank-*'):
        if directory.is_dir() and directory.name not in declared:
            raise RuntimeError('finite-bank job directory has no authoritative matched start/completion: '+directory.name)
    physical=[event for i,event in enumerate(events) if i in covered and event.get('kind')=='reservation' and event.get('stage')=='E2E3']
    rhs=sum(e.get('rhs',0) for e in physical);fe=sum(e.get('factorizations',0) for e in physical)
    if rhs!=11800 or fe!=2950:
        raise RuntimeError(f'physical bank reservations differ from declared complete bank: rhs={rhs}, FE={fe}; no inferred or duplicated cost')
    rss_by_unit={}
    for row in jobs:
        if row['children_maxrss'] is not None:
            unit=row['rss_unit'];rss_by_unit[unit]=max(rss_by_unit.get(unit,0),row['children_maxrss'])
    output=dict(schema='FINITE_BANK_FINAL_JOB_COST_V1',complete=True,prior_id='FINITE_PRIOR_ENUM4_V1',stage='E2E3',
        measured_process_cpu_seconds=sum(row['measured_process_cpu_seconds'] for row in jobs),
        summed_job_wall_seconds=sum(row['wall_seconds'] for row in jobs),
        physical_bank=dict(planned_rhs=11800,planned_factorizations=2950,actual_rhs=rhs,actual_factorizations=fe),
        jobs=jobs,failed_jobs_included=[row['label'] for row in jobs if row['exit_code']!=0],
        raw_wrapper_receipts=raw_receipts,max_child_rss_by_unit=rss_by_unit,
        cpu_source='Only authoritative finite-bank-* job_completion measured CPU, each event counted once',
        wall_source='Only those measured wrapper walls, each once; no solver timers or estimates added',
        original_in_child_receipts_preserved=['online/COST_RECEIPT.json','CHECKPOINT.json','COMPLETE.json'],
        in_child_CPU_warning='Raw bank receipts omit their own in-flight wrapper CPU; use this terminal setup receipt for attribution',
        no_additional_ledger_charge=True,whole_campaign_cost_authority='results/budget_ledger.json',
        finalizer_own_cost='Charged separately by its E2E3 wrapper; not part of bank setup and never manually charged again')
    write_json(bank/'online/JOB_COST_FINAL.json',output)
    print(json.dumps(dict(complete=True,path=str(bank/'online/JOB_COST_FINAL.json'),measured_process_cpu_seconds=output['measured_process_cpu_seconds'],rhs=rhs,factorizations=fe)))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',default=str(ROOT));parser.add_argument('--bank')
    run(parser.parse_args(argv))


if __name__=='__main__':main()
