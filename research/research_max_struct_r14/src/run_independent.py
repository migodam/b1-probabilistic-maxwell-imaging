"""Dispatch the frozen pair or historical C2 runner in this same CPU-budgeted process."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import runpy
import sys,time

BASELINES=('cie','dual_cie','global_dual','g19','dr','direct')
PAIR_METHODS=tuple('pair'+p for p in ('00','01','10','11'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--legacy-repo',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--task-id',required=True)
    parser.add_argument('--cpu',type=int)
    args=parser.parse_args();root=args.root.resolve()
    manifest=json.loads(args.manifest.read_text())
    matches=[t for t in manifest['tasks'] if t['task_id']==args.task_id]
    if len(matches)!=1:raise ValueError('Exactly one manifest task must match')
    task=matches[0];method=task.get('algorithm',task['method'])
    if task.get('metric_schedule') is not None:
        runner=root/'src/run_control.py';family='R14_FROZEN_STAGE_METRIC_SCHEDULE'
        schedules={'pair01_to11':['01','01','01','11'],'pair11_to01':['11','11','11','01']}
        if method not in schedules or task.get('pair')!=schedules[method][0] or task['metric_schedule']!=schedules[method]:
            raise ValueError('Unregistered metric schedule')
    elif method in PAIR_METHODS:
        runner=root/'src/run_experiment.py';family='R14_STATEFUL_DUALMETRIC'
        if task.get('pair')!=method.removeprefix('pair'):raise ValueError('Pair task fields disagree')
    elif method in BASELINES:
        runner=root/'runtime/equal_budget_C2_runner.py';family='R13_C2_NUMERICS_WITH_SHARED_TOTAL120_BUDGET'
    else:raise ValueError('Unregistered independent method: '+str(method))
    if not runner.is_file():raise FileNotFoundError('Required frozen runner absent: '+str(runner))
    for key in ('campaign','task_id'):
        value=task[key]
        if not value or Path(value).name!=value or any(c in value for c in '/\\') or value in ('.','..'):
            raise ValueError('Unsafe manifest '+key)
    output=root/'RESULTS'/task['campaign']/task['task_id']
    if output.exists():raise FileExistsError('Never overwrite predictions: '+str(output))
    reference=root/'runtime/frozen_C2_runner.py'
    match=None
    if family=='R13_C2_NUMERICS_WITH_SHARED_TOTAL120_BUDGET' and reference.is_file():
        original=reference.read_text(encoding='utf8')
        expected=original.replace('deadline=min(cap,time.process_time()+cap*(1-proposal_fraction))','deadline=cap')
        match=runner.read_text(encoding='utf8')==expected
        if not match:raise ValueError('Baseline differs from declared one-line shared-budget change')
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),task_id=task['task_id'],
        case_id=task['case_id'],campaign=task['campaign'],method=method,runner_family=family,
        dispatch='runpy same process; no subprocess/restarted CPU clock',runner=str(runner),
        runner_bytes=runner.stat().st_size,historical_reference=str(reference),
        historical_reference_present=reference.is_file(),historical_byte_identity=None,
        declared_budget_only_source_match=match,
        budget_policy='All methods proposal<=72, remaining total120 may be used by original corrector; unchanged numerical solver',
        manifest=str(args.manifest.resolve()),whole_process_CPU_at_dispatch=time.process_time(),
        argv_in=list(sys.argv),truth_access='NONE',wrapper_imports='standard library only',
        no_optimizer_or_backend_imported=True)
    directory=root/'audit/dispatch_provenance'/task['campaign'];directory.mkdir(parents=True,exist_ok=True)
    destination=directory/(task['task_id']+'.json')
    with destination.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    # Preserve all original arguments; only argv[0] names the selected program.
    previous=sys.argv;sys.argv=[str(runner),*previous[1:]]
    try:
        runpy.run_path(str(runner),run_name='__main__')
    finally:
        sys.argv=previous
        record['whole_process_CPU_after_runner']=time.process_time()
        if output.is_dir():
            with (output/'WRAPPER_PROVENANCE.json').open('x') as stream:
                json.dump(record,stream,indent=2);stream.write('\n')


if __name__=='__main__':main()
