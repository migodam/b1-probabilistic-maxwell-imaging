"""Account for every paid online branch and actual interval replay separately."""

import json

from search import HERE


def read(path):
    return json.loads(path.read_text())


def time_seconds(value):
    minutes, seconds = value.rsplit(':',1)
    if ':' in minutes:
        hours, minutes = minutes.split(':')
    else:
        hours=0
    return 3600*float(hours)+60*float(minutes)+float(seconds)


def main():
    branches=[]
    for folder in ('runs','cie','gn','gn_short','continuation','wri'):
        for path in sorted((HERE/folder).glob('*.json')):
            payload=read(path)
            branches.append(dict(branch=folder,observation=payload['observation'],
                                 receipt=str(path.relative_to(HERE)),
                                 wall_seconds=payload['wall_seconds'],cpu_seconds=payload['cpu_seconds'],
                                 operations=payload.get('counts',payload.get('native_counts',{})),
                                 inherited_cie_operations=payload.get('cie_counts'),
                                 projected_state_operations=payload.get('projected_counts'),
                                 numerical_result=payload.get('selected_numeric_residual',
                                     payload.get('seed_residual',payload.get('native_residual')))))
    verifier_runs=[]
    for path in sorted(HERE.glob('verified*.json')):
        data=read(path)
        verifier_runs.append(dict(receipt=path.name,wall_seconds=data['wall_seconds'],
                                  cpu_seconds=data['cpu_seconds'],shared_setup=data['geometry_prep'],
                                  statuses={k:[r['status'] for r in v] for k,v in data['cases'].items()}))
    interrupted=[]
    cut=HERE/'runs/interrupted_redundant_second_pass.txt'
    if cut.exists():
        for line in cut.read_text().splitlines():
            if ' --object optfresh_' not in line:
                continue
            pid,wall,cpu,command = line.strip().split(maxsplit=3)
            name=command.split(' --object ',1)[1].split()[0]
            completed=next(branch for branch in branches if branch['branch']=='gn' and
                           branch['observation']==f'{name}__pw360_10x20__pair1pct')
            interrupted.append(dict(pid=int(pid),object=name,
                                    process_cpu_at_interrupt=time_seconds(cpu),
                                    paid_cpu_outside_completed_first_case=max(0.,time_seconds(cpu)-completed['cpu_seconds']),
                                    process_wall_at_interrupt=time_seconds(wall),
                                    status='Redundant longer 2% search superseded by a completed independent 2% branch.'))
    out=dict(platform='Apple M4; single thread per branch or proof; multiple independent branches overlap',
             wall_scope='Per-view wall includes actual parallel contention. Do not sum per-view wall and call it total study wall.',
             online_candidates=branches,actual_interval_runs=verifier_runs,
             online_branch_cpu_sum=sum(r['cpu_seconds'] for r in branches),
             interval_run_cpu_sum=sum(r['cpu_seconds'] for r in verifier_runs),
             interrupted_with_paid_cpu=interrupted,
             supplementary_notes=['No waveform truth was read by any paid online branch.',
                 'The exploratory projected WRI geometry preparation had 256 native FFT columns; its paid shared timing is in runs/wri02.log.',
                 'The first continuation initialization of 05 2% stopped before a state solve because the expanded support did not include all previously active cells. The corrected initialization keeps every selected positive cell. Its failed startup time was not separately retained.',
                 'Archived lower-bound proposal and proof times are included per row in lower_bounds.json and were not rerun as part of this inversion.'])
    (HERE/'cost_ledger.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
