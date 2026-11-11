"""Combine completed-child accounting, preserving the failed controller record."""
import argparse,json,time
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
    q=r/'queues/development_v1';p2=r/'queues/development_v1_recovery/RECEIPT.json'
    second=json.loads(p2.read_text());first=json.loads((q/'RECEIPT.tmp').read_text())
    if second['pending'] or second['running']:raise RuntimeError('Recovery still active')
    jobs=first['jobs']+second['jobs'];out=r/'RESULTS/development_v1'
    for j in jobs:
        j['child_status']=json.loads((out/j['task_id']/'RECEIPT.json').read_text())['status']
    old=q/'RECEIPT.json'
    if old.exists():old.rename(q/'FAILED_CONTROLLER_LAST_RECEIPT.json')
    complete=dict(status='COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY',campaign='development_v1',pending=0,running=[],jobs=jobs,
        finished_child_CPU_seconds=sum(j['external']['child_process_CPU_seconds'] for j in jobs),
        controller_CPU_seconds=first['controller_CPU_seconds']+second['controller_CPU_seconds'],
        wall_seconds=first['wall_seconds']+second['wall_seconds'],
        wall_scope='sum of two active-controller windows; excludes intervention gap; not child-wall sum',
        end_to_end_wall_seconds_including_intervention=time.time()-q.stat().st_ctime,
        unknown_aborted_child_attempts=4,aborted_child_CPU_seconds=None,
        uncertainty='Four children terminated on failed controller/SSH exit before final receipts; unmeasured costs retained, never set to zero.',
        recovery_source='../development_v1_recovery/RECEIPT.json',original_source='RECEIPT.tmp',
        old_successful_runs_repeated=0)
    old.write_text(json.dumps(complete,indent=2)+'\n');print(json.dumps({k:complete[k] for k in ('status','finished_child_CPU_seconds','unknown_aborted_child_attempts')}))


if __name__=='__main__':main()
