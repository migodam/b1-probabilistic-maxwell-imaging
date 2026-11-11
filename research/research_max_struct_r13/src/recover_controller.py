"""Explicit recovery after the recorded Windows status-replace failure.

Preserves every completed run, failed partial directory and original controller
receipt. It is not an automatic scientific retry policy.
"""
from pathlib import Path
import argparse,json,shutil,time


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
    old=root/'queues/development_v1';snap=json.loads((old/'RECEIPT.tmp').read_text())
    config=json.loads((root/'configs/development_v1.json').read_text());pending=[];partial=[]
    for task in config['tasks']:
        out=root/'RESULTS'/task['campaign']/task['task_id']
        if (out/'RECEIPT.json').exists():continue
        if out.exists():
            dest=root/'audit/aborted_controller_children'/task['task_id'];dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.move(str(out),str(dest));partial.append(dict(task_id=task['task_id'],saved_partial=str(dest.relative_to(root)),
                CPU_seconds=None,reason='child terminated on SSH/controller exit before its final receipt; no recoverable process handle'))
        pending.append(task)
    config['campaign']='development_v1_recovery';config['tasks']=pending
    (root/'configs/development_v1_recovery.json').write_text(json.dumps(config,indent=2)+'\n')
    record=dict(status='EXPLICIT_RECOVERY_PREPARED',original_error='Windows WinError5 replacing open status JSON',
        original_completed=len(snap['jobs']),recovered_tasks=len(pending),partial_attempts=partial,
        original_completed_child_CPU_seconds=snap['finished_child_CPU_seconds'],original_controller_snapshot=snap,
        existing_complete_outputs_unchanged=True,scientific_retry_not_automatic=True)
    (root/'audit/CONTROLLER_FAILURE_AND_RECOVERY.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(complete=record['original_completed'],pending=len(pending),partial=len(partial))))


if __name__=='__main__':main()
