"""Resume the reviewed zero-RHS offline guard stop without resetting any book.

All previous physical work and the failed scoring receipt remain untouched.
The additional point baseline was registered before P1B data generation.
This finite serial queue has no daemon, scheduler, or extra authorization.
"""
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from wr.accounting import write


def read(path):return json.loads(path.read_text())


def main():
    if (ROOT/'CLOSED.json').exists():raise RuntimeError('Campaign already closed')
    with (ROOT/'REVIEWED_RESUME.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if read(ROOT/'results/P1A/reference/GATE.json').get('probability_reference_qualified') is not True:
            raise RuntimeError('P1A real reference gate not qualified')
        failed=ROOT/'results/jobs/p1b-independent-offline-risk'
        receipt=read(failed/'RECEIPT.json')
        if receipt.get('status')!='FAILED' or receipt.get('returncode')!=1:
            raise RuntimeError('Expected preserved guard failure receipt differs')
        if 'OFFLINE truth access blocked' not in (failed/'OUTPUT.log').read_text():
            raise RuntimeError('Unexpected previous failure, requires separate review')
        registration='P1B_STRONG_POINT_CONTROL_ADDENDUM.json'
        snapshot=ROOT/'results/jobs/p1b-independent-controlled-observations/CONFIG_SNAPSHOT'/registration
        if snapshot.read_bytes()!=(ROOT/'configs'/registration).read_bytes():
            raise RuntimeError('Continuous point control was not fixed before P1B observations')
        audit=ROOT/'AUDIT/point_guard_recovery';audit.mkdir(exist_ok=True)
        preserved=audit/'ORIGINAL_P1B_QUEUE_CHECKPOINT.json'
        if not preserved.exists():shutil.copyfile(ROOT/'P1B_QUEUE_CHECKPOINT.json',preserved)
        record=audit/'REVIEW.json'
        if not record.exists():
            if (ROOT/'results/P1B/offline_scores/TRUTH_ACCESS.json').exists():
                raise RuntimeError('Offline truth was accessed before completion of the point control')
            write(record,dict(reviewed_utc=datetime.now(timezone.utc).isoformat(),
                prior_failure='Explicit pre-GT guard rejected scoring before additional point control was complete',
                prior_failed_job=str(failed.relative_to(ROOT)),failed_receipt_preserved=True,
                guard_failure_new_RHS=0,GT_access_marker_absent=True,
                addendum_literal_copy_present_in_pre_observation_config_snapshot=True,
                budget_or_run_reset=False,old_RHS_and_CPU_all_remain_charged=True,
                authorization='Existing USER_RESEARCH_W_P0R_TO_P1_20261023 only'))
        jobs=[
            ('p1b-continuous-Full-validated-point-control','wr.continuous_point_control'),
            ('p1b-independent-offline-risk-after-guard-review','wr.progressive_offline'),
            ('p1b-full-Woodbury-ROM-GMRES-control','wr.green_reuse_control')]
        completed=[]
        for label,module in jobs:
            dest=ROOT/'results/jobs'/label
            previous=dest/'RECEIPT.json'
            if previous.exists():
                old=read(previous)
                if old.get('status')!='COMPLETE' or old.get('returncode')!=0:
                    raise RuntimeError('Failed recovery job preserved; no automatic retry: '+label)
                completed.append(label);continue
            if dest.exists():raise RuntimeError('Partial recovery job, parent review required: '+label)
            write(ROOT/'P1B_QUEUE_CHECKPOINT.json',dict(status='RUNNING',label=label,
                completed=completed,previous_guard_stop_preserved=str(preserved.relative_to(ROOT))))
            print('START '+label,flush=True)
            code=subprocess.run([sys.executable,str(ROOT/'tools/run_job.py'),'P1B',label,
                                 sys.executable,'-m',module],cwd=ROOT.parents[1]).returncode
            if code:
                write(ROOT/'P1B_QUEUE_CHECKPOINT.json',dict(status='HALTED_REQUIRES_REVIEW',
                    label=label,returncode=code,completed=completed,prior_failure_preserved=True))
                raise SystemExit(code)
            completed.append(label);print('COMPLETE '+label,flush=True)
        write(ROOT/'P1B_QUEUE_CHECKPOINT.json',dict(
            status='P1B_NUMERICAL_WORK_FINISHED_PENDING_SCIENTIFIC_REVIEW',
            reviewed_recovery_completed=completed,old_completed_jobs_preserved=True))
        code=subprocess.run([sys.executable,str(ROOT/'tools/finish_b2_after_p1b.py')],
                            cwd=ROOT.parents[1]).returncode
        if code:raise SystemExit(code)


if __name__=='__main__':main()
