"""Finite serial queue for the already authorized, frozen W campaign.

No service/daemon, no budget refresh. Every paid subprocess still uses the
same stage/scoped append-only Book. Unexpected errors stop the queue.
"""
import sys,json,time,subprocess
from pathlib import Path

root=Path(__file__).resolve().parents[1];repo=root.parents[1]
python=Path(sys.executable).absolute()
sys.path.insert(0,str(root/'src'))
from wr.accounting import write

jobs=[
 ('C','c-fixed-pair-band-aperture-loss','wr.phase_diagram',['configs/FIRST_JOINT_NOMINATION.json','--include-loss']),
 ('C','c-wide-reoptimized-pair','wr.wide_reopt',['configs/FIRST_JOINT_NOMINATION.json']),
 ('B2','b2-task-preserving-2026102202','wr.task_preserving',['2026102202','narrow']),
 ('B2','b2-native-cavity-2026102202','wr.native_edit',['2026102202','narrow','cavity']),
 ('B2','b2-native-cal-2026102201','wr.native_cal',['2026102201','narrow']),
 ('B2','b2-native-joint-2026102201','wr.native_edit',['2026102201','narrow','joint']),
 ('P1A','p1a-independent-controlled-observations','wr.prob_data',['P1A']),
 ('P1A','p1a-full-Maxwell-quadrature-bank','wr.prob_bank',[]),
 ('P1A','p1a-independent-probability-reference','wr.prob_reference',[]),
]

def main():
    # This prior is published before any new data; do not create or edit it here.
    if not (root/'configs/W_CONTROLLED_PRIOR.json').exists():raise RuntimeError('controlled family not frozen')
    previous=root/'results/jobs/b2-joint-stress-and-ablation/RECEIPT.json'
    while not previous.exists():time.sleep(5)
    if json.loads(previous.read_text())['returncode']!=0:raise RuntimeError('preceding physics job failed; inspect before continuation')
    completed=[]
    for stage,label,module,args in jobs:
        receipt=root/'results/jobs'/label/'RECEIPT.json'
        if receipt.exists():
            if json.loads(receipt.read_text())['returncode']!=0:raise RuntimeError('prior failed job requires explicit review: '+label)
            completed.append(label);continue
        write(root/'QUEUE_CHECKPOINT.json',dict(status='RUNNING',active_label=label,stage=stage,completed=completed,future=[r[1] for r in jobs if r[1] not in completed and r[1]!=label]))
        print('START '+label,flush=True)
        code=subprocess.run([str(python),str(root/'tools/run_job.py'),stage,label,str(python),'-m',module,*args],cwd=repo).returncode
        if code:
            write(root/'QUEUE_CHECKPOINT.json',dict(status='HALTED_REQUIRES_REVIEW',label=label,returncode=code,completed=completed))
            raise SystemExit(code)
        completed.append(label);print('COMPLETE '+label,flush=True)
    write(root/'QUEUE_CHECKPOINT.json',dict(status='P1A_REFERENCE_FINISHED_REVIEW_GATE_BEFORE_P1B',completed=completed))

if __name__=='__main__':main()
