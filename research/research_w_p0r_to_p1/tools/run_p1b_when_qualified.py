"""Continue the authorized independent W benchmark after the REAL P1A gate.

No new authorization, daemon or budget. This finite serial coordinator waits
for the already running P1A job receipt, then reuses the same append-only book.
If P1A fails, it records NOT_RUN and performs no fresh P1B physical query.
"""
import sys,json,time,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[1];repo=root.parents[1]
sys.path.insert(0,str(root/'src'))
from wr.accounting import write

def main():
    predecessor=root/'results/jobs/p1a-independent-probability-reference/RECEIPT.json'
    while not predecessor.exists():
        state=root/'QUEUE_CHECKPOINT.json'
        if state.exists() and json.loads(state.read_text()).get('status')=='HALTED_REQUIRES_REVIEW':
            write(root/'P1B_QUEUE_CHECKPOINT.json',dict(status='WAITING_FOR_PREDECESSOR_REVIEW',no_P1B_queries=True));return
        time.sleep(15)
    record=json.loads(predecessor.read_text())
    gate=root/'results/P1A/reference/GATE.json'
    if record['returncode']!=0 or not gate.exists() or json.loads(gate.read_text()).get('probability_reference_qualified') is not True:
        write(root/'P1B_QUEUE_CHECKPOINT.json',dict(status='NOT_RUN_P1A_REFERENCE_GATE_CLOSED',no_P1B_queries=True));return
    jobs=[
      ('p1b-actual-Full-material-rule','wr.progressive_full_bank',[]),
      ('p1b-independent-controlled-observations','wr.prob_data',['P1B']),
      ('p1b-four-arm-inference','wr.progressive_inference',[]),
      ('p1b-point-Full-validation','wr.progressive_point_validation',[]),
      ('p1b-fixed-first-observation-handoff','wr.prob_handoff',[]),
      ('p1b-independent-offline-risk','wr.progressive_offline',[]),
      ('p1b-full-Woodbury-ROM-GMRES-control','wr.green_reuse_control',[]),
    ]
    completed=[]
    for label,module,args in jobs:
        receipt=root/'results/jobs'/label/'RECEIPT.json'
        if receipt.exists():
            if json.loads(receipt.read_text())['returncode']!=0:raise RuntimeError('Failed P1B job needs parent review: '+label)
            completed.append(label);continue
        write(root/'P1B_QUEUE_CHECKPOINT.json',dict(status='RUNNING',label=label,completed=completed,probability_gate='P1A qualified'))
        print('START '+label,flush=True)
        code=subprocess.run([sys.executable,str(root/'tools/run_job.py'),'P1B',label,sys.executable,'-m',module,*args],cwd=repo).returncode
        if code:
            write(root/'P1B_QUEUE_CHECKPOINT.json',dict(status='HALTED_REQUIRES_REVIEW',label=label,returncode=code,completed=completed));raise SystemExit(code)
        completed.append(label);print('COMPLETE '+label,flush=True)
    write(root/'P1B_QUEUE_CHECKPOINT.json',dict(status='P1B_NUMERICAL_WORK_FINISHED_PENDING_SCIENTIFIC_REVIEW',completed=completed))

if __name__=='__main__':main()
