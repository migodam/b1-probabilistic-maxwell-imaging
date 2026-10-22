"""Explicit journal reconciliation after the user's memory-pressure stop.

No physics or job restart. Paid reservations remain charged. A conservative
upper bound on interrupted controlled wall is charged and marked as such.
"""
from pathlib import Path
import sys,json,time,os,fcntl,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write

def main():
    processes=subprocess.check_output(['ps','-axo','pid,args'],text=True)
    if any(('tools/run_p3l_' in line or 'tools/generate_p3l_' in line) and 'recover_p3l_' not in line for line in processes.splitlines()):
        raise RuntimeError('campaign processes still running; do not reconcile active jobs')
    snapshot=json.loads((OUT/'RESOURCE_LEDGER.json').read_text());active=snapshot['active_wall']
    if not active:raise RuntimeError('no interrupted jobs; reconciliation must not repeat')
    stats={job:dict(start=None,last=None,stage=None,RHS=0,states=0,LU=0,physics=0.) for job in active}
    with (OUT/'RESOURCE_LEDGER.jsonl').open() as f:
        for line in f:
            row=json.loads(line);job=row['job']
            if job not in stats:continue
            s=stats[job];s['stage']=row['stage'];s['last']=row['utc_unix']
            if row['event']=='JOB_START':s['start']=row['utc_unix']
            for k in ['RHS','states','LU']:s[k]+=row[k]
            s['physics']+=row['physics_wall_seconds']
    now=time.time();receipts=[]
    for job,limit in active.items():
        s=stats[job];jobdir=OUT/s['stage']/'jobs'/job
        if (jobdir/'CLOSED.json').exists():raise RuntimeError('unexpected closure requires manual audit')
        lock=(jobdir/'job.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if s['start'] is None:raise RuntimeError('missing paid job start')
        lower=max(0.,s['last']-s['start']);upper=min(float(limit),max(lower,now-s['start']))
        book=Book.__new__(Book);book.stage=s['stage'];book.job=job;book.path=OUT/'RESOURCE_LEDGER.jsonl';book.snapshot=OUT/'RESOURCE_LEDGER.json';book.control=(OUT/'budget.lock').open('a');book.context=dict(reason='user reported memory overflow; live process inventory confirms no campaign process')
        receipt=dict(status='INTERRUPTED_MEMORY_PRESSURE',reason='upstream strong-reference transfer cache retained old dense states; all physical reservations preserved',new_physics_RHS=0,paid_RHS=s['RHS'],paid_state_attempts=s['states'],paid_LU_attempts=s['LU'],recorded_physics_wall_seconds=s['physics'],controlled_wall_seconds=upper,controlled_wall_is_conservative_upper_bound=True,observed_wall_lower_bound=lower,unobserved_tail_may_be_overcharged=True,original_reservation_seconds=limit,available_endpoint_files=len(list((jobdir/'endpoints').glob('*.npz'))),attempt_files=len(list((jobdir/'attempts').glob('*.npz'))),no_automatic_replay=True)
        book._update('INTERRUPTED_JOB_RECONCILIATION',controlled=upper,release_wall=True,recovery=receipt)
        write(jobdir/'CLOSED.json',receipt);book.control.close();fcntl.flock(lock,fcntl.LOCK_UN);lock.close();receipts.append(dict(job=job,stage=s['stage'],**receipt))
    write(OUT/'AUDIT/MEMORY_INTERRUPTION_RECOVERY.json',dict(recovered_at_unix=now,new_physics_RHS=0,original_ledgers_not_reset=True,receipts=receipts,implementation_change='campaign-only receiver-transfer cache capacity4, max one paid worker on resume'))
    print(json.dumps(dict(jobs=len(receipts),new_physics_RHS=0,charged_RHS_preserved=sum(v['paid_RHS'] for v in receipts))))

if __name__=='__main__':main()
