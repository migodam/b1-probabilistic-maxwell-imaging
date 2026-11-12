"""One finite dependency-only recovery after the C3 control queue releases CPUs.

The original C2 source, data, successful cases and archive are never changed.
This driver adds only the two missing, unchanged inherited R1 modules to C4.
It neither reads truth nor changes a numerical algorithm or budget.
"""
from pathlib import Path
import argparse,json,subprocess,sys,time


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--legacy-repo',type=Path,required=True)
    p.add_argument('--controls-root',type=Path,required=True);a=p.parse_args()
    dependency=a.controls_root/'queues/confirmation_online_controls16/RECEIPT.json'
    print('Waiting for the finite C3 control queue to release its P cores.',flush=True)
    while not dependency.exists():time.sleep(1)
    previous=json.loads(dependency.read_text())
    if previous.get('pending') or previous.get('running'):
        raise RuntimeError('Control queue is not terminal')
    started=time.perf_counter();logs=[]
    for script,arguments in [
        ('finite_queue.py',['--root',a.root,'--legacy-repo',a.legacy_repo,
            '--manifest',a.root/'configs/physics_recovery_v1.json','--cores','0,2,4,6']),
        ('archive_campaign.py',['--root',a.root,'--campaign','physics_recovery_v1'])]:
        t=time.perf_counter();code=subprocess.call([sys.executable,str(a.root/'src'/script),*map(str,arguments)],cwd=a.root)
        logs.append(dict(script=script,returncode=code,wall_seconds=time.perf_counter()-t))
        if code:raise RuntimeError('Technical recovery subprocess failed: '+script)
    receipt=dict(status='FINITE_RECOVERY_SEQUENCE_COMPLETED',subprocesses=logs,
        wall_seconds=time.perf_counter()-started,P_cores=[0,2,4,6],BLAS_threads=1,
        reason='unchanged inherited R1 modules were absent from original remote deployment',
        algorithm_changed=False,truth_access='NONE')
    (a.root/'FINITE_RECOVERY_SEQUENCE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
