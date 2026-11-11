"""One finite, foreground research sequence; exits when all queues finish.

Waits for the separate development process to release its P cores. Never
interrupts it or starts another task's workload. No scientific quality gates.
"""
from pathlib import Path
import argparse,json,subprocess,sys,time


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True)
    p.add_argument('--development-root',type=Path,required=True);a=p.parse_args();r=a.root
    dependency=a.development_root/'queues/representation_v1/RECEIPT.json'
    print('Waiting for the finite development queue to release its four P cores.',flush=True)
    while not dependency.exists():time.sleep(1)
    previous=json.loads(dependency.read_text())
    if previous['pending'] or previous['running']:raise RuntimeError('Dependency has not terminated')
    cores='0,2,4,6,8,10,12,14';logs=[];wall=time.perf_counter()
    def run(script,*arguments):
        argv=[sys.executable,str(r/'src'/script),*map(str,arguments)]
        start=time.perf_counter();code=subprocess.call(argv,cwd=r)
        logs.append(dict(script=script,arguments=list(map(str,arguments)),returncode=code,wall_seconds=time.perf_counter()-start))
        if code!=0:raise RuntimeError(f'Technical process failure: {script}, return code {code}; evidence retained for explicit diagnosis/recovery.')
    def campaign(name):
        run('finite_queue.py','--root',r,'--legacy-repo',a.legacy_repo,'--manifest',r/'configs'/(name+'.json'),'--cores',cores)
        run('archive_campaign.py','--root',r,'--campaign',name)
    campaign('confirmation_core48')
    campaign('confirmation_mechanism16')
    run('prepare_conditional_rescue.py','--root',r)
    campaign('confirmation_rescue48')
    campaign('physics_v1')
    summary=dict(status='FINITE_SEQUENCE_COMPLETED',truth_access='NONE',wall_seconds=time.perf_counter()-wall,subprocesses=logs,
                 P_cores=list(map(int,cores.split(','))),no_new_daemon=True,no_GTX_scoring=True)
    (r/'FINITE_SEQUENCE_RECEIPT.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
