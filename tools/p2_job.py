"""Meter one external zero-physics or training job into the P2 host budget."""
import argparse,json,time,subprocess,resource,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.b1_joint.accounting import Ledger,dump
p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
if not cmd:raise SystemExit('command required')
root=ROOT/'research/joint_p2/results';ledger=Ledger(root);remaining=7200-ledger.charged_wall()
folder=root/'jobs'/a.id;folder.mkdir(parents=True,exist_ok=False)
r0=resource.getrusage(resource.RUSAGE_CHILDREN)
try:
    with ledger.job(a.id):
        with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
            r=subprocess.run(cmd,stdout=out,stderr=err,timeout=max(1,remaining))
        if r.returncode:raise RuntimeError(f'external job exit {r.returncode}')
finally:
    r1=resource.getrusage(resource.RUSAGE_CHILDREN)
    dump(folder/'CHILD_CPU.json',dict(child_cpu_seconds=r1.ru_utime+r1.ru_stime-r0.ru_utime-r0.ru_stime))
    print((folder/'stdout.log').read_text()[-3500:]);print((folder/'stderr.log').read_text()[-1500:])
