"""Record actual host process wall/CPU and physics reservation delta."""
import argparse,json,time,subprocess,resource,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
root=Path('results');root.mkdir(exist_ok=True); path=root/'JOB_LEDGER.jsonl'
cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
if not cmd: raise SystemExit('command required')
if path.exists() and any(json.loads(l).get('id')==a.id for l in path.read_text().splitlines()): raise SystemExit('unique job id required')
def count():
 f=root/'PHYSICS_LEDGER.jsonl';return sum(json.loads(l)['event']=='BEGIN' for l in f.read_text().splitlines()) if f.exists() else 0
before=count();t=time.perf_counter(); r0=resource.getrusage(resource.RUSAGE_CHILDREN)
job=root/'jobs'/a.id;job.mkdir(parents=True,exist_ok=False)
with (job/'stdout.log').open('w') as out,(job/'stderr.log').open('w') as err:
 r=subprocess.run(cmd,stdout=out,stderr=err)
r1=resource.getrusage(resource.RUSAGE_CHILDREN)
row=dict(id=a.id,command=cmd,returncode=r.returncode,wall_seconds=time.perf_counter()-t,cpu_seconds=r1.ru_utime+r1.ru_stime-r0.ru_utime-r0.ru_stime,physics_evaluations=count()-before)
with path.open('a') as f:f.write(json.dumps(row)+'\n')
print(json.dumps(row));print((job/'stdout.log').read_text()[-2500:]);print((job/'stderr.log').read_text()[-1500:]);sys.exit(r.returncode)
