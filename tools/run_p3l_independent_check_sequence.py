"""Independent native-coordinate ESS and learned full-MH after fixed training."""
from pathlib import Path
import sys,subprocess,time,json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,write

def main():
    tasks=[(f'full_ess_c{c}','A',f'reference_DEV_030_ess_mode_c{c}',
            ['run_p3l_job.py','reference','DEV_030','--method','ess_mode','--chain',str(c),
             '--warmup','1024','--production','4096']) for c in range(4)]
    tasks.append(('audit_full_ess','A','audit_DEV_030_ess_mode',['audit_p3l_ess.py']))
    tasks.append(('learned_full_correction','D','learned_full_correction',['run_p3l_learned_correction.py']))
    folder=OUT/'schedules/independent_sampler_and_learned_correction';folder.mkdir(parents=True,exist_ok=True)
    manifest=[dict(name=n,stage=s,job=j,args=a) for n,s,j,a in tasks]
    p=folder/'TASKS.json'
    if p.exists() and json.loads(p.read_text())!=manifest:raise RuntimeError('schedule already frozen differently')
    if not p.exists():write(p,manifest)
    dependency=OUT/'schedules/nonlinear_property_revision/CLOSED.json'
    start=time.monotonic()
    while not dependency.exists():
        if time.monotonic()-start>26000:raise RuntimeError('registered predecessor not complete; inspect, no retry')
        time.sleep(15)
    if json.loads(dependency.read_text())['status']!='COMPLETE':raise RuntimeError('predecessor not completed')
    done=[]
    for name,stage,job,args in tasks:
        c=OUT/stage/'jobs'/job/'CLOSED.json'
        if c.exists():
            if json.loads(c.read_text())['status']!='COMPLETE':raise RuntimeError('failed job: do not replay')
            done.append(dict(name=name,previously_completed=True));continue
        with (folder/(name+'.log')).open('x') as log:
            child=subprocess.Popen([sys.executable,str(ROOT/'tools'/args[0]),*args[1:]],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            write(folder/'PROGRESS.json',dict(done=done,running=name,PID=child.pid))
            print('START',name,child.pid,flush=True);code=child.wait()
        done.append(dict(name=name,exit_code=code));print('END',name,code,flush=True)
        write(folder/'PROGRESS.json',dict(done=done,running=None))
        if code:raise RuntimeError('preserve failure and stop sequence for scientific review')
    write(folder/'CLOSED.json',dict(status='COMPLETE',done=done))

if __name__=='__main__':main()
