"""Finish independent checks before the one frozen, new locked risk panel."""
from pathlib import Path
import sys,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,write

def main():
    tasks=[('simple_full_MH','D','simple_correction_exposed_DEV',['run_p3l_simple_correction.py']),
           ('locked_truth','D','dataset_locked',['generate_p3l_observations.py','locked']),
           ('locked_online','D','locked_imaging_frozen_panel',['run_p3l_locked_imaging.py']),
           ('locked_full_IS','D','locked_importance_frozen_128',['run_p3l_locked_importance.py']),
           ('locked_offline_scoring','D','evaluate_locked_imaging',['run_p3l_accounted_analysis.py','locked_imaging'])]
    folder=OUT/'schedules/frozen_locked_risk';folder.mkdir(parents=True,exist_ok=True)
    manifest=[dict(name=n,stage=s,job=j,args=a) for n,s,j,a in tasks]
    p=folder/'TASKS.json'
    if p.exists() and json.loads(p.read_text())!=manifest:raise RuntimeError('frozen sequence differs')
    if not p.exists():write(p,manifest)
    dependency=OUT/'schedules/independent_sampler_and_learned_correction/CLOSED.json'
    start=time.monotonic()
    while not dependency.exists():
        if time.monotonic()-start>26000:raise RuntimeError('existing checks unfinished: inspect, do not restart physics')
        time.sleep(15)
    if json.loads(dependency.read_text())['status']!='COMPLETE':raise RuntimeError('independent sequence not complete')
    if not (OUT/'METHODS_FROZEN_FOR_LOCKED.json').exists():raise RuntimeError('methods freeze absent')
    version=json.loads((OUT/'configs/LOCKED_CODE_VERSION.json').read_text())
    done=[]
    for name,stage,job,args in tasks:
        c=OUT/stage/'jobs'/job/'CLOSED.json'
        if c.exists():
            if json.loads(c.read_text())['status']!='COMPLETE':raise RuntimeError('failure retained, no reset')
            done.append(dict(name=name,previously_completed=True));continue
        # Runtime results may be dirty; estimator/kernel source must remain
        # exactly at the committed version recorded before observations.
        dirty=subprocess.check_output(['git','diff',version['git_commit'],'--name-only','--','src/p3l','src/p1_material','tools/run_p3l_locked_imaging.py','tools/run_p3l_locked_importance.py','tools/generate_p3l_observations.py','tools/evaluate_p3l_locked_imaging.py','tools/run_p3l_joint_gaussian.py'],cwd=ROOT,text=True)
        if dirty.strip():raise RuntimeError('frozen estimator/target source changed; no held-out execution')
        with (folder/(name+'.log')).open('x') as log:
            child=subprocess.Popen([sys.executable,str(ROOT/'tools'/args[0]),*args[1:]],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            write(folder/'PROGRESS.json',dict(done=done,running=name,PID=child.pid))
            print('START',name,child.pid,flush=True);code=child.wait()
        done.append(dict(name=name,exit_code=code));print('END',name,code,flush=True)
        write(folder/'PROGRESS.json',dict(done=done,running=None))
        if code:raise RuntimeError('stop after failed job; preserve all states and costs')
    write(folder/'CLOSED.json',dict(status='COMPLETE',done=done))

if __name__=='__main__':main()
