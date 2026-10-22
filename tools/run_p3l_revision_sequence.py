"""Registered serial revision; no concurrent paid workers or automatic retries."""
from pathlib import Path
import sys
import subprocess
import json
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,write

def main():
    tasks=[('material_metric_bracket','C','material_bracket_saved_arrays',['analyze_p3l_material_brackets.py']),
           ('property_spine_build','B','build_property_spine_DEV030',['run_p3l_spine.py','--build'])]
    tasks += [(f'spine_hmc_c{c}','B',f'spine_hmc_DEV030_pilot_c{c}',['run_p3l_spine.py','--chain',str(c)]) for c in range(4)]
    tasks += [('audit_spine','B','audit_spine_DEV030_pilot',['audit_p3l_spine.py'])]
    for model in ['mlp_diag','mlp_full','mlp_mdn','mlp_flow','cnn_full']:
        for seed in [2026102270,2026102271,2026102272]:
            tasks.append((f'train_{model}_{seed}','D',f'train_{model}_s{seed}',
                          ['train_p3l_density.py',model,'--seed',str(seed)]))
    dest=OUT/'schedules/nonlinear_property_revision';dest.mkdir(parents=True,exist_ok=True)
    manifest=[dict(name=n,stage=s,job=j,args=a) for n,s,j,a in tasks]
    path=dest/'TASKS.json'
    if path.exists() and json.loads(path.read_text())!=manifest:
        raise RuntimeError('frozen serial revision differs; cannot rewrite history')
    if not path.exists():write(path,manifest)
    dependency=OUT/'B/jobs/joint_probit_gaussian/CLOSED.json'
    started=time.monotonic()
    while not dependency.exists():
        write(dest/'PROGRESS.json',dict(waiting_for='joint_probit_gaussian',no_paid_work_launched=True))
        if time.monotonic()-started>7200:
            raise RuntimeError('existing job not completed; inspect before any new paid worker')
        time.sleep(10)
    if json.loads(dependency.read_text())['status']!='COMPLETE':
        raise RuntimeError('dependency failed; inspect partial artifacts, do not replay')
    done=[]
    for name,stage,job,args in tasks:
        closed=OUT/stage/'jobs'/job/'CLOSED.json'
        if closed.exists():
            if json.loads(closed.read_text())['status']!='COMPLETE':
                raise RuntimeError('failed job requires explicit analysis, never a reset')
            done.append(dict(name=name,previously_completed=True));continue
        with (dest/(name+'.log')).open('x') as log:
            child=subprocess.Popen([sys.executable,str(ROOT/'tools'/args[0]),*args[1:]],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            write(dest/'PROGRESS.json',dict(done=done,running=name,PID=child.pid))
            print('START',name,child.pid,flush=True)
            code=child.wait()
        done.append(dict(name=name,exit_code=code))
        print('END',name,code,flush=True)
        write(dest/'PROGRESS.json',dict(done=done,running=None))
        if code:raise RuntimeError('serial revision stopped on failure; inspect, no implicit retry')
    write(dest/'CLOSED.json',dict(status='COMPLETE',done=done))

if __name__=='__main__':main()
