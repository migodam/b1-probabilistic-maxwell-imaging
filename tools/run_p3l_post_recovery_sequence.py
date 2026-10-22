"""One paid process at a time after the recorded memory interruption."""
from pathlib import Path
import sys,subprocess,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,write

def main():
    tasks=[('conditional_departure','B','conditional_GN_departure_audit',['analyze_p3l_conditionals.py']),
           ('gaussian_image_diagnostic','B','evaluate_gaussian_reference_images',['run_p3l_accounted_analysis.py','gaussian_images']),
           ('true_projector_geometry','C','geometry_projector_DEV000_DEV028',['run_p3l_geometry.py']),
           ('slice3_quadrature','B','slice3_reference_DEV000',['run_p3l_slice_reference.py'])]
    tasks += [(f'hmc_pilot_c{c}','A',f'hmc_DEV_030_pilot_c{c}',['run_p3l_hmc.py','--chain',str(c)]) for c in range(4)]
    dest=OUT/'schedules/post_memory_recovery';dest.mkdir(parents=True,exist_ok=True)
    path=dest/'TASKS.json';manifest=[dict(name=n,stage=s,job=j,args=a) for n,s,j,a in tasks]
    if path.exists() and json.loads(path.read_text())!=manifest:raise RuntimeError('registered serial sequence differs')
    write(path,manifest);done=[]
    for name,stage,job,args in tasks:
        closed=OUT/stage/'jobs'/job/'CLOSED.json'
        if closed.exists():
            if json.loads(closed.read_text())['status']!='COMPLETE':raise RuntimeError('failed job requires explicit analysis, no restart')
            done.append(dict(name=name,previously_completed=True));continue
        log=(dest/(name+'.log')).open('a')
        child=subprocess.Popen([sys.executable,str(ROOT/'tools'/args[0]),*args[1:]],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        write(dest/'PROGRESS.json',dict(done=done,running=name,PID=child.pid));print('START',name,child.pid,flush=True)
        code=child.wait();log.close();done.append(dict(name=name,exit_code=code));print('END',name,code,flush=True)
        write(dest/'PROGRESS.json',dict(done=done,running=None))
        if code:raise RuntimeError('paid sequence stopped after failure; do not reset/replay')
    write(dest/'CLOSED.json',dict(status='COMPLETE',done=done))

if __name__=='__main__':main()
