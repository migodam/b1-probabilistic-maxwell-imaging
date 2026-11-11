"""Finite R12 queue; binds only its own children to distinct P cores."""
from pathlib import Path
import os,sys,subprocess,time,json,concurrent.futures
import psutil
ROOT=Path(__file__).resolve().parents[1]
TASKS=[(c,m) for c in ['r10_twod_0001','r10_twod_0004','r10_twod_0022']
       for m in ['direct','multistart','alternating_projection','douglas_rachford']]

def lane(tasks,core):
    result=[]
    for cid,method in tasks:
        out=ROOT/'audit/remote_logs';out.mkdir(parents=True,exist_ok=True)
        log=out/f'{cid}__{method}.log'
        code='import psutil,runpy,sys; psutil.Process().cpu_affinity(['+str(core)+']); sys.argv='+repr(['basin_escape.py','--case',cid,'--method',method])+'; runpy.run_path('+repr(str(ROOT/'experiments/basin_escape.py'))+',run_name="__main__")'
        env=dict(os.environ,MAXWELL_REUSE_REPO='D:/AI/MAX_STRUCT_R11_DEV_20261107',
                 PYTHONPATH='D:/AI/MAX_STRUCT_R12_20261111;D:/AI/MAX_STRUCT_R11_DEV_20261107;D:/AI/MAX_STRUCT_R11_DEV_20261107/vendor',
                 OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        started=time.time()
        with log.open('w') as f:p=subprocess.run([sys.executable,'-c',code],stdout=f,stderr=subprocess.STDOUT,env=env)
        row=dict(case=cid,method=method,core=core,returncode=p.returncode,start=started,wall=time.time()-started)
        result.append(row);(out/f'{cid}__{method}.launch.json').write_text(json.dumps(row,indent=2))
    return result

if __name__=='__main__':
    t=time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        jobs=[pool.submit(lane,TASKS[i::4],2*i) for i in range(4)]
        rows=[r for j in jobs for r in j.result()]
    (ROOT/'audit/REMOTE_QUEUE.json').write_text(json.dumps(dict(tasks=rows,wall=time.time()-t),indent=2))
