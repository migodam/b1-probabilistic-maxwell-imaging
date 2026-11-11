"""Uniform pre-score extra physical-corrector budget; original endpoints immutable."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,time,json,resource,subprocess,argparse
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).absolute().parent
sys.path.insert(0,str(HERE))
from run_vector import model,CASES,METHODS,dump
CONTRACT=dict(status='FROZEN_BEFORE_ANY_GT_SCORE',reason='Root requested common physical corrector budget diagnostic after all observed max120 stops, before any material scoring',maxiter=1000,CPU_cap_seconds=60,ftol=1e-13,gtol=1e-9,maxls=30,start='Each arm selected original120 endpoint',selection='minimum original physical Phi over initial and accepted callbacks; final returned point also audited',preserve='All original120 predictions and costs immutable',truth_access='NONE')
class BudgetReached(BaseException):pass
def run(cid,method):
    cfg=json.loads((HERE/'CONTINUATION_CONTRACT.json').read_text());base=HERE/'RESULTS'/cid/method
    out=base/'long_corrector';out.mkdir(exist_ok=False);startcpu=time.process_time();wall=time.perf_counter();m=model();Y=np.load(HERE/'data_public'/(cid+'.npz'))['Y'];z=np.load(base/'prediction.npz')['prediction'].real
    best=z.copy();f,df,J,E=m.physical(best,Y,False);bestphi=f;initialphi=f;curve=[];trials=[];solvercpu=time.process_time()
    def check():
        if time.process_time()-solvercpu>=cfg['CPU_cap_seconds']:raise BudgetReached()
    def fun(x):check();return m.physical(x,Y)
    def callback(x):
        nonlocal best,bestphi
        f,df,J,E=m.physical(x,Y,False);trials.append(x.copy());curve.append(dict(cpu_seconds=time.process_time()-solvercpu,phi=f,physical_LU=m.counts['physical_LU']))
        if f<bestphi:bestphi=f;best=x.copy()
        check()
    try:
        res=minimize(fun,z,jac=True,method='L-BFGS-B',bounds=[(0,4)]*m.N,callback=callback,options=dict(maxiter=cfg['maxiter'],ftol=cfg['ftol'],gtol=cfg['gtol'],maxls=cfg['maxls']))
        callback(res.x);stop=dict(nit=int(res.nit),nfev=int(res.nfev),success=bool(res.success),message=str(res.message))
    except BudgetReached:stop=dict(success=False,message='COMMON_ADDITIONAL_CPU_CAP',accepted_callbacks=len(trials))
    f,df,J,E=m.physical(best,Y,False)
    np.savez_compressed(out/'prediction.npz',initial=z,prediction=best,J=J,E=E,trials=np.asarray(trials))
    dump(out/'RECEIPT.json',dict(status='COMPLETE',case_id=cid,method=method,initial_phi=float(initialphi),physical_phi=float(f),physical_data_residual=float(np.sqrt(df)),CPU_seconds=time.process_time()-startcpu,wall_seconds=time.perf_counter()-wall,process_CPU_near_exit=time.process_time(),counts=m.counts,optimizer=stop,curve=curve,truth_access='NONE'))
    print(json.dumps(dict(case=cid,method=method,phi=f,CPU_seconds=time.process_time()-startcpu,stop=stop)),flush=True)
def queue():
    cp=HERE/'CONTINUATION_CONTRACT.json'
    if not cp.exists():dump(cp,CONTRACT)
    logs=[]
    for c in CASES:
        for k in METHODS:
            path=HERE/'RESULTS'/c/k/'RECEIPT.json'
            while not path.exists():time.sleep(2)
            if (path.parent/'long_corrector/RECEIPT.json').exists():continue
            t=time.time();before=resource.getrusage(resource.RUSAGE_CHILDREN);p=subprocess.run([sys.executable,__file__,'--case',c,'--method',k],capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
            logs.append(dict(case_id=c,method=k,returncode=p.returncode,external_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall_seconds=time.time()-t,stdout=p.stdout,stderr=p.stderr))
            dump(HERE/'CONTINUATION_QUEUE.json',dict(status='RUNNING' if not p.returncode else 'FAILED',rows=logs));print(p.stdout or p.stderr,flush=True)
            if p.returncode:sys.exit(p.returncode)
    dump(HERE/'CONTINUATION_QUEUE.json',dict(status='COMPLETE',rows=logs))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case');p.add_argument('--method');a=p.parse_args()
    if a.case:run(a.case,a.method)
    else:queue()
