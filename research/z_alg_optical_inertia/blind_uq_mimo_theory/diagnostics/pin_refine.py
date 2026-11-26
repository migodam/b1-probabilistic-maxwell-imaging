"""Allow exterior compensation after an exact single-site update reaches its limit.
The array and ROI are fixed by the previous blind reconstruction and physical profile.
"""
from pathlib import Path
import argparse,json,os,sys,time
for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'): os.environ[k]='1'
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from verified_wave import P,VerifiedWave,up,rd,fro
sys.path.insert(0,str(P/'imaging'))
from shared import Model

class Feasible(Exception):pass
def work(case,root,limit=96):
    m=Model.read(case)
    wd=P/'results'/root/case
    vfile=P/'ports/public/operators'/f'{m.operator_id}.npz'
    with np.load(vfile) as a:g=VerifiedWave(a)
    for f in sorted(wd.glob('roi_??.npz')):
        inner=f.with_name(f.stem+'_inner.json'); out=f.with_name(f.stem+'_pin.json')
        if not inner.exists() or out.exists():continue
        report=json.loads(inner.read_text());run=time.perf_counter();cpu=time.process_time()
        if report['attainable_peak_to']==4.5:continue
        with np.load(f) as d: x=d['material'].copy();mask=d['roi']
        idx=report['selected']['pixel'] if report['selected'] else int(np.flatnonzero(mask)[0])
        x[idx]=4.5
        lo=np.zeros_like(x);lo[idx]=4.5;hi=np.full_like(x,4.5)
        calls=0;zbest=x.copy();best=np.inf
        def fg(z):
            nonlocal calls,best,zbest
            loss,grad=m.data_grad(z);ratio=m.forward(z)['residual_ratio']
            calls+=1
            if ratio<best:best=ratio;zbest=z.copy()
            if ratio<.82:raise Feasible
            return loss,grad
        try:minimize(fg,x,method='L-BFGS-B',jac=True,bounds=np.column_stack((lo,hi)),options=dict(maxiter=limit,ftol=1e-14,gtol=1e-10,maxls=25))
        except Feasible:pass
        c,proof=g.inverse(zbest);state=m.forward(zbest);cur=state['j'];b=zbest[:,None]*g.B
        e=g.solve_error(zbest,cur,b,up(rd(1)*max(zbest)*fro(g.B)),c)
        pred=g.S@cur
        ru=up((fro(pred-m.Y)+g.snorm*e+rd(g.n)*g.snorm*fro(cur)+rd(2)*(fro(pred)+fro(m.Y)))/m.radius)
        np.savez_compressed(f.with_name(f.stem+'_pin.npz'),material=zbest,roi=mask,current=cur,pixel=idx)
        out.write_text(json.dumps(dict(case=case,roi=report['roi'],pixel=idx,peak=4.5,peak_other_endpoint=report['attainable_peak_from'],
            original_residual_ratio_upper=ru,verified=bool(ru<=1 and zbest[idx]==4.5),
            endpoint_hull_only=True,exact_continuous_path=False,inverse_replay=proof,loss_and_adjoint_calls=calls,
            physics=dict(m.cost),wall_s=time.perf_counter()-run,cpu_s=time.process_time()-cpu),indent=2))
        print(case,report['roi'],ru,round(time.perf_counter()-run,2),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('cases',nargs='+');a.add_argument('--kind',default='profiles');a.add_argument('--limit',type=int,default=96);q=a.parse_args()
    with threadpool_limits(1):
        for case in q.cases:work(case,q.kind,q.limit)
