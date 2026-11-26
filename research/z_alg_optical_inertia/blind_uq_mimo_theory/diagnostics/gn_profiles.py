"""Low-dimensional data normal equations; every material pixel remains free.
Use only the original training readout, a frozen blind image and its fitted maps.
The projected nonlinear iterations produce physical candidates. Verification is separate.
"""
import argparse,json,os,sys,time
from pathlib import Path
for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.linalg import solve
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'imaging'))
from shared import Model
def fit(case):
    m=Model.read(case);base=P/'imaging/results'/case
    stationary=base/'stationary'; dest=P/'results/gn_profiles'/case;dest.mkdir(exist_ok=True,parents=True)
    with np.load(stationary/'baseline.npz') as a: masks=a['rois'];img=a['chi'].ravel()
    for ri,rm in enumerate(masks,1):
        out=dest/f'roi_{ri:02}.json'
        if out.exists():continue
        tic=time.perf_counter();cpu=time.process_time();mask=rm.ravel();hi=np.full(m.n,4.5);hi[mask]=0.;lo=np.zeros(m.n)
        seeds=[img]; rec=[]
        for stage in (stationary,base):
            for nm,key in [('b1_predictions',f'x_roi_{ri:02}'),('b1_controls','physical_chi')]:
                f=stage/(nm+'.npz')
                if f.exists():
                    with np.load(f) as a:
                        if key in a:seeds.append(a[key].ravel())
        physical=P/'results/stable_profiles'/case/f'roi_{ri:02}.npz'
        if physical.exists():
            with np.load(physical) as a:seeds.append(a['material'].ravel())
        seeds=[np.minimum(hi,np.maximum(lo,v)) for v in seeds];x=min(seeds,key=lambda v:m.forward(v)['residual_ratio'])
        lam=.002;target=.60; initial=m.forward(x)['residual_ratio'];calls=0
        for it in range(72):
            st=m.forward(x);old=float(st['residual_ratio'])
            if old<target:break
            receiver=m.solve(x,m.SH,adjoint=True).conj().T
            J=(receiver[:,None,:]*st['field'].T[None,:,:]).reshape(m.Y.size,m.n)
            A=np.r_[J.real,J.imag];err=(m.Y-st['pred']).ravel();b=np.r_[err.real,err.imag]
            corr=A.T@b; active=(hi>lo)&(((x>lo+1e-8)&(x<hi-1e-8))|((x<=lo+1e-8)&(corr>0))|((x>=hi-1e-8)&(corr<0)))
            if not np.any(active):break
            scale=max(float(np.linalg.norm(A[:,active],ord='fro')**2/min(A.shape[0],active.sum())),1e-17)
            for inner in range(5):
                inds=np.flatnonzero(active)
                if not len(inds):break
                H=A[:,inds];d=H.T@solve(H@H.T+lam*scale*np.eye(len(b)),b,assume_a='pos'); bad=((x[inds]<=lo[inds]+1e-8)&(d<0))|((x[inds]>=hi[inds]-1e-8)&(d>0))
                if not bad.any() or inner==4:break
                active[inds[bad]]=False
            if not len(inds):break
            step=np.zeros_like(x);step[inds]=d; step*=min(1.,1.0/max(abs(d).max(),1e-16))
            ok=False
            for fraction in (1.,.5,.25,.125,.0625):
                z=np.clip(x+fraction*step,lo,hi);ratio=float(m.forward(z)['residual_ratio']);calls+=1
                if ratio<old*(1-1e-7):
                    x=z;ok=True
                    lam=max(1e-6,lam*.65) if fraction>=.5 else min(.5,lam*1.8);break
            if not ok:lam=min(1e3,lam*8)
            rec.append(dict(iteration=it+1,starting_ratio=old,trial_ratio=ratio,damping_after=lam,accepted=ok,free=int(len(inds))))
            if not ok and lam>=1e3:break
        st=m.forward(x);val=float(st['residual_ratio'])
        np.savez_compressed(out.with_suffix('.npz'),material=x,roi=mask,current=st['j'])
        out.write_text(json.dumps(dict(case=case,roi=ri,cells=int(mask.sum()),initial_ratio=initial,fit_ratio=val,cap_or_attained_peak=0 if val<1 else None,strict_check='see independent inner receipt',n_trials=calls,optimization=rec,costs=dict(m.cost),seconds_wall=time.perf_counter()-tic,seconds_cpu=time.process_time()-cpu),indent=2))
        print(case,ri,initial,val,len(rec),round(time.perf_counter()-tic,2),flush=True)
if __name__=='__main__':
    arg=argparse.ArgumentParser();arg.add_argument('cases',nargs='+');a=arg.parse_args()
    with threadpool_limits(1):
        for c in a.cases:fit(c)
