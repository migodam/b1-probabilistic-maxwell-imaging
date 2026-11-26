"""Bounded Gauss--Newton witness search on the native dressed operator.
The hierarchy restricts the proposals only. Reported data residuals are full wave.
"""
import argparse, json, os, sys, time
from pathlib import Path
for key in ("VECLIB_MAXIMUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from scipy.linalg import solve,svdvals
from scipy.ndimage import zoom
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/"imaging"))
from shared import Model

def solve_case(case, target=.72, rounds=38, mask_path=None, output=None):
    t=time.perf_counter(); cpu=time.process_time()
    m=Model.read(case); side=m.side; radius=m.radius
    yy,xx=np.indices((side,side)); n=m.n
    zero=np.zeros(n,bool) if mask_path is None else np.load(mask_path).reshape(-1).astype(bool)
    levels=[8,16,32] if side==32 else [8,16,32,side]
    store=P/'results'/'fits';store.mkdir(parents=True,exist_ok=True)
    out=store/(output or case); x=np.zeros(n); records=[]
    for level in levels:
        if level>side:continue
        labels=((yy*level//side)*level+xx*level//side).ravel()
        rX=np.bincount(labels,weights=x,minlength=level**2)/np.bincount(labels,minlength=level**2)
        lower=np.zeros_like(rX); upper=np.full_like(rX,4.5)
        trust=.9; damping=.003
        for it in range(rounds):
            x=rX[labels]; x[zero]=0
            J=m.solve(x,x[:,None]*m.B);E=m.B+m.green.apply(J);Y0=m.S@J
            R=m.Y-Y0; err=np.linalg.norm(R); ratio=err/radius
            entry=dict(level=level,it=it,ratio=float(ratio),wall_s=time.perf_counter()-t);records.append(entry)
            if ratio <= target:break
            H=m.solve(x,m.SH,adjoint=True).conj().T
            jac=(H[:,None,:]*E.T[None,:,:]).reshape(R.size,n);jac[:,zero]=0
            # Aggregate a real parameter update per tile in the unchanged EM model.
            order=np.argsort(labels,kind='stable')
            starts=np.flatnonzero(np.r_[True,np.diff(labels[order])!=0])
            jac=np.add.reduceat(jac[:,order],starts,axis=1)
            G=np.r_[jac.real,jac.imag]; y=np.r_[R.real.ravel(),R.imag.ravel()]
            # Simple active-set Levenberg step in data (or parameter) coordinates.
            grad=G.T@y
            active=((rX<1e-10)&(grad<=0))|((rX>4.5-1e-10)&(grad>=0))
            free=~active
            Gn=G[:,free];norm2=np.linalg.norm(Gn,ord='fro')**2/max(1,min(Gn.shape));alpha=max(1e-16,norm2*damping)
            if not np.any(free):break
            if Gn.shape[1]<=Gn.shape[0]:
                d=solve(Gn.T@Gn+alpha*np.eye(Gn.shape[1]),Gn.T@y,assume_a='pos',check_finite=False)
            else:d=Gn.T@solve(Gn@Gn.T+alpha*np.eye(Gn.shape[0]),y,assume_a='pos',check_finite=False)
            delta=np.zeros_like(rX);delta[free]=d
            if max(abs(delta))>trust:delta*=trust/max(abs(delta))
            accepted=False
            for a in (1.,.5,.25,.125,.0625):
                proposal=np.clip(rX+a*delta,lower,upper);trial=proposal[labels];trial[zero]=0
                JT=m.solve(trial,trial[:,None]*m.B)
                val=np.linalg.norm(m.Y-m.S@JT)
                if val < err*(1-1e-8):
                    rX=proposal; trust=min(2,trust*1.22);damping=max(1e-5,damping*.7);accepted=True;break
            if not accepted:
                trust*=.72;damping*=8
                if trust<.009:break
        x=rX[labels];x[zero]=0
        cur=m.solve(x,x[:,None]*m.B);p=m.S@cur;rr=np.linalg.norm(m.Y-p)/radius
        print(case,level,it,rr,round(time.perf_counter()-t,2),flush=True)
        np.savez_compressed(out.with_suffix('.npz'),material=x.reshape(side,side),current=cur,Yhat=p,r=radius)
        if rr<=target:break
    out.with_suffix('.json').write_text(json.dumps(dict(case=case,target=target,method="native full-wave bounded GN from air with successively released spatial tiles",roi_held_zero=int(zero.sum()),residual_over_r=float(rr),wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,cost=dict(m.cost),history=records),indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('cases',nargs='+');ap.add_argument('--rounds',type=int,default=38);ap.add_argument('--target',type=float,default=.72);ap.add_argument('--mask');ap.add_argument('--output')
    args=ap.parse_args()
    with threadpool_limits(1):
        for case in args.cases:solve_case(case,args.target,args.rounds,args.mask,args.output)
