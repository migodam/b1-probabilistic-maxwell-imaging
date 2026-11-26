"""Whole-domain two-sided forward verification of data-only native images.
Original FFT symbol; receiver-oriented primal and dual residuals.
"""
import argparse,json,sys,time
from pathlib import Path
from verified_wave import P,VerifiedWave,fro,up,rd,opnorm
import numpy as np
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(P/'imaging'))
from shared import Model
def one(g,m,x,j,allow_cached=None):
    inv,info=allow_cached if allow_cached else g.inverse(x)
    Y=m.Y;r=m.radius;sup=float(x.max())
    # A posteriori use S A^{-1}=D^* + (S-D^* A) A^{-1}.
    D=m.solve(x,m.SH,adjoint=True,rtol=2e-12)
    # Exact S rows, A^*D = D - G^* X D. Reuse same verified unitary transform.
    XD=x[:,None]*D
    arr=np.zeros((XD.shape[1],g.p,g.p),complex)
    arr[:,:g.side,:g.side]=XD.T.reshape(-1,g.side,g.side)
    e=g.transform(g.transform(arr)*g.sym.conj(),True)[:,:g.side,:g.side].reshape(XD.shape[1],-1).T
    GHd=e+g.R.conj().T@(g.L.conj().T@XD)
    Er=D-GHd-g.S.conj().T
    adjerr=up(fro(Er)+g.full_green_err*fro(XD)+rd(4)*(fro(D)+g.gnorm*sup*fro(D)+g.snorm))
    gain=up(opnorm(D)+inv*adjerr)
    J=m.solve(x,x[:,None]*g.B,rtol=2e-12)
    gj=g.apply(J,True); b=x[:,None]*g.B
    res=J-x[:,None]*gj-b
    primal=up(fro(res)+sup*g.full_green_err*fro(J)+rd(5)*(fro(J)+sup*fro(gj)+sup*fro(g.B)))
    F=g.S@J; yf=up(gain*primal+rd(g.n)*g.snorm*fro(J)+rd(2)*(fro(F)+fro(Y)))
    ratio=up((fro(F-Y)+yf)/r)
    return dict(residual_upper_over_r=ratio,numeric=float(np.linalg.norm(F-Y)/r),material_min=float(x.min()),
                material_max=float(x.max()),inverse=info,observation_gain_upper=gain,adjoint_residual_upper=adjerr,primal_residual_upper=primal,
                data_rounding_solve_upper=yf),J
def run(case,kind):
    dest=P/'results/native_replay'/case;dest.mkdir(exist_ok=True,parents=True)
    m=Model.read(case);src=P/'imaging/results'/case; masks=np.load(src/'stationary/baseline.npz')['rois']
    with np.load(P/'ports/public/operators'/f'{m.operator_id}.npz') as a:g=VerifiedWave(a)
    if kind=='controls':
        files=[(str(b.relative_to(src)),b,'physical_chi','physical_j') for b in (src/'stationary/b1_controls.npz',src/'b1_controls.npz') if b.exists()]
    else:
        wd=P/'results'/kind/case
        files=[(p.stem,p,'material','current') for p in sorted(wd.glob('roi_??.npz'))]
    for label,path,key,jk in files:
        out=dest/(kind+'_'+label.replace('/','_').replace('.npz','')+'.json')
        if out.exists():continue
        tic=time.perf_counter();cpu=time.process_time()
        with np.load(path) as a:x=a[key].ravel();j=a[jk]
        rec,_=one(g,m,x,j)
        rec.update(case=case,file=str(path.relative_to(P)),key=key,precision=128,
                   exact_model='binary64 input Fourier symbol, fixed loading, S,B,Y and radius',
                   roi_peaks=[float(x[mask.ravel()].max()) for mask in masks],wall_s=time.perf_counter()-tic,cpu_s=time.process_time()-cpu)
        out.write_text(json.dumps(rec,indent=2));print(out.name,rec['residual_upper_over_r'],rec['roi_peaks'],round(rec['wall_s'],2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('cases',nargs='+');ap.add_argument('--kind',default='controls');q=ap.parse_args()
    with threadpool_limits(1):
        for c in q.cases:run(c,q.kind)
