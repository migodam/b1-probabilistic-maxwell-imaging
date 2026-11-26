"""A local Gaussian imaging comparator, with uncertainty scale fixed before scoring."""
from pathlib import Path
import argparse,sys,json,time,os
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.linalg import solve
from scipy.fft import fft2,ifft2
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'imaging'))
from shared import Model
def run(case):
    d=P/'imaging/results'/case/'stationary'; f=d/'baseline.npz'
    if not f.exists():return
    out=P/'results/gaussian';out.mkdir(parents=True,exist_ok=True)
    previous=out/(case+'.npz')
    if previous.exists():
        with np.load(previous) as old,np.load(f) as current:
            if 'roi_masks' in old and np.array_equal(old['roi_masks'],current['rois']):return
        past=out/'before_final_mask_audit';past.mkdir(exist_ok=True)
        for ext in ('.npz','.json'):
            src=out/(case+ext);target=past/(case+ext)
            if src.exists() and not target.exists():src.rename(target)
    t=time.perf_counter();cpu=time.process_time();m=Model.read(case);a=np.load(f);x=a['chi'].ravel();masks=a['rois']
    st=m.forward(x);p=m.solve(x,m.SH,adjoint=True).conj().T
    jac=(p[:,None,:]*st['field'].T[None,:,:]).reshape(m.Y.size,m.n)
    A=np.r_[jac.real,jac.imag];res=(m.Y-st['pred']).ravel();b=np.r_[res.real,res.imag]
    sigma=m.radius/np.sqrt(len(b))
    freq=np.fft.fftfreq(m.side,d=2/m.side)*2*np.pi
    k2=freq[:,None]**2+freq[None,:]**2
    # Unit pixel prior variance, a fixed physical 0.12 m smoothing length.
    sym=np.exp(-.5*(.12**2)*k2);sym*=.75/np.sqrt(np.mean(sym**2))
    def filt(z):
        z=np.asarray(z); shape=z.shape
        return ifft2(fft2(z.reshape(-1,m.side,m.side),axes=(1,2))*sym,axes=(1,2)).real.reshape(shape)
    AL=filt(A);C=AL@AL.T+sigma*sigma*np.eye(len(b))
    correction=filt(solve(C,AL,assume_a='pos',check_finite=False))
    rng=np.random.default_rng(231771);stats=[]
    for _ in range(16):
        z=rng.standard_normal((32,m.n));eps=sigma*rng.standard_normal((32,len(b)))
        draw=np.clip(x+filt(z)+(b-z@AL.T-eps)@correction,0,4.5)
        stats.append(np.column_stack([np.max(draw[:,mk.ravel()],axis=1) for mk in masks]))
    vals=np.vstack(stats);quant=np.quantile(vals,[.025,.5,.975],axis=0).T
    rec=dict(case=case,intervals=quant.tolist(),target='peak contrast of each stationary blind ROI',
             prior_marginal_sd=.75,prior_smoothing_m=.12,draws=len(vals),clipped_to=[0,4.5],
             scope='local linear Gaussian model at CIE endpoint; finite-draw credible interval, not a simultaneous deterministic outer bound',
             training_residual_over_r=st['residual_ratio'],seconds_wall=time.perf_counter()-t,
             seconds_cpu=time.process_time()-cpu,physics=dict(m.cost))
    (out/(case+'.json')).write_text(json.dumps(rec,indent=2))
    np.savez_compressed(out/(case+'.npz'),peak_draws=vals,quantiles=quant,roi_masks=masks)
    print(case,len(masks),quant.tolist(),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('cases',nargs='+');q=ap.parse_args()
    with threadpool_limits(1):
        for case in q.cases:run(case)
