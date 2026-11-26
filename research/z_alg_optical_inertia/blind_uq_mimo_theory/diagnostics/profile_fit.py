"""Trace attainable ROI peaks using only public data and frozen blind ROIs."""
from pathlib import Path
import argparse, json, os, sys, time
for name in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(P/'imaging'))
from shared import Model
from run import physical_fit
U=4.5

def fit(case,ids=None,stationary=False):
    folder=P/'imaging/results'/case
    m=Model.read(case); images=folder/'stationary' if stationary else folder
    if not (images/'baseline.npz').is_file():return
    x0=np.load(images/'baseline.npz')
    masks=x0['rois']; originals=x0['chi'].ravel()
    rescue=folder/'b1_predictions.npz'; initial_candidates=[originals]
    full=folder/'b1_controls.npz'
    if (images/'b1_controls.npz').exists():full=images/'b1_controls.npz'
    if full.exists(): initial_candidates.append(np.load(full)['physical_chi'].ravel())
    allrows=[]
    for ri,mask in enumerate(masks,1):
        if ids and ri not in ids: continue
        t=time.perf_counter(); cpu=time.process_time()
        dest=P/'results'/('stable_profiles' if stationary else 'profiles')/case; dest.mkdir(parents=True,exist_ok=True)
        meta=dest/f'roi_{ri:02}.json'
        if meta.exists(): continue
        mask=mask.ravel()
        constrained=originals.copy()
        if rescue.exists() and not stationary:
            name=f'x_roi_{ri:02}'
            with np.load(rescue) as inp:
                if name in inp: constrained=inp[name].ravel()
        best=constrained; hist=[]; least=None; accepted=None
        # Direct removal with physical adjustable exterior. Warm-start from the best
        # clipped data-only proposal; never use held readings to choose an image.
        upper=np.full(m.n,U); upper[mask]=0
        probes=[np.minimum(v,upper) for v in [constrained,*initial_candidates]]
        best=min(probes,key=lambda v:m.forward(v)['residual_ratio'])
        base=m.forward(best)['residual_ratio']
        if base>=.45:
            best,cur,rec=physical_fit(m,best,upper,max_updates=210,target=.45)
        else: rec={'residual_over_radius':float(base),'accepted_updates':0}
        hist.append(dict(cap=0,**{k:rec[k] for k in ('residual_over_radius','accepted_updates')}))
        probe=best.copy();cur=m.forward(probe)['j'];res=m.forward(probe)['residual_ratio']
        if res<.985: least=0.;accepted=(probe.copy(),cur.copy(),res)
        # Continuous cap relaxation gives an attainable *inner* endpoint if zero
        # remains unresolved. Each successful step starts from the last feasible map.
        if accepted is None:
            seed=min(initial_candidates,key=lambda v:m.forward(v)['residual_ratio'])
            if m.forward(seed)['residual_ratio']>.93:
                seed,_,_=physical_fit(m,seed,np.full(m.n,U),max_updates=160,target=.75)
            for cap in (2.25,1.0,.5,.2,0.):
                upper[mask]=cap
                z,j,rec=physical_fit(m,np.minimum(seed,upper),upper,max_updates=95,target=.65)
                val=rec['residual_over_radius']
                hist.append(dict(cap=cap,residual_over_radius=val,accepted_updates=rec['accepted_updates']))
                if val<.99:
                    least=float(z[mask].max());accepted=(z.copy(),j.copy(),val);seed=z
                else:
                    if accepted is None or cap>0:break
        if accepted is None and m.forward(initial_candidates[-1])['residual_ratio']<.99:
            x=initial_candidates[-1];st=m.forward(x)
            accepted=(x.copy(),st['j'].copy(),st['residual_ratio']);least=float(x[mask].max())
        if accepted is not None:
            x,j,ratio=accepted
            np.savez_compressed(dest/f'roi_{ri:02}.npz',material=x,roi=mask,current=j)
        else:ratio=None
        rec=dict(case=case,roi=ri,cells=int(mask.sum()),cap_or_attained_peak=least,fit_ratio=ratio,
                 costs=dict(m.cost),optimization=hist,seconds_wall=time.perf_counter()-t,seconds_cpu=time.process_time()-cpu)
        meta.write_text(json.dumps(rec,indent=2));print(case,ri,least,ratio,round(rec['seconds_wall'],2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('cases',nargs='+');ap.add_argument('--roi',type=int,action='append')
    ap.add_argument('--stationary',action='store_true');args=ap.parse_args()
    with threadpool_limits(1):
        for case in args.cases:fit(case,args.roi,args.stationary)
