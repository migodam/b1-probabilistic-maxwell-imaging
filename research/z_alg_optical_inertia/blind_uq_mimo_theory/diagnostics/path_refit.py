"""Continuation inside the same nonlinear material box, with an archived physical endpoint.
Lower the ROI cap in finite steps. Warm starts are measured-data fits, never true material.
"""
import argparse,json,sys,time,os
from pathlib import Path
for key in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'imaging'))
from shared import Model
from run import physical_fit
def run(case,indices):
    m=Model.read(case); src=P/'imaging/results'/case; stable=src/'stationary'
    dst=P/'results/path_profiles'/case;dst.mkdir(parents=True,exist_ok=True)
    a=np.load(stable/'baseline.npz'); masks=a['rois']
    for i in indices:
        if i>len(masks):continue
        out=dst/f'roi_{i:02}.json';mask=masks[i-1].ravel();seeds=[];trace=[]
        if out.exists():continue
        tic=time.perf_counter();cpu=time.process_time()
        for b in [stable,src]:
            ctl=b/'b1_controls.npz'
            if ctl.exists():seeds.append(np.load(ctl)['physical_chi'].ravel())
        ref=P/'results/stable_profiles'/case/f'roi_{i:02}.npz'
        if ref.exists():seeds.append(np.load(ref)['material'].ravel())
        if m.side==32:seeds.extend([np.zeros(m.n),np.full(m.n,.3)])
        x=min(seeds,key=lambda z:m.forward(z)['residual_ratio']);st=m.forward(x);upper=np.full(m.n,4.5)
        if st['residual_ratio']>.6:
            rc=[]
            for z in seeds:
                z,_,a=physical_fit(m,z,upper,max_updates=300 if m.side==32 else 100,target=.56)
                rc.append((a['residual_over_radius'],z));trace.append(dict(stage='find unrestricted',rho=a['residual_over_radius'],updates=a['accepted_updates']))
                if a['residual_over_radius']<.8:break
            x=min(rc,key=lambda a:a[0])[1]
        bestx=x.copy();st=m.forward(x);bestval=st['residual_ratio'];lastcap=float(x[mask].max())
        if bestval<.98:
            # Save this already-physical anchor even if later continuation is interrupted.
            np.savez_compressed(out.with_suffix('.npz'),material=x,roi=mask,current=st['j'])
        for cap in (3.5,2.75,2.,1.5,1.,.65,.35,.15,0.):
            if cap>lastcap:continue
            upper[mask]=cap;candidate,j,rec=physical_fit(m,np.minimum(x,upper),upper,max_updates=160 if m.side==32 else 68,target=.67)
            ratio=rec['residual_over_radius'];trace.append(dict(cap=cap,rho=ratio,updates=rec['accepted_updates']))
            if ratio<.985:
                bestx=candidate.copy();bestval=ratio;lastcap=float(candidate[mask].max())
                np.savez_compressed(out.with_suffix('.npz'),material=candidate,roi=mask,current=j)
            elif ratio>max(2.,bestval*2):break
            x=candidate
        out.write_text(json.dumps(dict(case=case,roi=i,cells=int(mask.sum()),fit_ratio=bestval,cap_or_attained_peak=lastcap if bestval<.985 else None,trace=trace,
                 costs=dict(m.cost),seconds_wall=time.perf_counter()-tic,seconds_cpu=time.process_time()-cpu),indent=2))
        print(case,i,lastcap,bestval,round(time.perf_counter()-tic,1),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('rois',type=int,nargs='+');a=p.parse_args()
    with threadpool_limits(1):run(a.case,a.rois)
