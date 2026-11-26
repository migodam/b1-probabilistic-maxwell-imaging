"""Continue the one primary 32-grid case where even an unrestricted public fit failed."""
from pathlib import Path
import json,sys,time
from gn_profiles import P,Model,np,threadpool_limits
from run import physical_fit
def run():
    c='n32_c07';init=P/'imaging/results'/c/'stationary';dest=P/'results/gap_profiles'/c;dest.mkdir(parents=True,exist_ok=True)
    a=np.load(init/'baseline.npz'); masks=a['rois'];m=Model.read(c);hi=np.full(m.n,4.5)
    x=np.load(init/'b1_controls.npz')['physical_chi'].ravel();trace=[];tic=time.perf_counter();cpu=time.process_time()
    for k in range(3):
        x,j,rec=physical_fit(m,x,hi,max_updates=450,target=.55)
        trace.append(dict(iteration=k+1,rho=rec['residual_over_radius'],n=rec['accepted_updates']))
        np.savez_compressed(dest/'unrestricted.npz',material=x,current=j)
        if rec['residual_over_radius']<.95:break
    (dest/'unrestricted.json').write_text(json.dumps(dict(trace=trace,cost=dict(m.cost),wall=time.perf_counter()-tic,cpu=time.process_time()-cpu),indent=2))
    for i,mask in enumerate(masks,1):
        start=x.copy();mask=mask.ravel();table=[];tic=time.perf_counter();cpu=time.process_time();best=None;val=None
        for cap in (4.5,3.,2.,1.,.5,.2,0.):
            up=hi.copy();up[mask]=cap
            z,j,rec=physical_fit(m,np.minimum(start,up),up,max_updates=300,target=.65)
            rho=rec['residual_over_radius'];table.append(dict(cap=cap,rho=rho,n=rec['accepted_updates']))
            if rho<.99:
                best=float(z[mask].max());val=rho;np.savez_compressed(dest/f'roi_{i:02}.npz',material=z,roi=mask,current=j)
            if rho>2:break
            start=z
        (dest/f'roi_{i:02}.json').write_text(json.dumps(dict(case=c,roi=i,cells=int(mask.sum()),fit_ratio=val,cap_or_attained_peak=best,trace=table,seconds_wall=time.perf_counter()-tic,seconds_cpu=time.process_time()-cpu),indent=2))
        print(c,i,best,val,flush=True)
if __name__=='__main__':
    with threadpool_limits(1):run()
