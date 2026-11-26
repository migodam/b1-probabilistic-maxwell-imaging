"""Try a spatially extended alternative in the same blind region as a zero witness."""
import argparse,heapq,json,os,sys,time
from pathlib import Path
for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from verified_wave import P,VerifiedWave,up,rd,fro
sys.path.insert(0,str(P/'imaging'))
from shared import Model
def patch_of(mask):
    ds=distance_transform_edt(mask);n=mask.shape[0];i=int(np.argmax(ds));heap=[(-ds.flat[i],i)]
    visited={i};picked=[]
    while heap and len(picked)<max(1,int(np.ceil(mask.sum()/4))):
        _,j=heapq.heappop(heap);picked.append(j);y,x=divmod(j,n)
        for v,w in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
            z=v*n+w
            if 0<=v<n and 0<=w<n and mask[v,w] and z not in visited:
                visited.add(z);heapq.heappush(heap,(-ds[v,w],z))
    return np.array(picked)
class Enough(Exception):pass
def run(case,kind='stable_profiles'):
    wd=P/'results'/kind/case
    avail=[]
    for p in wd.glob('roi_??.json'):
        q=json.loads(p.read_text())
        if q['cap_or_attained_peak']==0 and q['fit_ratio']<.98:avail.append((q['cells'],p))
    if not avail:return
    _,f=max(avail);D=P/'results/extended';D.mkdir(exist_ok=True);out=D/(case+'.json')
    if out.exists():return
    t=time.perf_counter();cpu=time.process_time()
    m=Model.read(case);a=np.load(f.with_suffix('.npz'));x=a['material'];roi=a['roi'].reshape(m.shape)
    patch=patch_of(roi);lo=np.zeros(m.n);lo[patch]=2.;hi=np.full(m.n,4.5)
    seeds=[x,np.load(P/'imaging/results'/case/'stationary/baseline.npz')['chi'].ravel()]
    ctl=P/'imaging/results'/case/'b1_controls.npz'
    if ctl.exists():seeds.append(np.load(ctl)['physical_chi'].ravel())
    x=min([np.clip(z,lo,hi) for z in seeds],key=lambda v:m.forward(v)['residual_ratio'])
    best=np.inf;zbest=x.copy(); calls=0
    def fg(z):
        nonlocal calls,best,zbest
        loss,grad=m.data_grad(z);r=m.forward(z)['residual_ratio'];calls+=1
        if r<best:best=r;zbest=z.copy()
        if r<.83:raise Enough
        return loss,grad
    try:
        minimize(fg,x,method='L-BFGS-B',jac=True,bounds=np.c_[lo,hi],options=dict(maxiter=155,ftol=1e-14,gtol=1e-10,maxls=20))
    except Enough:pass
    mat=m.forward(zbest);cur=mat['j'];proof={};ru=None
    if best<1:
        with np.load(P/'ports/public/operators'/f'{m.operator_id}.npz') as op:g=VerifiedWave(op)
        inv,proof=g.inverse(zbest);ej=g.solve_error(zbest,cur,zbest[:,None]*g.B,up(rd(1)*max(zbest)*fro(g.B)),inv)
        pred=g.S@cur;ru=up((fro(pred-m.Y)+g.snorm*ej+rd(g.n)*g.snorm*fro(cur)+rd(2)*(fro(pred)+fro(m.Y)))/m.radius)
    np.savez_compressed(D/(case+'.npz'),material=zbest,roi=roi,patch=patch,current=cur)
    res=dict(case=case,roi_index=int(f.stem[-2:]),cells=int(roi.sum()),changed_patch_cells=len(patch),peak_floor=2.,
             patch_area_cm2=len(patch)*(200/m.side)**2,roi_average=float(zbest[roi.ravel()].mean()),
             minimum_on_patch=float(min(zbest[patch])),residual_over_r_numeric=best,upper_residual_over_r=ru,
             zero_file=str(f.relative_to(P)),online_selection='largest verified-or-pending zero ROI by cells; connected 25% interior patch',
             inverse=proof,cost=dict(m.cost),evaluations=calls,wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu)
    out.write_text(json.dumps(res,indent=2));print(case,res['roi_index'],len(patch),ru,best,flush=True)
if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('cases',nargs='+');pa.add_argument('--kind',default='stable_profiles');args=pa.parse_args()
    with threadpool_limits(1):
        for c in args.cases:run(c,args.kind)
