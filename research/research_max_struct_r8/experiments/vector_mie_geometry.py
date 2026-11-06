from pathlib import Path
import sys,os,json,time,resource
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r7.src.vector_inverse import VectorInverse
from research.research_max_struct_r1.experiments.qualify_backend import mie
OUT=ROOT/'RESULTS/vector3d/mie_geometry';OUT.mkdir(exist_ok=True,parents=True)
def main():
    t=time.perf_counter();c=time.process_time();radius=.145;k=2*np.pi;theta=np.linspace(0,np.pi,181);dirs=np.c_[np.sin(theta),theta*0,np.cos(theta)];_,s2,_=mie(radius,k,3,np.cos(theta));ref=1j/k*s2[:,None]*np.c_[np.cos(theta),theta*0,-np.sin(theta)];rows=[]
    # Fixed physical change: average electrostatic sphere polarizability, not chi.
    # MG interface mixing: chi_eff=3 f chi/(3+(1-f)chi). No fitted parameters.
    a=(np.arange(24)+.5)/24-.5;sub=np.stack(np.meshgrid(a,a,a,indexing='ij'),-1).reshape(-1,3)
    for ns in [5,6,8]:
        m=VectorInverse(ns,.4/ns);frac=np.array([np.mean(np.sum((p+m.h*sub)**2,axis=1)<=radius**2) for p in m.points]);chi=6*frac/(5-2*frac);m.reset_counts();_,J,_,lin=m.physical(chi);phase=np.exp(-1j*k*dirs@m.points.T);dip=np.einsum('rn,ncp->rcp',phase,m.h**3*J.reshape(-1,3,6));F=k*k/(4*np.pi)*(dip-dirs[:,:,None]*np.einsum('rc,rcp->rp',dirs,dip)[:,None,:]);far=F[:,:,0];row=dict(nside=ns,N=m.N,h=m.h,effective_volume=float(frac.sum()*m.h**3),exact_volume=4*np.pi*radius**3/3,mie_far_vector_relative_error=float(la.norm(far-ref)/la.norm(ref)),linear_residual=lin,counts=m.counts,setup=m.setup);rows.append(row);np.savez_compressed(OUT/f'grid{ns}.npz',points=m.points,fraction=frac,chi=chi,far=far,mie_far=ref,theta=theta)
    rec=dict(status='COMPLETED',rule='24³ fixed midpoint fractional geometry with electrostatic polarizability-preserving Maxwell-Garnett interface mixing chi_eff=3fchi/(3+(1-f)chi); no volume rescale/no fit/no scan',scope='independent analytic Mie geometry qualification only; not exact dielectric cut-cell homogenization; inversion model unchanged',rows=rows,wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);(OUT/'RECEIPT.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
