"""One physically motivated geometry correction: fixed24³ subcell volume fractions.
Validation only; inverse data, support, bounds, and predictions remain untouched.
Volume averaged chi=f*(epssphere-1) is a pulse/interface approximation, not an
exact cut-cell Maxwell interface treatment or continuum certificate.
"""
from pathlib import Path
import os,sys,json,time,resource
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r7.src.vector_inverse import VectorInverse
from research.research_max_struct_r1.experiments.qualify_backend import mie
OUT=ROOT/'RESULTS/vector_fractional_validation';OUT.mkdir(exist_ok=True)
with threadpool_limits(limits=1):
    tick=time.perf_counter();cpu=time.process_time();theta=np.linspace(0,np.pi,181);dirs=np.column_stack([np.sin(theta),np.zeros_like(theta),np.cos(theta)])
    _,s2,sigma=mie(.145,2*np.pi,3,np.cos(theta));ref=1j/(2*np.pi)*s2[:,None]*np.column_stack([np.cos(theta),np.zeros_like(theta),-np.sin(theta)])
    a=(np.arange(24)+.5)/24-.5;sub=np.stack(np.meshgrid(a,a,a,indexing='ij'),axis=-1).reshape(-1,3);rows=[]
    for ns in [5,6]:
        model=VectorInverse(ns,.4/ns);f=np.mean(np.sum((model.points[:,None,:]+model.h*sub[None,:,:])**2,axis=-1)<=.145**2,axis=1);chi=2*f
        model.reset_counts();Y,J,E,lin=model.physical(chi);p=model.h**3*J.reshape(-1,3,6)
        phase=np.exp(-1j*model.backend.k*dirs@model.points.T);dip=np.einsum('rn,ncp->rcp',phase,p)
        F=model.backend.k**2/(4*np.pi)*(dip-dirs[:,:,None]*np.einsum('rc,rcp->rp',dirs,dip)[:,None,:]);far=F[:,:,0]
        row=dict(nside=ns,N=model.N,h=model.h,fractional_cells=int(np.sum((f>0)&(f<1))),effective_volume=float(np.sum(f)*model.h**3),analytic_volume=4*np.pi*.145**3/3,mie_far_vector_relative_error=float(la.norm(far-ref)/la.norm(ref)),linear_relative_residual=lin,physical_counts=model.counts,setup=model.setup)
        rows.append(row);np.savez_compressed(OUT/('grid%d.npz'%ns),points=model.points,fraction=f,chi=chi,far=far,mie_far=ref,theta=theta)
    rec=dict(status='COMPLETED',subcell_rule='fixed24³ equally spaced midpoint samples in every cell; same rule both grids; no scan',dielectric_rule='volume averaged contrast chi=f*(3-1), isotropic CM+RR response retained',inverse_changes='NONE',rows=rows,refinement_error_decreased=rows[1]['mie_far_vector_relative_error']<rows[0]['mie_far_vector_relative_error'],scope='geometry aliasing diagnostic only; interface approximation remains, no continuum imaging advantage certificate',wall_s=time.perf_counter()-tick,cpu_s=time.process_time()-cpu,process_peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    (OUT/'RECEIPT.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec))
