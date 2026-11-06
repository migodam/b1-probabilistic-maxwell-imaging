"""One root-authorized grid16 refinement; no material/scale fitting.

Run only after the finite inversion queue exits (max two numerical children).
Existing TensorDDA arbitrary occupied-cell API removes exactly chi=0 rows.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
WALL,CPU=time.perf_counter(),time.process_time()
import sys,json,resource,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy import linalg as la
from scipy.special import spherical_jn,spherical_yn
from threadpoolctl import threadpool_limits
from collections import Counter
from research.research_max_struct_r9.src.vector_extension import paid_density_state
from research.research_max_struct_r1.src.vector_backend import TensorDDA,plane_waves
from research.research_max_struct_r1.experiments.qualify_backend import mie
OUT=ROOT/'RESULTS/vector_extension/validation_refinement16'

def diagnostics(far,reference,theta):
    alpha=np.vdot(reference,far)/np.vdot(reference,reference)
    return dict(relative_complex_vector_error=float(la.norm(far-reference)/la.norm(reference)),
        diagnostic_unapplied_best_complex_gain=[float(alpha.real),float(alpha.imag)],
        diagnostic_gain_amplitude=float(abs(alpha)),diagnostic_gain_phase_degrees=float(np.angle(alpha,deg=True)),
        relative_shape_residual_after_unapplied_gain=float(la.norm(far-alpha*reference)/la.norm(reference)),
        theta_degrees=(theta[::30]*180/np.pi).tolist(),
        theta_component_samples=[[float(z.real),float(z.imag)] for z in np.sum(far*np.c_[np.cos(theta),theta*0,-np.sin(theta)],axis=1)[::30]],
        reference_theta_component_samples=[[float(z.real),float(z.imag)] for z in np.sum(reference*np.c_[np.cos(theta),theta*0,-np.sin(theta)],axis=1)[::30]])

def run():
    OUT.mkdir(exist_ok=False);receipt=dict(status='RUNNING',authorization='root: one grid16 MG same sphere after large grid8/10 errors; no sweep',
        epsilon_r=9,chi=8,radius=.145,domain_width=.4,nside=16,subcells_per_axis=24,
        scope='same finite CM+RR DDA, only exactly zero-contrast background eliminated through existing arbitrary-point API; not a new backend',
        no_material_or_amplitude_fit=True,no_inversion_change=True,phaseI_modified=False,counts={})
    cost=Counter()
    try:
        old=[]
        for n in (8,10):
            with np.load(ROOT/f'RESULTS/vector_extension/validation/grid{n}.npz') as z:
                old.append(dict(nside=n,**diagnostics(z['far'],z['Mie_far'],z['theta'])))
        receipt['unchanged_8_10_angular_diagnostics']=old
        n=16;width=.4;h=width/n;k=2*np.pi;r=.145;eps=9;chi=eps-1
        ijk=np.indices((n,)*3).reshape(3,-1).T;points=(ijk-(n-1)/2)*h
        nearest=np.maximum(abs(points)-h/2,0);possible=np.sum(nearest**2,axis=1)<=r*r
        a=(np.arange(24)+.5)/24-.5;sub=np.stack(np.meshgrid(a,a,a,indexing='ij'),-1).reshape(-1,3)
        candidates=points[possible]
        fraction=np.array([np.mean(np.sum((p+h*sub)**2,axis=1)<=r*r) for p in candidates]);positive=fraction>0
        points=candidates[positive];fraction=fraction[positive];M=3*len(points)
        receipt['preassembly_memory_estimate']=dict(full16_cells=4096,active_cells=len(points),vector_dimension=M,
            complex_dense_matrix_bytes=16*M*M,conservative_eight_matrix_bytes=8*16*M*M,
            full16_single_matrix_bytes=16*(3*4096)**2,
            basis='Goff/G/K/LU plus conservative assembly temporaries; estimate not peak measurement')
        receipt['geometry_volume']=float(h**3*fraction.sum());receipt['analytic_volume']=4*np.pi*r**3/3
        # Same 24³ midpoint fraction, MG static polarizability mixing as grids8/10.
        x=3*fraction*chi/(3+(1-fraction)*chi)
        backend=TensorDDA(points,h**3);inc=plane_waves(points,k,[[0,0,1]],[[1,0,0]])
        state=paid_density_state(backend.G,inc,x,None,cost)
        theta=np.linspace(0,np.pi,181);directions=np.c_[np.sin(theta),theta*0,np.cos(theta)]
        s1,s2,sigma=mie(r,k,eps,np.cos(theta));ref=1j/k*s2[:,None]*np.c_[np.cos(theta),theta*0,-np.sin(theta)]
        def far_at(d):
            dip=np.exp(-1j*k*d@points.T)@(h**3*state['J'].reshape(-1,3))
            return k*k/(4*np.pi)*(dip-d*np.sum(d*dip,axis=1)[:,None])
        far=far_at(directions);receipt['grid16_angular_diagnostics']=diagnostics(far,ref,theta)
        mu,w=np.polynomial.legendre.leggauss(24);phi=np.arange(32)*2*np.pi/32;mm,pp=np.meshgrid(mu,phi,indexing='ij')
        directions_quad=np.stack([np.sqrt(1-mm*mm)*np.cos(pp),np.sqrt(1-mm*mm)*np.sin(pp),mm],axis=-1).reshape(-1,3)
        quad=far_at(directions_quad).reshape(24,32,3);cross=float(np.sum(w[:,None]*np.sum(abs(quad)**2,axis=2))*2*np.pi/32)
        receipt.update(scattering_cross_section=cross,mie_scattering_cross_section=sigma,
            scattering_cross_section_relative_error=abs(cross-sigma)/sigma,linear_relative_residual=state['linear_relative_residual'],
            internal_wavelength=1/np.sqrt(eps),cells_per_internal_wavelength=1/np.sqrt(eps)/h,
            sphere_external_size_parameter=k*r,sphere_internal_size_parameter=k*r*np.sqrt(eps))
        # Diagnostic exact Mie partial-wave weights, not fitted dipole parameters.
        ns=np.arange(1,13);xx=k*r;m=np.sqrt(complex(eps))
        def ric(z):
            j=spherical_jn(ns,z);jp=spherical_jn(ns,z,True);y=spherical_yn(ns,z);yp=spherical_yn(ns,z,True)
            return z*j,j+z*jp,z*(j+1j*y),j+1j*y+z*(jp+1j*yp)
        px,dpx,hx,dhx=ric(xx);pm,dpm,_,_=ric(m*xx)
        aa=(m*pm*dpx-px*dpm)/(m*pm*dhx-hx*dpm);bb=(pm*dpx-m*px*dpm)/(pm*dhx-m*hx*dpm)
        weight=(2*ns+1)*(abs(aa)**2+abs(bb)**2)
        receipt['exact_mie_partial_wave_diagnostic']=[dict(n=int(i),electric_abs=float(abs(av)),magnetic_abs=float(abs(bv)),
            scattering_fraction=float(ww/weight.sum())) for i,av,bv,ww in zip(ns[:4],aa[:4],bb[:4],weight[:4])]
        sensitivity=[]
        for factor in (.95,.975,1.,1.025,1.05):
            kk=k*factor;xx=kk*r;px,dpx,hx,dhx=ric(xx);pm,dpm,_,_=ric(m*xx)
            an=(m*pm*dpx-px*dpm)/(m*pm*dhx-hx*dpm);bn=(pm*dpx-m*px*dpm)/(pm*dhx-m*hx*dpm)
            _,_,total=mie(r,kk,eps,np.array([0.]))
            sensitivity.append(dict(k_factor=factor,k=kk,external_size_parameter=xx,internal_size_parameter=float(xx*m.real),
                magnetic_dipole_abs=float(abs(bn[0])),electric_dipole_abs=float(abs(an[0])),mie_scattering_cross_section=total,
                scope='five predeclared analytic-only frequencies, no DDA solves/fit or change to inversion data'))
        receipt['analytic_only_five_frequency_sensitivity']=sensitivity
        np.savez_compressed(OUT/'grid16.npz',points=points,fraction=fraction,chi=x,theta=theta,far=far,Mie_far=ref)
        receipt['status']='COMPLETED';receipt['continuum_accuracy']='NOT_CERTIFIED; interpret measured refinement and all prior failures'
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
    receipt.update(counts=dict(cost),CPU_seconds=time.process_time()-CPU,wall_seconds=time.perf_counter()-WALL,
        process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (OUT/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt),flush=True)
    return 0 if receipt['status']=='COMPLETED' else 1

if __name__=='__main__':
    with threadpool_limits(limits=1):raise SystemExit(run())
