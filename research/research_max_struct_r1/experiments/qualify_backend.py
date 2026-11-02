"""One frozen sphere and one refinement; independent analytic Mie control."""
import json,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn,spherical_yn
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from vector_backend import TensorDDA,sphere_grid


def mie(radius,k,epsilon,mu):
    """Analytic Riccati-Bessel sphere, h_n^(1), not dipole discretization.

    Input mu is cos(scattering angle); incident k=z, E=x.
    Return conventional dimensionless S1,S2 and exact scattering cross section.
    """
    x=k*radius;m=np.sqrt(complex(epsilon));ns=np.arange(1,13)
    def ric(z):
        j=spherical_jn(ns,z);jp=spherical_jn(ns,z,True)
        y=spherical_yn(ns,z);yp=spherical_yn(ns,z,True)
        return z*j,j+z*jp,z*(j+1j*y),j+1j*y+z*(jp+1j*yp)
    px,dpx,xx,dxx=ric(x);pm,dpm,_,_=ric(m*x)
    a=(m*pm*dpx-px*dpm)/(m*pm*dxx-xx*dpm)
    b=(pm*dpx-m*px*dpm)/(pm*dxx-m*xx*dpm)
    S1=np.zeros_like(mu,dtype=complex);S2=S1.copy()
    prev=np.zeros_like(mu);pn=np.ones_like(mu)
    for n,an,bn in zip(ns,a,b):
        tau=n*mu*pn-(n+1)*prev;c=(2*n+1)/(n*(n+1))
        S1+=c*(an*pn+bn*tau);S2+=c*(an*tau+bn*pn)
        nxt=((2*n+1)*mu*pn-(n+1)*prev)/n
        prev,pn=pn,nxt
    sigma=2*np.pi/k**2*np.sum((2*ns+1)*(abs(a)**2+abs(b)**2))
    return S1,S2,float(sigma)


def main():
    out=ROOT/'results/backend_qualification';out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();cpu=time.process_time();rows=[]
    summary=dict(status='RUNNING',radius_lambda=.12,epsilon=3,mu=1,k=2*np.pi,
                 grid_h_lambda=[.04,.03],LU=0,RHS=0,arbitrary_acceptance_gate=False,
                 M_Hcurl_status='NOT_RUN',failure_costs=[])
    try:
        k=2*np.pi;theta=np.linspace(0,np.pi,181)
        dirs=np.column_stack([np.sin(theta),np.zeros_like(theta),np.cos(theta)])
        s1,s2,sigma=mie(.12,k,3,np.cos(theta))
        etheta=np.column_stack([np.cos(theta),np.zeros_like(theta),-np.sin(theta)])
        ref=1j/k*s2[:,None]*etheta
        mu,w=np.polynomial.legendre.leggauss(24);phi=np.arange(32)*2*np.pi/32
        mm,pp=np.meshgrid(mu,phi,indexing='ij')
        quad=np.stack([np.sqrt(1-mm*mm)*np.cos(pp),np.sqrt(1-mm*mm)*np.sin(pp),mm],axis=-1).reshape(-1,3)
        for h in (.04,.03):
            p,v=sphere_grid(.12,h);model=TensorDDA(p,v,k)
            eps=np.broadcast_to(3*np.eye(3),(len(p),3,3)).copy()
            state=model.solve(eps,[[0,0,1],[0,0,1]],[[1,0,0],[0,1,0]],out/f'h_{h:.2f}_solve_receipt.json')
            summary['LU']+=state.receipt['LU'];summary['RHS']+=state.receipt['RHS']
            f=state.far(dirs)['E'].reshape(-1,3,2)[:,:,0]
            q=state.far(quad)['E'].reshape(24,32,3,2)
            cross=float(np.sum(w[:,None]*np.sum(abs(q[:,:,:,0])**2,axis=2))*2*np.pi/32)
            # A near receiver converges to the separately evaluated far operator.
            near=state.near(1000*dirs[::30])
            far_recovered=near['E_scattered'].reshape(-1,3,2)[:,:,0]*1000*np.exp(-1j*k*1000)
            ff=f[::30]
            row=dict(h=h,N=len(p),volume=len(p)*v,analytic_volume=4*np.pi*.12**3/3,
                abs_mie_far_vector_relative_error=float(la.norm(f-ref)/la.norm(ref)),
                scattering_cross_section=cross,mie_scattering_cross_section=sigma,
                scattering_cross_section_relative_error=float(abs(cross-sigma)/sigma),
                near_at_1000lambda_vs_far_relative_difference=float(la.norm(far_recovered-ff)/la.norm(ff)),
                far_transverse_relative_error=float(la.norm(np.einsum('rc,rc->r',dirs,f))/la.norm(f)),
                solve_receipt=state.receipt,
                H_status='analytic_curl_point_dipoles_with_symmetric_cell_zero_self',
                continuum_Hcurl_qualification='NOT_RUN')
            rows.append(row)
            np.savez_compressed(out/f'h_{h:.2f}_fields.npz',points=p,E=state.E,H=state.H,J=state.J,
                 p=state.p,j=state.j,E_exc=state.E_exc,theta=theta,far=f,mie_far=ref)
            print(json.dumps(row),flush=True)
        summary.update(status='COMPLETED',rows=rows,refinement_far_error_decreased=rows[1]['abs_mie_far_vector_relative_error']<rows[0]['abs_mie_far_vector_relative_error'],
                       refinement_cross_section_error_decreased=rows[1]['scattering_cross_section_relative_error']<rows[0]['scattering_cross_section_relative_error'],
                       qualification='BOUNDED_DDA_QUALIFICATION_ONLY_NO_CONTINUUM_CERTIFICATE')
    except Exception as exc:
        summary.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),rows=rows)
        summary['failure_costs'].append(dict(error=repr(exc),wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu))
        raise
    finally:
        summary.update(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
        (out/'QUALIFICATION.json').write_text(json.dumps(summary,indent=2,allow_nan=False))


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
