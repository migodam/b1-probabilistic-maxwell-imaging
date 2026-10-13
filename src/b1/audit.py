from __future__ import annotations
import json,time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from .physics import Forward,a13
from .budget import Ledger


def run(root,config):
    root=Path(root); led=Ledger(root/'results',config['physics_new_cap']); f=Forward(config,led,root/'results/forward_cache')
    checks={}; elapsed=time.perf_counter()
    def record(name,ok,**kw): checks[name]=dict(passed=bool(ok),**kw)
    saved={}
    for H in range(3):
        main=f.run(H,2.5,.10,kind='AUDIT_MAIN',cache=False); reference=f.run(H,2.5,.10,kind='AUDIT_RECURSIVE',implementation='recursive',cache=False)
        diff=np.linalg.norm(main[0]-reference[0])/max(np.linalg.norm(main[0]),1e-30)
        saved[H]=main
        record(f'independent_layer_matching_{H}',diff<1e-10,relative_field_error=float(diff),boundary_backward_residual=main[3]['boundary_backward_residual'])
        record(f'passive_power_{H}',main[3]['Qabs']>=-1e-12,**main[3])
    hom=f.run(0,2.5,.10,kind='AUDIT_A13_HOMOGENEOUS',implementation='a13_homogeneous',cache=False)
    diff=np.linalg.norm(saved[0][0]-hom[0])/np.linalg.norm(hom[0]); record('unmodified_A13_coefficients',diff<1e-11,relative_field_error=float(diff))
    for H in range(3):
        vacuum=f.run(H,1.,0.,kind='AUDIT_VACUUM',cache=False)
        record(f'vacuum_{H}',np.linalg.norm(vacuum[0])<1e-12,field_norm=float(np.linalg.norm(vacuum[0])))
    same=f.run(0,2.5,.10,kind='AUDIT_REDUNDANT_INTERFACE',cache=False,radii=[.35,.65,1.],eps_layers=[2.5+.1j]*3)
    diff=np.linalg.norm(same[0]-hom[0])/np.linalg.norm(hom[0]);record('redundant_interfaces_reduce_to_homogeneous',diff<1e-10,relative_error=float(diff))
    for H in [1,2]:
        high=f.run(H,2.5,.10,kind='AUDIT_TRUNCATION',nmax=14,cache=False)
        diff=np.linalg.norm(high[0]-saved[H][0])/np.linalg.norm(high[0]);record(f'multipole_refinement_{H}',diff<1e-10,relative_field_error=float(diff))
    for H in [0,1]:
        small=f.run(H,2.5,.10,kind='AUDIT_RAYLEIGH',k=.001,nmax=4,cache=False)
        e=2.5+.1j
        if H==0: K=(e-1)/(e+2)
        else:
            v=.5**3; ec=1.
            K=((e-1)*(ec+2*e)+v*(ec-e)*(1+2*e))/((e+2)*(ec+2*e)+2*v*(e-1)*(ec-e))
        theory=-1j*2/3*.001**3*K;err=abs(small[1][0]-theory)/abs(theory)
        record(f'quasistatic_solid_or_shell_{H}',err<1e-5,relative_a1_error=float(err))
    lossless=f.run(2,2.5,0.,kind='AUDIT_LOSSLESS',cache=False)
    record('lossless_energy_balance',abs(lossless[3]['Qabs'])<1e-10,Qabs=lossless[3]['Qabs'])
    nodes,weights=leggauss(64);a,b=saved[2][1:3]
    angular=np.array([a13.mie_s1_s2(mu,a,b) for mu in nodes]);cs=np.pi/config['wavenumber']**2*np.sum(weights*np.sum(abs(angular)**2,axis=1))
    exact=saved[2][3]['Qsca']*np.pi*config['outer_radius']**2
    record('vector_angular_power_integral',abs(cs-exact)/exact<1e-11,relative_error=float(abs(cs-exact)/exact))
    s1,s2=a13.mie_s1_s2(1.,a,b);ext=4*np.pi/config['wavenumber']**2*np.real(s1)
    record('optical_theorem',abs(ext-saved[2][3]['Qext']*np.pi)<1e-11,absolute_error=float(abs(ext-saved[2][3]['Qext']*np.pi)))
    u=np.array([0,0,1.]);pol=np.array([1.,0,0]);n=np.array([.6,0,.8]);field=a13.mie_farfield_vector(u,pol,n,a,b,config['wavenumber'],10)
    Q=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
    rotated=a13.mie_farfield_vector(Q@u,Q@pol,Q@n,a,b,config['wavenumber'],10)
    record('vector_transversality_rotation',abs(n@field)<1e-12 and np.linalg.norm(rotated-Q@field)<1e-12,transversality=float(abs(n@field)),rotation_error=float(np.linalg.norm(rotated-Q@field)))
    out=dict(status='PASS' if all(c['passed'] for c in checks.values()) else 'FAIL',checks=checks,wall_seconds=time.perf_counter()-elapsed,resources=led.summary(),scope='3D vector far-field of concentric nonmagnetic spheres only, floating-point audits, not interval certification')
    (root/'results/FORWARD_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],checks=len(checks),resources=out['resources'])));
    if out['status']!='PASS': raise RuntimeError('forward audit failed; no inference')
    return out
