#!/usr/bin/env python3
"""Small CPU checks for Research A. Not an imaging or speed benchmark.
Dependencies: numpy, scipy. Run: python verify_algebra.py --output results.json
Uses eight electric dipoles, six incident plane waves, and a fixed grid.
"""
from __future__ import annotations
import argparse, json, platform, time
from pathlib import Path
import numpy as np
from scipy.special import ndtr


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))


def dyadic(r, k):
    distance=np.linalg.norm(r)
    if distance == 0: raise ValueError('Self interaction excluded explicitly.')
    u=r/distance; kr=k*distance
    return np.exp(1j*kr)/(4*np.pi*distance**3)*((kr**2+1j*kr-1)*np.eye(3)+(-kr**2-3j*kr+3)*np.outer(u,u))


def alpha(chi, volume, k):
    a0=3*volume*chi/(chi+3)
    return a0/(1-1j*k**3/(6*np.pi)*a0)


def main(output):
    started=time.perf_counter(); rng=np.random.default_rng(261014)
    xyz=np.array([[x,y,z] for x in (-.22,.22) for y in (-.22,.22) for z in (-.22,.22)])
    k=1.6; volume=.012; N=len(xyz); n=3*N
    G=np.zeros((n,n),complex)
    for i in range(N):
        for j in range(N):
            if i!=j: G[3*i:3*i+3,3*j:3*j+3]=dyadic(xyz[i]-xyz[j],k)
    directions=np.vstack([np.eye(3),-np.eye(3)])
    Einc=np.empty((n,6),complex)
    for s,d in enumerate(directions):
        pol=np.eye(3)[(s%3+1)%3]
        Einc[:,s]=(np.exp(1j*k*(xyz@d))[:,None]*pol).reshape(-1)
    # Six physical far-field scalar channels (three directions, two polarizations).
    S=[]
    for d in (np.array([1.,2.,3.]),np.array([-2.,1.,2.]),np.array([2.,-3.,1.])):
        d=d/np.linalg.norm(d); p=np.cross(d,[0.,0.,1.]); p=p/np.linalg.norm(p)
        for pol in (p,np.cross(d,p)):
            S.append((np.exp(-1j*k*(xyz@d))[:,None]*pol[None,:]).reshape(-1))
    S=np.asarray(S)*k*k/(4*np.pi)
    chi0=np.full(N,1.5+.10j); a0=alpha(chi0,volume,k)
    A0=np.diag(np.repeat(a0,3)); L0=np.eye(n)-A0@G
    R0=np.linalg.solve(L0,np.eye(n)); j0=R0@A0@Einc; E0=Einc+G@j0
    chi1=chi0.copy(); chi1[0]=0.; chi1[5]=3.2+.22j
    a1=alpha(chi1,volume,k); DA=np.diag(np.repeat(a1-a0,3)); A1=A0+DA
    L1=np.eye(n)-A1@G; j1=np.linalg.solve(L1,A1@Einc); delta=j1-j0
    X=R0@DA@E0; T=R0@DA@G
    relative=np.linalg.solve(np.eye(n)-T,X)
    idx=np.r_[np.arange(0,3),np.arange(15,18)]
    U=np.eye(n)[:,idx]; D=DA[np.ix_(idx,idx)]
    Q=R0@U; core=np.eye(len(idx))-D@(U.T@G@Q)
    localized=Q@np.linalg.solve(core,D@U.T@E0)
    # Exact shared-material least-squares consistency, over all components/sources.
    E1=Einc+G@j1
    def consistency(j,E):
        J=j.reshape(N,3,6); EF=E.reshape(N,3,6)
        den=np.sum(np.abs(EF)**2,axis=(1,2))
        aa=np.sum(EF.conj()*J,axis=(1,2))/den
        return float(np.linalg.norm(J-aa[:,None,None]*EF)/np.linalg.norm(J)),aa
    exact_defect, a_recovered=consistency(j1,E1)
    Pdark=np.eye(n)-np.linalg.pinv(S)@S
    v=Pdark@(rng.normal(size=j0.shape)+1j*rng.normal(size=j0.shape))
    v*=.2*np.linalg.norm(j0)/np.linalg.norm(v)
    dark_j=j0+v; dark_E=Einc+G@dark_j
    dark_defect,_=consistency(dark_j,dark_E)
    # At the tangent level require L0*v_s = delta_a*E0_s, same scalar at each cell.
    tangent_defect,_=consistency(L0@v,E0)
    # Exact dual residual identity for an arbitrary non-Galerkin approximation.
    jhat=j1+.025*(rng.normal(size=j1.shape)+1j*rng.normal(size=j1.shape))
    Z=np.linalg.solve(L1.conj().T,S.conj().T)
    Zhat=Z+.04*(rng.normal(size=Z.shape)+1j*rng.normal(size=Z.shape))
    rp=A1@Einc-L1@jhat; rd=S.conj().T-L1.conj().T@Zhat
    ycorr=S@jhat+Zhat.conj().T@rp
    dual_remainder=rd.conj().T@np.linalg.solve(L1,rp)
    # Neumann truncation bound and quadratic/cubic response slopes in alpha coordinates.
    h=np.repeat(rng.normal(size=N)+.2j*rng.normal(size=N),3)
    H=np.diag(h)*.006
    amps=np.array([.0125,.025,.05,.1,.2]); first=[]; second=[]; etas=[]
    for t in amps:
        xt=R0@(t*H)@E0; tt=R0@(t*H)@G
        exact=np.linalg.solve(np.eye(n)-tt,xt)
        first.append(np.linalg.norm(S@(exact-xt)))
        second.append(np.linalg.norm(S@(exact-xt-tt@xt)))
        etas.append(float(np.linalg.norm(tt,2)))
    eta=float(np.linalg.norm(T,2))
    bound=float(np.linalg.norm(S,2)*eta**2/(1-eta)*np.linalg.norm(X)) if eta<1 else None
    remainder=float(np.linalg.norm(S@(delta-X-T@X)))
    # Pure algebra counterexample: a receiver-dark change that no common material realizes.
    # Two cells; receiver sees only cell 1. Cell 2 has fields 1 and 2 for the two sources.
    e=np.array([1.,2.]); requested=np.array([1.,1.]); fitted=(e@requested)/(e@e)
    inconsistent_res=float(np.linalg.norm(requested-fitted*e))
    # A single global material parameter changes a diagonal matrix at every cell.
    aglobal=alpha(chi0+.1,volume,k); rank_global=int(np.linalg.matrix_rank(np.diag(np.repeat(aglobal-a0,3))))
    results={
      'scope':'8-dipole CPU algebra verification; not reconstruction, topology coverage, runtime comparison, or a continuum certificate',
      'seed':261014,'n_cells':N,'n_complex_current':n,'n_sources':6,'n_receiver_channels_per_source':len(S),
      'backend':'self-contained electric dyadic dipole matrix; CM plus radiation reaction; diagonal Green excluded',
      'physical_convention':'exp(-i omega t), scalar isotropic passive chi in real free-space background',
      'finite_edit':'one cell removed to vacuum, one contrast changed; fixed Green matrix',
      'condition_L0':float(np.linalg.cond(L0)),'condition_L1':float(np.linalg.cond(L1)),
      'relative_feedback_norm':eta,
      'finite_relative_identity_relative_error':rel(relative,delta),
      'localized_schur_relative_error':rel(localized,delta),
      'localized_output_relative_error':rel(S@localized,S@delta),
      'exact_current_material_consistency_relative_defect':exact_defect,
      'recovered_alpha_relative_error':rel(a_recovered,a1),
      'dark_current_receiver_relative_error':float(np.linalg.norm(S@v)/(np.linalg.norm(S)*np.linalg.norm(v))),
      'dark_current_material_consistency_relative_defect':dark_defect,
      'dark_tangent_common_material_relative_defect':tangent_defect,
      'dual_corrected_residual_identity_relative_error':rel(S@j1-ycorr,dual_remainder),
      'neumann_remainder_observed':remainder,'neumann_remainder_upper_bound_if_eta_lt_1':bound,
      'linear_response_error_slope':float(np.polyfit(np.log(amps),np.log(first),1)[0]),
      'quadratic_response_error_slope':float(np.polyfit(np.log(amps),np.log(second),1)[0]),
      'two_illumination_inconsistent_current_best_scalar_residual':inconsistent_res,
      'one_parameter_global_matrix_update_rank':rank_global,
      'two_fixed_Gaussian_hypotheses_equal_prior_errors':{str(d):float(ndtr(-d/2)) for d in (1.529652,1.295221,.900170,2.467944,8.84)},
      'sphere_cavity_r_0_25_relative_volume':.25**3,
      'runtime_seconds':time.perf_counter()-started,
      'versions':{'python':platform.python_version(),'numpy':np.__version__}
    }
    assert results['finite_relative_identity_relative_error']<1e-11
    assert results['localized_schur_relative_error']<1e-11
    assert exact_defect<1e-12 and dark_defect>1e-3
    assert results['dark_current_receiver_relative_error']<1e-12
    assert results['dual_corrected_residual_identity_relative_error']<1e-10
    assert bound is not None and remainder<=bound*(1+1e-8)
    assert 1.95<results['linear_response_error_slope']<2.05
    assert 2.90<results['quadratic_response_error_slope']<3.10
    results['all_assertions_passed']=True
    Path(output).write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',default='results.json')
    args=parser.parse_args(); main(args.output)
