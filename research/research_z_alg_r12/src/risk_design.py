"""Predict the probability of a surviving competing minimum before drawing new data."""
from functools import lru_cache
import json
from time import perf_counter

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr

from freq_fold import candidate_basis, exact_profile
from freq_model import Fit, FrequencyTM, ri
from mirror_spectrum import Mirror, CENTER, V, STRONG, ROOT
from spectral_layouts import one

def conditional_risk(c,d,N,mu,var,side,sigma,P=2.):
    # Known calibrated gain. c and K are in an orthonormal normal coordinate system.
    A=c@c; K=2*side*c+mu*(N.T@d); B=c@K
    cmean=(K@K+var*(d@d))/A
    b=B/A
    csd=2*sigma*np.linalg.norm(c)/(np.sqrt(P)*A)
    tcov=-2*sigma**2*(c@K)/(P*A**2)
    tsd2=sigma**2*(K@K+var*(d@d))/(P*A**2)
    cond_sd=np.sqrt(max(0.,tsd2-(tcov/csd)**2))
    def density(z):
        C=cmean+csd*z
        p=C/2-3*b*b/4
        if p>=0:return 0.
        center=(b**3-b*C)/2
        width=4*(-p/3)**1.5
        Tmean=tcov*z/csd
        probability=(ndtr((center+width-Tmean)/cond_sd)-ndtr((center-width-Tmean)/cond_sd)
                     if cond_sd else float(abs(Tmean-center)<width))
        return np.exp(-z*z/2)*probability/np.sqrt(2*np.pi)
    end=min(12.,(1.5*b*b-cmean)/csd)
    return quad(density,-12.,end,epsabs=5e-9)[0] if end> -12. else 0.

def study(label,model,basis,pair,reps,seed,noise_ratio=.45):
    center,v,W,a,_=basis
    _,j,h=model.response(center,hessian=True)
    N=np.linalg.svd(ri(j)@W,full_matrices=True)[0][:,W.shape[1]:]
    c=N.T@ri(np.einsum('kij,i,j->k',h,v,v))/2
    d=ri(model.response(center,derivative=True)[3]@v)
    theta=np.linalg.norm(c)/np.linalg.norm(d)
    sigma=noise_ratio*np.linalg.norm(c)
    def risk(mult,side):
        return conditional_risk(c,d,N,0.,(theta*mult)**2/2,side,sigma)
    chosen=brentq(lambda t:max(risk(t,b) for b in (-1.,1.))-.10,.5,4.)
    dk=[('noise_free_plus_one_percent',1.01),('ten_percent_risk',chosen),('double_band',2.)]
    risks={name:[risk(t,side) for side in (-1,1)] for name,t in dk}
    record=[];tic=perf_counter();cost=model.factorizations
    for rep in range(reps):
        side=(-1,1)[rep%2]
        truth,wrong=(pair[0],pair[1]) if side==-1 else (pair[1],pair[0])
        nobs=model.response(truth)[0].size
        rng=np.random.default_rng(seed+rep)
        noise=(rng.standard_normal((3,nobs))+1j*rng.standard_normal((3,nobs)))*a*a*sigma/np.sqrt([1.,.5,.5])[:,None]
        for name,t in dk:
            off=a*theta*t*np.array([0.,-1.,1.])
            obs=[model.response(truth,b)[0]+e for b,e in zip(off,noise)]
            f=Fit(model,off,obs,np.sqrt([1.,.5,.5]))
            @lru_cache(maxsize=110)
            def gh(x):
                _,g,h,_,stable=exact_profile(f,basis,x)
                return g,h,stable
            grid=np.linspace(-4.,4.,17)
            turning=[brentq(lambda x:gh(x)[1],l,r) for l,r in zip(grid,grid[1:]) if gh(l)[1]*gh(r)[1]<0]
            points=[-4.,*turning,4.]
            roots=[brentq(lambda x:gh(x)[0],l,r) for l,r in zip(points,points[1:]) if gh(l)[0]*gh(r)[0]<0]
            rec,nlu=f.solve(wrong)
            record.append(dict(rep=rep,side=side,method=name,physical_halfband=float(max(off)),
                               known_wrong_start_RMS=float(np.linalg.norm(rec-truth)/np.sqrt(v.size)),
                               nearer_generating_branch=bool(np.linalg.norm(rec-truth)<np.linalg.norm(rec-wrong)),
                               weak_minima=sum(bool(gh(x)[1]>0) for x in roots),critical_x=roots,
                               smallest_stable_eigenvalue=min(gh(x)[2] for x in roots) if roots else None,
                               TRF_factorizations=nlu))
    return dict(label=label,geometry=getattr(model,'geometry',None),pair=np.asarray(pair).tolist(),center=center.tolist(),v=v.tolist(),
                amplitude=a,noise_ratio=noise_ratio,unit_power_sigma_real=float(a*a*sigma),risk=risks,chosen_scale=chosen,
                samples=record,wall_s=perf_counter()-tic,profile_and_fit_factorizations=model.factorizations-cost)

if __name__=='__main__':
    records=[]
    mirror=Mirror([1,2]);a=.2
    records.append(study('mirror',mirror,(CENTER,V,STRONG,a,True),[CENTER-a*V,CENTER+a*V],40,871001))
    out=ROOT/'results/risk_design.json'
    out.write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    # Next independent physical layout after the six earlier layouts, with fresh noise.
    layout=one(6)
    if layout['included']:
        model=FrequencyTM(layout['geometry']); pair=np.asarray(layout['pair']);basis=candidate_basis(model,pair)
        data=study('next_balanced_layout',model,basis,pair,40,983711)
        data['geometry_generation']=layout['fold_search']
        records.append(data)
    else:
        records.append(layout)
    out.write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
