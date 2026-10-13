"""Generic tempered SMC with prior-preserving moves and full-Mie DA correction.

The proxy only screens proposals; accepted moves retain the specified numerical
full-Mie marginal target. Failure to reach beta=1 is reported, never called a
posterior. Finite-particle coverage is not guaranteed by ESS.
"""
from __future__ import annotations
import numpy as np
from scipy.special import ndtr,ndtri,logsumexp
from .material import BOUNDS,WIDTH,prior_sample,event_ids
from .likelihood import log_likelihood

def systematic(w,rng):return np.searchsorted(np.cumsum(w),(np.arange(len(w))+rng.random())/len(w))

def run_smc(response,forward,y,sigma,seed,n=64,full_cap=340):
    rng=np.random.default_rng(seed);ledger=forward.backend.ledger;before=ledger.count();x=prior_sample(rng,n)
    full=np.array([forward.predict(xx,1.6,kind='SMC_INITIAL_FULL')[0] for xx in x]);ll=log_likelihood(full,y,sigma);la=log_likelihood(response.predict(x),y,sigma)
    logw=np.full(n,-np.log(n));beta=0.;logz=0.;trace=[];status='FULL_TARGET_REACHED';stage2_total=0;stage2_reject=0
    for stage in range(64):
        def cess(delta):return n*np.exp(2*logsumexp(logw+delta*ll)-logsumexp(logw+2*delta*ll))
        hi=1-beta;lo=0.
        if cess(hi)<.8*n:
            for _ in range(40):
                mid=(lo+hi)/2
                if cess(mid)<.8*n:hi=mid
                else:lo=mid
            delta=(lo+hi)/2
        else:delta=hi
        nextbeta=min(1.,beta+delta);increment=logsumexp(logw+delta*ll);logz+=increment;logw=logw+delta*ll-increment;w=np.exp(logw);ess=1/np.sum(w*w);resampled=False
        if ess<.65*n:
            ids=systematic(w,rng);x=x[ids].copy();ll=ll[ids].copy();la=la[ids].copy();logw.fill(-np.log(n));resampled=True
        accepted=0;screened=0;full_used=0
        for i in range(n):
            if ledger.count()-before>=full_cap:status='FULL_QUERY_SUBBUDGET_REACHED';break
            t=ndtri((x[i]-BOUNDS[:,0])/WIDTH)
            proposal=(np.sqrt(1-.35**2)*t+.35*rng.normal(size=4)) if rng.random()>.15 else rng.normal(size=4)
            xp=BOUNDS[:,0]+WIDTH*ndtr(proposal)
            if np.any(xp<=BOUNDS[:,0]) or np.any(xp>=BOUNDS[:,1]):continue
            lap=log_likelihood(response.predict(xp),y,sigma)[0]
            if np.log(rng.random())>min(0.,nextbeta*(lap-la[i])):continue
            screened+=1;fp=forward.predict(xp,1.6,kind='SMC_DELAYED_ACCEPTANCE_FULL')[0];lp=float(log_likelihood(fp,y,sigma));full_used+=1;stage2_total+=1
            if np.log(rng.random())<=min(0.,nextbeta*(lp-ll[i]+la[i]-lap)):
                x[i]=xp;ll[i]=lp;la[i]=lap;accepted+=1
            else:stage2_reject+=1
        trace.append(dict(beta=nextbeta,conditional_ESS_target=.8*n,ESS_before_resampling=ess,resampled=resampled,screened=screened,accepted=accepted,new_full=full_used))
        beta=nextbeta
        if beta>=1-1e-12:break
        if status!='FULL_TARGET_REACHED':break
    if beta<1-1e-12:status='PARTIAL_TEMPERATURE_NOT_POSTERIOR'
    w=np.exp(logw);p=np.bincount(event_ids(x),weights=w,minlength=12)
    return dict(status=status,beta=beta,seed=seed,particles=n,event_mass=p,ESS=float(1/np.sum(w*w)),log_evidence=float(logz),trace=trace,actual_full_evaluations=ledger.count()-before,stage2_rejections=stage2_reject,stage2_trials=stage2_total,target='full-vector-Mie-matching-n14:C:bounded-shared-gain',finite_particle_guarantee=False),dict(latent=x,weights=w,log_likelihood=ll)
