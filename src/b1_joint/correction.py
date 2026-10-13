"""Stratified continuous importance residuals correct the full-Mie target.

Z_full(A) = Z_proxy(A) + E_qA[(L_full-L_proxy) p0/qA].
qA is a 95% computable adaptive-cell proposal + 5% conditional proper prior.
All likelihood evaluations use the same observation and bounded shared gain.
Numerical cubature error and Monte Carlo residual error remain separate.
"""
from __future__ import annotations
import numpy as np
from scipy.special import logsumexp
from .material import BOUNDS,WIDTH,SCALES,geometry,event_intervals,prior_event_masses,prior_logpdf
from .likelihood import log_likelihood

def coefficient_box(box):
    result=[]
    for j,b in enumerate(box):
        if b is None:result.append((-1.,1.))
        elif j<2:
            lo,hi=b;result.append((-(1-(hi/SCALES[j])**2)**-2,-(1-(lo/SCALES[j])**2)**-2))
        else:result.append(tuple(b))
    return np.array(result)

class EventProposal:
    def __init__(self,leaves,event):
        selected=[a for a in leaves if a['event_id']==event];self.boxes=np.array([coefficient_box(a['box']) for a in selected]);self.volumes=np.prod(np.diff(self.boxes,axis=2).squeeze(-1)/WIDTH,axis=1)
        self.prob=np.array([a['scaled_mass'] for a in selected]);self.prob/=self.prob.sum();self.event=event
        c,s=divmod(event,3);intervals=event_intervals();self.priorbox=np.array([intervals[0][c] or (-1.,1.),intervals[1][s] or (-1.,1.),BOUNDS[2],BOUNDS[3]])
        self.prior_mass=prior_event_masses()[event]
    def sample(self,rng,n):
        which=rng.choice(len(self.boxes),n,p=self.prob);box=self.boxes[which].copy()
        use_prior=rng.random(n)<.05;box[use_prior]=self.priorbox
        return rng.uniform(box[:,:,0],box[:,:,1])
    def logpdf(self,x):
        x=np.atleast_2d(x);a=(x[:,None,:]>=self.boxes[None,:,:,0]).all(-1)&(x[:,None,:]<=self.boxes[None,:,:,1]).all(-1)
        # Random continuous draws have zero probability of leaf boundaries.
        if not np.all(a.sum(1)==1):raise ValueError('proposal partition did not cover exactly once')
        factor=.95*(a@(self.prob/self.volumes))+.05/self.prior_mass
        return -np.log(WIDTH).sum()+np.log(factor)

def corrected_reference(result,adaptive,response,forward,n_per_event,rng,boot=2000):
    shift=float(np.max(result['loglik']));base_total=np.exp(result['log_evidence']-shift);base=result['event_mass']*base_total
    deltas=[];is_weights=[];samples=[];fullfields=[];proxyfields=[];logqs=[];features=[]
    for j in range(12):
        proposal=EventProposal(adaptive['leaves'],j);x=proposal.sample(rng,n_per_event);q=proposal.logpdf(x)
        full=np.array([forward.predict(xx,1.6,kind='EVENT_FULL_LIKELIHOOD_RESIDUAL_CORRECTION')[0] for xx in x]);proxy=response.predict(x)
        lf=log_likelihood(full,result['y'],result['sigma']);lp=log_likelihood(proxy,result['y'],result['sigma']);factor=np.exp(prior_logpdf(x)-q)
        residual=np.exp(lp-shift)*np.expm1(lf-lp)*factor
        samples.append(x);fullfields.append(full);proxyfields.append(proxy);logqs.append(q);deltas.append(residual)
        is_weights.append(np.exp(lf-shift)*factor)
        features.append(np.column_stack((geometry(x),x[:,2:])))
    deltas=np.array(deltas);is_weights=np.array(is_weights);x=np.array(samples);v=np.array(features)
    z=base+deltas.mean(axis=1)
    if np.any(z<=0):raise FloatingPointError('negative corrected event integral; no clipping allowed')
    total=z.sum();p=z/total
    boot_z=np.stack([base[j]+rng.choice(deltas[j],(boot,n_per_event),replace=True).mean(1) for j in range(12)],axis=1)
    if np.any(boot_z<=0):raise FloatingPointError('correction interval unresolved: negative bootstrap evidence')
    bp=boot_z/boot_z.sum(1,keepdims=True);interval=np.quantile(bp,[.025,.975],axis=0).T
    conditional=[]
    for j in range(12):
        mean0=np.asarray(result['conditional'][j]['mean']);cov0=np.asarray(result['conditional'][j]['covariance'])
        first=base[j]*mean0+np.mean(deltas[j,:,None]*v[j],axis=0)
        second=base[j]*(cov0+np.outer(mean0,mean0))+np.mean(deltas[j,:,None,None]*(v[j,:, :,None]*v[j,:,None,:]),axis=0)
        mean=first/z[j];cov=second/z[j]-np.outer(mean,mean)
        conditional.append(dict(event_id=j,mean=mean,SD=np.sqrt(np.maximum(np.diag(cov),0)),covariance=cov))
    # Weighted conditional full-Mie samples; event masses use the lower-variance
    # residual-control-variate estimator, disclosed as different approximations.
    conditional_weights=is_weights/is_weights.sum(1,keepdims=True)
    weights=(p[:,None]*conditional_weights).ravel()
    summary=dict(event_mass=p,event_mass_interval95=interval,event_mass_MC_se=bp.std(axis=0,ddof=1),proxy_event_mass=result['event_mass'],max_event_change=float(abs(p-result['event_mass']).max()),log_evidence=float(np.log(total)+shift),conditional=conditional,n_per_event=n_per_event,physical_calls=12*n_per_event,bootstrap_replicates=boot,proposal='0.95 adaptive-cell density +0.05 prior conditional on each event; exact density q in latent dx',estimator='base continuous proxy integral + independent full-minus-proxy importance residual',same_prior_likelihood_and_noise=True,numerical_status='EMPIRICAL_FULL_MIE_RESIDUAL_REFERENCE; no uniform/continuum certificate',numerical_base_error='tensor17 vs13 + independent physical-coordinate cubature separately reported',mass_interval_scope='Monte Carlo residual only; does not include quadrature or global untested solver error')
    arrays=dict(latent=x.reshape(-1,4),weights=weights,event_ids=np.repeat(np.arange(12),n_per_event),F_full=np.array(fullfields).reshape(-1,32),F_proxy=np.array(proxyfields).reshape(-1,32),log_q_conditional=np.array(logqs).ravel(),residual_weights=deltas,full_is_weights=is_weights,base_event_integrals=base,likelihood_shift=shift,y=result['y'],sigma=result['sigma'],event_mass=p,event_mass_interval95=interval)
    return summary,arrays
