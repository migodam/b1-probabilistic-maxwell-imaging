"""Conditional predictive decision design. This module never imports a solver.

Particles approximate the stated continuous target; they are not a new prior.
Candidate-only distances are set heuristics and have no posterior meaning.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import logsumexp
from scipy.stats import truncnorm
from src.b1_joint.material import geometry, event_ids, RMAX
from src.b1_joint.likelihood import log_likelihood, conditional_gain, log_cdf_difference

@dataclass
class Action:
    id: str
    fields: np.ndarray
    sigma: np.ndarray
    blocks: float=1.
    kind: str='frequency'

def field_l1(x,y):
    """Exact volume-weighted componentwise epsilon L1, s_ref=1, R=1."""
    x=np.atleast_2d(x);y=np.atleast_2d(y);gx=geometry(x);gy=geometry(y)
    def ring(h):return np.where(h>0,(.65+h)**3-(.65-h)**3,0.)
    vx=gx[:,0]**3+ring(gx[:,1]);vy=gy[:,0]**3+ring(gy[:,1])
    overlap=np.minimum(gx[:,0,None],gy[None,:,0])**3+ring(np.minimum(gx[:,1,None],gy[None,:,1]))
    both=1-vx[:,None]-vy[None,:]+overlap
    diff=np.abs(x[:,None,2]-y[None,:,2])+np.abs(x[:,None,3]-y[None,:,3])
    ex=np.abs(x[:,2]-1)+np.abs(x[:,3]);ey=np.abs(y[:,2]-1)+np.abs(y[:,3])
    return np.maximum(both*diff+(vy[None,:]-overlap)*ex[:,None]+(vx[:,None]-overlap)*ey[None,:],0.)

def normalized(logw):
    return np.exp(logw-logsumexp(logw,axis=-1,keepdims=True))

def resample_indices(weights,n,seed):
    rng=np.random.default_rng(seed);u=(np.arange(n)+rng.random(n))/n
    c=np.cumsum(weights);c[-1]=1
    return np.searchsorted(c,u)

def posterior_weights(prior, fields, y, sigma):
    return normalized(np.log(prior)+log_likelihood(fields,y,sigma))

def conditional_log_increment(f0,y0,s0,fa,ya,sa,model='shared'):
    """All old/new blocks share a single persistent bounded gain per world.

    Returns draws-by-material conditional log likelihood up to a draw-only
    constant, sufficient for normalized posterior and EIG/risk. The determinant
    and truncation interval terms are retained.
    """
    f0=np.asarray(f0);fa=np.asarray(fa);ya=np.atleast_2d(ya)
    if model=='independent':
        return np.array([log_likelihood(fa,z,sa) for z in ya])
    if model=='fixed': return -.5*np.sum(((fa[None]-ya[:,None])/sa)**2,axis=-1)
    # Integrate the conditional truncated normal for gain analytically.
    m,v=conditional_gain(f0,y0,s0); A=1/v;B=m/v
    af=fa/sa; z=ya/sa
    addA=np.sum(af*af,axis=1); newA=A+addA;newB=B[None]+z@af.T
    post_m=newB/newA;post_s=1/np.sqrt(newA)
    logZ=log_cdf_difference((.7-m)/np.sqrt(v),(1.3-m)/np.sqrt(v))
    logZnew=log_cdf_difference((.7-post_m)/post_s,(1.3-post_m)/post_s)
    return .5*(newB**2/newA-B[None]**2/A)-.5*np.log(newA/A)[None]+logZnew-logZ[None]

def predictive_draws(w,f0,y0,s0,action,uniforms,normals,gain_uniforms):
    c=np.cumsum(w);c[-1]=1.;ids=np.searchsorted(c,uniforms)
    m,v=conditional_gain(f0[ids],y0,s0);sd=np.sqrt(v)
    g=truncnorm.ppf(gain_uniforms,(.7-m)/sd,(1.3-m)/sd,loc=m,scale=sd)
    return g[:,None]*action.fields[ids]+normals[:,:action.fields.shape[1]]*action.sigma

def expected_utility(x,w,f0,y0,s0,action,decisions,seed=26101530,outer=256,model='shared'):
    rng=np.random.default_rng(seed);u=rng.random(outer);ug=rng.random(outer);z=rng.normal(size=(outer,64))
    obs=predictive_draws(w,f0,y0,s0,action,u,z,ug)
    inc=conditional_log_increment(f0,y0,s0,action.fields,obs,action.sigma,model)
    post=normalized(np.log(np.maximum(w,1e-300))[None]+inc)
    loss=field_l1(decisions,x);risk0=float(np.min(loss@w));risks=np.min(post@loss.T,axis=1)
    events=event_ids(x);onehot=np.eye(12)[events];p=w@onehot;ep=post@onehot
    H=lambda a: -np.sum(np.where(a>0,a*np.log(np.maximum(a,1e-300)),0),axis=-1)
    gains=risk0-risks
    return dict(current_risk=risk0,expected_future_risk=float(risks.mean()),risk_reduction=float(gains.mean()),
        mc_se=float(gains.std(ddof=1)/np.sqrt(outer)),event_EIG=float(H(p)-H(ep).mean()),
        expected_event_error=float((1-ep.max(1)).mean()),event_quantiles=np.quantile(ep,[.05,.5,.95],axis=0),
        posterior_mass_martingale_max=float(abs(ep.mean(0)-p).max()),outer=outer,finite_bank=True,
        prediction_model=model,risk_draws=risks)

def joint_gain_pair(f0a,f0b,faa,fab,s0,sa,y0):
    """Bounded per-world gains; shared within a world, independent across worlds.

    Compare each candidate's gain to its old-data fitted gain with old-data
    curvature, then obtain the new pair's minimal distance. This conditional
    regularization prevents declaring full prior gain freedom after conditioning.
    The finite candidate-pair minimum is a screening score, never a global bound.
    """
    from scipy.optimize import lsq_linear
    ma,va=conditional_gain(np.atleast_2d(f0a),y0,s0);mb,vb=conditional_gain(np.atleast_2d(f0b),y0,s0)
    M=np.column_stack((faa/sa,-fab/sa));P=np.diag(1/np.sqrt([va[0],vb[0]]));q=P@np.clip([ma[0],mb[0]],.7,1.3)
    fit=lsq_linear(np.vstack((M,P)),np.r_[np.zeros(M.shape[0]),q],bounds=(.7,1.3),tol=1e-10)
    return float(np.sum((M@fit.x)**2)),fit.x

def set_scores(x,f0,y0,s0,action,max_pairs=128):
    """No candidate probabilities: material-separated composite pair screening."""
    md=field_l1(x,x);event=event_ids(x);pairs=np.array(np.where(np.triu((event[:,None]!=event[None])&(md>=.005),1))).T
    if len(pairs)>max_pairs:
        # Include hardest current-data pairs, with deterministic tie handling.
        m,v=conditional_gain(f0,y0,s0);gg=np.clip(m,.7,1.3)
        scores=np.sum(((gg[pairs[:,0],None]*f0[pairs[:,0]]-gg[pairs[:,1],None]*f0[pairs[:,1]])/s0)**2,axis=1)
        pairs=pairs[np.argsort(scores,kind='stable')[:max_pairs]]
    if not len(pairs):return dict(raw=0.,robust=0.,conflict=0.,pair_count=0)
    raw=[];robust=[];conflict=[]
    for a,b in pairs:
        raw.append(np.sum(((action.fields[a]-action.fields[b])/action.sigma)**2))
        score,g=joint_gain_pair(f0[a],f0[b],action.fields[a],action.fields[b],s0,action.sigma,y0);robust.append(score)
        # Existing two-world secant compensation fits: one scale maps b into a.
        # Stack first, optimize once, then compare to two independent scales.
        f,g0=f0[b]/s0,f0[a]/s0;fa,ga=action.fields[b]/action.sigma,action.fields[a]/action.sigma
        fit=lambda u,v: np.clip(np.dot(u,v)/max(np.dot(u,u),1e-300),.7/1.3,1.3/.7)
        c0=fit(f,g0);ca=fit(fa,ga);cj=fit(np.r_[f,fa],np.r_[g0,ga])
        sep=np.sum((c0*f-g0)**2)+np.sum((ca*fa-ga)**2)
        joint=np.sum((cj*np.r_[f,fa]-np.r_[g0,ga])**2)
        conflict.append(max(float(joint-sep),0.))
    # Material relevance weighting does not turn pair scores into a posterior.
    return dict(raw=float(np.mean(raw)),robust=float(np.min(robust)),
        conflict=float(np.mean(np.array(conflict)*md[pairs[:,0],pairs[:,1]])),pair_count=len(pairs),
        min_pair_robust=float(np.min(robust)),median_pair_robust=float(np.median(robust)))
