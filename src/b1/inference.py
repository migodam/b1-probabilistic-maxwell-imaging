from __future__ import annotations
import time
import numpy as np
from scipy.special import logsumexp
from scipy.stats import chi2
from numpy.polynomial.legendre import leggauss


def quadrature_nodes(order,config):
    x,w=leggauss(order); z=(x+1)/2; ww=w/2; prior=config['material_prior']
    lo=np.array([prior['eps_real'][0],prior['eps_imag'][0]]); hi=np.array([prior['eps_real'][1],prior['eps_imag'][1]])
    pairs=np.asarray([(a,b) for a in z for b in z]); weights=np.asarray([a*b for a in ww for b in ww])
    return lo+pairs*(hi-lo),weights


def log_likelihood(F,y,sigma,gain_sd=0.):
    """Real Gaussian with optional marginalized g~N(1,gain_sd^2).
    C(theta)=sigma^2 I + gain_sd^2 F F^T; its determinant is retained.
    sigma is experiment-fixed, never normalized by each candidate.
    """
    F=np.asarray(F);r=np.asarray(y)-F;m=F.shape[-1]
    ff=np.sum(F*F,axis=-1); rf=np.sum(r*F,axis=-1)
    a=gain_sd**2/sigma**2; denom=1+a*ff
    q=(np.sum(r*r,axis=-1)-a*rf*rf/denom)/sigma**2
    ld=m*np.log(sigma*sigma)+np.log(denom)
    return -.5*(q+ld+m*np.log(2*np.pi)),q


def weighted_quantiles(nodes,weights,prior_node_weights,bounds):
    # Piecewise-constant reconstruction on prior-mass bins. Uniform likelihood
    # reproduces uniform-prior quantiles; quadrature CDFs are still approximate.
    unique=np.unique(nodes); mass=np.array([weights[nodes==x].sum() for x in unique]); mass/=mass.sum()
    pm=np.array([prior_node_weights[nodes==x].sum() for x in unique]); pm/=pm.sum()
    edges=bounds[0]+(bounds[1]-bounds[0])*np.r_[0.,np.cumsum(pm)]; edges[-1]=bounds[1]
    return np.interp([.05,.5,.95],np.r_[0.,np.cumsum(mass)],edges).tolist()


def infer(bank,y,sigma,prior,gain_sd=0.):
    start=time.perf_counter();H=bank['H'];F=bank['F'];weights=bank['weights']; theta=bank['theta']
    logl,mah=log_likelihood(F,y,sigma,gain_sd)
    logzh=np.asarray([logsumexp(logl[H==h]+np.log(weights[H==h])) for h in range(3)])
    total=logsumexp(np.log(prior)+logzh); joint=np.log(np.asarray(prior)[H])+logl+np.log(weights)-total
    W=np.exp(joint);p=np.exp(np.log(prior)+logzh-total)
    reps=[];summ=[]
    for h in range(3):
        ids=np.flatnonzero(H==h); cond=np.exp(logl[ids]+np.log(weights[ids])-logzh[h])
        # Representative = conditional maximum likelihood under flat continuous prior.
        # Its event mass is p[h], NOT its probability as an exact image.
        i=ids[np.argmax(logl[ids])]; reps.append(int(i));mu=cond@theta[ids];cov=((theta[ids]-mu).T*cond)@(theta[ids]-mu)
        qs=[]
        for j in range(2):
            bounds=bank.get('bounds',[(theta[:,0].min(),theta[:,0].max()),(theta[:,1].min(),theta[:,1].max())])[j]
            qs.append(weighted_quantiles(theta[ids,j],cond,weights[ids],bounds))
        summ.append(dict(H=h,mean=mu.tolist(),covariance=cov.tolist(),quantiles_05_50_95=qs,effective_quadrature_nodes=float(1/np.sum(cond**2))))
    return dict(topology_probabilities=p.tolist(),log_evidence=logzh.tolist(),log_total_evidence=float(total),
                log_bayes_factors_vs_solid=(logzh-logzh[0]).tolist(),weights=W,representative_ids=reps,conditional_material=summ,
                map_index=int(np.argmax(np.log(np.asarray(prior)[H])+logl)),
                posterior_predictive_tail_probability=float(W@chi2.sf(mah,F.shape[1])),
                posterior_expected_whitened_residual_per_component=float(W@mah/F.shape[1]),
                wall_seconds=time.perf_counter()-start)


def refinement(coarse,fine,config):
    p=np.abs(np.array(coarse['topology_probabilities'])-fine['topology_probabilities'])
    bf=np.abs(np.array(coarse['log_bayes_factors_vs_solid'])-fine['log_bayes_factors_vs_solid'])
    ok=float(p.max())<=config['probability_ref_tolerance'] and float(bf.max())<=config['log_bayes_factor_ref_tolerance']
    return dict(status='NUMERICAL_REFERENCE_PASS' if ok else 'REFERENCE_UNRESOLVED',max_probability_difference=float(p.max()),max_log_BF_difference=float(bf.max()),
                probability_tolerance=config['probability_ref_tolerance'],log_BF_tolerance=config['log_bayes_factor_ref_tolerance'],certified_bound=False)


def top_k(result,bank,K=2):
    p=np.asarray(result['topology_probabilities']);order=np.argsort(-p,kind='stable');items=[]
    for h in order[:K]:
        i=result['representative_ids'][h];items.append(dict(H=int(h),posterior_mass=float(p[h]),eps_real=float(bank['theta'][i,0]),eps_imag=float(bank['theta'][i,1]),representative_index=i))
    return dict(items=items,remaining_modeled_mass=float(p[order[K:]].sum()),omitted_topologies=[int(h) for h in order[K:]],
                outside_adopted_prior_support_mass=0.0,unmodeled_geometry_mass=None,unmodeled_geometry_status='Broader physical model not defined: mass unquantified; outside this 3-template prior is excluded by prior, not ruled out by data',
                exact_image_probability=0.,event='each discrete fixed-geometry topology with its whole continuous complex-permittivity region',probability_status='approximate numerical quadrature')


def material_image(H,theta,config,grid=151):
    z=np.linspace(-1.1,1.1,grid);X,Y=np.meshgrid(z,z);r=np.hypot(X,Y);inside=r<=config['outer_radius'];void=np.zeros_like(inside)
    for a,b in config['topologies'][int(H)]['void_intervals']:void|=(r>=a)&(r<=b)&inside
    eps=np.ones(r.shape,complex);eps[inside&~void]=complex(*theta)
    return z,eps,void


def field_summaries(result,bank,config,grid=151):
    ims=[];voids=[];cond=[]
    for h in range(3):
        z,im,void=material_image(h,result['conditional_material'][h]['mean'],config,grid);ims.append(im);voids.append(void)
        _,cov= result['conditional_material'][h]['mean'], np.array(result['conditional_material'][h]['covariance'])
        inside=(abs(im-1)>1e-12);cond.append(np.sum(np.diag(cov))*inside)
    p=np.asarray(result['topology_probabilities']);ims=np.asarray(ims);mu=np.sum(p[:,None,None]*ims,axis=0)
    variance=np.sum(p[:,None,None]*(abs(ims-mu)**2+np.asarray(cond)),axis=0)
    return dict(axis=z,posterior_mean=mu,material_complex_sd=np.sqrt(variance),void_probability=np.sum(p[:,None,None]*np.asarray(voids),axis=0),
                cavity_center_probability='fixed point at origin, not inferred; annular void has no unique point center')
