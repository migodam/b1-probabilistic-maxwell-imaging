"""Finite quadrature probability algebra only. No PDE, Book, data or truth I/O.

Default complex observation shape(5,4,12), sigma=.01 for EACH real quadrature.
One shared gamma=a exp(i phi), a~U[.99,1.01],phi~U[-.01,.01]. Nuisance
integration is normalized Gauss-Legendre quadrature, never profiling.
Synthetic checks cannot establish the P1/Maxwell qualification gate.
"""
from __future__ import annotations
import numpy as np

OBSERVATION_SHAPE=(5,4,12)
AMPLITUDE_BOUNDS=(.99,1.01)
PHASE_BOUNDS=(-.01,.01)


def _logsumexp(values,axis=-1):
    a=np.asarray(values,dtype=np.longdouble)
    top=np.max(a,axis=axis,keepdims=True)
    safe=np.where(np.isfinite(top),top,0.)
    total=np.sum(np.exp(a-safe),axis=axis,keepdims=True)
    with np.errstate(divide='ignore'):
        out=safe+np.log(total)
    return np.squeeze(out,axis=axis)


def _shape(shape):
    s=tuple(shape)
    if len(s)!=3 or any(isinstance(z,bool) or int(z)!=z or z<=0 for z in s):
        raise ValueError('observation_shape must be three positive integer axes(frequency,source,receiver)')
    return tuple(int(z) for z in s)


def _observations(f,y,shape):
    shape=_shape(shape)
    f,y=np.asarray(f),np.asarray(y)
    if f.shape[-3:]!=shape or y.shape[-3:]!=shape:
        raise ValueError(f'f and y must end with observation shape{shape}')
    if not np.isfinite(f).all() or not np.isfinite(y).all():
        raise ValueError('f,y must be finite')
    return np.broadcast_arrays(f,y)


def nuisance_quadrature(n=11):
    """Normalized uniform-prior quadrature; weights sum1, not rectangle area."""
    if isinstance(n,bool) or int(n)!=n or n<1:
        raise ValueError('n must be a positive integer; frozen comparisons use11/21')
    z,w=np.polynomial.legendre.leggauss(int(n))
    a=(AMPLITUDE_BOUNDS[0]+AMPLITUDE_BOUNDS[1])/2+(AMPLITUDE_BOUNDS[1]-AMPLITUDE_BOUNDS[0])*z/2
    p=(PHASE_BOUNDS[0]+PHASE_BOUNDS[1])/2+(PHASE_BOUNDS[1]-PHASE_BOUNDS[0])*z/2
    aa,pp=np.meshgrid(a,p,indexing='ij')
    weights=(w[:,None]*w[None,:]/4).ravel()
    return dict(amplitude=aa.ravel().astype(np.longdouble),phase=pp.ravel().astype(np.longdouble),
        weights=weights.astype(np.longdouble),n=int(n),node_count=int(n)**2)


def _blocks(blocks,ns):
    if blocks is None:return [tuple(range(ns))]
    if blocks is True:return [(j,) for j in range(ns)]
    result=[]
    for block in blocks:
        b=tuple(block)
        if not b or any(isinstance(j,bool) or int(j)!=j or not 0<=j<ns for j in b):
            raise ValueError('source blocks must be nonempty lists of valid source indices')
        result.append(tuple(int(j) for j in b))
    if sorted(j for b in result for j in b)!=list(range(ns)):
        raise ValueError('source blocks must partition all sources exactly once')
    return result


def sufficient_statistics(f,y,*,observation_shape=OBSERVATION_SHAPE,source_blocks=None):
    """F=sum|f|²,C=sum f conj(y),Y=sum|y|²; block sums before nuisance integral."""
    shape=_shape(observation_shape);f,y=_observations(f,y,shape)
    blocks=_blocks(source_blocks,shape[1]);rows=[]
    for ids in blocks:
        fb,yb=f[..., :,list(ids),:],y[..., :,list(ids),:]
        fr,fi=np.asarray(fb.real,dtype=np.longdouble),np.asarray(fb.imag,dtype=np.longdouble)
        yr,yi=np.asarray(yb.real,dtype=np.longdouble),np.asarray(yb.imag,dtype=np.longdouble)
        axes=(-3,-2,-1)
        rows.append(dict(source_indices=ids,complex_count=shape[0]*len(ids)*shape[2],
            F=np.sum(fr*fr+fi*fi,axis=axes),Y=np.sum(yr*yr+yi*yi,axis=axes),
            C_real=np.sum(fr*yr+fi*yi,axis=axes),C_imag=np.sum(fi*yr-fr*yi,axis=axes)))
    totals={key:sum(row[key] for row in rows) for key in ('F','Y','C_real','C_imag')}
    return dict(**totals,complex_count=int(np.prod(shape)),source_blocks=rows,
        shared_nuisance_rule='sum block statistics, form ONE likelihood, integrate gamma ONCE')


def log_marginal_likelihood(f,y,*,sigma=.01,n=11,observation_shape=OBSERVATION_SHAPE,
                            source_blocks=None,return_node_weights=False):
    """Absolute density integral p(y|f), broadcast over f/y batch axes.

    Phi=(a² F-2a Re(exp(i phi)C)+Y)/(2 sigma²). Normalization is
    -(2M)log(sigma sqrt(2pi)); no channel-count normalization. Source blocks
    are only an additive evaluation device, not independent nuisance fits.
    Return absolute logZ and posterior nuisance moments under quadrature.
    All arrays are finite quadrature approximations, not continuous proofs.
    """
    if not np.isfinite(sigma) or sigma<=0:raise ValueError('sigma must be finite positive')
    stat=sufficient_statistics(f,y,observation_shape=observation_shape,source_blocks=source_blocks)
    grid=nuisance_quadrature(n);a,p=grid['amplitude'],grid['phase']
    F,Y,cr,ci=(stat[k][...,None] for k in ('F','Y','C_real','C_imag'))
    cross=np.cos(p)*cr-np.sin(p)*ci
    residual_sq=a*a*F-2*a*cross+Y
    roundoff=64*np.finfo(np.longdouble).eps*(a*a*F+2*np.abs(a*cross)+Y)
    if np.any(residual_sq < -roundoff):
        raise FloatingPointError('Sufficient-statistic residual became materially negative')
    clipped=int(np.count_nonzero(residual_sq<0))
    residual_sq=np.maximum(residual_sq,0.)
    constant=-2*stat['complex_count']*np.log(np.longdouble(sigma)*np.sqrt(np.longdouble(2)*np.pi))
    minimum_square=np.min(residual_sq,axis=-1,keepdims=True)
    relative_nodes=-(residual_sq-minimum_square)/(2*np.longdouble(sigma)**2)+np.log(grid['weights'])
    relative_logz=_logsumexp(relative_nodes)
    logz=constant-np.squeeze(minimum_square,axis=-1)/(2*np.longdouble(sigma)**2)+relative_logz
    pw=np.exp(relative_nodes-relative_logz[...,None])
    pw/=np.sum(pw,axis=-1,keepdims=True)
    am=np.sum(pw*a,axis=-1);pm=np.sum(pw*p,axis=-1)
    gamma_real=a*np.cos(p);gamma_imag=a*np.sin(p)
    gr=np.sum(pw*gamma_real,axis=-1);gi=np.sum(pw*gamma_imag,axis=-1)
    moments=dict(amplitude_mean=am,amplitude_variance=np.sum(pw*(a-am[...,None])**2,axis=-1),
        phase_mean=pm,phase_variance=np.sum(pw*(p-pm[...,None])**2,axis=-1),
        amplitude_phase_covariance=np.sum(pw*(a-am[...,None])*(p-pm[...,None]),axis=-1),
        gamma_mean=gr+1j*gi,gamma_real_mean=gr,gamma_imag_mean=gi,
        gamma_modulus_squared_mean=np.sum(pw*a*a,axis=-1),
        gamma_real_imag_covariance=np.sum(pw*(gamma_real-gr[...,None])*(gamma_imag-gi[...,None]),axis=-1),
        gamma_real_variance=np.sum(pw*(gamma_real-gr[...,None])**2,axis=-1),
        gamma_imag_variance=np.sum(pw*(gamma_imag-gi[...,None])**2,axis=-1))
    result=dict(log_normalizing_constant=constant,logZ=logz,posterior_nuisance_moments=moments,
        quadrature_n=int(n),quadrature_node_count=int(n)**2,prior_quadrature_weight_sum=np.sum(grid['weights']),
        complex_observation_count=stat['complex_count'],real_observation_count=2*stat['complex_count'],
        sigma_per_real_quadrature=float(sigma),sufficient_statistics=stat,
        residual_roundoff_clipped_node_count=clipped,
        maximum_residual_square_roundoff_scale=np.max(roundoff),
        roundoff_diagnostic_status='heuristic64eps scale; not a certified numerical error bound',
        accumulation_dtype=str(np.dtype(np.longdouble)),
        accumulation_epsilon=np.finfo(np.longdouble).eps,
        approximation='finite Gauss-Legendre integration of one shared nuisance; convergence not certified',
        P1_gate='NOT_EVALUATED_ALGEBRA_ONLY')
    if return_node_weights:result.update(nuisance_nodes=grid,posterior_node_weights=pw)
    return result


def finite_grid_posterior(log_likelihood,prior_density,cell_weights):
    """Proper density*volume posterior on a FINITE grid; never node counting.

    Raw quadrature evidence uses supplied prior density. The integral of that
    prior quadrature may differ from1, so expose it and also return evidence
    for the separately normalized finite prior. No silent evidence correction.
    """
    ll,p,v=(np.asarray(z,dtype=np.longdouble) for z in (log_likelihood,prior_density,cell_weights))
    if ll.ndim!=1 or p.shape!=ll.shape or v.shape!=ll.shape or not len(ll):
        raise ValueError('Likelihood,prior_density,cell_weights must be nonempty matching1D vectors')
    if np.any(np.isnan(ll)) or np.any(ll==np.inf) or not np.isfinite(p).all() or not np.isfinite(v).all() or np.any(p<0) or np.any(v<0):
        raise ValueError('Finite nonnegative density/weights; log likelihood may include-inf only')
    if not np.any((p>0)&(v>0)):raise ValueError('Prior quadrature must have positive mass')
    with np.errstate(divide='ignore'):logmass=np.log(p)+np.log(v)
    logprior=_logsumexp(logmass)
    anchor=np.max(ll)
    if not np.isfinite(anchor):raise ValueError('No positive posterior mass on finite nodes')
    relative=ll-anchor+logmass
    logrelative=_logsumexp(relative);loge=anchor+logrelative
    if not np.isfinite(loge):raise ValueError('No positive posterior mass on finite nodes')
    posterior=np.exp(relative-logrelative);posterior/=np.sum(posterior)
    return dict(posterior_weights=posterior,prior_quadrature_mass=np.exp(logprior),
        log_prior_quadrature_mass=logprior,
        log_evidence_raw_quadrature=loge,
        log_evidence_normalized_finite_prior=loge-logprior,
        interpretation='finite quadrature approximation; no continuous posterior/CDF certification',
        P1_gate='NOT_EVALUATED_ALGEBRA_ONLY')


def _weights(weights,n):
    w=np.asarray(weights,dtype=np.longdouble)
    if w.shape!=(n,) or not np.isfinite(w).all() or np.any(w<0) or abs(np.sum(w)-1)>1e-10:
        raise ValueError('weights must be finite nonnegative normalized1D probabilities')
    return w


def finite_grid_summary(weights,*,branch_ids=None,event_masks=None,target_positions=None):
    w=np.asarray(weights);w=_weights(w,len(w));out={'interpretation':'finite weighted quadrature only'}
    if branch_ids is not None:
        b=np.asarray(branch_ids)
        if b.shape!=w.shape:raise ValueError('branch_ids must match weights')
        out['branch_weights']=[dict(branch=label.item() if isinstance(label,np.generic) else label,
            weight=np.sum(w[b==label])) for label in np.unique(b)]
    if event_masks is not None:
        out['event_weights']={}
        for name,mask in event_masks.items():
            m=np.asarray(mask)
            if m.dtype!=bool or m.shape!=w.shape:raise ValueError('Each event mask must be boolean and match weights')
            out['event_weights'][name]=np.sum(w[m])
    if target_positions is not None:
        x=np.asarray(target_positions,dtype=np.longdouble)
        if x.ndim==1:x=x[:,None]
        if x.ndim!=2 or len(x)!=len(w) or not np.isfinite(x).all():raise ValueError('Positions must have finite shape(N,d)')
        mean=np.sum(w[:,None]*x,axis=0);center=x-mean
        cov=np.einsum('i,ij,ik->jk',w,center,center)
        out.update(target_position_mean=mean,target_position_variance=np.diag(cov),target_position_covariance=cov)
    return out


def weighted_cdf(nodes,weights,thresholds):
    """CDF of weighted finite node atoms, not a node-bin continuous CDF."""
    x=np.asarray(nodes,dtype=np.longdouble);t=np.asarray(thresholds,dtype=np.longdouble)
    if x.ndim!=1 or not np.isfinite(x).all() or np.isnan(t).any():raise ValueError('Finite1D nodes and non-NaN thresholds required')
    w=_weights(weights,len(x));flat=t.ravel()
    values=np.array([np.sum(w[x<=z]) for z in flat],dtype=np.longdouble).reshape(t.shape)
    return dict(values=values,interpretation='CDF of finite quadrature-node atoms; not certified continuous CDF or continuous node-bin CDF')


def weighted_crps(nodes,weights,outcome):
    """CRPS for the finite atomic forecast, O(NlogN); no continuous claim."""
    x=np.asarray(nodes,dtype=np.longdouble)
    if x.ndim!=1 or not np.isfinite(x).all() or not np.isfinite(outcome):raise ValueError('Finite1D nodes and scalar outcome required')
    w=_weights(weights,len(x));order=np.argsort(x,kind='stable');x,w=x[order],w[order]
    before_w=np.cumsum(w)-w;before_x=np.cumsum(w*x)-w*x
    half_pairwise=np.sum(w*(x*before_w-before_x))
    return np.sum(w*np.abs(x-outcome))-half_pairwise


def brier_score(event_probability,event_outcome):
    p=np.asarray(event_probability,dtype=np.longdouble);o=np.asarray(event_outcome)
    if not np.isfinite(p).all() or np.any((p<0)|(p>1)) or np.any((o!=0)&(o!=1)):
        raise ValueError('Probability in[0,1] and binary observed outcome required')
    return (p-o)**2
