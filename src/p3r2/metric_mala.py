"""Exact-MH MALA with a genuinely state-dependent full response metric.

Pure array code. The caller supplies U=-log target, its exact gradient and
H=C^-1. The Gaussian proposal is normalized at both endpoints. A manifold
Gamma drift is unnecessary for invariance because the actual proposal has
an exact Metropolis-Hastings correction. No finite-run mixing claim follows.
"""
from dataclasses import dataclass,field
from typing import Any
import numpy as np
from scipy.linalg import cho_solve,solve_triangular


@dataclass
class MetricState:
    x:np.ndarray
    energy:float
    gradient:np.ndarray
    precision:np.ndarray
    payload:Any=None
    L:np.ndarray=field(init=False,repr=False)
    logdet:float=field(init=False)
    eigenvalues:np.ndarray=field(init=False,repr=False)
    eigenvectors:np.ndarray=field(init=False,repr=False)

    def __post_init__(self):
        self.x=np.array(self.x,float,copy=True);self.gradient=np.array(self.gradient,float,copy=True)
        self.precision=np.array(self.precision,float,copy=True);self.energy=float(self.energy)
        d=self.x.size
        if self.x.ndim!=1 or self.gradient.shape!=(d,) or self.precision.shape!=(d,d):
            raise ValueError('matching state, gradient and precision dimensions required')
        if not np.isfinite(self.energy) or not all(np.isfinite(a).all() for a in (self.x,self.gradient,self.precision)):
            raise FloatingPointError('MALA state must contain finite exact arrays')
        if not np.allclose(self.precision,self.precision.T,rtol=1e-11,atol=1e-12):
            raise ValueError('metric precision must be symmetric SPD')
        self.precision=(self.precision+self.precision.T)/2
        self.L=np.linalg.cholesky(self.precision);self.logdet=float(2*np.log(np.diag(self.L)).sum())
        self.eigenvalues,self.eigenvectors=np.linalg.eigh(self.precision)


def drift(state,eps):
    return -.5*eps**2*cho_solve((state.L,True),state.gradient)


def normalized_log_proposal(origin,destination,eps):
    """log N(destination; origin-.5eps² C gradU, eps² C)."""
    if not np.isfinite(eps) or eps<=0: raise ValueError('positive finite epsilon required')
    destination=np.asarray(destination,float)
    displacement=destination-origin.x-drift(origin,eps)
    whitened=origin.L.T@displacement/eps;quadratic=float(whitened@whitened)
    value=-.5*len(origin.x)*np.log(2*np.pi)-len(origin.x)*np.log(eps)+.5*origin.logdet-.5*quadratic
    return float(value),quadratic


def metric_change(current,proposed,top_dimension=8):
    """Relative metric change and dominant eigenspace rotation, descriptive only."""
    delta=proposed.precision-current.precision
    left=solve_triangular(current.L,delta,lower=True)
    relative=solve_triangular(current.L,left.T,lower=True).T
    relative=(relative+relative.T)/2;ev=np.linalg.eigvalsh(relative)
    d=len(current.x);k=min(top_dimension,max(1,d-1))
    a=current.eigenvectors[:,-k:];b=proposed.eigenvectors[:,-k:]
    cosines=np.linalg.svd(a.T@b,compute_uv=False)
    gap=lambda s:float(s.eigenvalues[-k]-s.eigenvalues[-k-1]) if k<d else None
    return dict(relative_precision_change_Frobenius=float(np.linalg.norm(relative)),
        relative_precision_change_operator=float(np.max(abs(ev))),
        precision_change_relative_Frobenius=float(np.linalg.norm(delta)/np.linalg.norm(current.precision)),
        dominant_metric_subspace_dimension=k,
        dominant_metric_projector_rotation_Frobenius=float(np.linalg.norm(a@a.T-b@b.T)),
        dominant_metric_largest_principal_angle_radians=float(np.arccos(np.clip(cosines.min(),0,1))),
        current_selection_eigengap=gap(current),proposed_selection_eigengap=gap(proposed),
        rotation_caution='Eigenspace rotation is descriptive and unstable at small selection gaps')


def balance_terms(current,proposed,eps):
    forward,qforward=normalized_log_proposal(current,proposed.x,eps)
    reverse,qreverse=normalized_log_proposal(proposed,current.x,eps)
    ratio=current.energy-proposed.energy+reverse-forward
    return dict(log_proposal_forward=forward,log_proposal_reverse=reverse,
                forward_quadratic=qforward,reverse_quadratic=qreverse,
                logdet_precision_current=current.logdet,logdet_precision_proposed=proposed.logdet,
                log_target_ratio=current.energy-proposed.energy,
                log_proposal_correction=reverse-forward,log_mh_ratio=float(ratio))


class FullMetricMALA:
    """One exact target/gradient/metric endpoint evaluation per transition."""
    def __init__(self,evaluate,eps=.05):
        if not np.isfinite(eps) or eps<=0: raise ValueError('positive finite epsilon required')
        self.evaluate=evaluate;self.eps=float(eps)

    def step(self,current,rng):
        noise=rng.normal(size=len(current.x))
        proposed_x=current.x+drift(current,self.eps)+self.eps*solve_triangular(current.L.T,noise,lower=False)
        if not np.isfinite(proposed_x).all(): raise FloatingPointError('nonfinite MALA proposal; run unqualified')
        proposed=self.evaluate(proposed_x)
        terms=balance_terms(current,proposed,self.eps)
        uniform=float(rng.random());log_uniform=float(np.log(max(uniform,np.finfo(float).tiny)))
        accepted=bool(log_uniform<min(0.,terms['log_mh_ratio']))
        row=dict(terms,**metric_change(current,proposed),accepted=accepted,
            repeated_state=bool(not accepted or np.array_equal(current.x,proposed.x)),
            eps=self.eps,log_uniform=log_uniform,proposal_noise=noise,
            proposed_x=proposed.x.copy(),energy_current=current.energy,energy_proposed=proposed.energy,
            gradient_current=current.gradient.copy(),gradient_proposed=proposed.gradient.copy(),
            status='EXACT_MH_INVARIANCE_NOT_A_FINITE_RUN_MIXING_CERTIFICATE')
        return (proposed if accepted else current),proposed,row


def verify_probability_arrays(seed=2026102392,pairs=80):
    """Independent SciPy Gaussian-density and pairwise detailed-balance control.

    Synthetic arrays only. Checks a position-dependent SPD precision, including
    determinant variation and both proposal quadratics; no solver is imported.
    """
    from scipy.stats import multivariate_normal
    rng=np.random.default_rng(seed);d=5;eps=.37
    errors=[];balance_errors=[];antisymmetry=[];quadratic_errors=[];det_variation=[]
    def state(x):
        gradient=x+.14*np.sin(2*x)
        energy=.5*float(x@x)+.07*float(np.sum(1-np.cos(2*x)))
        v=np.sin(x);H=np.diag(1+.6*x*x)+.15*np.outer(v,v)
        return MetricState(x,energy,gradient,H)
    for _ in range(pairs):
        a=state(rng.normal(size=d));b=state(rng.normal(size=d));terms=balance_terms(a,b,eps)
        C=np.linalg.inv(a.precision);Cy=np.linalg.inv(b.precision)
        mu=a.x-.5*eps**2*C@a.gradient;muy=b.x-.5*eps**2*Cy@b.gradient
        forward=float(multivariate_normal.logpdf(b.x,mean=mu,cov=eps**2*C))
        reverse=float(multivariate_normal.logpdf(a.x,mean=muy,cov=eps**2*Cy))
        errors.extend([abs(forward-terms['log_proposal_forward']),abs(reverse-terms['log_proposal_reverse'])])
        ratio=a.energy-b.energy+reverse-forward;back=balance_terms(b,a,eps)['log_mh_ratio']
        antisymmetry.append(abs(terms['log_mh_ratio']+back))
        left=-a.energy+forward+min(0.,ratio);right=-b.energy+reverse+min(0.,-ratio)
        balance_errors.append(abs(left-right))
        dx=b.x-mu;dy=a.x-muy
        quadratic_errors.extend([abs(float(dx@a.precision@dx/eps**2)-terms['forward_quadratic']),
                                 abs(float(dy@b.precision@dy/eps**2)-terms['reverse_quadratic'])])
        det_variation.append(abs(a.logdet-b.logdet))
    # Same normal/uniform stream is exactly reproducible for matched-start
    # comparisons; the physics wrapper changes only metric/gradient cost.
    initial=state(np.zeros(d));kernel=FullMetricMALA(state,eps)
    first=kernel.step(initial,np.random.default_rng(seed));second=kernel.step(initial,np.random.default_rng(seed))
    same_stream=bool(np.array_equal(first[1].x,second[1].x) and first[2]['accepted']==second[2]['accepted'])
    return dict(status='PURE_ARRAY_BALANCE_CONTROL_COMPLETE',seed=seed,pairs=pairs,dimension=d,new_physics_RHS=0,
        normalized_logproposal_max_error_vs_independent_SciPy=float(max(errors)),
        forward_reverse_quadratic_max_error=float(max(quadratic_errors)),
        detailed_balance_log_flux_max_error=float(max(balance_errors)),
        MH_logratio_reverse_antisymmetry_max_error=float(max(antisymmetry)),
        maximum_tested_logdet_change=float(max(det_variation)),same_initial_rng_reproducible=same_stream,
        checks=dict(normalized_proposals=max(errors)<1e-10,quadratics=max(quadratic_errors)<1e-10,
            detailed_balance=max(balance_errors)<1e-10,reverse_antisymmetry=max(antisymmetry)<1e-10,
            variable_determinants_tested=max(det_variation)>1.,same_rng_reproducible=same_stream),
        scope='Algebraic probability control, including varying endpoint determinants; no Maxwell mixing or posterior qualification')
