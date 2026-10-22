"""Saved-matrix primal/dual decision bounds; no physics backend or book.

The algebraic enclosures are conditional on the supplied stability lower bound.
Floating-point inflation uses a standard gamma_k model for complex dot products
and absolute products. Assumptions: IEEE float64 round-to-nearest arithmetic,
ordinary stable dot-product evaluation (including BLAS), finite intermediates,
and no material underflow/overflow. These assumptions are not verified by a
rounding-mode-controlled interval library. Results therefore explicitly say
CONDITIONAL_FLOATING_POINT_BOUND, never globally interval-certified.

A numeric matrix may be stable when a sufficient Neumann bound is unavailable.
No full reference prediction is accepted by these evaluators to choose fallback.
"""
import numpy as np


UNIT_ROUNDOFF = np.finfo(float).eps/2
BOUND_LABEL = 'CONDITIONAL_FLOATING_POINT_BOUND_NOT_INTERVAL_CERTIFIED'


def _gamma(k):
    t = int(k)*UNIT_ROUNDOFF
    return t/(1-t) if t < 1 else np.inf


def _up(value, operations=8):
    with np.errstate(over='ignore',invalid='ignore'):
        return float(np.nextafter(float(value)*(1+_gamma(operations)),np.inf))


def _down_positive(value, operations=8):
    return max(0.,float(np.nextafter(float(value)/(1+_gamma(operations)),0.)))


def _matrix(value, name):
    value = np.asarray(value,complex)
    if value.ndim!=2 or not np.all(np.isfinite(value)):
        raise ValueError(f'{name} must be a finite matrix')
    return value


def _counts():
    return dict(matrix_products=0,absolute_matrix_products=0,
                induced_norm_upper_evaluations=0,frobenius_norm_upper_evaluations=0,
                projected_solve_rhs_columns=0,projected_linear_solves=0,
                residual_center_norm_evaluations=0,matrix_actions=[],
                new_maxwell_solves=0,new_physical_rhs=0,new_physical_factorizations=0)


def _abs_upper(value):
    # Complex absolute values involve squares/square-root. Multiplicative
    # gamma inflation assumes normal finite arithmetic, as stated above.
    return np.abs(value)*(1+_gamma(4))


def _product(left,right,counts,label):
    k = left.shape[1]
    product = left@right
    # gamma_(8k+4) conservatively accounts for real operations in a complex
    # dot product; the positive product is inflated for its own evaluation.
    positive = _abs_upper(left)@_abs_upper(right)
    with np.errstate(over='ignore',invalid='ignore'):
        error = _gamma(8*k+4)*positive/(1-_gamma(2*k+4))
    counts['matrix_products']+=1
    counts['absolute_matrix_products']+=1
    counts['matrix_actions'].append(dict(label=label,left_shape=list(left.shape),
                                         rhs_columns=right.shape[1],inner_dimension=k))
    return product,error


def _difference(left,right,right_error=0.):
    computed = left-right
    error = right_error+_gamma(4)*(_abs_upper(left)+_abs_upper(right))
    return computed,error


def _induced_upper(value,counts,entry_error=0.):
    """||A||2 <= sqrt(||A||1 ||A||inf), with entrywise error inflation."""
    counts['induced_norm_upper_evaluations']+=1
    positive = _abs_upper(value)+entry_error
    one = _up(np.max(np.sum(positive,axis=0)),value.shape[0]+8)
    infinity = _up(np.max(np.sum(positive,axis=1)),value.shape[1]+8)
    # sqrt(a)*sqrt(b) avoids overflowing a*b when a,b are representable.
    return _up(np.sqrt(one)*np.sqrt(infinity))


def _frobenius_upper(value,counts,entry_error=0.):
    counts['frobenius_norm_upper_evaluations']+=1
    positive = _abs_upper(value)+entry_error
    scale = float(np.max(positive)) if positive.size else 0.
    if scale==0:
        return 0.
    if not np.isfinite(scale):
        return np.inf
    return _up(scale*np.sqrt(np.sum((positive/scale)**2)),2*positive.size+8)


def _result_unavailable(reason,counts,**values):
    return dict(status='BOUND_UNAVAILABLE',reason=reason,bound_label=BOUND_LABEL,
                beta=None,counts=counts,**values)


def stability_from_saved_inverse(L0,R0,Lnew):
    """Conditional lower bound for sigma_min(Lnew) from a saved inverse.

    rho0 bounds ||I-R0 L0||2. If rho0<1 then B0=||R0||2/(1-rho0)
    bounds ||L0^-1||2. kappa=B0||Lnew-L0||2<1 gives
    beta=(1-kappa)/B0. No inverse or solve is performed here. 'Unavailable'
    only means the sufficient conditions were not established.
    """
    L0,R0,Lnew = (_matrix(a,n) for a,n in [(L0,'L0'),(R0,'R0'),(Lnew,'Lnew')])
    if L0.shape[0]!=L0.shape[1] or not L0.shape[0] or R0.shape!=L0.shape or Lnew.shape!=L0.shape:
        raise ValueError('matching nonempty square matrices required')
    counts = _counts()
    product,product_error = _product(R0,L0,counts,'SAVED_INVERSE_RESIDUAL_R0_L0')
    residual,residual_error = _difference(np.eye(len(L0)),product,product_error)
    rho0 = _induced_upper(residual,counts,residual_error)
    r0norm = _induced_upper(R0,counts)
    delta,delta_error = _difference(Lnew,L0)
    delta_norm = _induced_upper(delta,counts,delta_error)
    values = dict(rho0_upper=rho0,saved_inverse_norm_upper=r0norm,
                  delta_norm_upper=delta_norm,kappa=None,
                  inverse0_norm_upper=None,inverse_new_norm_upper=None,
                  numerical_assumptions='gamma_k complex dot products; finite normal float64 arithmetic; no directed interval implementation')
    if not np.isfinite(rho0) or rho0>=1 or not np.isfinite(r0norm) or r0norm<=0:
        return _result_unavailable('SAVED_INVERSE_RESIDUAL_SUFFICIENT_CONDITION_FAILED_NOT_SINGULARITY',counts,**values)
    denominator = _down_positive(1-rho0)
    inverse0 = _up(r0norm/denominator) if denominator>0 else np.inf
    kappa = _up(inverse0*delta_norm)
    values.update(inverse0_norm_upper=inverse0,kappa=kappa)
    if not np.isfinite(kappa) or kappa>=1:
        return _result_unavailable('NEUMANN_SUFFICIENT_CONDITION_FAILED_NOT_SINGULARITY',counts,**values)
    denominator = _down_positive(1-kappa)
    inverse_new = _up(inverse0/denominator) if denominator>0 else np.inf
    beta = _down_positive(1/inverse_new) if np.isfinite(inverse_new) and inverse_new>0 else 0.
    values['inverse_new_norm_upper']=inverse_new
    if not np.isfinite(beta) or beta<=0:
        return _result_unavailable('FINITE_POSITIVE_STABILITY_BOUND_NOT_ESTABLISHED',counts,**values)
    return dict(status='SUFFICIENT_STABILITY_BOUND_AVAILABLE',beta=beta,
                bound_label=BOUND_LABEL,counts=counts,**values)


def dual_corrected_prediction(L,rhs,S,V,Z0,beta,*,observation=None,sigma=None):
    """Projected state plus primal/dual output-error algebra on saved arrays.

    Shapes: L(n,n), rhs(n,m), S(p,n), V(n,r), Z0(n,p). Galerkin state is
    V(V*LV)^-1 V*rhs; its numeric solve error needs no separate assumption
    because the full primal residual is bounded afterwards. Z0 is an arbitrary
    supplied approximate receiver adjoint. beta must be an independently
    justified positive lower bound on sigma_min(L), or a stability-result dict.

    Returns fc=S j+Z0*rp and eps on the Frobenius complex output error.
    Roundoff in residuals and in fc is included. observation and sigma must
    be supplied together to form a real-whitened likelihood interval.
    """
    L,rhs,S,V,Z0 = (_matrix(a,n) for a,n in [(L,'L'),(rhs,'rhs'),(S,'S'),(V,'V'),(Z0,'Z0')])
    n=L.shape[0]
    if (L.shape!=(n,n) or n==0 or rhs.shape[0]!=n or S.shape[1]!=n
            or V.shape[0]!=n or not 0<V.shape[1]<=n or Z0.shape!=(n,S.shape[0])
            or rhs.shape[1]==0 or S.shape[0]==0):
        raise ValueError('incompatible L/rhs/S/V/Z0 dimensions')
    if (observation is None)!=(sigma is None):
        raise ValueError('observation and sigma must be supplied together')
    if sigma is not None and (not np.isfinite(sigma) or sigma<=0):
        raise ValueError('positive finite sigma required')
    y = _matrix(observation,'observation') if observation is not None else None
    if y is not None and y.shape!=(S.shape[0],rhs.shape[1]):
        raise ValueError('observation output layout mismatch')
    counts=_counts()
    supplied_status=beta.get('status') if isinstance(beta,dict) else 'CALLER_SUPPLIED'
    beta=beta.get('beta') if isinstance(beta,dict) else beta
    if beta is None or not np.isfinite(beta) or beta<=0:
        return _result_unavailable('NO_POSITIVE_SUPPLIED_STABILITY_LOWER_BOUND',counts,
                                  prediction=None,eps_complex=None,eps_white=None,likelihood_interval=None)
    LV,_=_product(L,V,counts,'L_ON_TRIAL_BASIS')
    reduced,_=_product(V.conj().T,LV,counts,'GALERKIN_OPERATOR')
    reduced_rhs,_=_product(V.conj().T,rhs,counts,'GALERKIN_RHS')
    try:
        coefficients=np.linalg.solve(reduced,reduced_rhs)
        counts['projected_solve_rhs_columns']=rhs.shape[1]
        counts['projected_linear_solves']=1
    except np.linalg.LinAlgError:
        return _result_unavailable('PROJECTED_OPERATOR_SOLVE_FAILED_NOT_FULL_MATRIX_SINGULARITY',counts,
                                  prediction=None,eps_complex=None,eps_white=None,likelihood_interval=None)
    state,_=_product(V,coefficients,counts,'GALERKIN_STATE')
    Lj,Lj_error=_product(L,state,counts,'PRIMAL_RESIDUAL_ACTION')
    rp,rp_error=_difference(rhs,Lj,Lj_error)
    LstarZ,LstarZ_error=_product(L.conj().T,Z0,counts,'DUAL_RESIDUAL_ACTION')
    Rd,Rd_error=_difference(S.conj().T,LstarZ,LstarZ_error)
    primal_norm=_frobenius_upper(rp,counts,rp_error)
    dual_norm=_frobenius_upper(Rd,counts,Rd_error)
    primal_roundoff=_frobenius_upper(np.zeros_like(rp),counts,rp_error)
    output,output_error=_product(S,state,counts,'PROJECTED_OUTPUT')
    correction,correction_error=_product(Z0.conj().T,rp,counts,'DUAL_OUTPUT_CORRECTION')
    fc=output+correction
    fc_error=output_error+correction_error+_gamma(4)*(_abs_upper(output)+_abs_upper(correction))
    fc_arithmetic=_frobenius_upper(np.zeros_like(fc),counts,fc_error)
    adjoint_norm=_induced_upper(Z0.conj().T,counts)
    eps_residual=_up(dual_norm*primal_norm/float(beta))
    eps_roundoff=_up(fc_arithmetic+_up(adjoint_norm*primal_roundoff))
    eps=_up(eps_residual+eps_roundoff)
    if not np.isfinite(eps) or not np.all(np.isfinite(fc)):
        return _result_unavailable('NONFINITE_ERROR_ARITHMETIC',counts,
                                  prediction=None,eps_complex=None,eps_white=None,likelihood_interval=None)
    result=dict(status='OUTPUT_BOUND_AVAILABLE',bound_label=BOUND_LABEL,
                beta=float(beta),supplied_stability_status=supplied_status,
                prediction=fc,projected_state=state,primal_residual=rp,dual_residual=Rd,
                primal_residual_norm_upper=primal_norm,dual_residual_norm_upper=dual_norm,
                eps_complex=eps,eps_residual_term=eps_residual,eps_roundoff_term=eps_roundoff,
                eps_white=None,likelihood_interval=None,counts=counts)
    if y is not None:
        delta,delta_error=_difference(fc,y)
        # Computed residual-center norm has its own arithmetic uncertainty;
        # include it as extra output radius before whitening.
        residual_center_roundoff=_frobenius_upper(np.zeros_like(delta),counts,delta_error)
        residual_center=float(np.linalg.norm(delta))
        counts['residual_center_norm_evaluations']+=1
        norm_arithmetic=_gamma(2*delta.size+8)*residual_center
        effective_eps=_up(eps+residual_center_roundoff+norm_arithmetic)
        scale=_up(np.sqrt(2)/float(sigma))
        eps_white=_up(scale*effective_eps)
        center_white=scale*residual_center
        # Inflate center multiplication uncertainty too; ensure lower and
        # upper endpoints move outward under the stated arithmetic model.
        eps_white=_up(eps_white+_gamma(16)*abs(center_white))
        lo=_down_positive(max(0.,center_white-eps_white))
        hi=_up(center_white+eps_white)
        lower=_down_positive(.5*lo*lo)
        upper=_up(.5*hi*hi)
        if not np.isfinite(upper):
            result.update(status='BOUND_UNAVAILABLE',reason='NONFINITE_LIKELIHOOD_INTERVAL')
        else:
            result.update(eps_white=eps_white,white_residual_center_norm=center_white,
                          likelihood_interval=(lower,upper),sigma=float(sigma))
    return result


def same_uniform_decision(current_interval,candidate_interval,logprior_proposal_ratio,logu):
    """Compare the SAME supplied log uniform against a bounded MH ratio.

    Caller must retain proposal and logu across REFINE calls. This function
    neither generates randomness nor changes the proposal. Unavailable
    intervals yield REFINE. Boundary equality is conservatively unresolved.
    """
    c=float(logprior_proposal_ratio);logu=float(logu)
    if not np.isfinite(c) or not np.isfinite(logu) or logu>=0:
        raise ValueError('finite log prior/proposal ratio and finite logu<0 required')
    base=dict(logu=logu,logprior_proposal_ratio=c,retain_same_uniform=True,
              bound_label=BOUND_LABEL)
    if current_interval is None or candidate_interval is None:
        return dict(action='REFINE',reason='LIKELIHOOD_BOUND_UNAVAILABLE',log_ratio_interval=None,**base)
    intervals=[]
    for interval in [current_interval,candidate_interval]:
        a=np.asarray(interval,float)
        if a.shape!=(2,) or not np.all(np.isfinite(a)) or a[0]<0 or a[0]>a[1]:
            raise ValueError('likelihood interval must have finite 0<=lower<=upper')
        intervals.append(a)
    current,candidate=intervals
    raw_lower=c+current[0]-candidate[1]
    raw_upper=c+current[1]-candidate[0]
    arithmetic=_gamma(8)*(abs(c)+np.max(np.abs(current))+np.max(np.abs(candidate)))
    lower=float(np.nextafter(raw_lower-arithmetic,-np.inf))
    upper=float(np.nextafter(raw_upper+arithmetic,np.inf))
    if not np.isfinite(lower) or not np.isfinite(upper):
        return dict(action='REFINE',reason='NONFINITE_LOG_RATIO_BOUND',log_ratio_interval=None,**base)
    action='ACCEPT' if logu<lower else 'REJECT' if logu>upper else 'REFINE'
    return dict(action=action,reason='SHARED_UNIFORM_INTERVAL_COMPARISON',
                log_ratio_interval=(lower,upper),**base)
