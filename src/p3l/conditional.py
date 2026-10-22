"""Frozen Gaussian64 conditional GN references and explicit likelihood correction.

No solver is instantiated here. Every Target.evaluate/jacobian invocation is
visible to its caller-owned accounting. All properties remain native uniforms.
Frames alter coordinates, never the original Gaussian prior or its measure.
"""
import numpy as np
from scipy.linalg import eigh, solve_triangular, cho_solve
from scipy.optimize import least_squares
from scipy.special import logsumexp

from .moments import gaussian_material_moments

LOG_2PI = float(np.log(2*np.pi))


def _native(q):
    q = np.asarray(q, float)
    if q.shape != (66,) or not np.all(np.isfinite(q)) or np.any(q[64:] < 0) or np.any(q[64:] > 1):
        raise ValueError("finite Gaussian64 and native box properties required")
    return q


def make_frame(q0, jacobian, weak_dimension, *, kind="prior_relative",
               material_derivative=None, material_weights=None, material_ridge=None):
    """Array-only fixed frame construction; no target calls or adaptive refresh.

    prior_relative selects smallest eigenvalues of Jz.T Jz.
    material_generalized selects smallest generalized eigenvalues relative to
    a regularized material Gram matrix, then orthonormalizes the chosen span.
    posterior_material_risk selects leading C^.5 M C^.5 eigenvectors mapped by
    C^.5, C=(I+Jz.T Jz)^-1, then orthonormalizes. These are declared frozen
    heuristics; no optimality, posterior validity, or novelty is asserted.
    """
    q0 = _native(q0).copy()
    J = np.asarray(jacobian,float)
    if J.ndim != 2 or J.shape[1] not in (64,66) or not np.all(np.isfinite(J)):
        raise ValueError("finite measurement Jacobian with 64 or 66 columns required")
    if not isinstance(weak_dimension,int) or not 0 <= weak_dimension < 64:
        raise ValueError("weak_dimension must be an integer in [0,63]")
    J = J[:,:64]
    Hdata = J.T@J
    Hdata = .5*(Hdata+Hdata.T)
    eigenvalues, U = eigh(Hdata)
    ridge_used = 0.
    if kind == "prior_relative":
        selected = U[:,:weak_dimension]
        selected_values = eigenvalues[:weak_dimension]
    elif kind in ("material_generalized", "posterior_material_risk"):
        D = np.asarray(material_derivative)
        if D.ndim != 2 or D.shape[1] not in (64,66) or not np.all(np.isfinite(D)):
            raise ValueError("finite complex material derivative with 64 or 66 columns required")
        D = D[:,:64]
        w = (np.full(D.shape[0],1/D.shape[0]) if material_weights is None
             else np.asarray(material_weights,float))
        if w.shape != (D.shape[0],) or not np.all(np.isfinite(w)) or np.any(w < 0):
            raise ValueError("nonnegative finite per-material-cell weights required")
        M = np.real(D.conj().T@(w[:,None]*D))
        M = .5*(M+M.T)
        if kind == "material_generalized":
            ridge_used = max(float(np.trace(M))/64*1e-8,1e-12) if material_ridge is None else float(material_ridge)
            if not np.isfinite(ridge_used) or ridge_used <= 0:
                raise ValueError("positive generalized-material metric ridge required")
            vals, vectors = eigh(Hdata,M+ridge_used*np.eye(64))
            selected, selected_values = vectors[:,:weak_dimension], vals[:weak_dimension]
        else:
            # Negative Hdata eigenvalues at ordinary roundoff are clipped.
            Croot = (U*(1/np.sqrt(1+np.maximum(eigenvalues,0))))@U.T
            risk = Croot@M@Croot
            vals, vectors = eigh(.5*(risk+risk.T))
            selected = Croot@vectors[:,-weak_dimension:][:,::-1] if weak_dimension else np.zeros((64,0))
            selected_values = vals[-weak_dimension:][::-1] if weak_dimension else np.zeros(0)
    else:
        raise ValueError("unknown fixed-frame kind")
    # Complete QR gives a unit-Jacobian chart, even for generalized vectors.
    Q = np.eye(64) if weak_dimension == 0 else np.linalg.qr(selected,mode="complete")[0]
    return dict(q0=q0, Vw=Q[:,:weak_dimension], Vs=Q[:,weak_dimension:],
        weak_dimension=weak_dimension, strong_dimension=64-weak_dimension,
        kind=kind, selection_eigenvalues=selected_values,
        generalized_metric_ridge=ridge_used,
        orthogonality_error=float(np.linalg.norm(Q.T@Q-np.eye(64))),
        gaussian_prior_weak_center=-Q[:,:weak_dimension].T@q0[:64],
        gaussian_prior_strong_center=-Q[:,weak_dimension:].T@q0[:64],
        reference_measure="da db du; absolute chart determinant magnitude 1",
        status="FROZEN_ARRAY_DERIVED_FRAME_NOT_A_POSTERIOR_CERTIFICATE")


def _parts(frame,a,u):
    Vw, Vs = np.asarray(frame["Vw"],float), np.asarray(frame["Vs"],float)
    q0 = _native(frame["q0"])
    a, u = np.asarray(a,float), np.asarray(u,float)
    if a.shape != (Vw.shape[1],) or u.shape != (2,) or not np.all(np.isfinite(a)) or not np.all(np.isfinite(u)) or np.any(u < 0) or np.any(u > 1):
        raise ValueError("finite weak coordinates and two native properties required")
    if Vw.shape[0] != 64 or Vs.shape != (64,64-Vw.shape[1]):
        raise ValueError("invalid frame dimensions")
    Q = np.c_[Vw,Vs]
    if not np.all(np.isfinite(Q)) or not np.allclose(Q.T@Q,np.eye(64),atol=1e-10,rtol=1e-10):
        raise ValueError("fixed Gaussian frame must be Euclidean orthonormal")
    return q0[:64]+Vw@a, Vs, u


def compose(frame,a,u,b):
    base, Vs, u = _parts(frame,a,u)
    b = np.asarray(b,float)
    if b.shape != (Vs.shape[1],) or not np.all(np.isfinite(b)):
        raise ValueError("finite strong coordinates required")
    return np.r_[base+Vs@b,u]


def _loggamma(target,point,y,region=None):
    if region is not None and not bool(region(point.q.copy())):
        return -np.inf
    return -float(target.phi(point,y))-.5*float(point.q[:64]@point.q[:64])-32*LOG_2PI


def _event(callback,event,**fields):
    if callback is not None:
        callback(dict(event=event,**fields))


def evaluate_slice(target,y,frame,a,u,*,b=None,point=None,region=None,callback=None):
    """A fixed-strong slice value, explicitly not the original weak marginal."""
    if b is None:
        b = np.zeros(frame["strong_dimension"])
    q = compose(frame,a,u,b)
    calls = 0
    if point is None:
        _event(callback,"EVALUATION_ATTEMPT",q=q,purpose="fixed_strong_slice")
        point = target.evaluate(q,purpose="p3l_fixed_strong_slice")
        calls = 1
    elif not np.array_equal(point.q,q):
        raise ValueError("explicit slice point does not match supplied coordinates")
    lg = _loggamma(target,point,y,region)
    _event(callback,"SLICE",q=q,endpoint=point.key,
           log_gamma=float(lg) if np.isfinite(lg) else None,zero_target_density=not np.isfinite(lg))
    return dict(point=point,q=q,log_gamma=lg,target_family="SLICE_PROFILE",
        target_subtype="SLICE",strong_direction_treatment="fixed",
        reference_measure="da du at explicitly fixed b", approximation_status="EXACT_DENSITY_EVALUATION_ON_ARTIFICIAL_SLICE",
        target_evaluations=calls,original_weak_marginal=False)


def _factor(H):
    H = .5*(H+H.T)
    fallback = None
    try:
        L = np.linalg.cholesky(H)
    except np.linalg.LinAlgError:
        vals, vecs = eigh(H)
        floor = max(1.,64*np.finfo(float).eps*max(1.,float(np.max(np.abs(vals)))))
        H = (vecs*np.maximum(vals,floor))@vecs.T
        H = .5*(H+H.T)
        L = np.linalg.cholesky(H)
        fallback = dict(reason="GN_CHOLESKY_FAILED",eigenvalue_floor=floor,
                        original_minimum_eigenvalue=float(vals.min()))
    return H,L,fallback


def reference_logpdf(reference,b):
    b = np.asarray(b,float)
    mean = np.asarray(reference["mean_b"],float)
    if b.shape[-1:] != mean.shape or not np.all(np.isfinite(b)):
        raise ValueError("finite strong coordinates with matching final axis required")
    delta = b-mean
    whitened = delta@reference["precision_cholesky"]
    return .5*reference["logdet_precision"]-.5*mean.size*LOG_2PI-.5*np.sum(whitened**2,axis=-1)


def sample_reference(reference,samples,rng):
    if not isinstance(samples,int) or samples < 1:
        raise ValueError("positive integer samples required")
    mean = np.asarray(reference["mean_b"],float)
    noise = rng.normal(size=(samples,mean.size))
    return mean+solve_triangular(reference["precision_cholesky"].T,noise.T,lower=False).T


def conditional_material_moments(prior,frame,a,u,reference):
    """Exact decoder marginal moments under the specified Gaussian reference."""
    base, Vs, u = _parts(frame,a,u)
    zmean = base+Vs@reference["mean_b"]
    A = np.asarray(prior.A,float)
    if A.ndim != 2 or A.shape[1] != 64:
        raise ValueError("original prior.A must map Gaussian64 to field cells")
    field_mean = np.asarray(prior.config["field_mean"],float)+A@zmean
    spread = solve_triangular(reference["precision_cholesky"],(A@Vs).T,lower=True)
    field_sd = np.linalg.norm(spread,axis=0)
    real_interval, imag_interval = np.asarray(prior.real_interval), np.asarray(prior.imag_interval)
    delta = complex(real_interval[0]+np.diff(real_interval)[0]*u[0],
                    imag_interval[0]+np.diff(imag_interval)[0]*u[1])
    moments = gaussian_material_moments(field_mean,field_sd,delta,host=prior.host,
        mask=prior.mask,halfwidth=prior.temperature)
    return dict(mean=moments.mean,variance=moments.variance,
        second_abs_moment=moments.second_abs_moment,fraction_mean=moments.fraction.mean,
        field_mean=field_mean,field_sd=field_sd,properties=u.copy(),
        target_family="GAUSSIAN_CONDITIONAL_REFERENCE",
        reference_measure="normalized Gaussian db; native fixed properties",
        approximation_status="DECODER_MOMENTS_EXACT_ONLY_UNDER_DECLARED_GAUSSIAN_REFERENCE",
        spatial_covariance_returned=False,local_truncation_accounted=False)


def conditional_reference(target,y,frame,a,u,*,max_evaluations=6,initial_b=None,
                          initial_point=None,callback=None,target_family="FULL",region=None):
    """At most 4--8 fresh full states for a conditional GN mode/reference.

    The optimizer uses the FULL conditional energy even when a LOCAL_C region
    is provided for subsequent correction. Its untruncated GN volume never
    estimates local retained probability mass. Production parameters are fixed.
    """
    if not isinstance(max_evaluations,int) or not 4 <= max_evaluations <= 8:
        raise ValueError("conditional max_evaluations must be an integer in [4,8]")
    if (target_family not in ("FULL","LOCAL_C") or
            (target_family == "LOCAL_C" and region is None) or
            (target_family == "FULL" and region is not None)):
        raise ValueError("FULL or explicitly region-defined LOCAL_C target required")
    base, Vs, u = _parts(frame,a,u)
    s = Vs.shape[1]
    b0 = np.zeros(s) if initial_b is None else np.asarray(initial_b,float)
    qinitial = compose(frame,a,u,b0)
    cache = {}
    best = [None,np.inf,None]
    counts = dict(target_evaluations=0,jacobian_evaluations=0,fun_calls=0,jac_calls=0)
    def remember(p,b):
        r = np.r_[target.white(p.state.data-y),p.q[:64]]
        if not np.all(np.isfinite(r)):
            raise FloatingPointError("nonfinite conditional residual")
        cost = .5*float(r@r)
        cache[np.asarray(b,float).tobytes()] = (p,r)
        if cost < best[1]:
            best[:] = [p,cost,np.asarray(b,float).copy()]
        return p,r
    if initial_point is not None:
        if not np.array_equal(initial_point.q,qinitial):
            raise ValueError("initial point must exactly match declared (a,u,b)")
        remember(initial_point,b0)
    def get(b):
        key = np.asarray(b,float).tobytes()
        if key not in cache:
            if counts["target_evaluations"] >= max_evaluations:
                raise RuntimeError("conditional fresh-state cap")
            q = compose(frame,a,u,b)
            counts["target_evaluations"] += 1
            _event(callback,"EVALUATION_ATTEMPT",q=q,purpose="conditional_mode")
            p = target.evaluate(q,purpose="p3l_conditional_GN_mode")
            remember(p,b)
            actual_residual = cache[key][1]
            _event(callback,"EVALUATION",q=q,endpoint=p.key,
                   cost=.5*float(actual_residual@actual_residual),best_cost=best[1])
        return cache[key]
    def fun(b):
        counts["fun_calls"] += 1
        return get(b)[1]
    def jac(b):
        counts["jac_calls"] += 1
        p,_ = get(b)
        if p.J is None:
            counts["jacobian_evaluations"] += 1
        J = np.asarray(target.jacobian(p),float)
        return np.vstack((J[:,:64]@Vs,Vs))
    def iteration(intermediate_result):
        _event(callback,"ITERATION",b=intermediate_result.x.copy(),cost=float(intermediate_result.cost))
    result = least_squares(fun,b0,jac=jac,method="trf",tr_solver="exact",
        max_nfev=max_evaluations,x_scale="jac",ftol=1e-8,xtol=1e-8,gtol=1e-6,
        callback=iteration)
    mode,cost,bmode = best
    if mode.J is None:
        counts["jacobian_evaluations"] += 1
    Jb = np.asarray(target.jacobian(mode),float)[:,:64]@Vs
    H,L,fallback = _factor(Jb.T@Jb+np.eye(s))
    C = cho_solve((L,True),np.eye(s))
    logdet = float(2*np.log(np.diag(L)).sum())
    U_mode = cost+32*LOG_2PI
    logZ = -U_mode+.5*s*LOG_2PI-.5*logdet
    gradient = Jb.T@target.white(mode.state.data-y)+Vs.T@mode.q[:64]
    reference = dict(mean_b=bmode,precision_b=H,covariance_b=C,
        precision_cholesky=L,logdet_precision=logdet,logZ_GN=logZ,
        U_mode=U_mode,objective_mode=cost,original_gaussian_log_normalizer=-32*LOG_2PI,
        reference_measure="normalized Gaussian density with respect to db",
        hessian_definition="Jb.T Jb + I; Gauss-Newton, NOT exact U_bb",
        approximation_status="UNTRUNCATED_GN_GAUSSIAN_REFERENCE_NOT_POSTERIOR",
        numerical_fallback=fallback)
    profile = dict(target_family="SLICE_PROFILE",target_subtype="PROFILE",
        strong_direction_treatment="conditional_optimization",log_profile=-U_mode,
        reference_measure="da du; optimized b section, not weak marginal",
        approximation_status="FINITE_OPTIMIZATION_PROFILE_NOT_A_MARGINAL")
    material = conditional_material_moments(target.prior,frame,a,u,reference)
    summary = dict(point=mode,q=mode.q.copy(),a=np.asarray(a).copy(),u=u.copy(),
        reference=reference,profile=profile,material_moments=material,
        target_family=target_family,reference_measure="da db du; original Gaussian64 x uniform2",
        strong_direction_treatment="conditional_GN_reference",
        approximation_status="APPROXIMATE_TARGET_OR_IMPORTANCE_REFERENCE",
        counts=counts,scipy_status=int(result.status),scipy_success=bool(result.success),
        scipy_message=str(result.message),nfev=int(result.nfev),njev=int(result.njev),
        strong_gradient_inf=float(np.max(np.abs(gradient))),
        conditional_mode_target="FULL_CONDITIONAL_ENERGY",
        local_region_contains_mode=None if region is None else bool(region(mode.q.copy())),
        GN_logZ_ignores_local_support=region is not None,
        conditional_mode_converged=bool(np.max(np.abs(gradient)) < 1e-4),
        caution="Finite mode error, GN curvature error and support truncation remain; no posterior-reference qualification.")
    _event(callback,"CONDITIONAL_REFERENCE",q=mode.q.copy(),endpoint=mode.key,
        logZ_GN=logZ,counts=counts,scipy_status=int(result.status),numerical_fallback=fallback)
    return summary


def importance_correct(target,y,frame,a,u,reference,*,samples,rng,region=None,
                       target_family="FULL",callback=None):
    """K exact-target evaluations under Gaussian r(b), returning finite-K SNIS.

    logweights = log[p0(q) exp(-Phi(q)) 1_C] - log r(b). The estimated
    normalizer is nonnegative and unbiased under ordinary reference draws;
    its logarithm and self-normalized material estimates are generally biased.
    This function does not implement a pseudo-marginal Markov chain.
    """
    if (target_family not in ("FULL","LOCAL_C") or
            (target_family == "LOCAL_C" and region is None) or
            (target_family == "FULL" and region is not None)):
        raise ValueError("FULL or explicitly region-defined LOCAL_C target required")
    bs = sample_reference(reference,samples,rng)
    log_r = reference_logpdf(reference,bs)
    logweights, materials, points = [], [], []
    for i,b in enumerate(bs):
        q = compose(frame,a,u,b)
        _event(callback,"EVALUATION_ATTEMPT",q=q,purpose="conditional_importance",sample=i)
        p = target.evaluate(q,purpose="p3l_exact_conditional_importance")
        lw = _loggamma(target,p,y,region)-log_r[i]
        if np.isnan(lw) or np.isposinf(lw):
            raise FloatingPointError("invalid exact importance logweight")
        logweights.append(lw)
        materials.append(target.prior.decode(q)[0])
        points.append(p)
        _event(callback,"IMPORTANCE_SAMPLE",q=q,endpoint=p.key,sample=i,
               logweight=float(lw) if np.isfinite(lw) else None,zero_weight=not np.isfinite(lw))
    lw = np.asarray(logweights)
    positive = bool(np.isfinite(lw).any())
    if positive:
        lsum = float(logsumexp(lw))
        weights = np.exp(lw-lsum)
        images = np.asarray(materials)
        mean = np.einsum("i,i...->...",weights,images)
        variance = np.einsum("i,i...->...",weights,np.abs(images-mean)**2)
        ess = float(1/(weights@weights))
        logZ_estimate = lsum-np.log(samples)
    else:
        weights = np.zeros(samples)
        mean,variance,ess,logZ_estimate = None,None,0.,-np.inf
    result = dict(target_family=target_family,reference_measure="da db du; original Gaussian64 x uniform2",
        strong_direction_treatment="exact_likelihood_importance_correction",
        approximation_status="FINITE_K_SELF_NORMALIZED_IMPORTANCE_BIASED_MOMENT_ESTIMATE",
        status="FINITE_K_UNQUALIFIED" if positive else "NO_POSITIVE_WEIGHTS",
        target_evaluations=samples,jacobian_evaluations=0,points=points,
        strong_samples=bs,log_reference_density=log_r,logweights=lw,weights=weights,
        importance_ess=ess,material_mean=mean,material_variance=variance,
        logZ_estimate=float(logZ_estimate),logZ_estimate_is_unbiased=False,
        normalizer_estimator_is_nonnegative_unbiased=True,
        material_estimate_is_unbiased=False,posterior_samples_returned=False,
        maximum_normalized_weight=float(weights.max()),
        caution="Importance ESS does not prove support coverage; no MCSE or posterior teacher qualification inferred.")
    _event(callback,"IMPORTANCE_SUMMARY",samples=samples,importance_ess=ess,
        logZ_estimate=float(logZ_estimate) if np.isfinite(logZ_estimate) else None,status=result["status"])
    return result
