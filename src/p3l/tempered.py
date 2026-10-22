"""Generic adaptive tempered SMC; no data, solver, Book or training imports.

The caller supplies normalized-prior draws and log densities in one fixed
coordinate system. Original Gaussian x uniform coordinates and exact probit
coordinates are both supported without changing any prior inside this module.
All function invocation counts are exposed; parent wrappers own real costs.
"""
import numpy as np
from scipy.linalg import cho_solve, solve_triangular
from scipy.special import logsumexp


def _logvalue(value,kind):
    value = float(value)
    if np.isnan(value) or np.isposinf(value):
        raise ValueError(f"{kind} returned NaN or +inf")
    return value


def _multiply(delta,values):
    # beta=0 must not form 0*(-inf) for valid zero-likelihood particles.
    return np.zeros_like(values) if delta == 0 else delta*values


def conditional_ess(logweights,loglikelihood,delta_beta):
    """N (E_w[G])²/E_w[G²], G=exp(delta_beta*loglikelihood).

    logweights must be normalized. This conditional ESS is distinct from
    ordinary importance ESS; it starts at N even with unequal incoming weights.
    """
    lw,ll = np.asarray(logweights,float),np.asarray(loglikelihood,float)
    if lw.ndim != 1 or ll.shape != lw.shape or not np.isfinite(logsumexp(lw)):
        raise ValueError("matching nonempty log-weight/likelihood vectors required")
    if np.isnan(ll).any() or np.isposinf(ll).any() or not np.isfinite(delta_beta) or delta_beta < 0:
        raise ValueError("invalid likelihoods or negative temperature increment")
    lw = lw-logsumexp(lw)
    if delta_beta == 0:
        return float(lw.size)
    active = np.isfinite(lw)&np.isfinite(ll)
    if not active.any():
        return 0.
    # A common likelihood offset cancels from CESS. Remove it before doubling
    # log increments, including very large constant negative log likelihoods.
    inc = np.full_like(ll,-np.inf)
    with np.errstate(over="ignore"):
        inc[active] = delta_beta*(ll[active]-np.max(ll[active]))
    first,second = logsumexp(lw+inc),logsumexp(lw+2*inc)
    if not np.isfinite(first) or not np.isfinite(second):
        return 0.
    return float(np.clip(np.exp(np.log(lw.size)+2*first-second),0.,lw.size))


def choose_beta(beta,logweights,loglikelihood,*,target_fraction=.8,
                min_delta_beta=1e-8,bisection_steps=64):
    """Largest next temperature meeting conditional ESS, or explicit failure.

    A too-small safe increment fails instead of forcing a floor that violates
    the ESS contract. A final increment to beta=1 may be smaller than the floor.
    """
    if not 0 <= beta < 1 or not 0 < target_fraction < 1 or not np.isfinite(min_delta_beta) or min_delta_beta <= 0:
        raise ValueError("invalid beta or conditional-ESS configuration")
    target = target_fraction*np.asarray(logweights).size
    remaining = 1.-beta
    if conditional_ess(logweights,loglikelihood,remaining) >= target:
        return dict(status="OK",beta=1.,delta_beta=remaining,
                    conditional_ess=conditional_ess(logweights,loglikelihood,remaining))
    lo,hi = 0.,remaining
    for _ in range(bisection_steps):
        mid = .5*(lo+hi)
        if conditional_ess(logweights,loglikelihood,mid) >= target:
            lo = mid
        else:
            hi = mid
    if lo < min_delta_beta or beta+lo <= beta:
        return dict(status="TEMPERING_STALLED_BELOW_BETA_FLOOR",beta=beta,
                    delta_beta=0.,safe_increment=lo,conditional_ess=None)
    return dict(status="OK",beta=beta+lo,delta_beta=lo,
                conditional_ess=conditional_ess(logweights,loglikelihood,lo))


def systematic_resample(weights,rng):
    """Return parent particle indices, with exactly one random offset."""
    w = np.asarray(weights,float)
    if w.ndim != 1 or not w.size or not np.all(np.isfinite(w)) or np.any(w < 0) or w.sum() <= 0:
        raise ValueError("finite nonnegative weights of positive sum required")
    cumulative = np.cumsum(w/w.sum())
    cumulative[-1] = 1.
    positions = (np.arange(w.size)+rng.random())/w.size
    return np.searchsorted(cumulative,positions,side="right")


def tempered_smc(sampler_prior,logprior,loglikelihood,*,rng,N=128,mutation_steps=2,
                 mutation="random_walk",proposal_scale=.2,step_size=.1,
                 gaussian_dimension=None,gradloglikelihood=None,gradlogprior=None,
                 precision=None,cess_target=.8,resample_threshold=.5,
                 min_delta_beta=1e-8,max_stages=256,callback=None):
    """Adaptive conditional-ESS tempering with exact MH mutation correction.

    sampler_prior(rng,N)->(N,d); logprior/loglikelihood(x)->scalar.
    Optional gradients are with respect to exactly those d coordinates.
    precision(x,beta)->SPD(d,d) is a deterministic frozen function of state
    and temperature. Its inverse is the MALA covariance preconditioner.

    Mutations: symmetric random_walk; Gaussian-reference pcn; or normalized
    two-endpoint mala. gaussian_dimension explicitly declares N(0,I) in the
    first g coordinates and (when g=d-2) independent Uniform[0,1]^2 properties.
    With g=d this includes exact probit Gaussian coordinates. Without this
    declaration, MALA needs an explicit gradlogprior. Native-box MALA proposals
    are rejected outside support, never reflected or clipped. random_walk/pcn
    may reflect the explicitly declared uniform property block symmetrically.

    The output always includes its actual terminal beta and failure status.
    Adaptive finite-N schedules, evidence and moments have finite-particle
    error/bias; weight ESS is never a global exploration certificate.
    """
    if not isinstance(N,int) or N < 2 or not isinstance(mutation_steps,int) or mutation_steps < 0:
        raise ValueError("N >= 2 and nonnegative integer mutation_steps required")
    if mutation not in ("random_walk","pcn","mala"):
        raise ValueError("mutation must be random_walk, pcn or mala")
    if not np.isfinite(proposal_scale) or proposal_scale <= 0 or not np.isfinite(step_size) or step_size <= 0:
        raise ValueError("positive finite proposal scales required")
    if not 0 < cess_target < 1 or not 0 <= resample_threshold <= 1:
        raise ValueError("invalid ESS or resampling threshold")
    if not np.isfinite(min_delta_beta) or min_delta_beta <= 0 or not isinstance(max_stages,int) or max_stages < 1:
        raise ValueError("positive beta floor and integer maximum stages required")
    x = np.array(sampler_prior(rng,N),dtype=float,copy=True)
    if x.ndim != 2 or x.shape[0] != N or x.shape[1] < 1 or not np.all(np.isfinite(x)):
        raise ValueError("prior sampler must return finite (N,d) states")
    d = x.shape[1]
    g = gaussian_dimension
    if g is not None and (not isinstance(g,int) or g < 1 or g not in (d,d-2)):
        raise ValueError("Gaussian reference dimension must equal d or d-2")
    if g == d-2 and (np.any(x[:,g:] < 0) or np.any(x[:,g:] > 1)):
        raise ValueError("declared uniform prior samples violate [0,1]^2")
    if mutation == "pcn" and (g is None or proposal_scale > 1):
        raise ValueError("pcn needs a declared Gaussian reference and scale <= 1")
    if mutation == "mala" and (gradloglikelihood is None or (gradlogprior is None and g is None)):
        raise ValueError("MALA requires likelihood and prior gradients (or declared Gaussian product prior)")
    counters = dict(prior_draws=N,logprior_evaluations=0,likelihood_evaluations=0,
                    likelihood_gradient_evaluations=0,prior_gradient_evaluations=0,precision_evaluations=0)
    configuration = dict(N=N,dimension=d,mutation_steps=mutation_steps,mutation=mutation,
        proposal_scale=proposal_scale,step_size=step_size,gaussian_dimension=g,
        cess_target=cess_target,resample_threshold=resample_threshold,
        min_delta_beta=min_delta_beta,max_stages=max_stages)
    beta,log_evidence = 0.,0.
    history = []
    roots = np.arange(N)
    logw = np.full(N,-np.log(N))
    lp,ll = np.empty(N),np.empty(N)
    stage = 0

    def emit(event,**fields):
        if callback is not None:
            callback(dict(event=event,stage=stage,beta=beta,counters=dict(counters),**fields))

    def prior(q,particle):
        counters["logprior_evaluations"] += 1
        value = _logvalue(logprior(q.copy()),"logprior")
        return value

    def likelihood(q,particle):
        counters["likelihood_evaluations"] += 1
        emit("LIKELIHOOD_EVALUATION_ATTEMPT",particle=particle,x=q.copy())
        value = _logvalue(loglikelihood(q.copy()),"loglikelihood")
        emit("LIKELIHOOD_EVALUATED",particle=particle,
             loglikelihood=value if np.isfinite(value) else None,zero_likelihood=not np.isfinite(value))
        return value

    def gradient_likelihood(q,particle):
        counters["likelihood_gradient_evaluations"] += 1
        emit("LIKELIHOOD_GRADIENT_ATTEMPT",particle=particle,x=q.copy())
        value = np.asarray(gradloglikelihood(q.copy()),float)
        if value.shape != (d,) or not np.all(np.isfinite(value)):
            raise ValueError("likelihood gradient must be finite shape (d,)")
        return value

    def gradient_prior(q):
        if gradlogprior is None:
            return np.r_[-q[:g],np.zeros(d-g)]
        counters["prior_gradient_evaluations"] += 1
        value = np.asarray(gradlogprior(q.copy()),float)
        if value.shape != (d,) or not np.all(np.isfinite(value)):
            raise ValueError("prior gradient must be finite shape (d,)")
        return value

    def metric(q,particle):
        if precision is None:
            return np.eye(d),0.
        counters["precision_evaluations"] += 1
        emit("PRECISION_EVALUATION_ATTEMPT",particle=particle,x=q.copy())
        H = np.array(precision(q.copy(),beta),dtype=float,copy=True)
        if H.shape != (d,d) or not np.all(np.isfinite(H)) or not np.allclose(H,H.T,atol=1e-12,rtol=1e-12):
            raise ValueError("precision must be finite symmetric (d,d)")
        L = np.linalg.cholesky(.5*(H+H.T))
        return L,float(2*np.log(np.diag(L)).sum())

    def reflect_property(q):
        if g == d-2:
            q[g:] = 1-np.abs(1-np.remainder(q[g:],2.))
        return q

    def mutate():
        accepted,attempted = 0,0
        for i in range(N):
            current = x[i].copy()
            current_lp,current_ll = lp[i],ll[i]
            # Local caches retain only one particle's gradient/precision, never
            # a giant per-particle physical Point or factorization store.
            current_gl,current_L,current_logdet = None,None,None
            for _ in range(mutation_steps):
                attempted += 1
                proposal = current.copy()
                correction = 0.
                if mutation == "pcn":
                    proposal[:g] = np.sqrt(1-proposal_scale**2)*current[:g]+proposal_scale*rng.normal(size=g)
                    if g < d:
                        proposal[g:] += proposal_scale*rng.normal(size=d-g)
                    reflect_property(proposal)
                    correction = .5*(proposal[:g]@proposal[:g]-current[:g]@current[:g])
                elif mutation == "random_walk":
                    proposal += proposal_scale*rng.normal(size=d)
                    reflect_property(proposal)
                else:
                    if not np.isfinite(current_ll):
                        # Zero-target particles must not use undefined log
                        # gradients. A symmetric RW can escape at beta>0.
                        proposal += proposal_scale*rng.normal(size=d)
                        reflect_property(proposal)
                    else:
                        if current_gl is None:
                            current_gl = gradient_likelihood(current,i)
                        if current_L is None:
                            current_L,current_logdet = metric(current,i)
                        gradient = gradient_prior(current)+beta*current_gl
                        drift = .5*step_size*cho_solve((current_L,True),gradient)
                        delta = drift+np.sqrt(step_size)*solve_triangular(current_L.T,rng.normal(size=d),lower=False)
                        proposal += delta
                proposed_lp = prior(proposal,i)
                if np.isneginf(proposed_lp):
                    emit("MUTATION_REJECTED_PRIOR_SUPPORT",particle=i)
                    continue
                proposed_ll = likelihood(proposal,i)
                proposed_gl,proposed_L,proposed_logdet = None,None,None
                if mutation == "mala" and np.isfinite(current_ll) and np.isfinite(proposed_ll):
                    proposed_gl = gradient_likelihood(proposal,i)
                    proposed_L,proposed_logdet = metric(proposal,i)
                    reverse_gradient = gradient_prior(proposal)+beta*proposed_gl
                    reverse_drift = .5*step_size*cho_solve((proposed_L,True),reverse_gradient)
                    forward_white = current_L.T@(delta-drift)
                    reverse_white = proposed_L.T@(-delta-reverse_drift)
                    # Normalization constants -d/2 log(2*pi*h) cancel, but
                    # both varying determinants and quadratic terms remain.
                    correction = (.5*(proposed_logdet-current_logdet)
                        +(forward_white@forward_white-reverse_white@reverse_white)/(2*step_size))
                current_logtarget = current_lp+float(_multiply(beta,np.array(current_ll)))
                proposed_logtarget = proposed_lp+float(_multiply(beta,np.array(proposed_ll)))
                if np.isneginf(current_logtarget):
                    ratio = np.inf if np.isfinite(proposed_logtarget) else -np.inf
                else:
                    ratio = proposed_logtarget-current_logtarget+correction
                yes = bool(np.log(max(rng.random(),np.finfo(float).tiny)) < min(0.,ratio))
                if yes:
                    current,current_lp,current_ll = proposal,proposed_lp,proposed_ll
                    current_gl,current_L,current_logdet = proposed_gl,proposed_L,proposed_logdet
                    accepted += 1
                emit("MUTATION",particle=i,accepted=yes,
                     log_mh_ratio=float(ratio) if np.isfinite(ratio) else None,
                     zero_target_rejection=bool(np.isneginf(ratio)))
            x[i],lp[i],ll[i] = current,current_lp,current_ll
        return dict(accepted=accepted,attempted=attempted,
                    acceptance=accepted/attempted if attempted else None)

    for i in range(N):
        lp[i] = prior(x[i],i)
        if not np.isfinite(lp[i]):
            raise ValueError("prior sampler returned a point outside supplied prior support")
        ll[i] = likelihood(x[i],i)
    status = "RUNNING"
    emit("INITIALIZED",particles=x.copy(),ancestor_ids=roots.copy(),
         logweights=logw.copy(),loglikelihood=ll.copy(),logprior=lp.copy(),configuration=configuration)
    while beta < 1. and stage < max_stages:
        stage += 1
        selection = choose_beta(beta,logw,ll,target_fraction=cess_target,min_delta_beta=min_delta_beta)
        if selection["status"] != "OK":
            status = selection["status"]
            emit("TEMPERING_FAILED",selection=selection)
            break
        previous_beta = beta
        increment = _multiply(selection["delta_beta"],ll)
        log_increment = float(logsumexp(logw+increment))
        if not np.isfinite(log_increment):
            status = "NO_POSITIVE_INCREMENTAL_WEIGHTS"
            emit("TEMPERING_FAILED",selection=selection)
            break
        logw = logw+increment-log_increment
        log_evidence += log_increment
        beta = selection["beta"]
        weights = np.exp(logw)
        ess_before = float(1/(weights@weights))
        resampled = ess_before < resample_threshold*N
        ancestors = systematic_resample(weights,rng) if resampled else np.arange(N)
        if resampled:
            x,lp,ll = x[ancestors].copy(),lp[ancestors].copy(),ll[ancestors].copy()
            roots = roots[ancestors]
            logw = np.full(N,-np.log(N))
        before_mutation = dict(counters)
        mutation_record = mutate()
        record = dict(stage=stage,beta_before=previous_beta,beta=beta,
            conditional_ess=selection["conditional_ess"],weight_ess_before_resampling=ess_before,
            resampled=resampled,parent_indices=ancestors.copy(),ancestor_ids=roots.copy(),
            distinct_prior_ancestors=int(np.unique(roots).size),log_evidence_increment=log_increment,
            log_evidence=log_evidence,mutation=mutation_record,
            mutation_call_counts={k:counters[k]-before_mutation[k] for k in counters},
            cumulative_call_counts=dict(counters))
        history.append(record)
        emit("STAGE_COMPLETE",history=record,particles=x.copy(),logweights=logw.copy(),
             loglikelihood=ll.copy(),logprior=lp.copy())
    if beta == 1.:
        status = "COMPLETED_FINITE_N_UNQUALIFIED"
    elif status == "RUNNING":
        status = "MAXIMUM_TEMPERING_STAGES_REACHED"
    weights = np.exp(logw)
    result = dict(status=status,beta=beta,particles=x,logprior=lp,loglikelihood=ll,
        weights=weights,logweights=logw,log_evidence=log_evidence,
        ancestor_ids=roots,history=history,counters=counters,
        configuration=configuration,completed_stages=len(history),attempted_stages=stage,
        weight_ess=float(1/(weights@weights)),mutation=mutation,
        target_family="CALLER_DECLARED_PRIOR_TIMES_LIKELIHOOD_POWER",
        reference_measure="caller coordinate measure; no implicit transform",
        approximation_status="ADAPTIVE_FINITE_N_SMC_NOT_EXACT_POSTERIOR_DRAWS",
        normalized_prior_sampling_assumed=True,log_likelihood_constants_as_supplied=True,
        adaptive_evidence_unbiasedness_claim=False,posterior_qualification=False,
        caution="Finite-N adaptive tempering may be biased and miss modes; weight ESS is not global mixing evidence.")
    emit("FINISHED",status=status,terminal_beta=beta,weight_ess=result["weight_ess"],log_evidence=log_evidence)
    return result
