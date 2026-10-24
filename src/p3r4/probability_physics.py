"""Pure diagnostics for paid frozen-sample likelihood evaluations; no solver."""
import numpy as np
from scipy.special import gammainc,logsumexp


def stable_importance(phi,x,logq):
    phi=np.asarray(phi,float);x=np.asarray(x,float);logq=np.asarray(logq,float)
    if phi.ndim!=1 or x.shape!=(len(phi),66) or logq.shape!=phi.shape:raise ValueError('paired sample dimensions')
    if not np.isfinite(phi).all() or not np.isfinite(x).all() or not np.isfinite(logq).all():
        return dict(status='NOT_EVALUATED_INCOMPLETE',planned_samples=len(phi),finite_samples=int(np.isfinite(phi).sum()))
    logw=-phi-.5*np.sum(x*x,axis=1)-logq
    maximum=float(np.max(logw));centered=logw-maximum
    lognorm=logsumexp(centered);normalized_log= centered-lognorm
    weights=np.exp(normalized_log);z=maximum+lognorm
    return dict(status='EVALUATED',log_weights=logw,normalized_weights=weights,
                ESS=float(np.exp(-logsumexp(2*normalized_log))),maximum_weight=float(weights.max()),
                sample_count=len(phi),log_mean_unnormalized_weight=float(z-np.log(len(phi))),
                interpretation='algebraic importance-weight concentration, not independent-chain ESS or FULL qualification; duplicate/rejected states retained')


def likelihood_rank(phi_samples,phi_truth,seed):
    samples=np.asarray(phi_samples,float)
    if not np.isfinite(samples).all() or not np.isfinite(phi_truth):return dict(status='NOT_EVALUATED_INCOMPLETE',planned_samples=len(samples),finite_samples=int(np.isfinite(samples).sum()))
    lower=int(np.count_nonzero(samples<phi_truth));ties=int(np.count_nonzero(samples==phi_truth));u=float(np.random.default_rng(seed).uniform())
    return dict(status='EVALUATED',truth_phi=float(phi_truth),samples=len(samples),lower_count=lower,tie_count=ties,
                randomized_inclusive_rank=float((lower+u*(ties+1))/(len(samples)+1)),
                statistic='Phi(theta;y), data-dependent statistic; marginal rank diagnostic, not FULL posterior certificate')


def gamma256_summary(phi):
    a=np.asarray(phi,float);finite=np.isfinite(a);pits=gammainc(256,a[finite])
    return dict(status='EVALUATED' if finite.all() else 'INCOMPLETE_ALL_FAILURES_RETAINED',planned_parents=len(a),finite_parents=int(finite.sum()),failed_parents=int((~finite).sum()),
                mean_all_parents=float(a.mean()) if finite.all() and len(a) else None,
                finite_only_mean_descriptive=float(a[finite].mean()) if finite.any() else None,
                Gamma_shape=256,Gamma_scale=1,fixed_sample_index=0,
                PIT_finite_mean=float(pits.mean()) if len(pits) else None,
                PIT_below_025=int(np.sum(pits<.025)),PIT_above_975=int(np.sum(pits>.975)),
                interpretation='prior-predictive fixed-draw residual benchmark; incomplete parents are not removed from planned denominator; no FULL qualification')
