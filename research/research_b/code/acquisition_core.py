"""Research B: finite-bank acquisition reference core (no Maxwell solver).

All probabilities are conditional on the supplied, explicitly declared bank.
A particle must include persistent nuisance, or the caller must provide the
correct conditional predictive kernel. Atomic banks do not certify continuous
posterior accuracy or unknown-mode coverage. No external services are used.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np
from scipy.special import logsumexp, ndtr
from scipy.stats import truncnorm


def probability_vector(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or len(p) == 0 or not np.isfinite(p).all() or (p < 0).any():
        raise ValueError('Expected a finite, nonnegative, nonempty vector.')
    if not np.isclose(p.sum(), 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError('Weights must sum to one; do not silently discard remainder mass.')
    return p


def entropy(p: np.ndarray) -> float:
    p = probability_vector(p)
    x = p[p > 0]
    return float(-np.dot(x, np.log(x)))


def bayes_risk(weights: np.ndarray, loss: np.ndarray) -> float:
    """loss[d, i]: loss of decision d when particle i is the material world."""
    w = probability_vector(weights)
    loss = np.asarray(loss, dtype=float)
    if loss.ndim != 2 or loss.shape[1] != len(w) or not np.isfinite(loss).all():
        raise ValueError('loss must be a finite decision-by-particle matrix.')
    if (loss < 0).any() or loss.shape[0] == 0:
        raise ValueError('Losses must be nonnegative and decisions nonempty.')
    return float(np.min(loss @ w))


def posterior_update(weights: np.ndarray, log_likelihood_increment: np.ndarray) -> np.ndarray:
    """Increment MUST condition on old data and persistent nuisance.

    With nuisance integrated out, use log p(old,new|material)-log p(old|material).
    Do not integrate a new independent gain prior for each observation block.
    """
    w = probability_vector(weights)
    ell = np.asarray(log_likelihood_increment, dtype=float)
    if ell.shape != w.shape or np.isnan(ell).any() or np.isposinf(ell).any():
        raise ValueError('Invalid log-likelihood increment.')
    lw = np.full_like(w, -np.inf)
    positive = w > 0
    lw[positive] = np.log(w[positive]) + ell[positive]
    z = logsumexp(lw)
    if not np.isfinite(z):
        raise ValueError('Zero predictive support: expand/re-evaluate support, do not reset weights.')
    return np.exp(lw - z)


def event_masses(weights: np.ndarray, labels: np.ndarray, n_events: int) -> np.ndarray:
    w = probability_vector(weights)
    e = np.asarray(labels)
    if e.shape != w.shape or not np.issubdtype(e.dtype, np.integer):
        raise ValueError('Integer event labels required for every particle.')
    if n_events < 1 or (e < 0).any() or (e >= n_events).any():
        raise ValueError('Labels must index the complete modeled event partition.')
    return np.bincount(e, weights=w, minlength=n_events)


def discrete_utility(weights: np.ndarray, channel: np.ndarray, loss: np.ndarray,
                     labels: np.ndarray | None = None, n_events: int | None = None) -> dict:
    """Exact finite observation experiment: channel[i,z] = p(z|world i)."""
    w = probability_vector(weights)
    c = np.asarray(channel, float)
    if c.ndim != 2 or c.shape[0] != len(w) or not np.isfinite(c).all() or (c < 0).any():
        raise ValueError('Invalid channel.')
    if not np.allclose(c.sum(axis=1), 1.0):
        raise ValueError('Every channel row must sum to one.')
    joint = w[:, None] * c
    pred = joint.sum(axis=0)
    good = pred > 0
    post = np.zeros_like(joint)
    post[:, good] = joint[:, good] / pred[None, good]
    future = float(sum(pred[z] * bayes_risk(post[:, z], loss) for z in np.where(good)[0]))
    r0 = bayes_risk(w, loss)
    out = {'current_risk': r0, 'expected_future_risk': future,
           'risk_reduction': r0 - future, 'predictive': pred,
           'posterior': post, 'posterior_mean_weights': post @ pred}
    if labels is not None:
        if n_events is None:
            raise ValueError('n_events is required with labels.')
        he = entropy(event_masses(w, labels, n_events))
        hfuture = sum(pred[z] * entropy(event_masses(post[:, z], labels, n_events))
                      for z in np.where(good)[0])
        out['event_eig_nats'] = he - hfuture
    return out


@dataclass(frozen=True)
class GaussianUtility:
    current_risk: float
    expected_future_risk: float
    risk_reduction: float
    mc_standard_error: float
    hoeffding_halfwidth: float
    event_eig_nats: float
    event_posterior_quantiles: list[list[float]]
    modeled_bank_only: bool = True


def gaussian_bank_utility(weights: np.ndarray, means: np.ndarray, covariance: np.ndarray,
                          loss: np.ndarray, labels: np.ndarray, n_events: int,
                          uniforms: np.ndarray, normal_draws: np.ndarray,
                          failure_probability: float = 0.01,
                          simultaneous_actions: int = 1) -> GaussianUtility:
    """Nested exact-bank inner sums; iid outer predictive samples.

    Use the SAME uniforms/normal_draws for candidate actions (common random
    numbers). Per-action dimensions must agree, or caller must construct a
    registered coupling. Bounds cover outer MC only, not bank/proxy error.
    The decision set supplied in loss is fixed. A medoid-only set evaluates
    constrained material selection risk, not unrestricted field MAE risk.
    """
    w = probability_vector(weights)
    mu = np.asarray(means, float)
    cov = np.asarray(covariance, float)
    u = np.asarray(uniforms, float)
    z = np.asarray(normal_draws, float)
    ell = np.asarray(loss, float)
    if mu.ndim != 2 or mu.shape[0] != len(w) or not np.isfinite(mu).all():
        raise ValueError('Invalid particle means.')
    if cov.shape != (mu.shape[1], mu.shape[1]) or not np.allclose(cov, cov.T):
        raise ValueError('Covariance must be symmetric and match the channels.')
    chol = np.linalg.cholesky(cov)
    if u.ndim != 1 or len(u) < 2 or (u < 0).any() or (u >= 1).any():
        raise ValueError('At least two uniforms in [0,1) are required.')
    if z.shape != (len(u), mu.shape[1]) or not np.isfinite(z).all():
        raise ValueError('Invalid normal_draws shape.')
    if not 0 < failure_probability < 1 or simultaneous_actions < 1:
        raise ValueError('Invalid confidence parameters.')
    r0 = bayes_risk(w, ell)
    e0 = event_masses(w, labels, n_events)
    cdf = np.cumsum(w)
    cdf[-1] = 1.0
    idx = np.searchsorted(cdf, u, side='right')
    observations = mu[idx] + z @ chol.T
    white_mu = np.linalg.solve(chol, mu.T).T
    white_obs = np.linalg.solve(chol, observations.T).T
    risks = np.empty(len(u))
    entropies = np.empty(len(u))
    masses = np.empty((len(u), n_events))
    for k, obs in enumerate(white_obs):
        logl = -0.5 * np.sum((white_mu - obs[None, :]) ** 2, axis=1)
        post = posterior_update(w, logl)
        risks[k] = bayes_risk(post, ell)
        masses[k] = event_masses(post, labels, n_events)
        entropies[k] = entropy(masses[k])
    lmax = float(np.max(ell))
    half = lmax * np.sqrt(np.log(2 * simultaneous_actions / failure_probability) / (2 * len(u)))
    return GaussianUtility(r0, float(risks.mean()), float(r0 - risks.mean()),
                           float(risks.std(ddof=1) / np.sqrt(len(u))), float(half),
                           float(entropy(e0) - entropies.mean()),
                           np.quantile(masses, [0.05, 0.5, 0.95], axis=0).tolist())


def effective_information(material_jacobian: np.ndarray, nuisance_jacobian: np.ndarray,
                          nuisance_precision: np.ndarray | None = None) -> np.ndarray:
    """Already jointly whitened real Jacobians; stack blocks BEFORE this call.

    Proper prior precision enters ONCE. This is a local linear-Gaussian
    diagnostic, not marginal event evidence or a nonlinear branch certificate.
    """
    b, n = np.asarray(material_jacobian, float), np.asarray(nuisance_jacobian, float)
    if b.ndim != 2 or n.ndim != 2 or b.shape[0] != n.shape[0]:
        raise ValueError('Incompatible Jacobians.')
    h = n.T @ n
    if nuisance_precision is None:
        inv = np.linalg.pinv(h, hermitian=True)
    else:
        p = np.asarray(nuisance_precision, float)
        if p.shape != h.shape or not np.allclose(p, p.T) or np.linalg.eigvalsh(p).min() < -1e-12:
            raise ValueError('Prior precision must be symmetric PSD.')
        inv = np.linalg.pinv(h + p, hermitian=True)
    result = b.T @ b - b.T @ n @ inv @ n.T @ b
    return (result + result.T) / 2


def gaussian_equal_covariance_tv(white_mean_distance: np.ndarray | float) -> np.ndarray:
    d = np.asarray(white_mean_distance, float)
    if (d < 0).any() or not np.isfinite(d).all():
        raise ValueError('Distances must be finite and nonnegative.')
    return 2 * ndtr(d / 2) - 1


def shared_gain_loglik(y: np.ndarray, field: np.ndarray, sigma: float,
                       prior_mean: float = 1.0, prior_sd: float = 0.15,
                       lower: float = 0.7, upper: float = 1.3,
                       order: int = 96) -> float:
    """Proper truncated-real-gain quadrature; all old+new channels share ONE gain.

    Real channel variance sigma**2. For the conditional increment subtract the
    old-data result from the concatenated-data result. Quadrature order must be
    converged for the actual likelihood; default order is not a certificate.
    """
    y, f = np.asarray(y, float).ravel(), np.asarray(field, float).ravel()
    if y.shape != f.shape or sigma <= 0 or prior_sd <= 0 or lower >= upper or order < 2:
        raise ValueError('Invalid gain likelihood inputs.')
    if not np.isfinite(y).all() or not np.isfinite(f).all():
        raise ValueError('Non-finite data.')
    nodes, qw = np.polynomial.legendre.leggauss(order)
    gain = lower + (nodes + 1) * (upper - lower) / 2
    logprior = truncnorm.logpdf(gain, (lower-prior_mean)/prior_sd,
                               (upper-prior_mean)/prior_sd, loc=prior_mean, scale=prior_sd)
    residual = (y[None, :] - gain[:, None] * f[None, :]) / sigma
    loglike = -0.5 * np.sum(residual**2, axis=1) - len(y) * np.log(sigma*np.sqrt(2*np.pi))
    return float(logsumexp(np.log(qw*(upper-lower)/2) + logprior + loglike))


def certified_interval_winner(estimates: Sequence[float], error_radii: Sequence[float]) -> int | None:
    """Return an action only when its lower bound exceeds all competitors' upper bounds.

    'Certified' is conditional on caller-provided simultaneous valid bounds.
    Empirical holdout maxima or bootstrap alone do not establish such bounds.
    """
    u, b = np.asarray(estimates, float), np.asarray(error_radii, float)
    if u.ndim != 1 or b.shape != u.shape or len(u) == 0 or (b < 0).any():
        raise ValueError('Invalid intervals.')
    if not np.isfinite(u).all() or not np.isfinite(b).all():
        raise ValueError('Finite intervals required.')
    top = int(np.argmax(u-b))
    others = np.arange(len(u)) != top
    return top if not others.any() or u[top]-b[top] > np.max((u+b)[others]) else None
