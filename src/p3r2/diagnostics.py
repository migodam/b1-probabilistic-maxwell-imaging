"""Rank split Rhat, Geyer autocorrelation ESS and material-function MCSE.

Inputs use (chain, draw, *function_shape). Complex material values are split
into real/imaginary components, adding a final axis of size 2. Warmup must be
removed by the caller. These diagnostics cannot certify global exploration.
"""
import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata


def _array(values):
    a = np.asarray(values)
    if np.iscomplexobj(a):
        a = np.stack((a.real, a.imag), axis=-1)
    a = np.asarray(a, dtype=float)
    if a.ndim < 2 or a.shape[0] < 2 or a.shape[1] < 4:
        raise ValueError("require at least two chains and four draws per chain")
    if not np.all(np.isfinite(a)):
        raise ValueError("nonfinite samples: diagnostics are undefined")
    return a


def _split(x):
    n = x.shape[1]//2
    return np.concatenate((x[:,:n], x[:,-n:]), axis=0)


def _rank(x):
    r = rankdata(x, method="average").reshape(x.shape)
    return ndtri((r-3/8)/(x.size+1/4))


def _rhat(x):
    m, n = x.shape
    if np.all(x == x[:, :1]):
        return np.nan if np.all(x == x.flat[0]) else np.inf
    W = np.mean(np.var(x, axis=1, ddof=1))
    B = n*np.var(np.mean(x, axis=1), ddof=1)
    if W == 0:
        return np.nan if B == 0 else np.inf
    return float(np.sqrt(((n-1)*W/n+B/n)/W))


def _ess(x):
    """Multi-chain variance plus initial-positive, initial-monotone pair sums."""
    m, n = x.shape
    if np.all(x == x.flat[0]):
        return np.nan
    if np.all(x == x[:, :1]):
        return 1.0
    centered = x-x.mean(axis=1, keepdims=True)
    fft_size = 1 << (2*n-1).bit_length()
    f = np.fft.rfft(centered, n=fft_size, axis=1)
    # Biased autocovariance / n, as used by Stan/ArviZ's ESS estimator.
    acov = np.fft.irfft(f*np.conjugate(f), n=fft_size, axis=1)[:,:n]/n
    W = np.mean(acov[:,0])*n/(n-1)
    variance = W*(n-1)/n + np.var(x.mean(axis=1), ddof=1)
    if variance <= 0 or W == 0:
        return np.nan if variance <= 0 else 1.0
    rho = 1-(W-acov.mean(axis=0))/variance
    rho[0] = 1.
    pairs = []
    for t in range(0, n-1, 2):
        pair = rho[t]+rho[t+1]
        if pair < 0:
            break
        pairs.append(pair if not pairs else min(pair, pairs[-1]))
    if not pairs:
        return float(m*n)
    tau = max(-1+2*np.sum(pairs), 1/np.log10(m*n))
    # Antithetic chains can legitimately have ESS > number of stored draws.
    return float(min(m*n/tau, m*n*np.log10(m*n)))


def _each(a, function):
    shape = a.shape[2:]
    flat = a.reshape(a.shape[0], a.shape[1], -1)
    result = [function(flat[:,:,k]) for k in range(flat.shape[-1])]
    return np.asarray(result).reshape(shape)


def rank_normalized_split_rhat(values):
    """Maximum of rank-normalized and folded rank-normalized split Rhat.

    Constants across all chains return NaN; chains constant at different
    values return inf. Neither is silently reported as converged.
    """
    a = _array(values)
    def calc(x):
        s = _split(x)
        if np.all(s == s.flat[0]):
            return np.nan
        rank_rhat = _rhat(_rank(s))
        folded_rhat = _rhat(_rank(np.abs(s-np.median(s))))
        return float(np.nanmax([rank_rhat, folded_rhat]))
    return _each(a, calc)


def effective_sample_size(values, method="bulk", probability=.05):
    """Geyer ESS: raw/mean, bulk ranks, or minimum 5%/95% tail indicator ESS."""
    a = _array(values)
    if method not in ("bulk", "tail", "mean", "raw"):
        raise ValueError("method must be bulk, tail, mean or raw")
    if not 0 < probability < .5:
        raise ValueError("tail probability must lie in (0,.5)")
    def calc(x):
        s = _split(x)
        if np.all(s == s.flat[0]):
            return np.nan
        if method == "bulk":
            return _ess(_rank(s))
        if method == "tail":
            low, high = np.quantile(s, [probability, 1-probability])
            lo = _ess((s <= low).astype(float))
            hi = _ess((s <= high).astype(float))
            # Quantile ties can make an indicator constant: it is undefined,
            # not a reason to ignore the failed tail and report the other.
            return min(lo, hi) if np.isfinite(lo) and np.isfinite(hi) else np.nan
        return _ess(s)
    return _each(a, calc)


def material_diagnostics(values):
    """Diagnostics for every material observable, including variance MCSE.

    Variance MCSE uses ESS of (h-pooled_mean)^2. It is an asymptotic plug-in
    estimate, not an exact confidence bound. Globally constant observables
    have MCSE=0 but undefined ESS/Rhat and status CONSTANT_UNASSESSED.
    """
    a = _array(values)
    mean = a.mean(axis=(0,1))
    sd = a.std(axis=(0,1), ddof=1)
    const = np.all(a == a[0,0], axis=(0,1))
    stuck = np.all(a == a[:, :1], axis=(0,1)) & ~const
    raw_ess = effective_sample_size(a, "mean")
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_mcse = sd/np.sqrt(raw_ess)
    centered_square = (a-mean)**2
    variance_ess = effective_sample_size(centered_square, "mean")
    with np.errstate(divide="ignore", invalid="ignore"):
        variance_mcse = centered_square.std(axis=(0,1), ddof=1)/np.sqrt(variance_ess)
    mean_mcse = np.where(const, 0., np.where(stuck, np.nan, mean_mcse))
    variance_mcse = np.where(const, 0., np.where(stuck, np.nan, variance_mcse))
    status = np.where(const, "CONSTANT_UNASSESSED", np.where(stuck,
                      "CHAINS_STUCK_AT_DIFFERENT_VALUES", "DIAGNOSTICS_ONLY"))
    return dict(mean=mean, variance=sd**2, rhat=rank_normalized_split_rhat(a),
        bulk_ess=effective_sample_size(a, "bulk"), tail_ess=effective_sample_size(a, "tail"),
        mean_ess=raw_ess, mean_mcse=mean_mcse, variance_ess=variance_ess,
        variance_mcse=variance_mcse, status=status,
        original_chains=a.shape[0], draws_per_chain=a.shape[1],
        complex_components="real,imag on final axis" if np.iscomplexobj(values) else None,
        caution="No global mixing certificate; MCSE assumes stationary adequately explored chains.")


def crosscheck_chains(first, second, *, independent=False, tolerance=None):
    """Compare two mechanisms' material means, using independent-run MCSE.

    Set independent=True only with independently seeded runs. A tolerance
    must be frozen in material units; absent one, no agreement pass is made.
    This is per-observable, with no multiple-comparison correction.
    """
    a, b = material_diagnostics(first), material_diagnostics(second)
    if a["mean"].shape != b["mean"].shape:
        raise ValueError("material function shapes differ")
    difference = a["mean"]-b["mean"]
    combined = np.sqrt(a["mean_mcse"]**2+b["mean_mcse"]**2) if independent else np.full_like(difference, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        standardized = np.abs(difference)/combined
    valid = ((a["status"] == "DIAGNOSTICS_ONLY") & (b["status"] == "DIAGNOSTICS_ONLY")
             & np.isfinite(combined) & np.isfinite(a["rhat"]) & np.isfinite(b["rhat"]))
    agreement = None
    if tolerance is not None:
        tol = np.broadcast_to(np.asarray(tolerance, float), difference.shape)
        if not np.all(np.isfinite(tol)) or np.any(tol <= 0):
            raise ValueError("positive frozen tolerances required")
        # Conservative compatibility criterion, not equivalence of targets:
        # entire approximate 95% difference interval lies inside tolerance.
        agreement = valid & (np.abs(difference)+1.96*combined < tol)
    return dict(difference=difference, independent_mcse=combined,
        standardized_difference=standardized, valid=valid,
        agreement_within_frozen_tolerance=agreement, first=a, second=b,
        caution="Agreement does not exclude shared bias or missed modes; finite diagnostics are not acceptance gates.")


# Readable alias for callers accustomed to the conventional short name.
split_rhat = rank_normalized_split_rhat
