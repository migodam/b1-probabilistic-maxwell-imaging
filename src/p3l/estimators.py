"""Material decisions from explicit samples; no physics or truth access.

All summaries condition on the samples supplied by the caller. They neither
qualify a posterior reference nor estimate mixing, coverage or escape mass.
Complex variance means E|X-E X|², not E[(X-E X)²].
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


N_PIXELS = 576


def _samples(samples):
    x = np.asarray(samples, dtype=complex)
    if x.ndim != 2 or x.shape[0] == 0 or x.shape[1] != N_PIXELS:
        raise ValueError('samples must have nonempty shape (n,576)')
    if not np.all(np.isfinite(x)):
        raise ValueError('samples must be finite')
    return x


def posterior_mean_variance(samples, *, ddof=0):
    """Unweighted sample mean and pointwise complex modulus variance.

    ddof=0 is the empirical posterior variance; ddof=1 has its usual iid
    sample-variance interpretation only when supplied draws are independent.
    Rejected MCMC states must remain in the input; no deduplication is done.
    """
    x = _samples(samples)
    if not isinstance(ddof, (int, np.integer)) or not 0 <= ddof < len(x):
        raise ValueError('ddof must be an integer in [0,n)')
    mean = np.mean(x, axis=0)
    variance = np.sum(np.abs(x-mean)**2, axis=0)/(len(x)-ddof)
    return dict(mean=mean, variance=variance, sample_count=len(x), ddof=int(ddof),
                covariance_scope='PIXEL_MODULUS_VARIANCE_ONLY')


@dataclass(frozen=True)
class GeometricMedian:
    value: np.ndarray
    converged: np.ndarray
    iterations: np.ndarray
    mean_absolute_distance: np.ndarray
    normalized_subgradient_residual: np.ndarray


def complex_geometric_median(samples, *, tolerance=1e-10, max_iterations=2000):
    """Pointwise minimizer of sum |x_i-t| using modified Weiszfeld.

    A coincident iterate is treated with the multiplicity-aware subgradient
    rule; repeated samples are retained with their original mass. A median
    need not be unique (e.g. two distinct samples). Failed convergence is
    returned per pixel, never silently upgraded to a certified optimum.
    """
    x = _samples(samples)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('positive finite tolerance required')
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError('positive integer max_iterations required')
    values = np.empty(N_PIXELS, complex)
    converged = np.zeros(N_PIXELS, bool)
    iterations = np.zeros(N_PIXELS, int)
    objectives = np.empty(N_PIXELS)
    residuals = np.empty(N_PIXELS)
    for j in range(N_PIXELS):
        points = x[:, j]
        y = complex(np.median(points.real), np.median(points.imag))
        scale = max(1., float(np.max(np.abs(points))))
        coincidence = 8*np.finfo(float).eps*scale
        for k in range(1, max_iterations+1):
            delta = points-y
            distance = np.abs(delta)
            nonzero = distance > coincidence
            multiplicity = len(points)-np.count_nonzero(nonzero)
            if not np.any(nonzero):
                converged[j] = True
                break
            directions = delta[nonzero]/distance[nonzero]
            norm_r = abs(np.sum(directions))
            if norm_r <= multiplicity+tolerance*len(points):
                converged[j] = True
                break
            weights = 1/distance[nonzero]
            # Centered form avoids cancellation for a large common offset.
            step = np.sum(weights*delta[nonzero])/np.sum(weights)
            y_new = y+(1-min(1., multiplicity/norm_r))*step
            if abs(y_new-y) <= tolerance*scale:
                y = y_new
                break
            y = y_new
        delta = points-y
        distance = np.abs(delta)
        nonzero = distance > coincidence
        multiplicity = len(points)-np.count_nonzero(nonzero)
        r = abs(np.sum(delta[nonzero]/distance[nonzero])) if np.any(nonzero) else 0.
        residual = max(0., r-multiplicity)/len(points)
        # This is a numerical stationarity diagnostic, not an interval proof.
        converged[j] = residual <= max(tolerance, 32*np.finfo(float).eps)
        values[j], iterations[j] = y, k
        objectives[j], residuals[j] = np.mean(distance), residual
    return GeometricMedian(values, converged, iterations, objectives, residuals)


def component_quantiles(samples, levels=(.025, .5, .975)):
    """Marginal Re/Im quantiles; not simultaneous or joint complex regions."""
    x = _samples(samples)
    levels = np.asarray(levels, float)
    if (levels.ndim != 1 or len(levels) == 0 or not np.all(np.isfinite(levels))
            or np.any((levels < 0) | (levels > 1))):
        raise ValueError('levels must be a nonempty vector in [0,1]')
    return dict(levels=levels.copy(), real=np.quantile(x.real, levels, axis=0),
                imag=np.quantile(x.imag, levels, axis=0),
                interpretation='COMPONENTWISE_MARGINAL_EMPIRICAL_QUANTILES')


def material_diagnostics(samples):
    """Truth-free empirical spread diagnostics; no risk or coverage claim."""
    x = _samples(samples)
    moments = posterior_mean_variance(x)
    mean, variance = moments['mean'], moments['variance']
    rms = np.sqrt(np.mean(np.abs(x-mean)**2, axis=1))
    return dict(sample_count=len(x), mean_material_variance=float(np.mean(variance)),
                material_spread_rms=float(np.sqrt(np.mean(variance))),
                sample_distance_from_mean_rms=rms,
                sample_mean_real=np.mean(x.real, axis=1),
                sample_mean_imag=np.mean(x.imag, axis=1),
                maximum_pixel_variance=float(np.max(variance)),
                interpretation='EMPIRICAL_SPREAD_NOT_GT_RISK_OR_MCSE')


def _decode(prior, q):
    decoded = prior.decode(q)
    value = decoded[0] if isinstance(decoded, tuple) else decoded
    value = np.asarray(value, complex)
    if value.shape != (N_PIXELS,) or not np.all(np.isfinite(value)):
        raise ValueError('decoder must return a finite complex material (576,)')
    return value


def project_legal_bayes_action(posterior_mean, prior, initial_q, *,
                               posterior_latent_samples=None, maximum_starts=4,
                               max_nfev=1000, tolerance=1e-9):
    """Local decoder projection of an explicit material posterior mean.

    Minimize ||g(q)-mean||²/576 with no prior penalty. The supplied MaterialPrior
    must declare dimension and gaussian_dimension; remaining coordinates have
    bounds [0,1]. Starts are the supplied q and deterministic farthest-point
    selections of supplied latent draws in fixed Euclidean latent distance.
    This selection uses neither decoder fit nor truth. No random restart or
    physics occurs. The returned best local result is not a global argmin.
    """
    mean = np.asarray(posterior_mean, complex)
    if mean.shape != (N_PIXELS,) or not np.all(np.isfinite(mean)):
        raise ValueError('posterior_mean must be finite and have shape (576,)')
    dimension, gaussian = int(prior.dimension), int(prior.gaussian_dimension)
    if not 0 <= gaussian <= dimension:
        raise ValueError('invalid prior coordinate dimensions')
    if (not isinstance(maximum_starts, (int, np.integer)) or maximum_starts < 1
            or not isinstance(max_nfev, (int, np.integer)) or max_nfev < 1
            or not np.isfinite(tolerance) or tolerance <= np.finfo(float).eps):
        raise ValueError('positive integer caps and tolerance above machine epsilon required')
    lower, upper = np.full(dimension, -np.inf), np.full(dimension, np.inf)
    lower[gaussian:], upper[gaussian:] = 0., 1.
    def validate(q):
        q = np.asarray(q, float)
        if (q.shape != (dimension,) or not np.all(np.isfinite(q))
                or np.any(q < lower) or np.any(q > upper)):
            raise ValueError('latent start/draw is nonfinite or outside prior bounds')
        return q.copy()
    starts = [validate(initial_q)]
    if posterior_latent_samples is not None:
        candidates = np.asarray(posterior_latent_samples, float)
        if candidates.ndim != 2 or candidates.shape[1] != dimension:
            raise ValueError('posterior_latent_samples must have shape (n,dimension)')
        candidates = [validate(q) for q in candidates]
        if candidates:
            candidates = np.array(candidates)
            distances = np.sum((candidates-starts[0])**2, axis=1)
            for _ in range(min(maximum_starts-1, len(candidates))):
                index = int(np.argmax(distances))
                if distances[index] == 0:
                    break
                starts.append(candidates[index].copy())
                distances = np.minimum(distances, np.sum((candidates-starts[-1])**2, axis=1))
    def residual(q):
        delta = _decode(prior, q)-mean
        return np.r_[delta.real, delta.imag]/np.sqrt(N_PIXELS)
    def jacobian(q):
        derivative = np.asarray(prior.derivative(q), complex)
        if derivative.shape != (N_PIXELS, dimension) or not np.all(np.isfinite(derivative)):
            raise ValueError('decoder derivative must have shape (576,dimension) and be finite')
        return np.r_[derivative.real, derivative.imag]/np.sqrt(N_PIXELS)
    trials = []
    best = None
    for index, start in enumerate(starts):
        initial_objective = float(np.sum(residual(start)**2))
        fit = least_squares(residual, start, jac=jacobian if callable(getattr(prior, 'derivative', None)) else '2-point',
                            bounds=(lower, upper), max_nfev=max_nfev,
                            ftol=tolerance, xtol=tolerance, gtol=tolerance)
        objective = float(2*fit.cost)
        gradient = 2*np.asarray(fit.grad)
        projected_gradient = fit.x-np.clip(fit.x-gradient, lower, upper)
        trial = dict(start_index=index, initial_q=start, q=fit.x.copy(),
                     initial_objective=initial_objective, objective=objective,
                     gradient=gradient, gradient_inf=float(np.max(np.abs(gradient))),
                     projected_gradient_inf=float(np.max(np.abs(projected_gradient))),
                     scipy_optimality=float(fit.optimality), active_mask=fit.active_mask.copy(),
                     success=bool(fit.success), status=int(fit.status), message=str(fit.message),
                     nfev=int(fit.nfev), njev=None if fit.njev is None else int(fit.njev))
        trials.append(trial)
        # Retain a start if numerical optimization happens to worsen it.
        choice_q = fit.x if objective <= initial_objective else start
        choice_objective = min(objective, initial_objective)
        if best is None or choice_objective < best[0]:
            best = (choice_objective, choice_q.copy(), index, objective <= initial_objective)
    objective, q, index, optimized = best
    return dict(q=q, material=_decode(prior, q), projection_mse=objective,
                projection_rms=float(np.sqrt(objective)), selected_start=index,
                selected_optimized_result=bool(optimized),
                optimization_success=bool(optimized and trials[index]['success']),
                gradient=trials[index]['gradient'].copy() if optimized else None,
                projected_gradient_inf=trials[index]['projected_gradient_inf'] if optimized else None,
                optimization_status=trials[index]['status'] if optimized else None,
                trials=trials,
                start_selection='FIXED_EUCLIDEAN_LATENT_FARTHEST_POINT_NO_GT',
                status='LOCAL_NUMERICAL_PROJECTION_NOT_GLOBAL_ARGMIN',
                objective='FULL_GRID_COMPLEX_SQUARED_DISTANCE_TO_SUPPLIED_MEAN_NO_PRIOR')
