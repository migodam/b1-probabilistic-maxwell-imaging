"""Finite-host FFT-TSOM and full-wave projected Gauss–Newton."""

from dataclasses import dataclass
from functools import lru_cache
from time import perf_counter

import numpy as np
from scipy.fft import ifft2
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve, svd
from scipy.optimize import minimize

from .ndt import BOUNDS, Host, host


def energy(a):
    a = np.asarray(a)
    return float(np.vdot(a, a).real)


@lru_cache(maxsize=None)
def _adjacent(n):
    cells = np.arange(n * n).reshape(n, n)
    first = np.r_[cells[:, :-1].ravel(), cells[:-1, :].ravel()]
    second = np.r_[cells[:, 1:].ravel(), cells[1:, :].ravel()]
    return first, second


@dataclass(frozen=True)
class MaterialPenalty:
    """Smoothed anisotropic TV and a smoothed L1 penalty about the intact plate."""

    n: int = 24
    tv: float = 0.022
    amplitude: float = 0.014
    smooth: float = 0.04

    def evaluate(self, x, hessian=None):
        a, b = _adjacent(self.n)
        differences = x[a] - x[b]
        size_x = np.hypot(x, self.smooth)
        size_edge = np.hypot(differences, self.smooth)
        wx = self.amplitude / x.size
        we = self.tv / differences.size
        value = wx * (size_x - self.smooth).sum() + we * (size_edge - self.smooth).sum()
        gradient = wx * x / size_x
        edge_grad = we * differences / size_edge
        np.add.at(gradient, a, edge_grad)
        np.add.at(gradient, b, -edge_grad)
        if hessian is not None:
            # Lagged diffusivity is a positive majorizer of the exact TV curvature.
            diag = np.diag_indices(x.size)
            weight = we / size_edge
            hessian[diag] += wx / size_x
            np.add.at(hessian, (a, a), weight)
            np.add.at(hessian, (b, b), weight)
            hessian[a, b] -= weight
            hessian[b, a] -= weight
        return float(value), gradient


@lru_cache(maxsize=64)
def acquisition_models(n, design, *, debye=False):
    if debye:
        from .ndt import host_debye
        host_function = host_debye
    else:
        host_function = host
    selected = []
    for frequency in design.freqs:
        h = host_function(n, frequency)
        selected.append(Host(h.Eh[:, design.tx], h.Sh[design.rx, :], h.Gh, h.xy))
    return tuple(selected)


def fullwave(models, x, *, jacobian=False, material_factors=None):
    """Exact material forward map and analytic distorted-wave Jacobian."""
    data, derivatives = [], []
    x = np.asarray(x).reshape(-1)
    if material_factors is None:
        material_factors = np.ones(len(models))
    for h, multiplier in zip(models, material_factors, strict=True):
        if np.any(x):
            deviation = multiplier * x[:, None]
            state = -deviation * h.Gh
            state.flat[::x.size + 1] += 1
            factor = lu_factor(state, overwrite_a=True, check_finite=False)
            current = lu_solve(factor, deviation * h.Eh, check_finite=False)
        else:
            current = np.zeros_like(h.Eh)
        data.append(h.Sh @ current)
        if jacobian:
            if np.any(x):
                field = h.Eh + h.Gh @ current
                receiver = lu_solve(factor, h.Sh.T, trans=1, check_finite=False).T
            else:
                field, receiver = h.Eh, h.Sh
            jac = (receiver[:, None, :] * field.T[None, :, :]).reshape(-1, x.size)
            derivatives.append(multiplier * jac)
    if jacobian:
        return np.asarray(data), np.concatenate(derivatives)
    return np.asarray(data)


def _bound_step(H, gradient, x):
    """Solve the regularized GN quadratic on the box by a scaled bounded minimizer."""
    lower, upper = BOUNDS
    scale = np.sqrt(np.diag(H))
    free = ~(((x <= lower + 1e-10) & (gradient > 0)) |
             ((x >= upper - 1e-10) & (gradient < 0)))
    newton = np.zeros_like(x)
    factor = cho_factor(H[np.ix_(free, free)], check_finite=False)
    newton[free] = cho_solve(factor, -gradient[free], check_finite=False)
    start = (np.clip(x + newton, lower, upper) - x) * scale
    bounds = list(zip((lower - x) * scale, (upper - x) * scale))

    def quadratic(z):
        step = z / scale
        curvature = H @ step
        return float(step @ (gradient + 0.5 * curvature)), (gradient + curvature) / scale

    solution = minimize(
        quadratic, start, jac=True, bounds=bounds, method="L-BFGS-B",
        options={"maxiter": 45, "ftol": 1e-10, "gtol": 1e-7, "maxcor": 12},
    )
    return solution.x / scale


def projected_gn(models, y, reference, penalty, *, iterations=18, material_factors=None):
    """Recompute the full wavefield, then solve a box-constrained GN model."""
    begin = perf_counter()
    x = np.zeros(penalty.n ** 2)
    weight = 1 / (reference * np.sqrt(y.size))
    history, snapshots = [], []
    for iteration in range(iterations):
        predicted, jac = fullwave(models, x, jacobian=True, material_factors=material_factors)
        residual = (predicted - y).ravel() * weight
        jac_real = np.ascontiguousarray(np.concatenate((jac.real, jac.imag)) * weight)
        residual_real = np.r_[residual.real, residual.imag]
        H = jac_real.T @ jac_real
        reg, reg_grad = penalty.evaluate(x, H)
        objective = 0.5 * energy(residual_real) + reg
        gradient = jac_real.T @ residual_real + reg_grad
        step = _bound_step(H, gradient, x)
        if np.max(np.abs(step)) < 5e-5:
            break
        directional = gradient @ step
        for power in range(9):
            alpha = 0.5 ** power
            candidate = np.clip(x + alpha * step, *BOUNDS)
            new_prediction = fullwave(models, candidate, material_factors=material_factors)
            new_data = 0.5 * energy((new_prediction - y) * weight)
            new_obj = new_data + penalty.evaluate(candidate)[0]
            if new_obj <= objective + 1e-4 * alpha * directional:
                break
        else:
            break
        x = candidate
        runtime = perf_counter() - begin
        history.append((iteration + 1, runtime, new_obj, new_data, 0.0))
        snapshots.append(x.copy())
        if objective - new_obj <= 1e-8 * max(1., objective):
            break
    return x, history, snapshots


@lru_cache(maxsize=None)
def fourier_dictionary(n, halfwidth):
    """Orthogonal IFFT columns for the four low-frequency corner blocks."""
    indices = np.r_[np.arange(halfwidth), np.arange(n - halfwidth, n)]
    rows, cols = np.meshgrid(indices, indices, indexing="ij")
    impulses = np.zeros((n, n, rows.size), dtype=complex)
    impulses[rows.ravel(), cols.ravel(), np.arange(rows.size)] = 1
    return ifft2(impulses, axes=(0, 1), norm="ortho").reshape(n * n, -1)


def _tsom_subspaces(models, y, sigma, relative_cut):
    subspaces = []
    for h, data in zip(models, y):
        left, singular, adj = svd(h.Sh, full_matrices=False, check_finite=False)
        projected = left.conj().T @ data
        signal = np.sqrt(np.mean(abs(projected) ** 2, axis=1))
        keep = np.flatnonzero((singular >= relative_cut * singular[0]) & (signal > 2.75 * sigma))
        # Keep a leading subspace. Adjacent weaker modes may lie below the
        # noise floor, so average across illuminations before selecting a prefix.
        rank = int(keep[-1] + 1) if keep.size else 0
        bright = adj[:rank].conj().T
        deterministic = bright @ (projected[:rank] / singular[:rank, None])
        subspaces.append((bright, deterministic, float(singular[0])))
    return subspaces


def fft_tsom(models, y, reference, sigma, penalty, *, stages=((3, 14), (5, 30)),
             relative_cut=0.13, state_tradeoff=0.32, material_factors=None):
    """Solve the projected-current FFT-TSOM cost by variable projection.

    Each complex current minimizer uses the LS normal matrix of the data and
    state equations. This is an exact update of the current block at fixed
    real material. The envelope gradient updates all material cells at once.
    """
    begin = perf_counter()
    x = np.zeros(penalty.n ** 2)
    data_weight = 1 / (reference * np.sqrt(y.size))
    subspaces = _tsom_subspaces(models, y, sigma, relative_cut)
    if material_factors is None:
        material_factors = np.ones(len(models))
    history, snapshots = [], []
    evaluations = 0
    for halfwidth, iterations in stages:
        fft_columns = fourier_dictionary(penalty.n, halfwidth)
        systems = []
        for h, data, (bright, det, sh_norm), multiplier in zip(
            models, y, subspaces, material_factors, strict=True
        ):
            # Explicit products of the small IFFT dictionary and the finite
            # host are precomputed once per stage. The host is not circulant.
            state_weight = state_tradeoff * sh_norm * data_weight
            projected = fft_columns - bright @ (bright.conj().T @ fft_columns)
            # Compress redundant coordinates, e.g. at the full Fourier span.
            # This preserves the range of (I - V_L V_L*) F_low.
            basis, singular, _ = svd(projected, full_matrices=False, check_finite=False)
            basis = basis[:, singular > 1e-9]
            response = multiplier * (h.Gh @ basis)
            field_det = multiplier * (h.Eh + h.Gh @ det)
            sensor = data_weight * (h.Sh @ basis)
            rhs_sensor = data_weight * (data - h.Sh @ det)
            systems.append((state_weight, det, basis, response, field_det, sensor, rhs_sensor,
                            sensor.conj().T @ sensor, sensor.conj().T @ rhs_sensor))

        last_data = last_state = 0.0

        def reduced(material):
            nonlocal evaluations, last_data, last_state
            evaluations += 1
            reg, gradient = penalty.evaluate(material)
            material = material.reshape(-1, 1)
            data_loss = state_loss = 0.0
            for (ws, det, basis, response, field_det, sensor, rhs_sensor, gram_sensor,
                 cross_sensor) in systems:
                state = ws * (basis - material * response)
                rhs = ws * (material * field_det - det)
                gram = gram_sensor + state.conj().T @ state
                factor = cho_factor(gram, overwrite_a=True, check_finite=False)
                coefficients = cho_solve(
                    factor, cross_sensor + state.conj().T @ rhs, check_finite=False
                )
                residual_state = state @ coefficients - rhs
                residual_data = sensor @ coefficients - rhs_sensor
                total_field = field_det + response @ coefficients
                gradient -= ws * np.sum((total_field.conj() * residual_state).real, axis=1)
                data_loss += 0.5 * energy(residual_data)
                state_loss += 0.5 * energy(residual_state)
            last_data, last_state = data_loss, state_loss
            return reg + data_loss + state_loss, gradient

        def record(intermediate_result):
            runtime = perf_counter() - begin
            history.append((evaluations, runtime, intermediate_result.fun, last_data, last_state))
            snapshots.append(intermediate_result.x.copy())

        solution = minimize(
            reduced, x, method="L-BFGS-B", jac=True, bounds=[BOUNDS] * x.size, callback=record,
            options={"maxiter": iterations, "ftol": 1e-10, "gtol": 5e-7, "maxcor": 12},
        )
        x = solution.x
    ranks = [bright.shape[1] for bright, _, _ in subspaces]
    return x, history, snapshots, ranks
