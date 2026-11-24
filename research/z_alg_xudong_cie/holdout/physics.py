"""Historical scalar-TM VIE wiring on the stored crop/pad Fourier operator."""

from collections import Counter
from pathlib import Path
import re
import sys
import warnings

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from scipy.special import hankel1, jv

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.z_alg_xudong_cie.holdout.data import K, RESEARCH, WIDTH, elapsed, timers

R8 = RESEARCH / "delegated/max_struct_r1/r8_worktree"
SRAP = RESEARCH / "srap_cie_array_study_20261121"


class NativeForwardError(ArithmeticError):
    def __init__(self, message, cost):
        super().__init__(message)
        self.cost = cost


def historical():
    for path in (R8, SRAP):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    import data2d
    return data2d


def centers(side):
    axis = (np.arange(side) + .5) * (WIDTH / side) - WIDTH / 2
    y, x = np.meshgrid(axis, axis, indexing="ij")
    return np.column_stack((x.ravel(), y.ravel()))


def independent_ports(side, tx_deg, rx_deg):
    x = centers(side)
    tx = np.deg2rad(tx_deg)
    rx = np.deg2rad(rx_deg)
    directions = np.c_[np.cos(tx), np.sin(tx)]
    receiver_positions = 3.75 * np.c_[np.cos(rx), np.sin(rx)]
    radius = WIDTH / side / np.sqrt(np.pi)
    amplitude = .5j * np.pi * K * radius * jv(1, K * radius)
    s = amplitude * hankel1(0, K * np.linalg.norm(receiver_positions[:, None] - x[None, :], axis=2))
    b = np.exp(1j * K * (x @ directions.T))
    return s, b


class NativeForward:
    def __init__(self, spectrum, side):
        self.spectrum, self.side = spectrum, side
        if min(spectrum.shape) < 2 * side - 1:
            raise ValueError("The saved Fourier embedding is too small.")
        self.kernel = ifft2(spectrum)
        self.diagonal = self.kernel[0, 0]
        self.cost = Counter()

    def apply(self, value):
        side = self.side
        ncols = value.reshape(side**2, -1).shape[1]
        self.cost["green_columns"] += ncols
        fields = value.reshape(side, side, ncols)
        padded = fft2(fields, s=self.spectrum.shape, axes=(0, 1))
        result = ifft2(self.spectrum[:, :, None] * padded, axes=(0, 1))[:side, :side]
        return result.reshape(value.shape)

    def stable_dense(self, m, b, s, *, rtol):
        active = np.flatnonzero(m)
        n = len(active)
        y, x = np.divmod(active, self.side)
        p, q = self.spectrum.shape
        local_m = m[active]
        build = timers()
        matrix = np.empty((n, n), dtype=complex, order="F")
        for first in range(0, n, 128):
            stop = min(first + 128, n)
            rows = slice(first, stop)
            dy = (y[rows, None] - y[None, :]) % p
            dx = (x[rows, None] - x[None, :]) % q
            matrix[rows] = -local_m[rows, None] * self.kernel[dy, dx]
        np.fill_diagonal(matrix, 1 - local_m * self.diagonal)
        matrix_build = elapsed(build)
        factor_start = timers()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", LinAlgWarning)
                lu = lu_factor(matrix, check_finite=False, overwrite_a=True)
        except LinAlgWarning as exc:
            raise NativeForwardError(str(exc), dict(fine_active_cells=n, dense_LU_attempts=1,
                                                     matrix_build=matrix_build,
                                                     factorization=elapsed(factor_start))) from exc
        factorization = elapsed(factor_start)
        load = local_m[:, None] * b[active]
        effective_norms = np.maximum(np.linalg.norm(load, axis=0), 1e-300)
        solve_start = timers()
        state = lu_solve(lu, load, check_finite=False)
        full = np.zeros(b.shape, dtype=complex)
        for refinement in range(3):
            full[active] = state
            defect = load - state + local_m[:, None] * self.apply(full)[active]
            residual = float(np.max(np.linalg.norm(defect, axis=0) / effective_norms))
            if residual <= 5 * rtol:
                break
            if refinement < 2:
                state += lu_solve(lu, defect, check_finite=False)
        counts = dict(fine_active_cells=n, forward_rhs=b.shape[1],
                      dense_factor_rhs=b.shape[1] * (1 + refinement), iterative_refinement=refinement,
                      fine_state_residual=residual if np.isfinite(residual) else None,
                      method="active_dense_LU", dense_LU=1, matrix_build=matrix_build,
                      factorization=factorization, solve_and_refinement=elapsed(solve_start))
        if not residual <= 5 * rtol:
            raise NativeForwardError("Dense fallback did not meet its saved-Fourier state defect contract.",
                                     counts)
        return s[:, active] @ state, counts

    def forward(self, m, b, s, *, rtol=2e-9):
        started = timers()
        original_count = self.cost["green_columns"]
        m = np.asarray(m).reshape(-1)
        if (m.shape != (self.side**2,) or b.shape[0] != len(m) or s.shape[1] != len(m) or
                np.iscomplexobj(m) or not np.isfinite(m).all() or np.min(m) < 0):
            raise ValueError("A real nonnegative, native-grid material is required.")
        if not np.any(m):
            response = np.zeros((s.shape[0], b.shape[1]), complex)
            cost = dict(fine_active_cells=0, forward_rhs=0, fine_state_residual=0.0,
                        method="exact_zero_current", dense_LU=0)
        else:
            iterative_start = timers()
            try:
                response, cost = historical().history.fine_solve(
                    m, b, s, self.apply, self.diagonal, rtol=rtol)
                cost.update(method="archived_active_FFT_GMRES", dense_LU=0)
            except RuntimeError as exc:
                incomplete = re.search(r"fine_solve src (\d+):", str(exc))
                failed = dict(reason=str(exc),
                              attempted_rhs=int(incomplete.group(1)) + 1 if incomplete else None,
                              native_fft_green_columns=self.cost["green_columns"] - original_count,
                              **elapsed(iterative_start))
                try:
                    response, cost = self.stable_dense(m, b, s, rtol=rtol)
                except (NativeForwardError, np.linalg.LinAlgError, MemoryError) as dense_error:
                    failed_cost = dict(getattr(dense_error, "cost", {}))
                    failed_cost.update(
                        initial_gmres_failure=failed,
                        native_fft_green_columns=self.cost["green_columns"] - original_count,
                        **elapsed(started))
                    raise NativeForwardError(str(dense_error), failed_cost) from dense_error
                cost.update(initial_gmres_failure=failed)
        return response, {**cost, "native_fft_green_columns": self.cost["green_columns"] - original_count,
                          **elapsed(started)}


def nearest_lift(material, source_side=32, native_side=96):
    factor, rem = divmod(native_side, source_side)
    if rem:
        raise ValueError("The inverse grid does not divide the physical grid.")
    image = np.asarray(material).reshape(source_side, source_side)
    return np.repeat(np.repeat(image, factor, axis=0), factor, axis=1).ravel()
