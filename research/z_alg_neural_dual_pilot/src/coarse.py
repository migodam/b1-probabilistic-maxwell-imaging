"""Full-grid binary64 proposal objective for the 32² pilot.

The returned scores are numerical. The independent 96² verifier decides
which proposals can make a physical claim.
"""

from collections import OrderedDict
from copy import copy
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh
from scipy.special import expit
from threadpoolctl import threadpool_limits

from chart import MU_SCALE, coordinates, profile, pullback


PAPER = Path(__file__).resolve().parents[2] / "z_alg_optical_inertia"
GEOMETRY = PAPER / "synthetic" / "public"


def _load_paper_inertia():
    spec = spec_from_file_location("_paper1_inertia", PAPER / "core" / "inertia.py")
    module = module_from_spec(spec)
    previous = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


_inertia = _load_paper_inertia()


class Coarse:
    """Full-grid score and gradients with respect to physical block multipliers.

    Y is complex with shape (20 receivers, 10 sources). u is the nonnegative
    scalar material cap. mu and lam are real scalars or equal-length arrays of
    1, 16 or 64 block weights. mu must be nonnegative and lam is signed.
    Block index is y * blocks_per_axis + x.

    smooth softens positive-part clipping with tau=.01*||Y||_F²/10.
    The inertia and cutoff index still use hard counting. The hard squared
    distance in the stored data units is independently available as d2_exact.
    When no source mode is retained, smooth uses the negative value of the
    s-th smallest eigenvalue of K to guide the search back across the boundary.
    """

    def __init__(self, geometry=GEOMETRY):
        started = _inertia.clock()
        with np.load(Path(geometry) / "coarse_green.npz", allow_pickle=False) as green:
            self.G = green["G"]
            self.centers = green["centers"]
        with np.load(Path(geometry) / "pw360_10x20.npz", allow_pickle=False) as ports:
            self.S = ports["S"]
            self.B = ports["B"]
        self.s = self.B.shape[1]
        self.r = self.S.shape[0]
        if self.B.shape != (32**2, 10) or self.S.shape != (20, 32**2):
            raise ValueError("Expected the frozen 32², 20-receiver, 10-source geometry.")
        self.S_star = _inertia.adjoint(self.S)
        with threadpool_limits(limits=1, user_api="blas"):
            self._geometry = _inertia.Prepared(
                self.G, self.S, self.B, np.zeros((self.r, self.s), dtype=complex)
            )

        coord_lo = self.centers.min(axis=0)
        relative = (self.centers - coord_lo) / np.ptp(self.centers, axis=0)
        self.block_ids = {}
        for width in (1, 4, 8):
            axes = np.minimum((relative * width).astype(int), width - 1)
            self.block_ids[width**2] = width * axes[:, 1] + axes[:, 0]
        self._rescue_cache = OrderedDict()
        self.setup_cost = _inertia.timing(started)

    def eval(self, Y, u, mu, lam, smooth=False):
        started = _inertia.clock()
        Y = np.asarray(Y, dtype=complex)
        if Y.shape != (self.r, self.s) or not np.isfinite(Y).all():
            raise ValueError("Y must be a finite matrix matching the frozen receiver and source order.")
        if np.ndim(u) != 0:
            raise ValueError("The material cap u must be a scalar.")
        mu, lam = np.atleast_1d(mu), np.atleast_1d(lam)
        if mu.shape != lam.shape or mu.ndim != 1 or len(mu) not in self.block_ids:
            raise ValueError("mu and lambda must have the same length of 1, 16, or 64.")

        cells = self.block_ids[len(mu)]
        with threadpool_limits(limits=1, user_api="blas"):
            problem = copy(self._geometry)
            problem.operations = self._geometry.operations.copy()
            initial_ops = problem.operations.copy()
            problem.Y = Y
            problem.SY = self.S_star @ Y
            problem.Y2 = _inertia.hermitian(_inertia.adjoint(Y) @ Y)
            problem.operations.update(receiver_observation_products=1, observation_grams=1)
            row = problem.evaluate_ldl(u, mu[cells], lam[cells], return_solution=True)
            result = self._score(problem, u, cells, len(mu), row, smooth)
            if smooth and row["nu"] is not None and row["nu"] >= self.s:
                self._rescue(result, problem, u, mu, lam, cells)

        timings = _inertia.timing(started)
        result.update(wall_s=timings["wall_seconds"], cpu_s=timings["cpu_seconds"])
        result.update(operations=dict(problem.operations - initial_ops))
        return result

    def _score(self, problem, u, cells, n_blocks, row, smooth):
        nu = row["nu"]
        currents = row.get("solution_current")
        result = dict(
            d2=0.0, d2_exact=0.0, grad_mu=np.zeros(n_blocks), grad_lam=np.zeros(n_blocks),
            grad_exact_mu=np.zeros(n_blocks), grad_exact_lam=np.zeros(n_blocks),
            nu=nu, eigenvalues=row.get("H_eigenvalues_descending"), active=0,
            status=row["status"], valid=nu is not None, source_cutoff_gap=None,
            source_zero_gap=None, temperature=None,
            KX_minus_L_Frobenius=row.get("KX_minus_L_Frobenius"),
        )
        if currents is None:
            return result

        linear = row["linear_term"]
        schur = _inertia.hermitian(problem.Y2 - _inertia.adjoint(linear) @ currents)
        eig, vectors = eigh(schur, check_finite=False)
        problem.operations["H_eigendecompositions"] += 1
        eig, vectors = eig[::-1], vectors[:, ::-1]
        retained = eig[nu:]
        active = retained > 0.0
        hard_score = float(retained[active].sum())
        result.update(
            d2=hard_score, d2_exact=hard_score, eigenvalues=eig, active=int(active.sum()),
            source_cutoff_gap=float(eig[nu - 1] - eig[nu]) if nu else None,
            source_zero_gap=float(np.min(np.abs(retained))),
        )
        if smooth:
            tau = max(0.01 * float(problem.Y2.trace().real) / self.s, 1.0e-14)
            weights = expit(retained / tau)
            score = np.maximum(retained, 0) + tau * np.logaddexp(0, -np.abs(retained) / tau)
            result.update(d2=float(score.sum()), temperature=tau)
        else:
            weights = active.astype(float)

        # The source spectral weights turn tr(P dH) into weighted cell sums.
        local_modes = vectors[:, nu:]
        J = currents @ local_modes
        E = (self.B + self.G @ currents) @ local_modes
        mixed = J.conj() * E
        real_sector = np.abs(J)**2 - u * mixed.real
        result["grad_mu"] = np.bincount(cells, weights=real_sector @ weights, minlength=n_blocks)
        result["grad_lam"] = np.bincount(cells, weights=mixed.imag @ weights, minlength=n_blocks)
        exact_weights = active.astype(float)
        result["grad_exact_mu"] = np.bincount(
            cells, weights=real_sector @ exact_weights, minlength=n_blocks,
        )
        result["grad_exact_lam"] = np.bincount(
            cells, weights=mixed.imag @ exact_weights, minlength=n_blocks,
        )
        return result

    def _rescue(self, result, problem, u, mu, lam, cells):
        """Directional guidance in the region where negative modes exhaust sources."""
        key = (float(u), tuple(mu), tuple(lam))
        if key in self._rescue_cache:
            self._rescue_cache.move_to_end(key)
            edge, gap, grad_mu, grad_lam = self._rescue_cache[key]
            problem.operations["K_rescue_cache_hits"] += 1
        else:
            K, _ = problem.matrices(u, mu[cells], lam[cells])
            eigenvalues, vectors = eigh(
                K, subset_by_index=(self.s - 2, self.s), driver="evr", check_finite=False,
                overwrite_a=True,
            )
            problem.operations["dense_K_rescue_eigendecompositions"] += 1
            edge = float(eigenvalues[1])
            gap = float(np.min(np.diff(eigenvalues)))
            mode = vectors[:, 1]
            mixed = mode.conj() * (self.G @ mode)
            grad_mu = np.bincount(
                cells, weights=np.abs(mode)**2 - u * mixed.real, minlength=len(mu),
            )
            grad_lam = np.bincount(cells, weights=mixed.imag, minlength=len(mu))
            self._rescue_cache[key] = (edge, gap, grad_mu, grad_lam)
            if len(self._rescue_cache) > 32:
                self._rescue_cache.popitem(last=False)

        # Scale to the mean squared source-data norm. Tanh keeps the contribution
        # bounded when a profile is far from admitting any retained source modes.
        energy = float(problem.Y2.trace().real) / self.s
        clipped = np.tanh(edge / MU_SCALE)
        factor = (energy / MU_SCALE) * (1.0 - clipped**2)
        result.update(
            d2=float(energy * clipped), grad_mu=factor * grad_mu,
            grad_lam=factor * grad_lam, rescue_K_eigenvalue=edge, rescue_K_gap=gap,
        )

    def optimize(self, Y, u, blocks=16, *, warm_start=None, time_limits=(1.0, 5.0),
                 max_evals=None):
        """Paid deterministic ascent, with best hard-score profiles at wall-clock prefixes.

        Warm starts are (mu, lam) in physical units and are paid initial evaluations.
        Each accepted trial uses a strictly higher hard score or, on a zero plateau,
        a strictly higher smooth guide. Returned candidates maximize the hard
        score. On an all-zero plateau the guide breaks ties among valid profiles.
        """
        wall_start, cpu_start = perf_counter(), process_time()
        limits = sorted(set(map(float, time_limits)))
        if any(not np.isfinite(limit) or limit <= 0 for limit in limits):
            raise ValueError("Time limits must be finite positive seconds.")
        if max_evals is not None and (int(max_evals) != max_evals or max_evals < 1):
            raise ValueError("max_evals must include at least the paid start.")
        if not limits and max_evals is None:
            raise ValueError("An optimization requires a time or evaluation limit.")
        if blocks not in self.block_ids:
            raise ValueError("Unknown block count.")

        if warm_start is None:
            a = np.zeros(blocks)
            b = np.zeros(blocks)
            mu, lam = profile(a, b)
        else:
            mu, lam = (np.atleast_1d(x).astype(float) for x in warm_start)
            if mu.shape != (blocks,) or lam.shape != (blocks,):
                raise ValueError("Warm start must have the requested number of blocks.")
            a, b = coordinates(mu, lam)

        def elapsed():
            return perf_counter() - wall_start, process_time() - cpu_start

        def capture(row, mu, lam, step, accepted):
            wall, cpu = elapsed()
            return dict(
                mu=mu.tolist(), lam=lam.tolist(), d2=row["d2_exact"],
                guide=row["d2"], nu=row["nu"], active=row["active"],
                status=row["status"], valid=row["valid"], accepted=accepted,
                trust_step=step, elapsed_wall_s=wall, elapsed_cpu_s=cpu,
                eval_wall_s=row["wall_s"], eval_cpu_s=row["cpu_s"],
                operations=row["operations"], rescue_K_eigenvalue=row.get("rescue_K_eigenvalue"),
            )

        current = self.eval(Y, u, mu, lam, smooth=True)
        history = [capture(current, mu, lam, 0.0, True)]
        radius = 0.5
        stop = "time_limit"
        maximum_time = limits[-1] if limits else np.inf
        while elapsed()[0] < maximum_time:
            if max_evals is not None and len(history) >= max_evals:
                stop = "evaluation_limit"
                break
            hard = current["d2_exact"] > 0
            grad_mu, grad_lam = pullback(
                a,
                current["grad_exact_mu" if hard else "grad_mu"],
                current["grad_exact_lam" if hard else "grad_lam"],
            )
            scale = max(np.max(np.abs(grad_mu)), np.max(np.abs(grad_lam)))
            if not np.isfinite(scale) or scale == 0:
                stop = "no_usable_gradient"
                break

            trial_a = a + radius * (grad_mu / scale)
            trial_b = b + radius * (grad_lam / scale)
            trial_mu, trial_lam = profile(trial_a, trial_b)
            trial = self.eval(Y, u, trial_mu, trial_lam, smooth=True)
            accepted = trial["valid"] and (
                not current["valid"]
                or trial["d2_exact"] > current["d2_exact"]
                or (not hard and trial["d2_exact"] == 0 and trial["d2"] > current["d2"])
            )
            history.append(capture(trial, trial_mu, trial_lam, radius, accepted))
            if accepted:
                a, b, current = trial_a, trial_b, trial
                radius = min(2.0, 1.35 * radius)
            else:
                radius *= 0.5
                if radius < 2.0**-14:
                    stop = "trust_radius_floor"
                    break

        final_wall, final_cpu = elapsed()
        limit_summaries = {
            limit: _prefix(history, limit, final_wall) for limit in limits
        }
        final = _prefix(history, maximum_time, final_wall)
        return dict(
            **final, prefixes=limit_summaries, paid_wall_s=final_wall, paid_cpu_s=final_cpu,
            eval_count_total=len(history), trials_total=len(history) - 1,
            accepted_updates_total=sum(item["accepted"] for item in history[1:]),
            over_budget_evals=len(history) - final["eval_count"],
            stop_reason=stop, trace=history,
        )

    def refine(self, Y, u, warm_start, *, steps=3):
        """Charge one warm-start score and up to three trial updates."""
        blocks = len(np.atleast_1d(warm_start[0]))
        return self.optimize(
            Y, u, blocks, warm_start=warm_start, time_limits=(), max_evals=1 + steps,
        )


def _prefix(history, limit, final_wall):
    completed = [item for item in history if item["elapsed_wall_s"] <= limit]
    valid = [item for item in completed if item["valid"]]
    first_best = max(valid, default=None, key=lambda item: item["d2"])
    best = max(
        valid, default=None, key=lambda item: (item["d2"], item["guide"] if item["d2"] == 0 else 0),
    )
    return dict(
        best=best, eval_count=len(completed),
        zero_plateau_tiebreak_used=best is not first_best,
        accepted_updates=sum(item["accepted"] for item in completed[1:]),
        completed_wall_s=completed[-1]["elapsed_wall_s"] if completed else 0.0,
        completed_cpu_s=completed[-1]["elapsed_cpu_s"] if completed else 0.0,
        wall_budget_used_s=min(limit, final_wall),
    )
