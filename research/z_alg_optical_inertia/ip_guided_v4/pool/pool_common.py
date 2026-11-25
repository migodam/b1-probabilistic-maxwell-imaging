"""Common physical solver, branch-constrained refiner and atomic evidence writer."""

from collections import Counter
import json
import os
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FLOW = ROOT / "Flow Matching" / "FlowMultiDecode_Local"
RESULTS = Path(os.environ.get("GUIDED_POOL_OUTPUT_ROOT", HERE))


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return {"wall_seconds": perf_counter() - start[0],
            "cpu_seconds": process_time() - start[1]}


def write_json(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def write_arrays(path, **arrays):
    temporary = path.with_suffix(".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)


class FullWave:
    """Exact active restriction of the stored VIE. One same-matrix LU for all RHS."""

    def __init__(self, green, sensor, incident, observed):
        self.green = green
        self.sensor = sensor
        self.incident = incident
        self.observed = observed
        self.scale = np.linalg.norm(observed) ** 2
        self.counts = Counter()

    def evaluate(self, material, gradient=False, incident=None, sensor=None):
        if incident is None:
            incident = self.incident
        if sensor is None:
            sensor = self.sensor
        system = -material[:, None] * self.green
        system.flat[::len(material) + 1] += 1
        factor = lu_factor(system, check_finite=False)
        state = lu_solve(factor, material[:, None] * incident, check_finite=False)
        self.counts.update(lu=1, forward_rhs=incident.shape[1])
        reading = sensor @ state
        if not gradient:
            return reading, state, factor
        residual = reading - self.observed
        adjoint = lu_solve(factor, sensor.conj().T @ residual, trans=2, check_finite=False)
        self.counts.update(adjoint_rhs=incident.shape[1])
        field = incident + self.green @ state
        derivative = np.real(np.sum(adjoint.conj() * field, axis=1))
        loss = np.vdot(residual, residual).real
        return .5 * loss / self.scale, derivative / self.scale

    def refine(self, start, lower, upper, max_calls):
        """Same native, box-constrained objective for every screening arm."""
        if not (np.all(lower <= start) and np.all(start <= upper)):
            raise ValueError("Starting material is outside the branch box.")
        begin = clock()
        counts_before = self.counts.copy()
        best = [np.inf, None]
        calls = 0

        class BudgetReached(Exception):
            pass

        def obj(material):
            nonlocal calls
            if calls >= max_calls:
                raise BudgetReached
            calls += 1
            value, gradient = self.evaluate(material, gradient=True)
            if value < best[0]:
                best[:] = [value, material.copy()]
            return value, gradient

        try:
            result = minimize(obj, start, jac=True, method="L-BFGS-B",
                              bounds=list(zip(lower, upper, strict=True)),
                              options=dict(ftol=1e-12, gtol=1e-10, maxcor=8))
            reason = str(result.message)
        except BudgetReached:
            reason = "full-wave objective evaluation limit"
        material = best[1]
        data, state, _ = self.evaluate(material)
        residual = float(np.linalg.norm(data - self.observed))
        defect = np.linalg.norm(state - material[:, None] * (self.incident + self.green @ state))
        if not np.isfinite(residual) or not np.all((lower <= material) & (material <= upper)):
            raise ArithmeticError("Refinement produced an invalid endpoint.")
        return material, data, dict(residual=residual, state_equation_defect=float(defect),
                                    counts=dict(self.counts - counts_before),
                                    stop_reason=reason, **elapsed(begin))
