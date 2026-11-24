"""Propose one legal 12x12-block contrast and its coherent wave from geometry alone."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from verify_physical_current import (CAP, N, NumericGeometry, blocks, elapsed, geometry_data,
                                     select_direction, solve)

TILE_SIDE = 12
SEED = 8789
NUMERICAL_TOLERANCE = 8e-10


def run(output, receipt):
    started, cpu = perf_counter(), process_time()
    sensor, incident, symbol = geometry_data()
    tile_indices = blocks(TILE_SIDE)
    numeric = NumericGeometry(sensor, incident, symbol, tile_indices)
    label = np.empty(N, dtype=int)
    for i, pixels in enumerate(tile_indices):
        label[pixels] = i
    uniform = np.full(N, CAP)
    selection = select_direction(numeric, uniform, rtol=NUMERICAL_TOLERANCE)
    b = incident @ selection
    current = adjoint = None

    def objective(parameters):
        nonlocal current, adjoint
        m = parameters[label]
        current = solve(numeric, m, m * b, x0=current, rtol=NUMERICAL_TOLERANCE)
        adjoint = solve(numeric, m, current, adjoint=True, x0=adjoint,
                        rtol=NUMERICAL_TOLERANCE)
        field = b + numeric.green(current)
        energy = float(np.vdot(current, current).real)
        gradient = 2 * np.real(adjoint.conj() * field)
        block_gradient = np.bincount(label, weights=gradient)
        return -energy / 1000, -block_gradient / 1000

    initial = np.random.default_rng(SEED).uniform(.4, CAP, TILE_SIDE**2)
    result = minimize(objective, initial, jac=True, method="L-BFGS-B",
                      bounds=[(0.0, CAP)] * TILE_SIDE**2,
                      options=dict(maxiter=40, ftol=1e-9, gtol=1e-7))
    material = result.x[label]
    illumination = select_direction(numeric, material, rtol=NUMERICAL_TOLERANCE)
    np.savez_compressed(output, m=material, a=illumination, params=result.x)
    report = dict(
        scope="geometry-only numerical proposal; certify using verify_physical_current.py",
        block_side=TILE_SIDE, material_cap=CAP, seed=SEED,
        iterations=int(result.nit), forward_adjoint_objective_evaluations=int(result.nfev),
        stop_status=int(result.status), stop_message=str(result.message),
        numerical_green_columns=numeric.green_columns,
        state_and_adjoint_rhs=2 * (incident.shape[1] + int(result.nfev)),
        **elapsed(started, cpu))
    receipt.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output, args.receipt)
