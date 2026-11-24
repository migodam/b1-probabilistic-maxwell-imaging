"""Search ring weights for two and three required high pixels at the archived radius."""

import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from area_certificate import (HERE, OBSERVATION_RADII, RING_ID, RING_THRESHOLD, clock,
                              read_json, ring_data, save_json, timing)
from inertia import adjoint, hermitian

SCALE = 1e-5
BASE_PROFILES = (
    (4.571871609296927e-5, -6.845644029318437e-5),
    (5e-6, -2e-5),
)
PARTITION_SIDE = 4


def prepare(model, positions, q, radius):
    side = PARTITION_SIDE
    frac = (positions - positions.min(axis=0)) / np.ptp(positions, axis=0)
    block = np.minimum((frac * side).astype(int), side - 1)
    order = side * block[:, 1] + block[:, 0]
    blocks = side**2
    records = dict(evaluations=0, inertia_rejected=0)
    best = dict(squared_tail=-1.0, parameters=None)

    def objective(parameters):
        records["evaluations"] += 1
        sector = SCALE * parameters[:blocks][order]
        optical = SCALE * parameters[blocks:][order]
        stationary = model.evaluate_ldl(RING_THRESHOLD, sector, optical, return_solution=True)
        if stationary["nu"] != 0 or "solution_current" not in stationary:
            records["inertia_rejected"] += 1
            return 100.0, np.zeros(2 * blocks)
        x = stationary["solution_current"]
        rhs = stationary["linear_term"]
        values, vectors = eigh(hermitian(model.Y2 - adjoint(rhs) @ x), check_finite=False)
        values, vectors = values[::-1], vectors[:, ::-1]
        retained = (np.arange(model.s) >= q) & (values > 0)
        squared = float(np.sum(values[retained]))
        if squared > best["squared_tail"]:
            best.update(squared_tail=squared, parameters=parameters.copy())
        active = vectors[:, retained]
        j = x @ active
        total = (model.B + model.G @ x) @ active
        z = j.conj() * total
        sector_derivative = np.sum((j.conj() * j).real - RING_THRESHOLD * z.real, axis=1)
        optical_derivative = np.sum(z.imag, axis=1)
        gradient = np.r_[np.bincount(order, weights=sector_derivative, minlength=blocks),
                         np.bincount(order, weights=optical_derivative, minlength=blocks)]
        return -squared / radius**2, -SCALE * gradient / radius**2

    return objective, best, records


def run(output):
    start = clock()
    model, positions = ring_data()
    radius = read_json(OBSERVATION_RADII)[RING_ID]["frozen_observation_ball"]
    blocks = PARTITION_SIDE**2
    tracks, profiles = [], []
    for q in (1, 2):
        objective, best, counts = prepare(model, positions, q, radius)
        results = []
        for sector, optical in BASE_PROFILES:
            initial = np.r_[np.full(blocks, sector / SCALE), np.full(blocks, optical / SCALE)]
            result = minimize(objective, initial, jac=True, method="L-BFGS-B",
                              bounds=[(.002, 80.0)] * blocks + [(-120.0, 120.0)] * blocks,
                              options=dict(maxiter=35, maxls=16, ftol=1e-9, gtol=1e-8))
            results.append(dict(nit=int(result.nit), nfev=int(result.nfev), status=int(result.status),
                                message=str(result.message)))
        if best["parameters"] is None:
            raise ArithmeticError("No fixed-positive-inertia proposal was found.")
        params = best["parameters"] * SCALE
        profiles.append(dict(id=f"spatial_q{q}", partition_side=PARTITION_SIDE,
                             mu=params[:blocks].tolist(), **{"lambda": params[blocks:].tolist()}))
        bound = float(np.sqrt(best["squared_tail"]))
        tracks.append(dict(q=q, numerical_tail_residual=bound, trials=results, **counts))
        print(q, "numerical residual", bound, "archived radius", radius, flush=True)
    save_json(output, dict(scope="observed case only, floating point proposals", radius=radius,
                           tracks=tracks, profiles=profiles, operations=dict(model.operations),
                           **timing(start)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "ring_spatial_proposals.json")
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output)
