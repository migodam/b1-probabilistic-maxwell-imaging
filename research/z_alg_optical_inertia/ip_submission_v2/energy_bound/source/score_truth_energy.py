"""Offline enclosure of the coherent-current norm for the saved true media."""

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from scipy.linalg import eigh, lu_solve
from threadpoolctl import threadpool_limits

from verify_energy_bound import (STUDY, CAP, N, SIDE, ExactGeometry, ball, blocks, elapsed,
                                 gershgorin, h, hi, load, lo, mid, scale_rows, star)
from certificate import frobenius_squared

sys.path.insert(0, str(STUDY / "synthetic"))
from nonvacuity import ExactFourierGeometry, StableActiveModel

HERE = Path(__file__).resolve().parent
PRECISION = 128
ARB_THREADS = 1


def score_case(model, exact, b, material):
    started, cpu = perf_counter(), process_time()
    inverse_bound, inverse_proof = model.inverse_bound()
    solve_start, solve_cpu = perf_counter(), process_time()
    current = np.zeros_like(b)
    active = model.active
    current[active, :] = lu_solve(model.lu, material[active, None] * b[active, :], check_finite=False)
    solve_cost = elapsed(solve_start, solve_cpu)

    audit_start, audit_cpu = perf_counter(), process_time()
    estimate = ball(current)
    m = [arb(float(k)) for k in material]
    defect = scale_rows(exact.b + exact.green(estimate), m) - estimate
    state_defect = frobenius_squared(defect).sqrt().upper()
    error_upper = (state_defect * inverse_bound).upper()

    gram = star(estimate) * estimate
    _, vectors = eigh(h(mid(gram)), check_finite=False)
    basis = ball(vectors)
    basis_gram = star(basis) * basis
    basis_floor = gershgorin(basis_gram, upper=False)
    if not basis_floor > 0:
        raise ArithmeticError("The approximate source eigenspace is not certified invertible.")
    top = ball(vectors[:, -1:])
    numerical_lower = (star(top) * gram * top)[0, 0].real / (star(top) * top)[0, 0].real
    numerical_upper = gershgorin(star(basis) * gram * basis, upper=True) / basis_floor
    norm_lower = max(arb(0), (numerical_lower.sqrt() - error_upper).lower())
    norm_upper = (numerical_upper.sqrt() + error_upper).upper()
    result = dict(
        true_material_peak=float(material.max()), truth_above_tested_cap=bool(material.max() > CAP),
        nonzero_true_fine_cells=model.n, coherent_energy_lower=lo(norm_lower**2),
        coherent_energy_upper=hi(norm_upper**2),
        native_numerical_current_energy=float(np.linalg.norm(current, 2)**2),
        per_stored_source_numerical_energies=np.sum(abs(current)**2, axis=0).tolist(),
        unweighted_state_defect_frobenius_upper=hi(state_defect),
        unweighted_state_error_frobenius_upper=hi(error_upper), inverse_proof=inverse_proof,
        new_true_forward_rhs=b.shape[1], inverse_proposal_rhs=model.n,
        arb_green_columns=b.shape[1], native_factorization_count=1,
        solve_cost=solve_cost, audit_cost=elapsed(audit_start, audit_cpu),
        **elapsed(started, cpu))
    return result


def run(output):
    started, cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = PRECISION, ARB_THREADS
    s, b, symbol, observations, _ = load()
    archived = ExactFourierGeometry(STUDY / "synthetic/public/fine_geometry.npz")
    exact = ExactGeometry(s, b, symbol, blocks(4))
    setup = elapsed(started, cpu)
    records = {}
    for key in observations:
        case = key.split("__")[0]
        path = STUDY / "synthetic/evaluation" / f"{case}.npz"
        with np.load(path, allow_pickle=False) as source:
            material = source["m_fine"].reshape(N)
        factor_start, factor_cpu = perf_counter(), process_time()
        model = StableActiveModel(archived, material)
        factor_cost = elapsed(factor_start, factor_cpu)
        record = score_case(model, exact, b, material)
        records[key] = dict(truth_source=str(path.relative_to(STUDY)),
                            factorization_cost=factor_cost, **record)
        print(key, record["coherent_energy_lower"], record["coherent_energy_upper"], flush=True)
    receipt = dict(
        scope="offline score only; none of these truth arrays or energies enter the competitor cap",
        model="exact saved binary64 96x96 scalar TM arrays, 400 MHz, original full ring",
        energy_units="unweighted sum_i |j_i|^2 in the archived cell contrast-current coordinates",
        statistic="max_{||a||_2=1} ||J_true a||_2^2",
        native_side=SIDE, area_weight_factor_m2=(2 / SIDE)**2,
        interval_precision_bits=ctx.prec, interval_threads=ctx.threads, numeric_threads=1,
        setup=setup, cases=records, **elapsed(started, cpu))
    output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE.parent / "receipts/truth_energy_offline.json")
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output)
