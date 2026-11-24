"""Enclose native currents attained by explicit materials inside the contrast box."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, lgmres
from threadpoolctl import threadpool_limits

from verify_energy_bound import (STUDY, SINGLE, CAP, N, SIDE, ExactGeometry, NumericGeometry,
                                 ball, blocks, elapsed, expand_block, frobenius_squared,
                                 h, hi, lo, scale_rows)

HERE = Path(__file__).resolve().parent
RECEIPTS = HERE.parent / "receipts"
GEOMETRY = STUDY / "synthetic/public/fine_geometry.npz"
PROFILE_INDEX = 2


def geometry_data():
    with np.load(GEOMETRY, allow_pickle=False) as f:
        return tuple(f[key] for key in ("PW360_10x20__S", "PW360_10x20__B", "spectrum"))


def solve(numeric, material, rhs, *, adjoint=False, x0=None, rtol=2e-12):
    if adjoint:
        adj_symbol = numeric.symbol.conj()

        def action(x):
            return x - numeric.green(material * x, adj_symbol)
    else:
        def action(x):
            return x - material * numeric.green(x)

    operator = LinearOperator((N, N), matvec=action, dtype=complex)
    answer, info = lgmres(operator, rhs, x0=x0, rtol=rtol, atol=1e-12, maxiter=45, inner_m=35)
    if info:
        raise ArithmeticError("Numerical native physical solve failed to propose a current.")
    return answer


def select_direction(numeric, m, *, rtol=2e-12):
    response = np.column_stack([solve(numeric, m, m * rhs, rtol=rtol) for rhs in numeric.b.T])
    _, vectors = eigh(h(response.conj().T @ response), check_finite=False)
    return vectors[:, -1]


def resolve_coercive_form():
    transfer = json.loads((RECEIPTS / "radiation_transfer.json").read_text())
    profiles = json.loads((SINGLE / "coherent_general_profiles.json").read_text())
    prior = json.loads((SINGLE / "coherent_general_verified.json").read_text())["cases"]
    proposal = profiles[PROFILE_INDEX]
    comparison = next(p for p in prior.values() if p["profile_index"] == PROFILE_INDEX)
    if any(proposal[key] != comparison["multipliers"][key] for key in ("mu", "sector", "optical")):
        raise ValueError("The frozen radiation-transfer and receiver coercivity profiles differ.")
    transfer_floor = arb(transfer["profiles"][str(PROFILE_INDEX)]["normalized_state_floor_lower"])
    if not (transfer_floor > 0 and all(v >= 0 for row in proposal["sector"] for v in row)):
        raise ArithmeticError("A positive, receiver-free state floor and sector weights are required.")
    t_add = arb(transfer["radiation_coefficient"]) / arb(proposal["mu"])
    w, t = (expand_block(proposal[key]) for key in ("sector", "optical"))
    z = [acb(arb(CAP) * arb(float(ww)), arb(float(tt)) + t_add) / 2
         for ww, tt in zip(w, t, strict=True)]
    return z, transfer_floor


def verify_case(exact, numeric, m, a, z, floor):
    if (m.shape != (N,) or a.shape != (numeric.b.shape[1],) or
            not np.isfinite(m).all() or not np.isfinite(a).all() or (m < 0).any() or (m > CAP).any()):
        raise ValueError("The trial material and source must satisfy the native box and channel layout.")
    started, cpu = perf_counter(), process_time()
    numerical_start = numeric.green_columns
    numerical_incident = numeric.b @ a
    numerical_current = solve(numeric, m, m * numerical_incident)
    field = ball((numerical_incident + numeric.green(numerical_current))[:, None])
    x = scale_rows(field, [arb(float(v)) for v in m])
    wave = ball(a[:, None])
    wave_norm = frobenius_squared(wave).sqrt()
    if not wave_norm > 0:
        raise ArithmeticError("The incident coefficient vector must have positive norm.")
    active_weights = [zz if mi > 0 else acb(0) for zz, mi in zip(z, m, strict=True)]
    defect = scale_rows(field - exact.green(x) - exact.b * wave, active_weights)
    error = (2 * frobenius_squared(defect).sqrt() / floor).upper()
    xnorm = frobenius_squared(x).sqrt()
    lower = max(arb(0), ((xnorm - error) / wave_norm).lower())
    upper = ((xnorm + error) / wave_norm).upper()
    return dict(
        native_material_min=float(m.min()), native_material_max=float(m.max()),
        native_material_mean=float(m.mean()), active_cells=int(np.count_nonzero(m)),
        fixed_a_real=a.real.tolist(), fixed_a_imaginary=a.imag.tolist(),
        coherent_energy_lower=lo(lower**2), coherent_energy_upper=hi(upper**2),
        numerical_energy=float(np.linalg.norm(numerical_current / np.linalg.norm(a))**2),
        weighted_constitutive_defect_frobenius_upper=hi(frobenius_squared(defect).sqrt()),
        current_norm_error_upper=hi(error),
        numerical_green_columns=numeric.green_columns - numerical_start,
        arb_green_columns=1, physical_rhs_for_fixed_direction=1, **elapsed(started, cpu))


def run(output, proposals):
    started, cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = 128, 1
    s, b, spectrum = geometry_data()
    exact = ExactGeometry(s, b, spectrum, blocks(4))
    numeric = NumericGeometry(s, b, spectrum, blocks(4))
    z, floor = resolve_coercive_form()
    setup = elapsed(started, cpu)
    m = np.full(N, CAP)
    selection_start, selection_cpu = perf_counter(), process_time()
    a = select_direction(numeric, m)
    selection = dict(numerical_green_columns=numeric.green_columns,
                     physical_rhs_for_direction_selection=b.shape[1],
                     **elapsed(selection_start, selection_cpu))
    records = {"uniform_0.75": dict(direction_selection=selection,
                                    **verify_case(exact, numeric, m, a, z, floor))}
    for path in proposals:
        with np.load(path, allow_pickle=False) as f:
            m, a = f["m"], f["a"]
        records[path.stem] = verify_case(exact, numeric, m, a, z, floor)
    for key, record in records.items():
        print(key, record["coherent_energy_lower"], record["coherent_energy_upper"], flush=True)
    receipt = dict(
        scope="whole box scale, fixed explicit legal media and fixed unit waves; data and truth unused",
        energy_units="unweighted sum_i |j_i|^2 on 96x96 native cells",
        native_side=SIDE, material_cap=CAP, source_count=b.shape[1],
        no_receiver_coercivity_receipt="radiation_transfer.json",
        coercivity_profile_index=PROFILE_INDEX, coercive_floor=lo(floor),
        interval_precision_bits=ctx.prec, arb_threads=ctx.threads, numeric_threads=1,
        setup=setup, cases=records, **elapsed(started, cpu))
    output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RECEIPTS / "physical_current.json")
    parser.add_argument("--proposal", type=Path, action="append", default=[])
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output, args.proposal)
