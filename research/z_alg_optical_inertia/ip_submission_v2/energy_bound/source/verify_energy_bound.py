"""Verify a simultaneous, noise-ball-dependent cap on coherent native currents."""

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, cg
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]
SINGLE = STUDY / "submission/single"
sys.path.insert(0, str(SINGLE))

from coherent_dual import expand_block
from general_fourier_verify import (GAP, PROOFS, conjugate_values, exact_block_diagonal,
                                    exact_compact_matrix, prove_inertia, scale_rows)
from local_moment import (CAP, N, SIDE, ExactGeometry, NumericGeometry, ball, blocks,
                          elapsed, gershgorin, h, hi, load, lo, mid, star)
from certificate import frobenius_squared

SHIFT = 9 / 1024
TARGET_GAP = 12 / 1024
GAP_PROFILES = {
    "01": (16, 1 / 16, 278 / 1024),
    "05": (16, 1 / 16, 220 / 1024),
}
PRECISION = 128
ARB_THREADS = 4


def compact_gap(exact, proposal, configuration, gap):
    """Repeat the frozen Fourier-envelope construction with a larger target gap."""
    _, obj, _, _, _, zr, zi, contrast = configuration
    cutoff, beta_value, gamma_value = GAP_PROFILES[obj]
    center, delta = acb(zr, zi), arb(contrast)
    beta, gamma = arb(beta_value), arb(gamma_value)
    _, _, _, multiplier, diagonal = exact_block_diagonal(proposal, gamma, beta, center, delta)
    diagonal = [d - (arb(gap) - arb(GAP)) for d in diagonal]
    if not all(d > 0 for d in diagonal):
        raise ArithmeticError("The local diagonal after shifting must remain positive.")
    compact, nnegative, omitted = exact_compact_matrix(
        exact, multiplier, diagonal, len(exact.symbol), cutoff, gamma, beta, center, delta,
        proposal["mu"])
    return dict(target_gap=gap, compact_dimension=compact.nrows(),
                cutoff=cutoff, beta=beta_value, gamma=gamma_value, center=[zr, zi],
                contrast=contrast,
                omitted_Schur_slack_lower=omitted,
                diagonal_minimum_lower=lo(min(diagonal)), **prove_inertia(compact, nnegative))


def stationary_proposal(numeric, y, proposal, shift):
    mu = proposal["mu"]
    sector = expand_block(proposal["sector"])
    optical = expand_block(proposal["optical"])
    z = (CAP * sector + 1j * optical) / 2
    receiver, adjoint = numeric.s, numeric.s.conj().T
    attenuation = numeric.sigma**2 / (mu + numeric.sigma**2)
    symbol_adj = numeric.symbol.conj()

    def action(x):
        return ((sector - shift) * x + adjoint @ (receiver @ x) / mu
                - z * numeric.green(x) - numeric.green(z.conj() * x, symbol_adj))

    def precondition(x):
        return x - numeric.v @ (attenuation * (numeric.v.conj().T @ x))

    operator = LinearOperator((N, N), matvec=action, dtype=complex)
    pre = LinearOperator((N, N), matvec=precondition, dtype=complex)
    rhs = adjoint @ y / mu + z[:, None] * numeric.b
    columns, iteration_counts = [], []
    for column in rhs.T:
        counts = []
        estimate, info = cg(operator, column, M=pre, rtol=2e-13, atol=0.0, maxiter=600,
                            callback=lambda _: counts.append(1))
        if info:
            raise ArithmeticError("The shifted stationary proposal did not converge.")
        columns.append(estimate)
        iteration_counts.append(len(counts))
    return ball(np.column_stack(columns)), iteration_counts


def exact_form(exact, y, proposal, shift, gap_shifted, estimate):
    """Lower enclosure of the source matrix from an arbitrary native trial."""
    mu = arb(proposal["mu"])
    sector = [arb(float(w)) for w in expand_block(proposal["sector"])]
    optical = [arb(float(t)) for t in expand_block(proposal["optical"])]
    z = [acb(arb(CAP) * w, t) / 2 for w, t in zip(sector, optical, strict=True)]
    rotated = scale_rows(estimate, conjugate_values(z))
    applied = (scale_rows(estimate, [w - arb(shift) for w in sector])
               + exact.s_star * (exact.s * estimate) / mu
               - scale_rows(exact.green(estimate), z) - exact.green(rotated, adjoint=True))
    rhs = exact.s_star * y / mu + scale_rows(exact.b, z)
    defect = frobenius_squared(applied - rhs).upper()
    correction = (mu * defect / gap_shifted).upper()
    form = star(y) * y + mu * (star(estimate) * applied - star(estimate) * rhs
                               - star(rhs) * estimate)
    for i in range(form.nrows()):
        form[i, i] -= correction
    return form, defect.sqrt(), correction


def uniform_negative_bound(form):
    """Return M>=0 such that the enclosed Hermitian form is >= -M I."""
    values, vectors = eigh(h(mid(form)), check_finite=False)
    basis = ball(vectors)
    gram_floor = gershgorin(star(basis) * basis, upper=False)
    if not gram_floor > 0:
        raise ArithmeticError("The proposed source basis could not be validated.")
    congruence_floor = gershgorin(star(basis) * form * basis, upper=False)
    bound = (max(arb(0), -congruence_floor) / gram_floor).upper()
    return bound, lo(gram_floor), lo(congruence_floor), float(values[0])


def verify_case(exact, numeric, y, radius, proposal, gap, shift):
    started, cpu = perf_counter(), process_time()
    exact_start, numeric_start = exact.green_columns, numeric.green_columns
    shifted_gap = arb(gap) - arb(shift)
    if not (arb(shift) > 0 and shifted_gap > 0):
        raise ArithmeticError("Shift must be strictly between 0 and the proved unshifted gap.")
    estimate, counts = stationary_proposal(numeric, y, proposal, shift)
    form, defect, correction = exact_form(exact, ball(y), proposal, shift, shifted_gap, estimate)
    source_negative, gram_floor, congruence_floor, value = uniform_negative_bound(form)
    cap = (arb(radius)**2 + source_negative) / (arb(proposal["mu"]) * arb(shift))
    return dict(energy_upper=hi(cap), epsilon_noise_ball_radius=radius,
                unshifted_scaled_K_gap_lower=gap, scaled_energy_shift=shift,
                shifted_scaled_gap_lower=lo(shifted_gap),
                stationary_defect_frobenius_upper=hi(defect),
                source_matrix_lower=-hi(source_negative),
                source_trial_congruence_floor=congruence_floor,
                source_basis_gram_floor=gram_floor,
                nominal_smallest_source_eigenvalue=value,
                defect_source_correction_upper=hi(correction),
                cg_iterations_by_source=counts,
                arb_green_columns=exact.green_columns - exact_start,
                numerical_green_columns=numeric.green_columns - numeric_start,
                **elapsed(started, cpu))


def run(output, replay_gap, target_gap, shift):
    started, cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = PRECISION, ARB_THREADS
    s, b, symbol, observations, baseline = load()
    frozen = json.loads((SINGLE / "coherent_general_verified.json").read_text())["cases"]
    radii = {row["id"]: row["inflated_radius"] for row in baseline["observations"]}
    exact = ExactGeometry(s, b, symbol, blocks(4))
    numeric = NumericGeometry(s, b, symbol, blocks(4))
    setup = elapsed(started, cpu)
    records = {}
    for configuration in PROOFS:
        index, obj, *_ = configuration
        case = next(key for key in observations if f"_{obj}__" in key)
        prior = frozen[case]
        if prior["profile_index"] != index or not prior["proved_K_scaled_gap"] == GAP:
            raise ValueError("The frozen gap receipt is inconsistent with the selected profile.")
        proposal = prior["multipliers"]
        proof = None
        if replay_gap:
            gap_start, gap_cpu = perf_counter(), process_time()
            proof = compact_gap(exact, proposal, configuration, target_gap)
            proof.update(elapsed(gap_start, gap_cpu))
            gap = target_gap
        else:
            gap = prior["proved_K_scaled_gap"]
        record = verify_case(exact, numeric, observations[case], radii[case], proposal, gap, shift)
        records[case] = dict(profile_index=index, multipliers=proposal, gap_replay=proof, **record)
        print(case, "certified energy cap", record["energy_upper"], flush=True)
    receipt = dict(
        model="exact saved binary64 arrays, 96x96 scalar TM on the original 10x20 full ring",
        data="synthetic/public plus the original padded radii in results/theory/fine_feasibility.json",
        hypothesis="real 0 <= m_i <= 0.75, unit complex coherent a, ||Sj-Ya||_2 <= r",
        scope="each coherent current in the specified noise ball; independent of truth arrays",
        energy_units="unweighted sum_i |j_i|^2 in the archived cell contrast-current coordinates",
        area_weight_factor_m2=(2 / SIDE)**2,
        material_peak_cap=CAP, native_side=SIDE, sources=b.shape[1], receivers=s.shape[0],
        interval_precision_bits=ctx.prec, interval_threads=ctx.threads, numeric_threads=1,
        inherited_gap_receipt="submission/single/coherent_general_verified.json",
        full_state_equation_solves=0, setup=setup, cases=records, **elapsed(started, cpu))
    output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE.parent / "receipts/energy_bound.json")
    parser.add_argument("--inherited-gap", action="store_true")
    parser.add_argument("--gap", type=float, default=TARGET_GAP)
    parser.add_argument("--shift", type=float)
    args = parser.parse_args()
    default_shift = GAP / 2 if args.inherited_gap else SHIFT
    shift = args.shift if args.shift is not None else default_shift
    with threadpool_limits(1):
        run(args.output, not args.inherited_gap, args.gap, shift)
