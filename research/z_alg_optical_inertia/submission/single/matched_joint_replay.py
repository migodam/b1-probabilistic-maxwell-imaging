"""Verify the number of positive illumination directions at matched 4x4 weights."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, cg
from threadpoolctl import threadpool_limits

from coherent_dual import expand_block
from general_fourier_verify import conjugate_values, scale_rows
from local_moment import (CAP, HERE, N, ExactGeometry, NumericGeometry, ball, blocks, elapsed,
                          gershgorin, h, hi, load, lo, mid, norm_upper, star)


def replay(exact, numeric, proposal, y):
    start, cpu = perf_counter(), process_time()
    mu = proposal["mu"]
    w = expand_block(proposal["sector"])
    t = expand_block(proposal["optical"])
    z = (CAP*w+1j*t)/2
    s, s_star = numeric.s, numeric.s.conj().T
    receive_factor = numeric.sigma**2/(mu+numeric.sigma**2)

    def action(x):
        return (s_star@(s@x)/mu + w*x - z*numeric.green(x)
                - numeric.green(z.conj()*x, numeric.symbol.conj()))

    def precondition(x):
        return x - numeric.v@(receive_factor*(numeric.v.conj().T@x))

    operator = LinearOperator((N, N), matvec=action, dtype=complex)
    pre = LinearOperator((N, N), matvec=precondition, dtype=complex)
    rhs = s_star@y/mu + z[:, None]*numeric.b
    solutions, solves = [], []
    for column in rhs.T:
        counts = []
        x, info = cg(operator, column, M=pre, rtol=2.e-13, atol=0., maxiter=600,
                     callback=lambda _: counts.append(1))
        if info:
            raise ArithmeticError("Matched joint stationary solve failed to converge.")
        solutions.append(x)
        solves.append(len(counts))

    x, yb = ball(np.column_stack(solutions)), ball(y)
    weight = [arb(float(value)) for value in w]
    zb = [acb(arb(CAP)*ww, arb(float(tt)))/2 for ww, tt in zip(weight, t, strict=True)]
    ax = (scale_rows(x, weight) + exact.s_star*(exact.s*x)/arb(mu)
          - scale_rows(exact.green(x), zb)
          - exact.green(scale_rows(x, conjugate_values(zb)), adjoint=True))
    l = exact.s_star*yb/arb(mu) + scale_rows(exact.b, zb)
    upper = star(yb)*yb + arb(mu)*(star(x)*ax - star(x)*l - star(l)*x)
    vectors = eigh(h(mid(upper)), check_finite=False)[1]
    vneg = ball(vectors[:, :-1])
    nupper = gershgorin(star(vneg)*upper*vneg, upper=True)
    if not nupper < 0:
        raise ArithmeticError("The matched coefficient matrix has unresolved extra positive directions.")
    v = ball(vectors)
    gram_min = gershgorin(star(v)*v, upper=False)
    maximum_upper = gershgorin(star(v)*upper*v, upper=True)
    if not (gram_min > 0 and maximum_upper > 0):
        raise ArithmeticError("The upper eigenvalue bound is unresolved.")
    return dict(
        second_and_lower_source_direction_congruence_upper=hi(nupper),
        maximum_original_source_quadratic_upper=hi(max(upper[i, i].real.upper()
                                                       for i in range(upper.nrows()))),
        joint_trace_distance_upper=hi((maximum_upper/gram_min).sqrt()),
        exact_stationary_matrix_defect_upper=norm_upper(ax-l),
        cg_iterations=solves, **elapsed(start, cpu))


def run(output):
    start, cpu = perf_counter(), process_time()
    profiles = json.loads((HERE/"coherent_general_native.json").read_text())["rows"]
    verified = json.loads((HERE/"coherent_general_verified.json").read_text())["cases"]
    s, b, symbol, observations, _ = load()
    ctx.prec, ctx.threads = 104, 4
    exact = ExactGeometry(s, b, symbol, blocks(4))
    numeric = NumericGeometry(s, b, symbol, blocks(4))
    results = {}
    for case, proof in verified.items():
        if not proof["excludes_at_padded_radius"] or not proof["proved_K_scaled_gap"] > 0:
            raise ArithmeticError("Run the positive Fourier proof first.")
        exact0, numeric0 = exact.green_columns, numeric.green_columns
        profile = profiles[proof["profile_index"]]
        if any(profile[key] != proof["multipliers"][key] for key in ("mu", "sector", "optical")):
            raise ValueError("The native proposal and the verified profile have different multipliers.")
        result = replay(exact, numeric, profile, observations[case])
        result.update(profile_index=proof["profile_index"], positive_source_directions=1,
                      joint_trace_distance_lower=proof["coherent_lower"],
                      warm_exact_fft_columns=exact.green_columns-exact0,
                      warm_numeric_fft_columns=numeric.green_columns-numeric0)
        results[case] = result
        print(case, result["joint_trace_distance_lower"], result["joint_trace_distance_upper"], flush=True)
    report = dict(native_cells=N, partition=[4, 4], cap=CAP, cases=results,
                  proved_K_source="coherent_general_verified.json", full_state_equation_solves=0,
                  **elapsed(start, cpu))
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"matched_joint_replay.json")
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output)
