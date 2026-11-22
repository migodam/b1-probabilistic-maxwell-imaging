"""Outward native verifier for nonuniform sector and optical dual weights."""

import argparse
import json
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import cholesky, eigh
from scipy.sparse.linalg import LinearOperator, cg
from threadpoolctl import threadpool_limits

from coherent_dual import expand_block
from local_moment import (CAP, HERE, N, SIDE, ExactGeometry, NumericGeometry, ball, blocks, elapsed,
                          gershgorin, h, hi, load, lo, mid, norm_upper, star)

GAP = 1 / 1024
# index, object, cutoff, beta, omitted shift, real/imaginary center, contrast
PROOFS = (
    (0, "05", 13, 1/16, 418/1024, 103/256, 39/256, 55/256),
    (2, "01", 13, 1/32, 606/1024, 140/256, 46/256, 50/256),
)


def scale_rows(matrix, values):
    return acb_mat(matrix.nrows(), matrix.ncols(),
                   [v * matrix[i, j] for i, v in enumerate(values) for j in range(matrix.ncols())])


def conjugate_values(values):
    return [value.conjugate() for value in values]


def abs_squared(z):
    return z.real**2 + z.imag**2


def exact_block_diagonal(proposal, gamma, beta, center, contrast):
    w = [arb(v) for row in proposal["sector"] for v in row]
    t = [arb(v) for row in proposal["optical"] for v in row]
    z = [acb(arb(CAP)*a, b)/2 for a, b in zip(w, t, strict=True)]
    m = [(a-center) / contrast for a in z]
    diagonal = [a-gamma-beta*abs_squared(b)-arb(GAP) for a, b in zip(w, m, strict=True)]
    if not all(d > 0 for d in diagonal):
        raise ArithmeticError("The exact local identity term does not stay positive.")
    return w, t, z, m, diagonal


def exact_compact_matrix(exact, m, diagonal, shape, cutoff, gamma, beta, center, contrast, mu):
    padding, nrec = shape, exact.s.nrows()
    modes_1d = range(-cutoff, cutoff+1)
    qs = [(y, x) for y in modes_1d for x in modes_1d]
    keep = {(y % padding, x % padding) for y, x in qs}
    omitted_slack = arb(1000)
    for y in range(padding):
        for x in range(padding):
            if (y, x) in keep:
                continue
            g = exact.symbol[y][x]
            schur = (gamma - 2*(center*g).real - contrast**2 * abs_squared(g)/beta).lower()
            if not schur > 0:
                raise ArithmeticError("An omitted Fourier block is unresolved.")
            omitted_slack = min(omitted_slack, schur)

    roots = [acb(arb(2)*k/padding).exp_pi_i() for k in range(padding)]
    deltas = range(-2*cutoff, 2*cutoff+1)
    block_width = SIDE // 4
    sums = {delta: [sum((roots[delta*k % padding] for k in range(
                   j*block_width, (j+1)*block_width)), acb(0)) for j in range(4)] for delta in deltas}
    weights = [[1/d for d in diagonal],
               [value/d for value, d in zip(m, diagonal, strict=True)],
               [abs_squared(value)/d for value, d in zip(m, diagonal, strict=True)]]
    gram_table = []
    for ww in weights:
        gram_table.append({
            (dy, dx): sum((ww[4*j+k]*sums[dy][j]*sums[dx][k]
                           for j in range(4) for k in range(4)), acb(0))/padding**2
            for dy in deltas for dx in deltas
        })
    block_indices = np.empty(N, dtype=int)
    for i, where in enumerate(exact.indices):
        block_indices[where] = i
    sw = scale_rows(exact.s_star, [weights[0][i] for i in block_indices])
    sm = scale_rows(exact.s_star, [weights[1][i].conjugate() for i in block_indices])
    f = acb_mat(len(modes_1d), SIDE, [roots[-k*j % padding] for k in modes_1d for j in range(SIDE)])
    projections = []
    for field in (sw, sm):
        fields = []
        for j in range(nrec):
            column = acb_mat(SIDE, SIDE, [field[i, j] for i in range(N)])
            fields.append(((f * column * f.transpose())/padding).entries())
        projections.append(fields)
    gram_receivers = exact.s * sw
    rank = len(qs)
    small = acb_mat(2*rank+nrec, 2*rank+nrec)
    p0, p1, p2 = gram_table
    for i, (yi, xi) in enumerate(qs):
        for j, (yj, xj) in enumerate(qs):
            delta = (yj-yi, xj-xi)
            small[i, j] = p0[delta]
            small[i, rank+j] = p1[delta]
            small[rank+j, i] = p1[delta].conjugate()
            small[rank+i, rank+j] = p2[delta]
        for j in range(nrec):
            rec = 2*rank+j
            v, z = projections[0][j][i], projections[1][j][i]
            small[i, rec], small[rec, i] = v, v.conjugate()
            small[rank+i, rec], small[rec, rank+i] = z, z.conjugate()
    for i in range(nrec):
        for j in range(nrec):
            small[2*rank+i, 2*rank+j] = gram_receivers[i, j]
        small[2*rank+i, 2*rank+i] += arb(mu)
    nnegative = 0
    for i, (y, x) in enumerate(qs):
        g = exact.symbol[y % padding][x % padding]
        d00, d01 = gamma - 2*(center*g).real, -contrast*g.conjugate()
        det = beta*d00 - abs_squared(d01)
        if not (det > 0 or det < 0):
            raise ArithmeticError("A kept Fourier block has indeterminate determinant.")
        nnegative += int(det < 0)
        small[i, i] += beta/det
        small[rank+i, rank+i] += d00/det
        small[i, rank+i] -= d01/det
        small[rank+i, i] -= d01.conjugate()/det
    return small, nnegative, lo(omitted_slack)


def prove_inertia(small, nnegative):
    centered = h(mid(small))
    values, modes = eigh(centered, check_finite=False, driver="evd")
    if int((values < 0).sum()) != nnegative:
        raise ArithmeticError("Numerical inertia of the lower envelope is insufficient.")
    shift = float(np.min(abs(values))/2)
    v = modes[:, :nnegative]
    w = v * np.sqrt(-2*values[:nnegative])
    lifted = h(centered + w @ w.conj().T)
    lifted.flat[::len(lifted)+1] -= shift
    l = ball(cholesky(lifted, lower=True, check_finite=False))
    wb, vb = ball(w), ball(v)
    defect = small + wb * star(wb) - l * star(l)
    for i in range(small.nrows()):
        defect[i, i] -= arb(shift)
    bound = norm_upper(defect)
    if not bound < shift:
        raise ArithmeticError("The interval shifted-Cholesky defect is too large.")
    negative = gershgorin(star(vb)*(small*vb + arb(shift)*vb), upper=True)
    if not negative < 0:
        raise ArithmeticError("The proposed negative eigenspace is unresolved.")
    return dict(compact_inertia_negatives=nnegative, compact_defect_upper=bound,
                compact_shift=shift, compact_negative_congruence_upper=hi(negative))


def native_coherent_bound(exact, numeric, mu, z, proposal, a, y):
    wave = ball(a[:, None])
    yb = ball(y) * wave
    bb = exact.b * wave
    idx = np.empty(N, dtype=int)
    for i, where in enumerate(exact.indices):
        idx[where] = i
    z_native = [z[i] for i in idx]
    zd = (CAP*expand_block(proposal["sector"]) + 1j*expand_block(proposal["optical"]))/2
    w = expand_block(proposal["sector"])
    receiver = numeric.s
    receiver_star = receiver.conj().T

    def action(x):
        return (receiver_star@(receiver@x)/mu + w*x - zd*numeric.green(x)
                - numeric.green(zd.conj()*x, numeric.symbol.conj()))

    def precondition(x):
        w = numeric.sigma**2/(mu+numeric.sigma**2)
        return x - numeric.v@(w*(numeric.v.conj().T@x))

    op = LinearOperator((N, N), matvec=action, dtype=complex)
    pre = LinearOperator((N, N), matvec=precondition, dtype=complex)
    rhs = receiver_star@(y@a)/mu + zd*(numeric.b@a)
    counts = []
    estimate, info = cg(op, rhs, M=pre, rtol=2.e-13, atol=0., maxiter=600,
                        callback=lambda _: counts.append(1))
    if info:
        raise ArithmeticError("The native coherent stationary proposal did not converge.")
    x = ball(estimate[:, None])
    ax = (scale_rows(x, [arb(float(v)) for v in w]) + exact.s_star*(exact.s*x)/arb(mu)
          - scale_rows(exact.green(x), z_native)
          - exact.green(scale_rows(x, conjugate_values(z_native)), adjoint=True))
    rhsb = exact.s_star*yb/arb(mu) + scale_rows(bb, z_native)
    defect = norm_upper(ax-rhsb)
    profiled_value = (star(yb)*yb + arb(mu)*(star(x)*ax - star(x)*rhsb - star(rhsb)*x))[0, 0].real
    lower = (profiled_value-arb(mu)*arb(defect)**2/arb(GAP))/(star(wave)*wave)[0, 0].real
    if not lower > 0:
        raise ArithmeticError("No positive native coherent objective was verified.")
    return lower.sqrt(), dict(cg_iterations=len(counts), stationary_defect_upper=defect,
                              normalized_coherent_projected_data_norm_upper=norm_upper(yb))


def run(output, indices):
    start, cpu = perf_counter(), process_time()
    s, b, symbol, observations, baseline = load()
    profile_rows = json.loads((HERE/"coherent_general_native.json").read_text())["rows"]
    radii = {row["id"]: row["inflated_radius"] for row in baseline["observations"]
             if row["id"] in observations}
    ctx.prec, ctx.threads = 104, 4
    exact = ExactGeometry(s, b, symbol, blocks(4))
    numeric = NumericGeometry(s, b, symbol, blocks(4))
    setup = elapsed(start, cpu)
    results = {}
    for index, obj, cutoff, beta, gamma, zr, zi, contrast in PROOFS:
        if index not in indices:
            continue
        row = profile_rows[index]
        case = next(case for case in observations if f"_{obj}__" in case)
        a_data = row["cases"][case]
        a = np.array(a_data["top_coherent_real"]) + 1j*np.array(a_data["top_coherent_imaginary"])
        mu = row["mu"]
        t0, c0 = perf_counter(), process_time()
        center, delta = acb(zr, zi), arb(contrast)
        beta0, gamma0 = arb(beta), arb(gamma)
        _, _, z, m, d = exact_block_diagonal(row, gamma0, beta0, center, delta)
        small, nnegative, omitted = exact_compact_matrix(exact, m, d, len(symbol), cutoff,
                                                         gamma0, beta0, center, delta, mu)
        gram_cost = elapsed(t0, c0)
        t1, c1 = perf_counter(), process_time()
        proof = prove_inertia(small, nnegative)
        inertia_cost = elapsed(t1, c1)
        t2, c2 = perf_counter(), process_time()
        exact_start, numeric_start = exact.green_columns, numeric.green_columns
        lower, station = native_coherent_bound(exact, numeric, mu, z, row, a, observations[case])
        margin = lower - arb(radii[case])
        multipliers = {key: row[key] for key in ("mu", "sector", "optical")}
        result = dict(profile_index=index, multipliers=multipliers, cap=CAP, Fourier_cutoff=cutoff, gamma=gamma, beta=beta,
                      center=[zr, zi], contrast=contrast, compact_dimension=small.nrows(),
                      proved_K_scaled_gap=GAP, omitted_Schur_slack_lower=omitted,
                      diagonal_minimum=lo(min(d)), coherent_lower=lo(lower),
                      radius=radii[case], margin_lower=lo(margin),
                      excludes_at_padded_radius=bool(margin > 0), **proof, **station,
                      native_stationarity_cost=elapsed(t2, c2), exact_fft_columns=exact.green_columns-exact_start,
                      numeric_fft_columns=numeric.green_columns-numeric_start,
                      gram_cost=gram_cost, inertia_cost=inertia_cost, **elapsed(t0, c0))
        print(case, "outward bound", result["coherent_lower"], "margin", result["margin_lower"],
              "seconds", result["wall_seconds"], flush=True)
        results[case] = result
        report = dict(native_side=SIDE, partition=[4, 4], cases=results,
                      interval_precision_bits=ctx.prec, interval_threads=ctx.threads,
                      setup=setup, full_state_equation_solves=0, **elapsed(start, cpu))
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=lambda x: HERE/x, default=HERE/"coherent_general_verified.json")
    parser.add_argument("--profile", action="append", type=int, choices=(0, 2))
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output, args.profile or (0, 2))
