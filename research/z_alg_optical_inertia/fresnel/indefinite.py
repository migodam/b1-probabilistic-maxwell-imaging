"""Outward verification of useful grouped Fresnel inertia proposals."""
import argparse
import json

import numpy as np
from scipy.linalg import LinAlgError, cholesky, eigh, solve
from threadpoolctl import threadpool_limits

from physical import FREQUENCIES, RESULTS, adjoint, clock, elapsed, hermitian, load_case, operators
from scan import PROFILES, base_matrix, exact_measurement_groups
from verify import base_assembly_error, down, gamma, gemm_error, inflate, norm_upper, rayleigh_bounds, up


def bound_a_gap(a, assembly_error, negative):
    started = clock()
    n = len(a)
    values, modes = eigh(a, subset_by_index=(0, negative), check_finite=False)
    if not values[negative-1] < 0 < values[negative]:
        raise ArithmeticError("The base-proposal inertia is unresolved.")
    shift = float(min(-values[negative-1], values[negative]) / 2.)
    z = modes[:, :negative]
    w = z * np.sqrt(-2. * values[:negative])
    corrections = w @ adjoint(w)
    shift_a = a + corrections
    shift_a.flat[::n+1] -= shift
    f = np.tril(cholesky(shift_a, lower=True, check_finite=False))
    ff = f @ adjoint(f)
    defect = inflate(norm_upper(shift_a-ff) + assembly_error
                     + gemm_error(w, adjoint(w)) + gemm_error(f, adjoint(f))
                     + gamma(32)*(norm_upper(a) + norm_upper(corrections) + norm_upper(ff) + np.sqrt(n)*shift))

    az = a @ z + shift*z
    proj = adjoint(z) @ az
    rows = np.abs(proj.real) + np.abs(proj.imag)
    np.fill_diagonal(rows, 0.)
    row_end = np.sum(rows, axis=1)
    diagonal = proj.diagonal().real
    nr = norm_upper(z)
    neg_error = inflate(nr * (gemm_error(a, z) + assembly_error * nr
                        + gamma(24)*(norm_upper(a)+shift)*nr) + gemm_error(adjoint(z), az)
                        + gamma(8*negative+32)*float(np.max(abs(diagonal)+row_end)))
    negative_upper = up(up(np.max(diagonal+row_end)) + neg_error)
    alpha = down(shift - defect)
    if alpha <= 0. or negative_upper >= 0.:
        raise ArithmeticError(f"Stabilized-base rounds did not close. "
                              f"shift={shift:.6g}, defect={defect:.6g}, negative_upper={negative_upper:.6g}.")
    return alpha, dict(A_negative_verified=negative, A_gap_lower=alpha,
                       A_positive_defect_upper=defect, A_negative_shifted_upper=negative_upper,
                       **elapsed(started))


def corrected_dual_gram(a, rhs, x, a_error, rhs_error, a_gap):
    """A second-order stationary enclosure of R* A^-1 R."""
    na, nr, nx = map(norm_upper, (a, rhs, x))
    ax = a @ x
    ea = inflate(gemm_error(a, x, na, nx) + a_error*nx)
    error_r = inflate(norm_upper(ax-rhs) + ea + rhs_error + gamma(16)*(norm_upper(ax)+nr))
    rx = adjoint(rhs) @ x
    xax = adjoint(x) @ ax
    m = hermitian(rx + adjoint(rx) - xax)
    error_m = inflate(2*nx*rhs_error + nx*ea + error_r**2/a_gap + 2*gemm_error(adjoint(rhs), x, nr, nx)
                      + gemm_error(adjoint(x), ax, nx)
                      + gamma(32)*(2*norm_upper(rx)+norm_upper(xax)))
    return m, error_m, error_r


def receiver_gap(u, base_error):
    """Use a full small eigensystem only as a proposal for inertia of U."""
    d, v = eigh(u, check_finite=False)
    vd = v * d
    vdv = vd @ adjoint(v)
    gram = adjoint(v) @ v
    gram.flat[::len(v)+1] -= 1.
    nv = norm_upper(v)
    ortho = inflate(norm_upper(gram) + gemm_error(adjoint(v), v, nv, nv)
                    + gamma(16)*(nv**2+np.sqrt(len(v))))
    defect = inflate(norm_upper(u-vdv) + base_error + gemm_error(vd, adjoint(v), nb=nv)
                     + gamma(24)*(nv*norm_upper(vd)+norm_upper(u)))
    gap = down(down(float(np.min(abs(d))) * down(1.-ortho)) - defect)
    return d, v, int((d < 0).sum()), gap, defect


def verify_cohort(a, x_all, p, geometry, case, r, tx, mu, a_error, nz, a_gap, a_nu):
    s = geometry["S"][r]
    measured = case["Y"][np.ix_(r, tx)]
    nrows = len(r)
    rhs = np.c_[adjoint(s), p[:, tx]]
    x = np.c_[x_all[:, r], x_all[:, geometry["S"].shape[0]+tx]]
    rhs_error = inflate(gamma(64)*nz*norm_upper(case["B"][:, tx]))
    m, m_error, _ = corrected_dual_gram(a, rhs, x, a_error, rhs_error, a_gap)
    u0 = m[:nrows, :nrows].copy()
    u0.flat[::nrows+1] += mu
    d0 = measured - m[:nrows, nrows:]
    assembly = inflate(gamma(16)*(norm_upper(m)+norm_upper(measured)+np.sqrt(nrows)*mu))
    d, v, u_nu, u_gap, defect = receiver_gap(u0, inflate(m_error+assembly))
    k_nu = a_nu - u_nu
    condition = dict(receiver_update_gap_lower=u_gap, receiver_update_defect_upper=defect,
                     receiver_update_assembly_error_upper=assembly,
                     A_gram_error_upper=m_error, K_negative_numeric=k_nu, channels=measured.size)
    if u_gap <= 0. or k_nu < 0:
        return dict(outward_status="receiver_inertia_unresolved", **condition)
    if k_nu >= len(tx):
        return dict(outward_status="certified_sources_exhausted", K_negative_verified=k_nu,
                    lower_squared=0., **condition)

    q = v @ ((adjoint(v) @ d0) / d[:, None])
    nq, ns, ny = map(norm_upper, (q, s, measured))
    # Evaluate the source-size corrected Gram at Z=S*Q+P itself. This
    # preserves cancellations in S*Q before using a Frobenius error bound.
    small_rhs = adjoint(s) @ q + p[:, tx]
    small_x = x[:, :nrows] @ q + x[:, nrows:]
    small_rhs_error = inflate(gemm_error(adjoint(s), q, ns, nq) + rhs_error
                             + gamma(16)*(ns*nq+norm_upper(p[:, tx])))
    small_gram, small_gram_error, small_defect = corrected_dual_gram(
        a, small_rhs, small_x, a_error, small_rhs_error, a_gap)
    nsx = norm_upper(small_x)
    uq_defect = inflate(norm_upper(mu*q+s@small_x-measured) + ns*small_defect/a_gap
                        + gemm_error(s, small_x, ns, nsx) + gamma(32)*(mu*nq+ns*nsx+ny))
    qy = adjoint(q) @ measured
    qq = adjoint(q) @ q
    h = mu*hermitian(qy + adjoint(qy) - mu*qq - small_gram)
    h_error = inflate(mu*(small_gram_error + uq_defect**2/u_gap
                      + 2*gemm_error(adjoint(q), measured, nq, ny) + mu*gemm_error(adjoint(q), q, nq, nq)
                      + gamma(48)*(2*norm_upper(qy)+mu*norm_upper(qq)+norm_upper(small_gram))))
    _, terms = rayleigh_bounds(h, h_error, skip=k_nu)
    square = max(0., down(sum(terms) * (1.-gamma(len(terms)+24))))
    return dict(outward_status="verified", K_negative_verified=k_nu, lower_squared=square,
                H_error_upper=h_error, source_eigenvalue_floors=terms, **condition)


def verify_candidate(scan, cases, geometry):
    profiles = {p.name: p for p in PROFILES}
    proposals = []
    for record in scan["cases"]:
        delta = record["empty_calibration"]["rms_transferred_over_all_scatter"]
        for row in record["grid"]:
            if not row["best_indefinite"]:
                continue
            for partition in ("adjacent_cohorts", "adjacent_pairs"):
                groups = row["curves"][partition]
                for j, mu in enumerate(scan["mus"]):
                    if mu < 1.e-7:
                        continue
                    per_group = [g["rows"][j] for g in groups]
                    if any(r["status"] == "indefinite" and r["lower_bound"] > 0. for r in per_group):
                        total = float(np.linalg.norm([r["lower_bound"] or 0. for r in per_group]))
                        proposals.append((total, delta, record["name"], row["u"],
                                          row["profile"], partition, mu))
    if not proposals:
        return dict(status="no_positive_grouped_candidate_at_mu_at_least_1e-7")
    total, delta, name, cap, profile, partition, mu = max(proposals)
    started = clock()
    case = cases[name]
    a, z, weight = base_matrix(geometry, cap, profiles[profile])
    a_error, nz = base_assembly_error(geometry, weight, z)
    a_nu = next(r["A_negative"] for r in scan["base_records"]
                if r["u"] == cap and r["profile"] == profile)
    proposal = dict(case=name, u=cap, profile=profile, partition=partition, mu=mu,
                    numerical_bound=total, empirical_empty_RMS_proxy=delta)
    try:
        a_gap, proof = bound_a_gap(a, a_error, a_nu)
        p = z[:, None]*case["B"]
        rhs = np.c_[adjoint(geometry["S"]), p]
        x = solve(a.astype(complex, copy=False), rhs, assume_a="her", check_finite=False)
        rows = [verify_cohort(a, x, p, geometry, case, r, tx, mu, a_error, nz, a_gap, a_nu)
                for r, tx in exact_measurement_groups(case)[partition]]
        sum_sq = sum(r.get("lower_squared", 0.) for r in rows)
        lower = max(0., down(np.sqrt(max(0., down(sum_sq * (1.-gamma(len(rows)+32)))))))
        return dict(status="discrete_indefinite_outward", proposal=proposal, base=proof, groups=rows,
                    lower_outward=lower, empirical_empty_RMS_margin=lower-delta, **elapsed(started))
    except (ArithmeticError, LinAlgError) as exc:
        return dict(status="outward_unresolved", reason=str(exc), proposal=proposal, **elapsed(started))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequency", type=int, choices=FREQUENCIES, action="append")
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api="blas"):
        for frequency in args.frequency or FREQUENCIES:
            scan = json.loads((RESULTS / f"scan_{frequency}GHz.json").read_text())
            geometry = operators(frequency)
            cases = {r["name"]: load_case(r["name"], frequency, geometry, model_audit=False)
                     for r in scan["cases"]}
            result = verify_candidate(scan, cases, geometry)
            (RESULTS / f"indefinite_{frequency}GHz.json").write_text(json.dumps(
                result, indent=2, allow_nan=False)+"\n")
            print(frequency, result.get("status"), result.get("lower_outward"), flush=True)


if __name__ == "__main__":
    main()
