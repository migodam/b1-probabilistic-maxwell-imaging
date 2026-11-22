"""Conservative IEEE-754 outward bounds on selected discrete Fresnel certificates.

The represented binary64 G, S, B, and masked Y are the exact discrete
inputs.  Special-function accuracy and continuous-model error are outside
this numerical enclosure.  The scanner is used only to propose candidates.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import (LinAlgError, cho_factor, cho_solve, cholesky, eigh,
                          solve)
from threadpoolctl import threadpool_limits

from physical import (FREQUENCIES, RESULTS, adjoint, clock, elapsed, hermitian,
                      load_case, operators)
from scan import PROFILES, base_matrix, exact_measurement_groups

EPS = np.finfo(float).eps


def up(value):
    return float(np.nextafter(value, np.inf))


def down(value):
    return float(np.nextafter(value, -np.inf))


def gamma(ops):
    """Roundoff bound, using EPS (twice the unit roundoff) conservatively."""
    error = up(float(ops) * EPS)
    if error >= .01:
        raise ArithmeticError("Rounding bound leaves supported dimension.")
    return up(error / down(1. - error))


def inflate(sum_or_product):
    return up(sum_or_product * (1. + gamma(128)))


def norm_upper(a):
    """Sequential accumulation bound also covers pairwise numpy accumulation."""
    size = a.size
    magnitudes = np.sum(a.real * a.real)
    if np.iscomplexobj(a):
        magnitudes += np.sum(a.imag * a.imag)
    value = np.sqrt(magnitudes)
    if not np.isfinite(value):
        raise ArithmeticError("Nonfinite numerical enclosure.")
    return up(value * (1. + gamma(2 * size + 16)))


def norm_lower(a):
    a_eps = gamma(2 * a.size + 16)
    return max(0., down(norm_upper(a) * (1. - inflate(4. * a_eps))))


def gemm_error(a, b, na=None, nb=None):
    na = norm_upper(a) if na is None else na
    nb = norm_upper(b) if nb is None else nb
    # At most six real primitive operations per complex dot-product term.
    # Two additional operations per term cover accumulation and regrouping.
    return inflate(gamma(8 * a.shape[1] + 16) * na * nb)


def base_assembly_error(geometry, weight, z):
    ng = norm_upper(geometry["G"])
    nw = norm_upper(weight)
    nz = inflate(np.max(np.abs(z)))
    return inflate(gamma(64) * (nw + 2 * nz * ng)), nz


def cholesky_gap(a, error):
    """Prove A ≥ positive floor using a shifted approximate Cholesky."""
    started = clock()
    component_norm = np.abs(a.real)
    if np.iscomplexobj(a):
        component_norm += np.abs(a.imag)
    np.fill_diagonal(component_norm, 0.)
    row_sums = np.sum(component_norm, axis=1)
    diagonal = a.diagonal().real
    nominal = float(np.min(diagonal - row_sums))
    row_magnitude = float(np.max(np.abs(diagonal) + row_sums))
    row_error = inflate(gamma(4 * len(a) + 32) * row_magnitude)
    floor = down(nominal - inflate(row_error + error))
    if floor > .05:
        return floor, dict(method="outward_gershgorin", gap_lower=floor, **elapsed(started))
    a_norm = norm_upper(a)
    for shift in (.125, .02, .002, .0002):
        proposal = a.copy(order="F")
        proposal.flat[::len(a)+1] -= shift
        try:
            factor = np.tril(cholesky(proposal, lower=True, check_finite=False, overwrite_a=True))
        except LinAlgError:
            continue
        f_norm = norm_upper(factor)
        computed = factor @ adjoint(factor)
        computed -= a
        computed.flat[::len(a)+1] += shift
        defect = norm_upper(computed)
        evaluated = inflate(defect + error + gemm_error(factor, factor, f_norm, f_norm)
                             + gamma(24) * (a_norm + f_norm**2 + np.sqrt(len(a))*shift))
        if down(shift-evaluated) > 0:
            floor = down(shift-evaluated)
            return floor, dict(method="outward_shifted_cholesky", gap_lower=floor,
                               shift=shift, defect_upper=evaluated, **elapsed(started))
    raise ArithmeticError("Positive base gap unresolved.")


def rayleigh_bounds(h, error, skip=0):
    """Numerator Gershgorin and Gram denominator bound for each trial span."""
    eig, vec = eigh(hermitian(h), check_finite=False)
    eig, vec = eig[::-1], vec[:, ::-1]
    hn = norm_upper(h)
    terms = []
    for j in range(skip+1, len(eig)+1):
        if eig[j-1] <= 0:
            break
        q = vec[:, :j]
        qn = norm_upper(q)
        hq = h @ q
        projected = adjoint(q) @ hq
        err = inflate(qn * gemm_error(h, q, hn, qn)
                      + gemm_error(adjoint(q), hq, qn))
        projected_abs = np.abs(projected.real) + np.abs(projected.imag)
        np.fill_diagonal(projected_abs, 0.)
        projection_row = np.sum(projected_abs, axis=1)
        center = projected.diagonal().real
        projection_floor = down(float(np.min(center - projection_row)) - inflate(
            err + gamma(8*j+32) * float(np.max(np.abs(center) + projection_row))))
        gram = adjoint(q) @ q
        gram_rows = np.sum(np.abs(gram.real) + np.abs(gram.imag), axis=1)
        denominator = inflate(float(np.max(gram_rows)) * (1. + gamma(8*j+32))
                               + gemm_error(adjoint(q), q, qn, qn))
        term = down(down(projection_floor / denominator) - error) if projection_floor > 0 else 0.
        terms.append(max(0., term))
    return eig, terms


def verified_form(s, y, a, p, x, mu, ea, ep, gap):
    """Enclose H through the completed-square/linear-system residual identity."""
    ns, na, np_, nx = map(norm_upper, (s, a, p, x))
    obs = s @ x - y
    no = norm_upper(obs)
    err_obs = inflate(gemm_error(s, x, ns, nx) + gamma(24) * (ns*nx + norm_upper(y)))
    ax = a @ x
    nax = norm_upper(ax)
    err_ax = inflate(gemm_error(a, x, na, nx) + ea * nx)
    defect_numerical = adjoint(s) @ obs + mu * (ax-p)
    defect = inflate(norm_upper(defect_numerical) + ns*err_obs + mu*(err_ax + ep)
                     + gemm_error(adjoint(s), obs, ns, no)
                     + gamma(32) * (ns * no + mu * (nax + np_)))
    cross = adjoint(x) @ p
    h0_raw = adjoint(obs) @ obs + mu * (adjoint(x) @ ax - cross - adjoint(cross))
    h0 = hermitian(h0_raw)
    error_h0 = inflate(err_obs * (2 * no + err_obs) + mu * nx * (err_ax + 2 * ep)
                       + gemm_error(adjoint(obs), obs, no, no)
                       + mu * (gemm_error(adjoint(x), ax, nx, nax)
                               + 2 * gemm_error(adjoint(x), p, nx, np_))
                       + gamma(48) * (no**2 + mu * nx * (nax + 2 * np_)))
    error = inflate(error_h0 + inflate(defect**2 / gap))
    eig, certified = rayleigh_bounds(h0, error)
    return dict(H_truncated_eigenvalue_floors=certified, H0_numeric_eigenvalues=eig.tolist(),
                control_point_form_error_upper=error_h0,
                stationary_defect_upper=defect, H_total_error_upper=error,
                lower_squared=max(0., down(sum(certified) * (1. - gamma(len(certified) + 32)))))


def prove(case, geometry, profile, cap, proposal, cache):
    started = clock()
    if (profile, cap) not in cache:
        a, z, weight = base_matrix(geometry, cap, profile)
        formation_error, nz = base_assembly_error(geometry, weight, z)
        gap, gap_record = cholesky_gap(a, formation_error)
        factor = cho_factor(a, lower=True, check_finite=False)
        rhs = adjoint(geometry["S"])
        if np.isrealobj(a):
            parts = cho_solve(factor, np.c_[rhs.real, rhs.imag], check_finite=False)
            m = rhs.shape[1]
            s_inverted = parts[:, :m] + 1j * parts[:, m:]
        else:
            s_inverted = cho_solve(factor, rhs, check_finite=False)
        # Limit retention to one factor profile and frequency at a time.
        cache.clear()
        cache[(profile, cap)] = (a, z, formation_error, nz, gap, gap_record, factor, s_inverted)
    a, z, ea, nz, gap, gap_record, factor, s_inverted = cache[(profile, cap)]
    p = z[:, None] * case["B"]
    ep = inflate(gamma(64) * nz * norm_upper(case["B"]))
    if np.isrealobj(a):
        parts = cho_solve(factor, np.c_[p.real, p.imag], check_finite=False)
        m = p.shape[1]
        p_inverted = parts[:, :m] + 1j * parts[:, m:]
    else:
        p_inverted = cho_solve(factor, p, check_finite=False)
    mu = proposal["mu"]
    rows = []
    for rx, tx in exact_measurement_groups(case)[proposal["partition"]]:
        s = geometry["S"][rx]
        y = case["Y"][np.ix_(rx, tx)]
        x0 = p_inverted[:, tx]
        rec = s_inverted[:, rx]
        c = hermitian(s @ rec)
        c.flat[::len(c)+1] += mu
        q = solve(c, y - s @ x0, assume_a="pos", check_finite=False)
        x = x0 + rec @ q
        verified = verified_form(s, y, a, p[:, tx], x, mu, ea, ep, down(mu*gap))
        rows.append(dict(retained_channels=y.size, tx_one_based=(tx+1).tolist(),
                         verified_K_negative=0, verified_K_gap_lower=down(mu*gap), **verified))
    lower_sq = max(0., down(sum(row["lower_squared"] for row in rows)*(1.-gamma(len(rows)+16))))
    r_cap = max(0., down(np.sqrt(lower_sq)))
    return dict(status="strict_discrete_positive_semidefinite", u=cap, profile=profile.name,
                family=profile.family, proposal=proposal, A=gap_record, A_assembly_error_upper=ea,
                retained_channels=sum(row["retained_channels"] for row in rows),
                certificate_distance_lower=r_cap, groups=rows, **elapsed(started))


def select_proposals(scan, record):
    bases = {(base["u"], base["profile"]): base["A_negative"] for base in scan["base_records"]}
    grid = record["grid"]
    # This empirical threshold is a sensitivity scenario.
    target = max(.03, record['empty_calibration']['rms_transferred_over_all_scatter'])
    choices = []
    for family in ("real", "optical", "two_zone"):
        candidates = [row for row in grid if row["family"] == family
                      and bases[(row["u"], row["profile"])] == 0]
        largest = max((row for row in candidates if row["best"]["numerical_bound"] > target),
                      key=lambda row: (row["u"], row["best"]["numerical_bound"]), default=None)
        if largest is not None:
            choices.append(largest)
    if not choices:
        choices = [max(grid, key=lambda row: row["best"]["numerical_bound"])]
    return choices


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS)
    parser.add_argument("--frequency", type=int, choices=FREQUENCIES, action="append")
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api="blas"):
        for frequency in args.frequency or FREQUENCIES:
            started = clock()
            scan = json.loads((args.output / f"scan_{frequency}GHz.json").read_text())
            geometry = operators(frequency)
            cases = {}
            profiles = {profile.name: profile for profile in PROFILES}
            for row in scan["cases"]:
                case = load_case(row["name"], frequency, geometry, model_audit=False)
                air = norm_lower(case["Y"][case["mask"]])
                results = []
                for candidate in select_proposals(scan, row):
                    profile = profiles[candidate["profile"]]
                    try:
                        verified = prove(case, geometry, profile, candidate["u"], candidate["best"], {})
                    except (ArithmeticError, LinAlgError) as exc:
                        verified = dict(status="outward_unresolved", reason=str(exc), u=candidate["u"],
                                        profile=profile.name, family=profile.family)
                    results.append(verified)
                    print(row["name"], frequency, verified["u"], verified["profile"],
                          verified.get("certificate_distance_lower", verified.get("reason")), flush=True)
                cases[row["name"]] = dict(air_only_discrete_distance_outward_lower=air, results=results)
            (args.output / f"verified_{frequency}GHz.json").write_text(json.dumps(
                dict(cases=cases, numeric_contract="IEEE-754 fp64, no overflow or subnormal matrix terms",
                     **elapsed(started)), indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
