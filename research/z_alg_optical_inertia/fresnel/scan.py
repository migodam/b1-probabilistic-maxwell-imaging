"""Fresnel 2005 masked multi-source sector/optical inertia scan."""
import argparse
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, eigvalsh, get_lapack_funcs
from threadpoolctl import threadpool_limits

from physical import (FREQUENCIES, NAMES, RESULTS, adjoint, clock, elapsed, hermitian,
                      load_case, operators)

CAPS = (.005, .02, .05, .1, .25, .5, 1., 2., 2.3)
# Common scaling from the physical receive operator.  Bound remains valid
# for each multiplier individually, including after their finite search.
MU_RELATIVE = 10. ** np.linspace(-8, 4, 17)
# Stated scenario ratios, not measurement noise bounds.
RADII = (.01, .03, .05, .11, .15)


@dataclass(frozen=True)
class Profile:
    family: str
    name: str
    weight_center: float = 1.
    weight_outer: float = 1.
    optical_center: float = 0.
    optical_outer: float = 0.

    def parameters(self, centers, cap):
        central = (np.max(np.abs(centers), axis=1) < .05)
        scale = np.where(central, self.weight_center, self.weight_outer)
        ratio = cap * np.where(central, self.optical_center, self.optical_outer)
        return scale, ratio


PROFILES = (
    Profile("real", "uniform_real"),
    Profile("optical", "uniform_positive", optical_center=1., optical_outer=1.),
    Profile("optical", "uniform_negative", optical_center=-1., optical_outer=-1.),
    Profile("optical", "uniform_positive_strong", optical_center=4., optical_outer=4.),
    Profile("two_zone", "two_zone_center", 1.5, .75, 1.5, .25),
    Profile("two_zone", "two_zone_outer", .75, 1.5, .25, 1.5),
)
PARTITIONS = {"individual": {8: 1, 18: 1},
              "adjacent_pairs": {8: 2, 18: 2},
              "adjacent_cohorts": {8: 4, 18: 6}}


def exact_measurement_groups(case):
    mask = case["mask"]
    count = mask.shape[1]
    groups = {}
    for key, sizes in PARTITIONS.items():
        step = sizes[count]
        groups[key] = []
        for start in range(0, count, step):
            tx = np.arange(start, start + step)
            rx = np.flatnonzero(mask[:, tx].all(axis=1))
            if not rx.size:
                raise ValueError("An acquisition cohort has no common receiver.")
            groups[key].append((rx, tx))
    return groups


def factor_inertia(factor, pivots):
    """Count signs of the 1×1 and 2×2 blocks of D in LAPACK's LDL*."""
    diagonal = factor.diagonal().real
    values = diagonal.copy()
    index = 0
    while index < len(pivots):
        if pivots[index] > 0:
            index += 1
        else:
            coupling = factor[index + 1, index]
            values[index:index + 2] = eigvalsh(
                np.array([[diagonal[index], coupling.conj()],
                          [coupling, diagonal[index + 1]]]), check_finite=False)
            index += 2
    return int((values < 0).sum()), float(np.min(np.abs(values)))


def base_matrix(geometry, cap, profile):
    weight, ratio = profile.parameters(geometry["centers"], cap)
    z = (cap * weight + 1j * ratio) * .5
    g = geometry["G"]
    if np.ptp(weight) == 0 and np.ptp(ratio) == 0:
        a = np.array(-cap * weight[0] * g.real + ratio[0] * g.imag, order="F")
    else:
        a = np.asfortranarray(-2. * hermitian(z[:, None] * g))
    a.flat[::len(a) + 1] += weight
    return a, z, weight


def receiver_reduction(geometry, observations, cap, profile, cost):
    """Factor A exactly once for all receiver masks and all three measurements."""
    started = clock()
    a, z, _ = base_matrix(geometry, cap, profile)
    s = geometry["S"]
    b = np.concatenate([case["B"] for case in observations], axis=1)
    p = z[:, None] * b
    routine = ("sytrf", "sytrs") if np.isrealobj(a) else ("hetrf", "hetrs")
    get_factor, get_solve = get_lapack_funcs(routine, (a,))
    factor, pivots, info = get_factor(a, lower=1, lwork=64*len(a), overwrite_a=0)
    cost["A_LDL"] += 1
    if info:
        raise ArithmeticError(f"Singular multiplier base, LAPACK LDL info={info}")
    nu_base, d_gap = factor_inertia(factor, pivots)
    rhs = np.asfortranarray(np.c_[adjoint(s), p])
    if np.isrealobj(factor):
        planes = np.asfortranarray(np.c_[rhs.real, rhs.imag])
        solved, info = get_solve(factor, pivots, planes, lower=1, overwrite_b=1)
        half = solved.shape[1] // 2
        solved = solved[:, :half] + 1j * solved[:, half:]
        cost["A_real_rhs"] += planes.shape[1]
    else:
        solved, info = get_solve(factor, pivots, rhs, lower=1, overwrite_b=0)
        cost["A_complex_rhs"] += rhs.shape[1]
    if info:
        raise ArithmeticError(f"Base inverse did not complete, LAPACK solve info={info}")
    defect = np.linalg.norm(a @ solved - rhs) / np.linalg.norm(rhs)
    if defect > 1.e-7:
        raise ArithmeticError(f"Base inverse relative backward residual {defect}")
    r = s.shape[0]
    partitioned = np.cumsum([case["B"].shape[1] for case in observations])
    return (hermitian(s @ solved[:, :r]), np.split(s @ solved[:, r:], partitioned[:-1], axis=1),
            np.split(p, partitioned[:-1], axis=1), np.split(solved[:, r:], partitioned[:-1], axis=1),
            dict(A_negative=nu_base, D_min_abs_pivot=d_gap,
                 A_inverse_relative_defect=float(defect), **elapsed(started)))


def cohort_curves(case, reduction, predicted, p, x, nu_base, mus, cache, cost):
    y = case["Y"]
    curves = {}
    for label, groups in case["groups"].items():
        curves[label] = []
        for rx, tx in groups:
            key = tuple(rx)
            if key not in cache:
                cache[key] = eigh(reduction[np.ix_(rx, rx)], check_finite=False)
                cost["receiver_Schur_eigensystems"] += 1
            d, v = cache[key]
            pair = np.ix_(rx, tx)
            dport = adjoint(v) @ (y[pair] - predicted[pair])
            e = hermitian(adjoint(p[:, tx]) @ x[:, tx])
            rows = []
            for mu in mus:
                shift = d + mu
                tol = max(1.e-15, np.finfo(float).eps * 128 * len(rx) * np.max(abs(shift)))
                neg = nu_base - int(np.count_nonzero(shift < 0))
                row = dict(nu=neg, receiver_update_min_abs=float(np.min(abs(shift))))
                if row["receiver_update_min_abs"] <= tol or not (0 <= neg <= p.shape[0]):
                    row.update(status="numeric_unresolved", lower_bound=None, H_eigs=None)
                elif neg >= len(tx):
                    row.update(status="sources_exhausted", lower_bound=0., H_eigs=None)
                else:
                    h = mu * (adjoint(dport) @ (dport / shift[:, None]) - e)
                    eig = eigvalsh(hermitian(h), check_finite=False)[::-1]
                    cost["H_eigensystems"] += 1
                    row.update(status="pd" if neg == 0 else "indefinite",
                               lower_bound=float(np.sqrt(np.maximum(eig[neg:], 0.).sum())),
                               H_eigs=eig.tolist())
                rows.append(row)
            curves[label].append(dict(tx=(tx+1).tolist(), rx=(rx+1).tolist(),
                channels=int(len(tx)*len(rx)), norm=float(np.linalg.norm(y[pair])), rows=rows))
    return curves


def choose(curves, mus):
    options = []
    for partition, groups in curves.items():
        for j, mu in enumerate(mus):
            together = [g["rows"][j] for g in groups]
            isolated = [r["lower_bound"] if r["lower_bound"] is not None else 0. for r in together]
            options.append(dict(mu=mu, partition=partition, mu_index=j,
                                numerical_bound=float(np.linalg.norm(isolated)),
                                inertia=[r["nu"] for r in together],
                                group_statuses=[r["status"] for r in together],
                                retained_channels=sum(g["channels"] for g in groups),
                                informative_channels=sum(g["channels"] for g, bound in zip(groups, isolated)
                                                         if bound > 0.),
                                informative_indefinite_groups=sum(r["status"] == "indefinite"
                                    and bool(r["lower_bound"]) for r in together)))
    best = max(options, key=lambda row: row["numerical_bound"])
    best_indefinite = max((row for row in options if row["informative_indefinite_groups"]),
                          key=lambda row: row["numerical_bound"], default=None)
    return best, best_indefinite


def reductions(frequency, output, *, model_audit):
    start = clock()
    geometry = operators(frequency)
    geometry_cost = elapsed(start)
    observations = [load_case(name, frequency, geometry, model_audit=model_audit) for name in NAMES]
    s = geometry["S"]
    mu_scale = np.linalg.norm(s)**2 / s.shape[0]
    mus = (mu_scale * MU_RELATIVE).tolist()
    for item in observations:
        item["groups"] = exact_measurement_groups(item)
        record = item["record"]
        record["partitions"] = {key: [dict(source_ids_one_based=(tx+1).tolist(),
                                        receiver_ids_one_based=(rx+1).tolist(), channels=len(tx)*len(rx))
                                       for rx, tx in groups] for key, groups in item["groups"].items()}
        # The strictly air-only hypothesis has j=0. Outward evaluation of its
        # normalization is computed in verify.py using the represented Y.
        record["air_only_normalized"] = 1.
        record["grid"] = []
    cost, base_records = Counter(), []
    for cap in CAPS:
        for profile in PROFILES:
            print(f"Fresnel {frequency} GHz / cap {cap:g} / {profile.name}", flush=True)
            c, predictions, p, x, base = receiver_reduction(geometry, observations, cap, profile, cost)
            base_records.append(dict(u=cap, profile=profile.name, **base))
            cache = {}
            for item, prediction, pi, xi in zip(observations, predictions, p, x, strict=True):
                curves = cohort_curves(item, c, prediction, pi, xi, base["A_negative"], mus, cache, cost)
                best, best_indefinite = choose(curves, mus)
                item["record"]["grid"].append(dict(u=cap, family=profile.family, profile=profile.name,
                                                   best=best, best_indefinite=best_indefinite, curves=curves))
    records = [item["record"] for item in observations]
    result = dict(frequency_GHz=frequency, side=60, width_m=.20,
                  multiplier_profiles=[vars(profile) for profile in PROFILES], mus=mus,
                  radii_scenarios=RADII, costs=dict(cost), geometry_cost=geometry_cost,
                  base_records=base_records, cases=records, **elapsed(start))
    (output/f"scan_{frequency}GHz.json").write_text(json.dumps(result, allow_nan=False))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS)
    parser.add_argument("--frequency", type=int, choices=FREQUENCIES, action="append")
    parser.add_argument("--without-nominal-forward", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1, user_api="blas"):
        for freq in args.frequency or FREQUENCIES:
            result = reductions(freq, args.output, model_audit=not args.without_nominal_forward)
            print(f"Fresnel {freq} GHz done: {result['wall_seconds']:.2f}s", flush=True)


if __name__ == "__main__":
    main()
