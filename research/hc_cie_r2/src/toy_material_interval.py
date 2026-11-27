#!/usr/bin/env python3
"""Certify strict local minima for the two full-forward material objectives."""
import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import mpmath as mp
import numpy as np

from toy_replay import save


iv = mp.iv
CENTERS = {
    "hard5": ["7.7008094099423445", "20.077386285156923", "8.231818051737394",
              "7.601310290275903", "11.665529400535293"],
    "soft3": ["10.326565116600996", "0.0", "9.959316299438955"],
}
# These conservative decimal outputs must be strictly enclosed by the computed bounds.
CLAIMS = {
    "hard5": ("2.08e-7", "2.07e-15", "3.4246e-5", "1e-30"),
    "soft3": ("8.26e-5", "8.26e-13", "0.0042476", "1e-30"),
}


def lo(x):
    return mp.mpf(x.a)


def hi(x):
    return mp.mpf(x.b)


def view(x):
    return mp.nstr(x, 28)


def exact_decimal_matrix(record):
    return iv.matrix([[iv.mpc(iv.mpf(re), iv.mpf(im)) for re, im in zip(rr, ii)]
                      for rr, ii in zip(record["real"], record["imag"])])


def verified_inverse(operator):
    """Outward-rounded Gaussian elimination, with pivots bounded away from zero."""
    n = operator.rows
    matrix = [[operator[i, j] for j in range(n)]+[iv.mpc(int(i == j)) for j in range(n)]
              for i in range(n)]
    min_pivot = mp.inf
    for i in range(n):
        index = max(range(i, n), key=lambda j: lo(abs(matrix[j][i])))
        matrix[index], matrix[i] = matrix[i], matrix[index]
        pivot = matrix[i][i]
        bound = lo(abs(pivot))
        assert bound > 0, "interval pivot contains zero"
        min_pivot = min(min_pivot, bound)
        for j in range(i+1, 2*n):
            matrix[i][j] /= pivot
        matrix[i][i] = iv.mpc(1)
        for k in range(n):
            if k == i:
                continue
            multiplicand = matrix[k][i]
            for j in range(i+1, 2*n):
                matrix[k][j] -= multiplicand*matrix[i][j]
            matrix[k][i] = iv.mpc(0)
    return iv.matrix([row[n:] for row in matrix]), min_pivot


def hankel_zero_one(z):
    """Enclose H_0 and H_1 via the convergent logarithmic Bessel series on 0<z<20."""
    assert 0 < lo(z) and hi(z) < 20
    half = z/2
    a = half**2
    n_terms = 120
    assert hi(4*a/(n_terms+2)**2)<mp.mpf(".5")
    base, harmonic = iv.mpf(1), iv.mpf(0)
    sum_zero, sum_one, log_extra, derive_extra = iv.mpf(1), iv.mpf(1), iv.mpf(0), iv.mpf(0)
    for k in range(1, n_terms+2):
        base = -base*a/k**2
        harmonic += iv.mpf(1)/k
        log_term = -base*harmonic
        derivative_term = 2*k*log_term/z
        one_term = base/(k+1)
        if k <= n_terms:
            sum_zero += base
            sum_one += one_term
            log_extra += log_term
            derive_extra += derivative_term

    def enclosure(first_missing):
        # Ratios after the first missing term are bounded by 4*(z²/4)/(k+1)²
        # even for the differentiated harmonic series. All are less than 1/2.
        upper = hi(2*abs(first_missing))
        return iv.mpf([-upper, upper])

    j0 = sum_zero+enclosure(base)
    j1 = half*(sum_one+enclosure(one_term))
    common = iv.ln(half)+iv.euler
    y0 = 2/iv.pi*(common*j0+log_extra+enclosure(log_term))
    y1 = 2/iv.pi*(common*j1-j0/z-derive_extra-enclosure(derivative_term))
    return iv.mpc(j0, y0), iv.mpc(j1, y1), j1


def geometry(record):
    n = len(record["positions"])
    points = [[iv.mpf(coord) for coord in row] for row in record["positions"]]
    angles = [iv.mpf(t) for t in record["receiver_angles"]]
    direct = [iv.mpf(t) for t in record["incident_angles"]]
    incident = [[iv.cos(t), iv.sin(t)] for t in direct]
    receiver = [[2*iv.cos(t), 2*iv.sin(t)] for t in angles]
    k = 2*iv.pi
    radius = iv.sqrt(iv.mpf(record["area"])/iv.pi)

    def norm2(a, b):
        return iv.sqrt(sum(((x-y)**2 for x, y in zip(a, b)), iv.mpf(0)))

    _, h1, j1 = hankel_zero_one(k*radius)
    prefactor = iv.pi*k*radius/2
    disk = iv.mpc(0, prefactor*j1)
    self_term = iv.mpc(0, prefactor)*h1-1
    g = iv.matrix(n)
    for i in range(n):
        g[i, i] = self_term
        for j in range(i):
            g[i, j] = g[j, i] = disk*hankel_zero_one(k*norm2(points[i], points[j]))[0]
    s = iv.matrix([[disk*hankel_zero_one(k*norm2(r, p))[0] for p in points] for r in receiver])
    b = iv.matrix([[iv.exp(iv.mpc(0, k*sum((u*v for u, v in zip(r, p)), iv.mpf(0))))
                    for r in incident] for p in points])
    chi_truth = [iv.mpf(v) for v in record["chi"]]
    inverse, _ = verified_inverse(iv.eye(n)-iv.diag(chi_truth)*g)
    y = s*inverse*iv.diag(chi_truth)*b
    return g, s, b, y, chi_truth


def real_inner(a, b):
    return sum((a[i,j].real*b[i,j].real+a[i,j].imag*b[i,j].imag
                for i in range(a.rows) for j in range(a.cols)), iv.mpf(0))


def state(chi, g, s, b, y, derivative=True):
    n = len(chi)
    x = iv.diag(chi)
    inverse, pivot_lower = verified_inverse(iv.eye(n)-x*g)
    j = inverse*x*b
    f = s*j-y
    objective = real_inner(f, f)/2
    if not derivative:
        return objective

    e = b+g*j
    basis = [iv.diag([int(i == j) for j in range(n)]) for i in range(n)]
    d_j = [inverse*h*e for h in basis]
    d_gj = [g*t for t in d_j]
    df = [s*t for t in d_j]
    grad = [real_inner(f, t) for t in df]
    hessian = iv.matrix(n)
    for i in range(n):
        for j in range(i+1):
            d2 = s*inverse*(basis[i]*d_gj[j]+basis[j]*d_gj[i])
            v = real_inner(df[i], df[j])+real_inner(f, d2)
            hessian[i,j] = hessian[j,i] = v
    return objective, grad, hessian, pivot_lower


def hessian_lower(interval_hessian, center_hessian, free):
    naive = np.array([[float(lo(center_hessian[i,j])) for j in free] for i in free])
    _, eigenvectors = np.linalg.eigh(naive)
    correction = iv.matrix([[iv.mpf(repr(float(x))) for x in row] for row in eigenvectors])
    restricted_hessian = iv.matrix([[interval_hessian[i,j] for j in free] for i in free])
    rotated = correction.T*restricted_hessian*correction
    gram = correction.T*correction

    def gershgorin_lower(matrix):
        return min(lo(matrix[i,i]-sum((abs(matrix[i,j]) for j in range(matrix.cols) if j != i),
                                      iv.mpf(0))) for i in range(matrix.rows))

    gram_upper = max(hi(sum((abs(gram[i,j]) for j in range(gram.cols)), iv.mpf(0)))
                     for i in range(gram.rows))
    numerator = gershgorin_lower(rotated)
    assert gershgorin_lower(gram)>0 and numerator>0, "no positive Hessian enclosure"
    return lo(iv.mpf(numerator)/iv.mpf(gram_upper))


def certify(record, case, use_geometry):
    center_values = CENTERS[case]
    if use_geometry:
        operators = geometry(record)
    else:
        operators = (*[exact_decimal_matrix(record[k]) for k in ("G", "S", "B", "Y")],
                     [iv.mpf(v) for v in record["chi"]])
    g, s, b, y, true_chi = operators
    chi = [iv.mpf(v) for v in center_values]
    radius = iv.mpf("1e-8")
    free = [i for i, v in enumerate(chi) if lo(v) != 0]
    active = [i for i in range(len(chi)) if i not in free]
    assert all(lo(chi[i]-radius)>0 and hi(chi[i]+radius)<30 for i in free)
    full_box = [chi[i]+iv.mpf([-hi(radius), hi(radius)]) if i in free else iv.mpf([0, hi(radius)])
                for i in range(len(chi))]
    f0, center_grad, center_hessian, center_pivot = state(chi, g, s, b, y)
    _, box_grad, box_hessian, box_pivot = state(full_box, g, s, b, y)
    lower = hessian_lower(box_hessian, center_hessian, free)
    center_gradient_upper = hi(iv.sqrt(sum((center_grad[i]**2 for i in free), iv.mpf(0))))
    margin = lo(iv.mpf(lower)*radius-iv.mpf(center_gradient_upper))
    positive_active = [lo(box_grad[i]) for i in active]
    f_bad_lower = lo(f0-iv.mpf(center_gradient_upper)*radius)
    f_true_upper = hi(state(true_chi, g, s, b, y, derivative=False))
    assert margin>0 and f_bad_lower>f_true_upper and all(t>0 for t in positive_active), (
        "no strict non-global local minimum certificate", view(margin), view(f_bad_lower),
        view(f_true_upper), [view(x) for x in positive_active]
    )
    hessian_claim, radial_claim, bad_claim, truth_claim = CLAIMS[case]
    assert lo(iv.mpf(lower)-iv.mpf(hessian_claim))>0
    assert lo(iv.mpf(margin)-iv.mpf(radial_claim))>0
    assert lo(iv.mpf(f_bad_lower)-iv.mpf(bad_claim))>0
    assert lo(iv.mpf(truth_claim)-iv.mpf(f_true_upper))>0
    active_claim = "0.00272999" if positive_active else None
    assert all(lo(iv.mpf(t)-iv.mpf(active_claim))>0 for t in positive_active)
    return dict(
        certified_conservative_bounds=dict(
            hessian_free_lower=hessian_claim, outward_radial_derivative_lower=radial_claim,
            bad_ball_objective_lower=bad_claim, truth_objective_upper=truth_claim,
            active_gradient_lower=active_claim,
        ),
        domain=("Exact decimal geometry parameters; k=2π; Bessel/Hankel power-series tail enclosed; "
                "data generated with these same enclosed operators.") if use_geometry else
               "Exact decimal stored G,S,B,Y finite matrix problem.",
        radius="1e-8", center=[str(v) for v in center_values],
        gradient_center_upper=view(center_gradient_upper),
        hessian_free_lower=view(lower), outer_radial_derivative_lower=view(margin),
        positive_active_gradients_lower=[view(x) for x in positive_active],
        free_cells_1based=[i+1 for i in free], active_lower_cells_1based=[i+1 for i in active],
        bad_ball_objective_lower=view(f_bad_lower), truth_objective_upper=view(f_true_upper),
        center_pivot_lower=view(center_pivot), box_pivot_lower=view(box_pivot),
        strict_non_global_minimum_certified=True,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    t0, c0 = perf_counter(), process_time()
    mp.mp.dps = 115
    iv.dps = 80
    directory = args.source/"research_work"
    report = dict(interval_decimal_precision=iv.dps,
                  certificate="Positive Hessian on the free ball and strictly outward boundary "
                              "gradient locate a unique interior minimum on its face. A strict "
                              "positive lower bound for the active gradient extends the conclusion "
                              "to the physical box. Smaller objective at the truth excludes a global minimum.")
    for case, input_name in (
        ("hard5", "tm_hardfiber_counterexample.json"),
        ("soft3", "tm_disjoint_exact_integral_counterexample_verified.json"),
    ):
        record = json.loads((directory/input_name).read_text(), parse_float=str)
        report[case] = {kind: certify(record, case, use_geometry)
                        for kind, use_geometry in [("decimal", False), ("geometry", True)]}
        print(case, report[case]["geometry"]["hessian_free_lower"], flush=True)
    report.update(wall_seconds=perf_counter()-t0, cpu_seconds=process_time()-c0)
    args.out.mkdir(parents=True, exist_ok=True)
    save(args.out/"material_interval.json", report)


if __name__ == "__main__":
    main()
