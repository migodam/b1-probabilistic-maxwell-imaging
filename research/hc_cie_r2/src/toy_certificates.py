#!/usr/bin/env python3
"""Independent interval and exact-rational checks for the frozen witnesses."""
import argparse
import json
from pathlib import Path
import platform
from time import perf_counter, process_time

import mpmath as mp
import numpy as np
import sympy as sp

from toy_replay import save


def soft_interval(source):
    mp.mp.dps = 110
    mp.iv.dps = 80
    iv = mp.iv
    problem = json.loads(source.read_text(), parse_float=str)

    def read_complex(key):
        pairs = problem[key]
        return iv.matrix([[iv.mpc(iv.mpf(x), iv.mpf(y)) for x, y in zip(row_real, row_imag)]
                          for row_real, row_imag in zip(pairs["real"], pairs["imag"])])

    def wedge(a, b):
        return [a[i, 0]*b[i, 1]-a[i, 1]*b[i, 0] for i in range(a.rows)]

    def flatten(a):
        return [a[i, j] for i in range(a.rows) for j in range(a.cols)]

    def realify(v):
        return [x.real for x in v]+[x.imag for x in v]

    def dot(v, w):
        return sum((x*y for x, y in zip(v, w)), iv.mpf(0))

    def bound_lower(t):
        return mp.mpf(t.a)

    def bound_upper(t):
        return mp.mpf(t.b)

    g, b, s, y, j, j_ref = (read_complex(k) for k in ("G", "B", "S", "Y", "J_bad", "J_true"))
    n, ns = j.rows, j.cols
    weight = iv.mpf(problem["weight"])
    dim = 2*n*ns
    e = b+g*j
    raw_wedge = wedge(j, e)

    def f(jj):
        return realify(flatten(s*jj-y)+[weight*v for v in wedge(jj, b+g*jj)])

    fj = f(j)
    fj_ref = f(j_ref)
    directions, green_directions, df = [], [], []
    for k in range(dim):
        h = iv.matrix(n, ns)
        h[(k//ns) % n, k % ns] = 1 if k < n*ns else iv.mpc(0, 1)
        gh = g*h
        d1, d2 = wedge(h, e), wedge(j, gh)
        df.append(realify(flatten(s*h)+[weight*(x+y) for x, y in zip(d1, d2)]))
        directions.append(h)
        green_directions.append(gh)
    gradient = [dot(c, fj) for c in df]
    hessian = iv.matrix(dim)
    for a in range(dim):
        for c in range(dim):
            d1 = wedge(directions[a], green_directions[c])
            d2 = wedge(directions[c], green_directions[a])
            residual_hessian = sum((iv.mpc(r.real, -r.imag)*(x+y)
                                    for r, x, y in zip(raw_wedge, d1, d2)), iv.mpc(0))
            hessian[a, c] = dot(df[a], df[c])+weight**2*residual_hessian.real
    nominal = np.array([[float(bound_lower(hessian[a, c])) for c in range(dim)]
                        for a in range(dim)])
    _, eigenvectors = np.linalg.eigh(nominal)
    change = iv.matrix([[iv.mpf(repr(float(x))) for x in row] for row in eigenvectors])
    rotated = change.T*hessian*change
    gram = change.T*change

    def gershgorin_lower(matrix):
        return min(bound_lower(matrix[a, a]-sum((abs(matrix[a, c]) for c in range(dim)
                                                 if c != a), iv.mpf(0)))
                   for a in range(dim))

    def rowsum_upper(matrix):
        return max(bound_upper(sum((abs(matrix[a, c]) for c in range(dim)), iv.mpf(0)))
                   for a in range(dim))

    rotated_lower = gershgorin_lower(rotated)
    gram_lower = gershgorin_lower(gram)
    gram_upper = rowsum_upper(gram)
    gradient2 = dot(gradient, gradient)
    df_norm2 = sum((dot(c, c) for c in df), iv.mpf(0))
    offdiag_norm2 = sum((abs(g[a, c])**2 for a in range(n) for c in range(n)
                        if a != c), iv.mpf(0))
    quadratic_bound2 = weight**2*offdiag_norm2
    radius = iv.mpf("1e-7")
    point_lower = iv.mpf("2.39e-6")
    kappa_upper = iv.mpf(".088")
    df_upper = iv.mpf("1.419")
    grad_upper = iv.mpf("6e-16")
    lip = 6*kappa_upper*(df_upper+2*kappa_upper*radius)
    reported_lip_upper = iv.mpf(".75")
    ball_lower = point_lower-reported_lip_upper*radius
    outward = ball_lower*radius-grad_upper
    f_center = dot(fj, fj)/2
    f_upper = bound_upper(dot(fj_ref, fj_ref)/2)
    # Interval comparisons use rational-decimal thresholds and outward rounding.
    checks = dict(
        gram_change_invertible=gram_lower > 0,
        center_positive=rotated_lower > bound_upper(point_lower*iv.mpf(gram_upper)),
        gradient_bounded=bound_upper(gradient2) < bound_lower(grad_upper**2),
        derivative_bounded=bound_upper(df_norm2) < bound_lower(df_upper**2),
        quadratic_remainder_bounded=bound_upper(quadratic_bound2) < bound_lower(kappa_upper**2),
        hessian_lipschitz_bounded=bound_upper(lip) < bound_lower(reported_lip_upper),
        radius_boundary_gradient_outward=bound_lower(outward) > 0,
        objective_in_ball_gt_4e_4=bound_lower(f_center-grad_upper*radius)>mp.mpf(".0004"),
        objective_near_truth_lt_1e_25=f_upper < mp.mpf("1e-25"),
    )
    assert all(checks.values()), checks
    return dict(
        scope="Exact decimal inputs only, no enclosure of the Hankel functions or continuous PDE.",
        interval_precision=mp.iv.dps, checks=checks,
        hessian_lower_on_ball="2.315e-6", radius="1e-7",
        outward_normalized_gradient_lower="2.309e-13",
        bad_ball_objective_lower="4e-4", near_truth_objective_upper="1e-25",
        comparison=dict(
            center_objective=str(mp.nstr(bound_lower(f_center), 28)),
            near_truth_objective_upper_display=str(mp.nstr(f_upper, 28)),
            rotated_hessian_gershgorin_lower_display=str(mp.nstr(rotated_lower, 28)),
            gram_upper_display=str(mp.nstr(gram_upper, 28)),
        ),
    )


def exact_tail(source_dir):
    R, I = sp.Rational, sp.I
    v = sp.Matrix([1, R("0.8"), R("0.6")])
    u = (R("0.15")+I*R("0.1"))*v
    obs = sp.Matrix([[1, 0, R("0.35")], [0, 1, R("-0.2")]])
    source_field = sp.Matrix([
        [1, R("0.4")+I*R("0.7")], [R("0.6")+I*R("0.1"), 1],
        [R("0.3")+I*R("0.8"), R("-0.5")+I*R("0.6")],
    ])
    bridge = obs.col_join(v.T)
    true_x = sp.diag(3, 8, 13)
    tau_max = R(15, 7)

    def real_inner(x, y):
        return sp.re(sp.expand_complex(sp.trace(x.conjugate().T*y)))

    def fixed_q_model(epsilon):
        tail = I*epsilon*sp.eye(3)
        green = -R("0.4")*sp.eye(3)+u*v.T+tail
        current = sp.simplify((sp.eye(3)-true_x*green).inv()*true_x*source_field)
        feedback = sp.simplify(v.T*current)
        return tail, sp.simplify(source_field+u*feedback), sp.simplify(bridge*current)

    small_eps = R("0.0001")
    _, drive, observed = fixed_q_model(small_eps)
    bases = [sp.diag(*[int(i == j) for j in range(3)]) for i in range(3)]
    lipschitz_columns = [bridge*e*drive for e in bases]
    gram = sp.Matrix(3, 3, lambda i, j: real_inner(lipschitz_columns[i],
                                                  lipschitz_columns[j]))
    gamma = R("0.79")
    shifted = gram-gamma**2*sp.eye(3)
    minors = [shifted[:k, :k].det() for k in (1, 2, 3)]
    alpha = tau_max*small_eps
    alpha_upper = R(1, 4000)
    a_f, m_norm, b_norm = 4, 3, 10
    perturb = m_norm*a_f*alpha_upper*(1/(1-alpha_upper)+1/(1-alpha_upper)**2)
    c2 = 2*m_norm*small_eps*a_f/(1-alpha_upper)**3
    rho = m_norm*tau_max*a_f/(1-alpha_upper)+b_norm
    curvature_lower = sp.factor((gamma-perturb)**2-rho*c2)
    checks = dict(
        gram_shift_positive=all(x > 0 for x in minors),
        bridge_frobenius_below_3=bool(real_inner(bridge, bridge) < m_norm**2),
        drive_frobenius_below_4=bool(real_inner(drive, drive) < a_f**2),
        target_frobenius_below_10=bool(real_inner(observed, observed) < b_norm**2),
        alpha_bound=bool(alpha < alpha_upper < 1),
        jacobian_margin=bool(gamma > perturb),
        curvature_bound=bool(curvature_lower > R(1, 2)),
    )

    witness_eps = R("0.1")
    tail, drive, observed = fixed_q_model(witness_eps)
    tau = sp.diag(R("2.1"), R("2.1"), R("0.05"))
    h = sp.diag(R("-0.5"), R("0.1"), R("0.86"))
    # Differentiate the rational solution directly along a symbolic line.
    t = sp.Symbol("t", real=True)
    state = (sp.eye(3)-(tau+t*h)*tail).inv()*(tau+t*h)*drive
    moving = bridge*state
    f0 = moving.subs(t, 0)-observed
    f1 = moving.diff(t).subs(t, 0)
    f2 = moving.diff(t, 2).subs(t, 0)
    curvature_negative = sp.factor(real_inner(f1, f1)+real_inner(f0, f2))
    checks.update(
        contraction_witness=bool(tau_max*witness_eps < 1),
        witness_point_inside_box=all(0 < v < tau_max for v in tau.diagonal()),
        witness_curvature_negative=bool(curvature_negative < 0),
    )
    originals = [json.loads((source_dir/name).read_text()) for name in (
        "variables_exact_convex_certificate.json", "variables_exact_nonconvex_witness.json")]
    original_match = [curvature_lower == R(originals[0]["strong_convexity_modulus_lower_exact"]),
                      curvature_negative == R(originals[1]["second_directional_derivative_exact"])]
    assert all(checks.values()) and all(original_match), (checks, original_match)
    return dict(
        scope="Reciprocal constructed matrix; q fixed at q_true independently for each tail.",
        checks=checks, exact_fraction_agrees_with_source=[bool(x) for x in original_match],
        positive=dict(alpha=str(alpha), lower_exact=str(curvature_lower),
                      lower_display=str(sp.N(curvature_lower, 30)),
                      shifted_real_gram_leading_minors=[str(x) for x in minors]),
        negative=dict(alpha=str(tau_max*witness_eps), curvature_exact=str(curvature_negative),
                      curvature_display=str(sp.N(curvature_negative, 30)),
                      point=["2.1", "2.1", "0.05"], direction=["-0.5", "0.1", "0.86"]),
    )


def exact_physical_ambiguity(source):
    original = json.loads(source.read_text())
    g = sp.Matrix(original["G"]).applyfunc(sp.sympify)
    b = sp.Matrix(original["B"]).applyfunc(sp.sympify)
    s = sp.Matrix([original["S"]]).applyfunc(sp.sympify)
    ident = sp.eye(3)
    plants = []
    currents, predictions = [], []
    for chi in (10, 20):
        operator = ident-chi*g
        current = sp.simplify(operator.inv()*chi*b)
        total = b+g*current
        inv = s*operator.inv()
        jac = sp.Matrix(2, 3, lambda i, j: sp.simplify(inv[j]*total[j, i]))
        real_jac = jac.applyfunc(sp.re).col_join(jac.applyfunc(sp.im))
        det = sp.factor((real_jac.T*real_jac).det())
        spectrum = (operator.conjugate().T*operator).eigenvals()
        plants.append(dict(chi=chi, jacobian_gram_determinant=str(det),
                           jacobian_full_rank=bool(det>0),
                           forward_condition=str(sp.sqrt(max(spectrum)/min(spectrum)))))
        currents.append(current)
        predictions.append(sp.simplify(s*current))
    hermitian_im = (g-g.conjugate().T)/(2*sp.I)
    energy_eigenvalues = (hermitian_im-s.conjugate().T*s).eigenvals()
    energy_eigenvalue_lower = min(energy_eigenvalues)
    assert predictions[0] == predictions[1] and currents[0] != currents[1]
    assert g == g.T and energy_eigenvalue_lower == sp.Rational(29, 320)
    assert all(x["jacobian_full_rank"] for x in plants)
    return dict(scope="Engineered linear receiver, exact reciprocal passive finite matrix.",
                two_positive_materials=plants, shared_data=[str(x) for x in predictions[0]],
                energy_gap=str(energy_eigenvalue_lower))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    directory = args.source/"research_work"
    t0, c0 = perf_counter(), process_time()
    interval = soft_interval(directory/"tm_disjoint_exact_integral_counterexample_verified.json")
    tail = exact_tail(directory)
    ambiguity = exact_physical_ambiguity(directory/"exact_two_material_ambiguity.json")
    save(args.out/"certificates.json", dict(
        interval=interval, tail=tail, finite_distance_ambiguity=ambiguity,
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         sympy=sp.__version__, mpmath=mp.__version__,
                         wall_seconds=perf_counter()-t0, cpu_seconds=process_time()-c0),
    ))
    print("interval, rational tail, and rational ambiguity checks passed", flush=True)


if __name__ == "__main__":
    main()
