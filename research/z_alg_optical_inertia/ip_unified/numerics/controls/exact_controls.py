"""Exact algebraic controls for the fixed-Green lower bounds, Pro §§2–4."""
import json
from pathlib import Path
from time import perf_counter, process_time

from sympy import I, Matrix, Poly, Rational as Q, diag, expand, eye, simplify, sqrt, symbols, zeros

HERE = Path(__file__).resolve().parent


def hermitian_real(a):
    return (a + a.H)/2


def hermitian_imag(a):
    return (a - a.H)/(2*I)


def trace_norm_squared(a):
    return simplify((a.H*a).trace())


def source_spectral_controls():
    h = diag(10, 1)
    direction_budget = diag(0, 100)
    # The exact dual certificate coincides with an attainable primal value.
    dual_p = diag(1, 0)
    dual_a = dual_p
    optimum = (dual_p*h-dual_a*direction_budget).trace()
    assert optimum == 10 and min(h.eigenvals()) == 1
    q = diag(0, 1)
    d = diag(10, 0)
    assert min((direction_budget-q).eigenvals()) >= 0 and (d-h+q).is_zero_matrix
    # Sharp joint scalar caps with h=(8,4,2), q bounded by (5,3,0), tr q<=6.
    caps = (5, 3, 0)
    beta = Q(6)
    q_diag = tuple(beta/sum(caps)*cap for cap in caps)
    difference = diag(8, 4, 2)-diag(*q_diag)
    assert difference.trace() == 8 and min(difference.eigenvals()) > 0
    return dict(direction_h=str(h), directional_budget=str(direction_budget),
                rank_only_squared_bound=1, direction_primal_squared_value=10,
                direction_dual_squared_value=10, scalar_h=[8, 4, 2],
                scalar_eigenvalue_caps=list(caps), scalar_trace_cap=6,
                scalar_optimal_deducted_eigenvalues=[str(v) for v in q_diag],
                scalar_sharp_squared_distance=8)


def nonconvex_with_no_negative_inertia():
    g = I
    j_one = simplify(1/(1-g))
    midpoint = j_one/2
    implied_material = simplify(midpoint/(1+g*midpoint))
    assert implied_material == Q(2, 5)+I/5
    return dict(green='i', sensor=1, drive=1, state_matrix=1, negative_inertia=0,
                legal_box=[0, 1], physical_endpoint_current=str(j_one),
                current_midpoint=str(midpoint), implied_material=str(implied_material),
                physical_current_set_convex=False)


def compensation():
    g = Matrix([[2 + I/2, Q(1, 4)], [Q(1, 4), I/2]])
    s = Matrix([[sqrt(3)+1, sqrt(3)-1], [sqrt(3)-1, sqrt(3)+1]])/4
    a = simplify(s.H*s)
    assert g == g.T and hermitian_imag(g) == eye(2)/2
    jtrue = Matrix([1, -3 + 6*I])
    b = (eye(2)/2-g)*jtrue
    b = simplify(b)
    assert b == Matrix([-Q(3, 4)-2*I, Q(5, 4)+9*I/2])
    y = s*jtrue
    k = simplify(a + eye(2) - hermitian_real(g))
    l = simplify(s.H*y + b/2)
    j0 = simplify(k.inv()*l)
    e0 = simplify(b + g*j0)
    h = simplify(trace_norm_squared(y)-(l.H*j0)[0])
    assert k == diag(-Q(1, 2), Q(3, 2))
    assert j0 == Matrix([Q(5, 4)-I, -Q(5, 12)+7*I/2])
    assert e0[0] == Q(103, 48)-5*I/2 and h == Q(199, 48)
    x = symbols('x', real=True)
    strict_row = (eye(2)-x*g)[0, :]
    f = simplify(strict_row*Matrix([0, sqrt(Q(2, 3))]))[0]
    t0 = simplify(strict_row*Matrix([sqrt(2), 0]))[0]
    t = sqrt(24)*t0
    d = sqrt(24)*(x*e0[0]-j0[0])
    assert expand(24*f.conjugate()*f) == x**2
    kappa_squared = expand(t.conjugate()*t).subs(x, Q(8, 17))
    assert kappa_squared == Q(48, 17)
    m = Q(3, 2)
    upper_left = expand(t.conjugate()*t-1)
    determinant = expand(m*upper_left-d.conjugate()*d)
    p = Poly(96*determinant, x)
    assert upper_left.subs(x, Q(8, 17)) == Q(31, 17)
    assert p.as_expr() == 4367*x**2 - 3768*x + 864 and p.discriminant() == -894528
    assert p.LC() > 0
    d0 = expand(d.conjugate()*d)
    crude = d0.subs(x, 0) + d0.subs(x, 1)
    deduction = crude/(kappa_squared-1)
    new = h-m
    percentage = simplify(100*sqrt(new/trace_norm_squared(y)))
    assert new == Q(127, 48) and trace_norm_squared(y) == Q(43, 2)
    assert crude == Q(12937, 96) and deduction > h
    return dict(green=str(g), sensor_gram=str(a), drive=str(b), true_peak=2, cap=1,
                passive_radiation=str(hermitian_imag(g)), negative_inertia=1, source_count=1,
                stationary_current=str(j0), stationary_first_field=str(e0[0]), source_matrix=str(h),
                weighted_frame_squared=str(x**2), compensated_block_first_minor=str(upper_left),
                compensated_block_first_minor_min='31/17', compensated_block_determinant=str(determinant),
                determinant_num_discriminant=str(p.discriminant()), old_distance_squared='0',
                direct_source_deduction=str(m), new_distance_squared_lower=str(new),
                signal_frobenius_squared=str(trace_norm_squared(y)),
                lower_distance_percent_signal_exact=str(percentage),
                lower_distance_percent_signal_numeric=float(percentage),
                crude_endpoint_center_bound=str(crude), crude_kappa_squared=str(kappa_squared),
                crude_deduction=str(deduction), crude_deduction_numeric=float(deduction),
                crude_distance_squared_lower='0')


def sector_and_optical(j, g, b):
    je = j.H*(b+g*j)
    return j.H*j-hermitian_real(je), hermitian_imag(je)


def source_multipliers():
    s = Matrix([[1, Q(1, 2)], [0, sqrt(3)/2]])
    a = simplify(s.H*s)
    y = 2*s
    j = Matrix([[Q(19, 20), Q(1, 5)], [Q(1, 5), Q(19, 20)]])
    q = [sector_and_optical(j[i, :], zeros(1, 1), eye(2)[i, :]) for i in range(2)]
    assert [(x.trace(), w.trace()) for x, w in q] == [(-Q(3, 400), 0)]*2
    scalar_ceiling = trace_norm_squared(s*j-y)
    assert scalar_ceiling == Q(373, 200) > max(a.eigenvals())
    t = symbols('t', nonnegative=True)
    ks = [a+diag(2, t), a+diag(t, 2)]
    ls = [2*a[:, 0]+Matrix([1, 0]), 2*a[:, 1]+Matrix([0, 1])]
    for k in ks:
        assert k[0, 0].is_positive and k.det().is_positive
    bound = simplify(trace_norm_squared(y)-sum((l.H*k.inv()*l)[0] for l, k in zip(ls, ks)))
    assert simplify(bound-(2-6/(12*t+11))) == 0
    point = bound.subs(t, 3)
    assert point == Q(88, 47) and point-scalar_ceiling == Q(69, 9400)
    return dict(green='0', drive='I2', cap=1, true_peak=2,
                true_minimum_squared_residual='2', trial_current=str(j),
                per_pixel_trial_sector_trace='-3/400', per_pixel_trial_optical_trace='0',
                scalar_ceiling_if_at_most_one_positive_source_eigenvalue='3/2',
                scalar_ceiling_if_two_positive_source_eigenvalues=str(scalar_ceiling),
                matrix_bound_formula='2 - 6/(12*t + 11)', t=3,
                matrix_bound=str(point), strict_difference=str(point-scalar_ceiling))


def rank_obstruction():
    # One dark spatial coordinate and two independently illuminated visible rows.
    g = diag(0, 0, 2)
    b = eye(2).col_join(zeros(1, 2))
    s = b.T
    y = 2*eye(2)
    mean = 2*b
    j_true = mean
    witness_high = diag(2, 2, 0)
    assert j_true == witness_high*(b+g*j_true) and s*j_true == y
    q, o = sector_and_optical(mean, g, b)
    assert q == 2*eye(2) and o == zeros(2, 2)
    w = Matrix([0, 0, 1])
    assert s*w == zeros(2, 1)
    dr = (w.H*w-hermitian_real(w.H*g*w))[0]
    do = hermitian_imag(w.H*g*w)[0]
    currents = [mean+w*(sign*2*identity_row)
                for identity_row in (eye(2)[0, :], eye(2)[1, :]) for sign in (-1, 1)]
    covariance = sum(((j-mean).H*(j-mean)/4 for j in currents), zeros(2, 2))
    assert covariance == 2*eye(2)
    assert q+dr*covariance == zeros(2, 2) and o+do*covariance == zeros(2, 2)
    individual_ranks = [((j-mean).H*(j-mean)).rank() for j in currents]
    for j in currents:
        assert s*j == y and max(sector_and_optical(j, g, b)[0].eigenvals()) == 2
    mu = Q(2)
    k = s.H*s+mu*(eye(3)-hermitian_real(g))
    l = s.H*y+(mu/2)*b
    h = simplify(y.H*y-l.H*k.inv()*l)
    assert k == diag(3, 3, -2) and h == eye(2)
    return dict(green=str(g), drive=str(b), sensor=str(s), cap=1,
                true_material=str(witness_high), true_minimum_regular_squared_residual='2',
                receiver_exact_mean=str(mean), receiver_null_ray=str(w), covariance=str(covariance),
                discrete_support_size=4, support_probabilities='1/4',
                pointwise_sector_max_eigenvalue='2', average_sector='zero', average_optical='zero',
                individual_negative_update_ranks=individual_ranks,
                expected_negative_update_rank=covariance.rank(),
                moment_optimum_squared_distance='0', inertia_sector_weight=int(mu),
                inertia_state_matrix=str(k), inertia_source_matrix=str(h),
                original_inertia_squared_distance='1',
                scope='Passive reciprocal algebraic model. Regular solutions suffice, including the generating material. The box also contains a singular dark-only mode at its third coordinate 1/2.')


def strictly_passive_rank_obstruction():
    # Add two positive-stiffness dark cells whose mean optical transfer cancels
    # the covariance's radiation. Every real nonnegative box point is regular.
    epsilon = Q(1, 20)
    g = diag(I*epsilon, I*epsilon, 2+I*epsilon, I*epsilon, I*epsilon)
    b = ((1-2*I*epsilon)*eye(2)).col_join(zeros(1, 2)).col_join(-I*eye(2)/4)
    s = eye(2).row_join(zeros(2, 3))
    x_true = diag(2, 2, 0, 0, 0)
    true_current = simplify((eye(5)-x_true*g).inv()*x_true*b)
    y = 2*eye(2)
    assert s*true_current == y and hermitian_imag(g) == epsilon*eye(5)
    mean = y.col_join(zeros(1,2)).col_join(eye(2)/2)
    w = Matrix([0,0,1,0,0])
    q, o = sector_and_optical(mean, g, b)
    assert q == Q(9,4)*eye(2) and o == -Q(9,80)*eye(2)
    samples = [(mean,Q(1,2))]
    for a in [eye(2)[0,:], eye(2)[1,:]]:
        samples.extend([(mean+3*w*a,Q(1,8)), (mean-3*w*a,Q(1,8))])
    assert sum((j*p for j,p in samples),zeros(5,2)) == mean
    covariance = sum(((j-mean).H*(j-mean)*p for j,p in samples),zeros(2,2))
    assert covariance == Q(9,4)*eye(2)
    for j,_ in samples:
        assert s*j == y and max(sector_and_optical(j, g, b)[0].eigenvals()) == Q(9,4)
    qs, os = zip(*(sector_and_optical(j, g, b) for j,_ in samples))
    for matrices in [qs,os]:
        assert sum((v*p for v,(_,p) in zip(matrices,samples)),zeros(2,2)) == zeros(2,2)
    k = s.H*s+2*(eye(5)-hermitian_real(g))
    l = s.H*y+b
    h = simplify(y.H*y-l.H*k.inv()*l)
    assert k == diag(3,3,-2,2,2) and h == Q(2317,2400)*eye(2)
    return dict(green=str(g), drive=str(b), sensor=str(s), cap=1, true_material=str(x_true),
                radiation_spectral_floor=str(epsilon), entire_real_nonnegative_box_regular=True,
                original_inertia_squared_distance=str(h[0,0]), moment_optimum_squared_distance='0',
                true_minimum_squared_residual=str(2/(1+epsilon**2)),
                covariance=str(covariance), expected_negative_update_rank=2,
                max_pointwise_negative_update_rank=1, pointwise_sector_max_eigenvalue='9/4')


def lossy_negative_projection():
    # G=iI, alpha=1 and B=I give M=2I, Cb=I/8 and Jc=iI/4.
    w = Matrix([[1, 1]])
    gamma = (w*w.H/2)[0, 0]
    j_c = I*eye(2)/4
    completion_center = diag(Q(1, 2), Q(1, 4))
    z0 = w*(j_c-completion_center)
    v = Matrix([[Q(3, 5), 4*I/5]])
    contraction_lift = w.H*v/sqrt(2)
    assert (v*v.H)[0] == 1 and contraction_lift.H*contraction_lift == v.H*v
    j = simplify(j_c+contraction_lift/4)
    image = simplify(w*(j-completion_center))
    projected = simplify(z0+sqrt(gamma)*v*sqrt(Q(1, 8)))
    ellipsoid_slack = simplify(eye(2)/8-2*(j-j_c).H*(j-j_c))
    assert image == projected and set(ellipsoid_slack.eigenvals()) == {0, Q(1, 8)}
    assert trace_norm_squared(image) == trace_norm_squared(projected)
    return dict(known_alpha=1, negative_output_row=str(w), negative_output_gram=str(gamma),
                optical_center=str(j_c), completed_square_center=str(completion_center),
                center_in_negative_output=str(z0), boundary_contraction=str(v),
                full_ellipsoid_current_lift=str(j), negative_output_from_full_lift=str(image),
                negative_output_from_projected_ellipsoid=str(projected),
                full_ellipsoid_slack_eigenvalues=['0', '1/8'])


def receive_ellipsoid():
    # α=1 is legal only for the stated lossy model.
    g, b, s, alpha = I*eye(2), eye(2), Matrix([[1, 0]]), 1
    m = hermitian_imag(g)+alpha*eye(2)
    fc = I*s*m.inv()*b/2
    rs = (s*m.inv()*s.H)[0]
    cb = b.H*m.inv()*b/4
    coefficient = simplify(sqrt(rs)*sqrt(cb[0, 0]))
    assert coefficient == Q(1, 4) and fc == Matrix([[I/4, 0]])
    y = fc+Matrix([[0, Q(3, 4)]])
    z = Matrix([[0, 1]])
    nearest = fc+coefficient*z
    dual_value = (z.H*(y-fc)).trace()-coefficient
    assert trace_norm_squared(z) == 1 and dual_value == Q(1, 2)
    assert trace_norm_squared(y-nearest) == dual_value**2
    x = diag(I, 0)
    j = (eye(2)-x*g).inv()*x*b
    lossy_forward = s*j
    assert j == x*(b+g*j) and lossy_forward == Matrix([[I/2, 0]])
    assert trace_norm_squared(lossy_forward-fc) == coefficient**2
    return dict(alpha=alpha, green=str(g), drive=str(b), sensor=str(s), ellipsoid_center=str(fc),
                ellipsoid_ball_radius=str(coefficient), proposed_data=str(y),
                dual_witness=str(z), nearest_ellipsoid_point=str(nearest),
                primal_distance=str(dual_value), dual_distance=str(dual_value),
                feasible_lossy_material=str(x), feasible_lossy_forward=str(lossy_forward),
                physical_boundary_distance_from_center=str(coefficient))


def main():
    w, c = perf_counter(), process_time()
    rows = dict(arithmetic='SymPy exact rational and radical algebra',
                source_spectral_controls=source_spectral_controls(),
                nonconvex_with_no_negative_inertia=nonconvex_with_no_negative_inertia(),
                compensation=compensation(),
                source_multipliers=source_multipliers(), moment_rank_obstruction=rank_obstruction(),
                strictly_passive_moment_rank_obstruction=strictly_passive_rank_obstruction(),
                lossy_negative_projection=lossy_negative_projection(),
                lossy_receive_ellipsoid=receive_ellipsoid(), wall_seconds=perf_counter()-w,
                cpu_seconds=process_time()-c)
    (HERE/'exact_controls.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
