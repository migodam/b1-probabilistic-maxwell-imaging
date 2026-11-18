"""Exact reduced Hessian on a box face for the archived physical objective."""
from argparse import ArgumentParser
import json
import time

import numpy as np
from scipy.linalg import eigh, lu_solve

from high_endpoints import ENDPOINTS, R12, bounds, evaluate, gradient_mapping, load_endpoint, noise_radius, physics


def hessians(model, reg, material, target, free):
    state = model.forward(material)
    jac = model.jacobian(state)[:, free]
    n = int(free.sum())
    with model.operation('FACE_HESSIAN', rhs=n+model.nr, face_inverse_RHS=n,
                         face_receiver_transpose_RHS=model.nr):
        left = lu_solve(model.factor, model.S.T, trans=1, check_finite=False)[free].T
        transfer = model.G[free] @ lu_solve(model.factor, model.identity[:, free], check_finite=False)
        fields = state['E'][free]
        residual = state['y'] - target
        mixed = transfer * ((left.T @ residual.conj()) @ fields.T)
        yn2 = np.linalg.norm(target)**2
        h_data = (jac.T @ jac + (mixed + mixed.T).real)/yn2
        d = reg.D @ material
        q = np.sqrt(d*d + reg.epsilon**2)
        weights = (2*reg.smooth + reg.tv*reg.epsilon**2/q**3)/reg.edges
        D = reg.D[:, free]
        prior_hessian = (D.T @ D.multiply(weights[:, None])).toarray()
    return h_data, h_data + prior_hessian


def spectrum(h, scale):
    values, vectors = eigh(h, check_finite=False)
    return dict(min=float(values[0]), max=float(values[-1]), negative=int(np.count_nonzero(values < -scale)),
                nearly_zero=int(np.count_nonzero(abs(values) <= scale)), tolerance=scale), vectors[:, 0]


def run(label, validate):
    spec, public, m = load_endpoint(label)
    model, reg = physics(public)
    lo, hi = bounds(public)
    m = np.clip(m, lo, hi)
    target, delta = public['Y'], noise_radius(public)
    t0, c0 = time.perf_counter(), time.process_time()
    values, gd, gr = evaluate(model, reg, m, target)
    g = gd + gr
    lower_active = m <= lo+1e-8
    upper_active = m >= hi-1e-8
    free = ~(lower_active | upper_active)
    if not free.any():
        raise ValueError(f'{label} has no numerically free directions')
    h_data, h_phi = hessians(model, reg, m, target, free)
    active_gradients = np.r_[g[lower_active], -g[upper_active]]
    tolerance = max(np.linalg.norm(h_phi, np.inf)*1e-10, 1e-13)
    phi_spectrum, weakest = spectrum(h_phi, tolerance)
    data_spectrum, _ = spectrum(h_data, tolerance)
    counts0 = model.counts.copy()
    duration = time.perf_counter()-t0
    checks = None
    if validate:
        v = np.zeros(m.size)
        rng = np.random.default_rng(20261119)
        v[free] = x = rng.normal(size=free.sum())
        v /= np.linalg.norm(x)
        limits = np.minimum(m[free]-lo, hi-m[free])/np.maximum(abs(v[free]), 1e-20)
        # The derivative of a smooth objective may be computed just outside
        # the face if needed; this constrained case admits a finite legal ray.
        step = min(.001, .1*np.min(limits))
        if step > 1e-7:
            g_plus = sum(evaluate(model, reg, m+step*v, target)[1:])
            g_minus = sum(evaluate(model, reg, m-step*v, target)[1:])
            fd = (g_plus[free]-g_minus[free])/(2*step)
            exact = h_phi @ v[free]
            checks = dict(step=float(step), relative_error=float(np.linalg.norm(fd-exact)/np.linalg.norm(exact)),
                          counts=dict(model.counts-counts0))
    newton_cost = None
    if phi_spectrum['min'] > 0:
        step_free = -np.linalg.solve(h_phi, g[free])
        p = m.copy()
        p[free] = np.clip(p[free]+step_free, lo, hi)
        dx = p-m
        calls = model.counts.copy()
        newton, nd, nr = evaluate(model, reg, p, target)
        newton_cost = dict(clipped_face_newton_RMS=float(np.linalg.norm(dx)/np.sqrt(m.size)),
                           new_Phi=float(newton['phi']), new_residual_over_noise=float(newton['data_residual']/delta),
                           projected_gradient_inf=float(np.linalg.norm(gradient_mapping(p, nd+nr, lo, hi), np.inf)),
                           counts=dict(model.counts-calls))
    row = dict(case=label, case_id=spec.case_id, phi=values['phi'], model_data_residual_over_noise=values['data_residual']/delta,
               projected_gradient_inf=float(np.linalg.norm(gradient_mapping(m, g, lo, hi), np.inf)),
               free_gradient_inf=float(np.linalg.norm(g[free], np.inf)), active_lower=int(lower_active.sum()),
               active_upper=int(upper_active.sum()), free=int(free.sum()),
               signed_active_gradient_min=float(active_gradients.min()) if len(active_gradients) else None,
               active_gradient_wrong_sign=int(np.count_nonzero(active_gradients < 0)),
               phi_hessian_on_free_box_face=phi_spectrum, data_hessian_on_free_box_face=data_spectrum,
               finite_difference=checks, single_face_newton_step=newton_cost,
               face_hessian_counts=dict(counts0), face_hessian_wall_seconds=duration,
               all_probe_counts=dict(model.counts), wall_seconds=time.perf_counter()-t0, cpu_seconds=time.process_time()-c0)
    return row


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=list(ENDPOINTS), nargs='+')
    args = parser.parse_args()
    results = [run(label, validate=(i==0)) for i, label in enumerate(args.case)]
    (R12/'results/high_hessian.json').write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    for row in results:
        print(row['case'], row['phi_hessian_on_free_box_face'], 'active sign min', row['signed_active_gradient_min'],
              'wrong signs', row['active_gradient_wrong_sign'], 'check', row['finite_difference'])
