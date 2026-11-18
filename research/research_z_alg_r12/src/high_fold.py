"""Feasible R11 directional fold proposals with full VIE and equal-ceiling controls."""
from argparse import ArgumentParser
import json
import time

import numpy as np
from scipy.linalg import lu_solve, svd
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

from high_endpoints import ENDPOINTS, R12, around_endpoint, bounds, load_endpoint, noise_radius, physics


class BudgetReached(Exception):
    pass


def directional_hessian(model, state, direction):
    """D²F(m)[v,v], retaining the factorization from the matching VIE forward."""
    with model.operation('DIRECTIONAL_HESSIAN', rhs=2*model.ns, tangent_RHS=2*model.ns):
        d = direction[:, None]
        j1 = lu_solve(model.factor, d * state['E'], check_finite=False)
        j2 = lu_solve(model.factor, d * (model.G @ j1), check_finite=False)
        return model.pack(2 * (model.S @ j2))


def ray_range(material, direction, lo, hi):
    active = np.abs(direction) > 1e-12
    x, d = material[active], direction[active]
    a, b = (lo-x)/d, (hi-x)/d
    # Delimit large extrapolation irrespective of individual-cell availability.
    cap = 2 * np.sqrt(material.size)
    return max(np.max(np.minimum(a, b)), -cap), min(np.min(np.maximum(a, b)), cap)


def real_quadratic_roots(a, b, c):
    if a == 0 or b*b < 4*a*c:
        return []
    disc = np.sqrt(b*b - 4*a*c)
    t = -.5 * (b + np.copysign(disc, b))
    return [t/a, c/t] if t != 0 else []


def rms(x):
    return float(np.linalg.norm(x) / np.sqrt(x.size))


def initial_proposals(model, material, target, lo, hi, delta):
    t0 = time.perf_counter()
    state = model.forward(material)
    jacobian = model.jacobian(state)
    free = (material > lo+1e-6) & (material < hi-1e-6)
    u, sigma, vh = svd(jacobian[:, free], full_matrices=False, check_finite=False)
    if not len(sigma):
        return [], [], [], None
    eps = np.finfo(float).tiny
    residual = model.pack(state['y'] - target)
    indexes = list(dict.fromkeys([len(sigma)-1, max(0, len(sigma)-2),
                                 int(np.argmin(np.abs(np.log(np.maximum(sigma, eps)/(delta/np.sqrt(material.size))))))]))
    proposed, weak_options, diagnostics = [], [], []
    numerical_range = sigma > sigma[0]*max(jacobian.shape)*np.finfo(float).eps
    for rank in indexes:
        s, ell = sigma[rank], u[:, rank]
        v = np.zeros(material.size)
        v[free] = vh[rank]
        q = directional_hessian(model, state, v)
        kappa = float(ell @ q)
        q_norm = float(np.linalg.norm(q))
        strong = sigma > max(10*s, sigma[0]*1e-7)
        normal = q - u[:, strong] @ (u[:, strong].T @ q)
        orthogonal = normal - kappa * ell
        left_null = q - u[:, numerical_range] @ (u[:, numerical_range].T @ q)
        normal_norm = float(np.linalg.norm(normal))
        correction = np.zeros(material.size)
        correction[free] = -vh[strong].T @ ((u[:, strong].T @ q) / sigma[strong])
        lower_ray, upper_ray = ray_range(material, v, lo, hi)
        r_l = float(ell @ residual)
        row = dict(singular_index=int(rank+1), sigma=float(s), strong_directions=int(strong.sum()),
                   hessian_response_norm=q_norm, kappa=kappa, normal_response_norm=normal_norm,
                   left_null_response_norm=float(np.linalg.norm(left_null)),
                   curvature_normal_fraction=normal_norm/max(q_norm, eps),
                   weak_alignment_of_normal=abs(kappa)/max(normal_norm, eps),
                   weak_alignment_of_all_curvature=abs(kappa)/max(q_norm, eps),
                   orthogonal_normal_fraction=float(np.linalg.norm(orthogonal))/max(normal_norm, eps),
                   signed_weak_data_residual=r_l, feasible_ray=[float(lower_ray), float(upper_ray)],
                   scalar_candidates=[])
        steps = [('same_response', float(-2*s/kappa))] if kappa != 0 else []
        steps.extend(('target_projection', float(a)) for a in real_quadratic_roots(.5*kappa, s, r_l))
        for kind, step in steps:
            reason = None
            if not np.isfinite(step) or not lower_ray <= step <= upper_ray:
                reason = 'outside_feasible_ray'
            elif abs(step) < .01*np.sqrt(material.size):
                reason = 'sub_0.01_material_RMS'
            normal_defect = .5*step**2*np.linalg.norm(orthogonal) if np.isfinite(step) and abs(step) < 1e100 else None
            item = dict(kind=kind, step=step if np.isfinite(step) else None,
                        ray_material_RMS=abs(step)/np.sqrt(material.size) if np.isfinite(step) else None,
                        orthogonal_second_order_defect_over_delta=float(normal_defect/delta) if normal_defect is not None else None,
                        rejection=reason)
            row['scalar_candidates'].append(item)
            if reason:
                continue
            simple = material + step*v
            corrected = simple + .5*step**2*correction
            corrected_ok = np.min(corrected) >= lo and np.max(corrected) <= hi
            candidate = corrected if corrected_ok else simple
            predicted = residual + step*(jacobian @ v) + .5*step**2*q
            if corrected_ok:
                predicted = predicted + .5*step**2*(jacobian @ correction)
            item.update(strong_correction_used=bool(corrected_ok), material_RMS=rms(candidate-material),
                        predicted_residual_over_delta=float(np.linalg.norm(predicted)/delta))
            if all(rms(candidate-prev[1]) > .01 for prev in proposed):
                proposed.append((f'{rank+1}:{kind}', np.clip(candidate, lo, hi), item))
        weak_options.append((v, lower_ray, upper_ray, row))
        diagnostics.append(row)
    # The slope and Hessian agree with a central physical derivative along one
    # feasible ray. This compares the new second derivative with the legacy J.
    validation = None
    for v, a, b, row in weak_options:
        if a < -.003 and b > .003:
            plus = model.jacobian(model.forward(material + .003*v)) @ v
            minus = model.jacobian(model.forward(material - .003*v)) @ v
            estimate = (plus-minus)/.006
            # Rebuild the base once to evaluate the same direction at the anchor.
            q = directional_hessian(model, model.forward(material), v)
            validation = dict(rank=row['singular_index'], step=.003,
                              derivative_relative_error=float(np.linalg.norm(estimate-q)/max(np.linalg.norm(q), eps)))
            break
    return proposed, weak_options, diagnostics, dict(validation=validation, wall_seconds=time.perf_counter()-t0)


def control_seeds(material, proposed, weak_options, lo, hi, shape, maximum):
    weak = []
    for v, a, b, row in weak_options:
        sized = [rms(m-material)*np.sqrt(material.size) for _, m, item in proposed
                 if any(item is other for other in row['scalar_candidates'])]
        for radius in sized or [.5 * max(-a, b)]:
            for s in (-1., 1.):
                step = float(np.clip(s * radius, a, b))
                candidate = np.clip(material+step*v, lo, hi)
                if rms(candidate-material) >= .01 and all(rms(candidate-item[1])>.01 for item in weak):
                    weak.append((f'weak-{row["singular_index"]}:{step:+.4g}', candidate, {}))
    rng = np.random.default_rng(2026111903)
    # Generic perturbed starts share the RMS of accepted proposals when present.
    # With no eligible fold start, this arm still tests finite nonlocal starts.
    radii = [rms(m-material) for _, m, _ in proposed[:maximum]] or [.5, 1., 2.]
    random = []
    for index in range(maximum):
        direction = gaussian_filter(rng.standard_normal(shape), 1.25).ravel()
        radius = radii[index % len(radii)]
        if not proposed and index % 2:
            seed = np.full(material.size, .1)
        else:
            seed = material
        candidate = np.clip(seed+radius*direction/rms(direction), lo, hi)
        random.append((f'smooth-perturbation-{index+1}', candidate, {}))
    return weak[:maximum], random


def trf(model, reg, m, target, lo, hi, max_forwards, max_rhs):
    """Original Phi as a sum of smooth least squares, with an exact material J."""
    scale = np.linalg.norm(target)
    smooth = np.sqrt(2*reg.smooth/reg.edges)
    tv = np.sqrt(2*reg.tv/reg.edges)
    difference = reg.D.toarray()
    calls0 = model.counts.copy()
    t0, cpu0 = time.perf_counter(), time.process_time()
    best = float('inf'), m.copy()
    history = []

    def guard(fw=0, rhs=0):
        if model.counts['FULL_FORWARD']-calls0['FULL_FORWARD']+fw > max_forwards or \
           model.counts['actual_complex_RHS']-calls0['actual_complex_RHS']+rhs > max_rhs:
            raise BudgetReached

    def state(x):
        if model.state is None or not np.array_equal(model.state['m'], x):
            guard(1, model.ns)
            model.forward(x)
        return model.state

    def residual(x):
        nonlocal best
        output = model.pack(state(x)['y'] - target) / scale
        d = reg.D @ x
        qs = np.sqrt(d*d + reg.epsilon**2)
        v = np.r_[output, smooth*d, tv*d/np.sqrt(qs+reg.epsilon)]
        phi = .5*float(v @ v)
        if not history or not np.array_equal(x, history[-1][1]):
            history.append((phi, x.copy()))
        if phi < best[0]:
            best = phi, x.copy()
        return v*100

    def jacobian(x):
        base = state(x)
        if 'J' not in base:
            guard(0, model.nr)
        J = model.jacobian(base) / scale
        d = reg.D @ x
        qs = np.sqrt(d*d+reg.epsilon**2)
        weight = .5*tv*np.sqrt(qs+reg.epsilon)/qs
        return 100*np.vstack((J, smooth*difference, weight[:, None]*difference))

    try:
        solution = least_squares(residual, m, jac=jacobian, method='trf', tr_solver='lsmr',
                                 bounds=(lo, hi), ftol=1e-13, xtol=1e-13, gtol=1e-10,
                                 tr_options=dict(atol=1e-6, btol=1e-6, maxiter=32),
                                 max_nfev=max_forwards)
        status = str(solution.message)
    except BudgetReached:
        status = 'PHYSICAL_FORWARD_OR_RHS_CEILING'
    return best[1], dict(best_Phi=best[0], status=status, evaluations=len(history),
                         phi_trace=[x[0] for x in history],
                         cost=dict(wall_seconds=time.perf_counter()-t0, cpu_seconds=time.process_time()-cpu0,
                                   counts=dict(model.counts-calls0)))


def run(label, output, max_forwards, max_rhs, maximum):
    spec, public, material = load_endpoint(label)
    lo, hi = bounds(public)
    material = np.clip(material, lo, hi)
    target, delta = public['Y'], noise_radius(public)
    proposed_model, reg = physics(public)
    candidates, weak_args, modes, setup = initial_proposals(proposed_model, material, target, lo, hi, delta)
    candidates = candidates[:maximum]
    weak, generic = control_seeds(material, candidates, weak_args, lo, hi, public['grid_shape'], maximum)
    setup_cost = dict(proposed_model.counts)
    # Diagnostic derivative validation is a shared measured overhead. The fold
    # discovery arm pays the actual initial forward, receiver-J and q solves.
    derivative_check_rhs = 3*proposed_model.ns + 2*proposed_model.nr + 2*proposed_model.ns if setup and setup['validation'] else 0
    derivative_check_fw = 3 if derivative_check_rhs else 0
    arms = {}
    arrays = {'initial': material}
    for arm, seeds, setup_rhs, setup_fw in (
        ('fold', candidates, setup_cost.get('actual_complex_RHS', 0)-derivative_check_rhs,
         setup_cost.get('FULL_FORWARD', 0)-derivative_check_fw),
        ('weak', weak, proposed_model.ns+proposed_model.nr, 1),
        ('multistart', generic, 0, 0),
        ('continue_fullj', [('archived-incumbent', material, {})], 0, 0),
    ):
        model, reg = physics(public)
        t0 = time.perf_counter()
        cpu0 = time.process_time()
        rows = []
        champion, winner = material.copy(), None
        initial = around_endpoint(model, reg, material, target, lo, hi, delta, with_jacobian=False)
        paid0 = model.counts.copy()
        best_phi = initial['phi']
        if not seeds and arm == 'fold':
            seeds = [('fallback-to-incumbent', material, {})]
        elif not seeds:
            seeds = [('no-finite-weak-start', material, {})]
        for group_index, (seed_name, seed, proposal_info) in enumerate(seeds):
            model.forward(seed)
            seed_relative = float(np.linalg.norm(model.state['y']-target) / delta)
            remaining = len(seeds)-group_index
            budget_fw = max(2, (max_forwards - setup_fw - (model.counts['FULL_FORWARD']-paid0['FULL_FORWARD'])) // remaining)
            budget_rhs = max(model.nr+model.ns, (max_rhs - setup_rhs - (model.counts['actual_complex_RHS']-paid0['actual_complex_RHS'])) // remaining)
            answer, row = trf(model, reg, seed, target, lo, hi, budget_fw, budget_rhs)
            rows.append(dict(seed=seed_name, seed_RMS=rms(seed - material), seed_residual_over_delta=seed_relative,
                             proposal=proposal_info, **row))
            arrays[f'{arm}_seed{group_index}'] = seed
            arrays[f'{arm}_candidate{group_index}'] = answer
            if row['best_Phi'] < best_phi:
                best_phi, champion, winner = row['best_Phi'], answer, seed_name
        paid = dict(model.counts-paid0)
        output_diag = around_endpoint(model, reg, champion, target, lo, hi, delta, with_jacobian=False)
        arrays[arm] = champion
        arms[arm] = dict(winner=winner, initial=initial, result=output_diag, trials=rows,
                         resulting_RMS=rms(champion - material), charged_setup_rhs=setup_rhs,
                         charged_setup_forwards=setup_fw, physical_search_counts=paid,
                         matched_search_rhs=paid.get('actual_complex_RHS', 0)+setup_rhs,
                         matched_search_forwards=paid.get('FULL_FORWARD', 0)+setup_fw,
                         final_diagnostic_counts=dict(model.counts - paid0 - model.counts.__class__(paid)),
                         online_and_diagnostic_counts=dict(model.counts),
                         wall_seconds=time.perf_counter()-t0, cpu_seconds=time.process_time()-cpu0)
    row = dict(label=label, case_id=spec.case_id, archived_method=spec.method,
               scope='exposed archived scalar-TM R11 transfer; candidates use old data and declared bounds only',
               noise_radius=delta, proposed_modes=modes, eligible_finite_fold_seeds=len(candidates),
               proposal_setup=setup, proposal_setup_total_counts=setup_cost,
               common_max_forwards=max_forwards, common_max_rhs=max_rhs, arms=arms)
    np.savez_compressed(output/f'{label}.npz', **arrays)
    (output/f'{label}.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
    summary = {name: (arm['winner'], round(arm['result']['residual_over_noise_radius'], 3),
                      arm['matched_search_forwards'], arm['matched_search_rhs']) for name, arm in arms.items()}
    print(label, len(candidates), summary, flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('case', nargs='+', choices=list(ENDPOINTS))
    parser.add_argument('--max-forwards', type=int, default=150)
    parser.add_argument('--max-rhs', type=int, default=4800)
    parser.add_argument('--max-seeds', type=int, default=3)
    parser.add_argument('--output', default=R12/'results/high_fold', type=type(R12))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for case in args.case:
        run(case, args.output, args.max_forwards, args.max_rhs, args.max_seeds)
