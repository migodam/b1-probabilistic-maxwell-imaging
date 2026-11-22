#!/usr/bin/env python3
"""After the public scan, compare saved and full-wave GN endpoints and radii."""
from scan import CANONICAL, HIGH, LA, LEGACY, OUTPUT, R8, REGIONS, SRAP, WORKSPACE, observations, read, write

import argparse
from collections import Counter
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve

sys.path.insert(0, str(HIGH))
from replay import ScalarTM
sys.path.insert(0, str(SRAP))
import legacy_la
from data2d import CONFIRM_ARRAYS, NOISE_SEED
from core.inertia import clock, timing

NAMES = dict(zip(CANONICAL + LEGACY, ('Austria · full 20×40', 'SRAP core–shell · half',
    'SRAP core–shell · full', 'SRAP double ring · half', 'Eccentric shell · half',
    'Legacy core–shell · half', 'Legacy double ring · half'), strict=True))


def evaluate(m, model, Y):
    state = model.forward(m)
    return state, float(np.linalg.norm(state['y']-Y))


def gauss_newton(m, model, Y, upper, rounds=6):
    """Fixed projected damped GN, accepted using the actual Maxwell residual."""
    before, started = model.counts.copy(), clock()
    m = np.clip(m, 0., upper)
    state, residual = evaluate(m, model, Y)
    history = [dict(residual=residual, accepted_scale=None)]
    for iteration in range(rounds):
        if model.state is not state:
            state, residual = evaluate(m, model, Y)
        J = model.jacobian(state)
        F = model.pack(state['y']-Y)
        rho = .03 * np.vdot(J, J).real / min(J.shape)
        if J.shape[0] <= J.shape[1]:
            H = J @ J.T
            H.flat[::len(H)+1] += rho
            delta = -J.T @ cho_solve(cho_factor(H, check_finite=False), F, check_finite=False)
        else:
            H = J.T @ J
            H.flat[::len(H)+1] += rho
            delta = -cho_solve(cho_factor(H, check_finite=False), J.T @ F, check_finite=False)
        accepted = False
        for scale in (1., .5, .25, .125, .0625):
            trial = np.clip(m+scale*delta, 0., upper)
            p, gap = evaluate(trial, model, Y)
            if gap < residual:
                m, state, residual = trial, p, gap
                history.append(dict(residual=residual, accepted_scale=scale))
                accepted = True
                break
        if not accepted:
            break
    return m, dict(bounds_upper=upper, proposals=iteration+1, accepted=len(history)-1,
                   residual_history=history, operations=dict(model.counts-before), **timing(started))


def archived(case, model, Y):
    if case == CANONICAL[0]:
        p = HIGH / f'results/{case}.npz'
        with np.load(p, allow_pickle=False) as z:
            raw, light = z['history__cie_dual__m'], z['history__physical_lbfgs__m']
        path = HIGH / f'results/{case}__comparison_arrays.npz'
        opts = ('task_constrained_gauss_newton__minus__100', 'task_constrained_gauss_newton__plus__100')
        with np.load(path, allow_pickle=False) as z:
            selected = min(opts, key=lambda x: np.linalg.norm(z[x+'__prediction']-Y))
            gn = z[selected+'__material']
        provenance = dict(dual_raw={'path': str(p.relative_to(WORKSPACE)), 'array':'history__cie_dual__m'},
                          physical_lbfgs={'path': str(p.relative_to(WORKSPACE)), 'array':'history__physical_lbfgs__m'},
                          physical_gn={'path':str(path.relative_to(WORKSPACE)), 'array':selected+'__material'})
    else:
        p = SRAP / f'results/predictions/{case}__dual.npz'
        with np.load(p, allow_pickle=False) as z:
            raw, light = z['m_30_raw'], z['m_30_light']
            receipts = read_history(z['rows_json'].item())
        gn, trace = gauss_newton(light, model, Y, upper=20.)
        provenance = dict(dual_raw={'path':str(p.relative_to(WORKSPACE)), 'array':'m_30_raw', 'receipt':receipts['raw']},
                          physical_lbfgs={'path':str(p.relative_to(WORKSPACE)), 'array':'m_30_light', 'receipt':receipts['light']},
                          physical_gn={'path':'results/roi/endpoint_arrays.npz', 'array':f'{case}__physical_gn__m',
                                       'warm_start':'m_30_light', 'trace':trace})
    return {'dual_raw':raw, 'physical_lbfgs':light, 'physical_gn':gn}, provenance


def read_history(text):
    import json
    return {r['view']:r for r in json.loads(text) if r['budget']==30}


def truth_and_direction(case, problem):
    """Recover saved clean numerical field from frozen seeds after sealing the scan."""
    Y = problem.Y
    if case == CANONICAL[0]:
        root = R8 / 'research/research_max_struct_r8'
        path = root / f'data/paper_benchmarks_private/{case}.npz'
        with np.load(path, allow_pickle=False) as z:
            m = z['m'].ravel()
        contract = read(root / 'delegated/paper_baselines/PAPER_CONTRACTS.json')['NIE2016']['cases']
        seed = 8600 + next(i for i,v in enumerate(contract) if v['id']==case)
        radius = read(HIGH / 'observation_radii.json')[case]['original_R8_noise_ball']
    elif case in CANONICAL:
        with np.load(SRAP / f'data/truth/{case}.npz', allow_pickle=False) as z:
            m = z['m'].ravel()
        with np.load(SRAP / f'data/public/{case}.npz', allow_pickle=False) as z:
            radius = float(z['delta'])
        object_id, view = case.split('__')
        seed = np.random.SeedSequence([NOISE_SEED, int(object_id.split('_')[1]),
                                      tuple(a.lower() for a in CONFIRM_ARRAYS).index(view), 0])
    else:
        from json import loads
        with np.load(LA / f'.runs/public/{case}.npz', allow_pickle=False) as z:
            metadata = loads(z['case_json'].item())
            m = legacy_la.coarse_truth(metadata, z['centers'])
            radius = float(z['delta'])
        seed = metadata['noise_seed']
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal(Y.shape) + 1j*rng.standard_normal(Y.shape)
    clean = Y - Z * (radius/np.linalg.norm(Z))
    return m, clean, radius


def metrics(m, truth, residual, Y, radii):
    difference = m-truth
    reference = np.linalg.norm(truth)
    return dict(port_residual=residual, port_relative_residual=residual/np.linalg.norm(Y),
                numerical_inside_radii={k: bool(residual <= r) for k,r in radii.items()},
                material_relative_L2=np.linalg.norm(difference)/reference,
                target_relative_L2=np.linalg.norm(difference[truth>0])/reference,
                background_relative_L2=np.linalg.norm(difference[truth==0])/reference,
                peak=float(m.max()), over_U4=int(np.count_nonzero(m>4.)), over_U20=int(np.count_nonzero(m>20.)))


def endpoint_score(output):
    scan = read(output / 'scan_U20.json')
    sealed_targets = {case for g in scan['geometries'] for case in g['case_results']}
    if sealed_targets != set(CANONICAL+LEGACY):
        raise RuntimeError('Complete the entire historical public scan before scoring any arrays.')
    results, arrays = {}, {}
    for case, p, radii, path in observations():
        started = clock()
        model = ScalarTM({'G':p.G, 'E0':p.B, 'S':p.S})
        methods, provenance = archived(case, model, p.Y)
        truth, clean, actual_noise = truth_and_direction(case, p)
        proxy, proxy_residual = evaluate(truth, model, p.Y)
        for label, value in (('truth_coarse_m', truth), ('proxy_prediction', proxy['y']), ('clean_fine_Y', clean)):
            arrays[f'{case}__{label}'] = value
        discrepancy = np.linalg.norm(proxy['y']-clean)
        radii = {**radii, 'conditional_noise_plus_coarse_proxy':actual_noise+discrepancy}
        method_scores = {}
        for method, m in methods.items():
            state, residual = evaluate(m, model, p.Y)
            method_scores[method] = {**metrics(m, truth, residual, p.Y, radii), 'provenance':provenance[method]}
            arrays[f'{case}__{method}__m'] = m
            arrays[f'{case}__{method}__prediction'] = state['y']
        projections_before, projection_start = Counter(model.counts), clock()
        regions = []
        for index, inside in enumerate(REGIONS):
            region = {'row':index//4, 'col':index % 4,
                      'reference_peak':float(truth[inside].max()),
                      'reference_mean':float(truth[inside].mean()),
                      'endpoints':{name:{'peak':float(m[inside].max()), 'mean':float(m[inside].mean())}
                                   for name,m in methods.items()}}
            candidate = methods['physical_gn'].copy()
            candidate[inside] = 0.
            state, gap = evaluate(candidate, model, p.Y)
            if np.max(candidate)>20. or np.min(candidate)<0.:
                raise ArithmeticError('Projection outside the public material box.')
            region['literal_absence_projection'] = dict(residual=gap, box_U=20.,
                numerical_inside_radii={k:bool(gap<=r) for k,r in radii.items()},
                state_relative_residual=state['state_relative_residual'])
            arrays[f'{case}__roi_absent_{index:02d}__prediction'] = state['y']
            regions.append(region)
        results[case] = dict(label=NAMES[case], measurement=path, Y_norm=float(np.linalg.norm(p.Y)), radii=radii,
            coarse_proxy_discrepancy_fro=float(discrepancy), coarse_proxy_observation_residual=proxy_residual,
            coarse_proxy_relative_discrepancy=float(discrepancy/np.linalg.norm(clean)),
            coarse_proxy_contrast_max=float(truth.max()), regions=regions, methods=method_scores,
            projection_receipt={'operations':dict(model.counts-projections_before), **timing(projection_start)},
            all_operations=dict(model.counts), **timing(started))
        print(case, 'proxy mismatch', results[case]['coarse_proxy_relative_discrepancy'],
              'GN rel residual', method_scores['physical_gn']['port_relative_residual'], flush=True)
    write(output / 'endpoint_metrics.json', results)
    np.savez_compressed(output / 'endpoint_arrays.npz', **arrays)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    endpoint_score(args.output)


if __name__ == '__main__':
    main()
