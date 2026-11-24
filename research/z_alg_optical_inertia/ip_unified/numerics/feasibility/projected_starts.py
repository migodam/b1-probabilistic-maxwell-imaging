"""Use inexpensive coarse material proposals, then audit each on the original native VIE."""
import argparse
import resource

import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from search_native import (CAPS, HERE, Model, N, SIDE, box_fit, clocks, elapsed, inputs, stem, write_json)
from lowdim import ProjectedFit, projected_model


def search(geometry, projection, obs, cap, args, shared):
    started = clocks()
    index, green, sensor, incident, _ = projection
    native = Model(*geometry, obs.y)
    fit = ProjectedFit(green, sensor, incident, obs.y)
    yy, xx = (np.indices((args.side, args.side))+.5)*2/args.side-1.
    gaussian = np.exp(-.5*((xx/.45)**2+(yy/.45)**2)).ravel()
    target = .998*obs.padded_radius
    path, shortlisted = [], []

    def evaluate(m, row):
        try:
            r, p, j = native.full(m[index], gradient=False)
        except ArithmeticError as exc:
            row['native_solve_failure'] = str(exc)
            return
        row['native_residual'] = float(r)
        shortlisted.append((r, m[index], p, j))

    for name, start, schedule in [('vacuum', np.zeros_like(gaussian), (None,)),
                                  ('generic_gaussian', .7*cap*gaussian, (None,)),
                                  ('source_field_continuation', .7*cap*gaussian, (.0001, .003, .05, None))]:
        seed = start
        for weight in schedule:
            t = clocks()

            def loss(m):
                value, gradient, _ = fit.calculate(m, weight)
                return value, gradient

            calc = minimize(loss, seed, method='L-BFGS-B', jac=True,
                            bounds=[(0., cap)]*len(seed), options={'ftol':1.e-11, 'gtol':1.e-9,
                            'maxiter':args.coarse_calls, 'maxfun':args.coarse_calls, 'maxcor':10})
            seed = np.clip(calc.x, 0, cap)
            _, _, r = fit.calculate(seed, weight)
            row = dict(seed=name, stage='projected_seed', penalty=weight, low_dim_residual=r,
                       evaluations=int(calc.nfev), stop=str(calc.message))
            if weight is None:
                evaluate(seed, row)
            row.update(elapsed(t))
            path.append(row)
        if shortlisted and min(x[0] for x in shortlisted) <= target:
            break
    if not shortlisted:
        res, ports, currents = native.full(np.zeros(N), gradient=False)
        shortlisted.append((res, np.zeros(N), ports, currents))
    r, m, ports, current = min(shortlisted, key=lambda x: x[0])
    if r > target and args.full_calls:
        try:
            candidate, stage = box_fit(native, m, np.arange(N), SIDE, cap, target,
                                       args.full_calls, args.full_cpu, True)
        except ArithmeticError as exc:
            path.append(dict(stage='native_9216_corrector', stop=str(exc)))
        else:
            path.append(dict(stage='native_9216_corrector', **stage))
            if candidate['residual'] < r:
                r, m, ports, current = (candidate[key] for key in ('residual', 'm', 'ports', 'j'))
    if m.min() < 0 or m.max() > cap or not np.all(np.isfinite(m)):
        raise ArithmeticError('Selected material violates the frozen box.')
    import search as old_search
    old_rtol, old_search.RTOL = old_search.RTOL, 1.e-10
    tightened_solve_failure = None
    try:
        r, ports, current = native.full(m, (current,), gradient=False)
    except ArithmeticError as exc:
        tightened_solve_failure = str(exc)
    finally:
        old_search.RTOL = old_rtol
    defect = m[:,None]*(native.b+native.green(current))-current
    base = args.out/stem(obs, cap)
    suffix = base.suffix
    np.savez_compressed(base.with_suffix(suffix+'.npz'), material=m, Y=obs.y,
                        Y_prediction=ports, current_approx=current)
    receipt = dict(observation=obs.id, cap=cap, radius=obs.padded_radius,
                   status='numerical_interior_pending_enclosure' if r < obs.padded_radius
                   else 'numerical_upper_candidate_above_radius', numerical_residual=float(r),
                   residual_over_radius=float(r/obs.padded_radius), native_peak=float(m.max()),
                   state_defect_numeric_frobenius=float(np.linalg.norm(defect)),
                   positive_cells=int(np.count_nonzero(m)), projected_side=args.side,
                   stages=path, tighter_native_recheck_failure=tightened_solve_failure, counts=dict(native.count), projected_counts=dict(fit.count),
                   shared_setup=shared, peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   **elapsed(started))
    write_json(base.with_suffix(suffix+'.json'), receipt)
    print(obs.id, cap, receipt['residual_over_radius'], receipt['wall_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--objects', nargs='+')
    parser.add_argument('--caps', nargs='+', type=float, choices=CAPS, default=CAPS)
    parser.add_argument('--side', type=int, default=16)
    parser.add_argument('--coarse-calls', type=int, default=36)
    parser.add_argument('--full-calls', type=int, default=4)
    parser.add_argument('--full-cpu', type=float, default=15.)
    parser.add_argument('--out', type=type(HERE), default=HERE/'projected_runs')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    geo, cases = inputs()
    t = clocks()
    projected = projected_model(*geo, next(iter(cases.values())).y, args.side)
    shared = dict(counts=dict(projected[-1]), **elapsed(t))
    write_json(args.out/'setup.json', shared)
    for cap in args.caps:
        for name in args.objects or cases:
            obs = cases[name]
            if not (args.out/(stem(obs,cap)+'.json')).exists():
                search(geo, projected, obs, cap, args, shared)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
