"""Search the frozen 32-by-32 noise balls with the historical native TM solver."""
import argparse
import json
from pathlib import Path
import resource
import sys

import numpy as np
from threadpoolctl import threadpool_limits

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]
PUBLIC = STUDY / 'physical_revision/array/observations'
KERNEL = STUDY / 'synthetic/public/fine_geometry.npz'
HISTORICAL = STUDY.parent / 'z_alg_xudong_cie/minpeak_native'
sys.path.insert(0, str(HISTORICAL))
from gn_recovery import gauss_newton
from search import Model, N, SIDE, box_fit, clocks, elapsed
from curve_inputs import Observation

CAPS = (1., 1.25, 1.5, 2., 3.)


def inputs():
    with np.load(KERNEL, allow_pickle=False) as f:
        spectrum = f['spectrum']
    with np.load(PUBLIC / 'arrays.npz', allow_pickle=False) as f:
        sensor, incident = f['S32'], f['B32']
    rows = json.loads((PUBLIC / 'manifest.json').read_text())['observations']
    cases = {}
    for row in rows:
        if (row['sources'], row['receivers'], row['noise_pct']) != (32, 32, 2):
            continue
        with np.load(PUBLIC / row['file'], allow_pickle=False) as f:
            y = f['Y']
        case = Observation(row['id'], row['object_id'], row['noise_pct'],
                           row['sample_radius'], row['padded_radius'], y)
        cases[case.object_id] = case
    return (spectrum, sensor, incident), cases


def stem(obs, cap):
    return f'{obs.id}__u{cap:g}'


def write_json(path, value):
    temp = path.with_suffix('.part')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def search(geometry, obs, cap, args):
    start = clocks()
    prefix = args.out / stem(obs, cap)
    model = Model(*geometry, obs.y)
    history, best = [], None
    target = .998*obs.padded_radius

    def record(state, row, stage):
        nonlocal best
        if best is None or state['residual'] < best['residual']:
            best = state
        history.append(dict(stage=stage, **row))
        print(obs.id, cap, stage, f"r/radius={state['residual']/obs.padded_radius:.8g}",
              f"wall={row['wall_seconds']:.5g}", flush=True)

    def gn(material, side, steps, stage):
        try:
            state, row = gauss_newton(model, material, side, np.arange(side**2), cap,
                                       steps, obs.padded_radius)
        except ArithmeticError as exc:
            history.append(dict(stage=stage, stop=str(exc)))
            return
        record(state, row, stage)

    yy, xx = (np.indices((SIDE, SIDE)) + .5) / (SIDE/2) - 1.
    generic_start = np.exp(-.5*((xx/.38)**2+(yy/.38)**2)).ravel()
    zero = np.zeros(N)
    # Two deterministic seeds use no other acquisition, label, or target geometry.
    for label, m in [('zero', zero), ('generic_gaussian', .65*cap*generic_start)]:
        gn(m, args.coarse_side, args.coarse_steps, label)
        if best is not None and best['residual'] <= target:
            break
    if best is None:
        raise ArithmeticError('Both starts failed before the first physical evaluation.')
    if best['residual'] > target and args.middle_steps:
        gn(best['m'], args.middle_side, args.middle_steps, 'middle_gauss_newton')
    if best['residual'] > target and args.full_calls:
        try:
            state, row = box_fit(model, best['m'], np.arange(N), SIDE, cap, target,
                                 args.full_calls, args.full_cpu, True)
        except ArithmeticError as exc:
            history.append(dict(stage='full_9216_box_lbfgsb', stop=str(exc)))
        else:
            record(state, row, 'full_9216_box_lbfgsb')
    material = best['m']
    if not np.all(np.isfinite(material)) or material.min() < 0 or material.max() > cap:
        raise ArithmeticError('Final material failed the native box check.')
    # Warm start, but a new solve independently checks all frozen receive/source ports.
    from search import RTOL as old_rtol
    import search as history_search
    history_search.RTOL = 1.e-10
    tightened_solve_failure = None
    try:
        numeric_residual, ports, currents = model.full(material, (best['j'],), gradient=False)
    except ArithmeticError as exc:
        tightened_solve_failure = str(exc)
        numeric_residual, ports, currents = best['residual'], best['ports'], best['j']
    finally:
        history_search.RTOL = old_rtol
    defect = material[:,None]*(model.b+model.green(currents))-currents
    arrays = dict(material=material, Y=obs.y, Y_prediction=ports, current_approx=currents)
    np.savez_compressed(prefix.with_suffix(prefix.suffix + '.npz'), **arrays)
    result = dict(observation=obs.id, cap=cap, radius=obs.padded_radius,
                  status='numerical_interior_pending_enclosure' if numeric_residual < obs.padded_radius
                  else 'numerical_upper_candidate_above_radius',
                  numerical_residual=float(numeric_residual), residual_over_radius=float(numeric_residual/obs.padded_radius),
                  raw_search_best_residual=best['residual'],
                  tighter_native_recheck_failure=tightened_solve_failure, native_peak=float(material.max()),
                  positive_cells=int(np.count_nonzero(material)),
                  state_defect_numeric_frobenius=float(np.linalg.norm(defect)),
                  optimizer=dict(coarse_side=args.coarse_side, coarse_steps=args.coarse_steps,
                                 middle_side=args.middle_side, middle_steps=args.middle_steps,
                                 full_calls=args.full_calls, full_cpu_seconds=args.full_cpu),
                  counts=dict(model.count), stages=history, peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  **elapsed(start))
    write_json(prefix.with_suffix(prefix.suffix + '.json'), result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--objects', nargs='+', default=['optfresh_01', 'optfresh_05'])
    parser.add_argument('--caps', type=float, nargs='+', default=CAPS, choices=CAPS)
    parser.add_argument('--coarse-side', type=int, default=12)
    parser.add_argument('--coarse-steps', type=int, default=7)
    parser.add_argument('--middle-side', type=int, default=24)
    parser.add_argument('--middle-steps', type=int, default=8)
    parser.add_argument('--full-calls', type=int, default=30)
    parser.add_argument('--full-cpu', type=float, default=60.)
    parser.add_argument('--out', type=Path, default=HERE / 'search_runs')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    geometry, cases = inputs()
    # Cover the most consequential box 1 and highest requested box first.
    for cap in args.caps:
        for name in args.objects:
            obs = cases[name]
            if not (args.out / (stem(obs, cap) + '.json')).exists():
                result = search(geometry, obs, cap, args)
                print('DONE', result['observation'], cap, result['residual_over_radius'],
                      result['wall_seconds'], flush=True)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
