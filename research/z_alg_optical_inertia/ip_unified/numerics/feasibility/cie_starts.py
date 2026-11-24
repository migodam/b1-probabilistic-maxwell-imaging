"""Replay the historical Dual-CIE as a public-data initializer for original box forwards."""
import argparse
import resource

import numpy as np
from threadpoolctl import threadpool_limits

from search_native import CAPS, HERE, N, SIDE, box_fit, clocks, elapsed, inputs, stem, write_json
from cie_seed import Adapter, legacy, run_cie_arm


def material_proposals(geometry, obs, args):
    started = clocks()
    model = Adapter(*geometry, obs.y, args.seed_wall)
    setup = dict(public_schedule=dict(modes=[5, 8, 12], betas=[8., 2., .5], rank=args.rank))
    selected = []
    stop = 'historical CIE schedule completed'
    try:
        for material, step, stage, stage_end in run_cie_arm(
                model, obs.y, setup, 'dual', np.zeros(N), obs.padded_radius):
            if step in (1, 16, 32) or stage_end:
                selected.append((step, stage, material.copy()))
            if step >= args.seed_steps:
                stop = 'CIE proposal step limit'
                break
    except legacy.TimeLimit:
        stop = 'CIE proposal wall limit'
    shared = dict(schedule=setup, seed_results=[dict(step=s, stage=k, peak=float(m.max()))
                                                   for s,k,m in selected],
                  stop=stop, auxiliary_counts=dict(model.clock.counts),
                  proposal_counts=dict(model.native.count), route=model.trace_route,
                  **elapsed(started))
    write_json(args.out/(obs.id+'__cie_setup.json'), shared)
    if selected:
        return selected[-2:], model.native, shared
    return [(0,0,np.zeros(N))], model.native, shared


def one(geometry, obs, args):
    seeds, native, shared = material_proposals(geometry, obs, args)
    for cap in args.caps:
        started = clocks()
        native.count.clear()  # Final-path counts exclude shared seed generation.
        path = []
        best_residual, best_m = np.inf, None
        for step, stage, proposal in seeds:
            m = np.clip(proposal, 0., cap)
            row = dict(stage='original_forward_of_box_projected_CIE', cie_step=step, cie_stage=stage,
                       auxiliary_peak=float(proposal.max()), peak=float(m.max()))
            try:
                res, ports, j = native.full(m, gradient=False)
            except ArithmeticError as exp:
                row['native_failure'] = str(exp)
            else:
                row['native_residual'] = float(res)
                if res < best_residual:
                    best_residual, best_m, final_ports, final_j = res, m, ports, j
            path.append(row)
        if best_m is None:
            best_m = np.zeros(N)
            best_residual, final_ports, final_j = native.full(best_m, gradient=False)
        target = obs.padded_radius*.998
        if args.full_calls and best_residual > target:
            try:
                state, row = box_fit(native, best_m, np.arange(N), SIDE, cap, target,
                                     args.full_calls, 15., True)
            except ArithmeticError as exc:
                path.append(dict(stage='native_box_corrector', stop=str(exc)))
            else:
                path.append(dict(stage='native_box_corrector', **row))
                if state['residual'] < best_residual:
                    best_m, best_residual, final_ports, final_j = (state[k] for k in ('m','residual','ports','j'))
        if best_m.min() < 0 or best_m.max() > cap or not np.all(np.isfinite(best_m)):
            raise ArithmeticError('CIE physical output violates the material box.')
        import search as old_search
        old_rtol, old_search.RTOL = old_search.RTOL, 1.e-10
        tightened_solve_failure = None
        try:
            res, final_ports, final_j = native.full(best_m, (final_j,), gradient=False)
        except ArithmeticError as exc:
            res = best_residual
            tightened_solve_failure = str(exc)
        finally:
            old_search.RTOL = old_rtol
        path_base = args.out/stem(obs, cap)
        ending = path_base.suffix
        np.savez_compressed(path_base.with_suffix(ending+'.npz'), material=best_m,
                            Y=obs.y, Y_prediction=final_ports, current_approx=final_j)
        defect = best_m[:,None]*(native.b+native.green(final_j))-final_j
        receipt = dict(observation=obs.id, cap=cap, radius=obs.padded_radius,
                       status='numerical_interior_pending_enclosure' if res < obs.padded_radius else
                       'numerical_upper_candidate_above_radius', numerical_residual=float(res),
                       residual_over_radius=float(res/obs.padded_radius), native_peak=float(best_m.max()),
                       positive_cells=int(np.count_nonzero(best_m)), stages=path,
                       tighter_native_recheck_failure=tightened_solve_failure,
                       state_defect_numeric_frobenius=float(np.linalg.norm(defect)), counts=dict(native.count),
                       shared_setup_file=obs.id+'__cie_setup.json',
                       peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, **elapsed(started))
        write_json(path_base.with_suffix(ending+'.json'), receipt)
        print(obs.id, cap, receipt['residual_over_radius'], shared['stop'], receipt['wall_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--objects', nargs='+', default=['optfresh_01', 'optfresh_05'])
    parser.add_argument('--caps', nargs='+', type=float, choices=CAPS, default=[1., 3.])
    parser.add_argument('--seed-wall', type=float, default=18.)
    parser.add_argument('--rank', type=int, default=20)
    parser.add_argument('--seed-steps', type=int, default=32)
    parser.add_argument('--full-calls', type=int, default=4)
    parser.add_argument('--out', type=type(HERE), default=HERE/'cie_runs')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    geom, cases = inputs()
    for name in args.objects:
        one(geom, cases[name], args)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
