"""Certify native material-specific upper residual endpoints independently of the search."""
import argparse
import json
import resource
import sys
from pathlib import Path

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import lu_solve
from threadpoolctl import threadpool_limits

from search_native import HERE, KERNEL, N, SIDE, STUDY, clocks, elapsed, inputs, stem, write_json

sys.path.insert(0, str(STUDY / 'synthetic'))
sys.path.insert(0, str(STUDY / 'submission/single'))
from local_moment import ExactGeometry
from nonvacuity import ExactFourierGeometry, StableActiveModel, ball, hi, lo, norm


class Proof:
    """Exact binary spectrum, exact binary material and Arb residual arithmetic."""
    def __init__(self, spectrum, sensor, incident):
        self.dense = ExactFourierGeometry(KERNEL)
        self.dense.S, self.dense.B = sensor, incident
        self.fft = ExactGeometry(sensor, incident, spectrum, [])

    def forward(self, m, current=None):
        started = clocks()
        if not np.any(m):
            ports_shape = (self.dense.S.shape[0], self.dense.B.shape[1])
            zeros = np.zeros((N, ports_shape[1]), dtype=complex)
            return acb_mat(*ports_shape), arb(0), zeros, dict(
                active_cells=0, interval_green_columns=0, dense_lu=0, proof='vacuum identity',
                primal_defect_frobenius_upper=0., port_prediction_error_upper=0., **elapsed(started))
        model = StableActiveModel(self.dense, m)
        inv, inverse = model.inverse_bound()
        if current is None:
            current = np.zeros((N, model.B.shape[1]), dtype=np.complex128)
            current[model.active] = lu_solve(model.lu, model.m[:, None]*model.B, check_finite=False)
        discarded_vacuum_norm = float(np.linalg.norm(current[m == 0]))
        current = current.copy()
        current[m == 0] = 0  # Exact vacuum rows are J_i=0.
        jb, bb, sb = map(ball, (current[model.active], model.B, model.S))
        g = self.fft.green(ball(current))
        x = list(map(arb, map(float, model.m)))
        error = acb_mat(model.n, bb.ncols(),
                        [x[i]*(bb[i, j]+g[int(k), j])-jb[i, j]
                         for i, k in enumerate(model.active) for j in range(bb.ncols())])
        # A^-1 on active cells controls the original native system exactly, because
        # each vacuum row of the full current equation is J_i=0.
        error_upper = (norm(sb).upper()*inv*norm(error).upper()).upper()
        predicted = sb*jb
        return predicted, error_upper, current, dict(primal_defect_frobenius_upper=hi(norm(error)),
                   sensor_frobenius_upper=hi(norm(sb)), port_prediction_error_upper=hi(error_upper),
                   discarded_numeric_vacuum_current_frobenius=discarded_vacuum_norm,
                   active_cells=model.n, interval_green_columns=bb.ncols(), dense_lu=1,
                   strict_inverse=inverse, run_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   **elapsed(started))


def distance(proof_port, prediction_error, y, radius):
    nominal = norm(proof_port-ball(y))
    lower = max(0., lo(nominal.lower()-prediction_error))
    upper = hi(nominal.upper()+prediction_error)
    return dict(radius=radius, distance_lower_for_fixed_candidate=lower,
                distance_upper_for_fixed_candidate=upper, distance_upper_over_radius=upper/radius,
                status=('certified_physical_feasible_witness' if upper <= radius else
                        'certified_candidate_outside_radius' if lower > radius else
                        'candidate_enclosure_straddles_radius'))


def prove_search_results(proof, cases, args, setup):
    certified, skipped = 0, 0
    for obs in cases.values():
        for cap in args.caps:
            prefix = stem(obs, cap)
            options = []
            for folder in args.search_dir:
                item = folder / (prefix+'.json')
                if item.exists() and item.with_suffix('.npz').exists():
                    record = json.loads(item.read_text())
                    options.append((record['numerical_residual'], record['positive_cells'], item.with_suffix('.npz')))
            if not options:
                continue
            possible = [o for o in options if o[1] <= args.max_active]
            numeric_distance, _, source = min(possible or options)
            source_ref = f'{source.parent.name}/{source.name}'
            out = args.out / (prefix+'.json')
            if out.exists():
                prior = json.loads(out.read_text())
                if (prior.get('source', f'search_runs/{source.name}') == source_ref
                        and 'distance_upper_for_fixed_candidate' in prior):
                    skipped += 1
                    continue
            with np.load(source, allow_pickle=False) as f:
                m, current = f['material'], f['current_approx']
                if not np.array_equal(f['Y'], obs.y):
                    raise ArithmeticError('Candidate observation does not match the frozen input.')
            if m.shape != (N,) or not np.isrealobj(m) or not np.all(np.isfinite(m)) or m.min() < 0 or m.max() > cap:
                raise ArithmeticError('Candidate violates the unmodified material box.')
            start = clocks()
            n_active = int(np.count_nonzero(m))
            row = dict(observation=obs.id, cap=cap, peak=float(m.max()), source=source_ref,
                       numerical_residual=numeric_distance, shared_geometry_setup=setup)
            if n_active > args.max_active:
                row.update(status='proof_not_run_dense_inverse_memory_limit', active_cells=n_active,
                           max_active=args.max_active)
            else:
                try:
                    p, e, _, forward = proof.forward(m, current)
                    row.update(distance(p, e, obs.y, obs.padded_radius), forward=forward)
                except ArithmeticError as exc:
                    row.update(status='candidate_enclosure_failed', failure=str(exc))
            row.update(elapsed(start))
            write_json(out, row)
            certified += 1
            print(prefix, row['status'], row.get('distance_upper_over_radius'),
                  row['wall_seconds'], flush=True)
    return certified, skipped


CONTROL_ID = 'control_disk075__pw360_32x32__pair2pct'


def low_peak_control(proof, args, setup):
    yy, xx = (np.indices((SIDE, SIDE)) + .5) * (2/SIDE) - 1.
    m = .625 * ((xx-.19)**2+(yy+.12)**2 <= .24**2).ravel()
    p, e, current, forward = proof.forward(m)
    center = np.array([complex(v) for v in p.entries()]).reshape((32, 32))
    rng = np.random.default_rng(2026112501)
    noise = rng.standard_normal(center.shape)+1j*rng.standard_normal(center.shape)
    y = center + noise*(.02*np.linalg.norm(center)/np.linalg.norm(noise))
    sampled_radius = hi(norm(ball(y)-ball(center)))
    r = hi(arb(sampled_radius)+norm(ball(y))/10**6)
    row = dict(observation=CONTROL_ID, cap=.75, peak=float(m.max()), geometry='Frozen native 400 MHz 2 m 32×32',
               noise_seed=2026112501, sample_radius=sampled_radius, shared_geometry_setup=setup,
               nominal_noise_denominator='Frobenius norm of saved clean binary64 port center',
               observed_noise_ratio=sampled_radius/np.linalg.norm(center), forward=forward,
               generating_control=True, real_material=True, **distance(p, e, y, r))
    if row['status'] != 'certified_physical_feasible_witness':
        raise ArithmeticError('Ground-truth control failed the preregistered ball inclusion test.')
    np.savez_compressed(args.out/(CONTROL_ID+'.npz'), material=m, Y=y, Y_clean=center, current_approx=current)
    write_json(args.out/(CONTROL_ID+'.json'), row)
    print(CONTROL_ID, row['status'], row['distance_upper_over_radius'], flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--search-dir', type=Path, action='append')
    ap.add_argument('--out', type=Path, default=HERE/'strict_endpoints')
    ap.add_argument('--objects', nargs='+')
    ap.add_argument('--caps', type=float, nargs='+', default=[1., 1.25, 1.5, 2., 3.])
    ap.add_argument('--max-active', type=int, default=3000)
    ap.add_argument('--control-only', action='store_true')
    args = ap.parse_args()
    if args.search_dir is None:
        args.search_dir = [HERE/'search_runs', HERE/'projected_runs', HERE/'cie_runs']
    args.out.mkdir(parents=True, exist_ok=True)
    invocation = clocks()
    geo, cases = inputs()
    if args.objects:
        cases = {name: cases[name] for name in args.objects}
    ctx.prec, ctx.threads = 128, 1
    started = clocks()
    proof = Proof(*geo)
    setup = elapsed(started)
    low_peak_run = not (args.out/(CONTROL_ID+'.json')).exists()
    if low_peak_run:
        low_peak_control(proof, args, setup)
    tested, replay_skipped = (0, 0) if args.control_only else prove_search_results(proof, cases, args, setup)
    with (args.out/'audit_invocations.jsonl').open('a') as f:
        f.write(json.dumps(dict(cases_written=tested, existing_cases_skipped=replay_skipped,
                               separate_control_run=low_peak_run, shared_setup=setup,
                               interval_bits=ctx.prec, max_active_cells=args.max_active,
                               **elapsed(invocation)), allow_nan=False)+'\n')


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
