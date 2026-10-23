"""ORACLE known-reference fixed-edit compensation; not online inference.

Run only through the existing serial physical queue. Import performs no solves.
The native chart changes coefficient representation, not support/noise/budgets.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import time
import uuid

import numpy as np
from p0.geometry import probes

from .accounting import ROOT, ScopedBook, BudgetStop, write
from .native_world import LinearWorld
from .physical_witness import ModulatedWorld, mesh_for
from .physics import Forward, pack, packj
from .optimization import trust_region


IDENTITY = 'ORACLE_PHYSICAL_EXISTENCE_KNOWN_REFERENCE_NATIVE_FIXED_EDIT_NOT_ONLINE'


def _ledger(book):
    if not book.path.exists():
        return []
    with book.path.open() as stream:
        fcntl.flock(stream, fcntl.LOCK_SH)
        return [json.loads(line) for line in stream if line.strip()]


def _cost(book, offset, started_cpu, started_wall):
    rows = [r for r in _ledger(book)[offset:] if r.get('scope') == book.scope]
    reservations = [r for r in rows if r.get('event') == 'reserve']
    return dict(scope=book.scope, new_charged_complex_rhs=sum(r.get('rhs', 0) for r in reservations),
                new_charged_factorizations=sum(r.get('factorizations', 0) for r in reservations),
                confirmed_backend_rhs=sum(r.get('rhs', 0) for r in rows if r.get('event') == 'backend_solve'),
                reserve_ids=[r['id'] for r in reservations],
                cpu_seconds=time.process_time()-started_cpu,
                wall_seconds=time.perf_counter()-started_wall,
                cpu_billing='measured here; external campaign launcher bills process CPU once',
                ledger=str(book.path))


def _stationarity(x, gradient):
    if gradient is None:
        return None
    return float(np.linalg.norm(x-np.clip(x-gradient, -1., 1.), ord=np.inf))


def main(seed, band='narrow', kind='joint'):
    """Fixed +.1 x probe translation or fixed (-.85,.85), r=.06 air cavity.

    Both FULL_J and fixed-seed 0.1-scale RANDOM_MULTISTART use 14 TR iterations,
    34 evaluator calls, the same Forward exact-current cache and [-1,1] box.
    Existing joint/cavity named allowances are reused; no budget is reset.
    """
    seed = int(seed)
    if kind not in ('joint', 'cavity'):
        raise ValueError('kind must be joint or cavity')
    if band not in ('narrow', 'wide'):
        raise ValueError('band must be narrow or wide')
    book = ScopedBook('B2', f'{kind}_{seed}_{band}')
    cfg = json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())
    frequencies = cfg['frequencies'][band]
    if cfg['noise']['each_real_quadrature_std'] != .01:
        raise RuntimeError('shared whitening differs from unchanged original noise')
    folder = f'B2/native_{kind}_direct_{seed}_{band}'
    out = ROOT/'results'/folder
    out.mkdir(parents=True, exist_ok=False)
    original_probe = probes()[8]
    if kind == 'joint':
        edited_probe = ((original_probe[0][0]+.1, original_probe[0][1]),
                        original_probe[1], original_probe[2])
        cavities = ()
    else:
        edited_probe = None
        cavities = (((-.85, .85), .06),)
    model = LinearWorld(seed, edited_probe, cavities=cavities)
    if model.n != 108:
        raise RuntimeError('native fixed chart must have 108 parameters')
    reference_model = ModulatedWorld(seed, original_probe if kind == 'joint' else None)
    zero = np.zeros(model.n)
    reference_zero = np.zeros(reference_model.n)
    grid = np.stack(np.meshgrid(np.linspace(-2., 2., 201),
                                np.linspace(0., 3., 151)), axis=-1).reshape(-1, 2)
    reference_material = reference_model.coefficient(reference_zero, 1.)(grid)
    raw_material = model.coefficient(zero, 1.)(grid)
    run_id = uuid.uuid4().hex
    cpu_start, wall_start = time.process_time(), time.perf_counter()
    ledger_start = len(_ledger(book))
    rows, failure = [], None
    d0 = raw = None
    active_method = 'REFERENCE_AND_RAW'
    current_x, accepted_x = zero.copy(), zero.copy()
    checkpoint = None
    latest_fj = None
    fw = None
    common_cost = None
    write(out/'PLAN.json', dict(seed=seed, band=band, kind=kind, identity=IDENTITY,
          scope=book.scope, scope_rhs_cap=book.scope_caps['rhs'], run_id=run_id,
          iterations_per_method=14, evaluator_calls_per_method=34, parameter_box=[-1., 1.],
          methods=['FULL_J', 'RANDOM_MULTISTART'], random_seed=2026102321+seed,
          random_scale=.1, chart='LinearWorld 54 hats / 108 bounded affine parameters',
          reference_probe=original_probe if kind == 'joint' else None,
          fixed_candidate_probe=edited_probe, fixed_cavities=cavities,
          fixed_edit_cannot_be_removed_by_corrector=True, gain=1.,
          original_source_noise_material_prior_unchanged=True,
          probabilities=None))

    def save_current(status, error=None, converged=False):
        np.savez_compressed(out/'CURRENT_ITERATE.npz', current_trial_x=current_x,
                            accepted_x=accepted_x,
                            accepted_anchor_available=checkpoint is not None,
                            material_grid=grid,
                            current_trial_material=model.coefficient(current_x, 1.)(grid),
                            accepted_material=model.coefficient(accepted_x, 1.)(grid),
                            accepted_data=(checkpoint['data'] if checkpoint is not None
                                           else np.empty((0, len(cfg['tx']), len(cfg['rx'])), complex)))
        write(out/'CURRENT_STATUS.json', dict(method=active_method, status=status,
              error=error, accepted_anchor_available=checkpoint is not None,
              stationarity=(_stationarity(checkpoint['x'], checkpoint['gradient'])
                            if checkpoint is not None else None), converged=bool(converged)))

    try:
        # In the cavity branch, this zero-reference mesh includes the same empty
        # circle interface but the reference MATERIAL has no added cavity/probe.
        base_mesh = mesh_for(seed, original_probe) if kind == 'joint' else mesh_for(seed, None, cavities=cavities)
        base = Forward(reference_model, base_mesh, frequencies, book, folder+'/base',
                       phase=f'B2_native_reference_{run_id}', identity=IDENTITY,
                       h=.145, pml_strength=6)
        save_current('REFERENCE_PENDING')
        d0, _ = base.evaluate(reference_zero, jac=False)
        np.savez_compressed(out/'REFERENCE.npz', x=reference_zero, data=d0,
                            material_grid=grid, material=reference_material,
                            frequencies=frequencies)
        fw = Forward(model, mesh_for(seed, edited_probe, cavities=cavities), frequencies,
                     book, folder+'/edit', phase=f'B2_native_compensation_{run_id}',
                     identity=IDENTITY, h=.145, pml_strength=6, retain_current=True)
        raw, _ = fw.evaluate(zero, jac=False)
        np.savez_compressed(out/'RAW.npz', x=zero, data=raw, reference_data=d0,
                            material_grid=grid, material=raw_material,
                            reference_material=reference_material, frequencies=frequencies)
        common_cost = _cost(book, ledger_start, cpu_start, wall_start)
        write(out/'COMMON_SETUP_COST.json', common_cost)
        rng = np.random.default_rng(2026102321+seed)
        random_initial = rng.normal(size=model.n)*.1
        if np.any(abs(random_initial)>1):
            raise RuntimeError('frozen random initialization unexpectedly outside native box')
        for method, x0 in [('FULL_J', zero), ('RANDOM_MULTISTART', random_initial)]:
            active_method = method
            current_x = accepted_x = x0.copy()
            checkpoint, latest_fj = None, None
            history = []
            method_ledger = len(_ledger(book))
            method_cpu, method_wall = time.process_time(), time.perf_counter()
            method_raw = None

            def evaluate(x, jac=False):
                nonlocal current_x, latest_fj, method_raw
                current_x = np.asarray(x).copy()
                data, j = fw.evaluate(x, jac=jac)
                if method_raw is None and np.array_equal(x, x0):
                    method_raw = data.copy()
                    np.savez_compressed(out/f'{method}_RAW.npz', x=x0, data=data,
                                        reference_data=d0, material_grid=grid,
                                        material=model.coefficient(x0, 1.)(grid),
                                        reference_material=reference_material, frequencies=frequencies)
                residual = pack(data-d0)
                jacobian = None if j is None else packj(j)
                if jac:
                    latest_fj = dict(x=np.asarray(x).copy(), data=data.copy(),
                                     gradient=jacobian.T@residual)
                return residual, jacobian

            def callback(event):
                nonlocal current_x, accepted_x, checkpoint
                history.append(event)
                if 'x' in event:
                    current_x = np.asarray(event['x']).copy()
                if event.get('event') == 'anchor':
                    accepted_x = np.asarray(event['x']).copy()
                    if latest_fj is None or not np.array_equal(latest_fj['x'], accepted_x):
                        raise RuntimeError('accepted anchor missing exact full-FJ checkpoint')
                    checkpoint = {k:v.copy() for k,v in latest_fj.items()}
                    np.savez_compressed(out/f'{method}_ACCEPTED.npz', **checkpoint,
                                        material_grid=grid,
                                        material=model.coefficient(accepted_x, 1.)(grid))
                write(out/f'{method}_ITERATIONS.json', history)
                # Persist latent current/accepted states without any recovery solve.
                np.savez_compressed(out/'CURRENT_LATENT.npz', current_trial_x=current_x,
                                    accepted_x=accepted_x,
                                    accepted_anchor_available=checkpoint is not None)

            write(out/f'{method}_ITERATIONS.json', history)
            try:
                fit = trust_region(evaluate, x0, model.lower, model.upper,
                                   prior_lambda=0., max_iter=14, max_calls=34,
                                   radius=.2, max_radius=.7, l2_radius=False,
                                   callback=callback, fatal_exceptions=(Exception,))
                if checkpoint is None or fit.residual is None:
                    raise RuntimeError('optimizer returned without a complete accepted anchor')
                if not np.array_equal(checkpoint['x'], fit.x):
                    raise RuntimeError('returned iterate differs from accepted checkpoint')
                final_data = checkpoint['data']
                row = dict(seed=seed, band=band, kind=kind, method=method,
                           status='NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION',
                           converged=bool(fit.success), stopping=str(fit.message),
                           iterations=int(fit.nit), nfev=int(fit.nfev), njev=int(fit.njev),
                           projected_gradient=float(fit.projected_gradient),
                           objective=float(fit.fun),
                           raw_distance=float(np.linalg.norm(raw-d0)/.01),
                           method_initial_distance=float(np.linalg.norm(method_raw-d0)/.01),
                           final_distance=float(np.linalg.norm(final_data-d0)/.01))
                np.savez_compressed(out/f'{method}_FINAL.npz', x=fit.x, initial=x0,
                                    data=final_data, reference_data=d0, raw_data=raw,
                                    method_initial_data=method_raw,
                                    material_grid=grid, material=model.coefficient(fit.x, 1.)(grid),
                                    raw_material=raw_material, reference_material=reference_material,
                                    frequencies=frequencies, gradient=checkpoint['gradient'])
                save_current(row['status'], converged=fit.success)
            except Exception as exc:
                failure = dict(method=method, error_type=type(exc).__name__, error=str(exc))
                row = dict(seed=seed, band=band, kind=kind,
                           status='INCOMPLETE_BUDGET_STOP' if isinstance(exc, BudgetStop)
                           else 'INCOMPLETE_NUMERICAL_OR_INTERFACE_FAILURE',
                           converged=False, **failure,
                           projected_gradient=(_stationarity(checkpoint['x'], checkpoint['gradient'])
                                               if checkpoint is not None else None))
                save_current(row['status'], failure)
            row.update(identity=IDENTITY, target_shift=.1 if kind == 'joint' else 0.,
                       target_mass_relative_change=0. if kind == 'joint' else None,
                       fixed_edit='unchanged-radius/contrast target translation' if kind == 'joint'
                       else 'strict-air cavity centered (-.85,.85), diameter .12',
                       edit_cannot_be_erased_by_latent_corrector=True,
                       same_source_noise_prior_bounds=True, gain=1.,
                       actual_incremental_cost=_cost(book, method_ledger, method_cpu, method_wall))
            rows.append(row)
            write(out/f'{method}_RESULTS.json', row)
            write(out/'RESULTS.json', rows)
            print(json.dumps(row), flush=True)
            if failure is not None:
                break
    except Exception as exc:
        failure = dict(method=active_method, error_type=type(exc).__name__, error=str(exc))
        save_current('INCOMPLETE_BUDGET_STOP' if isinstance(exc, BudgetStop)
                     else 'INCOMPLETE_NUMERICAL_OR_INTERFACE_FAILURE', failure)
    finally:
        costs = _cost(book, ledger_start, cpu_start, wall_start)
        write(out/'RHS_CPU_WALL_RECEIPT.json', costs)
        summary = dict(seed=seed, band=band, kind=kind, identity=IDENTITY,
                       status='COMPLETED_TWO_NUMERICAL_METHODS' if len(rows)==2 and failure is None
                       else 'INCOMPLETE_NO_CONVERGENCE_CLAIM',
                       completed_methods=[r['method'] for r in rows if r.get('status','').startswith('NUMERICAL_METHOD_COMPLETED')],
                       not_run_methods=[m for m in ('FULL_J','RANDOM_MULTISTART')
                                        if m not in [r['method'] for r in rows]],
                       failure=failure, common_setup_cost=common_cost,
                       actual_total_cost=costs, B2_totals=book.totals()['B2'],
                       probabilities=None, scientific_status='REQUIRES_INDEPENDENT_NUMERICAL_AND_GEOMETRY_VALIDATION')
        write(out/'RESULTS.json', rows)
        write(out/'SUMMARY.json', summary)
        write(out/('COMPLETE.json' if failure is None and len(rows)==2 else 'INCOMPLETE.json'), summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed', type=int)
    parser.add_argument('band', nargs='?', default='narrow', choices=['narrow','wide'])
    parser.add_argument('kind', nargs='?', default='joint', choices=['joint','cavity'])
    args = parser.parse_args()
    result = main(args.seed, args.band, args.kind)
    if result['status'] != 'COMPLETED_TWO_NUMERICAL_METHODS':
        raise SystemExit(2)
