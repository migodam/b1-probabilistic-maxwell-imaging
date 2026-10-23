"""Measurement-only wide-band continuation driver; importing does no work.

Run only in the campaign's serial physical queue. No truth/oracle imports.
The fixed selected support is a restricted inverse class, not a geometry proof.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import math
import time
import uuid
from pathlib import Path

import numpy as np
from wq.inverse_support import PixelGrid, pixel_mesh

from .accounting import ROOT, Book, BudgetStop, write
from .online_inference import GroupPixels, born_initialize, selected_groups
from .optimization import trust_region
from .physics import Forward


STAGES = ((1, 8, .002), (2, 6, .002), (3, 6, .002), (5, 10, 0.))
IDENTITY = 'ONLINE_MEASUREMENT_ONLY_FREQUENCY_CONTINUATION'


def _ledger_rows(book):
    if not book.path.exists():
        return []
    with book.path.open() as stream:
        fcntl.flock(stream, fcntl.LOCK_SH)
        return [json.loads(line) for line in stream if line.strip()]


def _receipt(book, phases):
    """Independent phase-filtered ledger sums, never inferred from optimizer calls."""
    rows = [r for r in _ledger_rows(book)
            if r.get('stage') == 'R' and r.get('phase') in phases]
    reservations = [r for r in rows if r.get('event') == 'reserve']
    solves = [r for r in rows if r.get('event') == 'backend_solve']
    factors = [r for r in rows if r.get('event') == 'backend_factor']
    return dict(ledger=str(book.path), phases=list(phases),
                charged_complex_rhs=sum(int(r.get('rhs', 0)) for r in reservations),
                charged_factorizations=sum(int(r.get('factorizations', 0)) for r in reservations),
                confirmed_backend_solve_rhs=sum(int(r.get('rhs', 0)) for r in solves),
                confirmed_backend_factor_events=len(factors),
                reserve_ids=[r['id'] for r in reservations],
                budget_rejections=[r for r in rows if r.get('event') == 'budget_reject'],
                cost_status='real R-stage ledger reservations; confirmed solves reported separately')


def _forecast(book, frequency_count, max_calls, iterations, source_count, receiver_count):
    # Driver enforces at most 1+iterations full-J calls; every remaining
    # optimizer call is F-only. One potentially uncached final F is included.
    max_jacobian_calls = 1 + iterations
    forecast = dict(rhs=frequency_count * (max_calls * source_count
                                           + max_jacobian_calls * receiver_count
                                           + source_count),
                    factorizations=frequency_count * (max_calls + 1),
                    max_optimizer_calls=max_calls,
                    max_full_jacobian_calls=max_jacobian_calls,
                    final_forward_calls_reserved_in_forecast=1)
    totals = book.totals()['R']
    cap = book.config['stages']['R']
    available = dict(rhs=cap['rhs'] - totals['rhs'],
                     factorizations=cap['factorizations'] - totals['factorizations'],
                     cpu_seconds=cap['cpu_seconds'] - totals['cpu_seconds']
                                 - (time.process_time() - book.start))
    reasons = [key for key in ('rhs', 'factorizations')
               if forecast[key] > available[key]]
    if available['cpu_seconds'] <= 0:
        reasons.append('cpu_seconds_exhausted')
    return dict(required=forecast, available=available, admitted=not reasons,
                rejection_reasons=reasons,
                cpu_forecast='No time prediction; Book enforces remaining CPU at each reservation',
                concurrency='Requires external serial physical queue; forecast does not reserve or extend budget')


def main(seed, iteration_scale=1., folder=None):
    """Run fixed [0.8], [0.8,0.9], [0.8,0.9,1.0], all-five schedule.

    iteration_scale in (0,1] may only reduce each stage's iteration allocation.
    An unaffordable stage and every later stage are explicitly skipped as whole
    stages. No budget is created, transferred, or automatically increased.
    """
    if not np.isfinite(iteration_scale) or not 0 < iteration_scale <= 1:
        raise ValueError('iteration_scale must be finite in (0,1]')
    seed = int(seed)
    if not (ROOT / 'results/R2/COMPLETE.json').exists():
        raise RuntimeError('derivative audit missing')
    cfg = json.loads((ROOT / 'configs/ORIGINAL_P0.json').read_text())
    frequencies = cfg['frequencies']['wide']
    if not np.array_equal(np.asarray(frequencies), np.array([.8, .9, 1., 1.1, 1.2])):
        raise RuntimeError('original wide frequencies differ from frozen continuation schedule')
    if cfg['noise']['each_real_quadrature_std'] != .01:
        raise RuntimeError('original noise differs from unchanged shared whitening')
    with np.load(ROOT / f'inputs/mean_{seed}_wide.npz') as archive:
        y = archive['data'].copy()
    if y.shape != (5, len(cfg['tx']), len(cfg['rx'])) or not np.all(np.isfinite(y)):
        raise ValueError('unexpected wide measurement shape or nonfinite values')
    book = Book('R')
    run_id = uuid.uuid4().hex
    relative = Path(folder or f'R3/online_{seed}_wide_FREQUENCY_CONTINUATION')
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('folder must be a results-relative destination')
    out = ROOT / 'results' / relative
    out.mkdir(parents=True, exist_ok=False)
    grid = PixelGrid()
    # Identical measurement-only initializer/support selection to online_inference.
    born = born_initialize(y, frequencies, grid)
    groups = selected_groups(born, grid)
    if len(groups) == 0:
        raise RuntimeError('measurement-only selected support is empty')
    model = GroupPixels(grid, groups)
    initial = np.r_[model.G.T @ born[:grid.n_pixels] / 4,
                    model.G.T @ born[grid.n_pixels:] / 4]
    x = initial.copy()
    np.savez_compressed(out / 'INITIAL.npz', initial=initial, groups=groups)
    plan = [dict(stage=k, frequencies=frequencies[:count], frequency_count=count,
                 base_iterations=base, iterations=max(1, math.ceil(base * iteration_scale)),
                 prior_lambda=lam, prior_center='zero', support='same fixed measurement-selected GroupPixels')
            for k, (count, base, lam) in enumerate(STAGES, 1)]
    write(out / 'PLAN.json', dict(seed=seed, identity=IDENTITY, run_id=run_id,
          iteration_scale=iteration_scale, stages=plan,
          original_prior_and_noise='unchanged; stage 4 explicitly removes numerical weak prior',
          support_exchange=('not performed; any proposal requires a separate measurement-only, '
                            'bound-aware outer support-exchange study')))
    forecasts = [_forecast(book, st['frequency_count'], 2*st['iterations']+10,
                           st['iterations'], len(cfg['tx']), len(cfg['rx'])) for st in plan]
    sequence_required = {key: sum(fc['required'][key] for fc in forecasts)
                         for key in ('rhs', 'factorizations')}
    sequence_available = forecasts[0]['available']
    sequence_reasons = [key for key in sequence_required
                        if sequence_required[key] > sequence_available[key]]
    if sequence_available['cpu_seconds'] <= 0:
        sequence_reasons.append('cpu_seconds_exhausted')
    sequence_forecast = dict(required=sequence_required, available=sequence_available,
                             admitted=not sequence_reasons, rejection_reasons=sequence_reasons,
                             stages=forecasts, no_budget_extension=True)
    write(out / 'SEQUENCE_BUDGET_FORECAST.json', sequence_forecast)
    stage_rows, phases = [], []
    final_data, final_frequencies = np.empty((0, len(cfg['tx']), len(cfg['rx'])), complex), []
    all_wide_complete, halted, failure = False, not sequence_forecast['admitted'], None
    start_cpu = time.process_time()
    try:
        for stage in plan:
            k, count, iterations = stage['stage'], stage['frequency_count'], stage['iterations']
            stage_out = out / f'stage_{k:02d}'
            stage_out.mkdir(exist_ok=False)
            max_calls = 2 * iterations + 10
            forecast = _forecast(book, count, max_calls, iterations, len(cfg['tx']), len(cfg['rx']))
            write(stage_out / 'BUDGET_FORECAST.json', forecast)
            if halted or not forecast['admitted']:
                status = ('SKIPPED_WHOLE_SEQUENCE_INSUFFICIENT_R_BUDGET'
                          if not sequence_forecast['admitted'] else
                          'SKIPPED_AFTER_UNAFFORDABLE_OR_FAILED_STAGE' if halted
                          else 'SKIPPED_INSUFFICIENT_R_BUDGET')
                row = dict(**stage, status=status, forecast=forecast, new_rhs=0,
                           new_factorizations=0)
                write(stage_out / 'FWI_HISTORY.json', [])
                write(stage_out / 'RESULTS.json', row)
                stage_rows.append(row)
                halted = True
                continue
            phase = f'R3_frequency_continuation_{seed}_{run_id}_s{k}'
            phases.append(phase)
            history = []
            last_anchor = x.copy()
            def callback(event):
                nonlocal last_anchor
                history.append(event)
                if event.get('event') == 'anchor':
                    last_anchor = np.asarray(event['x'], dtype=float).copy()
                write(stage_out / 'FWI_HISTORY.json', history)

            write(stage_out / 'FWI_HISTORY.json', history)
            try:
                fw = Forward(model, pixel_mesh(grid), stage['frequencies'], book,
                             str(relative / f'stage_{k:02d}' / 'forward'),
                             phase=phase, identity=IDENTITY)
                residual = fw.residual(y[:count])
                full_jacobian_calls = 0

                def evaluate(point, jac=False):
                    nonlocal full_jacobian_calls
                    if jac:
                        if full_jacobian_calls >= 1 + iterations:
                            # Enforce the forecast before any physical call even
                            # if an unexpected repeated accepted-refresh occurs.
                            raise RuntimeError('stage full-J forecast cap reached')
                        full_jacobian_calls += 1
                    return residual(point, jac=jac)

                fit = trust_region(evaluate, x, np.zeros(model.n),
                                   np.ones(model.n), prior_lambda=stage['prior_lambda'],
                                   max_iter=iterations, radius=.15, max_radius=.5,
                                   l2_radius=False, max_calls=max_calls,
                                   callback=callback, fatal_exceptions=(Exception,))
                if fit.residual is None:
                    raise RuntimeError('stage has no fully evaluated accepted anchor: ' + str(fit.message))
                data, _ = fw.evaluate(fit.x, jac=False)
                x = fit.x.copy()
                final_data, final_frequencies = data.copy(), list(stage['frequencies'])
                distance = float(np.linalg.norm(data - y[:count]) / .01)
                np.savez_compressed(stage_out / 'FINAL.npz', x=x,
                                    full_parameters=model.full(x), data=data,
                                    frequencies=stage['frequencies'], initial=last_anchor
                                    if not fit.anchors else fit.anchors[0]['x'], groups=groups)
                row = dict(**stage, status='COMPLETED_NUMERICAL_STAGE',
                           success=bool(fit.success), stopping=str(fit.message),
                           objective=float(fit.fun), projected_gradient=float(fit.projected_gradient),
                           iterations_used=int(fit.nit), nfev=int(fit.nfev), njev=int(fit.njev),
                           actual_forward_endpoint_calls=fw.calls,
                           d_mean_subset=distance, rhs_receipt=_receipt(book, [phase]))
                all_wide_complete = count == 5
            except Exception as exc:
                # Preserve the last accepted stage anchor without requesting any
                # extra fields during failure/budget handling.
                np.savez_compressed(stage_out / 'RETAINED_ANCHOR.npz', x=last_anchor,
                                    full_parameters=model.full(last_anchor), groups=groups)
                row = dict(**stage, status='FAILED_STAGE_NO_AUTOMATIC_RETRY',
                           error_type=type(exc).__name__, error=str(exc),
                           rhs_receipt=_receipt(book, [phase]))
                failure, halted = row, True
            finally:
                write(stage_out / 'RHS_RECEIPT.json', _receipt(book, [phase]))
            row['forecast'] = forecast
            stage_rows.append(row)
            write(stage_out / 'RESULTS.json', row)
            write(out / 'STAGES.json', stage_rows)
            if failure is not None and failure['error_type'] == 'BudgetStop':
                # A budget stop does not trigger another physical evaluation.
                halted = True
    finally:
        receipt = _receipt(book, phases)
        write(out / 'RHS_RECEIPT.json', receipt)
        np.savez_compressed(out / 'FINAL.npz', x=x, full_parameters=model.full(x),
                            initial=initial, groups=groups, data=final_data,
                            data_frequencies=final_frequencies)
        summary = dict(seed=seed, band='wide', method='FREQUENCY_CONTINUATION_TR',
                       identity=IDENTITY, run_id=run_id,
                       status=('ALL_FOUR_NUMERICAL_STAGES_COMPLETED' if all_wide_complete
                               else 'INCOMPLETE_CONTINUATION'), stages=stage_rows,
                       all_wide_complete=all_wide_complete,
                       sequence_budget_forecast=sequence_forecast,
                       all_wide_d_mean=(float(np.linalg.norm(final_data-y)/.01)
                                        if all_wide_complete else None),
                       final_data_frequencies=final_frequencies,
                       group_count=len(groups), occupied_area_max=float(4*len(groups)*grid.cell_area),
                       geometry_status='restricted feature support; independent prospective validator required',
                       support_exchange='not performed; active support held fixed throughout all correctors',
                       rhs_receipt=receipt, R_totals=book.totals()['R'],
                       driver_cpu_seconds_unbilled=time.process_time()-start_cpu,
                       cpu_accounting='external campaign launcher must charge process CPU once',
                       failure=failure, probabilities=None)
        write(out / 'STAGES.json', stage_rows)
        write(out / 'RESULTS.json', summary)
    print(json.dumps(summary), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed', type=int)
    parser.add_argument('--iteration-scale', type=float, default=1.)
    parser.add_argument('--folder', default=None)
    args = parser.parse_args()
    result = main(args.seed, iteration_scale=args.iteration_scale, folder=args.folder)
    if not result['all_wide_complete']:
        raise SystemExit(2)
