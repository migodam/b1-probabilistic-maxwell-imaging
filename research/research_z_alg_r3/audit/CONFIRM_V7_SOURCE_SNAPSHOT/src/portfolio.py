"""V7 guarded portfolios around immutable V6 solvers.

One common low-material FullJ start is attempted first. The remaining outer
wall budget goes to either all-receiver variable projection or reciprocal
FullJ restarts. Selection always uses the original all-observation residual.
No forward, Jacobian, grading, or private-data implementation lives here.
"""
from collections import Counter
import time
from .solver import solve_case, load_public, json_safe


PORTFOLIO_DEFAULTS = dict(wall_seconds=3., guard_wall_seconds=.20,
                          guard_max_nfev=80, seed=2026103001,
                          start_policy='strong_constants_sobol',
                          stop_relative_residual=1e-10,
                          fullj_nstarts=1024, fullj_max_nfev=200,
                          globaldual_nstarts=1024, globaldual_maxiter=200,
                          globaldual_maxfun=1000,
                          globaldual_correction_max_nfev=120,
                          globaldual_ftol=1e-12, globaldual_gtol=1e-10,
                          fit_ftol=1e-11, fit_xtol=1e-11, fit_gtol=1e-11)


def solve_portfolio(public_npz_or_dict, method='guarded_globaldual', config=None):
    """Run a counted common guard followed by a frozen V6 continuation.

    ``wall_seconds`` is the outer total budget. A continuation receives only
    its measured remaining portion. V6 finalization and wrapper bookkeeping
    can overrun the search budget slightly; actual outer time and overrun are
    reported, never subtracted or silently hidden.
    """
    started = time.perf_counter()
    if method not in ('guarded_globaldual','guarded_reciprocal'):
        raise ValueError('Unknown guarded portfolio method')
    cfg = dict(PORTFOLIO_DEFAULTS)
    cfg.update(config or {})
    if cfg['wall_seconds'] <= 0 or cfg['guard_wall_seconds'] <= 0:
        raise ValueError('Positive total and guard wall budgets required')
    public = load_public(public_npz_or_dict)
    phases, phase_wall, errors = {}, [], []
    candidates, traces, fits, events = [], [], [], []
    counts, timings, rejections = Counter(), Counter(), Counter()
    full_lu = full_rhs = reduced_lu = reduced_rhs = 0
    global_best = float('inf')
    common = dict(seed=cfg['seed'],start_policy=cfg['start_policy'],
                  stop_relative_residual=cfg['stop_relative_residual'],
                  fit_ftol=cfg['fit_ftol'],fit_xtol=cfg['fit_xtol'],fit_gtol=cfg['fit_gtol'])

    def elapsed():
        return time.perf_counter()-started

    def dispatch(name, core_method, phase_config):
        nonlocal full_lu, full_rhs, reduced_lu, reduced_rhs, global_best
        phase_start = elapsed()
        before = dict(full_lu=full_lu,full_rhs=full_rhs,reduced_lu=reduced_lu,reduced_rhs=reduced_rhs)
        result = solve_case(public,core_method,phase_config)
        phase_end = elapsed()
        outside_core_timer = max(0.,phase_end-phase_start-result['wall_seconds'])
        phases[name] = result
        phase_wall.append(dict(phase=name,method=core_method,outer_started_seconds=phase_start,
                               outer_finished_seconds=phase_end,outer_wall_seconds=phase_end-phase_start,
                               assigned_search_wall_seconds=phase_config['wall_seconds'],
                               outside_core_timer_seconds=outside_core_timer,
                               core_reported_wall_seconds=result['wall_seconds'],status=result['status']))
        candidate_offset = len(candidates)
        for index, old in enumerate(result['candidates']):
            row = dict(old)
            row.update(phase=name,phase_candidate_index=index,
                       phase_relative_discovery_seconds=old['discovery_seconds'],
                       discovery_seconds=phase_start+outside_core_timer+old['discovery_seconds'],
                       discovery_time_convention='Conservative upper bound: all phase time outside the V6 timer is charged before discoveries.',
                       origin=name+':'+old['origin'])
            for key, offset in (('discovery_full_complex_rhs',before['full_rhs']),
                                ('discovery_full_lu',before['full_lu']),
                                ('discovery_reduced_complex_rhs',before['reduced_rhs']),
                                ('discovery_reduced_lu',before['reduced_lu'])):
                if key in row:
                    row[key] += offset
            candidates.append(row)
        for old in result['best_so_far']:
            row = dict(old)
            row.update(phase=name,seconds=phase_start+outside_core_timer+old['seconds'],origin=name+':'+old['origin'],
                       time_convention='Conservative upper bound: all phase time outside the V6 timer is charged before discoveries.')
            for key, offset in (('full_complex_rhs',before['full_rhs']),('full_lu',before['full_lu']),
                                ('reduced_complex_rhs',before['reduced_rhs']),('reduced_lu',before['reduced_lu'])):
                if key in row:
                    row[key] += offset
            if row['relative_residual'] < global_best:
                global_best = row['relative_residual']
                traces.append(row)
        for old in result['fits']:
            row = dict(old)
            row.update(phase=name,origin=name+':'+old['origin'])
            if 'started_seconds' in row:
                row['started_seconds'] += phase_start
            if isinstance(row.get('candidate'),int):
                row['candidate'] += candidate_offset
            fits.append(row)
        for old in result['algebraic_events']:
            row = dict(old)
            row['phase'] = name
            # Full unmodified phase receipts preserve nested local indices.
            row['phase_local_indices'] = True
            if 'started_seconds' in row:
                row['started_seconds'] += phase_start
            if 'seconds' in row:
                row['seconds'] += phase_start
            events.append(row)
        a = result['accounting']
        full_lu += a['full_factorizations']
        full_rhs += a['full_complex_rhs']
        reduced_lu += a['reduced_factorizations']
        reduced_rhs += a['reduced_complex_rhs']
        counts.update(a['counts'])
        timings.update(a['operation_seconds'])
        rejections.update(result.get('rejections',{}))
        if result['status']=='FAILED_PAID':
            errors.append(dict(phase=name,status=result['status'],reason=result['stop_reason']))
        return result

    remaining = max(0.,cfg['wall_seconds']-elapsed())
    guard_config = dict(common,wall_seconds=min(cfg['guard_wall_seconds'],remaining),
                        fullj_nstarts=1,fullj_max_nfev=cfg['guard_max_nfev'],parameterization='material')
    guard = dispatch('guard','fullj',guard_config)
    guard_goal = guard['top1'] is not None and guard['top1']['relative_residual'] <= cfg['stop_relative_residual']
    stop_reason = 'GUARD_DATA_RESIDUAL_GOAL' if guard_goal else None
    if not guard_goal:
        remaining = max(0.,cfg['wall_seconds']-elapsed())
        if remaining > 0:
            continuation = dict(common,wall_seconds=remaining,fullj_max_nfev=cfg['fullj_max_nfev'])
            if method=='guarded_globaldual':
                continuation.update(globaldual_nstarts=cfg['globaldual_nstarts'],
                                    globaldual_maxiter=cfg['globaldual_maxiter'],
                                    globaldual_maxfun=cfg['globaldual_maxfun'],
                                    globaldual_correction_max_nfev=cfg['globaldual_correction_max_nfev'],
                                    globaldual_ftol=cfg['globaldual_ftol'],globaldual_gtol=cfg['globaldual_gtol'],
                                    parameterization='material')
                last = dispatch('continuation','globaldual',continuation)
            else:
                continuation.update(fullj_nstarts=cfg['fullj_nstarts'],parameterization='reciprocal')
                last = dispatch('continuation','fullj',continuation)
            stop_reason = last['stop_reason']
        else:
            stop_reason = 'TOTAL_OUTER_WALL_EXHAUSTED_AFTER_GUARD'
    ranking = sorted(range(len(candidates)),key=lambda index:candidates[index]['relative_residual'])
    top1 = candidates[ranking[0]] if ranking else None
    goal = top1 is not None and top1['relative_residual'] <= cfg['stop_relative_residual']
    status = 'COMPLETED_NUMERICAL' if goal or all(p['status']=='COMPLETED_NUMERICAL' for p in phases.values()) else 'BUDGET_EXHAUSTED'
    if stop_reason=='TOTAL_OUTER_WALL_EXHAUSTED_AFTER_GUARD' and not goal:
        status = 'BUDGET_EXHAUSTED'
    if errors:
        status = 'COMPLETED_WITH_PHASE_FAILURE' if candidates else 'FAILED_PAID'
    duration = elapsed()
    return json_safe(dict(status=status,stop_reason=stop_reason,method=method,
                          solver_revision='R3_PUBLIC_V7_GUARDED_PORTFOLIO',config=cfg,
                          case_id=public.get('case_id'),material_dimension=guard['material_dimension'],
                          receiver_count=guard['receiver_count'],source_count=guard['source_count'],
                          all_materials_unknown=True,objective='unweighted_all_observation_complex_L2',
                          selection='minimum original all-observation relative residual across every paid phase; no GT',
                          candidates=candidates,ranking=ranking,top1_index=ranking[0] if ranking else None,top1=top1,
                          best_so_far=traces,fits=fits,algebraic_events=events,rejections=dict(rejections),
                          wall_seconds=duration,phase_wall=phase_wall,phase_receipts=phases,phase_failures=errors,
                          outer_budget=dict(requested_wall_seconds=cfg['wall_seconds'],actual_wall_seconds=duration,
                                            reported_overrun_seconds=max(0.,duration-cfg['wall_seconds']),
                                            convention='Continuation gets measured remaining total time; immutable V6 finalization and wrapper bookkeeping are reported at actual cost.'),
                          accounting=dict(full_factorizations=full_lu,full_complex_rhs=full_rhs,
                                          reduced_factorizations=reduced_lu,reduced_complex_rhs=reduced_rhs,
                                          counts=dict(counts),operation_seconds=dict(timings),
                                          convention='Exact sum of phase operation ledgers, including both cached-operator imports; no discarded guard/failure work.'),
                          blas_thread_pools=guard['blas_thread_pools'],
                          scope='Finite discrete scalar-TM VIE; guarded optimizer portfolio only. V6 core unchanged; V6 exposed holdouts are development for this V7 revision. No global optimality, material-accuracy or completeness claim.'))
