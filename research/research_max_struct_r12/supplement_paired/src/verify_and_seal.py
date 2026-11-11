"""Audit saved public evidence; authorize offline scoring only by explicit flag.

Imports no solver, forward model, private data, or scoring module. The default
invocation writes nothing. No scientific quality threshold is imposed.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np


METHODS = {'g19', 'dr', 'ap', 'direct'}
TERMINAL = {'COMPLETE', 'CPU_BUDGET_REACHED', 'FAILED'}
VIEWS = ('initial', 'raw', 'corrected', 'primary', 'secondary')


def read(path):
    return json.loads(Path(path).read_text())


def arrays(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


def roundoff_equal(a, b):
    """Arithmetic consistency only: 256 float64 ulps of the value scale."""
    a, b = float(a), float(b)
    scale = max(abs(a), abs(b), np.finfo(float).tiny)
    return bool(np.isfinite(a) and np.isfinite(b)
                and abs(a-b) <= 256*np.finfo(float).eps*scale)


def prior_value(m, shape):
    """Saved original R only, without importing any physical backend.

    Same finite differences, width=19.5, smooth=tv=1e-4, epsilon=.01 as
    research_z_alg_r4_diversity.src.physics.Regularizer.evaluate.
    """
    a = np.asarray(m).reshape(shape)
    w = 1/19.5
    horizontal = -w*a[:, :-1] + w*a[:, 1:]
    vertical = -w*a[:-1, :] + w*a[1:, :]
    d = np.concatenate((horizontal.ravel(), vertical.ravel()))
    return float(1e-4*np.mean(d*d) + 1e-4*np.mean(np.sqrt(d*d+.01**2)-.01))


def audit(root, controller):
    errors, warnings, task_reports, loaded, stamped = [], [], [], {}, set()

    def require(condition, message):
        if not condition:
            errors.append(message)

    def stamp(path):
        path = Path(path).resolve()
        require(path.is_relative_to(root), 'Output stamp outside local root: '+str(path))
        if path.is_relative_to(root) and path.is_file():
            stamped.add(path)

    cfg = read(root/'configs/FROZEN_PROTOCOL.json')
    tasks = read(root/'configs/TASKS.json')['tasks']
    cases = read(root/'inputs/manifest.json')['cases']
    task_ids = [t['task_id'] for t in tasks]
    require(len(tasks) == 44 and len(set(task_ids)) == 44, 'Require exactly 44 unique frozen tasks')
    require(Counter(t['cohort'] for t in tasks) == {'A': 12, 'B': 16, 'C': 16},
            'Frozen task cohort counts must be A12/B16/C16')
    case_map = {(c['cohort'], c['case_id']): c for c in cases}
    require(len(case_map) == 31 and Counter(c['cohort'] for c in cases) == {'A': 3, 'B': 24, 'C': 4},
            'Public case denominator must remain A3/B24/C4')
    require(cfg['G19'].get('initialize_physical_current') is True
            and cfg['G19'].get('base_arm') == 'a_plus'
            and cfg['G19'].get('readout') == 'full_material_readout'
            and cfg['G19'].get('prepare') == 'prepare_full'
            and cfg['G19'].get('current_anchor') is None,
            'Frozen G19 projected-current/readout contract mismatch')

    queue = read(controller)
    require(queue.get('status') in {'COMPLETED', 'COMPLETED_WITH_ERRORS'}, 'Queue is not terminal')
    require(queue.get('total') == 44 and queue.get('pending') == 0 and not queue.get('running'),
            'Queue did not terminate all 44 tasks')
    jobs = queue.get('jobs', [])
    require(len(jobs) == 44 and {j['task_id'] for j in jobs} == set(task_ids),
            'Queue jobs must cover the exact 44 tasks without duplicates')
    job_map = {j['task_id']: j for j in jobs}
    stamp(controller)
    for name in ('EVENTS.jsonl', 'RESOURCE_PROBE.json'):
        p = controller.parent/name
        if p.exists():
            stamp(p)

    triggered, retained, trigger_details, public_cache = [], [], [], {}
    for key, case in case_map.items():
        cohort, ident = key
        public = arrays(root/case['input_file'])
        public_cache[key] = public
        m, Y, E, S = public['m_init'], public['Y'], public['E0'], public['S']
        require(m.shape == (576,) and np.isfinite(m).all(), f'{cohort}/{ident}: invalid public initializer')
        require(E.shape[0] == 576 and S.shape[1] == 576 and Y.shape == (S.shape[0], E.shape[1]),
                f'{cohort}/{ident}: source/receiver dimensions changed')
        require(case['source_count'] == E.shape[1] and case['receiver_count'] == S.shape[0],
                f'{cohort}/{ident}: public dimension metadata mismatch')
        selected_tasks = [t for t in tasks if (t['cohort'], t['case_id']) == key]
        if cohort == 'B':
            sigma = float(public['sigma'])
            delta = sigma*float(public['noise_radius'])
            residual = float(case['parent_paid_state']['physical_relative_residual'])*float(np.linalg.norm(Y))
            trigger = residual > delta
            require(case.get('parent_secondary_not_used') is True, ident+': B must use original primary')
            require(bool(case.get('trigger_from_original_primary_metadata')) == trigger,
                    ident+': saved B trigger disagrees with original primary public residual')
            require(roundoff_equal(case['raw_CN99_radius'], delta), ident+': original CN99 radius mismatch')
            require(roundoff_equal(case['original_primary_raw_residual_from_metadata'], residual),
                    ident+': stored primary raw residual mismatch')
            (triggered if trigger else retained).append(ident)
            trigger_details.append(dict(case_id=ident, original_primary_raw_residual=residual,
                                        delta_CN99=delta, restart_requested=trigger,
                                        retained_original=not trigger))
            if not trigger:
                require(not selected_tasks, ident+': inside-ball primary must be retained without restart')
        else:
            trigger = True
        if trigger:
            require(len(selected_tasks) == 4 and {t['method'] for t in selected_tasks} == METHODS,
                    f'{cohort}/{ident}: require four paired methods')
    require(len(triggered) == 4 and len(retained) == 20, 'B must contain four triggered and twenty retained originals')

    for task in tasks:
        tid = task['task_id']
        before_errors = len(errors)
        directory = root/'RESULTS'/tid
        report = dict(task_id=tid, status='EVIDENCE_UNREADABLE', terminal=False)
        try:
            require(task['method'] in METHODS, tid+': unknown method')
            require(directory.resolve().is_relative_to((root/'RESULTS').resolve()), tid+': unsafe output path')
            receipt = read(directory/'RECEIPT.json')
            report.update(status=receipt.get('status'), terminal=receipt.get('status') in TERMINAL,
                          execution_error=receipt.get('error'), terminations=receipt.get('terminations'))
            require(report['terminal'], tid+': receipt is nonterminal')
            require(receipt.get('truth_access') == 'NONE', tid+': unexpected online truth declaration')
            require(receipt.get('task') == task, tid+': receipt task contract differs from frozen task')
            require(receipt.get('case_id') == task['case_id'] and receipt.get('method') == task['method'],
                    tid+': receipt identity mismatch')
            if receipt.get('status') == 'FAILED':
                warnings.append(tid+': terminal FAILED; retained in every denominator and reported explicitly')
            job = job_map[tid]
            external = job.get('external', {})
            external_cpu = external.get('child_process_CPU_seconds')
            require(isinstance(external_cpu, (int, float)) and np.isfinite(external_cpu) and external_cpu >= 0,
                    tid+': missing external full-exit CPU record')
            require(isinstance(job.get('wall_seconds'), (int, float)) and job['wall_seconds'] >= 0,
                    tid+': missing external wall time')
            require(job.get('child_status') == receipt.get('status'), tid+': controller/child status mismatch')
            require(isinstance(job.get('returncode'), int), tid+': absent terminal process return code')
            report['external'] = external
            report['returncode'] = job.get('returncode')
            if job.get('returncode') != 0:
                warnings.append(tid+': nonzero terminal return code '+str(job.get('returncode')))
            require(isinstance(receipt.get('physics_cost'), dict), tid+': missing physical cost record')
            require(isinstance(receipt.get('whole_child_CPU_seconds'), (int, float)), tid+': missing child CPU record')
            if isinstance(external_cpu, (int, float)):
                require(external_cpu >= receipt['whole_child_CPU_seconds'], tid+': external CPU precedes receipt CPU')
            require(receipt.get('scoring_status') == 'SEALED_NOT_SCORED', tid+': online scoring status changed')

            history = read(directory/'HISTORY.json')
            materials = arrays(directory/'CANDIDATES.npz')['m']
            require(len(history) == len(materials), tid+': paid material/history length mismatch')
            public = public_cache[(task['cohort'], task['case_id'])]
            m0, Y = public['m_init'], public['Y']
            shape = tuple(map(int, public['grid_shape']))
            limits = cfg['budgets'][task['cohort']]
            cap = float(limits['whole_child'])
            proposal_cap = float(limits['proposal'])
            delta = float(receipt['delta_CN99'])
            require(receipt['CPU_cap'] == cap and receipt['phase_limits'] == limits, tid+': CPU contract mismatch')
            require(roundoff_equal(delta, float(public['sigma'])*float(public['noise_radius'])),
                    tid+': saved CN99 radius differs from original public radius')
            last_cpu = -np.inf
            corrector_deadlines = set()
            for index, row in enumerate(history):
                require(row['candidate'] == index, tid+': candidate sequence mismatch')
                require(np.asarray(materials[index]).shape == m0.shape and np.isfinite(materials[index]).all(),
                        tid+f': invalid paid material {index}')
                cpu, deadline = float(row['CPU_seconds']), float(row['deadline'])
                require(np.isfinite(cpu) and cpu >= last_cpu, tid+': nonmonotone/invalid paid CPU')
                last_cpu = cpu
                expected_eligible = cpu <= deadline and cpu <= cap
                require(row['eligible'] == expected_eligible, tid+f': eligibility mismatch at {index}')
                require(deadline <= cap and deadline > 0, tid+': phase deadline exceeds whole-child cap')
                if row['phase'] in ('initial', 'proposal'):
                    require(deadline == (cap if task['method'] == 'direct' else proposal_cap),
                            tid+': initial/proposal deadline mismatch')
                elif row['phase'] == 'direct':
                    require(task['method'] == 'direct' and deadline == cap, tid+': direct did not receive full cap')
                elif row['phase'] == 'corrector':
                    corrector_deadlines.add(deadline)
                else:
                    require(False, tid+': unknown paid phase '+str(row['phase']))
                require(isinstance(row.get('physics_cost'), dict), tid+': paid state has no operation costs')
                R = prior_value(materials[index], shape)
                require(roundoff_equal(row['R'], R), tid+f': original prior mismatch at {index}')
                require(roundoff_equal(row['Phi'], row['D']+row['R']), tid+f': Phi decomposition mismatch at {index}')
                require(roundoff_equal(row['D'], .5*row['relative_residual']**2),
                        tid+f': normalized data term mismatch at {index}')
                require(roundoff_equal(row['residual_norm'], row['relative_residual']*float(np.linalg.norm(Y))),
                        tid+f': raw residual metadata mismatch at {index}')
                require(row['noise_feasible'] == (row['residual_norm'] <= delta),
                        tid+f': CN99 qualification mismatch at {index}')
                legal = bool(np.isfinite(materials[index]).all() and np.min(materials[index]) >= .02
                             and np.max(materials[index]) <= 20.)
                require(row['legal'] == legal, tid+f': material legality metadata mismatch at {index}')
            require(len(corrector_deadlines) <= 1, tid+': corrector deadline changes within solve')
            require(receipt['eligible_states'] == sum(r['eligible'] for r in history)
                    and receipt['excluded_states'] == sum(not r['eligible'] for r in history),
                    tid+': eligible/excluded receipt counts disagree')
            require(not history or last_cpu <= receipt['whole_child_CPU_seconds'], tid+': paid CPU exceeds receipt CPU')

            eligible = [r for r in history if r['eligible']]
            proposal = [r for r in eligible if r['phase'] in ('initial', 'proposal')]
            feasible = [r for r in eligible if r['legal'] and r['noise_feasible']]
            pick = lambda rows: min(rows, key=lambda r:r['Phi'])['candidate'] if rows else None
            secondary = min(feasible, key=lambda r:(r['R'], r['Phi']))['candidate'] if feasible else None
            expected = dict(initial=0 if history else None, primary=pick(eligible), corrected=pick(eligible),
                            raw=(0 if history else None) if task['method'] == 'direct' else pick(proposal),
                            secondary=secondary)
            view_report = {}
            for view in VIEWS:
                npz, meta_path = directory/(view+'.npz'), directory/(view+'.json')
                exists = npz.is_file()
                require(exists == meta_path.is_file(), tid+'/'+view+': unpaired snapshot files')
                selected = receipt.get(view+'_candidate') if view != 'initial' else expected['initial']
                if selected is not None or view not in ('raw', 'corrected') or receipt['status'] != 'FAILED':
                    require(selected == expected[view], tid+'/'+view+': selector index mismatch')
                require(exists == (selected is not None), tid+'/'+view+': missing/unexpected selected snapshot')
                if not exists:
                    view_report[view] = dict(available=False, candidate=None)
                    continue
                meta, saved = read(meta_path), arrays(npz)
                index = meta['candidate']
                require(index == selected and 0 <= index < len(history), tid+'/'+view+': candidate metadata mismatch')
                if not 0 <= index < len(history):
                    continue
                require(meta == history[index], tid+'/'+view+': selected metadata differs from paid row')
                require(np.array_equal(saved['m'], materials[index]), tid+'/'+view+': selected material differs from paid pool')
                require(saved['j'].shape == public['E0'].shape and saved['Q'].shape == public['E0'].shape
                        and saved['Y_prediction'].shape == Y.shape, tid+'/'+view+': saved fields lost source dimensions')
                residual = saved['Y_prediction']-Y
                D = .5*float(np.vdot(residual, residual).real/np.vdot(Y,Y).real)
                require(roundoff_equal(meta['D'], D), tid+'/'+view+': saved prediction data term differs')
                constitutive = float(np.linalg.norm(saved['j']-saved['m'][:,None]*saved['Q'])
                                     /max(np.linalg.norm(saved['j']), np.finfo(float).tiny))
                if not np.isfinite(constitutive) or constitutive > 1e-10:
                    warnings.append(tid+'/'+view+': saved constitutive relative defect '+str(constitutive)
                                    +' exceeds diagnostic 1e-10; report only, no scientific rejection')
                view_report[view] = dict(available=True, candidate=index,
                                         constitutive_relative_defect=constitutive if np.isfinite(constitutive) else None)
                if view == 'initial':
                    require(np.array_equal(saved['m'], m0), tid+': initial material bytes differ from public input')
                    loaded[tid] = saved
            report['views'] = view_report
            report['paid_states'] = len(history)
            report['excluded_atomic_states'] = [dict(candidate=r['candidate'], CPU_seconds=r['CPU_seconds'],
                                                     deadline=r['deadline'], phase=r['phase'])
                                                for r in history if not r['eligible']]
            prefixes = read(directory/'PREFIXES.json')
            require([p['CPU_cap'] for p in prefixes] == cfg['prefixes'][task['cohort']], tid+': missing/reordered CPU prefixes')
            for prefix in prefixes:
                rows = [r for r in eligible if r['CPU_seconds'] <= prefix['CPU_cap']]
                require(prefix['candidate'] == pick(rows), tid+': prefix does not choose lowest paid Phi')
            if task['method'] == 'g19':
                execution, defect = receipt.get('G19_execution', {}), receipt.get('initial_current_projection')
                # An earlier atomic setup failure may have no completed projection;
                # preserve that explicit terminal failure instead of inventing it.
                completed_ls = execution.get('operation_counts', {}).get('initial_current_projection_SVD', 0) > 0
                if completed_ls or execution.get('status') == 'COMPLETED':
                    require(isinstance(defect, dict) and all(k in defect for k in ('relative_current_loss', 'rank', 'constitutive_defect')),
                            tid+': G19 projected-current rank/loss report missing')
                elif defect is None:
                    warnings.append(tid+': no completed initial LS projection before terminal setup/cap exit; rank/loss unavailable')
                require(execution.get('status') != 'RUNNING', tid+': capped G19 execution still claims RUNNING')
                if execution:
                    require(execution.get('sourcewise_cg') is False, tid+': G19 must retain global PR')
                report['G19_initial_projection'] = defect
            elif task['method'] in ('ap', 'dr'):
                require(isinstance(receipt.get('projection_cost'), dict) or receipt['status'] == 'FAILED',
                        tid+': missing fixed-projector operation costs')
            for p in directory.rglob('*'):
                if p.is_file():
                    stamp(p)
        except Exception as exc:
            errors.append(tid+': audit could not read complete terminal evidence: '+type(exc).__name__+': '+str(exc))
        report['contract_errors'] = errors[before_errors:]
        task_reports.append(report)

    initial_comparisons = []
    for key in case_map:
        cohort, ident = key
        ap_id, dr_id = f'{cohort}__{ident}__ap', f'{cohort}__{ident}__dr'
        if ap_id in loaded and dr_id in loaded:
            equal_fields = {}
            for field in ('m', 'j', 'Q', 'Y_prediction'):
                equal_fields[field] = bool(np.array_equal(loaded[ap_id][field], loaded[dr_id][field]))
                require(equal_fields[field],
                        f'{cohort}/{ident}: AP/DR initial physical '+field+' differ')
            initial_comparisons.append(dict(cohort=cohort, case_id=ident, fields_exactly_equal=equal_fields))
        elif any(t['cohort'] == cohort and t['case_id'] == ident for t in tasks):
            initial_comparisons.append(dict(cohort=cohort, case_id=ident,
                                            status='INITIAL_FIELDS_UNAVAILABLE_IN_TERMINAL_EVIDENCE',
                                            AP_available=ap_id in loaded, DR_available=dr_id in loaded))

    files = []
    for path in sorted(stamped):
        stat = path.stat()
        files.append(dict(path=str(path.relative_to(root)), size=stat.st_size, mtime_ns=stat.st_mtime_ns))
    report = dict(status='PUBLIC_CONTRACT_AUDIT_PASSED' if not errors else 'PUBLIC_CONTRACT_AUDIT_FAILED',
                  utc=datetime.now(timezone.utc).isoformat(), root=str(root), truth_access='NONE',
                  new_forward=0, new_adjoint=0, new_optimization=0, new_training=0,
                  tasks_expected=44, tasks_read=len(task_reports), terminal_status_counts=dict(Counter(r['status'] for r in task_reports)),
                  queue_status=queue.get('status'), terminal_failed_tasks=[r['task_id'] for r in task_reports if r['status']=='FAILED'],
                  B=dict(total=24, triggered=triggered, retained_original=retained, public_trigger_details=trigger_details),
                  errors=errors, warnings=warnings, task_reports=task_reports,
                  AP_DR_initial_exact_comparisons=initial_comparisons,
                  constitutive_check='Saved j=mQ only, relative diagnostic 1e-10; exceedances reported without quality gating; no G/S action or full forward',
                  arithmetic_check='Saved public Phi/D/R metadata uses 256 float64 ulps for arithmetic-order roundoff only',
                  timing_scope='Complete paid states checked against recorded absolute phase/whole-child deadlines; external full-exit costs retained, nested timers never added',
                  secondary_fallback='Offline scoring derives secondary_deploy from sealed secondary availability, otherwise primary; no fallback prediction is synthesized or sealed here',
                  files=files)
    return report, task_ids


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--controller', type=Path, help='Defaults to ROOT/queues/main/RECEIPT.json')
    parser.add_argument('--authorize-scoring', action='store_true', help='Explicit root authorization; write acceptance and, only on success, complete scoring seal')
    args = parser.parse_args()
    root = args.root.resolve()
    controller = (args.controller or root/'queues/main/RECEIPT.json').resolve()
    try:
        report, task_ids = audit(root, controller)
    except Exception as exc:
        report = dict(status='PUBLIC_CONTRACT_AUDIT_FAILED', errors=[type(exc).__name__+': '+str(exc)],
                      truth_access='NONE', new_forward=0, new_adjoint=0, new_optimization=0)
        task_ids = []
    if args.authorize_scoring:
        acceptance = root/'audit/ACCEPTANCE_PUBLIC.json'
        seal_path = root/'ALL44_PREDICTIONS_SEALED.json'
        if acceptance.exists() or seal_path.exists():
            raise SystemExit('Preserve existing acceptance/seal; refusing overwrite')
        acceptance.parent.mkdir(parents=True, exist_ok=True)
        acceptance.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        if not report.get('errors'):
            seal = dict(status='ALL44_PREDICTIONS_SEALED', utc=datetime.now(timezone.utc).isoformat(),
                        task_ids=task_ids, tasks=44, root_authorized_offline_scoring=True,
                        acceptance=str(acceptance.relative_to(root)), truth_access='NONE',
                        terminal_failed_tasks=report['terminal_failed_tasks'], files=report['files'],
                        scope='Public contract acceptance only; no algorithm-quality gate or scientific success claim',
                        fallback='Derived only during offline scoring; no fallback output sealed')
            seal_path.write_text(json.dumps(seal, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if not report.get('errors') else 1


if __name__ == '__main__':
    sys.exit(main())
