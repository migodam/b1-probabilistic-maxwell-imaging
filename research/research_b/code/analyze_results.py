"""Parent-paired statistics for saved Research B outcomes; no physics imports.

Two noise repetitions are averaged within each of the 32 independent parents
before the fixed 2,000 parent bootstrap. E0 remains a descriptive exposed panel.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np

try:
    from .cost import Cost, ROOT, dump
except ImportError:
    from cost import Cost, ROOT, dump

BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_SEED = 26101569
BASELINES = ('fixed_2.4', 'fixed_3.2', 'none')
MODES = ('mode1_candidate_only', 'mode2_reference_selection')
REFERENCE_ACTION_POLICIES = frozenset(('event_EIG', 'expected_material_risk',
    'calibrate_then_risk', 'risk_cost_stop', 'two_step_bundle'))
METRICS = ('selected_MAE', 'best_of_K_MAE', 'Coverage_at_K_005',
    'event_decision_accuracy', 'selected_structure_accuracy', 'measurement_blocks',
    'K', 'reference_posterior_ESS', 'residual_sentinel_rate', 'candidate_profile_chi2',
    'candidate_reference_regret', 'finite_reference_risk_diagnostic',
    'acquisition_seconds_assigned', 'effective_generator_refresh_rate',
    'observed_step', 'stopped_rate', 'selected_MAE_change_from_initial')
LOWER_IS_BETTER = frozenset(('selected_MAE', 'best_of_K_MAE', 'measurement_blocks',
    'candidate_reference_regret', 'finite_reference_risk_diagnostic',
    'acquisition_seconds_assigned', 'selected_MAE_change_from_initial'))
HIGHER_IS_BETTER = frozenset(('Coverage_at_K_005', 'event_decision_accuracy',
                             'selected_structure_accuracy'))


def _read(path):
    before = path.stat()
    value = json.loads(path.read_text())
    after = path.stat()
    if (before.st_mtime_ns, before.st_size) != (after.st_mtime_ns, after.st_size):
        raise RuntimeError(f'Input changed while reading: {path.name}; wait for stage completion')
    return value


def _number(value):
    if value is None:
        return None
    result = float(value)
    if not np.isfinite(result):
        raise ValueError('Nonfinite saved metric; do not silently discard it')
    return result


def _complete_mean(values):
    return None if any(v is None for v in values) else float(np.mean(values))


def _action_source(policy):
    if policy in REFERENCE_ACTION_POLICIES:
        return 'REFERENCE_ASSISTED_ACQUISITION'
    if policy in ('raw', 'robust', 'robust_abstain', 'conflict', 'conflict_diagnostic',
                   'classical_Fisher', 'joint_projected_Fisher'):
        return 'CANDIDATE_SET_AND_LOCAL_PHYSICS_ACQUISITION'
    return 'PREDETERMINED_OR_RANDOM_ACQUISITION'


def _baseline_scope(stage, baseline):
    if baseline == 'fixed_3.2':
        return ('E1_POST_HOC_STRONG_FIXED_DESCRIPTIVE_BASELINE' if stage == 'E1' else
                'E2_FROZEN_AFTER_E1_EXPOSED_FOLLOWUP' if stage == 'E2' else 'NOT_AVAILABLE_E0')
    return 'PREDETERMINED_SECOND_FREQUENCY' if baseline == 'fixed_2.4' else 'NO_ADDITIONAL_MEASUREMENT'


def _load_main_stage(results, stage):
    status_path = results/f'{stage}_STATUS.json'
    if not status_path.exists():
        raise RuntimeError(f'{stage} completion receipt absent; do not summarize a growing campaign')
    status = _read(status_path)
    if status.get('parents') != 32 or status.get('noise_per_parent') != 2:
        raise ValueError(f'{stage} is not the frozen 32-parent/two-noise panel')
    rows = _read(results/f'{stage}_POLICY_RESULTS.json')
    episodes = defaultdict(list)
    seen = set()
    for row in rows:
        key = (row['parent'], row['noise'], row['policy'], row['step'])
        if key in seen:
            raise ValueError(f'Duplicate saved episode step: {key}')
        seen.add(key)
        episodes[key[:3]].append(row)
    expected = {(p, n, policy) for p in range(32) for n in range(2)
                for policy in status['policies']}
    if set(episodes) != expected:
        raise ValueError(f'{stage} has incomplete/extra episodes, not a complete 32-parent panel')
    for (parent, noise, policy), episode in episodes.items():
        episode.sort(key=lambda row: row['step'])
        steps = [row['step'] for row in episode]
        if steps != list(range(max(steps)+1)):
            raise ValueError(f'Noncontiguous history: {(parent,noise,policy)}')
        last = episode[-1]
        if (last['step'] < status['max_new_blocks'] and policy != 'none'
                and not last.get('stop_status')):
            raise ValueError(f'Unfinished trajectory without explicit stopping: {(parent,noise,policy)}')
    return rows, episodes, status


def _extract(row, mode, initial):
    value = row[mode]
    reference_mode = mode == MODES[1]
    return dict(selected_MAE=_number(value['selected_MAE']),
        best_of_K_MAE=_number(value['best_of_K_MAE']),
        Coverage_at_K_005=_number(value['Coverage_at_K_005']),
        event_decision_accuracy=_number(value.get('event_decision_correct')),
        selected_structure_accuracy=_number(value.get('selected_structure_correct')),
        measurement_blocks=_number(row['measurement_blocks']),
        K=_number(row.get('reference_K' if reference_mode else 'candidate_only_K')),
        reference_posterior_ESS=_number(row.get('posterior_ESS')),
        residual_sentinel_rate=_number(row.get('residual_sentinel')),
        candidate_profile_chi2=_number(row.get('candidate_profile_chi2')),
        candidate_reference_regret=None if reference_mode else _number(row.get('candidate_reference_regret')),
        finite_reference_risk_diagnostic=_number(row.get('finite_reference_risk')),
        acquisition_seconds_assigned=_number(row.get('compute_acquisition_seconds')),
        effective_generator_refresh_rate=float(bool(row.get('generator_refresh')) and reference_mode),
        observed_step=float(row['step']),
        stopped_rate=float(bool(row.get('stop_status'))),
        selected_MAE_change_from_initial=float(value['selected_MAE']-initial[mode]['selected_MAE']))


def _main_parent_rows(stage, episodes, status):
    output = []
    horizons = [('initial', 0)] + [(f'step{i}', i) for i in range(1, status['max_new_blocks']+1)] + [('final', None)]
    by_parent = defaultdict(list)
    for (parent, noise, policy), episode in sorted(episodes.items()):
        for horizon, step in horizons:
            eligible = episode if step is None else [row for row in episode if row['step'] <= step]
            selected = eligible[-1]
            for mode in MODES:
                by_parent[(mode, policy, horizon, parent)].append((noise, selected,
                    _extract(selected, mode, episode[0])))
    for (mode, policy, horizon, parent), items in sorted(by_parent.items()):
        if sorted(noise for noise, _, _ in items) != [0, 1]:
            raise ValueError('A parent mean requires both fixed noise repetitions')
        source = _action_source(policy)
        output.append(dict(stage=stage, mode=mode, policy=policy, horizon=horizon,
            parent=parent, n_noise=2,
            action_information_source=source,
            material_selection='BOUNDED_SHARED_GAIN_PROFILE_NO_CANDIDATE_POSTERIOR' if mode == MODES[0]
                else 'FINITE_REFERENCE_PARTICLE_EXPECTED_MATERIAL_LOSS',
            pure_candidate_only_pipeline=mode == MODES[0] and source != 'REFERENCE_ASSISTED_ACQUISITION',
            reference_diagnostics_scope='ESS and finite_reference_risk describe shared reference particles, not candidate probabilities',
            exposure=status.get('exposure', 'same E1 parents; exposed followup' if stage == 'E2' else 'frozen methods, fresh prior parents'),
            stop_statuses='|'.join(sorted(set(row.get('stop_status', '') for _, row, _ in items))),
            observed_noise_steps='|'.join(str(row['step']) for _, row, _ in items),
            **{metric: _complete_mean([value[metric] for _, _, value in items]) for metric in METRICS}))
    return output


def _e0_parent_rows(results):
    """Historical seven-case summaries; no claim of 32 independent parents."""
    reference = _read(results/'E0_POLICY_REALIZED.json')
    candidates = _read(results/'E0_CANDIDATE_ONLY_REALIZED.json')
    candidate_map = {(r['case'], r['policy']): r for r in candidates if r['candidate_deletion'] == 'none'}
    baseline = {r['case']: r for r in reference if r['policy'] == 'no_measurement'}
    rename = {'no_measurement': 'none', 'predetermined_second_frequency': 'fixed_2.4', 'uniform_one_block': 'uniform'}
    output = []
    for original in reference:
        key = (original['case'], original['policy'])
        policy = rename.get(original['policy'], original['policy'])
        candidate = candidate_map.get(key)
        for mode in MODES:
            ref_mode = mode == MODES[1]
            source = original if ref_mode else candidate
            if source is None:
                continue
            selected = source['selected_material_MAE'] if ref_mode else source['candidate_only_selected_MAE']
            initial = source['initial_selected_MAE'] if ref_mode else source['initial_candidate_only_MAE']
            event_correct = source.get('structure_correct') if ref_mode else source.get('selected_structure_correct')
            for horizon in ('initial', 'step1', 'final'):
                use_initial = horizon == 'initial'
                mae = initial if use_initial else selected
                best = source['best_of_K_MAE']
                accuracy = (baseline[original['case']].get('structure_correct') if ref_mode else
                    candidate_map[(original['case'], 'no_measurement')].get('selected_structure_correct')) if use_initial else event_correct
                data = {metric: None for metric in METRICS}
                data.update(selected_MAE=float(mae), best_of_K_MAE=float(best),
                    Coverage_at_K_005=float(best <= .05), event_decision_accuracy=_number(accuracy),
                    # E0's reference 'structure_correct' records event decision, not selected image event.
                    selected_structure_accuracy=None if ref_mode else _number(accuracy),
                    measurement_blocks=0. if use_initial else float(original['measurement_blocks']),
                    K=float(source['K']), effective_generator_refresh_rate=0.,
                    observed_step=0. if use_initial else float(original['measurement_blocks']),
                    stopped_rate=0., selected_MAE_change_from_initial=float(mae-initial))
                if not ref_mode and not use_initial:
                    data.update(residual_sentinel_rate=_number(source.get('residual_sentinel')),
                                candidate_profile_chi2=_number(source.get('min_joint_profile_chi2')))
                output.append(dict(stage='E0', mode=mode, policy=policy, horizon=horizon,
                    parent=original['case'], n_noise=1,
                    action_information_source=_action_source(policy),
                    material_selection='BOUNDED_SHARED_GAIN_PROFILE_NO_CANDIDATE_POSTERIOR' if not ref_mode else 'FINITE_REFERENCE_PARTICLE_EXPECTED_MATERIAL_LOSS',
                    pure_candidate_only_pipeline=not ref_mode and policy not in REFERENCE_ACTION_POLICIES,
                    reference_diagnostics_scope='historical exposed seven-case panel; absent metrics stay null',
                    exposure='EXPOSED_DEVELOPMENT_CASES_NOT_32_INDEPENDENT_PARENTS',
                    stop_statuses='', observed_noise_steps='0' if use_initial else str(original['measurement_blocks']), **data))
    return output


def _summary(parent_rows):
    groups = defaultdict(list)
    for row in parent_rows:
        groups[(row['stage'], row['mode'], row['policy'], row['horizon'])].append(row)
    rows, availability = [], []
    for key, items in sorted(groups.items()):
        stage, mode, policy, horizon = key
        complete = len(items) == 32 if stage != 'E0' else True
        if not complete:
            raise ValueError(f'Incomplete parent group: {key}')
        numeric = {metric: _complete_mean([row[metric] for row in items]) for metric in METRICS}
        selected = np.array([row['selected_MAE'] for row in items])
        row = dict(stage=stage, mode=mode, policy=policy, horizon=horizon,
            n_parents=len(items), n_noise_episodes=sum(item['n_noise'] for item in items),
            action_information_source=items[0]['action_information_source'],
            material_selection=items[0]['material_selection'],
            pure_candidate_only_pipeline=items[0]['pure_candidate_only_pipeline'],
            statistical_unit='independent_parent_mean_of_two_noise_repetitions' if stage != 'E0' else 'exposed_development_case_descriptive_only',
            exposure=items[0]['exposure'], **numeric,
            selected_MAE_SD_between_parents=float(selected.std(ddof=1)),
            selected_MAE_p95_parent=float(np.quantile(selected, .95)))
        if stage != 'E0':
            indices = np.random.default_rng(BOOTSTRAP_SEED).integers(0, len(items), size=(BOOTSTRAP_REPLICATES,len(items)))
            low, high = np.quantile(selected[indices].mean(axis=1), [.025,.975])
            row.update(selected_MAE_CI95_low=float(low), selected_MAE_CI95_high=float(high))
        else:
            row.update(selected_MAE_CI95_low=None, selected_MAE_CI95_high=None)
        rows.append(row)
        availability.append(dict(stage=stage, mode=mode, policy=policy, horizon=horizon,
            available_parents={metric: sum(item[metric] is not None for item in items) for metric in METRICS}))
    return rows, groups, availability


def _paired(groups):
    out = []
    for (stage, mode, policy, horizon), items in sorted(groups.items()):
        by_parent = {row['parent']: row for row in items}
        parent_ids = sorted(by_parent, key=lambda parent: (isinstance(parent,str), parent))
        indices = np.random.default_rng(BOOTSTRAP_SEED).integers(0, len(parent_ids),
            size=(BOOTSTRAP_REPLICATES,len(parent_ids)))
        for baseline in BASELINES:
            if policy == baseline:
                continue
            base_items = groups.get((stage, mode, baseline, horizon))
            context = dict(stage=stage, mode=mode, policy=policy, horizon=horizon,
                baseline=baseline, baseline_scope=_baseline_scope(stage,baseline),
                delta_definition='policy_minus_baseline', bootstrap_seed=BOOTSTRAP_SEED,
                bootstrap_replicates=BOOTSTRAP_REPLICATES if stage != 'E0' else 0)
            if base_items is None:
                out.append(dict(**context, metric=None, status='BASELINE_NOT_AVAILABLE'))
                continue
            base = {row['parent']: row for row in base_items}
            if set(base) != set(by_parent):
                raise ValueError(f'Parent identities differ in paired comparison: {context}')
            blocks = np.array([by_parent[p]['measurement_blocks'] for p in parent_ids])
            base_blocks = np.array([base[p]['measurement_blocks'] for p in parent_ids])
            for metric in METRICS:
                values = [by_parent[p][metric] for p in parent_ids]
                controls = [base[p][metric] for p in parent_ids]
                row = dict(**context, metric=metric, n_parents=len(parent_ids),
                    n_noise_episodes=sum(by_parent[p]['n_noise'] for p in parent_ids),
                    policy_mean_blocks=float(blocks.mean()), baseline_mean_blocks=float(base_blocks.mean()),
                    measurement_block_delta=float((blocks-base_blocks).mean()),
                    equal_blocks_for_all_parents=bool(np.array_equal(blocks,base_blocks)),
                    preferred_direction='lower' if metric in LOWER_IS_BETTER else 'higher' if metric in HIGHER_IS_BETTER else 'descriptive_only')
                if any(value is None for value in values+controls):
                    row.update(status='METRIC_UNAVAILABLE_NO_PARTIAL_PARENT_DROPPING',
                               policy_mean=None, baseline_mean=None, paired_delta=None,
                               CI95_low=None, CI95_high=None)
                else:
                    values, controls = np.array(values), np.array(controls)
                    difference = values-controls
                    row.update(status='PARENT_PAIRED_DESCRIPTIVE_BOOTSTRAP' if stage != 'E0' else 'EXPOSED_CASE_PAIRED_DESCRIPTION_ONLY',
                        policy_mean=float(values.mean()), baseline_mean=float(controls.mean()),
                        paired_delta=float(difference.mean()),
                        parent_improvement_fraction=float(np.mean(difference < 0)) if metric in LOWER_IS_BETTER else
                            float(np.mean(difference > 0)) if metric in HIGHER_IS_BETTER else None)
                    if stage != 'E0':
                        low, high = np.quantile(difference[indices].mean(axis=1), [.025,.975])
                        row.update(CI95_low=float(low), CI95_high=float(high))
                    else:
                        row.update(CI95_low=None, CI95_high=None)
                out.append(row)
    return out


def _auxiliary(results):
    output = {}
    path = results/'compensation/COMPENSATION_MECHANISMS.json'
    if path.exists():
        data = _read(path)
        mechanisms = []
        for row in data['records']:
            distance = lambda name: row.get(name, {}).get('distances', {}).get('total_distance')
            full = row.get('full_validation', {})
            mechanisms.append(dict(geometry=row['geometry'], new_frequency=row['new_frequency'], status=row['status'],
                raw_proxy_distance=row.get('raw', {}).get('total_distance'),
                old_only_proxy_distance=distance('old_only'), new_only_proxy_distance=distance('new_only'),
                joint_shared_proxy_distance=distance('joint_shared'),
                independent_gain_proxy_distance=distance('independent_gain_negative_control'),
                calibration_old_proxy_distance=distance('calibration_plus_old'),
                calibration_joint_proxy_distance=distance('calibration_plus_joint'),
                full_joint_endpoint_distance=full.get('shared', {}).get('total_distance'),
                full_same_epsilon_shared_gain_refit=full.get('shared_gain_refit_at_same_epsilon', {}).get('total_distance'),
                full_same_epsilon_independent_gain_distance=full.get('independent_gain_at_same_epsilon', {}).get('total_distance'),
                full_same_epsilon_gain_conflict_squared=full.get('gain_conflict_squared'),
                full_required_gain_old=full.get('required_gain_old'), full_required_gain_new=full.get('required_gain_new'),
                full_required_gain_difference=full.get('required_gain_difference'),
                full_calibration_same_epsilon_distance=full.get('calibration_at_same_epsilon', {}).get('total_distance'),
                material_L1=row.get('joint_shared', {}).get('material_L1'),
                active_bounds=row.get('joint_shared', {}).get('active_bounds')))
        output['mechanisms'] = dict(source=str(path), observed_rows=len(mechanisms), expected_rows=12,
            full_validated_rows=data.get('full_validated_rows'), callback_attempts=data.get('callback_attempts_including_restored'),
            failures=data.get('failed_or_unresolved_callback_keys'), rows=mechanisms,
            scope='fixed-reference bounded best-found witnesses; proxy controls and same-epsilon full gain diagnostics remain distinct; not global minima or realized MAE')
    path = results/'E1_FORWARD_AUDIT.json'
    if path.exists():
        data = _read(path)
        rows = data['rows']
        metrics = ('relative', 'whitened', 'boundary_residual', 'n14_n18_white', 'recursion_white')
        summaries = []
        for frequency in sorted(set(row['k'] for row in rows)):
            subset = [row for row in rows if row['k']==frequency]
            current = dict(frequency=frequency, independent_materials=len(subset))
            for metric in metrics:
                available = [row[metric] for row in subset if row.get(metric) is not None]
                current[metric] = dict(available=len(available), maximum=max(available) if available else None,
                    p95=float(np.quantile(available,.95)) if available else None)
            summaries.append(current)
        output['forward_audit'] = dict(source=str(path), rows=len(rows), frequencies=summaries,
            uniform_bound=data.get('uniform_bound'), scope=data.get('scope'))
    path = results/'E1_INITIAL_REFERENCE_AUDIT.json'
    if path.exists():
        rows = _read(path)
        output['initial_reference_audit'] = dict(source=str(path), observations=len(rows),
            maximum_event_mass_difference=max((row['event_mass_max'] for row in rows),default=None),
            maximum_mean_difference=max((row['mean_max'] for row in rows),default=None),
            scope='initial-history order13/order17 continuous proxy integration only; future sequential particle ESS is not an error certificate')
    return output


def _csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(stages, job_id, root=ROOT):
    """Always charged as an independent ANALYSIS job; never calls evaluate."""
    root = Path(root)
    cost = Cost(root)
    with cost.job(job_id, 'ANALYSIS'):
        results = root/'results'
        parents, inputs = [], []
        for stage in stages:
            if stage == 'E0':
                parents.extend(_e0_parent_rows(results))
                inputs.append(dict(stage='E0', scope='exposed seven-case descriptive panel; no parent bootstrap'))
            else:
                raw, episodes, status = _load_main_stage(results,stage)
                parents.extend(_main_parent_rows(stage,episodes,status))
                inputs.append(dict(stage=stage, input_rows=len(raw), episodes=len(episodes), status=status))
        summary, groups, available = _summary(parents)
        paired = _paired(groups)
        auxiliary = _auxiliary(results)
        statistics = dict(schema='research-b-parent-statistics-v1', stages=list(stages), analysis_job_id=job_id,
            analysis_ledger_stage='ANALYSIS', analysis_ledger='COMPUTE_LEDGER.jsonl; appended by Cost.job on exit',
            new_solver_calls_by_this_analysis=0, GPU_training=0,
            bootstrap=dict(replicates=BOOTSTRAP_REPLICATES, seed=BOOTSTRAP_SEED,
                unit='independent parent; first average both noise repetitions',
                interval='2.5/97.5 percentile of paired parent-mean differences',
                delta='policy minus baseline; negative selected-MAE difference favors policy',
                multiple_comparisons='descriptive intervals, no multiplicity adjustment or automatic superiority conclusion',
                E0='no bootstrap: exposed development cases are not the 32-parent panel'),
            horizons='initial, step1/2/3 up to declared stage budget, final. Explicit stopped/none policies carry their last row; unfinished trajectories are rejected.',
            candidate_only_scope='mode1 is candidate profile material selection; reference-assisted acquisition policies are explicitly labeled, so those pipelines are not wholly candidate-only',
            selection_accuracy_definition='selected_structure_correct is agreement in the 12 posthoc central-radius/annular-width event bins; event_decision_accuracy may use a different maximum-mass event',
            ESS_scope='reference-particle ESS diagnostic in both output modes; no candidate posterior probabilities and no all-mode certificate',
            risk_scope='finite_reference_risk is reference-selection risk, even when displayed with mode1; candidate_reference_regret is the offline reference diagnostic of selected candidate',
            proposal_scope='mode1 pool is fixed625; mode2 refresh appends candidates from existing finite2048 support, not independent discovery of new material families',
            acquisition_cost_scope='saved compute_acquisition_seconds is assigned/shared recommendation computation, not total end-to-end cost; never sum shared values across policies',
            baseline_scopes={stage:{b:_baseline_scope(stage,b) for b in BASELINES} for stage in stages},
            same_parent_panel_across_E1_E2=True, E2_exposure='E1 parents reused, explicitly exposed followup; not a new blind cohort',
            summary_groups=len(summary), parent_rows=len(parents), paired_rows=len(paired),
            inputs=inputs, metric_availability=available, summaries=summary,
            paired_comparisons=paired, **auxiliary)
        _csv(results/'SUMMARY_METRICS.csv',summary)
        _csv(results/'PARENT_METRICS.csv',parents)
        _csv(results/'PAIRED_COMPARISONS.csv',paired)
        dump(results/'STATISTICAL_SUMMARY.json',statistics)
        print(json.dumps(dict(stages=list(stages),summary_groups=len(summary),parent_rows=len(parents),
            paired_rows=len(paired),bootstrap_unit='parent',new_Maxwell=0,GPU=0)),flush=True)
    return statistics


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stages', nargs='+', choices=['E0','E1','E2'], default=['E0','E1','E2'])
    parser.add_argument('--job', required=True)
    arguments = parser.parse_args()
    run(arguments.stages,arguments.job)
