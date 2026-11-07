"""Read recorded R5 receipts once; no solves, inference, training, or GT access.

Exact phase/model snapshots only. Nested physics timing is a breakdown of run
wall, not an additional task; original GN deployment is reference-only.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import time

STUDY = Path(__file__).resolve().parents[1]
ARMS = ('pixel_closed_coupled', 'haar_closed_coupled', 'db2_closed_coupled',
        'haar_closed_no_fullfield', 'haar_free_coupled')
EXTRA = 'haar_closed_noiseball_coupled'
COUNT_KEYS = {
    'LU': ('LU',),
    'source_complex_RHS': ('source_complex_RHS', 'source_RHS'),
    'loss_VJP_adjoint_complex_RHS': ('adjoint_complex_RHS',),
    'receiver_transpose_complex_RHS': ('receiver_transpose_RHS',),
    'actual_complex_RHS': ('actual_complex_RHS',),
    'G_matmul_calls': ('G_matmul_calls',),
    'G_adjoint_matmul_calls': ('G_adjoint_matmul_calls',),
    'S_matmul_calls': ('S_matmul_calls',),
    'S_adjoint_matmul_calls': ('S_adjoint_matmul_calls',),
}


def count(counters, alternatives):
    return next((counters[k] for k in alternatives if k in counters), None)


def counts(counters):
    return {key: count(counters, aliases) for key, aliases in COUNT_KEYS.items()}


def sum_known(values):
    values = list(values)
    return sum(values) if values and all(x is not None for x in values) else None


class Collector:
    def __init__(self, study):
        self.study = study
        self.rows, self.missing, self.sources, self.cautions = [], [], {}, []

    def read(self, path, required=True):
        path = Path(path)
        if not path.exists():
            if required:
                self.missing.append(str(path.relative_to(self.study)))
            return None
        stat = path.stat()
        self.sources[str(path.relative_to(self.study))] = dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
        return json.loads(path.read_text())

    def row(self, phase, category, arm, source, *, counters=None, fixed=None,
            task_wall=None, elapsed=None, physics=None, reference=False, **extra):
        values = counts(counters or {})
        if fixed:
            values.update({k: v for k, v in counts(fixed).items() if v is not None})
        self.rows.append(dict(scope='INDEPENDENT_DEPLOYMENT_REFERENCE' if reference else 'ACTUAL_NEW_EXECUTION',
                              phase=phase, category=category, arm=arm, source=str(source),
                              task_wall_seconds=task_wall, recorded_run_elapsed_seconds=elapsed,
                              physics_wall_seconds_nested=physics, **{**values, **extra}))

    def bridge_run(self, phase, category, arm, path):
        d = self.read(path)
        if d is None:
            return
        bridge = d.get('physical_validation', d if 'counts' in d else {})
        if not bridge.get('counts'):
            self.cautions.append('No recognized physics count receipt: '+str(path.relative_to(self.study)))
        self.row(phase, category, arm, path.relative_to(self.study),
                 counters=bridge.get('counts', {}), fixed=d.get('fixed_operator_calls'),
                 task_wall=d.get('wall_seconds'), elapsed=d.get('wall_seconds'),
                 physics=bridge.get('total_wall_seconds'),
                 cases=d.get('cases'), NN_loop_seconds=d.get('NN_loop_seconds'),
                 setup_seconds_nested=bridge.get('setup_seconds'),
                 status=d.get('status', bridge.get('runtime_status')),
                 time_scope='Recorded complete run wall; physics and NN timings are nested components, not additive')

    def training(self):
        for phase, folder, arms in [('strict', 'models', ARMS), ('noiseaware', 'models_noiseaware', (*ARMS, EXTRA))]:
            for arm in arms:
                directory = self.study/folder/arm
                summary = self.read(directory/'SUMMARY.json')
                bridge = self.read(directory/'physics_cost.json')
                if summary is None or bridge is None:
                    continue
                if summary.get('physics_counts') != bridge.get('counts'):
                    self.cautions.append(f'{folder}/{arm}: SUMMARY physics counters differ; standalone physics_cost used, not added')
                self.row(phase, 'training_including_VAL_selection', arm,
                         str((directory/'SUMMARY.json').relative_to(self.study))+' + physics_cost.json',
                         counters=bridge.get('counts', {}), fixed=summary.get('fixed_operator_calls'),
                         task_wall=summary.get('wall_seconds'), elapsed=summary.get('wall_seconds'),
                         physics=bridge.get('total_wall_seconds'), setup_seconds_nested=bridge.get('setup_seconds'),
                         epochs=summary.get('epochs'), best_epoch=summary.get('best_epoch'), status=summary.get('status'),
                         no_field_scope='No-full-field arm skips TRAIN forward loss, but actual VAL checkpoint-selection forwards are included',
                         NN_backward_FLOPs='UNRECORDED_UNKNOWN',
                         time_scope='One recorded per-model training run; exact bridge timing is nested, no phase-global elapsed inferred')
        path = self.study/'audit/ACTUAL_TRAIN_WIRING.json'
        d = self.read(path)
        if d:
            bridge = d.get('physics', {})
            self.row('strict', 'separate_actual_train_wiring', d.get('arm'), path.relative_to(self.study),
                     counters=bridge.get('counts', {}), fixed=d.get('fixed_operator_calls'),
                     task_wall=d.get('metrics', {}).get('seconds'),
                     elapsed=d.get('metrics', {}).get('seconds'), physics=bridge.get('total_wall_seconds'),
                     time_scope='Separate discarded-weight wiring invocation; not included in later model SUMMARY',
                     NN_backward_FLOPs='UNRECORDED_UNKNOWN')

    def preparation(self):
        path = self.study/'audit/NOISEAWARE_PREPARATION_RECEIPT.json'
        d = self.read(path)
        if d is None:
            return
        generated = [r for r in d['rows'] if r['disposition'] == 'GENERATED']
        copied = [r for r in d['rows'] if r['disposition'] == 'COPIED_FROZEN_CACHE']
        # Sum row-level original counters to retain source and receiver-J RHS.
        combined = {k: sum_known(count(r.get('actual_new_physics_counts') or {}, aliases) for r in generated)
                    for k, aliases in COUNT_KEYS.items()}
        self.row('noiseaware', 'new_common_TRAIN_anchor_preparation', 'shared_all_arms',
                 path.relative_to(self.study), counters=d['actual_new_physics_counts'],
                 task_wall=d.get('actual_train_task_wall_seconds_sum'),
                 elapsed=d.get('actual_elapsed_wall_seconds'),
                 cases=d.get('train_generated'), status=d['status'],
                 time_scope='Taskwall sums parallel TRAIN worker tasks; elapsed is entire preparation including VAL copies, not sumtaskwall',
                 **{k: v for k, v in combined.items() if v is not None and k not in ('LU', 'actual_complex_RHS')})
        self.row('noiseaware', 'cached_VAL_copy_operation', 'shared_all_arms', path.relative_to(self.study),
                 counters={'LU': 0, 'actual_complex_RHS': 0},
                 task_wall=d.get('actual_copy_wall_seconds_sum'), cases=len(copied),
                 time_scope='New copy wall only; copied original physics counters excluded from new execution')
        keys = set(k for r in copied for k in r.get('original_deployment_counts', {}))
        ref = {k: sum(r.get('original_deployment_counts', {}).get(k, 0) for r in copied) for k in keys}
        self.row('noiseaware', 'cached_VAL_original_GN_deployment_reference', 'shared_inputs_156',
                 path.relative_to(self.study), reference=True, counters=ref,
                 task_wall=sum_known(r.get('original_deployment_anchor_wall_seconds') for r in copied), cases=len(copied),
                 time_scope='Historical necessary deployment reference, not new R5 generation or added new execution')

    def predictions_and_depth(self):
        for phase, relative, arms in [('strict', 'results/VAL_PREDICTIONS', ARMS),
                                      ('noiseaware', 'results/noiseaware/VAL_PREDICTIONS', (*ARMS, EXTRA))]:
            for arm in arms:
                self.bridge_run(phase, 'public_prediction_and_final_exact_forward', arm, self.study/relative/f'{arm}.json')
        # Explicit snapshots only: never recursively gather archived runs.
        for phase, relative, arms in [('strict', 'results/depth', ARMS),
                                      ('noiseaware', 'results/noiseaware/depth', ARMS),
                                      ('noiseaware_extra', 'results/noiseaware/depth_noiseball', (EXTRA,))]:
            directory = self.study/relative
            for arm in arms:
                self.bridge_run(phase, 'fixed_depth_public_inference_and_final_forward', arm, directory/f'{arm}.json')
            self.bridge_run(phase, 'separate_depth_stage_physical_trace', 'shared_trace', directory/'PHYSICAL_TRACE_COST.json')

    def corrections(self):
        for phase, resultdir, auditdir, tag in [
                ('strict', 'results', 'audit', ''),
                ('noiseaware', 'results/noiseaware', 'audit/noiseaware', ''),
                ('noiseaware_extra', 'results/noiseaware', 'audit/noiseaware', 'noiseball')]:
            suffix = '.'+tag if tag else ''
            plan_path = self.study/auditdir/f'CORRECTOR_GROUP_SELECTION{suffix}.json'
            seal_path = self.study/auditdir/f'CORRECTOR_ONLINE_SEALED{suffix}.json'
            plan, seal = self.read(plan_path), self.read(seal_path)
            if plan is None or seal is None:
                continue
            for arm in plan['arms']:
                receipt_paths = [self.study/resultdir/f'CORRECTED_ONLINE{suffix}'/arm/f"{r['obs_id']}.json" for r in plan['records']]
                receipts = [r for r in (self.read(p) for p in receipt_paths) if r is not None]
                def summed(key, reference=False):
                    source = 'necessary_anchor_deployment' if reference else 'costcounts'
                    return sum_known(r.get(source, {}).get(key) for r in receipts)
                counter_keys = set(k for r in receipts for k in r.get('costcounts', {}))
                c = {k: summed(k) for k in counter_keys}
                self.row(phase, 'physical_corrector', arm,
                         str((self.study/resultdir/f'CORRECTED_ONLINE{suffix}'/arm).relative_to(self.study))+'/<frozen obs_id>.json',
                         counters=c,
                         task_wall=sum_known(r.get('worker_total_seconds_including_io_and_setup') for r in receipts),
                         cases=len(receipts),
                         status_counts=dict(Counter(r['status'] for r in receipts)),
                         incremental_budgeted_wall_seconds=sum_known(r.get('wall_seconds') for r in receipts),
                         NN_preparation_allocated_seconds_NOT_ADDED=sum_known(r.get('preparation_cost', {}).get('wall_seconds') for r in receipts),
                         accounting_complete=sum(r.get('accounting_status') == 'COMPLETE' for r in receipts),
                         time_scope='Workerwall includes actual IO/setup/correction only. Budgetwall includes prior allocated NN time and is diagnostic, NOT added again',
                         exact_snapshot_plan=str(plan_path.relative_to(self.study)), exact_snapshot_seal=str(seal_path.relative_to(self.study)))
                references = [r.get('necessary_anchor_deployment') for r in receipts]
                if all(x is not None for x in references):
                    ref_c = {k: sum_known(x.get(k) for x in references) for k in ('LU', 'actual_complex_RHS')}
                    self.row(phase, 'per_method_necessary_GN_deployment_reference', arm,
                             str((self.study/resultdir/f'CORRECTED_ONLINE{suffix}'/arm).relative_to(self.study))+'/<frozen obs_id>.json',
                             reference=True, counters=ref_c, task_wall=sum_known(x.get('wall_seconds') for x in references),
                             cases=len(references), time_scope='Independent per-method deployment invoice; reused common anchors repeated across method views, NEVER summed as new execution')

    def emit(self):
        output = self.study/'results'
        keys = list(dict.fromkeys(k for r in self.rows for k in r))
        with (output/'COST_LEDGER.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=keys)
            writer.writeheader()
            for row in self.rows:
                writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})
        actual = [r for r in self.rows if r['scope'] == 'ACTUAL_NEW_EXECUTION']
        by_category = defaultdict(list)
        for row in actual:
            by_category[row['category']].append(row)
        def totals(rows):
            return {key: dict(recorded_sum=sum(r[key] for r in rows if r.get(key) is not None),
                              rows_with_recorded_value=sum(r.get(key) is not None for r in rows),
                              rows_without_recorded_value=sum(r.get(key) is None for r in rows))
                    for key in (*COUNT_KEYS, 'task_wall_seconds')}
        notes = dict(schema='r5-actual-cost-accounting-v1', collection_unix=time.time(),
                     status='COMPLETE_RECORD_COLLECTION' if not self.missing else 'INCOMPLETE_MISSING_RECEIPTS',
                     ledger_rows=len(self.rows), actual_new_execution_rows=len(actual),
                     independent_deployment_reference_rows=len(self.rows)-len(actual),
                     totals_actual_recorded=totals(actual),
                     category_totals={k: totals(v) for k, v in by_category.items()},
                     missing_expected_receipts=self.missing, cautions=self.cautions,
                     exact_snapshot_roots=['models/<primary5>', 'models_noiseaware/<primary5+noiseball>',
                         'audit/ACTUAL_TRAIN_WIRING.json', 'audit/NOISEAWARE_PREPARATION_RECEIPT.json',
                         'results/VAL_PREDICTIONS/<primary5>.json', 'results/noiseaware/VAL_PREDICTIONS/<primary5+noiseball>.json',
                         'results/depth/<primary5>.json + PHYSICAL_TRACE_COST.json',
                         'results/noiseaware/depth/<primary5>.json + PHYSICAL_TRACE_COST.json',
                         'results/noiseaware/depth_noiseball/<noiseball>.json + PHYSICAL_TRACE_COST.json',
                         'Correction receipt names fixed by corresponding GROUP_SELECTION/ONLINE_SEALED files'],
                     source_receipt_count=len(self.sources),
                     full_bridge_counts='physics_cost or physical_validation counts counted once; matching SUMMARY counts are duplicate views, not additional work',
                     source_vs_adjoint='Source RHS, retained-LU loss-VJP adjoint RHS and receiver-transpose Jacobian RHS separately labelled. Missing fields stay blank/unknown, not invented zero.',
                     fixed_operators='G/S matmul counts are NOT Maxwell complex RHS. NN backward FLOPs are not recorded and remain unknown.',
                     global_elapsed='Recorded per-run elapsed retained in ledger. Sum of worker/task wall is NOT end-to-end elapsed; concurrent execution prevents a total elapsed/speed claim. Physics wall and NN loop wall are nested breakdowns, not added to run wall.',
                     correction_time='Predicted NN inference wall charged once in public prediction run. Corrector budget includes allocated NN prep diagnostically; new actual task wall uses worker_total_seconds, avoiding duplicated preparation.',
                     original_anchors='Historical strict anchors/reference creation are outside new R5 execution; required deployment is not free. Noiseaware TRAIN generation charged once; cached VAL copies have zero new physics. All deployment-reference rows excluded new totals and are not mutually additive.',
                     no_field='No-field control TRAIN skips exact-F objective; actual VAL selection F calls remain in model bridge totals.',
                     unavailable_counts_policy='Recorded sums with missing row values are labelled partial, not silently completed. No theoretical counts or FLOPs substituted.',
                     no_new_numerical_work=True, truth_or_TEST_arrays_read=False)
        t = notes['totals_actual_recorded']
        notes['RHS_breakdown_unclassified_recorded_total'] = t['actual_complex_RHS']['recorded_sum']-sum(t[k]['recorded_sum'] for k in ('source_complex_RHS', 'loss_VJP_adjoint_complex_RHS', 'receiver_transpose_complex_RHS'))
        notes['RHS_breakdown_scope'] = 'Unclassified is arithmetic remainder of authoritative total after recorded named components; no inferred assignment to source or adjoint categories.'
        notes['standalone_physical_trace_wall_seconds'] = sum(r.get('physics_wall_seconds_nested') or 0. for r in actual if r['category'] == 'separate_depth_stage_physical_trace')
        (output/'COST_NOTES.json').write_text(json.dumps(notes, indent=2, allow_nan=False)+'\n')
        print(json.dumps(dict(status=notes['status'], rows=len(self.rows), source_receipts=len(self.sources),
                              missing=len(self.missing), actual_recorded_LU=notes['totals_actual_recorded']['LU'],
                              actual_recorded_complex_RHS=notes['totals_actual_recorded']['actual_complex_RHS'])), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=STUDY)
    args = parser.parse_args()
    c = Collector(args.study.resolve())
    c.training()
    c.preparation()
    c.predictions_and_depth()
    c.corrections()
    c.emit()


if __name__ == '__main__':
    main()
