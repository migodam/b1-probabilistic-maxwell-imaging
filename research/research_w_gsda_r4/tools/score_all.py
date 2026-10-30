#!/usr/bin/env python3
"""Score sealed endpoints and render actual nine-row digit panels. No PDE calls."""
import time
WRAPPER_WALL0 = time.perf_counter()
from pathlib import Path
import os
import resource
import shutil
import argparse
import csv
import json
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gsda.score import score_endpoint, plot_digit_rows, plot_digit_triptych, plot_case, plot_quality_cost, plot_heldout_cost, plot_all_method_comparison, plot_physics_evidence, plot_native_controls, plot_native_heldout, plot_tv_ablation
WRAPPER_META = dict(root=ROOT, job='score_all')


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def infer_case(root, epdir, seal):
    case = seal.get('case_id', seal.get('case'))
    if case is not None: return case
    candidates = [p for p in epdir.parts if (root / 'inputs' / 'private' / p / 'truth.npz').exists()]
    if len(candidates) != 1: raise ValueError('Cannot infer unique generic case ID from endpoint path/seal')
    return candidates[0]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--results', type=Path)
    p.add_argument('--no-plots', action='store_true')
    p.add_argument('--no-montage', action='store_true', help='Render selected per-case figures only')
    p.add_argument('--case', action='append', dest='cases', help='Score/render only these case IDs; reuse existing scores for other sealed cases')
    p.add_argument('--job', default='score_all')
    p.add_argument('--tv-ablation-only', action='store_true', help='Plot the actual sealed same-data c006 TV0 versus TV=.2 ablation')
    p.add_argument('--paired-b0-controls-only', action='store_true', help='Render complete paired controls with the three sealed B0 extensions')
    p.add_argument('--native-heldout-only', action='store_true', help='Refresh only native-control heldout figure and paired CSV fields')
    p.add_argument('--native-controls-only', action='store_true', help='Plot the actual complete paired 18 native control endpoints')
    p.add_argument('--physics-figures-only', action='store_true', help='Plot existing absorption/spectral/mesh evidence without physics')
    p.add_argument('--comparison-only', action='store_true', help='Render fixed truth|B0|B0+FullGN|B4 all-nine comparison from merged actual scores')
    p.add_argument('--heldout-cost-only', action='store_true', help='Refresh saved heldout columns and its new cost/quality figure only')
    p.add_argument('--quality-cost-only', action='store_true', help='Render final cost/quality plots from existing merged scores; requires actual B4 scores')
    p.add_argument('--endpoint', action='append', type=Path, dest='endpoints', help='Score/render this sealed endpoint directory (or endpoint.npz); retain older aggregate records')
    a = p.parse_args()
    root = a.root.resolve()
    result_root = (a.results or root / 'results').resolve()
    selected_epdirs = None if not a.endpoints else {(x.parent if x.name == 'endpoint.npz' else x).resolve() for x in a.endpoints}
    if selected_epdirs:
        for epdir in selected_epdirs:
            epdir.relative_to(root)
            if not (epdir / 'endpoint.npz').exists() or not (epdir / 'ONLINE_SEALED.json').exists():
                p.error(f'Explicit endpoint is not sealed and ready: {epdir}')
    WRAPPER_META.update(root=root, job=a.job, results=str(result_root), selected_cases=a.cases, selected_endpoints=None if selected_epdirs is None else [str(x) for x in sorted(selected_epdirs)])
    if a.physics_figures_only:
        result = plot_physics_evidence(root)
        WRAPPER_META.update(physics_figures_only=True)
        print(json.dumps(result))
        return
    if a.quality_cost_only or a.heldout_cost_only or a.comparison_only or a.native_controls_only or a.native_heldout_only or a.paired_b0_controls_only or a.tv_ablation_only:
        aggregate_path = root / 'evaluation' / 'SCORES_ALL.json'
        saved_records = json.loads(aggregate_path.read_text())['records']
        result = plot_tv_ablation(root, saved_records) if a.tv_ablation_only else (plot_native_controls(root, saved_records, include_b0=True) if a.paired_b0_controls_only else (plot_native_heldout(root, saved_records) if a.native_heldout_only else (plot_native_controls(root, saved_records) if a.native_controls_only else (plot_all_method_comparison(root, saved_records) if a.comparison_only else (plot_heldout_cost(root, saved_records) if a.heldout_cost_only else plot_quality_cost(root, saved_records))))))
        WRAPPER_META.update(quality_cost_only=a.quality_cost_only, heldout_cost_only=a.heldout_cost_only, comparison_only=a.comparison_only, native_controls_only=a.native_controls_only, native_heldout_only=a.native_heldout_only, paired_b0_controls_only=a.paired_b0_controls_only, tv_ablation_only=a.tv_ablation_only, actual_rows=result.get('actual_rows', result.get('actual_heldout_rows')))
        print(json.dumps(result))
        return
    records, failures, unsealed, blocked, newly_scored = [], [], [], [], []
    for epdir in sorted({x.parent for x in result_root.rglob('endpoint.npz')}):
        sealfile = epdir / 'ONLINE_SEALED.json'
        if not sealfile.exists():
            unsealed.append(str(epdir.relative_to(root)))
            continue
        try:
            seal = json.loads(sealfile.read_text())
            case = infer_case(root, epdir, seal)
            if (a.cases and case not in a.cases) or (selected_epdirs is not None and epdir not in selected_epdirs):
                existing = epdir / 'SCORES.json'
                if existing.exists(): records.append(json.loads(existing.read_text()))
                continue
            row = score_endpoint(root, epdir, case)
            records.append(row)
            newly_scored.append(str(epdir.relative_to(root)))
            dump(epdir / 'SCORES.json', row)
            if not a.no_plots: plot_case(root, row)
        except Exception as exc:
            if 'METHOD_FREEZE required' in str(exc):
                blocked.append(str(epdir.relative_to(root)))
                continue
            failures.append(dict(endpoint_dir=str(epdir.relative_to(root)), error=repr(exc)))
    # Persist all previous scored records when the caller processes one new case/campaign.
    aggregate = root / 'evaluation' / 'SCORES_ALL.json'
    older = json.loads(aggregate.read_text()).get('records', []) if aggregate.exists() else []
    merged = {r['endpoint_dir']: r for r in older}
    merged.update({r['endpoint_dir']: r for r in records})
    records = [merged[k] for k in sorted(merged)]
    figures = []
    if not a.no_plots and not a.no_montage:
        groups = sorted({(r['campaign'], r['method']) for r in records if r['cohort'] == 'digit'})
        for campaign, method in groups:
            selected = [r for r in records if r['cohort'] == 'digit' and r['method'] == method and r['campaign'] == campaign]
            ids = [r['case_id'] for r in selected]
            if len(set(ids)) != len(ids):
                failures.append(dict(method=method, campaign=campaign, error='Multiple sealed endpoints per case/method; explicit measurement-only endpoint selection required'))
                continue
            safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in campaign + '_' + method)
            figures.append(plot_digit_rows(root, selected, method, root / 'figures' / f'{safe}_DIGITS_TRUTH_RECON_ERROR'))
            for quantity in ('contrast', 'loss'):
                figures.append(plot_digit_triptych(root, selected, method, campaign,
                    root / 'figures' / f'{safe}_DIGITS_{quantity.upper()}_9X3', quantity=quantity))
            if campaign in ('B0_FULL_FWI', 'B4_GSDA', 'FULL_GN_TV_SUPPLEMENT') and len(selected) == 9:
                alias_tag = {'B0_FULL_FWI':'B0', 'B4_GSDA':'B4', 'FULL_GN_TV_SUPPLEMENT':'FULL_GN_B0_WARM'}[campaign]
                for source_suffix, alias in (('TRUTH_RECON_ERROR', 'TRUTH_RECON_ERROR'), ('CONTRAST_9X3', 'CONTRAST'), ('LOSS_9X3', 'LOSS')):
                    for extension in ('png', 'pdf'):
                        source = root / 'figures' / f'{safe}_DIGITS_{source_suffix}.{extension}'
                        shutil.copyfile(source, root / f'DIGITS_1_TO_9_{alias_tag}_{alias}.{extension}')
    out = root / 'evaluation'
    dump(out / 'SCORES_ALL.json', dict(status='POST_SEAL_ARRAY_SCORED' if records else 'NOT_RUN',
         records=records, failures=failures, unsealed_endpoints=unsealed, digit_endpoints_blocked_until_method_freeze=blocked,
         figures=figures, support_threshold=1., scoring_physics_calls=0,
         nine_digit_reconstruction_claim=any(len({r['case_id'] for r in records if r.get('cohort') == 'digit' and r.get('campaign') == c}) == 9 for c in {r.get('campaign') for r in records})))
    flat = []
    for r in records:
        row = {k: r.get(k) for k in ('case_id', 'cohort', 'shape', 'digit', 'method', 'campaign', 'status', 'online_status', 'endpoint_dir', 'objective', 'physical_data_objective', 'stop', 'converged', 'selected_scaled_projected_gradient_inf', 'full_objective_comparable', 'selection', 'total_objective_calls', 'scope', 'cost_basis', 'not_cold_start', 'independent_blind_claim')}
        if row['physical_data_objective'] is None and r['metrics'].get('observation_white_norm') is not None:
            row['physical_data_objective'] = .5 * r['metrics']['observation_white_norm'] ** 2
        row.update({k: v for k, v in r['metrics'].items() if not isinstance(v, (dict, list))})
        row.update({f'true_{k}': v for k, v in r['metrics']['true_topology'].items()})
        row.update({f'predicted_{k}': v for k, v in r['metrics']['predicted_topology'].items()})
        row.update(r.get('current_selection', {}))
        row.update({k:(json.dumps(v, separators=(',', ':')) if isinstance(v, (dict,list)) else v) for k,v in r.get('mechanism', {}).items()})
        row['starts_and_stops'] = json.dumps(r.get('starts'), separators=(',', ':'))
        row['error'] = json.dumps(r.get('error'), separators=(',', ':'))
        guard = r.get('guard') or {}
        row['guard_denied_reads'] = len(guard['denied']) if 'denied' in guard else None
        row['guard_online_opens'] = guard.get('online_opens')
        row['guard_scope'] = guard.get('scope')
        row.update({f'cost_{k}': r['cost'].get(k) for k in ('full_fe','rhs','wrapper_process_cpu_seconds','wall_seconds','peak_rss_bytes')})
        row['whole_child_CPU'] = (r.get('end_to_end_cost') or {}).get('child_total_cpu_seconds')
        row['whole_child_wall'] = (r.get('end_to_end_cost') or {}).get('external_wall_seconds')
        row['whole_child_not_additive_to_Book'] = (r.get('end_to_end_cost') or {}).get('do_not_double_add')
        row.update({f'additional_cost_{k}': (r.get('additional_cost') or {}).get(k) for k in ('full_fe','rhs','wrapper_process_cpu_seconds','wall_seconds')})
        row.update({f'heldout_cost_{k}': r['heldout_cost'].get(k) for k in ('full_fe','rhs','wrapper_process_cpu_seconds','cpu_seconds','wall_seconds')})
        flat.append(row)
    if flat:
        fields = sorted(set().union(*(r.keys() for r in flat)))
        with (out / 'SCORES_ALL.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(flat)
    actual_scientific = [r for r in flat if r.get('cohort') != 'development' and r.get('status') == 'POST_SEAL_ARRAY_SCORED']
    if actual_scientific:
        experiment_fields = sorted(set().union(*(r.keys() for r in actual_scientific)))
        with (root / 'EXPERIMENT_RESULTS.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=experiment_fields); writer.writeheader(); writer.writerows(actual_scientific)
    control_rows = []
    for family in ('B1', 'B1_RESET', 'B2', 'B3', 'B4', 'B5'):
        campaign = 'CONTROL_' + family
        for case_id in ('c010', 'c011', 'c012'):
            matched = [r for r in flat if r.get('campaign') == campaign and r.get('case_id') == case_id]
            if matched:
                control_rows.extend(matched)
                continue
            folder = root / 'results' / campaign / case_id
            status = 'SEALED_NOT_SCORED' if (folder / 'ONLINE_SEALED.json').exists() else ('FAILED_NOT_SCORED' if (folder / 'FAILURE.json').exists() else ('UNSEALED_NOT_SCORED' if folder.exists() else 'NOT_RUN'))
            control_rows.append(dict(campaign=campaign, method=family, case_id=case_id, cohort='mechanism_control', status=status))
    control_fields = sorted(set().union(*(r.keys() for r in control_rows)))
    with (out / 'NATIVE_CONTROL_METRICS.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=control_fields); writer.writeheader(); writer.writerows(control_rows)
    primary_campaigns = {'B0_FULL_FWI','FULL_GN_TV_SUPPLEMENT','B4_GSDA'}
    exposed_digit_ablations = [r for r in flat if r.get('cohort') == 'digit' and r.get('campaign') not in primary_campaigns]
    if exposed_digit_ablations:
        fields = sorted(set().union(*(r.keys() for r in exposed_digit_ablations)))
        with (out / 'DIGIT_EXPOSED_ABLATION_METRICS.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(exposed_digit_ablations)
    digit_flat = [r for r in flat if r.get('cohort') == 'digit' and r.get('campaign') in primary_campaigns]
    frozen = any((root / name).is_file() for name in ('METHOD_FREEZE', 'METHOD_FREEZE.json', 'METHOD_FREEZE.md'))
    if frozen:
        # Private mapping is used only after the orchestrator freezes the method.
        case_map_path = root / 'inputs/private/CASE_MAP.json'
        case_map = json.loads(case_map_path.read_text())['cases'] if case_map_path.exists() else []
        digit_ids = {r['case_id']: r['digit'] for r in case_map if r['cohort'] == 'digit'}
        seen = {r['endpoint_dir'] for r in digit_flat}
        for batchfile in result_root.rglob('BATCH_STATUS.json'):
            batch = json.loads(batchfile.read_text()); batchdir = batchfile.parent
            attempts = {r['case']: r for r in batch.get('attempts', [])}
            for case_id in batch.get('cases', []):
                if case_id not in digit_ids: continue
                epdir = batchdir / case_id; relative = str(epdir.relative_to(root))
                if relative in seen: continue
                att = attempts.get(case_id, {})
                receipt = json.loads((epdir / 'RECEIPT.json').read_text()) if (epdir / 'RECEIPT.json').exists() else {}
                failure = json.loads((epdir / 'FAILURE.json').read_text()) if (epdir / 'FAILURE.json').exists() else att.get('reason')
                digit_flat.append(dict(case_id=case_id, cohort='digit', digit=digit_ids[case_id],
                    method=batchdir.name, campaign=batchdir.name, endpoint_dir=relative,
                    status=('UNSEALED_NOT_SCORED' if epdir.exists() and not (epdir / 'ONLINE_SEALED.json').exists() else 'NOT_RUN') if not att else 'ATTEMPTED_NOT_SCORED', online_status=att.get('online_status', att.get('status')),
                    error=json.dumps(failure), returncode=att.get('returncode'),
                    cost_full_fe=receipt.get('full_fe'), cost_rhs=receipt.get('rhs'),
                    cost_wrapper_process_cpu_seconds=receipt.get('wrapper_process_cpu_seconds'),
                    cost_wall_seconds=receipt.get('wall_seconds')))
                seen.add(relative)
    if digit_flat:
        fields = sorted(set().union(*(r.keys() for r in digit_flat)))
        with (out / 'DIGITS_1_TO_9_METRICS.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(digit_flat)
        shutil.copyfile(out / 'DIGITS_1_TO_9_METRICS.csv', root / 'DIGITS_1_TO_9_METRICS.csv')
        for campaign, alias_tag in (('B0_FULL_FWI', 'B0'), ('B4_GSDA', 'B4'), ('FULL_GN_TV_SUPPLEMENT', 'FULL_GN_B0_WARM')):
            actual_method = [r for r in digit_flat if r.get('campaign') == campaign and r.get('status') == 'POST_SEAL_ARRAY_SCORED']
            if len({r['case_id'] for r in actual_method}) == 9:
                with (root / f'DIGITS_1_TO_9_{alias_tag}_METRICS.csv').open('w', newline='') as f:
                    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(actual_method)
    WRAPPER_META.update(aggregate_records=len(records), newly_scored_endpoints=newly_scored, newly_scored_count=len(newly_scored), failures=len(failures), unsealed=len(unsealed), rendered_montages=len(figures))
    print(json.dumps(dict(aggregate_records=len(records), newly_scored_endpoints=newly_scored, failures=failures, unsealed=len(unsealed), digit_scoring_blocked=len(blocked), figures=figures,
                          scoring_physics_calls=0)))
    if failures: raise SystemExit(1)

if __name__ == '__main__':
    status, error = 'FAILED', None
    try:
        main()
        status = 'SCORING_COMPLETED'
    except BaseException as exc:
        error = repr(exc)
        raise
    finally:
        receipt_root = Path(WRAPPER_META.pop('root'))
        job = ''.join(c if c.isalnum() or c in '-_' else '_' for c in WRAPPER_META.pop('job'))
        receipt = dict(job=job, pid=os.getpid(), status=status,
            wrapper_process_cpu_seconds=time.process_time(),
            wall_seconds=time.perf_counter() - WRAPPER_WALL0,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
            full_fe=0, rhs=0, error=error, scope='Post-seal array scoring and plotting; no physics', **WRAPPER_META)
        dump(receipt_root / 'receipts' / f'{job}_{os.getpid()}.json', receipt)

