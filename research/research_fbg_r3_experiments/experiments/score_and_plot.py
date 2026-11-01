#!/usr/bin/env python3
"""Post-closure material scoring. Never import this module in an online runner.

All private file access occurs in main(), after the explicit closure guard.
No optimization or new Maxwell solve is performed. Physical residuals come from
sealed predictions / paid endpoint ledgers. Candidate GT scores are oracle-only.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import gzip
import json
import math
import numpy as np

STUDY = Path(__file__).resolve().parents[1]
DEFAULT_PRIVATE = STUDY.parents[3] / 'fbg_r3e_private_20261101'


def clean(value):
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value): return None
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    return value


def otsu_mask(material):
    """256-bin histogram Otsu; uniform material has empty foreground."""
    material = np.asarray(material, float).ravel()
    if not np.all(np.isfinite(material)): raise ValueError('Nonfinite material')
    if np.ptp(material) <= 1e-12: return np.zeros(material.size, bool), float(material[0])
    hist, edges = np.histogram(material, bins=256)
    centers = (edges[1:] + edges[:-1]) / 2
    weight = np.cumsum(hist)
    moment = np.cumsum(hist * centers)
    denom = weight * (weight[-1] - weight)
    between = np.full(hist.size, -np.inf)
    good = denom > 0
    between[good] = (moment[-1] * weight[good] - moment[good] * weight[-1])**2 / denom[good]
    threshold = float(centers[int(np.argmax(between))])
    return material > threshold, threshold


def material_metrics(material, truth, amplitude, support):
    material, truth = np.asarray(material, float).ravel(), np.asarray(truth, float).ravel()
    if material.shape != truth.shape: raise ValueError('Prediction / truth dimensions differ')
    fg, threshold = otsu_mask(material)
    support = np.asarray(support, bool).ravel()
    union = np.count_nonzero(fg | support)
    iou = np.count_nonzero(fg & support) / union if union else 1.
    nrmse = float(np.linalg.norm(material - truth) / np.sqrt(truth.size) / amplitude)
    return dict(amplitude_normalized_nrmse=nrmse, otsu_iou=float(iou),
                otsu_threshold=threshold, material_nrmse_pass=bool(nrmse <= .15),
                image_pass=bool(nrmse <= .15 and iou >= .7))


def npz_dict(path):
    if not path.exists(): return {}
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def write_csv(path, rows):
    if not rows: return
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as file:
        writer = csv.DictWriter(file, columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(clean(v), sort_keys=True) if isinstance(v, (dict, list)) else clean(v)
                             for k, v in row.items()})


def group_id(row):
    parts = [str(row['method'])]
    for key in ('k', 'budget_label', 'budget', 'wall_budget', 'noise_label', 'noise'):
        if key in row: parts.append(f'{key}={row[key]}')
    return ' / '.join(parts)


def endpoint_metadata(prefix):
    """Read already-paid candidate physical scores; never recompute forward maps."""
    path = Path(str(prefix) + '_ledger.jsonl.gz')
    out = {}
    if path.exists():
        with gzip.open(path, 'rt') as file:
            for line in file:
                run = json.loads(line)
                for i, endpoint in enumerate(run.get('corrected', [])):
                    out[run['run_id'] + f'_candidate{i}'] = {
                        'relative_residual': endpoint.get('relative_residual'),
                        'provenance': endpoint.get('provenance', {})}
    return out


def cluster_summary(rows):
    by_cluster = defaultdict(list)
    for row in rows: by_cluster[row['object_id']].append(row)
    clusters = []
    for oid, selected in sorted(by_cluster.items()):
        good = [r for r in selected if r.get('amplitude_normalized_nrmse') is not None]
        clusters.append(dict(object_id=oid, run_count=len(selected),
                             data_pass_fraction=float(np.mean([r['data_pass'] for r in selected])),
                             joint_pass_fraction=float(np.mean([r['joint_pass'] for r in selected])),
                             median_nrmse=float(np.median([r['amplitude_normalized_nrmse'] for r in good])) if good else None,
                             median_iou=float(np.median([r['otsu_iou'] for r in good])) if good else None))
    # Intervals describe finite-object resampling only, not population inference.
    values = np.array([r['joint_pass_fraction'] for r in clusters])
    rng = np.random.default_rng(20261101)
    interval = np.quantile(rng.choice(values, (10000, len(values)), replace=True).mean(1), [.025, .975]) if len(values) else [None, None]
    return dict(object_count=len(clusters), clusters=clusters,
                equal_object_mean_joint_pass=float(np.mean(values)) if len(values) else None,
                object_bootstrap_95_interval=clean(interval),
                inference_scope='Descriptive finite-object resampling; shared shapes/amplitudes/grids are not a population sample.')


def geometry_cluster_summary(rows):
    families = defaultdict(list)
    for row in rows: families[row['shape']].append(row)
    records = [dict(shape=shape, material_contracts=len({r['object_id'] for r in selected}),
                    run_count=len(selected), joint_pass_fraction=float(np.mean([r['joint_pass'] for r in selected])),
                    material_data_pass_fraction=float(np.mean([r.get('material_data_pass', False) for r in selected])))
               for shape, selected in sorted(families.items())]
    values = np.array([r['joint_pass_fraction'] for r in records])
    rng = np.random.default_rng(20261101)
    interval = np.quantile(rng.choice(values, (10000,len(values)),replace=True).mean(1), [.025,.975]) if len(values) else [None,None]
    return dict(fixed_geometry_families=len(records), material_contracts=len({r['object_id'] for r in rows}),
                families=records, equal_geometry_mean_joint_pass=float(np.mean(values)) if len(values) else None,
                geometry_bootstrap_95_interval=clean(interval),
                inference_scope='Primary descriptive resampling across fixed shape families; only three geometries in E1, very low confidence, not independent morphology draws.')


def summarize(rows):
    groups = defaultdict(list)
    for row in rows: groups[row['group']].append(row)
    summary = []
    for group, selected in sorted(groups.items()):
        good = [r for r in selected if r.get('amplitude_normalized_nrmse') is not None]
        oracle = [r for r in selected if r.get('oracle_best_image_nrmse') is not None]
        summary.append(dict(group=group, runs=len(selected), failures=sum(r['score_status'] != 'SCORED' for r in selected),
                            stopping_counts=dict(Counter(str(r.get('stop','NOT_RECORDED')) for r in selected)),
                            physical_budget_cap_stops=sum(r.get('stop') in ('FULL_LU_CAP','FULL_COMPLEX_RHS_CAP') for r in selected),
                            maximum_scored_unit_LU=max((r.get('model_cost',{}).get('LU',0) for r in selected), default=0),
                            maximum_scored_unit_RHS=max((r.get('model_cost',{}).get('actual_complex_RHS',0) for r in selected), default=0),
                            physical_caps_recorded_runs=sum(bool(r.get('physical_caps')) for r in selected),
                            physical_cap_overshoots=sum(bool(r.get('physical_caps')) and (r.get('model_cost',{}).get('LU',0)>(r['physical_caps'].get('max_lu') or np.inf) or r.get('model_cost',{}).get('actual_complex_RHS',0)>(r['physical_caps'].get('max_rhs') or np.inf)) for r in selected),
                            execution_failures=sum('FAIL' in str(r.get('run_status', '')) or 'ERROR' in str(r.get('run_status', '')) for r in selected),
                            data_passes=sum(r['data_pass'] for r in selected), image_passes=sum(r['image_pass'] for r in selected),
                            joint_passes=sum(r['joint_pass'] for r in selected),
                            material_data_passes=sum(r.get('material_data_pass', False) for r in selected),
                            material_data_pass_but_iou_fail=sum(r.get('material_data_pass', False) and not r['image_pass'] for r in selected),
                            data_pass_image_fail=sum(r['data_pass'] and not r['image_pass'] for r in selected),
                            median_nrmse=float(np.median([r['amplitude_normalized_nrmse'] for r in good])) if good else None,
                            median_iou=float(np.median([r['otsu_iou'] for r in good])) if good else None,
                            oracle_image_passes=sum(r.get('oracle_any_image_pass', False) for r in selected),
                            oracle_joint_passes=sum(r.get('oracle_any_joint_pass', False) for r in selected),
                            oracle_material_data_passes=sum(r.get('oracle_any_material_data_pass',False) for r in selected),
                            oracle_joint_coverage_runs=sum(r.get('oracle_joint_coverage_complete', False) for r in selected),
                            oracle_joint_passes_status='EXACT' if all(r.get('oracle_joint_coverage_complete', False) for r in selected) else 'LOWER_BOUND_MISSING_PHYSICAL_RESIDUAL',
                            oracle_median_best_nrmse=float(np.median([r['oracle_best_image_nrmse'] for r in oracle])) if oracle else None,
                            wall_seconds=sum(r.get('wall_seconds', 0.) or 0. for r in selected),
                            model_cost=dict(sum((Counter(r.get('model_cost', {})) for r in selected), Counter())),
                            actual_corrected_seed_counts=dict(Counter(str(r.get('actual_corrected_seeds', 'NOT_RECORDED')) for r in selected)),
                            object_cluster=cluster_summary(selected),
                            geometry_family_cluster=geometry_cluster_summary(selected)))
    return summary


def cost_rows(batch, rows):
    output = []
    for row in rows:
        counts = row.get('model_cost', {})
        aux = row.get('kernel_algebraic_counts', row.get('auxiliary_counts', {}))
        if not aux and isinstance(row.get('solver_cost'), dict):
            aux = row['solver_cost'].get('auxiliary', row['solver_cost'])
        auxiliary_status='RECORDED'
        if set(aux).issubset({'physical_caps'}) and str(row.get('method','')).startswith(('n1_','cie','tsom','csi')):
            auxiliary_status='NOT_RECORDED_AFTER_BUDGET_ABORT'
        def aux_count(key):return aux.get(key,0) if auxiliary_status=='RECORDED' else None
        qr=(aux.get('pivoted_QR',0)+aux.get('rank_revealing_qr',0)) if auxiliary_status=='RECORDED' else None
        output.append(dict(batch=batch, run_id=row['run_id'], group=group_id(row),
                           status=row.get('status'), stop=row.get('stop'), physical_caps=row.get('physical_caps',{}), execution_source=row.get('execution_source'), wall_seconds=row.get('wall_seconds'),
                           LU=counts.get('LU', 0), actual_complex_RHS=counts.get('actual_complex_RHS', 0),
                           source_RHS=counts.get('source_RHS', 0), receiver_transpose_RHS=counts.get('receiver_transpose_RHS', 0),
                           support_setup_RHS=counts.get('support_setup_RHS', 0),
                           auxiliary_accounting_status=auxiliary_status,QR=qr,
                           green_dense_matrix_products=aux_count('green_dense_matrix_products'),
                           green_vector_actions=aux_count('green_vector_actions'),
                           receiver_vector_actions=aux_count('receiver_vector_actions'),
                           actual_finite_seed_insertions=aux_count('finite_seed_insertions'),
                           rejected_finite_seed_branches=aux_count('seed_branches_rejected'),
                           green_basis_columns=aux_count('G_columns'),receiver_basis_columns=aux_count('S_columns'),
                           archived_endpoint_count=row.get('archived_endpoint_count'),
                           current_capacity_records=row.get('current_capacity_records', []),
                           seed_count_scope=row.get('seed_count_scope', 'Runner corrected_count'),
                           actual_proposed_seeds=row.get('proposal_count'), actual_corrected_seeds=row.get('corrected_count'),
                           phase_wall=row.get('phase_wall', {}), model_counts=counts, auxiliary_counts=aux,
                           setup_scope='Run costs as recorded; shared operator assembly is separate, not repeated per run.'))
    return output


def recorded_audit_cost_rows():
    """Recorded setup/failed-development charges only; never rerun an audit."""
    output = []
    def add(path,rid,wall,model,aux):
        output.append(dict(batch='recorded_audit_setup',run_id=rid,execution_source=str(path.resolve()),
                           group='recorded DEV/audit charge',status='RECORDED',wall_seconds=wall,
                           LU=model.get('LU',0),actual_complex_RHS=model.get('actual_complex_RHS',0),
                           model_counts=model,auxiliary_counts=aux,
                           setup_scope='Recorded receipt only; prior overwritten/unrecorded DEV attempts are not reconstructed. No new computation.'))
    root=STUDY/'audit'
    path=root/'literature_baseline_execution.json'
    if path.exists():
        record=json.loads(path.read_text())
        for i,row in enumerate(record.get('rows',[])):
            cost=row.get('cost',{})
            add(path,f'baseline_dev_{i}',cost.get('wall_seconds'),cost.get('model',{}),{k:v for k,v in cost.items() if k not in ('model','wall_seconds')})
    path=root/'literature_calibration_timing.json'
    if path.exists():
        record=json.loads(path.read_text())
        for field in ('timing_50','full_schedule'):
            for i,row in enumerate(record.get(field,[])):
                cost=row.get('cost',{})
                add(path,f'{field}_{i}',row.get('wall',cost.get('wall_seconds')),cost.get('model',{}),{k:v for k,v in cost.items() if k not in ('model','wall_seconds')})
    path=root/'E0_PRODUCTION_WIRING.json'
    if path.exists():
        record=json.loads(path.read_text())
        add(path,'E0_production_wiring',record.get('wall'),record.get('cost',{}),{})
    path=root/'N1_INITIAL_REPRESENTABILITY.json'
    if path.exists():
        record=json.loads(path.read_text())
        counts=dict(sum((Counter(r.get('cost',{})) for r in record.get('rows',[])),Counter()))
        add(path,'N1_initial_representability_all_cases',record.get('wall'),counts,{})
    output.append(dict(batch='recorded_audit_setup',run_id='N1_representability_first_loop_serialization_failed',
                       execution_source=str((root/'N1_INITIAL_REPRESENTABILITY.json').resolve()),
                       group='failed DEV serialization attempt',status='FAILED_SERIALIZATION_AFTER_PAID_PHYSICS',
                       wall_seconds=None,LU=24,actual_complex_RHS=144,model_counts=dict(LU=24,source_RHS=144,actual_complex_RHS=144),
                       auxiliary_counts={},count_provenance='Parent audited completed 24-forward-state loop before np.int64 JSON serialization failed',
                       setup_scope='Physics wall NOT_MEASURED; whole tool .243s includes imports and is not substituted. Subsequent proper diagnostic separately charged. No rerun.'))
    output.extend([
        dict(batch='recorded_native_DEV',run_id='owner_native_affine_DEV_trial',execution_source=str((root/'NATIVE_TSOM_FREEZE.json').resolve()),
             group='native-affine owner DEV',status='RECORDED_PARENT_RECEIPT',wall_seconds=.815661584,LU=2,actual_complex_RHS=30,
             source_RHS=12,tangent_RHS=18,model_counts=dict(LU=2,source_RHS=12,tangent_RHS=18,actual_complex_RHS=30),
             auxiliary_counts=None,auxiliary_accounting_status='NOT_RECORDED_IN_PARENT_COST_NOTE',
             setup_scope='One paid owner DEV trial, parent supplied measured wall and physical counters; no repeat.'),
        dict(batch='recorded_native_DEV',run_id='owner_native_weighted_wiring_physical_reference',execution_source=str((root/'NATIVE_TSOM_FREEZE.json').resolve()),
             group='native-affine weighted wiring',status='RECORDED_PARENT_RECEIPT',wall_seconds=None,LU=1,actual_complex_RHS=6,
             source_RHS=6,model_counts=dict(LU=1,source_RHS=6,actual_complex_RHS=6),auxiliary_counts=None,
             auxiliary_accounting_status='NOT_RECORDED_IN_PARENT_COST_NOTE',setup_scope='Separate paid physical reference; wall NOT_MEASURED; no repeat.')])
    native_gen=STUDY/'data/NATIVE_CONFIRMATION_FREEZE.json'
    if native_gen.exists():
        output.append(dict(batch='shared_setup',run_id='native_confirmation_generation',execution_source=str(native_gen.resolve()),
             generation_receipt=str((root/'NATIVE_CONFIRMATION_GENERATION.json').resolve()),
             group='native-affine independent confirmation generation',status='RECORDED_GENERATION_COUNTS_AND_WALL',
             wall_seconds=.008774708025157452,LU=12,actual_complex_RHS=72,source_RHS=72,
             model_counts=dict(LU=12,source_RHS=72,actual_complex_RHS=72),auxiliary_counts=None,
             setup_scope='Parent supplied generation charge: 6 objects ×2 layouts, clean/noisy share the generated clean field. No second generation or new physics run.'))
    return output


def recorded_external_setup_rows():
    output=[]
    for path in (STUDY/'results').glob('*_scores.json'):
        record=json.loads(path.read_text())
        for source in record.get('result_sources',[]):
            if not source.get('shared_setup'): continue
            counts=Counter()
            for entry in source['shared_setup']:
                counts.update({k:v for k,v in entry.get('counts',{}).items() if isinstance(v,(int,float)) and k!='wall_seconds'})
            rhs=sum(counts[k] for k in ('source_RHS','adjoint_RHS','receiver_RHS','tangent_RHS','support_RHS'))
            output.append(dict(batch='external_shared_setup',run_id='shared_operators_and_CIE_basis',
                               execution_source=str(Path(source['path']).resolve()),group='E3 shared setup',status='RECORDED',
                               wall_seconds=source.get('setup_wall_seconds'),LU=counts.get('LU',0),actual_complex_RHS=rhs,
                               model_counts={k:counts.get(k,0) for k in ('LU','source_RHS')},
                               auxiliary_counts=dict(counts),source_receipts=source['shared_setup'],
                               setup_scope='Setup charged once per finite queue; recorded RHS sum, 580 G/S basis columns and SVDs retained explicitly.'))
    failed_dir=STUDY/'results/e3_limited_diagnostic_20261101'
    if failed_dir.exists():
        output.append(dict(batch='unmeasured_failed_launches',run_id='E3_failed_launches_aggregate_UNMEASURED',
                           execution_source=str(failed_dir.resolve()),group='E3 failed launch overhead',
                           status='FAILED_LAUNCH_COST_UNMEASURED',wall_seconds=None,LU=None,actual_complex_RHS=None,
                           setup_scope='Parent disclosed failed launches; exact wall/physical cost not measured in downloaded receipts. Unknown, not zero; aggregate does not imply a measured attempt count.'))
    return output


def deduplicate_cost_ledger(rows):
    """Count reused sealed executions once, retaining genuinely repeated batches."""
    source_by_tag_run = {}
    raw_by_execution = {}
    for path in (STUDY/'results').glob('*_scores.json'):
        record = json.loads(path.read_text())
        tag = record.get('tag',record.get('batch'))
        sources = record.get('result_sources', [])
        if not sources and record.get('result_source'): sources=[dict(path=record['result_source'])]
        for source in sources:
            source_path = Path(source['path'])
            if not source_path.is_absolute(): source_path=Path.cwd()/source_path
            if not source_path.exists(): continue
            artifact=json.loads(source_path.read_text())
            for raw in artifact.get('rows',[]):
                source_dir=str(source_path.parent.resolve())
                source_by_tag_run[(tag,raw['run_id'])]=source_dir
                raw_by_execution[(source_dir,raw['run_id'])]=raw
    unique = {}
    for row in rows:
        row=dict(row)
        source=row.get('execution_source') or source_by_tag_run.get((row['batch'],row['run_id']))
        row['execution_source']=source
        key=(source or row['batch'],row['run_id'])
        tags=row.get('included_score_tags',[])
        if isinstance(tags,str):
            try: tags=json.loads(tags)
            except ValueError: tags=[]
        tags=set(tags+[row['batch']])
        if key in unique:
            old=unique[key]
            tags.update(old['included_score_tags'])
            first_batch=old['batch']
            old.update(row)
            old['batch']=first_batch
            old['included_score_tags']=sorted(tags)
        else:
            row['included_score_tags']=sorted(tags)
            unique[key]=row
    for key,row in unique.items():
        raw=raw_by_execution.get(key,{})
        if str(raw.get('method','')).startswith(('n1_','cie','tsom','csi')):
            cost=raw.get('solver_cost')
            if isinstance(cost,dict) and set(cost).issubset({'physical_caps'}):
                row['auxiliary_accounting_status']='NOT_RECORDED_AFTER_BUDGET_OR_EXECUTION_ABORT'
                for field in ('QR','green_dense_matrix_products','green_vector_actions','receiver_vector_actions','green_basis_columns','receiver_basis_columns','actual_finite_seed_insertions','rejected_finite_seed_branches'):
                    row[field]=None
                row['auxiliary_counts']=cost
    return list(unique.values())


def multistart_curves(scored,report_k=None):
    groups = defaultdict(list)
    for row in scored: groups[row['group']].append(row)
    outputs = []
    for group, rows in sorted(groups.items()):
        by_case = defaultdict(list)
        for row in rows: by_case[row['case_id']].append(row)
        starts = sorted({r.get('start',0) for r in rows})
        for k in range(1,len(starts)+1):
            if report_k and k not in report_k:continue
            selected_rows = []
            for case, case_rows in sorted(by_case.items()):
                pool = [r for r in case_rows if r.get('start',0) in starts[:k]]
                if len(pool) != k: continue
                def score(r):
                    v = r.get('online_physical_score')
                    return v if v is not None else np.inf
                online = min(pool,key=lambda r:(score(r),r['run_id']))
                result = dict(online)
                result.update(K=k, prefix_starts=starts[:k], selected_online_run_id=online['run_id'],
                    cumulative_run_ids=[r['run_id'] for r in pool],
                    oracle_final_endpoint_best_nrmse=min((r.get('amplitude_normalized_nrmse',np.inf) for r in pool),default=None),
                    oracle_final_endpoint_any_joint_pass=any(r['joint_pass'] for r in pool),
                    oracle_final_endpoint_any_material_data_pass=any(r.get('material_data_pass',False) for r in pool),
                    oracle_all_candidates_any_material_data_pass=any(r.get('oracle_any_material_data_pass',False) for r in pool),
                    oracle_all_candidates_any_joint_pass=any(r.get('oracle_any_joint_pass',False) for r in pool),
                    actual_corrected_seeds=sum(r.get('actual_corrected_seeds',0) or 0 for r in pool),
                    wall_seconds=sum(r.get('wall_seconds',0.) or 0. for r in pool),
                    model_cost=dict(sum((Counter(r.get('model_cost',{})) for r in pool),Counter())))
                if all(r.get('physical_caps') for r in pool):
                    result['physical_caps']={key:(sum(r['physical_caps'][key] for r in pool) if all(r['physical_caps'].get(key) is not None for r in pool) else None) for key in ('max_lu','max_rhs')}
                    result['physical_cap_scope']='Cumulative caps across paid starts in prefix'
                selected_rows.append(result)
            if not selected_rows: continue
            summary = summarize(selected_rows)[0]
            outputs.append(dict(group=group,K=k,starts=starts[:k],observations=len(selected_rows),
                                summary=summary,online_best_phi_rows=selected_rows,
                                oracle_final_endpoint_joint_passes=sum(r['oracle_final_endpoint_any_joint_pass'] for r in selected_rows),
                                oracle_all_candidates_joint_passes=sum(r['oracle_all_candidates_any_joint_pass'] for r in selected_rows),
                                oracle_final_endpoint_material_data_passes=sum(r['oracle_final_endpoint_any_material_data_pass'] for r in selected_rows),
                                oracle_all_candidates_material_data_passes=sum(r['oracle_all_candidates_any_material_data_pass'] for r in selected_rows)))
    return outputs


def multistart_plot(tag, curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    grouped = defaultdict(list)
    for row in curves: grouped[row['group']].append(row)
    fig, axes = plt.subplots(1,3,figsize=(16,5.5),layout='constrained')
    for i,(group,points) in enumerate(sorted(grouped.items())):
        points.sort(key=lambda r:r['K'])
        label = group.split(' / ')[0]
        color = plt.get_cmap('tab20')(i%20)
        yy = [r['summary']['joint_passes']/r['observations'] for r in points]
        axes[0].plot([r['K'] for r in points],yy,'o-',label=label,color=color,markersize=3)
        axes[1].plot([r['summary']['wall_seconds'] for r in points],yy,'o-',color=color,markersize=3)
        axes[2].plot([r['summary']['model_cost'].get('LU',0) for r in points],yy,'o-',color=color,markersize=3)
    axes[0].set(xlabel='Shared start count K (frozen prefix order)',ylabel='Online best-Φ joint PASS fraction')
    axes[1].set(xlabel='Actual sum inference wall over scored observations / s')
    axes[2].set(xlabel='Actual sum full physical LU over scored observations')
    for ax in axes: ax.set_ylim(-.02,1.02);ax.grid(alpha=.2)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=7,fontsize=8)
    actual_ks=sorted({r['K'] for r in curves});max_k=max(actual_ks,default=1)
    fig.suptitle(f'Actual multi-start prefixes K={actual_ks}; finite fixed shape families')
    fig.savefig(STUDY/'figures'/f'{tag}_K1_K{max_k}_cost_curves.png',dpi=165)
    plt.close(fig)


def short_method_label(group):
    method=group.split(' / ')[0]
    labels={'native':'Native','native_trf':'Native TRF','feedback':'Fixed feedback','feedback_beta:0.2':'Fixed feedback β=.2','cie_dual':'CIE dual','cie_single':'CIE single','cie_matched':'CIE matched mf2','material_seeds':'Material seeds','globaldual':'GlobalDual','reciprocal':'Reciprocal','data':'Data L-BFGS','tsom':'TSOM','csi':'CSI'}
    if method in labels:return labels[method]
    return method.replace('_chart_cie_native_tsom','').replace('_chart_cie','').replace('n1_','N1 ').replace('_',' ')


def plots(batch, scored, summaries, mapping, truths, materials, raw_rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    origin = 'upper' if batch.startswith('e3') else 'lower'
    raster_note = '; row0 at top, original raster order' if origin == 'upper' else ''
    figures = STUDY / 'figures'
    figures.mkdir(exist_ok=True)
    labels = [s['group'] for s in summaries]
    width = max(10, len(labels) * .9)
    fig, axes = plt.subplots(1, 3, figsize=(width, 4.4), layout='constrained')
    x = np.arange(len(summaries))
    for i, (field, title) in enumerate((('data_passes', 'Physical data PASS'), ('joint_passes', 'Joint image + data PASS'), ('median_nrmse', 'Median amplitude-normalized NRMSE'))):
        vals = [(s[field] / s['runs'] if i < 2 else (s[field] if s[field] is not None else np.nan)) for s in summaries]
        axes[i].bar(x, vals, color='#267c9e')
        axes[i].set_xticks(x, labels, rotation=65, ha='right', fontsize=7)
        axes[i].set_title(title)
        if i < 2: axes[i].set_ylim(0, 1)
        else: axes[i].axhline(.15, color='crimson', linestyle='--')
    fig.suptitle(f'{batch.upper()}: all sealed runs; finite object clusters, shared starts/layouts')
    fig.savefig(figures / f'{batch}_scores_summary.png', dpi=165)
    plt.close(fig)
    # All-object visual panels: selection uses only the online physical score.
    raw_lookup = {r['run_id']: r for r in raw_rows}
    all_groups = sorted({r.get('comparison_group',r['group']) for r in scored})
    cases = sorted({r['case_id'] for r in scored})
    panels=sorted({(mapping[c]['side'],mapping[c]['layout'],str(mapping[c].get('noise','clean'))) for c in cases})
    for side,layout,condition in panels:
        case_ids=sorted([c for c in cases if mapping[c]['side']==side and mapping[c]['layout']==layout and str(mapping[c].get('noise','clean'))==condition],
                        key=lambda c:(mapping[c]['shape'],mapping[c]['amplitude']))
        per_page=9 if batch=='e1_native_confirmation' else 8
        for page in range((len(all_groups)+per_page-1)//per_page):
            page_groups=all_groups[page*per_page:(page+1)*per_page]
            fig, axes = plt.subplots(len(case_ids), len(page_groups)+1,
                figsize=(2.45*(len(page_groups)+1), 2.55*len(case_ids)), squeeze=False, layout='constrained')
            for ri, case in enumerate(case_ids):
                meta = mapping[case]
                truth = truths[meta['object_id']]['m'].reshape(side, side)
                axes[ri, 0].imshow(truth, origin=origin, vmin=.02, vmax=20, cmap='viridis', interpolation='nearest')
                axes[ri, 0].set_title(f"GT: {meta['shape']}, amp={meta['amplitude']:g}\n{condition}; {case[-8:]}", fontsize=9)
                for ci, group in enumerate(page_groups, 1):
                    possible = [r for r in scored if r['case_id'] == case and r.get('comparison_group',r['group']) == group and r['run_id'] in materials]
                    if not possible:
                        axes[ri, ci].text(.5,.5,'FAILED / MISSING',ha='center',va='center',transform=axes[ri,ci].transAxes)
                        continue
                    def online_score(row):
                        raw = raw_lookup[row['run_id']]
                        value = raw.get('score', raw.get('online_score', np.inf))
                        return value if value is not None else np.inf
                    chosen = min(possible, key=lambda r: (online_score(r), r['run_id']))
                    im = axes[ri, ci].imshow(materials[chosen['run_id']].reshape(side, side), origin=origin,
                                            vmin=.02, vmax=20, cmap='viridis', interpolation='nearest')
                    axes[ri, ci].set_title(f"{short_method_label(group)}\ns={chosen.get('start','?')}; N={chosen['amplitude_normalized_nrmse']:.3f}; IoU={chosen['otsu_iou']:.2f}\nr={chosen['relative_residual']:.2g}; {'PASS' if chosen['joint_pass'] else ('IoU FAIL' if chosen.get('material_data_pass') else 'FAIL')}", fontsize=7)
                for ax in axes[ri]: ax.set_xticks([]); ax.set_yticks([])
            fig.suptitle(f'{batch.upper()} / {side}×{side} / {layout} / {condition}: lowest online physical score across shared starts; scale [.02,20]{raster_note}', fontsize=11)
            fig.savefig(figures / f'{batch}_materials_side{side}_{layout}_{condition}_p{page+1}.png', dpi=150)
            plt.close(fig)


def fullwave_publication_plots(scored,mapping,truths,materials):
    """The requested raster figures: values unchanged, row zero at top."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    methods=['native_trf','feedback_beta:0.2','cie_single','cie_dual']
    labels=['Native TRF','Fixed feedback β=.2','CIE single + corrector','CIE dual + corrector']
    cases=sorted({r['case_id'] for r in scored},key=lambda c:(mapping[c]['shape'],mapping[c]['layout'],c))
    def digit_order(c):
        text=mapping[c]['shape']
        return int(''.join(x for x in text if x.isdigit()) or 99)
    digits=sorted([c for c in cases if mapping[c]['shape'].startswith('digit') and 1<=digit_order(c)<=9],key=lambda c:(digit_order(c),mapping[c]['layout']))
    lookup={(r['case_id'],r['method']):r for r in scored}
    def draw(selected,path,errors=False):
        fig,axes=plt.subplots(len(selected),5,figsize=(12.8,2.45*len(selected)),squeeze=False,layout='constrained')
        material_image=error_image=None
        for ri,case in enumerate(selected):
            meta=mapping[case];side=meta['side'];gt=truths[meta['object_id']]['m'].reshape(side,side)
            material_image=axes[ri,0].imshow(gt,origin='upper',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest')
            axes[ri,0].set_title(f"GT {meta['shape']} / {meta['layout']}\n{case[-8:]}",fontsize=9)
            for ci,(method,label) in enumerate(zip(methods,labels),1):
                row=lookup.get((case,method))
                if row is None or row['run_id'] not in materials:
                    axes[ri,ci].text(.5,.5,'FAILED / MISSING',ha='center',va='center',transform=axes[ri,ci].transAxes)
                    continue
                mm=materials[row['run_id']].reshape(side,side)
                if errors:
                    error_image=axes[ri,ci].imshow(abs(mm-gt),origin='upper',vmin=0,vmax=20,cmap='magma',interpolation='nearest')
                else:
                    axes[ri,ci].imshow(mm,origin='upper',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest')
                axes[ri,ci].set_title(f"{label}\nN={row['amplitude_normalized_nrmse']:.3f}; IoU={row['otsu_iou']:.2f}\nr={row['relative_residual']:.2g}; {'PASS' if row['joint_pass'] else 'FAIL'}",fontsize=8)
            for ax in axes[ri]:ax.set_xticks([]);ax.set_yticks([])
        fig.colorbar(material_image,ax=axes[:,0],shrink=.2,label='Material m [.02,20]')
        if errors and error_image is not None:fig.colorbar(error_image,ax=axes[:,1:].ravel().tolist(),shrink=.2,label='Absolute material error [0,20]')
        fig.suptitle(f"Disclosed historical 24×24 targets; {'absolute-error panels' if errors else 'original material arrays'}; row0 at top, original raster order",fontsize=13)
        fig.savefig(STUDY/'figures'/path,dpi=165)
        plt.close(fig)
    draw(cases,'fullwave_all_shapes.png')
    if digits:
        draw(digits,'fullwave_digits_1to9.png')
        draw(digits,'fullwave_digit_errors.png',True)


def load_sealed_results(result_path):
    closed = json.loads(result_path.read_text())
    if closed.get('scope') == 'LIMITED_E3_MATURE_BASELINE_DIAGNOSTIC':
        tasks=closed.get('manifest',{}).get('tasks',[])
        if closed.get('tasks') != len(tasks) or not all((result_path.parent/(t['task_id']+'.json')).exists() for t in tasks):
            raise SystemExit('E3 finite queue is incomplete; no private scoring allowed')
        closed['original_status_field']=closed.get('status','ABSENT')
        closed['status']='PREDICTIONS_CLOSED_COMPLETE_FINITE_QUEUE_PARENT_AUTHORIZED'
        raw_rows=[]
        materials,predictions,candidates,endpoints={},{},{},{}
        for task in tasks:
            receipt=json.loads((result_path.parent/(task['task_id']+'.json')).read_text())
            result=receipt.get('result',{})
            rid=task['task_id']
            row=dict(run_id=rid,case_id=task.get('case_id',task.get('observation_id')),method=task['method'],
                     start=task.get('start',int(rid.rsplit('__s',1)[-1]) if '__s' in rid else 0),
                     side=int(np.sqrt(receipt.get('all_material_cells_unknown',576))),layout=task.get('layout','FROM_EVALUATION_METADATA'),
                     status=receipt.get('status'),relative_residual=result.get('relative_residual'),
                     score=result.get('ranking_objective'),wall_seconds=receipt.get('wall_seconds'),
                     wall_budget=receipt.get('budget_envelope_seconds',closed['manifest'].get('config',{}).get('wall_seconds')),
                     model_cost=receipt.get('physical_counts',receipt.get('partial_physical_counts',{})),
                     auxiliary_counts=receipt.get('current_costs',{}),execution_source=str(result_path.parent.resolve()),
                     selected_phase=result.get('selected_phase'),
                     current_schedule_status=(result.get('current') or {}).get('status'),
                     physical_residual_source='ALREADY_PAID_ORIGINAL_MAXWELL_RECEIPT_NO_NEW_SOLVE',
                     closed_candidate_material_keys=[])
            worker=result.get('corrector') or result
            row.update(worker_wrapper_status=worker.get('wrapper_status'),worker_wrapper_stop_reason=worker.get('wrapper_stop_reason'),
                       worker_phase_stop_reasons=[x.get('wrapper_stop_reason') for x in worker.get('phases',[])],
                       raw_declared_max_lu=closed['manifest'].get('config',{}).get('max_lu'),
                       raw_declared_max_rhs=closed['manifest'].get('config',{}).get('max_rhs'))
            if result.get('m') is not None: materials[rid]=np.asarray(result['m'],float)
            current=result.get('current') or {}
            if current:
                widths=[h['mf'] for h in current.get('history',[]) if h.get('mf') is not None]
                row.update(current_schedule_complete=current.get('status')=='SCHEDULE_COMPLETE',current_max_mf_reached=max(widths,default=None))
            endpoint_pool=list(current.get('candidates',[]))
            correction=result.get('corrector')
            if correction and correction.get('m') is not None: endpoint_pool.append(correction)
            if current.get('candidates'):
                aux_best=min(current['candidates'],key=lambda c:c['objective'])
                row['auxiliary_current_relative_residual']=aux_best.get('physical_relative_residual')
            for i,c in enumerate(endpoint_pool):
                key=rid+f'_candidate{i}'
                candidates[key]=np.asarray(c['m'],float)
                endpoints[key]=dict(relative_residual=c.get('physical_relative_residual',c.get('relative_residual')),
                                    provenance=dict(stage=c.get('stage'),reason=c.get('reason')))
                row['closed_candidate_material_keys'].append(key)
            row['corrected_count']=int(bool(correction))
            row['seed_count_scope']='One shared current incumbent correction for CIE; native/feedback main solves separate'
            raw_rows.append(row)
        return closed,raw_rows,materials,predictions,candidates,endpoints
    if 'CLOSED' not in closed.get('status', ''):
        raise SystemExit(f'Result artifact has no closed-predictions state: {result_path}')
    prefix = Path(str(result_path).removesuffix('_results.json'))
    raw_rows = closed['rows']
    for raw in raw_rows:
        raw['execution_source'] = str(result_path.parent.resolve())
    if raw_rows and (result_path.parent/(raw_rows[0]['run_id']+'.npz')).exists():
        materials, predictions, candidates, endpoints = {}, {}, {}, {}
        for raw in raw_rows:
            rid = raw['run_id']
            saved = npz_dict(result_path.parent / (rid+'.npz'))
            if 'm' in saved: materials[rid] = saved['m']
            if 'prediction' in saved: predictions[rid] = saved['prediction']
            detail_path = result_path.parent / (rid+'.json.gz')
            metadata = []
            if detail_path.exists():
                with gzip.open(detail_path, 'rt') as file:
                    detail = json.load(file)
                metadata = detail.get('candidates', [])
                raw['current_capacity_records'] = [
                    {k:h[k] for k in ('tag','event','stage','capacity_before','capacity_after','real_current_dof','apcs_rank','fourier_columns') if k in h}
                    for h in detail.get('history', []) if h.get('tag') == 'mode_selected' or h.get('event') == 'stage_end']
            common_seed_count = sum(str(c.get('tag', '')).startswith('common_corrector_') for c in metadata)
            raw.setdefault('corrected_count', common_seed_count)
            raw['archived_endpoint_count'] = len(saved.get('candidates', []))
            raw['seed_count_scope'] = 'Common auxiliary full-material corrector starts; native primary optimizer is separate'
            keys = []
            for i, cm in enumerate(saved.get('candidates', [])):
                key = rid+f'_candidate{i}'
                candidates[key] = cm
                keys.append(key)
                em = metadata[i] if i < len(metadata) else {}
                residual = em.get('relative_residual', em.get('physical_relative_residual'))
                if residual is None and str(em.get('tag', '')).startswith('common_corrector_'):
                    incumbent_records = [h for h in em.get('history', []) if h.get('physical_score') is not None and h.get('relative_residual') is not None]
                    if incumbent_records:
                        residual = min(incumbent_records, key=lambda h: h['physical_score'])['relative_residual']
                # The native corrector's incumbent archive stores physical score;
                # recover its already-paid residual from the matching history.
                if residual is None and em.get('physical_score') is not None:
                    for h in detail.get('history', []):
                        if h.get('physical_score') == em['physical_score']:
                            residual = h.get('relative_residual')
                            break
                endpoints[key] = dict(relative_residual=residual, provenance=em.get('tag', em.get('provenance')))
            raw['closed_candidate_material_keys'] = keys
            raw.setdefault('score', raw.get('physical_score'))
    else:
        materials = npz_dict(Path(str(prefix)+'_materials.npz'))
        predictions = npz_dict(Path(str(prefix)+'_predictions.npz'))
        candidates = npz_dict(Path(str(prefix)+'_candidate_materials.npz'))
        endpoints = endpoint_metadata(prefix)
    return closed, raw_rows, materials, predictions, candidates, endpoints


def native_affine_gain_plot(tag,scored,curves,mapping,truths,materials):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    chosen={x['group'].split(' / ')[0]:{a['case_id']:a for a in x['online_best_phi_rows']} for x in curves if x['K']==5}
    ordinary='n1_fourier_chart_cie_native_tsom'
    if ordinary not in chosen:return
    candidates=sorted({c for method,items in chosen.items() if method.startswith('n1_finite') for c,a in items.items()
                       if a['joint_pass'] and not chosen[ordinary][c]['joint_pass']})
    if not candidates:return
    case=candidates[0];meta=mapping[case];side=meta['side'];gt=truths[meta['object_id']]['m'].reshape(side,side)
    methods=[None,'cie_dual','cie_matched','feedback',ordinary,'n1_task_chart_cie_native_tsom',
             'n1_finite_chart_cie_native_tsom','n1_finite_fourier_chart_cie_native_tsom',
             'n1_finite_random_chart_cie_native_tsom','n1_finite_krylov_chart_cie_native_tsom']
    labels=['GT','Mature CIE dual','Matched CIE mf2','Fixed feedback','Ordinary Fourier GN','Task local GN',
            'Finite task','Finite Fourier','Finite random','Finite Krylov']
    fig,axes=plt.subplots(2,5,figsize=(13.5,6),layout='constrained')
    for ax,method,label in zip(axes.ravel(),methods,labels):
        if method is None:
            image=ax.imshow(gt,origin='lower',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest');ax.set_title(label)
        else:
            row=chosen[method][case]
            ax.imshow(materials[row['run_id']].reshape(side,side),origin='lower',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest')
            ax.set_title(f"{label}; s={row['start']}\nN={row['amplitude_normalized_nrmse']:.4f}; IoU={row['otsu_iou']:.2f}\nr={row['relative_residual']:.2g}; {'PASS' if row['joint_pass'] else 'FAIL'}",fontsize=9)
        ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(image,ax=axes.ravel().tolist(),shrink=.6,label='Physical material [.02,20]')
    fig.suptitle(f"Native-affine development: offline illustration of the shared finite-branch gain\n{meta['shape']}, amplitude={meta['amplitude']:g}, {meta['layout']}; {case}; all-case statistics reported separately",fontsize=11)
    fig.savefig(STUDY/'figures'/f'{tag}_finite_gain_case.png',dpi=165)
    plt.close(fig)


def native_confirmation_gain_plot(record,mapping,truths,materials):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    chosen={c['group'].split(' / ')[0]:{a['case_id']:a for a in c['online_best_phi_rows']}
            for c in record['multistart_K_prefix_curves'] if c['K']==3 and 'noise_label=clean' in c['group']}
    task='n1_task_chart_cie_native_tsom';fourier='n1_fourier_chart_cie_native_tsom'
    extra=sorted(c for c,a in chosen[task].items() if a['material_data_pass'] and not chosen[fourier][c]['material_data_pass'])
    if not extra:return
    case=extra[0];meta=mapping[case];side=meta['side'];gt=truths[meta['object_id']]['m'].reshape(side,side)
    methods=[None,'native','feedback','cie_dual',fourier,task,'n1_finite_chart_cie_native_tsom','n1_finite_random_chart_cie_native_tsom']
    labels=['GT','Native','Fixed feedback','Mature CIE dual','Ordinary Fourier','Task-selected mode','Finite task','Finite random']
    fig,axes=plt.subplots(2,4,figsize=(13.2,7.7),layout='constrained')
    self_iou=record['offline_gt_gate_diagnostics'][meta['object_id']]['gt_self_otsu_iou']
    for ax,method,label in zip(axes.ravel(),methods,labels):
        if method is None:
            image=ax.imshow(gt,origin='lower',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest')
            ax.set_title(f'GT curved, amplitude16\nGT self-Otsu IoU={self_iou:.3f}',fontsize=10)
        else:
            row=chosen[method][case];cost=row['model_cost']
            ax.imshow(materials[row['run_id']].reshape(side,side),origin='lower',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest')
            state=('material+data PASS; joint IoU FAIL' if row['material_data_pass'] and not row['joint_pass'] else ('joint PASS' if row['joint_pass'] else 'FAIL'))
            ax.set_title(f"{label}; chosen s={row['start']}\nN={row['amplitude_normalized_nrmse']:.5f}; IoU={row['otsu_iou']:.3f}; r={row['relative_residual']:.2g}\n{state}\nK3 sum: {row['wall_seconds']:.2f}s, {cost.get('LU',0)} LU, {cost.get('actual_complex_RHS',0)} RHS",fontsize=8.5)
        ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(image,ax=axes.ravel().tolist(),shrink=.65,label='Physical material [.02,20]')
    fig.suptitle(f'Fresh native-affine confirmation: one clean material+data gain; primary joint gate unchanged\n{case}; uniform; K3 starts0/1/4; per-start caps4s /120 LU /2400 RHS; N1 current DOF192\nOffline illustrative case; corresponding noisy observation not recovered by task',fontsize=11)
    fig.savefig(STUDY/'figures/e1_native_confirmation_task_material_gain.png',dpi=165)
    plt.close(fig)


def native_final_increment(records):
    index={r.get('tag',r.get('batch')):r for _,r in records}
    record=index.get('e1_native_confirmation')
    if not record:return None
    output=dict(evaluation='FRESH_BLIND_NATIVE_AFFINE_CONFIRMATION',primary_gate_changed=False,by_noise={},N1_E3='NOT_RUN',N2_E3='NOT_RUN')
    for noise in ('clean','noisy'):
        chosen={c['group'].split(' / ')[0]:c for c in record.get('multistart_K_prefix_curves',[])
                if c['K']==3 and 'noise_label='+noise in c['group']}
        sets={metric:{method:{a['case_id'] for a in c['online_best_phi_rows'] if a.get(metric)} for method,c in chosen.items()}
              for metric in ('data_pass','material_data_pass','joint_pass')}
        task='n1_task_chart_cie_native_tsom';fourier='n1_fourier_chart_cie_native_tsom'
        baseline_methods=('cie_dual','feedback','material_seeds','native')
        base=set().union(*(sets['material_data_pass'].get(m,set()) for m in baseline_methods))
        output['by_noise'][noise]=dict(
            counts={m:{metric:len(sets[metric][m]) for metric in sets} for m in chosen},
            task_primary_joint_gains_vs_fourier=sorted(sets['joint_pass'][task]-sets['joint_pass'][fourier]),
            task_primary_joint_losses_vs_fourier=sorted(sets['joint_pass'][fourier]-sets['joint_pass'][task]),
            task_secondary_material_data_gains_vs_fourier=sorted(sets['material_data_pass'][task]-sets['material_data_pass'][fourier]),
            task_secondary_material_data_gains_vs_strong_portfolio=sorted(sets['material_data_pass'][task]-base))
    post=index.get('e1_posthoc_mode_controls')
    if post:
        output['post_disclosure_component_scope']=post.get('classification')
        output['post_disclosure_clean_material_data_successes']=[dict(method=a['method'],case_id=a['case_id'],start=a.get('start'))
            for a in post.get('rows',[]) if a.get('noise_label')=='clean' and a.get('material_data_pass')]
    output['interpretation']='Primary joint coverage does not improve; one clean secondary material+data gain survives immutable metrics, but normal/soft also recover the case in post-disclosure diagnostics. Noise and cross-object robustness and A17-specific causal attribution are not established.'
    return output


def write_global_aggregate():
    """Export derived score tables; no truth reads, solves, or altered metrics."""
    results=STUDY/'results'
    records=[]
    for path in sorted(results.glob('*_scores.json')):
        r=json.loads(path.read_text())
        if r.get('status')=='POST_CLOSURE_SCORED':records.append((path,r))
    cohorts=[];by_case=[];rejected=[]
    for path,r in records:
        tag=r.get('tag',r.get('batch',path.stem))
        cohorts.append(dict(tag=tag,classification=r.get('classification'),score_file=str(path.relative_to(STUDY)),
                            scored_rows=len(r.get('rows',[])),material_contracts=len({a['object_id'] for a in r.get('rows',[])}),
                            summaries=r.get('summaries'),source_artifacts=r.get('result_sources'),
                            historical_subset_of='e3_all' if tag=='e3_pilot' else None,
                            truth_access_deviation=r.get('truth_access_deviation',r.get('truth_access_deviation_record'))))
        curves=r.get('multistart_K_prefix_curves',[])
        if curves:
            selected=[(c['K'],a) for c in curves if c['K'] in (1,5) or ('confirmation' in tag and c['K']==3)
                      for a in c['online_best_phi_rows']]
        else:selected=[(None,a) for a in r.get('rows',[])]
        for k,a in selected:
            by_case.append(dict(score_cohort=tag,classification=r.get('classification'),K=k,
                case_id=a['case_id'],object_id=a['object_id'],shape=a['shape'],amplitude=a['amplitude'],
                side=a.get('side'),layout=a.get('layout'),noise=a.get('noise_label','clean'),method=a['method'],
                selected_run_id=a['run_id'],selected_start=a.get('start'),run_status=a.get('run_status'),
                data_pass=a['data_pass'],material_data_pass=a.get('material_data_pass',a['data_pass'] and a.get('amplitude_normalized_nrmse',float('inf'))<=.15),
                joint_pass=a['joint_pass'],nrmse=a.get('amplitude_normalized_nrmse'),iou=a.get('otsu_iou'),
                relative_residual=a.get('relative_residual'),physical_data_gate=a.get('physical_data_gate'),
                original_online_physical_score=a.get('online_physical_score'),
                cumulative_worker_wall=a.get('wall_seconds'),LU=a.get('model_cost',{}).get('LU'),
                actual_complex_RHS=a.get('model_cost',{}).get('actual_complex_RHS'),
                oracle_final_endpoint_joint=a.get('oracle_final_endpoint_any_joint_pass'),
                oracle_available_candidates_joint=a.get('oracle_all_candidates_any_joint_pass',a.get('oracle_any_joint_pass'))))
        for c in r.get('candidate_scores',[]):
            if c.get('joint_pass'):continue
            if c.get('relative_residual') is None:reason='PHYSICAL_RESIDUAL_NOT_RECORDED'
            elif not c.get('data_pass'):reason='OFFLINE_PHYSICAL_DATA_GATE_FAIL'
            elif c.get('amplitude_normalized_nrmse',float('inf'))>.15:reason='OFFLINE_MATERIAL_GATE_FAIL'
            else:reason='OFFLINE_IOU_GATE_FAIL'
            rejected.append(dict(score_cohort=tag,run_id=c['run_id'],candidate_key=c['candidate_key'],
                category=reason,scope='Offline gate result; not automatically an online rejection',
                nrmse=c.get('amplitude_normalized_nrmse'),iou=c.get('otsu_iou'),relative_residual=c.get('relative_residual'),
                provenance=c.get('provenance')))
        for a in r.get('rows',[]):
            if not a['joint_pass']:
                rejected.append(dict(score_cohort=tag,run_id=a['run_id'],candidate_key='SELECTED_ONLINE_ENDPOINT',
                    category='SELECTED_ENDPOINT_JOINT_FAIL',scope='Selected sealed result, not a fabricated candidate',
                    nrmse=a.get('amplitude_normalized_nrmse'),iou=a.get('otsu_iou'),relative_residual=a.get('relative_residual'),
                    provenance=a.get('joint_failure_cause',a.get('score_status'))))
    # N2 root filters are different from post-corrector GT gates; preserve them.
    path=results/'e2_mechanism_ledger.jsonl.gz'
    if path.exists():
        with gzip.open(path,'rt') as file:
            for line in file:
                run=json.loads(line)
                for i,root in enumerate(run.get('roots',[])):
                    status=root.get('status')
                    if status=='LEGAL_SUBSET_ROOT':status='LEGAL_SUBSET_ROOT_FULL_DATA_FAIL'
                    rejected.append(dict(score_cohort='e2_root_filter',run_id=run['run_id'],candidate_key=f'root_{i}',
                        category=status,scope='Numerical conditional-root filter; legal subset roots are not full-data solutions',
                        relative_residual=root.get('all_data_relative_residual'),subset_relative_residual=root.get('subset_relative_residual'),
                        provenance={k:root.get(k) for k in ('reason','support','receiver','sources','certificate') if k in root}))
    with (results/'COST_LEDGER.csv').open() as file:ledger=list(csv.DictReader(file))
    setup_batches={'recorded_audit_setup','recorded_native_DEV','external_shared_setup','shared_setup','unmeasured_failed_launches'}
    inferred=[a for a in ledger if a['batch'] not in setup_batches]
    summary=dict(status='DERIVED_FROM_CLOSED_SCORED_ARTIFACTS',cohorts=cohorts,final_native_increment=native_final_increment(records),
        unique_paid_inference_executions=len(inferred),cost_ledger_rows=len(ledger),
        known_recorded_LU=sum(int(a['LU']) for a in ledger if a.get('LU')),
        known_recorded_complex_RHS=sum(int(a['actual_complex_RHS']) for a in ledger if a.get('actual_complex_RHS')),
        unknown_wall_records=[dict(batch=a['batch'],run_id=a['run_id']) for a in ledger if not a.get('wall_seconds')],
        scope='Cohort tables overlap through exact prediction reuse; only COST_LEDGER counts actual executions once. Mixed matrix dimensions and missing attempts prevent interpreting global counter sums as equal FLOPs or complete cost.',
        rejected_category_counts=dict(Counter(a['category'] for a in rejected)),
        gate=dict(nrmse=.15,iou=.7,data_clean=.001,data_noisy=.003),
        nrmse_definition='RMSE / declared per-object amplitude from sealed contract',
        native_affine_status=('FRESH_BLIND_NATIVE_AFFINE_CONFIRMATION_SCORED' if any(c['tag']=='e1_native_confirmation' for c in cohorts) else ('CONTROLLED_DEVELOPMENT_NO_NEW_FRESH_CONFIRMATION' if any('native_tsom' in c['tag'] for c in cohorts) else 'PENDING_CLOSED_SCORING')))
    (results/'RESULTS_SUMMARY.json').write_text(json.dumps(clean(summary),indent=2)+'\n')
    write_csv(results/'RESULTS_BY_CASE.csv',by_case)
    write_csv(results/'REJECTED_CANDIDATES.csv',rejected)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--batch', required=True, choices=['e1', 'e2', 'e3'])
    ap.add_argument('--results', required=True, type=Path, nargs='+')
    ap.add_argument('--private-dir', default=DEFAULT_PRIVATE, type=Path)
    ap.add_argument('--public-dir', default=STUDY/'data/e1_public', type=Path)
    ap.add_argument('--metadata-file', type=Path, help='Explicit canonical evaluation metadata; only opened after closure guards.')
    ap.add_argument('--authorize-closed', required=True, help='Must equal PREDICTIONS_CLOSED; parent authorization required.')
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--tag', help='Output/ledger tag; defaults to batch.')
    ap.add_argument('--curve-k',type=int,nargs='+',help='Only these actual prefix sizes are compared; default all.')
    ap.add_argument('--classification', default='POST_CLOSURE_CONTROLLED_STUDY')
    args = ap.parse_args()
    tag = args.tag or args.batch
    if args.authorize_closed != 'PREDICTIONS_CLOSED': raise SystemExit('Private scoring requires explicit parent closure authorization')
    raw_rows, materials, predictions, candidates, endpoints = [], {}, {}, {}, {}
    source_records = []
    run_ids = set()
    for result_path in args.results:
        closed, rows, mm, pp, cc, ee = load_sealed_results(result_path)
        overlap = run_ids.intersection(r['run_id'] for r in rows)
        if overlap: raise SystemExit(f'Duplicate run IDs across scoring inputs: {sorted(overlap)[:3]}')
        run_ids.update(r['run_id'] for r in rows)
        raw_rows.extend(rows)
        materials.update(mm); predictions.update(pp); candidates.update(cc); endpoints.update(ee)
        source_records.append(dict(path=str(result_path), status=closed['status'],
                                   batch_wall=closed.get('batch_wall', closed.get('elapsed_wall_seconds',closed.get('wall_seconds'))),
                                   workers=closed.get('workers',closed.get('manifest',{}).get('lanes',1)),
                                   setup_wall_seconds=closed.get('setup_wall_seconds'),shared_setup=closed.get('shared_setup')))
    # First private access. Do not move above both closure checks.
    mapping_rows = json.loads((args.metadata_file or (args.private_dir/'mapping.json')).read_text())
    mapping = {r['case_id']: r for r in mapping_rows}
    selected_case_ids = {r['case_id'] for r in raw_rows}
    selected_mapping = [r for r in mapping_rows if r['case_id'] in selected_case_ids]
    truths = {r['object_id']: npz_dict(args.private_dir/r['truth_file']) for r in selected_mapping}
    gate_diagnostics = {}
    diagnostic_path = STUDY/'results/gt_gate_diagnostic_scores.json'
    if diagnostic_path.exists():
        gate_diagnostics = json.loads(diagnostic_path.read_text()).get('objects', {})
    for oid, truth in truths.items():
        if oid not in gate_diagnostics:
            self_mask, threshold = otsu_mask(truth['m'])
            stored = np.asarray(truth['mask'], bool).ravel()
            union = np.count_nonzero(self_mask | stored)
            meta = next(r for r in selected_mapping if r['object_id'] == oid)
            gate_diagnostics[oid] = dict(object_id=oid, side=meta['side'], shape=meta['shape'], amplitude=meta['amplitude'],
                gt_self_otsu_iou=float(np.count_nonzero(self_mask & stored)/union) if union else 1.,
                gt_self_otsu_threshold=threshold, stored_support_rule='Unchanged stored truth mask; E1 and fresh confirmation use occupancy > .35, historical E3 uses its original mask')
    diagnostic_path.write_text(json.dumps(dict(status='OFFLINE_GATE_DIAGNOSTIC_AFTER_AUTHORIZED_CLOSURE',
        physical_calls=0, gates_changed=False, objects=gate_diagnostics), indent=2)+'\n')
    scored, candidate_scores = [], []
    for raw in raw_rows:
        meta = mapping[raw['case_id']]
        if 'noise_label' not in raw and 'noise' in meta:
            condition=meta['noise']
            raw['noise_label']=str(condition) if isinstance(condition,str) else ('clean' if condition==0 else f'noisy_{condition:g}')
            if isinstance(condition,(int,float)): raw['noise_level']=condition
        raw['side']=meta['side'];raw['layout']=meta['layout']
        truth = truths[meta['object_id']]
        row = {k: raw[k] for k in ('run_id', 'case_id', 'side', 'layout', 'start', 'method', 'k', 'budget', 'budget_label', 'noise', 'noise_label', 'wall_budget', 'physical_caps', 'base_method', 'chart', 'affine_mode', 'comparison_role', 'stop', 'worker_wrapper_status', 'worker_wrapper_stop_reason', 'worker_phase_stop_reasons', 'raw_declared_max_lu', 'raw_declared_max_rhs', 'selected_phase', 'current_schedule_status', 'current_schedule_complete', 'current_max_mf_reached', 'auxiliary_current_relative_residual', 'physical_residual_source', 'execution_source', 'current_capacity_records', 'seed_count_scope', 'archived_endpoint_count') if k in raw}
        representation=('NATIVE_AFFINE_DPCS_APCS' if '_native_tsom' in raw['method'] else ('ANCHORED_CURRENT_PROTOTYPE' if raw['method'].startswith('n1_') else ('CAPACITY_MATCHED_CIE_MF2' if raw['method']=='cie_matched' else 'MATURE_BASELINE_OR_PHYSICAL_OPTIMIZER')))
        row.update(representation=representation,object_id=meta['object_id'], shape=meta['shape'], amplitude=meta['amplitude'],
                   group=group_id(raw), comparison_group=group_id({k:v for k,v in raw.items() if k not in ('noise','noise_label','noise_level')}), run_status=raw.get('status'), score_status='MISSING_PREDICTION',
                   relative_residual=raw.get('relative_residual'), online_physical_score=raw.get('score', raw.get('physical_score')), wall_seconds=raw.get('wall_seconds'),
                   model_cost=raw.get('model_cost', {}), actual_corrected_seeds=raw.get('corrected_count'),
                   data_pass=False, image_pass=False, material_nrmse_pass=False, material_data_pass=False, joint_pass=False)
        noise = raw.get('noise_level', raw.get('noise', 0.))
        noisy = str(raw.get('noise_label', '')).lower() not in ('', 'clean', 'none', '0') or (isinstance(noise, (float, int)) and noise > 0)
        gate = .003 if noisy else .001
        row['physical_data_gate'] = gate
        if raw['run_id'] in materials:
            m = materials[raw['run_id']]
            try:
                row.update(material_metrics(m, truth['m'], meta['amplitude'], truth['mask']))
                row['score_status'] = 'SCORED'
                if raw['run_id'] in predictions:
                    obs = npz_dict(args.public_dir/f"{raw['case_id']}.npz")['Y']
                    measured = float(np.linalg.norm(predictions[raw['run_id']] - obs)/np.linalg.norm(obs))
                    row['saved_prediction_residual'] = measured
                    row['residual_record_difference'] = abs(measured - row['relative_residual']) if row['relative_residual'] is not None else None
                    row['relative_residual'] = measured
                residual = row['relative_residual']
                row['data_pass'] = bool(residual is not None and np.isfinite(residual) and residual <= gate)
                row['material_data_pass'] = bool(row['material_nrmse_pass'] and row['data_pass'])
                row['joint_pass'] = bool(row['image_pass'] and row['data_pass'])
                row['joint_failure_cause'] = 'PASS' if row['joint_pass'] else ('IOU_GATE_ONLY' if row['material_data_pass'] else 'MATERIAL_OR_DATA_GATE')
            except (ValueError, TypeError) as exc:
                row.update(score_status='INVALID_MATERIAL', scoring_error=str(exc), image_pass=False)
        pool = []
        pool_materials = []
        for key in raw.get('closed_candidate_material_keys', [k for k in candidates if k.startswith(raw['run_id']+'_candidate')]):
            if key not in candidates: continue
            cm = material_metrics(candidates[key], truth['m'], meta['amplitude'], truth['mask'])
            em = endpoints.get(key, {})
            residual = em.get('relative_residual')
            cm.update(run_id=raw['run_id'], candidate_key=key, relative_residual=residual,
                      data_pass=bool(residual is not None and residual <= gate), provenance=em.get('provenance'))
            cm['joint_pass'] = bool(cm['image_pass'] and cm['data_pass'])
            candidate_scores.append(cm)
            if not any(np.array_equal(candidates[key], old) for old in pool_materials):
                pool.append(cm)
                pool_materials.append(candidates[key])
        # Include selected online endpoint: it can be the retained anchor, outside corrected pool.
        corrected_distinct = len(pool)
        if row['score_status'] == 'SCORED':
            matches = [i for i, old in enumerate(pool_materials) if np.array_equal(materials[row['run_id']], old)]
            if matches:
                # An archived endpoint may omit its paid residual; the selected
                # incumbent's sealed prediction supplies that exact endpoint.
                for i in matches:
                    pool[i].update(relative_residual=row['relative_residual'], data_pass=row['data_pass'], joint_pass=row['joint_pass'])
            else:
                pool.append(row)
        row.update(oracle_pool_size_including_online=len(pool),
                   oracle_corrected_candidate_count=corrected_distinct,
                   oracle_best_image_nrmse=min((r['amplitude_normalized_nrmse'] for r in pool), default=None),
                   oracle_any_image_pass=any(r['image_pass'] for r in pool),
                   oracle_any_joint_pass=any(r['joint_pass'] for r in pool),
                   oracle_any_material_data_pass=any(r.get('material_nrmse_pass',False) and r['data_pass'] for r in pool),
                   oracle_joint_coverage_complete=all(r.get('relative_residual') is not None for r in pool))
        scored.append(row)
    summaries = summarize(scored)
    curves = multistart_curves(scored,args.curve_k) if any(r.get('start',0)>0 for r in scored) else []
    result = dict(status='POST_CLOSURE_SCORED', batch=args.batch, tag=tag, classification=args.classification, result_sources=source_records,
                  private_truth_opened_after_parent_closure=True,
                  private_truth_object_ids_loaded=sorted(truths),
                  offline_gt_gate_diagnostics={oid:gate_diagnostics[oid] for oid in truths},
                  gt_diagnostic_scope='Offline self-Otsu versus stored occupancy support only; no physical calls and no gate change.',
                  gates=dict(amplitude_normalized_nrmse=.15, otsu_iou=.7, clean_relative_data=.001, noisy_relative_data=.003),
                  definitions=dict(nrmse='sqrt(mean((m-m_gt)^2)) / declared per-object amplitude from sealed contract',
                      support='Unchanged stored truth mask (occupancy > .35 for E1/confirmation, original historical E3 mask); predicted foreground by 256-bin histogram Otsu',
                      oracle='Best image NRMSE / any joint pass among AVAILABLE sealed endpoints plus selected online incumbent; GT used only here. Joint-pass counts are known-pass lower bounds if physical residual coverage is incomplete.',
                      visualization='Lowest online physical score among shared starts; all material panels use [.02,20].',
                      independence='Material object IDs define contracts across amplitude/grid, not independent morphology draws. Primary descriptive resampling is by fixed shape family (three in E1); layouts and starts are repeated observations.'),
                  summaries=summaries, rows=scored, candidate_scores=candidate_scores,
                  multistart_K_prefix_curves=curves,
                  multistart_scope='Only actual frozen shared-start prefixes are eligible, with no extrapolation; best-Phi uses sealed online physical score, oracle uses GT only offline.',
                  failed_or_joint_failed_runs=[r['run_id'] for r in scored if not r['joint_pass']])
    if args.batch == 'e3':
        result['raw_worker_stop_counts']={method:dict(Counter(str(r.get('worker_wrapper_stop_reason')) for r in scored if r['method']==method)) for method in sorted({r['method'] for r in scored})}
        result.update(N1_E3='NOT_RUN',N2_E3='NOT_RUN',raster_convention='row0 at top, original raster order; unchanged material arrays')
    if args.batch == 'e1' and 'NOT_FRESH_BLIND' in args.classification:
        result['truth_access_deviation_record'] = 'audit/TRUTH_ACCESS_DEVIATION.md'
    if tag == 'e1_initial':
        result['initial_n1_capacity_failure'] = dict(classification='IMPLEMENTATION_CAPACITY_FAILURE',
            record='audit/N1_INITIAL_REPRESENTABILITY.json',
            inference='Initial N1 does not adequately evaluate normal-direction mechanism; retained, not a theory refutation.')
    output = STUDY/'results'/f'{tag}_scores'
    if args.batch == 'e2' and output.with_suffix('.json').exists():
        prior_record = json.loads(output.with_suffix('.json').read_text())
        if prior_record.get('truth_access_deviation'):
            result['truth_access_deviation'] = prior_record['truth_access_deviation']
    output.with_suffix('.json').write_text(json.dumps(clean(result), indent=2)+'\n')
    write_csv(output.with_suffix('.csv'), scored)
    ledger_path = STUDY/'results/COST_LEDGER.csv'
    previous = []
    if ledger_path.exists():
        with ledger_path.open() as file: previous = [r for r in csv.DictReader(file) if r['batch'] != tag]
    setup_rows = []
    if not any(r['batch'] == 'shared_setup' for r in previous):
        manifest_path = STUDY/'audit/SOURCE_MANIFEST.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            setup_rows.append(dict(batch='shared_setup', run_id='public_data_generation', group='shared setup',
                                   status='RECORDED', wall_seconds=manifest.get('generation_wall'),
                                   LU=None, actual_complex_RHS=None,
                                   setup_scope='Data generation wall includes original Green assembly and observation creation; physical counter totals not in SOURCE_MANIFEST. Shared once, excluded from inference sum.'))
    write_csv(ledger_path, deduplicate_cost_ledger(previous + setup_rows + cost_rows(tag, raw_rows) + recorded_audit_cost_rows() + recorded_external_setup_rows()))
    if not args.no_plots and args.batch in ('e1', 'e3'):
        plots(tag, scored, summaries, mapping, truths, materials, raw_rows)
        if curves: multistart_plot(tag,curves)
        if tag=='e1_native_tsom':native_affine_gain_plot(tag,scored,curves,mapping,truths,materials)
        if tag=='e1_native_confirmation':native_confirmation_gain_plot(result,mapping,truths,materials)
        if args.batch=='e3' and tag=='e3_all': fullwave_publication_plots(scored,mapping,truths,materials)
    write_global_aggregate()
    print(json.dumps(clean(dict(batch=args.batch, runs=len(scored), objects=len({r['object_id'] for r in scored}),
                               summaries=[{k:s[k] for k in ('group','runs','data_passes','joint_passes','median_nrmse','oracle_joint_passes')} for s in summaries],
                               output=str(output.with_suffix('.json')))), indent=2))


if __name__ == '__main__': main()
