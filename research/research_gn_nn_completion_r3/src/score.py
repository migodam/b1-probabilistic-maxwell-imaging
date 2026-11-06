"""Offline GT scoring. This module is never imported by public evaluation.

A size/mtime inventory seals all planned prediction receipts and arrays before
any truth file is opened. No hashes are used. Missing predictions stay NOT_RUN.
"""
import json
from pathlib import Path
import time
import numpy as np
from threadpoolctl import threadpool_limits
from research.research_gn_search_space_r2.src.score_search import otsu, iou


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    from research.research_gn_search_space_r2.src.search import safe
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe(value), indent=2, allow_nan=False)+'\n')


def npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {k:z[k] for k in z.files}


def stamp(path):
    p = Path(path)
    if not p.exists():
        return None
    s = p.stat()
    return dict(path=str(p.resolve()), size=s.st_size, mtime_ns=s.st_mtime_ns)


def seal_predictions(study, manifest_path, runs, *, splits=('val', 'test'), output=None, max_observations=None):
    """Close the supplied finite run x observation plan, preserving failures."""
    study = Path(study).resolve()
    manifest_path = Path(manifest_path).resolve()
    manifest = read_json(manifest_path)
    records = [r for r in manifest['records'] if r['split'] in splits]
    if max_observations is not None:
        records = records[:max_observations]
    receipts = []
    for run in runs:
        for record in records:
            base = study/'evaluation'/run/record['obs_id']
            result_json = base.with_suffix('.json')
            result_npz = base.with_suffix('.npz')
            j, z = stamp(result_json), stamp(result_npz)
            receipts.append(dict(run=run, obs_id=record['obs_id'],
                                 disposition='SEALED' if j and z else 'NOT_RUN', json=j, npz=z))
    output = Path(output) if output else study/'evaluation/PREDICTIONS_SEALED.json'
    snapshot_path=output.with_name(output.stem+'_MANIFEST.json')
    snapshot=dict(manifest,records=records)
    if 'parents' in snapshot:
        parents={r['parent_id']for r in records}
        snapshot['parents']=[p for p in snapshot['parents']if p['parent_id']in parents]
    write_json(snapshot_path,snapshot)
    value = dict(schema='gn-nn-r3-prediction-seal-v1', status='ALL_PLANNED_PREDICTIONS_SEALED',
                 sealed_unix=time.time(), manifest=stamp(snapshot_path), source_manifest=stamp(manifest_path), study=str(study),
                 data_root=str(manifest_path.parent.parent), splits=list(splits), runs=list(runs),
                 max_observations=max_observations, observation_ids=[r['obs_id'] for r in records],
                 truth_read=False, integrity_method='size and mtime inventory; no hashes', receipts=receipts)
    write_json(output, value)
    return output


def remap_path(value, original_study, current_study):
    path=Path(value)
    if not path.is_absolute():
        return Path(current_study)/path
    try:
        return Path(current_study)/path.relative_to(Path(original_study))
    except ValueError:
        return path


def validate_stamp(item, *, resolver=None, replay_frozen=False):
    if item is None:
        raise ValueError('Missing sealed artifact inventory')
    actual=stamp(resolver(item['path']) if resolver else item['path'])
    valid=actual is not None and (actual['size']==item['size'] if replay_frozen else actual==item)
    if not valid:
        raise ValueError('A sealed artifact changed or disappeared: '+str(item))


def material_metrics(m, truth, amplitude, alpha):
    m, truth = np.asarray(m).reshape(-1), np.asarray(truth).reshape(-1)
    mask = np.asarray(alpha).reshape(-1) >= .5
    rmse = float(np.sqrt(np.mean((m-truth)**2)))
    threshold = otsu(m)
    return dict(RMSE=rmse, NRMSE_amplitude=rmse/float(amplitude),
                NRMSE=rmse/float(amplitude), IoU_Otsu=iou(m>threshold, mask),
                otsu_threshold_reconstruction_only=threshold, truth_mask='alpha >= 0.5')


def _candidate_metrics(arrays, receipt, truth):
    amplitude = float(np.asarray(truth['amplitude']))
    rows = []
    for c, m in zip(receipt['candidates'], arrays['candidate_materials']):
        row = dict(c)
        row.update(material_metrics(m, truth['m'], amplitude, truth['alpha']))
        row['joint_success'] = bool(row['physical_data_feasible'] and row['NRMSE']<=.15 and row['IoU_Otsu']>=.7)
        rows.append(row)
    return rows


def _retention(candidates, arrays, truth):
    anchor = arrays['anchor']
    anchor_error = material_metrics(anchor, truth['m'], truth['amplitude'], truth['alpha'])['NRMSE']
    rows = []
    for alpha in (.5, 1.):
        raw = next((c for c in candidates if c['label']==f'raw_{alpha:g}'), None)
        corrected = next((c for c in candidates if c['label']==f'corrected_{alpha:g}'), None)
        if raw is None or corrected is None:
            continue
        p = arrays['candidate_materials'][raw['index']] - anchor
        h = arrays['candidate_materials'][corrected['index']] - anchor
        denom = float(p@p)
        gain = anchor_error - raw['NRMSE']
        rows.append(dict(alpha=alpha, displacement_projection_retention=float(h@p/denom) if denom>1e-30 else None,
                         displacement_norm_retention=float(np.linalg.norm(h)/np.sqrt(denom)) if denom>1e-30 else None,
                         raw_NRMSE=raw['NRMSE'], corrected_NRMSE=corrected['NRMSE'],
                         material_gain_retention=(anchor_error-corrected['NRMSE'])/gain if gain>1e-12 else None,
                         raw_had_positive_material_gain=bool(gain>1e-12)))
    return rows


def current_mechanism_metrics(arrays, truth, v48):
    """Offline current-target errors using the public cached receiver rowspace."""
    j=np.asarray(arrays['j_aux'])
    j0=np.asarray(arrays['anchor_j0'])
    target_j=np.asarray(truth['j'])
    delta=np.asarray(arrays.get('delta_current',j-j0))
    v=np.asarray(v48)
    if v.shape!=(576,48):
        raise ValueError('Cached V48 shape(576,48) required; no scorer SVD fallback')
    def energy(x):return float(np.vdot(x,x).real)
    gt_delta=target_j-j0
    null_target=gt_delta-v@(v.conj().T@gt_delta)
    error=j-target_j
    total,visible=energy(error),energy(v.conj().T@error)
    inc_total,inc_visible=energy(delta),energy(v.conj().T@delta)
    projected_error=energy(delta-null_target)
    target_energy=energy(null_target)
    gt_energy=energy(target_j)
    anchor_error=j0-target_j
    anchor_total,anchor_visible=energy(anchor_error),energy(v.conj().T@anchor_error)
    return dict(anchor_current_relative_L2_error=float(np.sqrt(anchor_total/max(gt_energy,1e-300))),
                anchor_current_error_total_energy=anchor_total,anchor_current_error_visible_energy=anchor_visible,
                anchor_current_error_null_energy=max(0.,anchor_total-anchor_visible),
                auxiliary_current_relative_L2_error=float(np.sqrt(total/max(gt_energy,1e-300))),
                projected_target_current_error_absolute=float(np.sqrt(projected_error)),
                projected_target_current_error_relative=float(np.sqrt(projected_error/max(target_energy,1e-300))),
                projected_target_error_normalized_by_GT_current=float(np.sqrt(projected_error/max(gt_energy,1e-300))),
                projected_target_energy=target_energy,
                current_error_total_energy=total,current_error_visible_energy=visible,
                current_error_null_energy=max(0.,total-visible),current_error_visible_fraction=visible/max(total,1e-300),
                current_increment_total_energy=inc_total,current_increment_visible_energy=inc_visible,
                current_increment_null_energy=max(0.,inc_total-inc_visible),
                current_increment_visible_fraction=inc_visible/max(inc_total,1e-300),
                saved_increment_consistency_L2=float(np.linalg.norm(delta-(j-j0))),
                target_definition='P0(jGT-j0), P0=I-V48V48*, all48 receiver directions retained',
                scope='OFFLINE_GT_CURRENT_DIAGNOSTIC_ONLY_NO_NEW_PHYSICS_OR_SVD')


def oracle_null_current(public, anchor, truth):
    """Offline natural current oracle; optimal current L2, not optimal material readout."""
    from research.research_z_alg_r4_diversity.src.physics import Model
    tick = time.perf_counter()
    model = Model(public)
    with threadpool_limits(limits=1):
        if 'V48' in public:
            v = np.asarray(public['V48'])
            projection_source = 'public cached full receiver rowspace V48'
        else:
            with model.operation('OFFLINE_RECEIVER_SVD', dense_SVD=1):
                _, s, vh = np.linalg.svd(model.S, full_matrices=False)
                if s[-1] <= np.finfo(float).eps*max(model.S.shape)*s[0]:
                    raise ValueError('Oracle requires full numerical receiver rank 48')
                v = vh.conj().T
            projection_source = 'offline paid full receiver SVD, no weak truncation'
        with model.operation('OFFLINE_ORACLE_NULL_PROJECTION', current_projection=1):
            delta = np.asarray(truth['j'])-np.asarray(anchor['j0'])
            j = np.asarray(anchor['j0'])+delta-v@(v.conj().T@delta)
            q = model.E0+model.G@j
            raw_m = np.real(np.sum(q.conj()*j, axis=1))/np.maximum(np.sum(np.abs(q)**2, axis=1),1e-30)
            m = np.clip(raw_m, .02, 20.)
            aux_resid = float(np.linalg.norm(model.S@j-public['Y'])/np.linalg.norm(public['Y']))
            defect = float(np.linalg.norm(j-m[:,None]*q)/max(np.linalg.norm(j),1e-300))
        state = model.forward(m, label='OFFLINE_oracle_null_real_material_readout')
        full_resid = float(np.linalg.norm(state['y']-public['Y'])/np.linalg.norm(public['Y']))
    out = material_metrics(m,truth['m'],truth['amplitude'],truth['alpha'])
    gate = .001 if float(np.asarray(public.get('noise',0.)))==0 else .003
    out.update(physical_data_relative_residual=full_resid, physical_data_feasible=full_resid<=gate,
               auxiliary_data_relative_residual=aux_resid, constitutive_relative_defect=defect,
               full_VIE_equation_defect=state['state_relative_residual'],
               projection_source=projection_source, clipped_cells=int(np.count_nonzero(raw_m!=m)),
               wall_seconds=time.perf_counter()-tick, costcounts=dict(model.counts),
               scope='OFFLINE_GT_NATURAL_CURRENT_ORACLE_ACHIEVABLE_REFERENCE_NOT_MATERIAL_OPTIMUM', truth_access='GT_CURRENT_AND_MATERIAL')
    out['joint_success'] = bool(out['physical_data_feasible'] and out['NRMSE']<=.15 and out['IoU_Otsu']>=.7)
    return out, dict(m_oracle=m, m_oracle_unclipped=raw_m, j_oracle=j)


def score_predictions(seal_path, *, output=None, include_oracle=False, study_override=None, replay_frozen=False):
    seal = read_json(seal_path)
    if seal.get('status') != 'ALL_PLANNED_PREDICTIONS_SEALED':
        raise ValueError('Prediction closure is required before GT scoring')
    original_study=Path(seal['study'])
    current_study=Path(study_override).resolve() if study_override else original_study
    if current_study!=original_study and not replay_frozen:
        raise ValueError('Relocated re-scoring requires explicit replay_frozen=True / --replay-frozen')
    resolver=lambda p:remap_path(p,original_study,current_study) if replay_frozen else Path(p)
    validate_stamp(seal['manifest'],resolver=resolver,replay_frozen=replay_frozen)
    # Validate the complete seal before the first truth open.
    for item in seal['receipts']:
        if item['disposition']=='SEALED':
            validate_stamp(item['json'],resolver=resolver,replay_frozen=replay_frozen); validate_stamp(item['npz'],resolver=resolver,replay_frozen=replay_frozen)
    manifest = read_json(resolver(seal['manifest']['path']))
    records = {r['obs_id']:r for r in manifest['records']}
    data_root, study = resolver(seal['data_root']), current_study
    rows, oracles = [], {}
    operators = None
    for item in seal['receipts']:
        record = records[item['obs_id']]
        row = {**record, 'run':item['run'], 'transport_status':item['disposition']}
        if item['disposition']!='SEALED':
            row.update(status='NOT_RUN', joint_success=False)
            rows.append(row)
            continue
        receipt = read_json(resolver(item['json']['path']))
        arrays = npz(resolver(item['npz']['path']))
        truth = npz(data_root/record['truth_path'])
        candidates = _candidate_metrics(arrays, receipt, truth)
        selected = receipt['top1']
        old = receipt.get('old_score_top1')
        row.update(status=receipt['status'], selection_status=receipt['selection_status'],
                   accounting_status=receipt['accounting_status'], wall_seconds=receipt['wall_seconds'],
                   cost_scope=receipt.get('cost_scope','INCREMENTAL_WARM_MODEL_POST_GN'),budget=receipt.get('budget',dict(wall_seconds=10.,max_LU=100,max_complex_RHS=6000)),
                   costcounts=receipt['costcounts'], auxiliary=receipt.get('auxiliary'),
                   anchor_gradient_diagnostic=receipt.get('anchor_gradient_diagnostic',dict(status='NOT_RECORDED',reason='Historical receipt has no paid-J anchor diagnostic')),
                   candidates=candidates, correction_retention=_retention(candidates,arrays,truth),
                   anchor=material_metrics(arrays['anchor'],truth['m'],truth['amplitude'],truth['alpha']),
                   top1=None if selected is None else candidates[selected['index']],
                   old_score_top1=None if old is None else candidates[old['index']])
        if 'j_aux' in arrays and 'anchor_j0' in arrays:
            if operators is None:
                operators=npz(data_root/manifest.get('operators_path','operators.npz'))
            if 'V48' not in operators:
                raise ValueError('Current mechanism scoring requires existing public V48')
            row['current_metrics']=current_mechanism_metrics(arrays,truth,operators['V48'])
            raw_full=next((c for c in candidates if c['label']=='raw_1'),None)
            row['current_metrics'].update(raw_full_material_NRMSE=None if raw_full is None else raw_full['NRMSE'],
                raw_full_physical_data_relative_residual=None if raw_full is None else raw_full['relative_residual'],
                selected_material_NRMSE=None if row['top1']is None else row['top1']['NRMSE'],
                constitutive_relative_defect=(receipt.get('auxiliary')or{}).get('constitutive_relative_defect'),
                auxiliary_data_relative_residual=(receipt.get('auxiliary')or{}).get('auxiliary_data_relative_residual'))
        materials=arrays['candidate_materials']
        groups=[]
        for index,m in enumerate(materials):
            match=next((g for g in groups if np.sqrt(np.mean((m-materials[g[0]])**2))<=1e-8),None)
            if match is None:groups.append([index])
            else:match.append(index)
        row.update(candidate_distinct_material_count=len(groups),candidate_distinct_material_RMS_tolerance=1e-8,
                   candidate_duplicate_slot_groups=groups,candidate_distinctness_scope='RMS-near-equal saved materials, not certified distinct nonlinear roots')
        all_best = min(candidates,key=lambda c:c['NRMSE']) if candidates else None
        qualified = [c for c in candidates if c['physical_data_feasible']]
        qualified_best = min(qualified,key=lambda c:c['NRMSE']) if qualified else None
        row.update(all_bestK_OFFLINE=all_best, qualified_bestK_OFFLINE=qualified_best,
                   all_bestK_count=len(candidates), qualified_bestK_count=len(qualified),
                   bestK_scope='finite five-point candidate set, not all possible nonlinear roots')
        row['joint_success'] = bool(selected is not None and candidates[selected['index']]['joint_success'] and receipt['accounting_status']=='COMPLETE')
        row['accounted_within_budget'] = bool(receipt['accounting_status']=='COMPLETE' and receipt['wall_seconds']<=row['budget']['wall_seconds'] and receipt['costcounts'].get('LU',0)<=row['budget']['max_LU'] and receipt['costcounts'].get('actual_complex_RHS',0)<=row['budget']['max_complex_RHS'])
        if include_oracle and record['obs_id'] not in oracles:
            if operators is None:
                operators=npz(data_root/manifest.get('operators_path','operators.npz'))
            public={**operators,**npz(data_root/record['public_path'])}
            anchor=npz(study/'anchors'/f'{record["obs_id"]}.npz')
            metric, oracle_arrays=oracle_null_current(public,anchor,truth)
            oracle_dir=study/('evaluation/offline_oracle_REPLAY'if replay_frozen else'evaluation/offline_oracle');oracle_dir.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(oracle_dir/f'{record["obs_id"]}.npz',**oracle_arrays)
            oracles[record['obs_id']]=dict(record,**metric)
        rows.append(row)
    groups={}
    for row in rows:
        key=(row['run'],row['split'],row.get('distribution','ID'),str(row['noise']),row.get('level','low'))
        groups.setdefault(key,[]).append(row)
    summary=[]
    for key,group in groups.items():
        valid=[r for r in group if r.get('top1') is not None]
        out=dict(zip(('run','split','distribution','noise','level'),key))
        out.update(n_planned=len(group),n_selected=len(valid),n_joint_success=sum(r['joint_success']for r in group),
                   joint_success_rate=sum(r['joint_success']for r in group)/len(group),
                   n_missing=sum(r['status']=='NOT_RUN'for r in group),
                   n_accounted_within_budget=sum(r.get('accounted_within_budget',False)for r in group))
        for metric in ('NRMSE','IoU_Otsu'):
            values=[r['top1'][metric]for r in valid]
            out[metric+'_mean']=float(np.mean(values))if values else None
            out[metric+'_median']=float(np.median(values))if values else None
        improvements=[r['anchor']['NRMSE']-r['top1']['NRMSE']for r in valid]
        out['NRMSE_improvement_vs_anchor_mean']=float(np.mean(improvements))if improvements else None
        summary.append(out)
    value=dict(schema='gn-nn-r3-scores-v1',seal_path=str(Path(seal_path).resolve()),
               truth_access='HISTORICAL_ALREADY_DISCLOSED_FROZEN_RESCORING'if replay_frozen else'OFFLINE_AFTER_COMPLETE_PREDICTION_CLOSURE',rows=rows,summary=summary,
               replay_frozen=bool(replay_frozen),original_study=str(original_study),current_study=str(current_study),
               replay_inventory_policy='existence and recorded file size only; unpacked mtimes are not blind-seal evidence'if replay_frozen else'original in-session size/mtime inventory checked',
               oracle_rows=list(oracles.values()),
               success_definition='NRMSE/amplitude <= .15 and IoU_Otsu >= .7 and clean data residual <= .001 / noisy <= .003; complete incremental cost receipt required',
               budget_overrun_policy='overruns and accounting status reported independently; n_accounted_within_budget is explicit',
               mask_definition='truth alpha >= .5', all_missing_preserved=True)
    output=Path(output)if output else Path(seal_path).parent/('HISTORICAL_REPLAY_SCORES.json'if replay_frozen else'SCORES.json')
    write_json(output,value)
    return output


def validation_gate(scores_path, *, output=None):
    """Summarize paired VAL raw/top1 evidence; never authorize expansion."""
    scores_path=Path(scores_path)
    scores=read_json(scores_path)
    if any(r['split']!='val' for r in scores['rows']):
        raise ValueError('Validation gate only accepts a VAL-only score artifact')
    groups={}
    for row in scores['rows']:
        groups.setdefault(row['run'],[]).append(row)
    summaries=[]
    for run,rows in groups.items():
        valid=[r for r in rows if r.get('top1')is not None]
        raw=[(r,next((c for c in r.get('candidates',[])if c['label']=='raw_1'),None))for r in rows]
        raw=[(r,c)for r,c in raw if c is not None]
        baseline=[r for r in rows if r.get('anchor')]
        base_mean=float(np.mean([r['anchor']['NRMSE']for r in baseline]))if baseline else None
        base_qualified=sum(any(c['kind']=='anchor'and c['physical_data_feasible']for c in r.get('candidates',[]))for r in rows)/len(rows)
        item=dict(run=run,n_planned=len(rows),n_scored=len(baseline),n_selected=len(valid),n_raw_full=len(raw),
                  anchor_mean_NRMSE_amplitude=base_mean,anchor_physical_qualification_rate=base_qualified)
        for label,pairs in [('raw_full',raw),('top1',[(r,r['top1'])for r in valid])]:
            means=[c['NRMSE']for _,c in pairs]
            paired_base=float(np.mean([r['anchor']['NRMSE']for r,_ in pairs]))if pairs else None
            mean=float(np.mean(means))if means else None
            qualified=None if label=='raw_full' and not pairs else sum(c['physical_data_feasible']for _,c in pairs)/len(rows)
            absolute=None if mean is None else paired_base-mean
            relative=None if mean is None or paired_base<=0 else absolute/paired_base
            item[label]=dict(mean_NRMSE_amplitude=mean,paired_anchor_mean_NRMSE_amplitude=paired_base,
                             absolute_NRMSE_improvement=absolute,relative_mean_NRMSE_improvement=relative,
                             physical_qualification_rate=qualified,qualification_drop_pp=None if qualified is None else 100*(base_qualified-qualified),
                             n_joint_success=sum(c['joint_success']for _,c in pairs),
                             pass_relative5pct=None if relative is None else relative>=.05,
                             pass_absolute_point05=None if absolute is None else absolute>=.05,
                             pass_qualification_drop5pp=None if qualified is None else base_qualified-qualified<=.05+1e-12)
        measured=[r for r in rows if 'wall_seconds'in r]
        walls=[r['wall_seconds']for r in measured]
        item['costs']=dict(wall_seconds_total=float(sum(walls)),wall_seconds_mean=float(np.mean(walls))if walls else None,
                           wall_seconds_max=max(walls)if walls else None,
                           LU=sum(r['costcounts'].get('LU',0)for r in measured),
                           complex_RHS=sum(r['costcounts'].get('actual_complex_RHS',0)for r in measured),
                           n_accounted_within_budget=sum(r.get('accounted_within_budget',False)for r in rows))
        main=item['top1']
        item['relative_gate_candidate']=bool(len(valid)==len(rows)and all(r.get('accounting_status')=='COMPLETE'for r in rows)
            and main['pass_relative5pct']and main['pass_qualification_drop5pp'])
        item['absolute_gate_candidate']=bool(len(valid)==len(rows)and all(r.get('accounting_status')=='COMPLETE'for r in rows)
            and main['pass_absolute_point05']and main['pass_qualification_drop5pp'])
        summaries.append(item)
    value=dict(schema='gn-nn-r3-validation-gate-v1',scores_path=str(scores_path.resolve()),rows=summaries,
               source_split='VAL_ONLY',relative_gate='(mean anchor amplitude-NRMSE minus mean corrected-top1 amplitude-NRMSE) / mean anchor >= .05, with physical qualification drop <= .05',
               absolute_gate_also_reported='absolute amplitude-NRMSE decrease >= .05, retained for unambiguous owner interpretation',
               decision='PARENT_SCIENTIFIC_REVIEW_REQUIRED_NO_AUTOMATIC_EXPANSION')
    output=Path(output)if output else scores_path.with_name(scores_path.stem+'_GATE.json')
    write_json(output,value)
    table=output.with_suffix('.md')
    lines=['# Paired validation evidence','',
           '| Run | N | Anchor NRMSE | Raw NRMSE | Corrected top1 NRMSE | Relative gain | Raw physical | Top1 physical | Qual. drop | Mean/max sec | LU / RHS |',
           '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    def f(x,d=5):return 'NA'if x is None else f'{x:.{d}f}'
    for row in summaries:
        raw,top,cost=row['raw_full'],row['top1'],row['costs']
        gain=top['relative_mean_NRMSE_improvement']
        lines.append(f'| {row["run"]} | {row["n_selected"]}/{row["n_planned"]} | {f(row["anchor_mean_NRMSE_amplitude"])} | {f(raw["mean_NRMSE_amplitude"])} | {f(top["mean_NRMSE_amplitude"])} | {f(None if gain is None else 100*gain,2)}% | {f(None if raw["physical_qualification_rate"]is None else 100*raw["physical_qualification_rate"],1)}% | {100*top["physical_qualification_rate"]:.1f}% | {top["qualification_drop_pp"]:.1f}pp | {f(cost["wall_seconds_mean"],3)}/{f(cost["wall_seconds_max"],3)} | {cost["LU"]} / {cost["complex_RHS"]} |')
    lines+=['','Relative gate: corrected top1 must improve paired mean amplitude-normalized NRMSE by at least5%, with physical qualification dropping by at most5 percentage points. JSON also reports the absolute .05 criterion. Raw and corrected qualification rates use the same clean .001 / noisy .003 full-VIE gates.','',
            'This table is validation evidence for the parent scientific decision. It does not authorize dataset expansion. Missing predictions, failures, incomplete cost receipts, and budget overruns remain explicit in the source scores.']
    table.write_text('\n'.join(lines)+'\n')
    return output,table
