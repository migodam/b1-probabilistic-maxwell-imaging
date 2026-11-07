"""R4 offline scoring, after finite public predictions are sealed."""
from pathlib import Path
import csv,json,time
import numpy as np
from research.research_gn_nn_completion_r3.src.score import (
    seal_predictions,score_predictions,read_json,write_json,npz,
    material_metrics,
)

def structure_metrics(m,truth,amplitude):
    x=np.asarray(m).reshape(24,24);t=np.asarray(truth).reshape(24,24)
    error=x-t
    edge=np.concatenate((np.diff(error,axis=0).ravel(),np.diff(error,axis=1).ravel()))
    return {'MAE':float(np.mean(abs(error))),
            'edge_RMSE_amplitude':float(np.sqrt(np.mean(edge**2))/amplitude),
            'continuous_relative_L2':float(np.linalg.norm(error)/np.linalg.norm(t))}

def csv_write(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v)if isinstance(v,(dict,list))else v for k,v in r.items()})

def score(seal_path,output):
    # The inherited scorer enforces the recorded seal before first GT arrays.
    protocol_seal=read_json(seal_path)
    if 'test' in protocol_seal['splits'] and not (Path(protocol_seal['study'])/'audit/METHOD_FREEZE.json').exists():
        raise RuntimeError('New TEST numerical scoring requires the final method freeze')
    path=score_predictions(seal_path,output=output,include_oracle=False)
    data=read_json(path);seal=read_json(seal_path);root=Path(seal['data_root']);study=Path(seal['study'])
    truth_cache={};truth_calibration={};anchor_cost_cache={};tables=[]
    for row in data['rows']:
        if row.get('transport_status')!='SEALED':
            table={k:row.get(k)for k in ('run','obs_id','parent_id','split','group','noise','family','digit','amplitude','seen_image','seen_parent','seen_class','seen_font','seen_rule')}
            table.update(status='NOT_RUN',joint_success=False,budgeted_joint_success=False,
                         accounted_within_budget=False,wall_seconds=None,anchor_NRMSE=None,
                         selected_NRMSE=None,selected_data_relative_residual=None,selected_IoU=None,
                         selected_edge_RMSE=None,LU=None,RHS=None)
            table.update(online_total_wall_seconds=None,online_total_LU=None,online_total_RHS=None,online_total_within55=False)
            row['budgeted_joint_success']=False;tables.append(table)
            continue
        key=row['truth_path']
        if key not in truth_cache:
            truth_cache[key]=npz(root/key)
            t=truth_cache[key];g=material_metrics(t['m'],t['m'],t['amplitude'],t['alpha'])
            truth_calibration[key]={'parent_id':row['parent_id'],'group':row.get('group'),'family':row.get('family'),
                                    'GT_self_IoU_Otsu':g['IoU_Otsu'],'GT_self_shape_gate':g['IoU_Otsu']>=.7}
        truth=truth_cache[key];a=npz(study/'evaluation'/row['run']/f"{row['obs_id']}.npz")
        amp=float(truth['amplitude'])
        row['actual_saved_candidate_count']=len(row['candidates'])
        row['bestK_scope']='OFFLINE over the actual finite saved candidate set; K is method-specific, not universally five; no posterior or global-root coverage'
        for candidate,material in zip(row['candidates'],a['candidate_materials']):
            candidate.update(structure_metrics(material,truth['m'],amp))
        # Copies of selected candidate dictionaries are from the same JSON model.
        for label in ('top1','old_score_top1','qualified_bestK_OFFLINE','all_bestK_OFFLINE'):
            c=row.get(label)
            if c is not None:c.update(structure_metrics(a['candidate_materials'][c['index']],truth['m'],amp))
        row['anchor'].update(structure_metrics(a['anchor'],truth['m'],amp))
        row['scene_identity']={k:row.get(k)for k in ('parent_id','split','group','family','digit','font','generator','seen_image','seen_parent','seen_class','seen_font','seen_rule')}
        row['GT_used_online']=False
        row['budgeted_joint_success']=bool(row['joint_success'] and row.get('accounted_within_budget',False))
        cold=row['cost_scope'].startswith(('INDEPENDENT','COLD','SEPARATE_FULL45'))
        own_receipt=read_json(study/'evaluation'/row['run']/f"{row['obs_id']}.json")
        # Old anchor preparation computes an NN feature even for ordinary
        # baselines. Keep it in execution ledgers, but do not charge an unused
        # feature as a necessary traditional deployment operation.
        unused_features=float(own_receipt.get('baseline_unused_feature_wall_seconds',0.))if cold else 0.
        row['actual_executed_method_wall_seconds']=row['wall_seconds']
        row['unused_NN_feature_wall_seconds_excluded_from_deployment']=unused_features
        if row['obs_id']not in anchor_cost_cache:
            anchor_cost_cache[row['obs_id']]=read_json(study/'anchors'/f"{row['obs_id']}.json")
        base=anchor_cost_cache[row['obs_id']]
        row['online_total_wall_seconds']=row['wall_seconds']-unused_features+(0. if cold else base['solve_wall_seconds'])
        row['online_total_LU']=row['costcounts'].get('LU',0)+(0 if cold else base['counts'].get('LU',0))
        row['online_total_RHS']=row['costcounts'].get('actual_complex_RHS',0)+(0 if cold else base['counts'].get('actual_complex_RHS',0))
        row['online_total_within55']=bool(row['accounting_status']=='COMPLETE'and row['online_total_wall_seconds']<=55 and row['online_total_LU']<=600 and row['online_total_RHS']<=36000)
        row['joint_success_total55']=bool(row['joint_success']and row['online_total_within55'])
        table={k:row.get(k)for k in ('run','obs_id','parent_id','split','group','noise','family','digit','amplitude','seen_image','seen_parent','seen_class','seen_font','seen_rule','status','joint_success','accounted_within_budget','wall_seconds')}
        table.update(anchor_NRMSE=row['anchor']['NRMSE'],
                     budgeted_joint_success=row['budgeted_joint_success'],
                     selected_NRMSE=None if row.get('top1')is None else row['top1']['NRMSE'],
                     selected_data_relative_residual=None if row.get('top1')is None else row['top1']['relative_residual'],
                     selected_IoU=None if row.get('top1')is None else row['top1']['IoU_Otsu'],
                     selected_edge_RMSE=None if row.get('top1')is None else row['top1']['edge_RMSE_amplitude'],
                     LU=row.get('costcounts',{}).get('LU'),RHS=row.get('costcounts',{}).get('actual_complex_RHS'))
        table.update({k:row[k]for k in ('online_total_wall_seconds','online_total_LU','online_total_RHS','online_total_within55','joint_success_total55')})
        table.update(actual_executed_method_wall_seconds=row['actual_executed_method_wall_seconds'],
                     unused_NN_feature_wall_seconds_excluded_from_deployment=unused_features)
        raw=next((c for c in row['candidates']if c['label']=='raw_1'),None)
        table.update(raw_NRMSE=None if raw is None else raw['NRMSE'],
                     raw_data_relative_residual=None if raw is None else raw['relative_residual'],
                     selected_label=None if row.get('top1')is None else row['top1']['label'],
                     anchor_fallback=bool(row.get('top1')and row['top1']['kind']=='anchor'),
                     selected_RMS_displacement=row.get('selected_RMS_displacement_from_anchor'),
                     old_score_NRMSE=None if row.get('old_score_top1')is None else row['old_score_top1']['NRMSE'],
                     distinct_material_candidates=row.get('candidate_distinct_material_count'),
                     qualified_bestK_NRMSE=None if row.get('qualified_bestK_OFFLINE')is None else row['qualified_bestK_OFFLINE']['NRMSE'],
                     all_bestK_NRMSE=None if row.get('all_bestK_OFFLINE')is None else row['all_bestK_OFFLINE']['NRMSE'])
        table['actual_saved_candidate_count']=row['actual_saved_candidate_count']
        tables.append(table)
    # Independent cold CIE has its own uniform starting material. Comparative
    # gains and failure-filled diagnostics use the shared strict GN, not that
    # easier cold reference. Endpoint quality itself remains unchanged.
    common={}
    ordered=sorted(data['rows'],key=lambda r:0 if r['run']=='strict_gn45_anchor'else 1)
    for row in ordered:
        if row.get('transport_status')!='SEALED'or row['run'].startswith('cie_'):continue
        candidate=next((c for c in row.get('candidates',[])if c['label']=='anchor'),None)
        if candidate is not None:
            common.setdefault(row['obs_id'],dict(row['anchor'],joint_success=bool(candidate['joint_success']),source_run=row['run']))
    by_table={(r['run'],r['obs_id']):r for r in tables}
    for row in data['rows']:
        reference=common.get(row['obs_id']);table=by_table[row['run'],row['obs_id']]
        table['method_anchor_NRMSE']=table['anchor_NRMSE']
        row['common_strict_GN_reference']=reference
        if reference is not None:
            table.update(anchor_NRMSE=reference['NRMSE'],common_strict_anchor_NRMSE=reference['NRMSE'],
                         common_strict_anchor_joint_success=reference['joint_success'])
            gained=bool(row.get('joint_success')and not reference['joint_success'])
            lost=bool(reference['joint_success']and not row.get('joint_success'))
            row.update(new_joint_success_over_strict_GN=gained,lost_joint_success_vs_strict_GN=lost)
            table.update(new_joint_success_over_strict_GN=gained,lost_joint_success_vs_strict_GN=lost)
    data.update(schema='gn-nn-r4-scored-v1',scored_unix=time.time(),
                truth_scope='VAL_MODEL_SELECTION' if set(seal['splits'])=={'val'} else 'NEW_TEST_AFTER_METHOD_AND_PREDICTION_FREEZE',
                statistical_unit='independent parent; clean/noisy and model seeds are correlated repetitions',
                current_energy_scope='src.score.current_metrics are explicitly auxiliary-head diagnostics; true physical current energies are computed separately by experiments/test_energy.py, whose runtime status must be reported')
    csv_write(Path(path).with_suffix('.csv'),tables)
    groups={}
    for row in tables:
        for grouping in ('ALL',row.get('group','unknown'),row.get('family','unknown')):
            groups.setdefault((row['run'],row['split'],row['noise'],grouping),[]).append(row)
    summaries=[]
    for (run,split,noise,group),rs in groups.items():
        valid=[r for r in rs if r['selected_NRMSE']is not None]
        summaries.append({'run':run,'split':split,'noise':noise,'grouping':group,'planned':len(rs),
                          'parents':len({r['parent_id']for r in rs}),'selected':len(valid),
                          'joint_success':sum(bool(r['joint_success'])for r in rs),
                          'budgeted_joint_success':sum(bool(r.get('budgeted_joint_success'))for r in rs),
                          'joint_success_total55':sum(bool(r.get('joint_success_total55'))for r in rs),
                          'anchor_fallback':sum(bool(r.get('anchor_fallback'))for r in rs),
                          'new_joint_success_over_strict_GN':sum(bool(r.get('new_joint_success_over_strict_GN'))for r in rs),
                          'lost_joint_success_vs_strict_GN':sum(bool(r.get('lost_joint_success_vs_strict_GN'))for r in rs),
                          'not_run':sum(r['status']=='NOT_RUN'for r in rs),
                          'NRMSE_mean':float(np.mean([r['selected_NRMSE']for r in valid]))if valid else None,
                          'anchor_NRMSE_mean':float(np.mean([r['anchor_NRMSE']for r in rs if r['anchor_NRMSE']is not None]))if any(r['anchor_NRMSE']is not None for r in rs)else None,
                          'IoU_mean':float(np.mean([r['selected_IoU']for r in valid]))if valid else None,
                          'wall_mean':float(np.mean([r['wall_seconds']for r in rs if r['wall_seconds']is not None]))if any(r['wall_seconds']is not None for r in rs)else None,
                          'LU_total':sum(r['LU']or 0 for r in rs),'RHS_total':sum(r['RHS']or 0 for r in rs)})
    csv_write(Path(path).with_name(Path(path).stem+'_SUMMARY.csv'),summaries)
    data['inherited_R3_summary']=data.pop('summary')
    data['summary']=summaries
    data['summary_scope']='R4 groups explicitly distinguish new_instance, unseen_class and unseen_style; inherited R3 ID defaults are not used for these claims'
    data['GT_shape_metric_calibration']={'parents':list(truth_calibration.values()),
        'warning':'Otsu on continuous GT can differ from alpha>=.5; this describes the frozen metric and is not a strict upper bound on reconstruction IoU. Gate unchanged.'}
    write_json(path,data)
    return path
