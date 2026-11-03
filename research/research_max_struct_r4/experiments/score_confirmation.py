"""Independent confirmation sealing/scoring. No algorithm or physics imports.

Preparation is not execution. `seal` uses public data only. Private truth is
opened solely by `score --authorize-confirmation-truth` after a complete seal.
Outputs are write-once. Missing tasks require explicit finalized dispositions.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import json
import math
import time
import numpy as np
from offline_score_common import (read,write_new,stamp,unchanged,public_metrics,
    prediction,physical,metrics,sum_known,write_csv_new,frozen_contract)

ROOT=Path(__file__).resolve().parents[1]
METHODS=('onehop','global','native_extended')
FINAL_OUTER={'OK','FAIL','FAILED','HARD_TIMEOUT','WORKER_FAILED','UNKNOWN','CANCELLED'}
FINAL_INNER={'SEALED_PREDICTION','FAILED','UNKNOWN','HARD_TIMEOUT','WORKER_FAILED','CANCELLED'}
PHYSICAL_KEYS=('LU','source_RHS','adjoint_RHS','actual_complex_RHS')


def utc():return datetime.now(timezone.utc).isoformat()


def specifications():
    public=ROOT/'data_public/confirmation'
    manifest=read(public/'manifest.json')
    cfg=read(ROOT/'configs/confirmation_frozen.json')
    execution=read(ROOT/'audit/confirmation/EXECUTION_MANIFEST.json')
    ids=[x['observation_id'] for x in manifest['records']]
    if len(ids)!=6 or len(set(ids))!=6 or cfg['primary_methods']!=list(METHODS):
        raise ValueError('Frozen six-object/three-method manifest mismatch')
    expected={ident+'__'+method for ident in ids for method in METHODS}
    tasks={x['task_id']:x for x in execution['tasks']}
    if set(tasks)!=expected or len(execution['tasks'])!=18:
        raise ValueError('Expected all18 frozen tasks exactly once')
    if any(x['backend']!='cpu' for x in tasks.values()):
        raise ValueError('Frozen confirmation is CPU; do not report legacy cuda_workers metadata as actual use')
    if not cfg.get('secondary_selection'):
        raise ValueError('Predeclared secondary Phi-selection contract missing')
    return public,manifest,cfg,tasks


def evidence_record(task,remote,explicit):
    ident=task['task_id']; obs,method=ident.split('__',1)
    receipt_path=remote/'results/confirmation'/(ident+'.json')
    sources=[];override=explicit.get(ident)
    outer=read(receipt_path) if receipt_path.exists() else {}
    if receipt_path.exists():sources.append(stamp(receipt_path))
    if not outer and override is None:
        raise RuntimeError('Task has no terminal receipt or explicit finalized disposition: '+ident)
    if override is not None and override.get('closed') is not True:
        raise RuntimeError('Explicit disposition must set closed=true: '+ident)
    outer_status=outer.get('status',override.get('status') if override else None)
    if outer_status not in FINAL_OUTER and override is None:
        raise RuntimeError('Outer task not closed: '+ident+' '+str(outer_status))
    inner=outer.get('result') or {}
    inner_status=inner.get('status','UNKNOWN')
    if inner_status not in FINAL_INNER and override is None:
        raise RuntimeError('Inner scientific disposition not closed: '+ident+' '+str(inner_status))
    status='UNKNOWN';reason=None
    if override:
        status=override['status'];reason=override.get('reason','Explicit finalized disposition')
        if status not in FINAL_OUTER|FINAL_INNER:
            raise ValueError('Unrecognized explicit terminal status: '+status)
    elif outer_status=='OK' and inner_status=='SEALED_PREDICTION':
        status='SEALED_PREDICTION'
    else:
        status=inner_status if inner_status!='UNKNOWN' else outer_status
        if status=='OK':status='UNKNOWN'
        reason=inner.get('error',outer.get('error','No sealed scientific prediction'))
    if outer.get('task',{}).get('task_id',ident)!=ident or inner.get('method',method)!=method:
        status='UNKNOWN';reason='Receipt task or method identity mismatch'
    prediction_path=remote/'predictions'/obs/(method+'.npz')
    if prediction_path.exists():sources.append(stamp(prediction_path))
    phi=None;validation_error=None
    if status=='SEALED_PREDICTION':
        try:
            m,yp=prediction(prediction_path,(24,24))
            public=public_metrics(ROOT/'data_public/confirmation'/(obs+'.npz'))
            phi=physical(m,yp,public,(24,24))['physical_objective']
            if not math.isfinite(phi):raise ValueError('Nonfinite original Phi')
        except Exception as exc:
            status='UNKNOWN';validation_error=repr(exc);reason='Sealed prediction unavailable/invalid'
    return dict(task_id=ident,observation_id=obs,method=method,outer_status=outer_status,
        inner_status=inner_status,disposition=status,closed=True,reason=reason,
        prediction_validation_error=validation_error,original_Phi_pretruth=phi,
        receipt=str(receipt_path),prediction=str(prediction_path),evidence=sources,
        explicit_disposition=override,outer=outer,
        frozen_configuration_matches=(inner.get('configuration')==task['config']) if inner.get('configuration') is not None else None)


def seal(args):
    start,cpu=time.perf_counter(),time.process_time()
    public,manifest,cfg,tasks=specifications()
    inputs=[stamp(public/'manifest.json'),stamp(ROOT/'configs/confirmation_frozen.json'),
            stamp(ROOT/'audit/confirmation/EXECUTION_MANIFEST.json'),stamp(Path(__file__)),
            stamp(Path(__file__).with_name('offline_score_common.py'))]
    explicit={}
    if args.dispositions:
        dispositions=read(args.dispositions)
        records=dispositions['records']
        explicit={x['task_id']:x for x in records}
        if len(explicit)!=len(records) or not set(explicit)<=set(tasks):
            raise ValueError('Duplicate or out-of-manifest explicit disposition')
        inputs.append(stamp(args.dispositions))
    entries=[evidence_record(tasks[key],args.remote_root,explicit) for key in sorted(tasks)]
    # No truth access above or below in this operation.
    secondary=[]
    for record in manifest['records']:
        ident=record['observation_id']
        inputs.append(stamp(public/record['path']))
        group=[x for x in entries if x['observation_id']==ident]
        available=[x for x in group if x['disposition']=='SEALED_PREDICTION']
        minimum=min((x['original_Phi_pretruth'] for x in available),default=None)
        secondary.append(dict(observation_id=ident,complete_three_arm_pool=len(available)==3,
            minimum_available_original_Phi=minimum,
            retained_available_methods=[x['method'] for x in available if x['original_Phi_pretruth']==minimum],
            selection_status='SEALED_COMPLETE_POOL' if len(available)==3 else 'PARTIAL_POOL_UNKNOWN_MINIMUM',
            tie_policy='strict original Phi minimum; retain all exact ties; no added numerical tolerance or GT tiebreak'))
    for path in (args.shared_setup,args.remote_root/'results/confirmation/runner_summary.json',
                 args.remote_root/'results/confirmation/worker_affinity.json',
                 ROOT/'audit/confirmation/REMOTE_RUN_RECEIPT.json',
                 ROOT/'audit/confirmation/WINDOWS_ZIP_PATH_NORMALIZATION.json'):
        if path.exists():inputs.append(stamp(path))
    result=dict(status='ALL18_DISPOSITIONS_SEALED',sealed_utc=utc(),expected_tasks=18,
        expected_objects=6,entries=entries,secondary_pretruth_selection=secondary,inputs=inputs,
        shared_setup=str(args.shared_setup.resolve()),remote_root=str(args.remote_root.resolve()),
        frozen_secondary_selection=cfg['secondary_selection'],
        denominator_policy='All6 objects for each method and secondary pool; FAILED/UNKNOWN retained, never discarded',
        truth_opened=False,hash_checks=False,new_physical_LU=0,new_complex_RHS=0,
        wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)
    write_new(args.seal,result)
    print(json.dumps(dict(status=result['status'],seal=str(args.seal),truth_opened=False,
                         dispositions={s:sum(x['disposition']==s for x in entries) for s in sorted({x['disposition'] for x in entries})})))


def task_cost(entry):
    outer=entry['outer'];inner=outer.get('result') or {}
    physical_cost=inner.get('physical_model_cost')
    complete=isinstance(physical_cost,dict) and entry['outer_status']=='OK'
    # Top-level total only. Never add setup/auxiliary/corrector nested costs.
    row=dict(observation_id=entry['observation_id'],method=entry['method'],
        disposition=entry['disposition'],physical_counts=physical_cost,
        physical_cost_complete=complete,task_worker_wall_seconds=outer.get('total_worker_s',outer.get('wall_observed_s')),
        adapter_wall_seconds=inner.get('wall_seconds'),adapter_cpu_seconds=inner.get('cpu_seconds'),
        worker_slot=outer.get('slot'),worker_pid=outer.get('pid'),
        solver_import_charge=inner.get('solver_import_charge'),memory=inner.get('memory'),
        auxiliary_operation_counts=inner.get('auxiliary_operation_counts'),
        failed_attempts=inner.get('failed_attempts'),outer_cost_status=outer.get('costs_status'),
        bootstrap_import_cost=None,actual_GPU_tasks=0,source_count=12)
    for key in PHYSICAL_KEYS:
        if not complete:
            row[key]=physical_cost.get(key) if isinstance(physical_cost,dict) else None
        else:
            row[key]=physical_cost.get(key,0)
    return row


def cost_summary(rows):
    fields=PHYSICAL_KEYS+('task_worker_wall_seconds','adapter_wall_seconds','adapter_cpu_seconds')
    out={key:sum_known(row.get(key) for row in rows) for key in fields}
    out['physical_counter_scope']='sum top-level physical_model_cost only; nested phase costs are diagnostic partitions'
    out['physical_cost_complete']=all(row['physical_cost_complete'] for row in rows)
    out['actual_GPU_tasks']=0;out['actual_GPU_workers']=0;out['source_count']=12
    return out


def build_costs(sealed):
    rows=[task_cost(x) for x in sealed['entries']]
    shared_path=Path(sealed['shared_setup'])
    shared=read(shared_path) if shared_path.exists() else None
    shared_counts=shared.get('physical_counts') if shared else None
    methods={m:cost_summary([r for r in rows if r['method']==m]) for m in METHODS}
    combined=cost_summary(rows)
    combined['shared_setup_charged_times']=1
    combined['shared_setup']=dict(source=str(shared_path),available=shared is not None,
        physical_counts=shared_counts,wall_seconds=shared.get('wall_seconds') if shared else None,
        cpu_seconds=shared.get('cpu_seconds') if shared else None,
        scope='Shared confirmation generation/operator preparation; charged once across all18 tasks, not once per method')
    combined['campaign_physical_totals_including_shared_once']={key:sum_known([
        combined[key]['total'],shared_counts.get(key,0) if shared_counts is not None else None]) for key in PHYSICAL_KEYS}
    imports=[]
    for row in rows:
        charge=row['solver_import_charge']
        if charge and charge.get('charged_once_per_worker'):
            imports.append(dict(pid=row['worker_pid'],slot=row['worker_slot'],**charge))
    combined['solver_import_charges_once_per_worker']=imports
    combined['solver_import_wall_seconds']=sum_known(x.get('wall_seconds') for x in imports)
    combined['solver_import_cpu_seconds']=sum_known(x.get('cpu_seconds') for x in imports)
    affinity_path=Path(sealed['remote_root'])/'results/confirmation/worker_affinity.json'
    affinity=read(affinity_path) if affinity_path.exists() else None
    expected_workers={x['pid'] for x in affinity} if isinstance(affinity,list) else None
    charged_workers={x['pid'] for x in imports}
    import_complete=expected_workers is not None and charged_workers==expected_workers and len(imports)==len(charged_workers)
    combined['worker_affinity']=affinity
    combined['solver_import_cost_complete']=import_complete
    if not import_complete:
        for key in ('solver_import_wall_seconds','solver_import_cpu_seconds'):
            combined[key].update(complete=False,total=None)
    combined['worker_time_not_makespan']=sum_known([combined['task_worker_wall_seconds']['total'],
        combined['solver_import_wall_seconds']['total']])
    run_path=ROOT/'audit/confirmation/REMOTE_RUN_RECEIPT.json'
    measured_run=read(run_path) if run_path.exists() else None
    combined['remote_campaign_elapsed_seconds']=measured_run.get('wall_seconds') if measured_run else None
    combined['remote_campaign_elapsed_source']=str(run_path) if measured_run else None
    combined['remote_campaign_elapsed_scope']=measured_run.get('scope') if measured_run else None
    combined['campaign_wall_seconds_for_cost']=combined['remote_campaign_elapsed_seconds']
    combined['currency_cost']=None
    combined['recorded_remote_CPU_seconds']=sum_known([combined['adapter_cpu_seconds']['total'],
        combined['solver_import_cpu_seconds']['total']])
    combined['recorded_CPU_plus_local_shared_generation_seconds']=sum_known([
        combined['recorded_remote_CPU_seconds']['total'],shared.get('cpu_seconds') if shared else None])
    combined['earlier_remote_bootstrap_import_seconds']=None
    combined['limitations']='Measured remote elapsed includes SSH/process startup and shutdown. Summed task/worker wall is work time, not concurrent elapsed; bootstrap cannot be separated, isolated per-job peak memory and currency tariff unavailable; do not infer GPU use from legacy cuda_workers=2'
    secondary=[dict(observation_id=ident,cost=cost_summary([r for r in rows if r['observation_id']==ident]))
               for ident in sorted({r['observation_id'] for r in rows})]
    return dict(task_rows=rows,primary_method_task_costs=methods,
        secondary_case_pools=secondary,all18_plus_shared_once=combined,
        secondary_pool_campaign_cost='Exactly all18_plus_shared_once; descriptive secondary is not another execution to add',
        primary_cost_comparison='Raw per-method task work shown. Shared generation may be shown once for each hypothetical standalone campaign, never tripled in actual aggregate.')


def score(args):
    if not args.authorize_confirmation_truth:
        raise PermissionError('Parent must explicitly execute score --authorize-confirmation-truth after sealing')
    paths=[args.output.with_suffix('.'+ext) for ext in ('json','csv','md')]
    if any(path.exists() for path in paths):raise FileExistsError('Preserve existing confirmation score; use a new output prefix')
    closed=read(args.seal)
    if closed.get('status')!='ALL18_DISPOSITIONS_SEALED' or len(closed['entries'])!=18:
        raise RuntimeError('All18 outer+inner terminal dispositions must be sealed first')
    for entry in closed['inputs']+[item for row in closed['entries'] for item in row['evidence']]:
        if not unchanged(entry):raise RuntimeError('Sealed evidence changed: '+entry['path'])
    # FIRST confirmation truth access occurs only here, after all18 are sealed.
    opened=utc();mapping=read(ROOT/'sealed_truth/confirmation/mapping_PRIVATE.json')
    families={r['observation_id']:r['family'] for r in mapping['records']}
    rows=[]
    for ident in sorted({entry['observation_id'] for entry in closed['entries']}):
        with np.load(ROOT/'sealed_truth/confirmation'/(ident+'.npz'),allow_pickle=False) as data:
            truth=data['m'].copy();support=data['mask'].ravel().astype(bool)
        public=public_metrics(ROOT/'data_public/confirmation'/(ident+'.npz'))
        for entry in [x for x in closed['entries'] if x['observation_id']==ident]:
            row=dict(observation_id=ident,family=families[ident],method=entry['method'],
                outer_status=entry['outer_status'],inner_status=entry['inner_status'],disposition=entry['disposition'],
                denominator_included=True,gate_evaluation='UNAVAILABLE',joint_gate=False,
                amplitude_normalized_nrmse=None,otsu_iou=None,relative_residual=None,
                whitened_residual=None,physical_objective=None,reason=entry['reason'])
            row['frozen_configuration_matches']=entry['frozen_configuration_matches']
            if entry['disposition']=='SEALED_PREDICTION':
                m,yp=prediction(entry['prediction'],(24,24))
                row.update(metrics(m,yp,public,truth,support,(24,24)),gate_evaluation='EVALUATED')
                inner=entry['outer'].get('result') or {}
                row['physical_objective_receipt_delta']=row['physical_objective']-inner['physical_objective'] if inner.get('physical_objective') is not None else None
                row['corrector_stop']=inner.get('corrector',{}).get('stop')
                row['complete_auxiliary_schedule']=inner.get('complete_auxiliary_schedule')
            rows.append(row)
    secondary=[]
    for selection in closed['secondary_pretruth_selection']:
        ident=selection['observation_id']
        candidates=[row for row in rows if row['observation_id']==ident and row['method'] in selection['retained_available_methods']]
        complete=selection['complete_three_arm_pool']
        item=dict(observation_id=ident,family=families[ident],method='secondary_min_original_Phi',
            denominator_included=True,selection_status=selection['selection_status'],
            selected_methods=selection['retained_available_methods'],complete_three_arm_pool=complete,
            minimum_original_Phi=selection['minimum_available_original_Phi'],
            candidate_metrics=candidates,
            joint_gate=bool(complete and candidates and all(r['joint_gate'] for r in candidates)),
            joint_gate_any_retained=bool(complete and candidates and any(r['joint_gate'] for r in candidates)),
            tie_policy=selection['tie_policy'],
            cost_policy='All3 attempted arms, including failures, plus shared campaign setup once')
        if complete and len(candidates)==1:
            item.update({key:candidates[0].get(key) for key in ('amplitude_normalized_nrmse','otsu_iou','relative_residual','whitened_residual','physical_objective')})
        else:
            item.update({key:None for key in ('amplitude_normalized_nrmse','otsu_iou','relative_residual','whitened_residual','physical_objective')})
        secondary.append(item)
    summaries=[]
    for method in METHODS+('secondary_min_original_Phi',):
        group=secondary if method=='secondary_min_original_Phi' else [r for r in rows if r['method']==method]
        summaries.append(dict(method=method,independent_objects=6,joint_pass=sum(r['joint_gate'] for r in group),
            denominator=6,available_scored_objects=sum(r.get('gate_evaluation')=='EVALUATED' for r in group) if method in METHODS else sum(r['complete_three_arm_pool'] for r in group),
            unavailable_or_incomplete=sum(r.get('gate_evaluation')!='EVALUATED' for r in group) if method in METHODS else sum(not r['complete_three_arm_pool'] for r in group),
            subgroup_counts=[dict(family=f,objects=sum(r['family']==f for r in group),joint_pass=sum(r['joint_gate'] for r in group if r['family']==f)) for f in sorted(set(families.values()))]))
    result=dict(status='CONFIRMATION_POST_SEAL_SCORED',first_confirmation_truth_opened_utc=opened,
        seal=str(args.seal),contract=frozen_contract(),rows=rows,secondary_rows=secondary,
        summaries=summaries,costs=build_costs(closed),
        scientific_scope='Six independent materials,18 tasks; descriptive counts only, no significance/population superiority claim',
        selection='Frozen original-Phi secondary set selected before truth, exact ties retained; no GT choice',
        new_physical_LU=0,new_complex_RHS=0,hash_checks=False)
    write_new(paths[0],result);write_csv_new(paths[1],rows+secondary)
    lines=['# R4 six-object confirmation scores','',
        'Post-seal independent scoring; six materials, three frozen arms per material. FAILED/UNKNOWN remain in each original six-object denominator. Descriptive counts only; no significance claim.','',
        '| Object | Family | Method | Disposition | RMSE/12 | IoU | Relative residual | Joint |',
        '|---|---|---|---|---:|---:|---:|---|']
    def f(x):return 'null' if x is None else f'{x:.9g}'
    for row in rows:
        lines.append(f"|{row['observation_id']}|{row['family']}|{row['method']}|{row['disposition']}|{f(row['amplitude_normalized_nrmse'])}|{f(row['otsu_iou'])}|{f(row['relative_residual'])}|{'PASS' if row['joint_gate'] else 'FAIL/UNAVAILABLE'}|")
    lines+=['','Secondary selection was sealed before GT. Exact-Phi ties are all retained, without GT tiebreak. A partial failed/UNKNOWN pool has unknown full-pool minimum; available endpoints remain visible and costs remain charged.','',
        '| Object | Retained available methods | Complete3 | Joint all retained |', '|---|---|---|---|']
    for row in secondary:lines.append(f"|{row['observation_id']}|{', '.join(row['selected_methods'])}|{row['complete_three_arm_pool']}|{row['joint_gate']}|")
    lines+=['','| Method | Joint pass / original objects | Unavailable/incomplete |','|---|---:|---:|']
    for row in summaries:lines.append(f"|{row['method']}|{row['joint_pass']} / 6|{row['unavailable_or_incomplete']}|")
    lines+=['','## Cost table','',
        'Top-level task physical totals only; nested setup/auxiliary/corrector are not added again. Each physical source solve has12 source RHS; small-block RHS are separate. Actual GPU tasks/workers=0 despite legacy runner metadata.','',
        '| Method task work | LU | Physical complex RHS | Worker wall sum s | Adapter CPU sum s |','|---|---:|---:|---:|---:|']
    for method,cost in result['costs']['primary_method_task_costs'].items():
        lines.append(f"|{method}|{f(cost['LU']['total'])}|{f(cost['actual_complex_RHS']['total'])}|{f(cost['task_worker_wall_seconds']['total'])}|{f(cost['adapter_cpu_seconds']['total'])}|")
    combined=result['costs']['all18_plus_shared_once']
    lines+=['', 'Actual total is all18 task work plus shared generation/setup **once**. The costed secondary pool is exactly that same total, not an additional run. Unknown costs retain null totals and known subtotals in JSON. Worker wall sums are not concurrent elapsed. Solver import charges are counted once per worker. Peak RSS has worker-lifetime scope.','',
        f"Measured remote elapsed used for time cost: {f(combined['remote_campaign_elapsed_seconds'])} s, including foreground SSH/process startup and shutdown; bootstrap cannot be separated. Recorded remote adapter+solver-import CPU: {f(combined['recorded_remote_CPU_seconds']['total'])} s; this excludes unmeasured bootstrap. Currency tariff/cost remains unknown.", '',
        f"Shared setup physical counts: `{json.dumps(combined['shared_setup']['physical_counts'])}`. Campaign physical totals including shared once: `{json.dumps(combined['campaign_physical_totals_including_shared_once'])}`.", '',
        'See JSON task rows for all failures, status pairs, per-task counters, local-operation counts, memory scope, source-RHS counts and incomplete costs. Scoring adds0 LU/0 RHS.']
    with paths[2].open('x') as stream:stream.write('\n'.join(lines)+'\n')
    for ident in sorted({r['observation_id'] for r in rows}):
        print(json.dumps(dict(observation_id=ident,rows=[r for r in rows if r['observation_id']==ident],secondary=[r for r in secondary if r['observation_id']==ident])))
    print(json.dumps(dict(status=result['status'],summaries=summaries,costs=combined)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['seal','score'])
    parser.add_argument('--remote-root',type=Path,default=ROOT/'results/confirmation_remote')
    parser.add_argument('--seal',type=Path,default=ROOT/'audit/confirmation/CONFIRMATION_PREDICTIONS_SEALED.json')
    parser.add_argument('--shared-setup',type=Path,default=ROOT/'audit/confirmation/GENERATION_RECEIPT.json')
    parser.add_argument('--dispositions',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'results/CONFIRMATION_SCORES')
    parser.add_argument('--authorize-confirmation-truth',action='store_true')
    args=parser.parse_args()
    if args.operation=='seal':seal(args)
    else:score(args)
