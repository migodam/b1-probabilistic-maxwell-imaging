"""Metadata-only R15 cost accounting. No arrays, GT, physics, hashes or science verdict."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import time

METHODS = ('g19','g19_s19','twofold','twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset','cie','dual_cie','dr')
MAIN='twofold_exposed24_v1'
FORK='fullstate_interventions_v1'


def numeric(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def write(path,data):
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def table(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list,tuple)) else v for k,v in row.items()})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path)
    # A later original-receipt copy can enrich costs without reading scientific arrays.
    p.add_argument('--packaging-receipt',type=Path);a=p.parse_args()
    root=a.root.resolve();out=(a.out or root/'analysis/cost_final_v1').resolve()
    if out.exists() or not out.is_relative_to(root/'analysis'):raise FileExistsError('New R15 analysis output required')
    start_cpu,start_wall=time.process_time(),time.perf_counter()
    read_files=[];processes=[];counter_rows=[];gaps=[];checks=[];tasks=[]
    corrector_rows=[];corrector_summary=[];termination_counts=[];render_attempts=[]
    def read(path):
        path=Path(path)
        if not path.is_absolute():path=root/path
        path=path.resolve()
        if not path.is_relative_to(root) or any(x.lower() in ('sealed_truth','truth','credentials','tokens','private_operational') for x in path.parts):raise ValueError('Only local admitted cost metadata')
        if path.suffix.lower()!='.json':raise ValueError('Only JSON metadata is read')
        read_files.append(str(path.relative_to(root)));return json.loads(path.read_text(encoding='utf-8'))
    def add_process(ident,category,source,data,host='LOCAL_APPLE_CPU',cpu=None,cpu_scope=None,body=None,wall=None,notes=None,status=None):
        if body is None:body=data.get('CPU_seconds',data.get('analysis_CPU_seconds',data.get('recorded_attempt_CPU_seconds')))
        if cpu_scope is None:
            preferences=[('whole_actual_OS_child_CPU_seconds','OS_WHOLE_PROCESS'),('actual_OS_child_CPU_seconds','OS_WHOLE_PROCESS'),('whole_process_CPU_seconds','PROCESS_CPU_SNAPSHOT'),
                ('whole_process_CPU_seconds_at_receipt','PROCESS_CPU_SNAPSHOT'),('whole_process_CPU_at_receipt_seconds','PROCESS_CPU_SNAPSHOT'),
                ('whole_process_CPU_at_receipt','PROCESS_CPU_SNAPSHOT')]
            for key,scope in preferences:
                if numeric(data.get(key)):cpu=data[key];cpu_scope=scope;break
            if cpu_scope is None:cpu=body;cpu_scope='BODY_CPU_ONLY' if numeric(body) else 'UNKNOWN'
        if wall is None:wall=data.get('whole_actual_wall_seconds',data.get('wall_seconds',data.get('analysis_wall_seconds',data.get('wall_since_imports'))))
        processes.append(dict(process_id=ident,category=category,host=host,source=str(source),status=status or data.get('status'),
            preferred_CPU_seconds=cpu if numeric(cpu) else None,CPU_scope=cpu_scope,
            body_CPU_seconds=body if numeric(body) else None,body_nested_not_additive=True,
            wall_seconds=wall if numeric(wall) else None,wall_intervals_not_added=True,
            CPU_null_reason=None if numeric(cpu) else 'Original process CPU not measured/present; not imputed0',
            notes=notes,additive_process_row=True))
    def counters(ident,category,namespace,values,source,scope='FRESH_EXECUTED'):
        if not isinstance(values,dict):return
        for key,value in values.items():
            if not numeric(value):continue
            if 'peak' in key or key.endswith('_rank') or key.endswith('_bytes'):aggregation='MAX_OR_SNAPSHOT_NOT_SUM'
            elif 'second' in key or 'microsecond' in key:aggregation='NESTED_TIME_NOT_ADDED_TO_CPU'
            elif 'dimension' in key:aggregation='SUM_DIMENSIONS_NOT_A_SINGLE_FACTOR_DIMENSION'
            else:aggregation='SUM_RECORDED_COUNTER'
            if namespace=='physical':
                family='MAXWELL_LU' if key=='LU' else 'MAXWELL_AGGREGATE_RHS' if key=='actual_complex_RHS' else 'PHYSICAL_ROLE_OVERLAPS_AGGREGATES'
            elif 'H_' in key or key=='H_cholesky':family='BACKGROUND_H_CHOLESKY_OR_TRIANGULAR_RHS'
            elif 'FFT' in key or 'IFFT' in key:family='AUX_FFT_TRANSFORMS_OR_CALLS'
            elif 'small_LU' in key or 'small_triangular' in key:family='AUX_SMALL_LU_OR_TRIANGULAR_RHS'
            elif 'setup' in key or 'SVD' in key or 'eigh' in key or 'Cholesky' in key:family='AUX_ACTUAL_SETUP'
            else:family='OTHER_AUX_COUNTER'
            counter_rows.append(dict(process_id=ident,category=category,namespace=namespace,counter=key,value=value,
                                     family=family,aggregation=aggregation,scope=scope,source=str(source)))
    def auxiliary(ident,category,data,source):
        if not isinstance(data,dict):return
        if isinstance(data.get('operation_counts'),dict):
            counters(ident,category,'auxiliary',data['operation_counts'],source)
            counters(ident,category,'adapter_diagnostic',data.get('adapter_diagnostic_cost'),source)
            counters(ident,category,'adapter_denominator_setup',data.get('denominator_setup_cost'),source)
        else:counters(ident,category,'auxiliary',data,source)

    ledger=read(f'analysis/{MAIN}/COST_LEDGER.json');queue_costs=read(f'analysis/{MAIN}/QUEUE_COSTS.json')
    qm=read(f'queues/{MAIN}/RECEIPT.json');qf=read(f'queues/{FORK}/RECEIPT.json')
    for campaign,q,count in [(MAIN,qm,264),(FORK,qf,56)]:
        if q.get('status') not in ('COMPLETED','COMPLETED_WITH_ERRORS') or q.get('running') or q.get('pending_tasks',0) or len(q.get('jobs',[]))!=count:raise ValueError('Fixed terminal queue required')
    expected={(f'r14_new_{i:04d}',m) for i in range(24) for m in METHODS}
    if len(ledger)!=264 or {(r['case_id'],r['solver']) for r in ledger}!=expected:raise ValueError('Full24x11 original cost ledger required')
    jobs={j['task_id']:j for j in qm['jobs']}
    for row in ledger:
        ident=row['task_id'];job=jobs[ident];external=job.get('external',{}).get('child_process_CPU_seconds')
        if external!=row.get('external_child_CPU_seconds'):raise ValueError('Main ledger/external receipt mismatch: '+ident)
        add_process(ident,'online_main',f'queues/{MAIN}/RECEIPT.json',{},'REMOTE_WINDOWS_CPU',external,'OS_WHOLE_PROCESS',
                    row.get('internal_final_process_CPU_seconds'),job.get('wall_seconds'),status=job.get('status'))
        physical=row.get('physical_cost_raw') or {};aux=row.get('auxiliary_cost') or {};projection=row.get('fixed_projection_cost') or {}
        counters(ident,'online_main','physical',physical,f'analysis/{MAIN}/COST_LEDGER.json')
        auxiliary(ident,'online_main',aux,f'analysis/{MAIN}/COST_LEDGER.json')
        counters(ident,'online_main','background_projection',projection,f'analysis/{MAIN}/COST_LEDGER.json')
        lu=physical.get('LU');ds=physical.get('full_factorization_dimension')
        tasks.append(dict(campaign=MAIN,task_id=ident,case_id=row['case_id'],method=row['solver'],queue_status=job.get('status'),
            run_status=row.get('run_status'),returncode=job.get('returncode'),external_CPU_seconds=external,
            internal_final_process_CPU_seconds=row.get('internal_final_process_CPU_seconds'),
            physical_LU=lu,physical_complex_RHS=physical.get('actual_complex_RHS'),
            factor_dimension_sum=ds,single_factor_dimension=ds/lu if numeric(ds) and numeric(lu) and lu>0 else None,
            termination_reason=row.get('termination_reason')))
        history_path=root/'RESULTS'/MAIN/ident/'HISTORY.json'
        audit_row=dict(task_id=ident,case_id=row['case_id'],method=row['solver'],history_path=str(history_path.relative_to(root)),
            first_completed_corrector_audit_CPU_seconds=None,first_completed_corrector_audit_candidate=None,
            first_completed_corrector_audit_eligible=None,corrector_audits=None,eligible_corrector_audits=None,
            nominal_remaining_to_120_after_first_audit=None,proposal_termination=None,first_audit_is_exact_start=False,
            raw_best_candidate_CPU_used_as_proposal_end=False)
        term=row.get('termination_reason')
        audit_row['proposal_termination']=term.get('proposal') if isinstance(term,dict) else None
        if history_path.is_file():
            history=read(history_path)
            completed=[h for h in history if h.get('phase')=='corrector' and numeric(h.get('CPU_seconds'))]
            audit_row.update(corrector_audits=len(completed),eligible_corrector_audits=sum(h.get('eligible') is True for h in completed),
                null_reason=None if completed else 'No saved completed corrector audit with finite CPU')
            if completed:
                first=min(completed,key=lambda h:h['CPU_seconds'])
                audit_row.update(first_completed_corrector_audit_CPU_seconds=first['CPU_seconds'],
                    first_completed_corrector_audit_candidate=first.get('candidate'),
                    first_completed_corrector_audit_eligible=first.get('eligible'),
                    nominal_remaining_to_120_after_first_audit=120-first['CPU_seconds'])
        else:audit_row['null_reason']='Original HISTORY missing; no first-audit time inferred'
        corrector_rows.append(audit_row)
    fork_task_records=[]
    for job in qf['jobs']:
        ident=job['task_id'];external=job.get('external',{}).get('child_process_CPU_seconds')
        add_process(ident,'online_fork',f'queues/{FORK}/RECEIPT.json',{},'REMOTE_WINDOWS_CPU',external,'OS_WHOLE_PROCESS',
                    job.get('child_internal_CPU_seconds'),job.get('wall_seconds'),status=job.get('status'))
        counters(ident,'online_fork','physical',job.get('physical_cost'),f'queues/{FORK}/RECEIPT.json')
        for prefix,namespace in [('auxiliary','auxiliary'),('fixed_projection','background_projection')]:
            inherited=job.get(prefix+'_inherited');cumulative=job.get(prefix+'_cumulative')
            if isinstance(inherited,dict) and isinstance(cumulative,dict):
                fresh={k:cumulative.get(k,0)-inherited.get(k,0) for k in set(cumulative)|set(inherited)
                       if numeric(cumulative.get(k,0)) and numeric(inherited.get(k,0))}
                checks.append(dict(check='fork_fresh_counter_matches_queue',task_id=ident,namespace=namespace,
                    matches=all(fresh.get(k,0)==v for k,v in (job.get(prefix+'_new') or {}).items())))
                counters(ident,'online_fork',namespace,fresh,f'queues/{FORK}/RECEIPT.json','CUMULATIVE_MINUS_INHERITED_ONLY')
            else:gaps.append(dict(id=ident+'_'+namespace,reason='Missing inherited/cumulative counter pair',value=None))
        fork_task_records.append(dict(task_id=ident,method=job.get('method'),case_id=job.get('case_id'),status=job.get('status'),
            returncode=job.get('returncode'),external_CPU_seconds=external,internal_final_process_CPU_seconds=job.get('child_internal_CPU_seconds')))
    for campaign,q in [(MAIN,qm),(FORK,qf)]:
        values=[j.get('external',{}).get('child_process_CPU_seconds') for j in q['jobs']]
        total=sum(v for v in values if numeric(v));checks.append(dict(check='external_CPU_recomputed_from_all_terminal_jobs',campaign=campaign,
            recomputed=total,reported=q.get('finished_child_CPU_seconds'),matches=total==q.get('finished_child_CPU_seconds'),missing=sum(not numeric(v) for v in values)))
        add_process(campaign+'_controller','online_controller',f'queues/{campaign}/RECEIPT.json',{},'REMOTE_WINDOWS_CPU',
            q.get('controller_CPU_seconds'),'CONTROLLER_BODY_CPU',wall=q.get('queue_wall_seconds'),notes='Disjoint from child CPU; queue wall is an interval, never summed with child wall')
    checks.append(dict(check='main_QUEUE_COSTS_controller_matches_original',matches=queue_costs[0]['controller_CPU_seconds']==qm['controller_CPU_seconds']))

    # Standalone validation attempts, never add nested total receipts again.
    for suffix in ['ATTEMPT01_RECEIPT.json','RECEIPT.json']:
        src='audit/denominator_wiring/'+suffix;d=read(src);ident='denominator_'+suffix
        add_process(ident,'validation',src,d)
        counters(ident,'validation','auxiliary',d.get('operation_cost'),src)
        counters(ident,'validation','setup_reassignment',d.get('counter_reassignment_setup_cost'),src)
    read('audit/denominator_wiring/TOTAL_COST.json')
    gaps.append(dict(id='denominator_attempt01_counter_reassignment',value=None,reason='Attempt01 preserves empty reassignment setup counter; no numeric setup count imputed; CPU still charged'))
    src='audit/integration_v1/RESULTS.json';d=read(src);add_process('integration_all7','validation',src,d,body=d.get('recorded_attempt_CPU_seconds'),wall=d.get('wall_since_imports'))
    counters('integration_all7','validation','physical',{'LU':d['physical_LU_total'],'actual_complex_RHS':d['physical_actual_complex_RHS_total']},src)
    for rp in sorted((root/'audit/integration_v1').glob('*/RECEIPT.json')):
        v=read(rp);ident='integration_'+rp.parent.name
        # The primary solver counter is fresh; restored copies retain history.
        aux=v.get('auxiliary_cost')
        if aux is None:aux=(v.get('extra') or {}).get('native_receipt')
        auxiliary(ident,'validation',aux,str(rp.relative_to(root)))
        counters(ident,'validation','background_projection',v.get('projector_cost'),str(rp.relative_to(root)))
        extra=v.get('extra') or {}
        if isinstance(extra.get('restore_auxiliary_cost'),dict):
            inherited=extra['cost_at17'];cum=extra['restore_auxiliary_cost']
            counters(ident+'_restore_increment','validation','auxiliary',
                {k:cum.get(k,0)-inherited.get(k,0) for k in set(cum)|set(inherited)},str(rp.relative_to(root)),'RESTORE_INCREMENT_NOT_CUMULATIVE')
        if extra.get('pending_cycle_restore'):
            replay=extra['pending_cycle_restore'];ref=replay['reference_auxiliary_cost'];cum=replay['replay_auxiliary_cost']
            # BEFORE reference not retained in receipt; AFTER subtraction would
            # erase actual replay. Preserve both as nonadditive provenance.
            counters(ident+'_pending_replay_cumulative','validation','auxiliary',cum,str(rp.relative_to(root)),'INHERITED_CUMULATIVE_NOT_ADDITIVE')
            gaps.append(dict(id='integration_pending_replay_fresh_aux',value=None,reason='No separately retained BEFORE counter in receipt; replay CPU and physical work charged in integration total; cumulative counts not added'))

    # Operator characterization groups were actually computed once per unique cache,
    # not once for every execution that reused exact-equal G/S.
    for campaign in ['rank_characterization/current_complex_v1','rank_characterization/real_material_v1','rank_all24/complex','rank_all24/real']:
        src='audit/'+campaign+'/RECEIPT.json';d=read(src);add_process(campaign,'setup',src,d,wall=d.get('wall_since_before_numpy_scipy_imports',d.get('wall_seconds')))
        for mp in sorted((root/'audit'/campaign).glob('operator_group_*/METADATA.json')):
            md=read(mp);counters(campaign+'_'+mp.parent.name,'setup','operator_setup',md.get('operation_counts'),str(mp.relative_to(root)))
    for name in ['old_endpoint_geometry_v1','readout_counterfactual_v1','paid_endpoint_diagnostic_v1']:
        src='audit/'+name+'/RECEIPT.json';d=read(src);add_process(name,'validation',src,d)
        phy=d.get('physical_cost')
        if phy is None:phy={'LU':d.get('physical_LU_total'),'actual_complex_RHS':d.get('actual_complex_RHS_total')}
        counters(name,'validation','physical',phy,src)
    cf=root/'audit/readout_counterfactual_v1'
    for rp in sorted(cf.glob('*/*_RECEIPT.json')):
        d=read(rp)
        counters('counterfactual_'+rp.parent.name+'_'+rp.stem,'validation','auxiliary',d.get('fresh_auxiliary_cost'),str(rp.relative_to(root)))
    src='audit/old_endpoint_geometry_v1/METADATA_ENRICHMENT_RECEIPT.json';d=read(src);add_process('old_endpoint_metadata_enrichment','validation_metadata',src,d)

    # Immediate analysis receipts only: archive copies and initial plot subfolders
    # are nested provenance, not independent process charges.
    for rp in sorted((root/'analysis').glob('*/*RECEIPT*.json')):
        if rp.parent==out:continue
        d=read(rp);name=rp.parent.name
        if name=='old_geometry_summary':category='plotting_and_analysis';host='LOCAL_APPLE_CPU'
        elif name=='fullstate_v1':category='plotting_and_analysis';host='REMOTE_WINDOWS_CPU'
        elif name==MAIN:category='scoring_and_plotting';host='LOCAL_APPLE_CPU'
        elif name=='intervention_scoring_v1':category='postseal_scoring';host='LOCAL_APPLE_CPU'
        elif name.startswith('cost_final'):category='cost_accounting_metadata';host='LOCAL_APPLE_CPU'
        else:category='offline_analysis';host='LOCAL_APPLE_CPU'
        add_process('analysis_'+name,category,str(rp.relative_to(root)),d,host)
        if name=='selection_geometry_v1':
            for record in d.get('model_construction_cost',[]):
                counters('selection_geometry_'+record['case_id'],'offline_analysis','model_constructor_cache',record.get('cost'),str(rp.relative_to(root)))
            counters('selection_geometry_regularizer','offline_analysis','regularizer',{'regularizer_evaluations':d.get('regularizer_evaluations')},str(rp.relative_to(root)))
    src='audit/analysis_attempts/main_scoring_attempt01/FAILED_RECEIPT.json';d=read(src);add_process('failed_main_scoring_attempt01','scoring_failed',src,d)
    for rp in sorted((root/'audit/analysis_attempts').glob('cost_summary_attempt*/RECEIPT.json')):
        d=read(rp);add_process(rp.parent.name,'cost_accounting_metadata',str(rp.relative_to(root)),d)
    src='audit/PACKAGING_REPAIR_V2.json';d=read(src)
    add_process('failed_raw_packaging_v1','publication_packaging',src,d,'REMOTE_WINDOWS_CPU',d.get('failed_external_CPU_seconds'),
        'OS_WHOLE_PROCESS',wall=d.get('failed_wall_seconds'),status='FAILED_BEFORE_ANY_ZIP',notes='Original failure preserved; zero ZIP output is not zero CPU')
    packaging_path=a.packaging_receipt
    if packaging_path is None:
        candidates=[root/'audit/remote_packaging_v2/PACKAGE_RECEIPT.json',root/'archives/r15_raw_v2/PACKAGE_RECEIPT.json',root/'audit/raw_packaging_v2/PACKAGE_RECEIPT.json']
        packaging_path=next((x for x in candidates if x.is_file()),None)
    if packaging_path:
        d=read(packaging_path);add_process('raw_packaging_v2','publication_packaging',str(packaging_path),d,'REMOTE_WINDOWS_CPU')
        checks.append(dict(check='raw_v2_original_receipt_present',matches=True,parts=len(d.get('parts',[])),
                           archive_bytes=d.get('actual_archive_bytes',sum(x.get('bytes',x.get('actual_bytes',0)) for x in d.get('parts',[])))))
        checks.append(dict(check='raw_v2_part_bytes_sum_matches_original',
            matches=sum(x.get('actual_bytes',x.get('bytes',0)) for x in d.get('parts',[]))==d.get('actual_archive_bytes'),
            parts=len(d.get('parts',[])),archive_bytes=d.get('actual_archive_bytes'),
            original_body_CPU_seconds=d.get('CPU_seconds'),process_snapshot_CPU_seconds=d.get('whole_process_CPU_at_receipt_seconds')))
    else:
        add_process('raw_packaging_v2_PENDING','publication_packaging','PENDING_ORIGINAL_REMOTE_PACKAGE_RECEIPT',{},'REMOTE_WINDOWS_CPU',
            cpu=None,cpu_scope='UNKNOWN',notes='Parent reports139.140625CPU /217.6562wall /16parts /7364794746bytes; raw receipt not local, excluded from verified numeric sums')
        gaps.append(dict(id='raw_packaging_v2_original_receipt',value=None,reason='Copy remote archives/r15_raw_v2/PACKAGE_RECEIPT.json, then use --packaging-receipt with fresh --out; no science recomputation',
            parent_reported_unverified=dict(CPU_seconds=139.140625,wall_seconds=217.6562,parts=16,bytes=7364794746)))
    for rp in sorted((root/'archives').glob('*/R15_REVIEW_EXPORT_V1_RECEIPT.json')):
        d=read(rp);add_process('review_export_v1','publication_IO',str(rp.relative_to(root)),d,'REMOTE_WINDOWS_CPU')
    remote=root/'audit/remote_packaging_v2'
    review_path=remote/'REVIEW_EXPORT_RECEIPT.json'
    if review_path.is_file():
        d=read(review_path)
        present=next((x for x in processes if x['process_id']=='review_export_v1'),None)
        if present:
            checks.append(dict(check='review_export_duplicate_receipt_CPU_matches',matches=present['body_CPU_seconds']==d.get('CPU_seconds')))
            present['duplicate_source_not_added']=str(review_path.relative_to(root))
        else:add_process('review_export_v1','publication_IO',str(review_path.relative_to(root)),d,'REMOTE_WINDOWS_CPU')
    pipeline_path=remote/'POST_MAIN_PIPELINE.json'
    if pipeline_path.is_file():
        pipe=read(pipeline_path)
        mapping={'fork56':FORK+'_controller','review_export':'review_export_v1','raw_package':'failed_raw_packaging_v1'}
        for stage in pipe.get('stages',[]):
            ident=mapping.get(stage['stage']);target=next((x for x in processes if x['process_id']==ident),None)
            external=stage.get('external',{}).get('child_process_CPU_seconds')
            if target is not None and numeric(external):
                if target['CPU_scope']=='OS_WHOLE_PROCESS':checks.append(dict(check='pipeline_external_duplicate_matches',process_id=ident,matches=target['preferred_CPU_seconds']==external))
                target['nested_previous_CPU_seconds']=target['preferred_CPU_seconds']
                target['nested_previous_CPU_scope']=target['CPU_scope']
                if ident==FORK+'_controller':target['body_CPU_seconds']=target['preferred_CPU_seconds']
                target.update(preferred_CPU_seconds=external,CPU_scope='OS_WHOLE_PROCESS',
                    external_source=str(pipeline_path.relative_to(root)),external_stage=stage['stage'],
                    external_returncode=stage.get('returncode'),wall_seconds=stage.get('wall_seconds'))
        add_process('post_main_pipeline_wrapper','publication_IO',str(pipeline_path.relative_to(root)),pipe,'REMOTE_WINDOWS_CPU',
            cpu=pipe.get('wrapper_CPU_seconds'),cpu_scope='WRAPPER_BODY_CPU',wall=pipe.get('wall_seconds'),
            notes='Wrapper CPU alone; stage child CPU represented in existing independent rows, never duplicated')
    failed_package_path=remote/'FAILED_PACKAGE_ATTEMPT01.json'
    if failed_package_path.is_file():
        d=read(failed_package_path);target=next(x for x in processes if x['process_id']=='failed_raw_packaging_v1')
        target['body_CPU_seconds']=d.get('CPU_seconds');target['nested_failure_process_snapshot_CPU_seconds']=d.get('whole_process_CPU_at_error_seconds')
        target['original_failed_package_receipt']=str(failed_package_path.relative_to(root))
    download=remote/'DOWNLOAD_RECEIPT.json'
    if download.is_file():
        d=read(download);add_process('remote_packaging_metadata_download','publication_IO',str(download.relative_to(root)),d)
        gaps.append(dict(id='remote_packaging_metadata_download_child_CPU',value=None,reason='Metadata download parent CPU measured; separate SSH/SCP child CPU not reported'))
    render=root/'FIGURES/final/RENDER_ATTEMPTS.jsonl'
    if render.is_file():
        read_files.append(str(render.relative_to(root)))
        for index,line in enumerate(render.read_text(encoding='utf-8').splitlines(),1):
            if not line.strip():continue
            d=json.loads(line);ident=f'final_render_attempt_{index:02d}'
            add_process(ident,'plotting',str(render.relative_to(root)),d)
            render_attempts.append(dict(attempt=index,status=d.get('status'),body_CPU_seconds=d.get('CPU_seconds'),
                whole_process_CPU_snapshot_seconds=d.get('whole_process_CPU_seconds_at_receipt'),wall_seconds=d.get('wall_seconds'),
                artifacts=len(d.get('artifacts',[])),started_utc=d.get('started_utc'),finished_utc=d.get('finished_utc'),
                nested_timer_policy='Use one process snapshot per render attempt, not body+snapshot'))
    # Recorded metadata/deployment/IO/static processes: retain actual measurements,
    # including unknown failed launches/partial transfers instead of inserting0.
    named=['LOCAL_ENVIRONMENT.json','MANIFEST_FREEZE_twofold_exposed24_v1.json','PUBLIC_ACQUISITION_SUMMARY.json',
        'RUNNER_IMPLEMENTATION.json','QUEUE_INTERVENTIONS_IMPLEMENTATION.json','QUEUE_INTERVENTIONS_SOURCE512_STATIC.json',
        'ANALYZE_FULLSTATE_IMPLEMENTATION.json','COMPARE_READOUTS_IMPLEMENTATION.json','PACKAGE_RAW_STATIC_CHECK.json',
        'PAID_ENDPOINT_DIAGNOSTIC_STATIC_CHECK.json','REVIEW_EXPORT_STATIC_CHECK.json','REMOTE_DEPLOY_V1.json',
        'REMOTE_FOLLOWUP_DEPLOY_V1.json','FULLSTATE_ANALYSIS_IMPORT.json','SCORING_DOWNLOAD_RECEIPT.json',
        'REMOTE_LAUNCH_FAILURE_V1.json','REVIEW_TRANSFER_REPAIR.json','PACKAGE_REPAIR_DEPLOY_V2.json']
    for filename in named:
        rp=root/'audit'/filename
        if not rp.is_file():gaps.append(dict(id=filename,value=None,reason='Expected metadata receipt absent'));continue
        d=read(rp);category='publication_IO' if any(x in filename for x in ['DEPLOY','DOWNLOAD','TRANSFER','IMPORT','LAUNCH']) else 'setup_implementation_metadata'
        add_process(filename,category,'audit/'+filename,d)
        if numeric(d.get('SCP_child_CPU_seconds')):
            add_process(filename+'_SCP_CHILD',category,'audit/'+filename,{},cpu=d['SCP_child_CPU_seconds'],cpu_scope='OS_CHILD_CPU',
                notes='Separate SCP subprocess CPU; parent process CPU is disjoint')
        if not numeric(processes[-1]['preferred_CPU_seconds']):gaps.append(dict(id=filename,value=None,reason=d.get('null_reason',d.get('transfer_CPU_null_reason','No original process CPU present'))))
        if category=='publication_IO' and numeric(d.get('CPU_seconds')) and not numeric(d.get('SCP_child_CPU_seconds')):
            gaps.append(dict(id=filename+'_unreported_child_IO_CPU',value=None,
                reason='Recorded parent/body CPU retained; separate SSH/SCP child CPU not reported by this metadata, not assumed included'))
    gaps += [dict(id='currency_cost',value=None,reason='No metered currency record'),
             dict(id='LLM_billing',value=None,reason='No token/billing record; not imputed0'),
             dict(id='publication_upload_CPU',value=None,reason='No original upload/build-site OS receipt in current evidence; do not claim zero or completion'),
             dict(id='plotting_separate_from_scoring',value=None,reason='Several actual receipts cover scoring+plotting together; no fabricated split'),
             dict(id='cross_host_speed',value=None,reason='Local Apple validation and remote Windows inverse CPU are different workloads/hosts; no absolute speed claim')]

    methods=[]
    for method in METHODS:
        rr=[r for r in tasks if r['method']==method];ext=[r['external_CPU_seconds'] for r in rr]
        methods.append(dict(method=method,expected_tasks=24,actual_terminal_rows=len(rr),queue_status_counts=dict(Counter(r['queue_status'] for r in rr)),
            run_status_counts=dict(Counter(r['run_status'] for r in rr)),external_CPU_seconds=sum(x for x in ext if numeric(x)),
            missing_external_CPU=sum(not numeric(x) for x in ext),internal_CPU_snapshot_seconds=sum(r['internal_final_process_CPU_seconds'] for r in rr if numeric(r['internal_final_process_CPU_seconds'])),
            physical_LU=sum(r['physical_LU'] for r in rr if numeric(r['physical_LU'])),
            actual_complex_RHS=sum(r['physical_complex_RHS'] for r in rr if numeric(r['physical_complex_RHS'])),
            factor_dimension_sum=sum(r['factor_dimension_sum'] for r in rr if numeric(r['factor_dimension_sum'])),
            single_factor_dimensions=sorted(set(r['single_factor_dimension'] for r in rr if numeric(r['single_factor_dimension'])))))
    def quantile(values,q):
        if not values:return None
        values=sorted(values);at=(len(values)-1)*q;lo=int(at);hi=min(lo+1,len(values)-1)
        return values[lo]+(values[hi]-values[lo])*(at-lo)
    for method in METHODS:
        rr=[r for r in corrector_rows if r['method']==method]
        values=[r['first_completed_corrector_audit_CPU_seconds'] for r in rr if numeric(r['first_completed_corrector_audit_CPU_seconds'])]
        reasons=Counter(json.dumps(r['proposal_termination'],sort_keys=True,ensure_ascii=False) for r in rr)
        corrector_summary.append(dict(method=method,expected_tasks=24,with_completed_corrector_audit=len(values),missing=24-len(values),
            min_CPU_seconds=min(values) if values else None,p25_CPU_seconds=quantile(values,.25),median_CPU_seconds=quantile(values,.5),
            mean_CPU_seconds=sum(values)/len(values) if values else None,p75_CPU_seconds=quantile(values,.75),max_CPU_seconds=max(values) if values else None,
            quantity='Minimum saved HISTORY CPU_seconds with phase=corrector; completed audit time, NOT exact corrector start/proposal end',
            nominal_remaining_to_120_median=quantile([120-v for v in values],.5)))
        for reason,count in sorted(reasons.items()):termination_counts.append(dict(method=method,expected_tasks=24,proposal_termination=json.loads(reason),count=count,unknown=json.loads(reason) is None))
    aggregates=[]
    for key in sorted({(r['category'],r['namespace'],r['counter'],r['family'],r['aggregation']) for r in counter_rows}):
        rr=[r for r in counter_rows if (r['category'],r['namespace'],r['counter'],r['family'],r['aggregation'])==key and r['scope']!='INHERITED_CUMULATIVE_NOT_ADDITIVE']
        if not rr:continue
        aggregation=key[-1];values=[r['value'] for r in rr]
        aggregates.append(dict(category=key[0],namespace=key[1],counter=key[2],family=key[3],aggregation=aggregation,
            recorded_rows=len(rr),value=max(values) if aggregation=='MAX_OR_SNAPSHOT_NOT_SUM' else sum(values),
            additive_to_CPU=False,additive_to_Maxwell_RHS=key[2]=='actual_complex_RHS' and key[1]=='physical'))
    method_counters=[]
    for method in METHODS:
        ids={r['task_id'] for r in tasks if r['method']==method};groups={}
        for row in counter_rows:
            if row['category']=='online_main' and row['process_id'] in ids:
                groups.setdefault((row['namespace'],row['counter'],row['family'],row['aggregation']),[]).append(row['value'])
        for key,values in sorted(groups.items()):
            method_counters.append(dict(method=method,expected_tasks=24,namespace=key[0],counter=key[1],family=key[2],aggregation=key[3],
                recorded_rows=len(values),value=max(values) if key[3]=='MAX_OR_SNAPSHOT_NOT_SUM' else sum(values),
                absent_counter_policy='Only recorded original counters aggregated; absent/unsupported counter not fabricated',additive_to_CPU=False))
    cpu_scopes=[]
    for key in sorted({(r['category'],r['host'],r['CPU_scope']) for r in processes}):
        rr=[r for r in processes if (r['category'],r['host'],r['CPU_scope'])==key];known=[r['preferred_CPU_seconds'] for r in rr if numeric(r['preferred_CPU_seconds'])]
        cpu_scopes.append(dict(category=key[0],host=key[1],CPU_scope=key[2],processes=len(rr),known_CPU_subtotal_seconds=sum(known),
            missing_CPU_processes=len(rr)-len(known),complete_CPU_total_seconds=sum(known) if len(known)==len(rr) else None))
    out.mkdir(parents=True,exist_ok=False)
    outputs={'TASK_COSTS264':tasks,'FORK_COSTS56':fork_task_records,'METHODS11':methods,'PROCESS_COSTS':processes,
        'CPU_BY_CATEGORY_SCOPE':cpu_scopes,'COUNTERS_RECORDED':counter_rows,'COUNTERS_BY_CATEGORY':aggregates,
        'METHOD_COUNTERS11':method_counters,'GAPS':gaps,'RECONCILIATION':checks,
        'FIRST_CORRECTOR_AUDIT264':corrector_rows,'FIRST_CORRECTOR_AUDIT_BY_METHOD':corrector_summary,
        'PROPOSAL_TERMINATION_COUNTS':termination_counts,'FINAL_RENDER_ATTEMPTS':render_attempts}
    for name,data in outputs.items():write(out/(name+'.json'),data);table(out/(name+'.csv'),data)
    totals=dict(main_external_CPU_recomputed=sum(j['external']['child_process_CPU_seconds'] for j in qm['jobs']),
        fork_external_CPU_recomputed=sum(j['external']['child_process_CPU_seconds'] for j in qf['jobs']),
        main_controller_CPU=qm['controller_CPU_seconds'],fork_controller_CPU=qf['controller_CPU_seconds'],
        main_physical_LU=sum(r['physical_LU'] for r in methods),main_actual_complex_RHS=sum(r['actual_complex_RHS'] for r in methods),
        fork_physical_LU=sum(j.get('physical_cost',{}).get('LU',0) for j in qf['jobs']),
        fork_actual_complex_RHS=sum(j.get('physical_cost',{}).get('actual_complex_RHS',0) for j in qf['jobs']),
        validation_physical_LU=sum(r['value'] for r in aggregates if r['category']=='validation' and r['namespace']=='physical' and r['counter']=='LU'),
        validation_actual_complex_RHS=sum(r['value'] for r in aggregates if r['category']=='validation' and r['namespace']=='physical' and r['counter']=='actual_complex_RHS'),
        all_cost_currency=None,LLM_billing=None,GPU_used=False,GPU_scope='Recorded CPU-only programs, single BLAS thread; no CUDA calls/installation; no GPU billing inferred',
        whole_study_CPU_total=None,whole_study_CPU_null_reason='Pending/unmeasured independent process costs and mixed body/snapshot/OS scopes; category subtotals preserved separately')
    by_id={r['process_id']:r for r in processes}
    totals['fork_controller_external_CPU']=by_id[FORK+'_controller']['preferred_CPU_seconds'] if by_id[FORK+'_controller']['CPU_scope']=='OS_WHOLE_PROCESS' else None
    totals['review_export_external_CPU']=by_id.get('review_export_v1',{}).get('preferred_CPU_seconds') if by_id.get('review_export_v1',{}).get('CPU_scope')=='OS_WHOLE_PROCESS' else None
    totals['raw_v2_body_CPU']=by_id.get('raw_packaging_v2',{}).get('body_CPU_seconds')
    totals['raw_v2_process_CPU_snapshot']=by_id.get('raw_packaging_v2',{}).get('preferred_CPU_seconds')
    totals['final_render_attempts']=len(render_attempts)
    totals['final_render_body_CPU']=sum(r['body_CPU_seconds'] for r in render_attempts if numeric(r['body_CPU_seconds']))
    totals['final_render_process_CPU_snapshots']=sum(r['whole_process_CPU_snapshot_seconds'] for r in render_attempts if numeric(r['whole_process_CPU_snapshot_seconds']))
    write(out/'TOTALS.json',totals)
    lines=['# R15最终机械费用汇总','',
        '仅读取本地JSON费用/队列元数据，不打开材料/GT数组，不运行physics、优化、评分或远端任务。费用不代表恢复质量或科学结论。', '',
        f'主264任务外部OS整进程CPU重算：{totals["main_external_CPU_recomputed"]:.6f}s；fork56：{totals["fork_external_CPU_recomputed"]:.6f}s。两组控制器记录CPU分别{totals["main_controller_CPU"]:.6f}s、{totals["fork_controller_CPU"]:.6f}s，与子进程分离。保留11方法×24完整分母、全部终态/失败/未知项。', '',
        f'Main实际Maxwell LU / aggregate complex RHS：{totals["main_physical_LU"]} / {totals["main_actual_complex_RHS"]}；fork：{totals["fork_physical_LU"]} / {totals["fork_actual_complex_RHS"]}。source_RHS、adjoint_RHS及operation角色是aggregate RHS子项，不能再相加。', '',
        f'已记录独立validation物理总费：{totals["validation_physical_LU"]} LU / {totals["validation_actual_complex_RHS"]} aggregate complex RHS。不叠加其case/arm子计数；分区/常规模型费用分别见原回执。', '',
        'full_factorization_dimension是原记录的累计维度之和；逐run除以LU所得实际单因子维度另列为576，不将sum当单矩阵维度。denominator_small_LU与其dimension_sum/triangular RHS、FFT调用/变换、background H Cholesky/setup/iteration RHS各保留独立counter命名空间，不折算为Maxwell LU/RHS。峰值bytes/rank仅取max/snapshot，不相加。', '',
        'Fork仅按cumulative−inherited计新auxiliary/background费用；继承setup不再收费。rank缓存每个实际operator_group的计算次数只记一次，24个复用execution不当成24次SVD。integration恢复历史counter不重复计；无法从原回执拆开的pending replay新aux保留unknown，但其CPU与物理费用已在真实进程总费中。', '',
        'CPU_BY_CATEGORY_SCOPE按host/category/OS整进程、process snapshot和body分开。PROCESS_COSTS每个独立进程仅选一个首选CPU，内部body只作嵌套参考；不把case/arm/plot子计时重复叠加。Wall是不同且可能并行的区间，不求单个总和；不比较Mac与Windows绝对速度。scoring+plotting无法单独拆分时保留组合类。', '',
        '失败评分attempt01与首次打包外部531.625s（0 ZIP）都保留。V2原PACKAGE_RECEIPT存在时按原16part文件/bytes合计核对，body和process snapshot分开；尚未copy时保留pending。POST_MAIN_PIPELINE给出的fork控制器/review-export/失败打包外部OS CPU只替代同一进程内部计时，不新增重复费用。停止的bulk SCP、失败detached launch及未计量部署/上传CPU不填0。currency与LLM费用未知保持null；GPU_used=false只限已记录CPU程序。', '',
        '| Method | Tasks | OS CPU seconds | Maxwell LU | Aggregate complex RHS |', '|---|---:|---:|---:|---:|']
    lines += [f'| {r["method"]} | {r["actual_terminal_rows"]}/24 | {r["external_CPU_seconds"]:.6f} | {r["physical_LU"]} | {r["actual_complex_RHS"]} |' for r in methods]
    lines += ['', '## 首个已完成corrector审计（非精确起点）', '',
        '只取每任务HISTORY中phase=corrector的最小有限CPU_seconds，不使用raw best candidate CPU推测proposal结束。已完成审计包含前置计算和该次正演/记录；它不是进入corrector的精确时刻。120减该数仅为该已记录时刻后的名义预算余额，不是corrector实际总时长。proposal终止标签另从原RECEIPT派生的COST_LEDGER读取；缺失保持null，不由时间猜测。', '',
        '| Method | Audits /24 | Min | P25 | Median | Mean | P75 | Max |', '|---|---:|---:|---:|---:|---:|---:|---:|']
    def fmt(value):return 'null' if value is None else f'{value:.6f}'
    lines += ['| '+r['method']+' | '+str(r['with_completed_corrector_audit'])+'/24 | '+' | '.join(fmt(r[k]) for k in
        ['min_CPU_seconds','p25_CPU_seconds','median_CPU_seconds','mean_CPU_seconds','p75_CPU_seconds','max_CPU_seconds'])+' |' for r in corrector_summary]
    lines += ['', 'proposal termination counts：', '', '| Method | Original proposal termination | Count /24 |', '|---|---|---:|']
    lines += [f'| {r["method"]} | {json.dumps(r["proposal_termination"],ensure_ascii=False)} | {r["count"]}/24 |' for r in termination_counts]
    lines += ['', '原始/派生counter：COUNTERS_RECORDED、COUNTERS_BY_CATEGORY、逐11方法METHOD_COUNTERS11；进程/CPU范围：PROCESS_COSTS、CPU_BY_CATEGORY_SCOPE；缺口：GAPS；源核对：RECONCILIATION。全部CSV并有JSON副本。', '']
    (out/'COMPACT.md').write_text('\n'.join(lines),encoding='utf-8')
    receipt=dict(status='COMPLETE_WITH_EXPLICIT_GAPS',source_files=read_files,metadata_files_read=len(read_files),
        expected_main=264,expected_methods=11,expected_forks=56,failed_reconciliation=sum(x.get('matches') is False for x in checks),
        gaps=len(gaps),totals=totals,corrector_timing_task_rows=len(corrector_rows),
        corrector_timing_definition='Min completed HISTORY phase=corrector CPU, NOT exact start; no rawbest CPU surrogate',
        CPU_seconds=time.process_time()-start_cpu,whole_process_CPU_at_receipt_seconds=time.process_time(),
        wall_seconds=time.perf_counter()-start_wall,GT_arrays_read=0,other_arrays_read=0,new_physics=0,new_optimization=0,hash_operations='NONE')
    write(out/'RECEIPT.json',receipt)
    print(json.dumps({k:receipt[k] for k in ['status','metadata_files_read','failed_reconciliation','gaps','totals','CPU_seconds','wall_seconds']}))


if __name__=='__main__':main()
