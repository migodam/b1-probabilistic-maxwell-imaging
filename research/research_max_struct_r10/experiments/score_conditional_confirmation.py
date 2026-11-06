"""Prepared conditional4 scorer. Requires separate root postseal authorization.

No private file is opened by importing this module. No forward solver is used.
Campaign cost counts each actual prefix/fallback once; deployment views overlap.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,csv,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics,dump,csvwrite,stamp
DEPLOYMENTS=('CIE_baseline','CIE_then_DualCIE','CIE_then_G19_PGP')
LABELS=dict(CIE_baseline='CIE baseline',CIE_then_DualCIE='Conditional CIE → Dual',CIE_then_G19_PGP='Conditional CIE → G19 PGP')
def js(p):return json.loads(p.read_text())
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def strict_sum(values):
    values=list(values);return sum(values) if all(v is not None for v in values) else None
def result_path(task,kind):
    return ROOT/'RESULTS/conditional_confirmation'/kind/task['task_id']

def gate(args):
    if args.authorize_postseal!='ROOT_EXPLICIT_CONDITIONAL4_POSTSEAL_AUTHORIZATION' or not args.authorization_note.strip():
        raise PermissionError('Separate root conditional4 postseal authorization required; private truth unopened')
    freeze=js(args.freeze);pipeline=js(args.pipeline_receipt)
    if freeze['status']!='FROZEN_BEFORE_NEW_DATA_GENERATION_AND_GT_ACCESS' or freeze['seed']!=2026110711 or len(freeze['cases'])!=4:raise ValueError('Unexpected frozen conditional cohort')
    if pipeline['status']!='ALL_ACTUAL_CONDITIONAL_PREDICTIONS_SEALED' or pipeline.get('truth_access')!='NONE' or not pipeline['shared_prefix_charged_once']:
        raise RuntimeError('All actual inverse outcomes and deployment selections must be sealed')
    common=js(ROOT/'delegated/twod/FREEZE_SCORING.json')
    if (common['primary_support']['threshold'],common['primary_success']['nrmse_threshold'],common['primary_success']['fixed_support_iou_threshold'])!=(1.,.15,.7):raise ValueError('Common frozen score definitions changed')
    manifest=js(args.prefix_manifest);public=js(ROOT/'data_public/conditional_confirmation/manifest.json')
    cases=[r['case_id'] for r in freeze['cases']];prefix={r['case_id']:r for r in manifest['tasks']}
    if len(prefix)!=4 or set(prefix)!=set(cases) or {r['case_id'] for r in public['records']}!=set(cases):raise ValueError('Four public identities mismatch')
    cp=args.pipeline_receipt.parent;trigger_doc=js(cp/'TRIGGER_DECISIONS.json');triggers={r['case_id']:r for r in trigger_doc['rows']}
    if set(triggers)!=set(cases):raise RuntimeError('Incomplete public trigger decisions')
    fallback_doc=js(cp/'FALLBACK_TASKS.json');fallbacks=fallback_doc['tasks'];jobs={};evidence=[stamp(args.freeze),stamp(args.pipeline_receipt),stamp(args.prefix_manifest),stamp(cp/'TRIGGER_DECISIONS.json'),stamp(cp/'FALLBACK_TASKS.json')]
    process_jobs=[];startup_jobs=[]
    for q in pipeline['queues']:
        if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS'):raise RuntimeError('A planned queue remains nonterminal')
        for j in q['jobs']:
            if j['status'] not in ('FINISHED','FAILED','LAUNCH_FAILED'):raise RuntimeError('Nonterminal inverse child')
            process_jobs.append(dict(j,queue_name=q['name']))
            if j['task_id'] in jobs:
                old=jobs[j['task_id']]
                if (q['name']!='startup_repair_queue' or j['status']!='FINISHED' or old['status']=='FINISHED'
                    or old.get('physical_model_cost') or 'All R10 outputs must stay inside' not in old.get('error_log_tail','')
                    or startup_jobs or not pipeline.get('original_pipeline_receipt')):
                    raise RuntimeError('Duplicate task is not the sole documented zero-physics startup repair')
                startup_jobs.append(old)
            jobs[j['task_id']]=j
    expected={r['task_id'] for r in prefix.values()}|{r['task_id'] for r in fallbacks}
    if (set(jobs)!=expected or len(jobs)!=pipeline['actual_inverse_runs'] or len(jobs)>12
        or len(process_jobs)!=pipeline.get('actual_worker_process_attempts',len(jobs))):raise RuntimeError('Actual inverse/worker receipt coverage mismatch')
    if startup_jobs:
        old=js(cp/pipeline['original_pipeline_receipt']);evidence.append(stamp(cp/pipeline['original_pipeline_receipt']))
        if old['status']!='ALL_ACTUAL_CONDITIONAL_PREDICTIONS_SEALED':raise RuntimeError('Missing archived first pipeline terminal receipt')
        declared=pipeline.get('engineering_startup_failures',[])
        if len(declared)!=1 or declared[0]['task_id']!=startup_jobs[0]['task_id'] or declared[0]['GetProcessTimes_CPU_seconds']!=startup_jobs[0].get('child_process_CPU_seconds') or declared[0]['physical_LU']!=0 or declared[0]['physical_complex_RHS']!=0:
            raise RuntimeError('Declared engineering-startup reconciliation mismatch')
    pipeline['_scoring_process_jobs']=process_jobs;pipeline['_scoring_startup_jobs']=startup_jobs
    attempts={};selectors={}
    for case in cases:
        task=prefix[case];d=result_path(task,'prefix');rp=d/'RECEIPT.json';r=js(rp) if rp.exists() else {};g=triggers[case]
        obs=ROOT/'data_public/conditional_confirmation'/(case+'.npz');z=load(obs);evidence.append(stamp(obs))
        valid=r.get('status')=='SEALED_PUBLIC_PREDICTION' and (d/'prediction.npz').exists()
        if valid:
            pred=load(d/'prediction.npz');yp=pred['Y_prediction'];sigma=float(z['sigma']);radius=float(z['noise_radius'])
            if not np.isfinite(yp).all() or sigma<=0 or not np.isfinite(radius):raise RuntimeError('Invalid prefix trigger inputs')
            statistic=float(np.linalg.norm(yp-z['Y'])/sigma)
            if g['status']!='EVALUATED' or bool(statistic>radius)!=g['trigger'] or not np.isclose(statistic,g['whitened_port_residual'],rtol=1e-12,atol=1e-12):raise RuntimeError('Sealed trigger differs from frozen public discrepancy')
            if jobs[task['task_id']]['status']!='FINISHED':raise RuntimeError('Sealed prefix child has unresolved execution failure')
            evidence.extend(stamp(d/n) for n in ('prediction.npz','prediction.json','RECEIPT.json','PUBLIC_TASK.json'))
        elif g['status']!='INVALID_PREFIX_NO_COMPUTABLE_TRIGGER' or g['trigger']:raise RuntimeError('Invalid prefix incorrectly treated as a computable trigger')
        attempts[task['task_id']]=dict(case_id=case,kind='prefix',task=task,directory=d,receipt=r,job=jobs[task['task_id']],valid=valid,
            startup_jobs=[s for s in startup_jobs if s['task_id']==task['task_id']])
        expected_fb=[case+'__conditional_fallback_'+label for label in ('G19_PGP','DualCIE')] if g['trigger'] else []
        actual_fb=[t for t in fallbacks if t['case_id']==case]
        if sorted(t['task_id'] for t in actual_fb)!=sorted(expected_fb):raise RuntimeError('Fallbacks do not exactly follow shared strict trigger')
        for ft in actual_fb:
            fd=result_path(ft,'fallback');fr=js(fd/'RECEIPT.json') if (fd/'RECEIPT.json').exists() else {};fv=fr.get('status')=='SEALED_PUBLIC_PREDICTION' and (fd/'prediction.npz').exists()
            if fv:
                evidence.extend(stamp(fd/n) for n in ('prediction.npz','prediction.json','RECEIPT.json','PUBLIC_TASK.json'))
                if ft['task_id'].endswith('G19_PGP'):
                    setup=js(fd/'DENOMINATOR_SETUP.json')
                    if (setup['projection'],setup.get('actual_denominator_rank',setup.get('green_strong_rank')),setup['matrix_inner_maxiter'],setup['basis'])!=('galerkin',19,8,'green_strong'):raise RuntimeError('Conditional G19 numerical contract mismatch')
                    evidence.append(stamp(fd/'DENOMINATOR_SETUP.json'))
            elif jobs[ft['task_id']]['status']=='FINISHED':raise RuntimeError('Fallback marked finished without a prediction seal')
            attempts[ft['task_id']]=dict(case_id=case,kind='fallback',task=ft,directory=fd,receipt=fr,job=jobs[ft['task_id']],valid=fv)
        for label in ('G19_PGP','DualCIE'):
            sd=ROOT/'RESULTS/conditional_confirmation/deployed'/(case+'__CIE_then_'+label);sp=sd/'SELECTION.json';sel=js(sp)
            if sel['case_id']!=case or sel['method']!='CIE_then_'+label or sel['trigger']!=g['trigger'] or sel['deployment_cost_additive'] is not False:raise RuntimeError('Invalid deployment selector seal')
            order=[task['task_id']]+([case+'__conditional_fallback_'+label] if g['trigger'] else [])
            eligible=[attempts[t] for t in order if attempts[t]['valid'] and attempts[t]['receipt'].get('physical_objective') is not None]
            expected_choice=min(eligible,key=lambda a:a['receipt']['physical_objective'])['task']['task_id'] if eligible else None
            chosen=sel.get('selected');actual_choice=chosen['task_id'] if chosen else None
            if actual_choice!=expected_choice:raise RuntimeError('Selector is not min originalPhi with prefix tie priority')
            expected_prefix_cpu=strict_sum([jobs[task['task_id']].get('child_process_CPU_seconds')]+[s.get('child_process_CPU_seconds') for s in startup_jobs if s['task_id']==task['task_id']])
            if sel.get('prefix_full_child_CPU_seconds')!=expected_prefix_cpu:raise RuntimeError('Selector prefix cost does not preserve known startup cost')
            expected_startup_cpu=strict_sum(s.get('child_process_CPU_seconds') for s in startup_jobs if s['task_id']==task['task_id'])
            if sel.get('prefix_engineering_startup_CPU_seconds',0.)!=expected_startup_cpu:raise RuntimeError('Selector engineering-startup provenance mismatch')
            evidence.append(stamp(sp))
            if chosen:
                dp=sd/'prediction.npz';cpred=load(attempts[actual_choice]['directory']/'prediction.npz');dout=load(dp)
                for key in ('m','Y_prediction'):
                    if not np.array_equal(cpred[key],dout[key]):raise RuntimeError('Deployed prediction differs from publicly selected attempt')
                evidence.extend(stamp(sd/n) for n in ('prediction.npz','prediction.json'))
            selectors[case,'CIE_then_'+label]=dict(selection=sel,directory=sd,selected_task_id=actual_choice)
    return freeze,common,pipeline,cases,prefix,triggers,attempts,selectors,evidence

def attempt_cost(a):
    r=a['receipt'];j=a['job'];c=r.get('physical_model_cost') or j.get('physical_model_cost') or {};ops=(r.get('auxiliary_execution_receipt') or {}).get('operation_counts') or {}
    launched=j['status']!='LAUNCH_FAILED';matrix=a['task'].get('material_readout',{}).get('label')=='conditional_G19_PGP'
    startup_cpu=strict_sum(s.get('child_process_CPU_seconds') for s in a.get('startup_jobs',[]))
    return dict(task_id=a['task']['task_id'],case_id=a['case_id'],kind=a['kind'],method_label=a['task'].get('material_readout',{}).get('label',a['task']['arm']),
        status=r.get('status','LAUNCH_FAILED' if not launched else'MISSING_RECEIPT'),queue_status=j['status'],full_child_CPU_seconds=strict_sum([j.get('child_process_CPU_seconds'),startup_cpu]),CPU_scope=j.get('CPU_scope'),
        actual_solve_child_CPU_seconds=j.get('child_process_CPU_seconds'),preserved_zero_physics_startup_CPU_seconds=startup_cpu,
        physical_LU=c.get('LU',0 if not launched else None),physical_source_RHS=c.get('source_RHS',0 if not launched else None),physical_adjoint_RHS=c.get('adjoint_RHS',0 if not launched else None),
        physical_complex_RHS=c.get('actual_complex_RHS',0 if not launched else None),physical_LU_dimension=576,
        small_denominator_LU=ops.get('denominator_small_LU',None if matrix and launched else 0),small_denominator_LU_dimension_sum=ops.get('denominator_small_LU_dimension_sum',None if matrix and launched else 0),
        small_denominator_complex_RHS=ops.get('denominator_small_triangular_complex_RHS',None if matrix and launched else 0),
        nested_R10_CPU_seconds=r.get('r10_total_cpu_seconds'),nested_R8_CPU_seconds=r.get('cpu_seconds'),nested_auxiliary_CPU_seconds=r.get('auxiliary_cpu_seconds'),nested_corrector_CPU_seconds=(r.get('corrector') or {}).get('cpu_seconds'),
        corrector_stop=(r.get('corrector') or {}).get('stop'),corrector_success=(r.get('corrector') or {}).get('success'),corrector_iterations=(r.get('corrector') or {}).get('iterations'))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--authorize-postseal',required=True);ap.add_argument('--authorization-note',required=True)
    ap.add_argument('--freeze',type=Path,default=ROOT/'delegated/twod/conditional_confirmation/FREEZE.json')
    ap.add_argument('--pipeline-receipt',type=Path,default=ROOT/'delegated/twod/conditional_confirmation/RECEIPT.json')
    ap.add_argument('--prefix-manifest',type=Path,default=ROOT/'configs/twod_conditional_prefix_tasks.json')
    ap.add_argument('--truth-dir',type=Path,default=ROOT/'sealed_truth/conditional_confirmation');ap.add_argument('--output',type=Path,default=ROOT/'delegated/conditional_confirmation_score')
    args=ap.parse_args();start=time.perf_counter();cpu=time.process_time()
    freeze,common,pipeline,cases,prefix,triggers,attempts,selectors,evidence=gate(args)
    if args.output.exists():raise FileExistsError('Write-once conditional4 score directory already exists')
    out=args.output;out.mkdir(parents=True)
    dump(out/'ALL_ACTUAL_AND_SELECTOR_SEAL.json',dict(status='ALL_ACTUAL_TERMINAL_AND_EIGHT_SELECTOR_SEALS_VERIFIED',root_authorization=args.authorization_note,
        actual_inverse_runs=len(attempts),triggered_cases=sum(g['trigger'] for g in triggers.values()),truth_opened_by_this_process=False,evidence=evidence))
    # FIRST new conditional GT access, after separate root authorization and public gate.
    sm=js(args.truth_dir/'score_metadata.json');meta={r['case_id']:r for r in sm['records']}
    if set(meta)!=set(cases):raise ValueError('Conditional score metadata identities differ from frozen4')
    truth={}
    for case in cases:
        with np.load(args.truth_dir/(case+'.npz'),allow_pickle=False) as z:truth[case]=z['m'].copy()
    scored_attempts={};stage_rows=[];costs={tid:attempt_cost(a) for tid,a in attempts.items()};pictures={}
    for tid,a in attempts.items():
        case=a['case_id'];m=meta[case];r=a['receipt'];d=a['directory'];c=costs[tid]
        row=dict(task_id=tid,case_id=case,kind=a['kind'],method_label=c['method_label'],status=c['status'],material_success=False,amplitude_NRMSE=None,fixed_support_IoU=None,
            family=m['family'],high_chi=m['high_chi'],amplitude=m['amplitude'],original_Phi=r.get('physical_objective'),auxiliary_selected_stage=r.get('auxiliary_selected_stage'),corrector_stop=c['corrector_stop'],corrector_success=c['corrector_success'],corrector_iterations=c['corrector_iterations'])
        if a['valid']:
            pred=load(d/'prediction.npz');obs=load(ROOT/'data_public/conditional_confirmation'/(case+'.npz'));sigma=float(obs['sigma']);radius=float(obs['noise_radius']);we=float(np.linalg.norm(pred['Y_prediction']-obs['Y'])/sigma)
            row.update(material_metrics(pred['m'],truth[case],m['amplitude'],m['high_chi']),whitened_residual=we,CN99_radius=radius,CN99_raw_noise_pass=bool(we<=radius),
                relative_physical_residual=float(np.linalg.norm(pred['Y_prediction']-obs['Y'])/np.linalg.norm(obs['Y'])))
            for p in sorted(d.glob('stage_*_physical.npz'))+[d/'selected_before_corrector.npz']:
                if p.exists():
                    state=load(p);sr=dict(task_id=tid,case_id=case,state=p.stem,**material_metrics(state['m'],truth[case],m['amplitude'],m['high_chi']))
                    if p.with_suffix('.json').exists():sr['original_Phi']=js(p.with_suffix('.json')).get('objective')
                    if 'Y_prediction' in state:
                        sr['whitened_residual']=float(np.linalg.norm(state['Y_prediction']-obs['Y'])/sigma);sr['CN99_raw_noise_pass']=bool(sr['whitened_residual']<=radius)
                    stage_rows.append(sr)
        else:row['failure_reason']=r.get('error',a['job'].get('error',a['job'].get('error_log_tail','Explicit terminal attempt failure')))
        scored_attempts[tid]=row
    opportunity_rows=[]
    for case in cases:
        opportunity_rows.append(scored_attempts[prefix[case]['task_id']])
        for label in ('DualCIE','G19_PGP'):
            tid=case+'__conditional_fallback_'+label
            if tid in scored_attempts:opportunity_rows.append(scored_attempts[tid])
            else:opportunity_rows.append(dict(case_id=case,task_id=tid,kind='fallback',method_label=label,status='NOT_RUN',
                reason='INVALID_PREFIX_NO_COMPUTABLE_TRIGGER' if triggers[case]['status']!='EVALUATED' else'PREFIX_DID_NOT_STRICTLY_EXCEED_NOISE_RADIUS',material_success=None,amplitude_NRMSE=None,fixed_support_IoU=None))
    deployments=[]
    for case in cases:
        ptid=prefix[case]['task_id'];pc=costs[ptid];g=triggers[case]
        for method in DEPLOYMENTS:
            selected=ptid if attempts[ptid]['valid'] else None;fbid=None;selection='Baseline prefix only';fallback_status='NOT_APPLICABLE'
            if method!='CIE_baseline':
                sel=selectors[case,method];selected=sel['selected_task_id'];selection=sel['selection']['selection'];label=method.removeprefix('CIE_then_');fbid=case+'__conditional_fallback_'+label if g['trigger'] else None
                fallback_status=costs[fbid]['status'] if fbid else'NOT_RUN'
            fc=costs[fbid] if fbid else None
            row=dict(case_id=case,method=method,family=meta[case]['family'],trigger=g['trigger'],trigger_status=g['status'],fallback_status=fallback_status,
                selected_task_id=selected,selected_kind=attempts[selected]['kind'] if selected else None,selection=selection,status='VALID_SELECTED_PREDICTION' if selected else'NO_VALID_SELECTED_PREDICTION',
                material_success=False,amplitude_NRMSE=None,fixed_support_IoU=None,deployment_cost_additive=False,
                prefix_full_child_CPU_seconds=pc['full_child_CPU_seconds'],own_fallback_full_child_CPU_seconds=fc['full_child_CPU_seconds'] if fc else 0.,
                deployment_full_child_CPU_seconds=strict_sum([pc['full_child_CPU_seconds'],fc['full_child_CPU_seconds'] if fc else 0.]),
                deployment_physical_LU=strict_sum([pc['physical_LU'],fc['physical_LU'] if fc else 0]),deployment_physical_complex_RHS=strict_sum([pc['physical_complex_RHS'],fc['physical_complex_RHS'] if fc else 0]),
                deployment_physical_source_RHS=strict_sum([pc['physical_source_RHS'],fc['physical_source_RHS'] if fc else 0]),
                deployment_physical_adjoint_RHS=strict_sum([pc['physical_adjoint_RHS'],fc['physical_adjoint_RHS'] if fc else 0]),deployment_physical_LU_dimension=576,
                deployment_small_denominator_LU=fc['small_denominator_LU'] if fc else 0,
                deployment_small_denominator_LU_dimension_sum=fc['small_denominator_LU_dimension_sum'] if fc else 0,
                deployment_small_denominator_complex_RHS=fc['small_denominator_complex_RHS'] if fc else 0,
                cost_scope='Prefix plus own fallback whether failed/unselected; alternative deployment views share prefix and must not be summed as campaign bill')
            if selected:
                s=scored_attempts[selected]
                for k,v in s.items():
                    if k not in ('task_id','case_id','kind','method_label','status','family'):row[k]=v
                pictures[case,method]=load(attempts[selected]['directory']/'prediction.npz')['m']
            deployments.append(row)
    summaries=[]
    for method in DEPLOYMENTS:
        z=[r for r in deployments if r['method']==method];valid=[r for r in z if r['amplitude_NRMSE'] is not None]
        summaries.append(dict(method=method,objects=4,valid_predictions=len(valid),material_successes=sum(r['material_success'] for r in z),CN99_passes=sum(bool(r.get('CN99_raw_noise_pass')) for r in z),
            mean_amplitude_NRMSE=float(np.mean([r['amplitude_NRMSE'] for r in valid])) if valid else None,median_amplitude_NRMSE=float(np.median([r['amplitude_NRMSE'] for r in valid])) if valid else None,
            mean_fixed_support_IoU=float(np.mean([r['fixed_support_IoU'] for r in valid])) if valid else None,
            deployment_full_child_CPU_seconds=strict_sum(r['deployment_full_child_CPU_seconds'] for r in z),deployment_physical_LU=strict_sum(r['deployment_physical_LU'] for r in z),
            deployment_physical_complex_RHS=strict_sum(r['deployment_physical_complex_RHS'] for r in z),
            deployment_small_denominator_LU=strict_sum(r['deployment_small_denominator_LU'] for r in z),
            deployment_small_denominator_LU_dimension_sum=strict_sum(r['deployment_small_denominator_LU_dimension_sum'] for r in z),
            deployment_small_denominator_complex_RHS=strict_sum(r['deployment_small_denominator_complex_RHS'] for r in z),deployment_cost_additive=False))
    csvwrite(out/'DEPLOYMENT_VIEWS.csv',deployments);csvwrite(out/'ALL_ATTEMPTS_WITH_NOT_RUN.csv',opportunity_rows);csvwrite(out/'ACTUAL_ATTEMPT_COSTS.csv',list(costs.values()));csvwrite(out/'ACTUAL_ATTEMPT_STAGES.csv',stage_rows);csvwrite(out/'DEPLOYMENT_SUMMARY.csv',summaries)
    csvwrite(out/'RAW_WORKER_PROCESS_COSTS_NOT_ADDITIVE.csv',[dict(task_id=j['task_id'],case_id=j['case_id'],queue_name=j['queue_name'],status=j['status'],
        full_child_CPU_seconds=j.get('child_process_CPU_seconds'),CPU_scope=j.get('CPU_scope'),physical_model_cost=j.get('physical_model_cost'),
        role='ZERO_PHYSICS_STARTUP_FAILURE' if j['status']!='FINISHED' and not j.get('physical_model_cost') and 'All R10 outputs must stay inside' in j.get('error_log_tail','') else'ACTUAL_SOLVE_OR_OTHER_TERMINAL_ATTEMPT') for j in pipeline['_scoring_process_jobs']])
    generation=js(args.freeze.parent/'GENERATION_RECEIPT.json');gen_cpu=generation.get('CPU_seconds');inverse_cpu=strict_sum(c['full_child_CPU_seconds'] for c in costs.values())
    if inverse_cpu is not None and not np.isclose(inverse_cpu,pipeline['actual_whole_inverse_child_CPU_seconds'],rtol=0,atol=1e-9):raise RuntimeError('Campaign unique child cost does not reconcile with producer raw processes')
    queue_cpu=strict_sum(q['external_queue_process_cost'].get('child_process_CPU_seconds') for q in pipeline['queues']);controller=pipeline.get('pipeline_controller_CPU_seconds')
    campaign=dict(shared_prefix_counted_once=True,unique_inverse_attempts=len(costs),prefix_runs=4,triggered_cases=sum(g['trigger'] for g in triggers.values()),fallback_runs=len(costs)-4,
        actual_worker_processes=len(pipeline['_scoring_process_jobs']),zero_physics_failed_startups=len(pipeline['_scoring_startup_jobs']),
        startup_CPU_already_in_prefix_cost=True,raw_worker_process_table_is_alternative_accounting_view=True,
        actual_unique_inverse_full_child_CPU_seconds=inverse_cpu,queue_external_own_process_CPU_seconds=queue_cpu,pipeline_recorded_own_CPU_seconds=controller,
        generation_CPU_seconds_charged_once=gen_cpu,known_campaign_CPU_sum=strict_sum([inverse_cpu,queue_cpu,controller,gen_cpu]),
        actual_inverse_physical_LU=strict_sum(c['physical_LU'] for c in costs.values()),actual_inverse_physical_complex_RHS=strict_sum(c['physical_complex_RHS'] for c in costs.values()),
        generation_physical_cost_separate=generation.get('physical_cost'),actual_small_denominator_LU=strict_sum(c['small_denominator_LU'] for c in costs.values()),
        generation_physical_factorization_dimension_per_call=576,generation_dimension_counter_note='Generator aggregate full_factorization_dimension sums four576 dimensions; it is not a2304 dimensional LU',
        actual_small_denominator_dimension_sum=strict_sum(c['small_denominator_LU_dimension_sum'] for c in costs.values()),actual_small_denominator_complex_RHS=strict_sum(c['small_denominator_complex_RHS'] for c in costs.values()),
        deployment_views_additive=False,nested_internal_timers_not_added=True,
        exclusions='Transport/SSH/packaging/archive and unmeasured pipeline lifetime tail are separately itemized engineering costs; no invented allocation to deployment views',
        physical_counter_scope='All actual inverse576LU/RHS once, including failed/unselected fallbacks; generation separate;19-dimensional systems never merged with576LU')
    dump(out/'CAMPAIGN_COSTS.json',campaign);dump(out/'TRIGGER_DECISIONS_USED.json',triggers);dump(out/'DEPLOYMENT_SUMMARY.json',summaries)
    render(out,cases,truth,pictures,deployments)
    review=['# Frozen conditional4 independent scoring','',
        'All four prespecified high objects are retained. Separate root postseal authorization and all actual terminal/selector seals precede the first conditional truth read. No new forward or inverse calculation was made.','',
        '| Deployment | Valid | Material success | CN99 pass | Mean NRMSE | Full child CPU s |',
        '|---|---:|---:|---:|---:|---:|']
    for s in summaries:review.append('| %s | %s/4 | %s/4 | %s/4 | %s | %s |'%(LABELS[s['method']],s['valid_predictions'],s['material_successes'],s['CN99_passes'],s['mean_amplitude_NRMSE'],s['deployment_full_child_CPU_seconds']))
    review+=['',f"Strict public-noise triggers: {campaign['triggered_cases']}/4; actual fallback attempts: {campaign['fallback_runs']}. Untriggered fallback opportunities remain NOT_RUN, not successes or failures.",'',
        '| Case | Deployment | Selected kind | NRMSE | IoU | Material pass | CN99 pass | Original Phi |',
        '|---|---|---|---:|---:|---|---|---:|']
    for r in deployments:review.append('| %s | %s | %s | %s | %s | %s | %s | %s |'%(r['case_id'],LABELS[r['method']],r['selected_kind'],r.get('amplitude_NRMSE'),r.get('fixed_support_IoU'),r.get('material_success'),r.get('CN99_raw_noise_pass'),r.get('original_Phi')))
    review+=['','The three deployment CPU columns overlap the same CIE prefixes and must not be added. CAMPAIGN_COSTS.json charges each unique actual child, all failed/unselected work, own queue/controller CPU and generation once. The known zero-physics startup rejection is included once in its prefix cost and separately visible in the raw worker-process view; these two accounting tables are not additive. Nested inverse timers are retained as provenance, not additional CPU. Small19 LU and physical576 LU remain separate.','',
        'This is a four-object independent conditional confirmation under the frozen trigger and minimum-original-Phi selector. Family and high contrast are confounded; no broad population or contrast-causal claim follows. Material failures, legal noise-compatible bias and bad data fits must be distinguished. A lower original Phi does not certify material recovery. Main24 and exposed right3 scores are unchanged.','',
        'Files: DEPLOYMENT_VIEWS.csv, ALL_ATTEMPTS_WITH_NOT_RUN.csv, ACTUAL_ATTEMPT_STAGES.csv, ACTUAL_ATTEMPT_COSTS.csv, RAW_WORKER_PROCESS_COSTS_NOT_ADDITIVE.csv, CAMPAIGN_COSTS.json and ALL4_DEPLOYMENT_MATERIALS.png/.pdf.']
    (out/'REVIEW.md').write_text('\n'.join(review)+'\n')
    dump(out/'RECEIPT.json',dict(status='COMPLETE',scope='INDEPENDENT_FROZEN_CONDITIONAL4_CONFIRMATION',root_authorization=args.authorization_note,
        objects=4,actual_inverse_attempts=len(attempts),triggered_cases=campaign['triggered_cases'],fallback_attempts=campaign['fallback_runs'],
        new_inverse_runs_by_scorer=0,new_forward_calls=0,new_physical_LU=0,new_physical_RHS=0,conditional_GT_access='POST_ALL_ACTUAL_TERMINAL_AND_EIGHT_SELECTOR_SEALS_ONLY',
        main24_modified=False,right3_modified=False,GT_used_for_trigger_or_selection=False,
        scope_limits='Four prespecified independent high objects; family/contrast confounded; no broad-population guarantee; no extra hard-case generation if no trigger',
        deployment_cost_warning='Three deployment columns overlap shared prefixes; charge campaign unique attempts once, never sum deployment columns',
        scoring_CPU_seconds=time.process_time()-cpu,scoring_wall_seconds=time.perf_counter()-start))
    print(json.dumps(dict(status='CONDITIONAL4_SCORING_COMPLETE',output=str(out),actual_inverse_attempts=len(attempts),triggered=campaign['triggered_cases'])))

def render(out,cases,truth,pictures,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(4,4,figsize=(11,10.5),constrained_layout=True)
    for i,case in enumerate(cases):
        for k,method in enumerate(('GT',)+DEPLOYMENTS):
            a=axs[i,k];m=truth[case] if method=='GT' else pictures.get((case,method));a.set_xticks([]);a.set_yticks([])
            if m is not None:im=a.imshow(m.reshape(24,24),origin='lower',cmap='viridis',vmin=.02,vmax=20)
            else:a.text(.5,.5,'NO VALID PREDICTION',ha='center',va='center',transform=a.transAxes,color='crimson')
            title=case.removeprefix('r10_conditional_')+' GT' if method=='GT' else LABELS[method]
            if method!='GT' and m is not None:
                r=next(r for r in rows if r['case_id']==case and r['method']==method);title+='\nN %.4f | IoU %.3f'%(r['amplitude_NRMSE'],r['fixed_support_IoU'])
            a.set_title(title,fontsize=8)
    fig.colorbar(im,ax=axs,shrink=.7,label='real contrast chi; common [.02,20]');fig.suptitle('Frozen conditional4 | all objects | shared CIE prefix | no GT trigger or selector',fontsize=11)
    fig.savefig(out/'ALL4_DEPLOYMENT_MATERIALS.png',dpi=180);fig.savefig(out/'ALL4_DEPLOYMENT_MATERIALS.pdf');plt.close(fig)

if __name__=='__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):main()
