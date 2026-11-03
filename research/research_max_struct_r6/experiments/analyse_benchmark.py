"""Offline R6 score/cost analysis, with no solver, prediction or GT imports.

Both completed independent SCORES.json files must exist before any analysis.
Only score JSON, original task/queue receipts, history JSON and stage sidecar
JSON are read. All outputs are new, write-once analysis/figure directories.
No scientific acceptance, novelty, convergence or superiority is decided here.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import argparse
import csv
import json
import math
import statistics
import time

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
ARMS=('cie_single','cie_dual','global_dual','a_material')
STRATA=('DEV9_DIGITS','DEV6_DIFFICULT','NEW6_CONFIRMATION')
TOTAL_FIELDS=('cpu_seconds','wall_seconds','LU','actual_complex_RHS','source_RHS','adjoint_RHS',
              'setup_cpu_seconds','setup_wall_seconds','auxiliary_cpu_seconds','auxiliary_wall_seconds',
              'corrector_cpu_seconds','corrector_wall_seconds','aux_plus_corrector_cpu_seconds',
              'aux_plus_corrector_wall_seconds','other_cpu_seconds','other_wall_seconds')
METRICS=('material_nrmse','otsu_iou','relative_residual','original_Phi','projected_gradient_inf')
NOTES=[
    'Each unique manifest task receipt is counted once. Scoring rows are not a second cost source.',
    'Preserved additional failed attempts are charged separately and never add material objects or final benchmark tasks. Outer queue process-wall measurements overlap each attempt receipt wall and are not added to it.',
    'Task CPU/wall totals are the original receipt scope, not OS-measured whole-process usage.',
    'Auxiliary and corrector timers are successive disjoint solver phases. Their sum excludes setup, interphase serialization and other receipt overhead.',
    'Other time = task total - auxiliary - corrector; it is not an independently measured timer and may be null or show small negative roundoff/coverage differences.',
    'Setup is inside task total and other time. Physics sub-operation timers and execution-receipt timers are nested; never add them to task or phase timers.',
    'Physical LU/RHS totals come only from physical_model_cost. Auxiliary operation counts and corrector evaluations are different quantities and are never added to LU/RHS.',
    'Four-worker queue elapsed is elapsed batch time; summed task wall is overlapping process time and is not batch elapsed. Their ratio is descriptive overlap only.',
    'Primary unique material object counts are independent-object denominators. Repeated arms, source columns and FFT duplicates add no independent material objects.',
    'Arm/stratum object counts are overlapping group denominators, not additive independent-object counts.',
    'All failed tasks stay in manifest denominators. Missing metrics/timers remain null and are not imputed as zeros, losses or successes.',
    'Material NRMSE comes only from the independent final score file. No intermediate material truth metrics are calculated.',
    'Original-Phi curves contain only already charged host/stage physical validations and corrector evaluations. Auxiliary per-iteration proxy objectives are not original Phi.',
    'Curve CPU is directly recorded: stage sidecars use auxiliary-start CPU; corrector evaluation history uses observer-creation CPU. The short intervening clock-origin offset was not logged, is not inferred, and is flagged in curve accounting.',
    'Sparse auxiliary physical endpoints do not establish a physical objective at every intermediate auxiliary update. Lines/steps only connect observed incumbents.',
    'No scientific conclusion or claim acceptance is supplied by this arithmetic/plotting script.'
]


def load(path):
    return json.loads(Path(path).read_text())


def resolve(path):
    path=Path(path)
    return path if path.is_absolute() else REPO/path


def numeric(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def stats(values):
    vals=[float(v) for v in values if numeric(v)]
    return dict(known=len(vals),missing=len(values)-len(vals),
                known_sum=math.fsum(vals) if vals else None,
                mean=statistics.mean(vals) if vals else None,
                median=statistics.median(vals) if vals else None,
                minimum=min(vals) if vals else None,maximum=max(vals) if vals else None)


def add_known(*values):
    return math.fsum(values) if all(numeric(v) for v in values) else None


def subtract(a,b):
    return a-b if numeric(a) and numeric(b) else None


def ratio(a,b):
    return a/b if numeric(a) and numeric(b) and b!=0 else None


def write_json(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def write_csv(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('x',newline='') as f:
        out=csv.DictWriter(f,fieldnames=fields);out.writeheader()
        for row in rows:
            out.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v
                          for k,v in row.items()})


def summarize_operations(rows,key):
    keys=sorted(set(k for row in rows for k in (row.get(key) or {})))
    return {k:stats([row[key].get(k,0) if isinstance(row.get(key),dict) else None for row in rows]) for k in keys}


def summarize(rows,batch,stratum,arm):
    ids=sorted({r['case_id'] for r in rows})
    result=dict(batch=batch,stratum=stratum,arm=arm,tasks=len(rows),unique_objects=len(ids),
        case_ids=ids,independent_material_objects_within_group=len(ids) if batch=='primary' else 0,
        sealed_predictions=sum(r['status']=='SEALED_PUBLIC_PREDICTION' for r in rows),
        statuses={s:sum(r['status']==s for r in rows) for s in sorted({r['status'] for r in rows})},
        failed_attempt_entries=sum(r['failed_attempt_count'] for r in rows),
        legacy_joint=dict(pass_count=sum(r.get('legacy_joint') is True for r in rows),
                          known=sum(isinstance(r.get('legacy_joint'),bool) for r in rows)),
        confirmation_joint=dict(pass_count=sum(r.get('confirmation_joint') is True for r in rows),
                                known=sum(isinstance(r.get('confirmation_joint'),bool) for r in rows)))
    result['statistics']={key:stats([r.get(key) for r in rows]) for key in TOTAL_FIELDS+METRICS}
    for key in ('physical_model_cost','auxiliary_physical_model_cost','corrector_physical_model_cost',
                'physical_unassigned_to_aux_or_corrector','auxiliary_operation_counts',
                'physical_operation_cpu','physical_operation_wall'):
        result[key]=summarize_operations(rows,key)
    return result


def primary_inputs():
    return dict(primary=ROOT/'RESULTS/scoring_primary/SCORES.json',
                fft=ROOT/'RESULTS/scoring_fft/SCORES.json')


def collect_batch(batch,score_path,manifest_path,queue_paths,seen_receipts):
    scores=load(score_path);manifest=load(manifest_path)
    rows=scores.get('rows')
    if not isinstance(rows,list):raise ValueError('Score file must contain rows: '+str(score_path))
    if scores.get('batch')!=batch:raise ValueError('Score batch mismatch: '+str(score_path))
    tasks=manifest['tasks'];index={r['task_id']:r for r in rows}
    if len(index)!=len(rows):raise ValueError('Duplicate score task ID in '+batch)
    if len({t['task_id'] for t in tasks})!=len(tasks):raise ValueError('Duplicate manifest task ID')
    if set(index)!={t['task_id'] for t in tasks}:raise ValueError('Score/manifest task coverage differs')
    receipts={};data=[];issues=[]
    for task in tasks:
        ident=task['task_id'];score=index[ident]
        if any(task[k]!=score[k] for k in ('case_id','arm')):raise ValueError('Score/manifest identity mismatch')
        directory=resolve(task['output_dir']);path=(directory/'RECEIPT.json').resolve()
        if path in seen_receipts:raise ValueError('A task receipt would be counted twice: '+str(path))
        seen_receipts.add(path);rec=load(path);receipts[ident]=rec
        if rec.get('status')!=score.get('status'):raise ValueError('Scoring/receipt status mismatch: '+ident)
        if rec.get('task_id') not in (None,ident):raise ValueError('Receipt task ID mismatch: '+ident)
        physical=rec.get('physical_model_cost')
        physical=physical if isinstance(physical,dict) else None
        auxiliary=rec.get('auxiliary_physical_model_cost')
        corrector=rec.get('corrector') or {};corrector_cost=corrector.get('cost')
        aux_cpu=rec.get('auxiliary_cpu_seconds');aux_wall=rec.get('auxiliary_wall_seconds')
        cor_cpu=corrector.get('cpu_seconds');cor_wall=corrector.get('wall_seconds')
        phase_cpu=add_known(aux_cpu,cor_cpu);phase_wall=add_known(aux_wall,cor_wall)
        row={key:value for key,value in score.items()
             if key not in ('auxiliary_execution_receipt','corrector_raw_summary')}
        row.update(batch=batch,receipt_path=str(path),output_dir=str(directory),
            receipt_cpu_scope=rec.get('cost_scope'),cpu_seconds=rec.get('cpu_seconds'),
            wall_seconds=rec.get('wall_seconds'),setup_cpu_seconds=rec.get('setup_cpu_seconds'),
            setup_wall_seconds=rec.get('setup_wall_seconds'),auxiliary_cpu_seconds=aux_cpu,
            auxiliary_wall_seconds=aux_wall,corrector_cpu_seconds=cor_cpu,corrector_wall_seconds=cor_wall,
            aux_plus_corrector_cpu_seconds=phase_cpu,aux_plus_corrector_wall_seconds=phase_wall,
            other_cpu_seconds=subtract(rec.get('cpu_seconds'),phase_cpu),
            other_wall_seconds=subtract(rec.get('wall_seconds'),phase_wall),
            physical_model_cost=physical,auxiliary_physical_model_cost=auxiliary,
            corrector_physical_model_cost=corrector_cost,
            auxiliary_operation_counts=(rec.get('auxiliary_execution_receipt') or {}).get('operation_counts'),
            physical_operation_cpu=rec.get('physical_operation_cpu'),
            physical_operation_wall=rec.get('physical_operation_wall'),
            failed_attempt_count=len(rec.get('failed_attempts') or []),
            process_peak_rss_bytes=rec.get('process_peak_rss_bytes'))
        for field in ('LU','actual_complex_RHS','source_RHS','adjoint_RHS'):
            row[field]=physical.get(field,0) if physical is not None else None
            if numeric(score.get(field)) and numeric(row[field]) and score[field]!=row[field]:
                issues.append(dict(task_id=ident,issue='Scoring/receipt cost mismatch',field=field,
                                   score=score[field],receipt=row[field]))
        # A missing operation key in a present physical Counter is an observed zero.
        if all(isinstance(x,dict) for x in (physical,auxiliary,corrector_cost)):
            row['physical_unassigned_to_aux_or_corrector']={
                k:physical.get(k,0)-auxiliary.get(k,0)-corrector_cost.get(k,0)
                for k in sorted(set(physical)|set(auxiliary)|set(corrector_cost))}
        else:row['physical_unassigned_to_aux_or_corrector']=None
        for key in ('other_cpu_seconds','other_wall_seconds'):
            if numeric(row[key]) and row[key]<-1e-6:
                issues.append(dict(task_id=ident,issue='Negative unassigned timing remainder',field=key,value=row[key]))
        data.append(row)
    queues=[];all_queue_ids=[]
    for queue_path in queue_paths:
        if not queue_path.is_file():
            issues.append(dict(issue='Queue segment receipt unavailable',path=str(queue_path)))
            continue
        queue=load(queue_path)
        qids=[r['task_id'] for r in queue.get('tasks',[])];all_queue_ids.extend(qids)
        if len(qids)!=len(set(qids)):
            issues.append(dict(issue='Duplicate task IDs inside one queue segment; original task receipt costs are not re-summed',path=str(queue_path)))
        queues.append(dict(path=str(queue_path),status=queue.get('status'),workers=queue.get('workers'),
            started_utc=queue.get('started_utc'),finished_utc=queue.get('finished_utc'),
            batch_wall_seconds=queue.get('batch_wall_seconds'),
            tasks=len(qids),task_ids=qids,failed_attempts=queue.get('failed_attempts') or [],
            process_wall_seconds=stats([r.get('process_wall_seconds') for r in queue.get('tasks',[])])))
    if set(all_queue_ids)!={r['task_id'] for r in tasks}:
        issues.append(dict(issue='Combined queue task-ID coverage differs from final manifest'))
    overlap=[]
    for left,right in zip(queues,queues[1:]):
        if left['finished_utc'] and right['started_utc']:
            end=datetime.fromisoformat(left['finished_utc']);start=datetime.fromisoformat(right['started_utc'])
            overlap.append(dict(left=left['path'],right=right['path'],gap_seconds=(start-end).total_seconds(),
                                sequential=start>=end))
    if any(not pair['sequential'] for pair in overlap):
        issues.append(dict(issue='Queue segments overlap; elapsed segment sum is not elapsed union'))
    elapsed=add_known(*(q['batch_wall_seconds'] for q in queues)) if queues else None
    batch_summary=dict(batch=batch,score_path=str(score_path),manifest_path=str(manifest_path),
        total_tasks=len(data),unique_material_objects=len({r['case_id'] for r in data}),
        additional_independent_objects=len({r['case_id'] for r in data}) if batch=='primary' else 0,
        classification='PRIMARY_SCIENTIFIC_COMPARISON' if batch=='primary' else 'ENGINEERING_DUPLICATE_OBSERVATIONS',
        queue_segments=queues,queue_segment_order_check=overlap,
        workers=sorted({q['workers'] for q in queues if q['workers'] is not None}),
        batch_elapsed_wall_seconds=elapsed,
        batch_elapsed_scope='Sum of recorded sequential segment elapsed times; gaps between segments excluded. Order check preserved.',
        task_wall=stats([r.get('wall_seconds') for r in data]),
        task_cpu=stats([r.get('cpu_seconds') for r in data]),
        physical_LU=stats([r.get('LU') for r in data]),
        physical_actual_complex_RHS=stats([r.get('actual_complex_RHS') for r in data]),
        queue_failed_attempt_entries=sum(len(q['failed_attempts']) for q in queues) if queues else None,
        queue_repeat_task_ids=sorted({ident for ident in all_queue_ids if all_queue_ids.count(ident)>1}),
        accounting_issues=issues)
    batch_summary['sum_task_wall_over_batch_elapsed']=ratio(batch_summary['task_wall']['known_sum'],batch_summary['batch_elapsed_wall_seconds'])
    return data,receipts,batch_summary


def collect_extra_attempts(path,seen_receipts):
    """Account preserved failed attempts separately from final task denominators."""
    if path is None:return [],[]
    entries=load(path)
    if isinstance(entries,dict):entries=entries.get('attempts',[])
    rows=[];issues=[]
    for entry in entries:
        receipt_path=resolve(entry['receipt_path']).resolve()
        if receipt_path in seen_receipts:
            raise ValueError('Extra attempt duplicates an already counted receipt: '+str(receipt_path))
        seen_receipts.add(receipt_path);rec=load(receipt_path)
        physical=rec.get('physical_model_cost');corrector=rec.get('corrector') or {}
        zero_evidence=None
        if entry.get('known_physical_zero_before_model'):
            evidence_path=resolve(entry['repair_evidence_path']);evidence=load(evidence_path)
            if not (physical is None and rec.get('status')=='FAILED'
                    and 'public=load_public(files)' in rec.get('traceback','')
                    and "Conflicting public value: material_bounds" in rec.get('traceback','')
                    and evidence.get('failed_task')==rec.get('task_id')
                    and evidence.get('failed_LU')==0 and evidence.get('failed_RHS')==0
                    and not rec.get('auxiliary_execution_receipt')):
                raise ValueError('Specific pre-model zero-cost evidence no longer matches preserved failed receipt')
            physical=dict(LU=0,actual_complex_RHS=0,source_RHS=0,adjoint_RHS=0)
            zero_evidence=dict(status='KNOWN_ZERO_BEFORE_MODEL_CONSTRUCTION',
                original_receipt_physical_model_cost=rec.get('physical_model_cost'),
                repair_evidence_path=str(evidence_path),classification=entry.get('zero_classification'),
                outer_queue_process_wall_seconds=evidence.get('failed_process_wall_seconds'))
        row=dict(batch=entry['batch'],task_id=entry.get('task_id',rec.get('task_id')),
            case_id=rec.get('case_id'),arm=rec.get('arm'),status=rec.get('status'),
            receipt_path=str(receipt_path),note=entry.get('note'),
            denominator_role='Additional charged attempt, not an additional material object or final task',
            cpu_seconds=rec.get('cpu_seconds'),wall_seconds=rec.get('wall_seconds'),
            auxiliary_cpu_seconds=rec.get('auxiliary_cpu_seconds'),
            auxiliary_wall_seconds=rec.get('auxiliary_wall_seconds'),
            corrector_cpu_seconds=corrector.get('cpu_seconds'),corrector_wall_seconds=corrector.get('wall_seconds'),
            physical_model_cost=physical,physical_zero_evidence=zero_evidence,error=rec.get('error'))
        for field in ('LU','actual_complex_RHS','source_RHS','adjoint_RHS'):
            row[field]=physical.get(field,0) if isinstance(physical,dict) else None
        if zero_evidence is not None:
            for key in ('auxiliary_cpu_seconds','auxiliary_wall_seconds','corrector_cpu_seconds','corrector_wall_seconds',
                        'aux_plus_corrector_cpu_seconds','aux_plus_corrector_wall_seconds'):
                row[key]=0.
            row['other_cpu_seconds']=row['cpu_seconds'];row['other_wall_seconds']=row['wall_seconds']
            row['phase_timing_zero_basis']='Both solver phases were not entered; specifically verified pre-model load failure'
        if physical is None:
            issues.append(dict(receipt_path=str(receipt_path),issue='Failed-attempt physical model cost unavailable; not imputed to zero'))
        rows.append(row)
    return rows,issues


def charged_cost_groups(rows,attempts):
    """Final task cost plus preserved extras; object denominators stay unique."""
    lookup={(r['batch'],r['case_id']):r for r in rows}
    for extra in attempts:
        base=lookup.get((extra['batch'],extra['case_id']))
        extra.update(stratum=base.get('stratum') if base else 'UNKNOWN',label=base.get('label') if base else None)
    result=[]
    for batch in ('primary','fft'):
        arms=[a for a in ARMS if any(r['batch']==batch and r['arm']==a for r in rows+attempts)]
        for stratum in ('ALL',)+STRATA:
            for arm in ('ALL',)+tuple(arms):
                def matches(r):
                    return r['batch']==batch and (stratum=='ALL' or r['stratum']==stratum) and (arm=='ALL' or r['arm']==arm)
                final=[r for r in rows if matches(r)];extra=[r for r in attempts if matches(r)]
                if not final and not extra:continue
                group=final+extra
                result.append(dict(batch=batch,stratum=stratum,arm=arm,final_tasks=len(final),
                    additional_attempts=len(extra),charged_attempts=len(group),
                    unique_material_objects=len({r['case_id'] for r in group}),
                    additional_independent_objects_from_extra_attempts=0,
                    statistics={k:stats([r.get(k) for r in group]) for k in TOTAL_FIELDS},
                    physical_model_cost=summarize_operations(group,'physical_model_cost')))
    return result


def pairwise(rows):
    pairs=[];summaries=[]
    for batch in ('primary','fft'):
        groups={}
        for row in rows:
            if row['batch']==batch:groups.setdefault(row['case_id'],{})[row['arm']]=row
        available={r['arm'] for r in rows if r['batch']==batch}
        for baseline in ('cie_dual','global_dual'):
            if baseline not in available:continue
            for ident,group in groups.items():
                a,b=group.get('a_material'),group.get(baseline)
                if a is None or b is None:continue
                row=dict(batch=batch,case_id=ident,stratum=a['stratum'],label=a['label'],
                    baseline=baseline,a_status=a['status'],baseline_status=b['status'])
                for metric in ('material_nrmse','original_Phi','relative_residual','cpu_seconds','wall_seconds','LU','actual_complex_RHS'):
                    delta=subtract(b.get(metric),a.get(metric))
                    row[metric+'_baseline_minus_A']=delta
                    row[metric+'_relative_reduction']=ratio(delta,b.get(metric))
                row['otsu_iou_A_minus_baseline']=subtract(a.get('otsu_iou'),b.get('otsu_iou'))
                row['cpu_baseline_over_A']=ratio(b.get('cpu_seconds'),a.get('cpu_seconds'))
                row['nrmse_A_lower']=(a['material_nrmse']<b['material_nrmse'] if numeric(a.get('material_nrmse')) and numeric(b.get('material_nrmse')) else None)
                gate='confirmation_joint' if a['stratum']=='NEW6_CONFIRMATION' else 'legacy_joint'
                row.update(joint_field=gate,a_joint=a.get(gate),baseline_joint=b.get(gate))
                pairs.append(row)
            selected=[r for r in pairs if r['batch']==batch and r['baseline']==baseline]
            for stratum in ('ALL',)+STRATA:
                group=[r for r in selected if stratum=='ALL' or r['stratum']==stratum]
                if not group:continue
                vals=[r['material_nrmse_baseline_minus_A'] for r in group]
                summaries.append(dict(batch=batch,baseline=baseline,stratum=stratum,
                    paired_objects=len(group),known_material_pairs=sum(numeric(v) for v in vals),
                    A_lower_nrmse=sum(numeric(v) and v>0 for v in vals),
                    A_higher_nrmse=sum(numeric(v) and v<0 for v in vals),
                    equal_nrmse=sum(numeric(v) and v==0 for v in vals),
                    missing_material_pairs=sum(not numeric(v) for v in vals),
                    A_joint_passes=sum(r['a_joint'] is True for r in group),
                    baseline_joint_passes=sum(r['baseline_joint'] is True for r in group),
                    known_joint_pairs=sum(isinstance(r['a_joint'],bool) and isinstance(r['baseline_joint'],bool) for r in group),
                    statistics={key:stats([r.get(key) for r in group]) for key in group[0]
                                if key.endswith(('_baseline_minus_A','_A_minus_baseline','_relative_reduction'))}))
    return pairs,summaries


def backend_pairs(rows):
    primary={(r['case_id'],r['arm']):r for r in rows if r['batch']=='primary'}
    result=[]
    for fft in (r for r in rows if r['batch']=='fft'):
        dense=primary.get((fft['case_id'],fft['arm']))
        if dense is None:continue
        row=dict(case_id=fft['case_id'],stratum=fft['stratum'],label=fft['label'],arm=fft['arm'],
                 independent_objects_added=0,dense_task_id=dense['task_id'],fft_task_id=fft['task_id'])
        for field in ('material_nrmse','original_Phi','relative_residual','cpu_seconds','wall_seconds',
                      'auxiliary_cpu_seconds','auxiliary_wall_seconds','corrector_cpu_seconds','LU','actual_complex_RHS'):
            row[field+'_dense']=dense.get(field);row[field+'_fft']=fft.get(field)
            row[field+'_dense_minus_fft']=subtract(dense.get(field),fft.get(field))
            row[field+'_dense_over_fft']=ratio(dense.get(field),fft.get(field))
        result.append(row)
    return result


def collect_curves(rows,receipts,stratum):
    selected=[r for r in rows if r['batch']=='primary' and r['stratum']==stratum]
    case_ids=sorted({r['case_id'] for r in selected})
    if len(case_ids)!=6:raise ValueError('Fixed curve stratum must contain all six objects')
    points=[];accounting=[]
    for row in selected:
        directory=Path(row['output_dir']);rec=receipts[row['task_id']]
        local=[];issues=[];history_path=directory/'HISTORY.json';corrector_path=directory/'CORRECTOR.json'
        # History is used for observed-cost consistency, never for surrogate Phi.
        history=load(history_path) if history_path.is_file() else []
        if not history:issues.append('Auxiliary history unavailable')
        stage_ends={r['stage']:r for r in history if r.get('event')=='stage_end'}
        incumbent=None
        for candidate in rec.get('auxiliary_candidate_scores') or []:
            stage=candidate['stage'];phi=candidate.get('objective')
            if numeric(phi):incumbent=phi if incumbent is None else min(incumbent,phi)
            sidecar=directory/f'stage_{stage:02d}_joint.json' if stage>=0 else None
            meta=load(sidecar) if sidecar is not None and sidecar.is_file() else {}
            cost=meta.get('physical_model_cost') or {}
            if stage>=0 and not meta:issues.append('Stage sidecar unavailable: '+str(stage))
            local.append(dict(batch='primary',task_id=row['task_id'],case_id=row['case_id'],arm=row['arm'],
                phase='initial' if stage<0 else 'auxiliary_stage',stage=stage,
                global_iteration=None if stage<0 else sum(1 for h in history if 'iteration' in h and h['stage']<=stage),
                evaluation=None,observed_Phi=phi,best_original_Phi=incumbent,
                cumulative_LU=cost.get('LU'),cumulative_actual_complex_RHS=cost.get('actual_complex_RHS'),
                cumulative_observed_cpu_seconds=meta.get('cpu_seconds'),
                cpu_origin='UNKNOWN_INITIAL_VALIDATION_TIME' if stage<0 else 'AUXILIARY_START',
                stop=stage_ends.get(stage,{}).get('status')))
        corrector=load(corrector_path) if corrector_path.is_file() else {}
        for evaluation in corrector.get('evaluation_history') or []:
            cost=evaluation.get('physical_model_cost') or {}
            local.append(dict(batch='primary',task_id=row['task_id'],case_id=row['case_id'],arm=row['arm'],
                phase='corrector_evaluation',stage=None,global_iteration=None,evaluation=evaluation['evaluation'],
                observed_Phi=evaluation.get('objective'),best_original_Phi=evaluation.get('best_objective'),
                cumulative_LU=cost.get('LU'),cumulative_actual_complex_RHS=cost.get('actual_complex_RHS'),
                cumulative_observed_cpu_seconds=evaluation.get('cpu_seconds'),cpu_origin='OBSERVER_CREATION',stop=None))
        if row['status']=='SEALED_PUBLIC_PREDICTION' and not corrector.get('evaluation_history'):
            issues.append('Paid corrector evaluation history unavailable')
        points.extend(local)
        accounting.append(dict(task_id=row['task_id'],case_id=row['case_id'],arm=row['arm'],
            selection='ALL six objects in preselected stratum, every fixed arm, no winner filter',
            stratum=stratum,physical_points=len(local),initial_point_has_cost_coordinates=False,
            cpu_origin_offset_seconds=None,cpu_origin_offset_status='NOT_LOGGED; no fabricated alignment',
            sources=[str(history_path),str(corrector_path)],issues=issues))
    return points,accounting,case_ids


def make_plots(rows,points,case_ids,figure_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    colors=dict(zip(ARMS,('#4477AA','#EE6677','#228833','#AA3377')))
    markers=dict(zip(STRATA,('o','s','^')))
    outputs=[]
    for batch in ('primary','fft'):
        fig,ax=plt.subplots(figsize=(8.5,5.2))
        group=[r for r in rows if r['batch']==batch]
        plotted=0
        for r in group:
            if not numeric(r.get('cpu_seconds')) or not numeric(r.get('material_nrmse')):continue
            ax.scatter(r['cpu_seconds'],r['material_nrmse'],s=45,alpha=.8,
                       marker=markers[r['stratum']],color=colors[r['arm']],edgecolors='white',linewidths=.45)
            plotted+=1
        ax.set_xlabel('Actual task receipt CPU (s)');ax.set_ylabel('Final material NRMSE')
        ax.set_title(('Primary: independent material objects' if batch=='primary' else
                      'FFT engineering duplicate observations')+f' ({plotted}/{len(group)} finite task points)')
        arm_handles=[Line2D([],[],marker='o',linestyle='',color=colors[a],label=a)
                     for a in ARMS if any(r['arm']==a for r in group)]
        stratum_handles=[Line2D([],[],marker=markers[s],linestyle='',color='0.35',label=s)
                         for s in STRATA if any(r['stratum']==s for r in group)]
        first=ax.legend(handles=arm_handles,loc='upper right',fontsize=8);ax.add_artist(first)
        ax.legend(handles=stratum_handles,loc='upper left',fontsize=7)
        ax.grid(alpha=.2);fig.tight_layout()
        path=figure_dir/f'{batch}_material_nrmse_vs_actual_cpu.png'
        fig.savefig(path,dpi=160);plt.close(fig);outputs.append(str(path))
    labels={r['case_id']:r['label'] for r in rows if r['batch']=='primary'}
    for xfield,xlabel,suffix in (
        ('cumulative_LU','Cumulative actual physical LU','LU'),
        ('cumulative_observed_cpu_seconds','Recorded cumulative CPU (s; see origin note)','CPU')):
        fig,axes=plt.subplots(2,3,figsize=(14,8),squeeze=False)
        for ident,ax in zip(case_ids,axes.ravel()):
            for arm in ARMS:
                p=[p for p in points if p['case_id']==ident and p['arm']==arm and
                   numeric(p.get(xfield)) and numeric(p.get('best_original_Phi'))]
                if not p:continue
                ax.step([r[xfield] for r in p],[r['best_original_Phi'] for r in p],where='post',
                        lw=1.35,color=colors[arm],label=arm,alpha=.9)
                stage=[r for r in p if r['phase']=='auxiliary_stage']
                ax.scatter([r[xfield] for r in stage],[r['best_original_Phi'] for r in stage],
                           color=colors[arm],s=15,zorder=3)
            ax.set_title(labels[ident]+'\n'+ident,fontsize=9)
            ax.set_yscale('log');ax.set_ylabel('Observed incumbent original Phi');ax.set_xlabel(xlabel,fontsize=8)
            ax.grid(alpha=.2)
        handles=[Line2D([],[],color=colors[a],label=a) for a in ARMS]
        fig.legend(handles=handles,loc='upper center',ncol=4,fontsize=9)
        fig.suptitle('Every object in the fixed six-object stratum; paid physical samples only',y=.965,fontsize=11)
        fig.text(.5,.012,'Auxiliary endpoints are sparse; connecting steps do not evaluate intervening physical Phi.'+
                 (' Stage CPU starts at auxiliary entry; corrector CPU starts at observer creation; offset not logged.' if suffix=='CPU' else ''),
                 ha='center',fontsize=7)
        fig.tight_layout(rect=(0,.04,1,.92))
        path=figure_dir/f'fixed_six_original_Phi_vs_{suffix}.png'
        fig.savefig(path,dpi=160);plt.close(fig);outputs.append(str(path))
    return outputs


def fmt(value):
    if value is None:return '—'
    if isinstance(value,bool):return 'PASS' if value else 'FAIL'
    if numeric(value):return f'{value:.5g}'
    return str(value).replace('|','\\|').replace('\n',' ')


def markdown(rows,batches,pairs,figures,output_dir,attempts):
    lines=['# R6 offline data and cost tables','',
           'This artifact reports arithmetic, measured costs and accounting gaps. Scientific interpretation is reserved for the parent review.','',
           '## Object and task denominators','',
           '| Batch | Unique observed objects | Additional independent objects | Tasks | Workers | Batch elapsed s | Sum task wall s | Sum task CPU s | Physical LU | Complex RHS |',
           '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for b in batches:
        lines.append('| '+' | '.join(fmt(v) for v in (b['batch'],b['unique_material_objects'],b['additional_independent_objects'],
            b['total_tasks'],b['workers'],b['batch_elapsed_wall_seconds'],b['task_wall']['known_sum'],
            b['task_cpu']['known_sum'],b['physical_LU']['known_sum'],b['physical_actual_complex_RHS']['known_sum']))+' |')
    lines+=['','## Charged attempts including preserved failures','',
            '| Batch | Final tasks | Extra attempts | All charged attempts | Charged CPU s | Charged task wall s | Charged LU | Charged RHS |',
            '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for b in batches:
        lines.append('| '+' | '.join(fmt(v) for v in (b['batch'],b['total_tasks'],b['additional_charged_attempts'],
            b['all_charged_attempts_including_final'],b['all_charged_task_cpu']['known_sum'],
            b['all_charged_task_wall']['known_sum'],b['all_charged_physical_LU']['known_sum'],
            b['all_charged_actual_complex_RHS']['known_sum']))+' |')
    if attempts:
        lines+=['','| Additional attempt | Status | Receipt CPU s | Receipt wall s | Outer process wall s (overlaps) | LU | RHS | Evidence |',
                '| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |']
        for r in attempts:
            e=r.get('physical_zero_evidence') or {}
            lines.append('| '+' | '.join(fmt(v) for v in (r['task_id'],r['status'],r['cpu_seconds'],r['wall_seconds'],
                e.get('outer_queue_process_wall_seconds'),r['LU'],r['actual_complex_RHS'],e.get('status'))) + ' |')
    lines+=['','## Pairwise A material arithmetic','',
            'Positive NRMSE delta means baseline NRMSE minus A NRMSE is positive. Missing pairs remain missing.','',
            '| Batch | Stratum | Baseline | Paired objects | Known NRMSE pairs | A lower | A higher | Equal | Median NRMSE delta |',
            '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for p in pairs:
        lines.append('| '+' | '.join(fmt(v) for v in (p['batch'],p['stratum'],p['baseline'],p['paired_objects'],
            p['known_material_pairs'],p['A_lower_nrmse'],p['A_higher_nrmse'],p['equal_nrmse'],
            p['statistics']['material_nrmse_baseline_minus_A']['median']))+' |')
    for batch in ('primary','fft'):
        lines+=['',f'## All {batch} case tasks','',
                '| Stratum | Object label | Case | Arm | Status | NRMSE | IoU | Original Phi | CPU s | Wall s | LU | RHS | Joint |',
                '| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |']
        for r in sorted((r for r in rows if r['batch']==batch),key=lambda r:(STRATA.index(r['stratum']),r['label'],r['case_id'],ARMS.index(r['arm']))):
            gate=r.get('confirmation_joint') if r['stratum']=='NEW6_CONFIRMATION' else r.get('legacy_joint')
            lines.append('| '+' | '.join(fmt(v) for v in (r['stratum'],r['label'],r['case_id'],r['arm'],r['status'],
                r.get('material_nrmse'),r.get('otsu_iou'),r.get('original_Phi'),r.get('cpu_seconds'),r.get('wall_seconds'),
                r.get('LU'),r.get('actual_complex_RHS'),gate))+' |')
    lines+=['','## Accounting notes','']+['- '+n for n in NOTES]
    lines+=['','## Figures','']
    import os
    for figure in figures:lines.append(f'![{Path(figure).stem}]({os.path.relpath(figure,output_dir)})')
    return '\n'.join(lines)+'\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--primary-scores',type=Path,default=primary_inputs()['primary'])
    ap.add_argument('--fft-scores',type=Path,default=primary_inputs()['fft'])
    ap.add_argument('--primary-manifest',type=Path,default=ROOT/'configs/primary_final_manifest.json')
    ap.add_argument('--fft-manifest',type=Path,default=ROOT/'configs/FFT_manifest.json')
    ap.add_argument('--primary-queue',type=Path,action='append',
                    help='Repeat for sequential queue segments; default original primary + completion queue')
    ap.add_argument('--fft-queue',type=Path,action='append')
    ap.add_argument('--attempts-manifest',type=Path,default=ROOT/'configs/extra_attempts.json',
                    help='Optional preserved additional charged attempts: [{batch,task_id,receipt_path,note}]')
    ap.add_argument('--output',type=Path,default=ROOT/'RESULTS/analysis')
    ap.add_argument('--figures',type=Path,default=ROOT/'FIGURES/analysis')
    ap.add_argument('--curve-stratum',choices=('NEW6_CONFIRMATION','DEV6_DIFFICULT'),default='NEW6_CONFIRMATION',
                    help='Plot every object in this fixed six-object stratum; never select wins')
    args=ap.parse_args()
    args.primary_queue=args.primary_queue or [ROOT/'audit/queue_primary/RECEIPT.json',ROOT/'audit/queue_primary_completion/RECEIPT.json']
    args.fft_queue=args.fft_queue or [ROOT/'audit/queue_fft/RECEIPT.json']
    # Hard precondition before opening either score file, receipts or histories.
    missing=[str(p) for p in (args.primary_scores,args.fft_scores) if not p.is_file()]
    if missing:ap.error('Both completed score files are required; missing: '+', '.join(missing))
    if args.output.exists() or args.figures.exists():
        ap.error('Write-once output/figure directories must not already exist; use a new explicit destination')
    wall,cpu=time.perf_counter(),time.process_time();seen=set();rows=[];receipts={};batches=[]
    for batch in ('primary','fft'):
        data,recs,summary=collect_batch(batch,getattr(args,batch+'_scores'),getattr(args,batch+'_manifest'),
                                      getattr(args,batch+'_queue'),seen)
        rows.extend(data);receipts.update(recs);batches.append(summary)
    final_receipt_count=len(seen)
    attempts,attempt_issues=collect_extra_attempts(args.attempts_manifest,seen)
    for summary in batches:
        extra=[r for r in attempts if r['batch']==summary['batch']]
        final=[r for r in rows if r['batch']==summary['batch']]
        summary['additional_charged_attempts']=len(extra)
        summary['all_charged_attempts_including_final']=len(extra)+len(final)
        summary['all_charged_task_cpu']=stats([r.get('cpu_seconds') for r in final+extra])
        summary['all_charged_task_wall']=stats([r.get('wall_seconds') for r in final+extra])
        summary['all_charged_physical_LU']=stats([r.get('LU') for r in final+extra])
        summary['all_charged_actual_complex_RHS']=stats([r.get('actual_complex_RHS') for r in final+extra])
        if (summary['queue_failed_attempt_entries'] or 0)>len(extra):
            summary['accounting_issues'].append(dict(issue='Queue failure entries exceed explicitly preserved extra attempts; inspect missing historical costs'))
    primary_ids={r['case_id'] for r in rows if r['batch']=='primary'}
    fft_ids={r['case_id'] for r in rows if r['batch']=='fft'}
    if not fft_ids<=primary_ids:raise ValueError('FFT duplicate set contains observations outside primary benchmark')
    summaries=[]
    for batch in ('primary','fft'):
        arms=[a for a in ARMS if any(r['batch']==batch and r['arm']==a for r in rows)]
        for stratum in ('ALL',)+STRATA:
            for arm in ('ALL',)+tuple(arms):
                group=[r for r in rows if r['batch']==batch and (stratum=='ALL' or r['stratum']==stratum)
                       and (arm=='ALL' or r['arm']==arm)]
                if group:summaries.append(summarize(group,batch,stratum,arm))
    pairs,pair_summaries=pairwise(rows);engineering=backend_pairs(rows)
    charged_groups=charged_cost_groups(rows,attempts)
    points,curve_accounting,case_ids=collect_curves(rows,receipts,args.curve_stratum)
    args.output.mkdir(parents=True,exist_ok=False);args.figures.mkdir(parents=True,exist_ok=False)
    write_json(args.output/'AGGREGATE_COSTS.json',dict(batch_summaries=batches,arm_stratum_summaries=summaries,
        primary_independent_objects=len(primary_ids),new_confirmation_independent_objects=len({r['case_id'] for r in rows if r['batch']=='primary' and r['stratum']=='NEW6_CONFIRMATION'}),
        fft_observation_duplicates=len(fft_ids),fft_additional_independent_objects=0,
        additional_attempt_accounting_issues=attempt_issues,notes=NOTES))
    write_json(args.output/'TASK_LEDGER.json',rows);write_csv(args.output/'TASK_LEDGER.csv',rows)
    write_json(args.output/'ADDITIONAL_ATTEMPTS.json',attempts);write_csv(args.output/'ADDITIONAL_ATTEMPTS.csv',attempts)
    write_json(args.output/'ALL_CHARGED_ARM_STRATUM_COSTS.json',charged_groups)
    write_json(args.output/'PAIRWISE_A_VS_BASELINES.json',dict(rows=pairs,summaries=pair_summaries))
    write_csv(args.output/'PAIRWISE_A_VS_BASELINES.csv',pairs)
    write_json(args.output/'FFT_ENGINEERING_MATCHED.json',engineering)
    write_csv(args.output/'FFT_ENGINEERING_MATCHED.csv',engineering)
    write_csv(args.output/'ORIGINAL_PHI_CURVES.csv',points)
    write_json(args.output/'CURVE_ACCOUNTING.json',dict(stratum=args.curve_stratum,case_ids=case_ids,tasks=curve_accounting,notes=NOTES[-4:-1]))
    figures=make_plots(rows,points,case_ids,args.figures)
    (args.output/'ALL_CASES.md').write_text(markdown(rows,batches,pair_summaries,figures,args.output,attempts))
    write_json(args.output/'ANALYSIS_RECEIPT.json',dict(status='COMPLETE_OFFLINE_ARITHMETIC_ONLY',
        created_utc=datetime.now(timezone.utc).isoformat(),primary_scores=str(args.primary_scores),fft_scores=str(args.fft_scores),
        unique_final_task_receipts_counted=final_receipt_count,unique_all_attempt_receipts_counted=len(seen),
        additional_attempts=len(attempts),attempts_manifest=str(args.attempts_manifest) if args.attempts_manifest else None,
        task_rows=len(rows),curve_stratum=args.curve_stratum,
        primary_curve_case_ids=case_ids,new_physical_calls=0,new_GT_files_opened=0,new_prediction_files_opened=0,
        figures=figures,output=str(args.output),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        scientific_review='NOT_PERFORMED; parent scientific interpretation required'))
    print(json.dumps(dict(status='COMPLETE_OFFLINE_ARITHMETIC_ONLY',tasks=len(rows),output=str(args.output),figures=str(args.figures))))


if __name__=='__main__':main()
