"""Metadata-only partial/final R14 cost and typed full-state catalog.

No numerical solver, optimizer, GT/material scoring, or NPZ numeric payload.
CPU and physical-counter inclusion are independent; nested/reused views stay
visible but cannot silently enter the additive totals twice.
"""
from __future__ import annotations
import argparse,ast,csv,json,math,re,struct,time,zipfile
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

ROWS=[];INPUTS=[];ISSUES=[];SEEN=set();QUEUE_CHILDREN={}


def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def read(path):
    path=Path(path).resolve()
    if any(p.lower() in ('sealed_truth','truth','data_private') for p in path.parts):raise ValueError('GT input forbidden')
    if not path.is_file():return None
    try:
        before=path.stat();data=json.loads(path.read_text(encoding='utf-8'));after=path.stat()
        INPUTS.append(dict(path=str(path),bytes=before.st_size,mtime_ns=before.st_mtime_ns,
            stable_during_read=(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)))
        return data
    except Exception as exc:ISSUES.append(dict(path=str(path),issue='UNREADABLE_JSON',error=repr(exc)));return None


def write(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def table(out,name,rows):
    write(out/(name+'.json'),rows);keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/(name+'.csv')).open('x',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})


def physical(cost):
    if not isinstance(cost,dict):return None,None
    lu=cost.get('physical_LU',cost.get('LU'))
    if 'actual_complex_RHS' in cost:rhs=cost['actual_complex_RHS']
    elif 'physical_complex_RHS' in cost:rhs=cost['physical_complex_RHS']
    elif 'source_RHS' in cost and 'adjoint_RHS' in cost:rhs=cost['source_RHS']+cost['adjoint_RHS']
    else:rhs=None
    return lu if finite(lu) else None,rhs if finite(rhs) else None


def row(ident,path,category,host,role,cpu=None,wall=None,counts=None,include_cpu=True,include_physics=True,
        parent=None,reused=False,status=None,scope=None,reason=None,**extra):
    if ident in SEEN:raise ValueError('Duplicate cost identity: '+ident)
    SEEN.add(ident);lu,rhs=physical(counts)
    item=dict(cost_id=ident,source_path=str(path),category=category,host=host,role=role,status=status,
        CPU_seconds=cpu if finite(cpu) else None,wall_seconds=wall if finite(wall) else None,
        CPU_null_reason=None if finite(cpu) else reason or 'No measured CPU field in source receipt;not zero',
        physical_LU=lu,physical_complex_RHS=rhs,physical_null_reason=None if lu is not None and rhs is not None else 'Missing original physical counter(s);operator columns/algebraic factors not substituted',
        include_CPU_in_additive_total=include_cpu,include_physics_in_additive_total=include_physics,
        parent_cost_id=parent,reused_evidence=reused,CPU_scope=scope,
        failed_or_incomplete=status in ('FAILED','LAUNCH_FAILED','COMPLETED_WITH_ERRORS','INCOMPLETE_SETUP') or extra.get('returncode') not in (None,0) or extra.get('analysis_failure') is True,
        raw_physical_counters=counts,**extra)
    ROWS.append(item);return item


def first_cpu(d):
    for key in ('OS_child_total_CPU_seconds','offline_whole_process_CPU_seconds','whole_process_CPU_seconds','whole_process_CPU_seconds_at_receipt','process_CPU_seconds_including_imports',
                'process_CPU_seconds','CPU_seconds','offline_CPU_seconds','cpu_seconds'):
        if finite(d.get(key)):return d[key],key
    return None,None


def queue_costs(root):
    for path in sorted((root/'queues').rglob('RECEIPT.json'),key=lambda p:(len(p.relative_to(root/'queues').parts),str(p))):
        d=read(path)
        if not isinstance(d,dict):continue
        queue_rel=path.parent.relative_to(root/'queues');name=d.get('campaign',queue_rel.parts[0]);jobs=d.get('jobs',[])
        external=[(j.get('external') or {}).get('child_process_CPU_seconds',j.get('child_process_CPU_seconds',j.get('child_CPU_seconds'))) for j in jobs]
        aggregate=d.get('solver_children_whole_OS_CPU_seconds',d.get('whole_children_OS_CPU_seconds',d.get('finished_child_CPU_seconds')))
        prefer_aggregate=finite(aggregate) and (not external or not all(finite(x) for x in external))
        host='Windows' if any('Windows' in str((j.get('external') or {}).get('CPU_scope')) for j in jobs) else 'Apple_local' if 'solver_children_whole_OS_CPU_seconds' in d or d.get('host')=='darwin' else 'UNKNOWN'
        host_note='Windows GetProcessTimes declared by jobs' if host=='Windows' else 'Local Apple controller receipt family;CPU accounting kept separate' if host=='Apple_local' else 'Queue receipt lacks explicit OS evidence'
        qid='queue:'+str(queue_rel)
        controller=d.get('controller_CPU_seconds');wall=d.get('queue_wall_seconds',d.get('controller_wall_seconds',d.get('controller_wall_including_source_freeze_seconds',d.get('wall_seconds'))))
        row(qid+':controller',path,'online',host,'controller',controller,wall,include_physics=False,
            status=d.get('status'),scope='Controller CPU only,disjoint from solver children;queue wall never added to CPU',host_qualification=host_note)
        row(qid+':children_aggregate',path,'online',host,'whole_child_CPU_aggregate',aggregate,wall,
            include_cpu=prefer_aggregate,include_physics=False,parent=None if prefer_aggregate else qid+':children',
            status='COMPLETED_WITH_ERRORS' if any(j.get('returncode') not in (None,0) for j in jobs) else d.get('status'),scope='Aggregate child OS CPU,not additional to per-job OS CPU',reused=True,queue_status=d.get('status'))
        all_subprocess=d.get('all_controller_subprocess_OS_CPU_seconds')
        if finite(all_subprocess) and finite(aggregate):
            row(qid+':controller_subprocess_support',path,'setup',host,'non_solver_controller_subprocess_remainder',all_subprocess-aggregate,
                include_physics=False,status=d.get('status'),scope='Same receipt allcontroller subprocessCPU minus solverchildrenCPU;not full aggregate twice')
        for j,value in zip(jobs,external):
            ident=j.get('task_id','UNKNOWN');cost_id=qid+':child:'+ident
            row(cost_id,path,'online',host,'solver_child_external',value,j.get('wall_seconds'),
                include_cpu=not prefer_aggregate,include_physics=False,parent=qid+':children_aggregate' if prefer_aggregate else None,
                status=j.get('status',j.get('child_status')),scope=(j.get('external') or {}).get('CPU_scope','External child field scope not explicitly recorded'),
                reason='Per-job external CPU absent;aggregate may exist but is not divided among tasks',
                returncode=j.get('returncode'),task_id=ident,campaign=name,remote_receipt_path=j.get('receipt'),log_path=j.get('log'),core=j.get('cpu'),
                child_internal_CPU_seconds=j.get('child_internal_CPU_seconds',j.get('child_reported_CPU_seconds')),child_receipt_status=j.get('child_status'),queue_attempt=str(queue_rel))
            QUEUE_CHILDREN.setdefault((name,ident),[]).append(cost_id)


def result_costs(root):
    vector_total_path=root/'analysis/vector4_v1/COST_TOTALS.json';vector=read(vector_total_path)
    vector_authority=isinstance(vector,dict)
    if vector_authority:
        children=vector.get('solver_children_whole_OS_CPU_seconds');all_subprocess=vector.get('all_controller_subprocess_OS_CPU_seconds')
        queue_vector=any(r['cost_id'].startswith('queue:vector4_v1') and r['role']=='whole_child_CPU_aggregate' and r['include_CPU_in_additive_total'] for r in ROWS)
        row('vector4:external_children',vector_total_path,'online','Apple_local','solver_children_external_aggregate',children,
            counts=dict(physical_LU=vector.get('physical_LU'),physical_complex_RHS=vector.get('physical_complex_RHS')),
            reused=True,include_cpu=not queue_vector,scope='Reused actual-cost summary;queue OS aggregate preferred when available;body timers are nested',
            original_controller_receipt_available=queue_vector)
        support=all_subprocess-children if finite(all_subprocess) and finite(children) else None
        row('vector4:controller_subprocess_support',vector_total_path,'setup','Apple_local','non_solver_controller_subprocess_remainder',support,
            include_cpu=not queue_vector,include_physics=False,reused=True,scope='all_controller_subprocess_OS_CPU minus solver_children_OS_CPU from same summary;not fullaggregate or queue remainder again')
        row('vector4:controller',vector_total_path,'online','Apple_local','controller',vector.get('controller_CPU_seconds'),vector.get('controller_wall_seconds'),
            include_cpu=not queue_vector,include_physics=False,reused=True,scope='Controller CPU disjoint from children;queue controller preferred when same source available')
    for path in sorted((root/'RESULTS').rglob('RECEIPT.json')):
        rel=path.relative_to(root/'RESULTS');campaign=rel.parts[0]
        d=read(path)
        if not isinstance(d,dict):continue
        task=d.get('task',{});ident=task.get('task_id',rel.parts[1] if len(rel.parts)>1 else str(rel))
        candidates=QUEUE_CHILDREN.get((campaign,ident),[]);vec=campaign=='vector4_v1' and vector_authority
        cpu,key=first_cpu(d)
        matches=[q for q in ROWS if q['cost_id'] in candidates and q.get('returncode')==0 and q.get('child_receipt_status')==d.get('status') and (not finite(q.get('child_internal_CPU_seconds')) or cpu is not None and abs(q['child_internal_CPU_seconds']-cpu)<1e-8)]
        parent=matches[-1]['cost_id'] if matches else None
        if candidates and parent is None:ISSUES.append(dict(path=str(path),issue='NO_MATCHING_CONTROLLER_ATTEMPT_FOR_CURRENT_RESULT',candidate_attempt_cost_ids=candidates,action='Current bodyCPU counted separately as fallback;prior failed/other attempts retained'))
        counts=d.get('physical_cost',d.get('operation_counts_all',d.get('physical_model_cost')))
        if counts is None and any(k in d for k in ('physical_LU','actual_complex_RHS')):counts=d
        host='Apple_local' if campaign=='vector4_v1' or campaign.startswith('shared_readout') else 'Windows' if parent and any(r['cost_id']==parent and r['host']=='Windows' for r in ROWS) else 'UNKNOWN'
        row('result:'+str(rel),path,'online',host,'solver_child_internal',cpu,d.get('wall_seconds'),counts,
            include_cpu=parent is None and not vec,include_physics=not vec,parent=parent or ('vector4:external_children' if vec else None),
            reused=False,status=d.get('status'),scope='Internal '+str(key)+';nested under external child when present;finalserialization/exit may be excluded',
            campaign=campaign,task_id=ident,auxiliary_cost=d.get('auxiliary_cost',d.get('auxiliary_receipt')),
            source_declaration=d.get('git_source'),CPU_authority='Fallback internal timing only;full lifetime unknown' if parent is None and not vec else 'Nested detail,not added')


def intervention_costs(root):
    directory=root/'audit/interventions';batch=read(directory/'BATCH_RECEIPT.json')
    jobs=(batch or {}).get('jobs',[]);lookup={}
    for j in jobs:
        label=f"{j.get('case_id')}_k{int(j.get('checkpoint',0)):05d}";ident='intervention:'+label+':external'
        lookup[label]=ident
        row(ident,directory/'BATCH_RECEIPT.json','validation','Apple_local','diagnostic_child_wrapper',j.get('child_CPU_seconds'),j.get('wall_seconds'),
            include_physics=False,status='COMPLETE' if j.get('returncode')==0 else 'FAILED',returncode=j.get('returncode'),
            scope='Batch-reported child_CPU_seconds;exact OS-collection code not recorded in this receipt;do not add body CPU')
    if batch is not None:row('intervention:controller',directory/'BATCH_RECEIPT.json','validation','Apple_local','batch_controller',
        wall=batch.get('wall_seconds'),include_physics=False,status='COMPLETE',reason='Only batch wall recorded;controller CPU not measured')
    for path in sorted(directory.glob('*/RECEIPT.json')):
        d=read(path)
        if not isinstance(d,dict):continue
        name=path.parent.name;cpu,key=first_cpu(d)
        row('intervention:'+name+':body',path,'validation','Apple_local','diagnostic_child_internal',cpu,d.get('diagnostic_wall_seconds'),
            d.get('physical_operation_counts'),include_cpu=name not in lookup,parent=lookup.get(name),status=d.get('status'),
            scope='Internal '+str(key)+';actual new intervention physics,not historical source-checkpoint counts',
            auxiliary_cost=d.get('diagnostic_auxiliary_operation_counts'),setup_counts=d.get('setup_operation_counts'),
            historical_checkpoint_counts_recharged=False)
        ledger=read(path.parent/'COST_LEDGER.json')
        if isinstance(ledger,list):
            for event in ledger:
                row('intervention:'+name+':event:'+str(event.get('event')),path.parent/'COST_LEDGER.json','validation','Apple_local','nested_event',
                    event.get('CPU_seconds'),event.get('wall_seconds'),event.get('physical_operation_delta'),False,False,
                    parent='intervention:'+name+':body',status=event.get('status'),scope='Nested paid event breakdown only',
                    event_label=event.get('label'),auxiliary_delta=event.get('auxiliary_operation_delta'))


def metadata_costs(root,out):
    handled=set()
    generation=root/'audit/GENERATION24_RECEIPT.json';g=read(generation)
    if isinstance(g,dict):
        row('generation24',generation,'setup','Apple_local','data_generation',g.get('process_CPU_seconds'),g.get('process_wall_seconds'),
            dict(LU=g.get('LU'),actual_complex_RHS=g.get('complex_RHS')),status=g.get('status'),
            scope='Entry through final generator receipt sample;final receiptwrite/exit excluded;case/setup timers nested',
            setup_details=g.get('setup'),failed_cases=g.get('failed_cases'),prior_probe=g.get('earlier_interpreter_probe'))
        for i,case in enumerate(g.get('case_costs',[])):
            row('generation24:case:'+str(i),generation,'setup','Apple_local','nested_generation_case',
                case.get('CPU_seconds_including_serialization',case.get('CPU_seconds')),
                case.get('wall_seconds_including_serialization',case.get('wall_seconds')),case.get('counts'),False,False,
                parent='generation24',status=case.get('status'),scope='Nested in full generator cost',case_id=case.get('case_id'))
        probe=g.get('earlier_interpreter_probe')
        if isinstance(probe,dict):row('generation24:prior_interpreter_probe',generation,'setup','Apple_local','failed_prior_import_probe',
            probe.get('CPU_seconds'),counts=dict(LU=probe.get('physical_LU'),actual_complex_RHS=probe.get('physical_RHS')),
            status='FAILED',scope=probe.get('scope'),reason='Earlier import probe CPU unmeasured;not included in generator process')
        handled.add(generation.resolve())
    candidates=[]
    for base_dir in [root/'audit',root/'analysis',root/'archives']:
        for p in base_dir.rglob('*.json'):
            if p.resolve().is_relative_to(out) or any(part.startswith('cost_ledger') for part in p.relative_to(root).parts):continue
            rel=p.relative_to(root);parts=rel.parts
            if 'interventions' in parts or 'trajectories' in parts:continue
            if 'RECEIPT' in p.name or 'DOWNLOAD_EXTRACT' in p.name or p.name.startswith('PLAN_FAILURE_') or p.name in ('SUMMARY.json','METRIC_QUALIFICATION.json','CAMPAIGN_PREPARATION.json','INDEPENDENT_PREPARATION.json','SELF_ONLY_PREPARATION.json'):
                # SUMMARY is included only for actual standalone local-closed-loop attempts.
                if p.name=='SUMMARY.json' and not p.parent.name.startswith('local_closed_loop'):continue
                candidates.append(p)
    for p in sorted(candidates):
        if p.resolve() in handled:continue
        d=read(p)
        if not isinstance(d,dict):continue
        rel=str(p.relative_to(root));cpu,key=first_cpu(d)
        category='scoring' if p.parent.is_relative_to(root/'analysis') or 'threshold_precision' in rel or 'TICK_FIX' in p.name else 'transport' if p.name.endswith('_ARCHIVE_RECEIPT.json') or 'DOWNLOAD_EXTRACT' in p.name else 'setup' if any(x in rel for x in ['PREPARATION','PSUTIL_RUNTIME']) else 'validation'
        host='Windows' if p.name.endswith('_ARCHIVE_RECEIPT.json') and 'dualmetric' in p.name else 'Apple_local'
        counts=d.get('physical_cost',d.get('physical_counts'))
        if p.parent.is_relative_to(root/'analysis'):
            counts=dict(LU=d.get('new_physical_LU',d.get('forward_LU')),actual_complex_RHS=d.get('new_RHS',d.get('complex_solver_RHS')))
        elif counts is None and any(k in d for k in ('forward_LU','new_forward_LU','new_physical_LU','LU')):
            counts=dict(LU=d.get('forward_LU',d.get('new_forward_LU',d.get('new_physical_LU',d.get('LU')))),actual_complex_RHS=d.get('complex_solver_RHS',d.get('solver_RHS',d.get('new_RHS',d.get('new_complex_RHS',d.get('complex_RHS'))))))
        elif d.get('new_forward')==0 and d.get('new_adjoint')==0 or d.get('new_physics')==0:counts=dict(LU=0,actual_complex_RHS=0)
        force_dir=root/'audit/material_force_decomposition_v1'
        force_parent=None
        if p.parent==force_dir and p.name=='RECEIPT.json' and (force_dir/'LAUNCH_RECEIPT.json').is_file():force_parent='metadata:audit/material_force_decomposition_v1/LAUNCH_RECEIPT.json'
        if p.parent==force_dir and p.name=='PLOT_REFRESH_RECEIPT.json' and (force_dir/'PLOT_LAUNCH_RECEIPT.json').is_file():force_parent='metadata:audit/material_force_decomposition_v1/PLOT_LAUNCH_RECEIPT.json'
        local_archive_root=root/'archives/local_fullstates_v1'
        archive_inner=(p.name=='PACKAGE_RECEIPT.json' and p.is_relative_to(local_archive_root) and (local_archive_root/'LOCAL_ARCHIVE_RECEIPT.json').is_file())
        if p.is_relative_to(local_archive_root):category='transport'
        row('metadata:'+rel,p,category,host,'nested_archive_phase' if archive_inner else 'standalone_process',cpu,d.get('wall_seconds',d.get('offline_wall_seconds')),counts,
            include_cpu=not archive_inner and force_parent is None,include_physics=not archive_inner,
            parent='metadata:archives/local_fullstates_v1/LOCAL_ARCHIVE_RECEIPT.json' if archive_inner else force_parent,
            status=d.get('status'),scope='Preferred most-complete recorded timer '+str(key)+';other body/all-attempt/child subtimers are nested and not added',
            raw_timer_fields={k:v for k,v in d.items() if ('CPU' in k or 'wall' in k) and isinstance(v,(int,float,str))},
            setup_counts=d.get('setup_ideal_counts',d.get('mechanism_setup_and_ideal_counts',d.get('setup_operation_counts'))),
            auxiliary_validation_counts=d.get('new_validation_operation_counts'),
            no_absolute_crosshost_speed_comparison=True,
            analysis_failure=p.parent.name=='exposed48_v1',
            analysis_failure_note='Parent-identified NATIVE export duplicatepair bug;preserve analysisCPU,not a physical prediction failure;v1-v2 numerical comparison pending' if p.parent.name=='exposed48_v1' else None)
        # Child trajectories execute inside these single-process diagnostic drivers.
        if p.name=='RECEIPT.json' and p.parent.name in ('finite_mode_trajectory_v1','coupled_mode_trajectory_v1'):
            for child in sorted((p.parent/'trajectories').glob('*/RECEIPT.json')):
                c=read(child)
                if not isinstance(c,dict):continue
                row('nested:'+str(child.relative_to(root)),child,'validation','Apple_local','nested_same_process_trajectory',c.get('CPU_seconds'),c.get('wall_seconds'),
                    c.get('physical_delta'),False,False,parent='metadata:'+rel,status=c.get('status'),
                    scope='Driver process already includes this trajectory;not an independent child process',auxiliary_cost=c.get('native_step_operation_counts'))
    for old in sorted((root/'audit').glob('cost_ledger*/RECEIPT.json')):
        if old.resolve().is_relative_to(out):continue
        d=read(old)
        if isinstance(d,dict):row('prior_ledger_build:'+str(old.relative_to(root)),old,'validation','Apple_local','prior_metadata_ledger_process',
            d.get('new_CPU_seconds'),d.get('new_wall_seconds'),dict(LU=d.get('new_physics'),actual_complex_RHS=d.get('new_physics')),
            status=d.get('status'),scope='Only this previous ledger build own measuredCPU;its copied experiment cost totals are not added again')
    qualification=root/'FAILURES.md'
    if qualification.is_file():
        text=qualification.read_text()
        if '0.409276' in text:
            row('documented_failure:metric_qualification_output01',qualification,'validation','Apple_local','documented_failed_metadata_formatting',
                wall=.409276,counts=dict(LU=0,actual_complex_RHS=0),status='FAILED',
                reason='FAILURES.md#5 explicitly records toolwall0.409276 anduncapturedprocessCPU',scope='Failed output-dict formatting after calculations;documentedzero forward/LU/RHS;no guessedCPU')
    reused=root/'audit/CAMPAIGN_PREPARATION.json'
    if reused.is_file():row('reused:historical48_public_inputs',reused,'reused','HISTORICAL_HOST_UNKNOWN','old_public_inputs',
        include_cpu=False,include_physics=False,reused=True,scope='R13public48/operator history reused;old generation/training cost unknown here andnot charged as newR14physics',
        reason='Historical preparation cost not reconstructed;new filesystem copy overhead separately represented')
    # Transfer files prove bytes exist,not their upload/download CPU or network bill.
    packages=[p for folder in (root/'private_operational',root/'archives') for p in folder.glob('*.zip')]
    for p in packages:
        row('unmeasured_transport:'+str(p.relative_to(root)),p,'transport','UNKNOWN','package_transfer_or_packaging_unmeasured',
            include_physics=False,status='FILE_PRESENT',reason='Artifact presence/bytes provide no actual upload/download CPU;pack CPU may overlap a separate measured archive receipt',
            bytes=p.stat().st_size,possible_overlap='Do not add a guessed package/transfer time to measured archive process')


def npy_header(stream):
    magic=stream.read(8)
    if not magic.startswith(b'\x93NUMPY'):raise ValueError('Not NPY')
    width=2 if magic[6]==1 else 4;length=struct.unpack('<H' if width==2 else '<I',stream.read(width))[0]
    return ast.literal_eval(stream.read(length).decode('latin1'))


def metadata_tree(archive):
    """Decode string metadata only;truncate huge history after preceding loop state."""
    if '__state_tree_json__.npy' not in archive.namelist():return {},'NO_TYPED_STATE_TREE'
    with archive.open('__state_tree_json__.npy') as f:
        head=npy_header(f);descr=head['descr']
        if 'U' not in descr:raise ValueError('Expected Unicode JSON metadata')
        chars=int(descr.split('U')[-1]);limit=min(chars*4,1048576)
        raw=f.read(limit);text=raw.decode('utf-32-be' if descr.startswith('>') else 'utf-32-le')
        marker=text.find(',["history",')
        if chars*4<=limit:tree=json.loads(text);scope='FULL_METADATA_JSON_ONLY'
        elif marker>=0:tree=json.loads(text[:marker]+']}');scope='PREFIX_THROUGH_STAGE_BEFORE_HISTORY;later numeric payloads/headernames only'
        else:return {},'METADATA_PREFIX_NOT_PARSEABLE_WITHOUT_LARGE_HISTORY;unknown scalars retained'
    def decode(v):
        if not isinstance(v,dict):return v
        tag=v.get('type')
        if tag=='dict':return {decode(k):decode(x) for k,x in v['items']}
        if tag in ('tuple','list'):return tuple(decode(x) for x in v['items'])
        if tag=='ndarray':return {'typed_array_key':v['key']}
        return v
    return decode(tree),scope


def checkpoint_catalog(root):
    rows=[];fields=[]
    targets=[]
    for folder in (root/'RESULTS',root/'audit'):
        for p in folder.rglob('*.npz'):
            if any(part.startswith('cost_ledger') for part in p.relative_to(root).parts):continue
            if any(s in p.name.upper() for s in ('CHECKPOINT','STATE','PROPOSAL_IN_FLIGHT')):targets.append(p)
    wanted=('m','A','J','Q','Vd','JD','GJD','SJD','Od','Os','Ls','old_grad','old_direction',
        'prstate__gam','prstate__w','prstate__direction','basis_mf1','basis_mf2','observations',
        'incident_field','receiver_operator','green_cache__spectrum','onehop_cache__indices','onehop_cache__goff','onehop_cache__diag_g')
    for path in sorted(set(targets)):
        try:
            with zipfile.ZipFile(path) as archive:
                names=set(archive.namelist());tree,scope=metadata_tree(archive);headers={}
                for name in wanted:
                    key=name+'.npy'
                    if key in names:
                        with archive.open(key) as stream:headers[name]=npy_header(stream)
                if not all(k in headers for k in ('m','A','J','Q')):
                    status='PARTIAL_OR_NON_PR_STATE';kind='UNKNOWN_OR_PARTIAL'
                else:
                    n=headers['m']['shape'][0];nr=headers['J']['shape'][0]
                    kind='VECTOR_CELL_MAJOR_XYZ' if nr==3*n else 'SCALAR_TM_NATIVE_FOURIER_APCS' if nr==n else 'UNKNOWN'
                    status='TYPED_AUXILIARY_STATE_CONTEXT_PRESENT' if '__state_tree_json__.npy' in names else 'ARRAYS_WITHOUT_LOOP_METADATA'
                pair=tree.get('pair')
                source='typed JSON' if pair is not None else None
                if pair is None:
                    hit=re.search(r'pair([01]{2})',str(path))
                    if hit:pair=hit.group(1);source='path label,not metadata-certified'
                n=headers.get('m',{}).get('shape',(None,))[0]
                jshape=headers.get('J',{}).get('shape',())
                ns=jshape[1] if len(jshape)>1 else None
                memory=(('old_grad' in headers and 'old_direction' in headers) or
                    isinstance(tree.get('prstate'),dict))
                absent_memory=(tree.get('old_grad','ABSENT') is None and tree.get('old_direction','ABSENT') is None)
                rowinfo=dict(path=str(path),bytes=path.stat().st_size,kind=kind,status=status,pair=pair,pair_source=source,
                    material_real_dimension=n,current_complex_shape=list(jshape),coefficient_complex_shape=list(headers.get('A',{}).get('shape',())),
                    sources=ns,stage=tree.get('stage'),stage_iteration=tree.get('iteration'),global_iteration=tree.get('iterations',tree.get('global_iteration')),
                    stage_pending=tree.get('stage_pending'),done=tree.get('done'),PR_memory_arrays_present=memory,
                    PR_memory_explicit_none=absent_memory,metadata_read_scope=scope,
                    cached_Green_present='green_cache__spectrum.npy' in names,actual4cell_block_cache_present='onehop_cache__goff.npy' in names,
                    source_metric_context_present=('Os' in headers and 'Od' in headers) or 'state_weights.npy' in names,
                    numeric_payloads_read=False,GT_access=False,
                    representation='Native potentially redundant Fourier coefficients with fixed receiver DPCS;ambient A count is not independent spatial current rank' if kind=='SCALAR_TM_NATIVE_FOURIER_APCS' else 'Cell-major xyz,shared real material;vector range basis andprstate require original compatible acquisition/model' if kind=='VECTOR_CELL_MAJOR_XYZ' else 'Insufficient typed context;do not call restart complete',
                    restore_scope='Auxiliary numerical state only;re-supply exact publicY/operators,regularizer/runtime callbacks. CPUclock/selection/controller/physicalLU cache not certified restored',
                    source_representation_fields=sorted(k[:-4] for k in names if any(t in k for t in ('basis_','prior__','E2','weights','frequencies','receiver_'))))
                rows.append(rowinfo)
                for name,h in headers.items():fields.append(dict(checkpoint_path=str(path),field=name,dtype=h['descr'],shape=list(h['shape']),fortran_order=h['fortran_order'],payload_read=False))
        except Exception as exc:
            ISSUES.append(dict(path=str(path),issue='CHECKPOINT_METADATA_ERROR',error=repr(exc)))
            rows.append(dict(path=str(path),status='UNREADABLE_METADATA',error=repr(exc),numeric_payloads_read=False))
    return rows,fields


def coverage(root):
    output=[]
    for path in sorted((root/'configs').glob('*.json')):
        cfg=read(path)
        if not isinstance(cfg,dict) or not isinstance(cfg.get('tasks'),list):continue
        name=cfg.get('campaign',path.stem);tasks=cfg['tasks'];found=0;failed=0;checkpoint_tasks=0;checkpoint_files=0
        for task in tasks:
            folder=root/'RESULTS'/name/task['task_id']
            states=list(folder.rglob('CHECKPOINT*.npz')) if folder.is_dir() else []
            checkpoint_tasks+=bool(states);checkpoint_files+=len(states)
            files=[folder/'RECEIPT.json']+list(folder.glob('attempt_*/RECEIPT.json'))
            exists=next((p for p in files if p.is_file()),None)
            if exists:
                found+=1;d=read(exists) or {};failed+=d.get('status')=='FAILED'
        q=read(root/'queues'/name/'RECEIPT.json')
        output.append(dict(campaign=name,manifest=str(path),expected_tasks=len(tasks),local_child_receipts=found,
            failed_receipts=failed,local_tasks_with_checkpoint=checkpoint_tasks,local_checkpoint_files=checkpoint_files,
            state_scope='Local availability only;absence may reflect archive scope or an algorithm without fullPR checkpoints,not a failed run',missing_local_child_receipts=len(tasks)-found,queue_status=(q or {}).get('status'),
            execution_role='MANIFEST_ONLY_SUPERSET_NOT_ADDITIONAL_RUN' if name=='dualmetric48_v1' else 'SEPARATE_CAMPAIGN',
            scope='Local evidence availability,not inference that missing jobs never ran;main48manifest overlaps pilot+remainder and isnot an additional executed cohort'))
    return output


def write_solver_root_report(root,out,verified,campaigns,dimensions,dr):
    def fmt(v):return f"{v:.6f}" if finite(v) else '未记录'
    def dr_sum(key):
        vals=[x.get(key) for x in dr]
        return sum(vals) if all(finite(v) for v in vals) else None
    with (root/'COST_LEDGER.md').open('w',encoding='utf-8') as f:
        f.write('# 本轮逆解阶段成本账本\n\n')
        f.write(f"已按真实终态回执核对 **576个成功solver child**，实际OS CPU合计 **{verified['successful_children_actual_OS_CPU_seconds']:.6f}秒**。其中Windows为{verified['Windows_success_CPU_seconds']:.6f}秒，Apple为{verified['Apple_success_CPU_seconds']:.6f}秒。另有shared旧失败attempt的16个child，实际OS CPU为{verified['failed_extra_actual_OS_CPU_seconds']:.6f}秒；成功与失败solver合计{verified['all_success_and_failed_solver_OS_CPU_seconds']:.6f}秒。该计数表示进程/费用闭合，不表示收敛、恢复或科学验收。\n\n")
        f.write('|Campaign|设备|成功child|实际child OS CPU（秒）|该队列wall（秒）|物理LU|物理复RHS|\n|---|---|---:|---:|---:|---:|---:|\n')
        for c in campaigns:f.write(f"|{c['campaign']}|{c['host']}|{c['successful_children']}|{fmt(c['actual_children_OS_CPU_seconds'])}|{fmt(c['queue_wall_seconds'])}|{c['physical_LU']}|{c['physical_complex_RHS']}|\n")
        f.write('\n每个wall只描述对应队列，**不相加为整体并行wall**，也不作跨设备绝对速度结论。Windows五队列已逐child核对GetProcessTimes与队列合计。Apple shared/vector只有批次OS总量；按method的OS CPU未独立采集，因此在细表保留null，不按body时间比例分摊。各method的实测body计时另列，且不与外层OS CPU重复相加。\n\n')
        f.write('## 物理求解与辅助动作分别计数\n\n|设备|物理矩阵维数|成功child|原始物理LU|原始物理复RHS|\n|---|---:|---:|---:|---:|\n')
        for d in dimensions:f.write(f"|{d['host']}|{d['physical_matrix_dimension']}×{d['physical_matrix_dimension']}|{d['successful_children']}|{d['physical_LU']}|{d['physical_complex_RHS']}|\n")
        f.write('\nscalar矩阵为576×576；vector矩阵为648×648，对应216个共享实材料cell、每cell三个场分量。两种物理与尺寸不混为同等LU工作量。\n\n')
        if dr:
            f.write('|DR任务数|固定背景H Cholesky|H setup RHS|H iteration RHS|\n|---:|---:|---:|---:|\n')
            f.write(f"|{len(dr)}|{dr_sum('H_cholesky')}|{dr_sum('H_setup_rhs')}|{dr_sum('H_iteration_rhs')}|\n")
        f.write('\nDR的固定背景H分解与迭代RHS单独保存，既不加到材料相关物理LU，也不加到物理Maxwell RHS；其CPU已包含在实际child计时。material-force-decomposition的56分支含152064个block复RHS，属于缓存4cell辅助动作，**不是物理solver RHS**；其LAUNCH外部CPU仅计一次，内层receipt与plot body为嵌套计时。\n\n')
        f.write('## 费用与来源边界\n\n当前RECEIPT/ACTUAL_CONFIG和终态队列足以核对576个solver费用。新independent264的完整NPZ原档仍在远端发布链路，本地缺NPZ只表示传输/compactreview范围，不是科学缺失。shared旧失败root与成功attempt_02分开保留；vector使用attempt_01。\n\n')
        f.write('exposed48_v1的分析输出重复pair错误已保留其20.179077秒分析CPU，不计为物理失败；v1/v2数值比较由主线程验收。初始import、输出格式失败及未测传输费用的CPU仍保留null，不补0。controller、generation、validation、scoring与archive/transport在详细账本单列；不把solver、wrapper、body与嵌套事件重复相加。最终上传/出版费用由主线程后续追加，本表不等待该阶段、也不把目前已知费用称为最终全项目总成本。\n\n')
        f.write('证据文件：['+out.name+'/REPORT.md]('+str(out/'REPORT.md')+') · [按campaign与method分组]('+str(out/'SOLVER_METHOD_TOTALS.csv')+') · [逐child费用]('+str(out/'SOLVER_TASK_COSTS.csv')+') · [DR背景费用]('+str(out/'DR_BACKGROUND_COSTS.csv')+') · [576进程核对]('+str(out/'SOLVER_COMPLETION_VERIFICATION.json')+')。本次只读取成本/状态metadata，未读GT数组、未运行新physics或改变算法。\n')


def solver_completion(root,out):
    """Actual successful solver-stage OSCPU only;not controller/body double totals."""
    specs=[('dualmetric_pilot12_v1',None,48,'Windows',576),
        ('dualmetric_remaining36_v1',None,144,'Windows',576),
        ('self_only16_v1',None,48,'Windows',576),('stage_control16_v1',None,32,'Windows',576),
        ('independent24_v1',None,264,'Windows',576),
        ('shared_readout4_v1','attempt_02',16,'Apple_local',576),
        ('vector4_v1','attempt_01',24,'Apple_local',648)]
    tasks=[];campaigns=[];dr=[];methods=[]
    for campaign,attempt,expected,host,dimension in specs:
        qp=root/'queues'/campaign
        if attempt:qp=qp/attempt
        qp=qp/'RECEIPT.json';q=read(qp)
        if not isinstance(q,dict):raise FileNotFoundError('Completed solver queue missing: '+str(qp))
        jobs=q.get('jobs',[])
        if len(jobs)!=expected or any(j.get('returncode')!=0 or j.get('child_status') in ('FAILED',None) for j in jobs):
            raise ValueError('Not all expected solver children are successful: '+campaign)
        total=q.get('finished_child_CPU_seconds',q.get('solver_children_whole_OS_CPU_seconds',q.get('whole_children_OS_CPU_seconds')))
        if not finite(total):raise ValueError('Actual completed childOSCPU required: '+campaign)
        collected=[]
        for j in jobs:
            folder=root/'RESULTS'/campaign/j['task_id'];rp=folder/attempt/'RECEIPT.json' if campaign=='vector4_v1' else folder/'RECEIPT.json'
            d=read(rp)
            if not isinstance(d,dict):raise FileNotFoundError('Local cost-only child receipt required,not NPZ: '+str(rp))
            task=d.get('task',{});method=task.get('method',j.get('method',j['task_id'].split('__')[1]))
            costs=d.get('physical_cost',d.get('operation_counts_all',{}));lu,rhs=physical(costs)
            external=(j.get('external') or {}).get('child_process_CPU_seconds',j.get('child_process_CPU_seconds'))
            body=d.get('CPU_seconds',d.get('whole_process_CPU_seconds'))
            item=dict(campaign=campaign,method=method,task_id=j['task_id'],case_id=task.get('case_id'),host=host,
                physical_matrix_dimension=dimension,attempt=attempt or 'root',queue_receipt=str(qp),child_receipt=str(rp),
                run_status=d.get('status'),returncode=j.get('returncode'),actual_child_OS_CPU_seconds=external if finite(external) else None,
                OS_CPU_null_reason=None if finite(external) else 'Only whole-campaign actual childOSCPU collected;per-task/per-method attribution unavailable andnot proportionallyallocated',
                internal_body_CPU_seconds=body,physical_LU=lu,physical_complex_RHS=rhs,
                projector_cost=d.get('projector_cost'),auxiliary_cost=d.get('auxiliary_cost',d.get('auxiliary_receipt')),
                cost_scope='Successful solverstage only;all setup/proposal/corrector/finalserialization in externalchildtime;noGT/NPZ needed',
                local_fullfields_available=any(folder.rglob('*.npz')),
                fullfield_scope='LocalNPZ absence is atransport/compactreview boundary,not scientific missingness')
            tasks.append(item);collected.append(item)
            if isinstance(d.get('projector_cost'),dict):
                dr.append(dict(campaign=campaign,method=method,task_id=j['task_id'],case_id=task.get('case_id'),host=host,
                    background_matrix_dimension=dimension,physical_LU=lu,physical_complex_RHS=rhs,
                    **d['projector_cost'],scope='Fixedbackground Hfactor/setup/iteration andpropagation actions;NOT physical Maxwell LU/RHS;already inside childCPU'))
        if host=='Windows':
            if not all(finite(t['actual_child_OS_CPU_seconds']) for t in collected) or abs(sum(t['actual_child_OS_CPU_seconds'] for t in collected)-total)>1e-6:
                raise ValueError('Windows aggregate andper-child GetProcessTimes mismatch: '+campaign)
        l=[t['physical_LU'] for t in collected];h=[t['physical_complex_RHS'] for t in collected]
        campaigns.append(dict(campaign=campaign,attempt=attempt or 'root',host=host,physical_matrix_dimension=dimension,
            successful_children=len(collected),expected_successful_children=expected,actual_children_OS_CPU_seconds=total,
            controller_CPU_seconds=q.get('controller_CPU_seconds'),
            queue_wall_seconds=q.get('queue_wall_seconds',q.get('controller_wall_seconds',q.get('controller_wall_including_source_freeze_seconds',q.get('wall_seconds')))),
            wall_scope='This queue only;never summed into an overall parallel-walltime',
            physical_LU=sum(l) if all(finite(v) for v in l) else None,
            physical_complex_RHS=sum(h) if all(finite(v) for v in h) else None,
            body_CPU_seconds_sum=sum(t['internal_body_CPU_seconds'] for t in collected) if all(finite(t['internal_body_CPU_seconds']) for t in collected) else None,
            per_method_actual_OS_available=host=='Windows',queue_receipt=str(qp)))
        for method in sorted({t['method'] for t in collected}):
            g=[t for t in collected if t['method']==method];osvals=[t['actual_child_OS_CPU_seconds'] for t in g]
            methods.append(dict(campaign=campaign,method=method,host=host,physical_matrix_dimension=dimension,successful_children=len(g),
                actual_children_OS_CPU_seconds=sum(osvals) if all(finite(v) for v in osvals) else None,
                OS_CPU_null_reason=None if all(finite(v) for v in osvals) else 'Apple batch OS total known;permethod OS not measured. BodyCPU retained separately,not substituted.',
                internal_body_CPU_seconds_sum=sum(t['internal_body_CPU_seconds'] for t in g),
                physical_LU=sum(t['physical_LU'] for t in g) if all(finite(t['physical_LU']) for t in g) else None,
                physical_complex_RHS=sum(t['physical_complex_RHS'] for t in g) if all(finite(t['physical_complex_RHS']) for t in g) else None))
    count=sum(c['successful_children'] for c in campaigns)
    if count!=576:raise ValueError('Expected576 successful solver processes')
    independent=next(c for c in campaigns if c['campaign']=='independent24_v1')
    shared=next(c for c in campaigns if c['campaign']=='shared_readout4_v1')
    vector=next(c for c in campaigns if c['campaign']=='vector4_v1')
    for observed,target in [(independent['actual_children_OS_CPU_seconds'],30369.53125),
        (shared['actual_children_OS_CPU_seconds'],952.745),(vector['actual_children_OS_CPU_seconds'],297.930588)]:
        if abs(observed-target)>1e-6:raise ValueError('Requested actual OSCPU checkpoint differs from source')
    failurepath=root/'queues/shared_readout4_v1/RECEIPT.json';failure=read(failurepath)
    if not failure or len(failure.get('jobs',[]))!=16 or any(j.get('returncode')==0 for j in failure['jobs']):raise ValueError('Original16failedattempts must be retained separately')
    failed_CPU=failure.get('solver_children_whole_OS_CPU_seconds')
    if not finite(failed_CPU) or abs(failed_CPU-.514267)>1e-8:raise ValueError('Failed attempt actualCPU differs')
    dimensions=[]
    for host,dim in sorted({(t['host'],t['physical_matrix_dimension']) for t in tasks}):
        g=[t for t in tasks if t['host']==host and t['physical_matrix_dimension']==dim]
        dimensions.append(dict(host=host,physical_matrix_dimension=dim,successful_children=len(g),
            physical_LU=sum(t['physical_LU'] for t in g),physical_complex_RHS=sum(t['physical_complex_RHS'] for t in g),
            scope='Original material-dependent physical Maxwell solves only;DRbackground andblockRHS not included'))
    actual=sum(c['actual_children_OS_CPU_seconds'] for c in campaigns)
    verified=dict(status='ALL576_SUCCESSFUL_SOLVERS_ACCOUNTED',successful_solver_children=576,
        failed_extra_attempts=16,successful_children_actual_OS_CPU_seconds=actual,
        Windows_success_CPU_seconds=sum(c['actual_children_OS_CPU_seconds'] for c in campaigns if c['host']=='Windows'),
        Apple_success_CPU_seconds=sum(c['actual_children_OS_CPU_seconds'] for c in campaigns if c['host']=='Apple_local'),
        failed_extra_actual_OS_CPU_seconds=failed_CPU,all_success_and_failed_solver_OS_CPU_seconds=actual+failed_CPU,
        independent264_OS_CPU_seconds=independent['actual_children_OS_CPU_seconds'],
        shared16_OS_CPU_seconds=shared['actual_children_OS_CPU_seconds'],vector24_OS_CPU_seconds=vector['actual_children_OS_CPU_seconds'],
        no_parallel_wall_sum=True,no_GT_read=True,no_new_physics=True,
        local_fullfields_not_required_for_cost_closure=True,
        publication_transport_cost_scope='Costs present locally atledger sampling only;final upload/publication cost to be appended by parent')
    for name,values in [('SOLVER_TASK_COSTS',tasks),('SOLVER_CAMPAIGN_TOTALS',campaigns),('SOLVER_METHOD_TOTALS',methods),
        ('SOLVER_PHYSICAL_DIMENSION_TOTALS',dimensions),('DR_BACKGROUND_COSTS',dr)]:table(out,name,values)
    write(out/'SOLVER_COMPLETION_VERIFICATION.json',verified)
    write_solver_root_report(root,out,verified,campaigns,dimensions,dr)
    return verified


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--solver-complete',action='store_true')
    args=parser.parse_args();root=args.root.resolve();out=args.out.resolve();cpu0,wall0=time.process_time(),time.perf_counter()
    if out.exists() or out==root or any(out.is_relative_to(root/n) for n in ('RESULTS','sealed_truth','data_public','src','configs','queues')):
        raise ValueError('New separate ledger output directory required')
    out.mkdir(parents=True,exist_ok=False);started=datetime.now(timezone.utc).isoformat()
    queue_costs(root);result_costs(root);intervention_costs(root);metadata_costs(root,out)
    availability=coverage(root);checkpoints,typed_fields=checkpoint_catalog(root)
    solver_verified=solver_completion(root,out) if args.solver_complete else None
    totals=[]
    for host,category in sorted({(r['host'],r['category']) for r in ROWS}):
        group=[r for r in ROWS if r['host']==host and r['category']==category]
        canonical=[r for r in group if r['include_CPU_in_additive_total']]
        physical_rows=[r for r in group if r['include_physics_in_additive_total']]
        totals.append(dict(host=host,category=category,known_CPU_subtotal=sum(r['CPU_seconds'] for r in canonical if finite(r['CPU_seconds'])),
            CPU_complete=bool(canonical) and all(finite(r['CPU_seconds']) for r in canonical),excluded_only_category=not canonical,unknown_CPU_rows=[r['cost_id'] for r in canonical if r['CPU_seconds'] is None],
            known_original_LU_subtotal=sum(r['physical_LU'] for r in physical_rows if finite(r['physical_LU'])),
            known_original_complex_RHS_subtotal=sum(r['physical_complex_RHS'] for r in physical_rows if finite(r['physical_complex_RHS'])),
            unknown_physics_rows=[r['cost_id'] for r in physical_rows if r['physical_LU'] is None or r['physical_complex_RHS'] is None],
            failed_or_incomplete_rows=sum(r['failed_or_incomplete'] for r in canonical),
            cost_scope='Known measured subtotal only;unknown not0,wall notCPU,crosshost notspeed ratio'))
    for name,values in [('COST_LEDGER',ROWS),('COST_TOTALS',totals),('COVERAGE',availability),('FULL_STATE_CATALOG',checkpoints),('TYPED_STATE_FIELDS',typed_fields),('SOURCE_INPUTS',INPUTS),('ISSUES',ISSUES)]:table(out,name,values)
    catalog_groups=defaultdict(list)
    for cp in checkpoints:
        rel=Path(cp['path']).relative_to(root);group='/'.join(rel.parts[:2]);catalog_groups[group].append(cp)
    with (out/'FULL_STATE_CATALOG.md').open('x') as f:
        f.write('# Full-state catalog for theory2 / Pro\n\nOnly typed JSONmetadata andNPYheaders were inspected. No m/J/Q/A payload,GT,forward or optimizer was read/run. CSV/JSON preserve every path andtyped field.\n\n|Local source group|Files|Kinds|Material dimensions|Source counts|Known stage rows|\n|---|---:|---|---|---|---:|\n')
        for group,items in sorted(catalog_groups.items()):
            kinds=sorted({str(x.get('kind')) for x in items});dims=sorted({str(x.get('material_real_dimension')) for x in items});sources=sorted({str(x.get('sources')) for x in items})
            f.write(f"|{group}|{len(items)}|{','.join(kinds)}|{','.join(dims)}|{','.join(sources)}|{sum(x.get('stage') is not None for x in items)}|\n")
        f.write('\nScalar state: m(realN),A(complexnativecoefficients×sources),J/Q(complexN×sources),oldprecurrentgradient/direction,DPCS/source metrics,actual4cell/Green caches andstage loopmetadata. A count is ambient,not independent current rank. Vector state is cell-majorXYZ withone real material percell,vectorrangebasis/prstate,andrequires thecompatibleexistingmodel.\n\nA full auxiliary iterator checkpoint isnot a whole-process restart: physicalLU cache,timers,controller/selector state andruntimeprior objects areoutside thisstructural certificate. PRmemory explicitNone may be a legalstageentry/reset;never silently substitute zero memory. Use SAME fullstate path for paired interventions;material-only endpoints areinsufficient. Cache presence isstructural,not a recomputed physics/identity check. Missing stage/pair metadata remainsunknown;path-derived pair labels areexplicitlyqualified.\n')
    receipt=dict(status='SOLVER_COMPLETE_PUBLICATION_COSTS_OPEN' if solver_verified else 'PARTIAL_LOCAL_EVIDENCE_LEDGER',solver_verification=solver_verified,started_utc=started,finished_utc=datetime.now(timezone.utc).isoformat(),
        cost_rows=len(ROWS),checkpoint_rows=len(checkpoints),typed_fields=len(typed_fields),issues=len(ISSUES),
        root=str(root),new_physics=0,new_optimizer=0,GT_arrays_read=False,NPZ_numeric_payloads_read=False,
        new_CPU_seconds=time.process_time()-cpu0,new_wall_seconds=time.perf_counter()-wall0,
        own_cost_scope='Ledger construction cost separate from input experiment costs;final report/receiptwrite/exit excluded',
        reuse='Copied old public inputs/checkpoints are not recharged as new generation;current replay/validation work is charged',
        no_SHA=True,finalization='Rerun with new --out after complete receipts arrive;never overwrite partial ledger')
    write(out/'RECEIPT.json',receipt)
    notes=['# R14 mechanical cost ledger — partial local evidence',
        'CPU totals are separated byApple/Windows andsetup/online/validation/scoring/transport. Only include_CPU_in_additive_total rows enter measured CPU subtotals. Original physical LU/RHS use a separate inclusion flag. Wall times are retained but never summed into process CPU. Missing cost fields remain null.',
        'Windows per-child GetProcessTimes dominates the corresponding solver internal timer;controllerCPU is disjoint. Where only an actual children aggregate is available,the aggregate is counted once andper-child body records remain nested. Modes andper-trajectory diagnostics inside one Python process are nested under theirdriver. Generation case/setup timers are nested under itsprocess total. Failed attempts remain billable.',
        'Vector OS CPU prefers its real queues/vector4_v1/attempt_01 controller receipt when local;the derived COST_TOTALS CPU remains nested while its physical totals are retained once. Shared-readout failed root attempt and successful attempt_02 are distinct cost authorities;matching uses receipt status and bodyCPU,not taskID alone. allsubprocess−solverchildren remainder is added only once from the preferred source. Intervention batch childCPU replaces bodyCPU without claiming an undocumentedOS collector.',
        'Uploaded/runtime ZIP presence establishesbytes only. Network/packCPU without a realreceipt remainsunknown,and mayoverlap a measuredarchive process. No guessed currency/GPU/peakmemory/time is added. Whole-process timer samples generallyexclude finalreceiptwrite/exit unless anexternal collector saysotherwise.',
        'COVERAGE recordslocallyavailable expectedchildreceipts,not a claim thatmissingjobs neverran. dualmetric48_v1 is a planning supersetofthepilot/remainder,not another192executedtasks. Partial/failure/unknown rows are not zero-cost or negative scientific findings.',
        'FULL_STATE_CATALOG andTYPED_STATE_FIELDS use onlyNPZ headers andtypedJSONmetadata. Huge historytrees aretruncated afterthe precedingstage/PRloopmetadata. No material/current arrayvalues were read. Full auxiliary state isnot a complete process-clock/selector/physicalLU restart,and a material-only endpoint isnot PRmemory.',
        '理论2/Pro使用：先按同一个checkpoint路径定位m/A/J/Q、旧precurrent梯度/方向、stage/iteration、DPCS/source metrics、实际Green/4cellcache。A的ambient维数不是独立current rank。scalar native Fourier和3D vector range坐标分开；PR memory显式None可能是合法stage入口，不能当成丢失后任意补零。正式扰动必须复用同一完整state，而不是只复用m；本目录不执行扰动、不作理论或科学验收。',
        'All costs aremechanical evidence. Parent retainsfinal scientific review.']
    with (out/'REPORT.md').open('x') as f:
        f.write('\n\n'.join(notes)+'\n\n|Host|Category|Known CPU s|CPU complete|\n|---|---|---:|---|\n')
        for r in totals:f.write(f"|{r['host']}|{r['category']}|{r['known_CPU_subtotal']:.6f}|{r['CPU_complete']}|\n")
        f.write(f"\nCheckpoint files:{len(checkpoints)};typed field rows:{len(typed_fields)};issues:{len(ISSUES)}.\n")
    print(json.dumps(receipt));return 0


if __name__=='__main__':raise SystemExit(main())
