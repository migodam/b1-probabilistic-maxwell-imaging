"""Metadata-only known-cost inventory. No score values, NPZ, truths or solvers.

Final queue jobs own external CPU; internal timers/cumulative analysis totals are
references. All unknowns and overlapping summaries remain explicit. Regenerate
into a NEW output directory after later queues/downloads complete.
"""
import argparse,csv,json,time
from pathlib import Path
from datetime import datetime,timezone


def load(path):return json.loads(Path(path).read_text())


def canonical(counter):
    lu=counter.get('physical_LU',counter.get('LU',counter.get('new_LU',counter.get('new_Maxwell_LU'))))
    if 'actual_complex_RHS' in counter:rhs=counter['actual_complex_RHS'];keys=['actual_complex_RHS']
    elif 'total_complex_RHS' in counter:rhs=counter['total_complex_RHS'];keys=['total_complex_RHS']
    elif 'physical_complex_RHS' in counter:rhs=counter['physical_complex_RHS'];keys=['physical_complex_RHS']
    elif any(k in counter for k in ['source_RHS','adjoint_RHS','tangent_RHS','receiver_transpose_RHS']):
        keys=[k for k in ['source_RHS','adjoint_RHS','tangent_RHS','receiver_transpose_RHS'] if k in counter];rhs=sum(counter[k] for k in keys)
    else:rhs=counter.get('RHS',counter.get('new_RHS'));keys=['RHS'] if 'RHS' in counter else ['new_RHS'] if 'new_RHS' in counter else []
    return lu,rhs,keys


def supplemental_receipts(root,add):
    """Read only explicit cost receipts, never result tables/material payloads."""
    zero=dict(LU=0,RHS=0)
    final_metadata=[(root/'audit/DOCUMENT_ASSEMBLY_RECEIPT.json','report_assembly'),(root/'audit/RUN_CONTRACTS_FINAL/BUILD_RECEIPT.json','metadata_contracts')]
    final_metadata += [(p,'publication_inventory') for p in (root/'audit').glob('GIT_PUBLICATION_*/GIT_ALLOWLIST.json')]
    final_metadata += [(p,'publication_link_check') for p in (root/'audit').glob('GIT_PUBLICATION_*/LINK_CHECK.json')]
    final_metadata += [(p,'publication_readback') for p in (root/'audit').glob('PUBLIC_READBACK*.json')]
    for f,phase in final_metadata:
        if not f.is_file():continue
        d=load(f);add('final_metadata:'+str(f.relative_to(root)),phase,d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=zero,scope='One actual metadata/document/publication verification process; nested records not charged again')
    extraction=list((root/'audit').glob('*EXTRACTION*RECEIPT.json'))
    extraction+=[root/'audit/MECHANISM_AND_RESCUE_INPUT_EXTRACTION.json']
    for f in extraction:
        if not f.is_file():continue
        d=load(f);items=d if isinstance(d,list) else[d]
        for i,item in enumerate(items):add('extraction:'+f.name+':'+str(i),'release_extraction',item.get('host','Apple_local'),f,f'[{i}]' if isinstance(d,list) else'',cpu=item.get('CPU_seconds'),wall=item.get('wall_seconds'),counter=zero,status=item.get('status'),reason='Batch entries each counted once; archive/file byte fields are not timing',scope='Actual extraction interval; no source/GT array payload opened by inventory')
    f=root/'audit/PHYSICS_DEPENDENCY_RECOVERY_PLAN.json'
    if f.is_file():
        d=load(f);add('dependency_runtime_pack','setup_packing',d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=zero,scope='Dependency-only runtime bundle pack; named bundle suppresses unknown-runtime duplicate')
        failed=d.get('first_metadata_preparation_attempt')
        if failed:add('dependency_pack_failed_metadata','setup_failed','Apple_local',f,'first_metadata_preparation_attempt',cpu=failed.get('CPU_seconds'),counter=zero,reason='Original failed metadata preparation timing unknown where null')
    f=root/'audit/PHYSICS_C4_PREPARATION_RECEIPT.json'
    if f.is_file():
        d=load(f)
        for field in ['controller_CPU_seconds','transport_child_CPU_seconds']:add('physics_C4_prepare:'+field,'setup_transport','Apple_local',f,field,cpu=d.get(field),wall=d.get('wall_seconds') if field=='controller_CPU_seconds' else None,counter=zero,status=d.get('status'),scope='Separate controller and transport-child CPU; parent wall recorded once, never summed child walls')
        rem=d.get('remote_preparation',{});add('physics_C4_remote_prepare','setup_packing','i9_remote',f,'remote_preparation',cpu=rem.get('CPU_seconds'),wall=rem.get('wall_seconds'),counter=zero,status=rem.get('status'),scope='Remote dependency/input preparation; separate host from transport controller')
        failed=d.get('initial_failed_transport')
        if failed:add('physics_C4_failed_transport','setup_failed','Apple_local',f,'initial_failed_transport',cpu=failed.get('CPU_seconds'),counter=zero,reason='Failed command-length transport before remote Python; original cost unknown, not zero')
    for f in (root/'audit').glob('*RECOVERY_COMPOSITION*.json'):
        d=load(f);add('recovery_composition:'+f.name,'release_composition',d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=zero,status=d.get('status'),scope='Byte-copy and receipt composition only; source attempts not charged again')
    for f in root.glob('release_preparation_*/PACKAGING_RECEIPT.json'):
        d=load(f);add('release_packaging:'+f.parent.name,'release_packaging',d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=zero,status=d.get('status'),scope='Whole packaging interval only; nested per-archive/split timers excluded')
    sources=[]
    for f in (root/'audit').glob('RELEASE_UPLOAD_COSTS*.jsonl'):
        for i,line in enumerate(f.read_text().splitlines()):
            if line.strip():sources.append((f,i,json.loads(line),'release_upload'))
    f=root/'audit/ARTIFACT_TRANSFER_COSTS.json'
    if f.is_file():
        d=load(f);items=d if isinstance(d,list) else d.get('transfers',[])
        for i,item in enumerate(items):sources.append((f,i,item,'release_download'))
    for f,i,item,phase in sources:
        host=item.get('host','Apple_local')
        for field in ['controller_CPU_seconds','transport_child_CPU_seconds']:add('transfer:'+f.name+':'+str(i)+':'+field,phase,host,f,f'[{i}].'+field,cpu=item.get(field),wall=item.get('wall_seconds') if field=='controller_CPU_seconds' else None,counter=zero,status=item.get('status'),reason='Only sanitized transfer timing; public download index is not a cost receipt',scope='Controller/transport child CPU separate; one observed transfer wall, no child-wall sum')
    for f in (root/'audit').rglob('AUDIT.json'):
        if not any('PUBLIC_ACCEPTANCE' in part for part in f.parts):continue
        d=load(f);cost=d.get('validation_cost',{})
        add('public_acceptance:'+str(f.relative_to(root)),'validation_public_interfaces',d.get('host','Apple_local'),f,'validation_cost',cpu=cost.get('CPU_seconds'),wall=cost.get('wall_seconds'),counter=dict(LU=cost.get('new_Maxwell_factorizations'),RHS=0 if cost.get('new_Maxwell_factorizations')==0 else None,matvec_actions_reference=cost.get('matvec_actions')),status=d.get('status'),reason='Every actual audit attempt batch charged, including failed WindowsPath/partial wiring audits; task rows not added',scope=cost.get('timer_scope','Outer public saved-state audit interval'))
    for name,phase in [('CROSS_CAMPAIGN_PUBLIC_INTERFACES.json','validation_public_interfaces'),('ARCHIVE_METADATA_REPACK_ATTEMPT.json','release_packaging_failed_attempt'),('PHYSICS_RECOVERY_PROMOTION_RECEIPT.json','release_composition')]:
        f=root/'audit'/name
        if not f.is_file():continue
        d=load(f);add('additional_receipt:'+name,phase,d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=zero,status=d.get('status'),scope='Actual independent metadata/interface/copy/packing attempt interval; later successful archive compression remains separate')
    # Receipt schema from evaluate_measured_holdout.py: totals own unique
    # material LU/RHS, not allocated per-view shares or setup subtimers.
    for f in (root/'analysis').rglob('*RECEIPT.json'):
        if any('fixture' in part.lower() or 'fake' in part.lower() for part in f.parts):continue
        if f.name in ['SCORING_RECEIPT.json','PLOT_RECEIPT.json','AGGREGATION_RECEIPT.json']:continue
        d=load(f)
        if 'unique_physical_LU' in d and 'total_evaluator_CPU_seconds' in d:
            add('heldout_evaluator:'+str(f.relative_to(root)),'offline_holdout_validation',d.get('host','Apple_local'),f,cpu=d.get('total_evaluator_CPU_seconds'),wall=d.get('total_evaluator_wall_seconds'),counter=dict(LU=d.get('unique_physical_LU'),RHS=d.get('unique_source_RHS')),dimension=576,status=d.get('status'),scope='Unique physical totals once; allocated view shares, forward/setup/import/guard subtimers and near-exit process timer all nested')
        elif 'physical_LU_known_sum' in d and 'physical_complex_RHS_known_sum' in d:
            add('stationarity:'+str(f.relative_to(root)),'offline_stationarity',d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=d.get('physical_LU_known_sum'),actual_complex_RHS=d.get('physical_complex_RHS_known_sum'),forward_RHS_reference=d.get('physical_forward_RHS_known_sum'),adjoint_RHS_reference=d.get('physical_adjoint_RHS_known_sum')),dimension=576,status=d.get('status'),scope='Deduplicated unique-material receipt sums once; requested400/evaluated rows and forward/adjoint subtotals are references, never additional solves')
        elif f.name=='REVIEW_TABLES_RECEIPT.json':
            add('review_table_serialization:'+str(f.relative_to(root)),'scoring_plotting',d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=zero,status=d.get('status'),scope='Compact table serialization outer CPU once; source aggregation and online costs not repeated')
        elif f.name=='PLOT_COST_RECEIPT.json':
            add('vector_plot:'+str(f.relative_to(root)),'scoring_plotting',d.get('host','Apple_local'),f,cpu=d.get('total_plot_CPU_seconds'),wall=d.get('total_plot_wall_seconds'),counter=dict(LU=d.get('new_LU'),RHS=d.get('new_RHS')),status=d.get('status'),scope='Whole vector plot process once; imports and per-figure CPU/wall sums are nested, not additive')
        elif f.name=='RECEIPT.json' and any(key in str(f.parent).lower() for key in ['vector','objective_truth','feasible','oracle','stationarity','portfolio_diagnostic']):
            no_forward=any(d.get(k)==0 for k in ['new_forward_solves','new_forward','new_Maxwell_LU','new_LU'])
            counter=d.get('physical_cost',d.get('physical_costs'))
            if counter is None:counter=zero if no_forward else None
            add('offline_diagnostic:'+str(f.relative_to(root)),'scoring_plotting',d.get('host','Apple_local'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=counter,status=d.get('status'),reason='Receipt outer timer only; absent CPU/physical counts remain unknown',scope='Offline figure/objective-truth/feasible-pair/oracle/stationarity receipt, no rowwise duplicates')


def main():
    started,wall0=time.process_time(),time.perf_counter();p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();root=a.root.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);rows=[];notes=[];seen=set()
    def add(identifier,phase,host,source,pointer='',cpu=None,wall=None,counter=None,dimension=None,cpu_include=True,physics_include=True,status=None,reason=None,scope=None):
        if identifier in seen:return
        seen.add(identifier);lu,rhs,keys=canonical(counter or{})
        rows.append(dict(cost_id=identifier,phase=phase,host=host,source=str(source.relative_to(root)) if isinstance(source,Path) and source.is_relative_to(root) else str(source),json_pointer=pointer,status=status,CPU_seconds=cpu,wall_seconds=wall,CPU_known_sum_include=cpu_include,physical_LU=lu,physical_complex_RHS=rhs,physical_known_sum_include=physics_include,system_dimension=dimension,physical_RHS_source_keys=keys,counter_metadata=counter,missing_or_overlap_reason=reason,timer_scope=scope))
    # Queues: merged final campaigns dominate component recovery receipts.
    queuedirs=list((root/'queues').iterdir()) if(root/'queues').is_dir() else[];selected={}
    for q in queuedirs:
        if not q.is_dir():continue
        final=q/'RECEIPT.json';candidates=sorted(q.glob('STATUS_*.json'))
        path=final if final.is_file() else candidates[-1] if candidates else None
        if path is None:continue
        try:d=load(path)
        except Exception as ex:notes.append(dict(source=str(path),issue='Unreadable snapshot',error=repr(ex)));continue
        selected[q.name]=(path,d)
    recovery_sources=set()
    for name,(path,d) in selected.items():
        if d.get('recovery_source'):recovery_sources.add((path.parent/str(d['recovery_source']).replace('\\','/')).resolve())
    for name,(path,d) in selected.items():
        host='Apple_local' if 'local' in name else'i9_remote';recovery_reference=path.resolve() in recovery_sources
        add('queue:'+name,'online',host,path,cpu=None,wall=d.get('wall_seconds'),cpu_include=False,physics_include=False,status=d.get('status'),reason='Merged parent queue already includes these jobs/controller window; reference only' if recovery_reference else'Queue elapsed is an independent interval; never sum child wall times.',scope=d.get('wall_scope','Parallel queue elapsed wall; not summed child wall'))
        add('controller:'+name,'online_controller',host,path,'controller_CPU_seconds',cpu=d.get('controller_CPU_seconds'),cpu_include=not recovery_reference,physics_include=False,status=d.get('status'),reason='Included in merged parent controller total' if recovery_reference else'Controller CPU separate from child CPU; partial status remains partial')
        for i,job in enumerate(d.get('jobs',[])):
            tid=job.get('task_id',f'job{i}');ext=job.get('external',{}) or{};cpu=ext.get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'));counter=None;dimension=None;receipt=None
            rp=root/'RESULTS'/name/tid/'RECEIPT.json'
            if rp.is_file():
                try:
                    receipt=load(rp);counter=receipt.get('physical_cost');components=receipt.get('components');dof=receipt.get('material_DOF');dimension=components*dof if isinstance(components,int) and isinstance(dof,int) else None
                except Exception as ex:notes.append(dict(source=str(rp),issue='Receipt cost unavailable',error=repr(ex)))
            add('job:'+name+':'+tid,'online',host,path,f'jobs[{i}].external',cpu=cpu,wall=job.get('wall_seconds'),counter=counter,dimension=dimension,cpu_include=not recovery_reference,physics_include=not recovery_reference,status=job.get('child_status'),reason='Already counted in merged parent' if recovery_reference else 'External CPU unknown; never use internal timer as exit CPU' if cpu is None else None,scope=ext.get('CPU_scope','Final queue child exit accounting'))
            if receipt is not None:
                add('internal:'+name+':'+tid,'online_reference',host,rp,'whole_process_CPU_seconds',cpu=receipt.get('whole_process_CPU_seconds'),counter=counter,dimension=dimension,cpu_include=False,physics_include=False,reason='Internal pre-exit timer/counters overlap authoritative job row; do not add')
                for field in ['projector_cost']:
                    if receipt.get(field):add('operator:'+name+':'+tid+':'+field,'online_operator_reference',host,rp,field,counter=receipt[field],dimension=dimension,cpu_include=False,physics_include=False,reason='Operator/H/port actions are distinct from physical solves; nested in job CPU')
        # Composition retains replaced technical failures as additional attempts.
        # Recovery jobs/controller are already in the merged parent; failures are not.
        for i,job in enumerate(d.get('failed_attempts',[])):
            ext=job.get('external') or {};tid=job.get('task_id',f'failed{i}');rp=root/str(job.get('result_receipt_relative',''))
            receipt=load(rp) if rp.is_file() else {};counter=receipt.get('physical_cost')
            add('failed_attempt:'+name+':'+str(i)+':'+tid,'online_failed_attempt',host,path,f'failed_attempts[{i}]',cpu=ext.get('child_process_CPU_seconds',job.get('child_process_CPU_seconds')),wall=job.get('wall_seconds'),counter=counter,cpu_include=not recovery_reference,physics_include=not recovery_reference,status=job.get('child_status'),reason='Separate original technical failure, not present in final primary jobs; missing timing/counters remain null',scope=ext.get('CPU_scope','Original child exit accounting; no failure cost inferred from successful retry'))
        for item in d.get('running',[]):add('running:'+name+':'+item['task_id'],'online_incomplete',host,path,'running',cpu=None,counter=None,status='RUNNING_AT_SNAPSHOT',reason='Unfinished child exit CPU/physical totals not available; not zero')
    # Expected manifests with no local queue evidence remain visible, including remote work.
    for f in (root/'configs').glob('*.json'):
        try:d=load(f)
        except Exception:continue
        name=d.get('campaign')
        if name and isinstance(d.get('tasks'),list) and name not in selected:
            add('unobserved_queue:'+name,'online_incomplete','Apple_local' if 'local' in name else'i9_remote',f,'tasks',cpu=None,counter=None,status='NO_LOCAL_QUEUE_RECEIPT',reason=f"Expected{len(d['tasks'])}tasks; local ledger cannot certify remote/current completion")
    # Analysis cost totals are a cross-check only, never another consumption total.
    for f in (root/'analysis').glob('*/COST_TOTALS.json'):
        d=load(f);name=f.parent.name
        add('analysis_cost_reference:'+name,'online_reference','Apple_local' if 'local' in name else'i9_remote',f,cpu=d.get('external_child_CPU_seconds_known_sum'),counter=dict(LU=d.get('physical_LU_known_sum'),actual_complex_RHS=d.get('physical_complex_RHS_known_sum')),cpu_include=False,physics_include=False,reason='Same jobs as queue; analysis totals cross-check, not additional work')
    recovery=root/'audit/CONTROLLER_FAILURE_AND_RECOVERY.json'
    if recovery.is_file():
        for i,item in enumerate(load(recovery).get('partial_attempts',[])):add('aborted_C1:'+item['task_id'],'online_failed','i9_remote',recovery,f'partial_attempts[{i}]',cpu=item.get('CPU_seconds'),counter=None,status='ABORTED_ATTEMPT',reason=item.get('reason','Lost process timing'))
    # Explicit missing-timer failures disclosed by parent; preserve rather than infer.
    add('local_psutil_import_failure','validation_failed','Apple_local','Parent orchestration disclosure: initial local psutil import failure',cpu=None,counter=dict(LU=0,RHS=0),status='FAILED_BEFORE_MODEL',reason='Original CPU/wall was not captured; no recoverable timing record')
    failed=root/'delegated/full_current/dual_pressure_diagnostic_v2/FAILED_ATTEMPT_ACCOUNTING.json'
    if failed.is_file():
        d=load(failed);add('dual_pressure_failed_serialization','validation_failed','Apple_local',failed,cpu=d.get('CPU_seconds'),counter=d.get('physical_cost_known_from_completed_calls'),dimension=576,status=d.get('status'),reason=d.get('CPU_scope'),scope='Known completed3LU/84RHS before JSON serialization failed; separate failed attempt')
    # Fresh generation and NN preprocessing: outer timer includes case/subsetup timings.
    for rel in ['delegated/evidence/generation/GENERATION_RECEIPT.json','delegated/evidence/nnprep/PREPROCESS_RECEIPT.json']:
        f=root/rel
        if not f.is_file():continue
        d=load(f);t=d.get('totals',{});lu=d.get('LU',t.get('LU'));rhs=d.get('complex_RHS',t.get('complex_solve_RHS'));counter=dict(LU=lu,RHS=rhs,**{k:v for k,v in t.items() if k not in ['LU','RHS']})
        add('setup:'+rel,'setup_generation' if 'generation/' in rel else'setup_NN','Apple_local',f,cpu=d.get('process_CPU_seconds'),wall=d.get('process_wall_seconds'),counter=counter,dimension='576 and1296 scalar systems' if 'generation/' in rel else None,status=d.get('status'),physics_include='generation/' not in rel,scope=d.get('cost_scope','Outer measured interval; case/setup timers already included; final metadata/exit outside interval'))
        if 'generation/' in rel:
            casefile=f.parent/'CASE_GENERATION_COSTS.json'
            for case in load(casefile) if casefile.is_file() else[]:
                counts=case['counts'];lu=counts.get('LU');rhs=counts.get('actual_complex_RHS')
                if lu not in [1,2] or rhs is None:continue
                for dimension in ([576,1296] if lu==2 else[576]):add('generation_solve:'+case['case_id']+':'+str(dimension),'setup_generation','Apple_local',casefile,case['case_id'],counter=dict(LU=1,RHS=rhs/lu),dimension=dimension,cpu_include=False,scope='Single forward per matrix dimension reconstructed from counted generator calls; case CPU included in outer process')
    launches=root/'delegated/evidence/generation/LAUNCH_ATTEMPTS.json'
    if launches.is_file():
        for item in load(launches).get('attempts',[]):
            if item.get('status')!='COMPLETE':add('generation_failed_launch:'+str(item['attempt']),'setup_failed','Apple_local',launches,cpu=item.get('CPU_seconds'),wall=item.get('wall_seconds'),counter=dict(LU=item.get('physical_LU'),RHS=0 if item.get('physical_LU')==0 else None),status=item.get('status'),reason='Failed before generation; original launch timing unknown')
    add('inherited_NN_training','inherited_prior','UNMEASURED_LEGACY','Existing R8 supervised physics weights; no new training',cpu=None,counter=None,cpu_include=False,physics_include=False,reason='Historical training CPU/GPU/forward cost unknown; not zero')
    known_bundles=set()
    for f in (root/'audit').glob('*BUNDLE.json'):
        d=load(f);known_bundles.add(d.get('bundle'));add('runtime_pack:'+f.name,'setup_packing','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=0,RHS=0),status=d.get('status','PACKED'),scope='Actual runtime packaging timer; no physics')
    dependency=root/'audit/PHYSICS_DEPENDENCY_RECOVERY_PLAN.json'
    if dependency.is_file():known_bundles.add(load(dependency).get('bundle'))
    known_bundles={Path(str(x).replace('\\','/')).name for x in known_bundles if x}
    for f in root.glob('runtime_*.zip'):
        if f.name not in known_bundles:add('runtime_pack_untimed:'+f.name,'setup_packing','Apple_local',f,cpu=None,counter=dict(LU=0,RHS=0),reason='Archive exists but no matching outer packaging timer receipt found; file bytes do not imply CPU')
    prep=root/'audit/DEVELOPMENT_INPUT_PREPARATION.json'
    if prep.is_file():
        d=load(prep);add('development_input_preparation','setup_packing','Apple_local',prep,cpu=d.get('CPU_seconds'),counter=dict(LU=0,RHS=0),reason='Wall not captured in receipt')
    # Exposed profile/forecast/pressure/witness diagnostics: only outer totals count.
    full=root/'delegated/full_current'
    for pattern in ['profile_diagnostic_v2/*/PROFILE.json','forecast_diagnostic_v2/*/FORECAST.json','dual_pressure_diagnostic_v2/*/DUAL_PRESSURE.json','same_shadow_witness/WITNESS.json']:
        for f in full.glob(pattern):
            d=load(f);counter=d.get('physical_model_cost',dict(LU=d.get('new_Maxwell_LU'),RHS=0 if d.get('new_Maxwell_LU')==0 else None));add('diagnostic:'+str(f.relative_to(full)),'validation_exposed','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=counter,dimension=576,status=d.get('status'),scope=d.get('timer_scope','Outer diagnostic CPU; nested preparation/route/forecast timers are included'))
            for field in ['profile_cost','projector_cost','total_dense_operator_cost']:
                if d.get(field):add('diagnostic_operators:'+str(f.relative_to(full))+':'+field,'validation_operator_reference','Apple_local',f,field,counter=d[field],dimension=576,cpu_include=False,physics_include=False,reason='Nested operator/QR/SVD/RHS counts retained separately; not additional physical solves or CPU')
    # Wiring execution leaves: parent VALIDATION summaries repeat these receipts.
    for sub in ['validation','equal_area_wiring_v2','nfft_equations_wiring']:
        for f in (full/sub).rglob('*EXECUTION.json'):
            d=load(f);counter=d.get('whole_adapter_physics_delta',d.get('physical_model_delta'));add('wiring:'+str(f.relative_to(full)),'validation_wiring','Apple_local',f,cpu=d.get('whole_adapter_CPU_seconds',d.get('CPU_seconds',d.get('cpu_seconds'))),wall=d.get('whole_adapter_wall_seconds',d.get('wall_seconds')),counter=counter,dimension='576 scalar' if sub!='equal_area_wiring_v2' else'original equal-area wiring grid; receipt scope',status=d.get('status'),scope=d.get('timer_scope','Adapter outer time; nested setup/source counters not added again'))
            for field in ['operation_counts','denominator_setup_cost','adapter_diagnostic_cost']:
                if d.get(field):add('wiring_operators:'+str(f.relative_to(full))+':'+field,'validation_operator_reference','Apple_local',f,field,counter=d[field],cpu_include=False,physics_include=False,reason='Nested adapter counts; no duplicate timing/physical RHS')
    nfft_validation=full/'nfft_equations_wiring/VALIDATION.json'
    if nfft_validation.is_file():
        dd=load(nfft_validation)
        for row in rows:
            if row['cost_id']=='wiring:nfft_equations_wiring/EXECUTION.json':row['CPU_known_sum_include']=False;row['missing_or_overlap_reason']='Included in measured outer NFFT validation interval'
        add('nfft_outer_validation','validation_wiring','Apple_local',nfft_validation,'validation_CPU_seconds',cpu=dd.get('validation_CPU_seconds'),wall=dd.get('validation_wall_seconds'),physics_include=False,scope='Outer wiring/derivative validation; complete EXECUTION CPU nested here, interrupted execution remains separate')
        if dd.get('derivative_check_extra_cost'):add('nfft_derivative_operator_counts','validation_operator_reference','Apple_local',nfft_validation,'derivative_check_extra_cost',counter=dd['derivative_check_extra_cost'],cpu_include=False,physics_include=False,reason='Operator calls, not Maxwell solve RHS; included in outer validation CPU')
    # Analytic physics qualification outer CPU plus dimension-specific disjoint solves.
    phys=root/'delegated/physics'
    for name in ['QUALIFICATION.json','VOLUME_MIE_QUALIFICATION.json','MIE_REFINEMENT_DIAGNOSIS.json']:
        f=phys/name
        if not f.is_file():continue
        d=load(f);counts=d.get('counts');add('physical_qualification:'+name,'validation_physics','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=counts,dimension=648 if counts else None,physics_include=counts is not None,status=d.get('status'),scope='Outer qualification timer; nested analytic solve times included')
        for key in ['analytic_Mie_rows','rows','staircase_control']:
            for i,item in enumerate(d.get(key,[])):
                solve=item.get('physical_solve')
                if solve:add('analytic_solve:'+name+':'+key+str(i),'validation_physics','Apple_local',f,f'{key}[{i}].physical_solve',cpu=solve.get('cpu_s'),wall=solve.get('wall_s'),counter=solve,dimension=3*solve['N'] if solve.get('N') else None,cpu_include=False,scope='Physical solve count is disjoint; timer already included in outer qualification')
    first=phys/'MIE_REFINEMENT_VOLUME_ONLY_FIRST_ATTEMPT.json'
    if first.is_file():
        d=load(first);add('Mie_first_summary_overlap','validation_reference','Apple_local',first,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),cpu_include=False,physics_include=False,reason='Same volume-solve receipts retained unchanged in final combined refinement summary; overlap unresolved, not counted twice')
    f=phys/'new3d/GENERATION_RECEIPT.json'
    if f.is_file():
        d=load(f);add('new3d_generation_timer','setup_generation','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),physics_include=False,scope='Outer generation timer includes inverse/fine setup; nested setup timers not added again')
        for i,item in enumerate(d.get('rows',[])):add('new3d_generation_solve:'+str(i),'setup_generation','Apple_local',f,f'rows[{i}].physical_costs',counter=item.get('physical_costs'),dimension=3*int(item['generation_side'])**3 if item.get('generation_side') else None,cpu_include=False,reason='Per-solve CPU not separately timed; included in outer generation interval')
    f=phys/'EXTERNAL_INPUT_RECEIPT.json'
    if f.is_file():
        d=load(f);add('external_input_timer','setup_validation','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),physics_include=False,scope='Outer input/calibration/verification timer; nested projector/verification CPU not added again')
        for i,item in enumerate(d.get('Fresnel',[])):add('external_Fresnel_checks:'+str(i),'setup_validation','Apple_local',f,f'Fresnel[{i}].physical_costs',counter=item.get('physical_costs'),dimension=item.get('G_shape',[None])[0],cpu_include=False,scope='Disjoint actual measured-input physical checks; not target reconstruction')
    for rel in ['INTEGRATION_CHECK.json','FRESNEL_PERIODIC_AUDIT_CHECK.json','classic_public/PACKAGING_RECEIPT.json','new3d/runner_public/PACKAGING_RECEIPT.json','fresnel_public/PACKAGING_RECEIPT.json']:
        f=phys/rel
        if not f.is_file():continue
        d=load(f);counter=d.get('counts') or(dict(LU=d.get('new_LU'),RHS=d.get('new_RHS')) if 'new_LU' in d else None)
        add('physical_pack_check:'+rel,'setup_validation','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=counter,status=d.get('status'),reason='Outer CPU or full solve counters absent where null; not reconstructed from successful comparison values',scope='Receipt-level totals only; auxiliary/current counters remain distinct')
        auxiliary=d.get('auxiliary',{})
        if isinstance(auxiliary,dict) and auxiliary.get('CPU_seconds') is not None:add('physical_check_partial_timer:'+rel,'setup_validation','Apple_local',f,'auxiliary.CPU_seconds',cpu=auxiliary['CPU_seconds'],wall=auxiliary.get('wall_seconds'),cpu_include=d.get('CPU_seconds') is None,physics_include=False,reason='Known inner interval; outer process timer absent' if d.get('CPU_seconds') is None else'Already included in outer timer')
    # Compression/download/release receipts: names are evidence, no invented transfers.
    for f in (root/'archives').glob('*ARCHIVE_RECEIPT.json'):
        d=load(f);add('compression:'+f.name,'release_compression',d.get('host','Apple_local' if 'local' in f.name.lower() else 'i9_remote'),f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=0,RHS=0),status=d.get('status','ARCHIVED'),scope='Archive compression process timer, separate from solve CPU')
    # Download costs are read only from audit/ARTIFACT_TRANSFER_COSTS.json below.
    if not (root/'audit/ARTIFACT_TRANSFER_COSTS.json').is_file():add('artifact_download_timing_unknown','release_download','Apple_local','Sanitized ARTIFACT_TRANSFER_COSTS receipt absent',cpu=None,counter=dict(LU=0,RHS=0),reason='Public DOWNLOADS.json is an index, not a timing receipt; bytes do not determine time')
    # Public selector: batch scan is additional CPU, task intervals nested inside.
    for f in (root/'predictions/min_change_paid_pool').glob('*/SELECTIONS.json') if(root/'predictions/min_change_paid_pool').is_dir() else[]:
        d=load(f);name=d.get('campaign',f.parent.name);selector_host=d.get('host','Apple_local')
        add('min_change_selector:'+name,'public_selection',selector_host,f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=d.get('new_Maxwell_LU',0),RHS=0),status=d.get('status'),scope=d.get('cost_policy','Additional paid-pool scan CPU; all original proposal/corrector costs retained'))
        for i,task in enumerate(d.get('tasks',[])):
            add('min_change_selector_task:'+name+':'+task['task_id'],'public_selection_reference',selector_host,f,f'tasks[{i}]',cpu=task.get('CPU_seconds'),wall=task.get('wall_seconds'),counter=dict(LU=0,RHS=0),cpu_include=False,physics_include=False,reason='Per-task scan interval nested in the authoritative selector batch CPU; never add twice')
    for name,host in [('ENVIRONMENT_APPLE.json','Apple_local'),('ENVIRONMENT_I9.json','i9_remote'),('ENVIRONMENT_READBACK_FAILED_ATTEMPT.json','unrecorded_metadata_host')]:
        f=root/'audit'/name
        if not f.is_file():continue
        d=load(f);add('environment_readback:'+name,'environment_metadata',host,f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=d.get('new_physical_LU',0),RHS=0),status=d.get('status','METADATA_READBACK'),reason=d.get('unknown_reason'),scope='Metadata query only; verified AppleM4/10cores/24GiB and i9 package versions. Failed threadpoolctl query CPU unknown,0physical calls.')
    # Scoring/plotting are real work, but most historic runs have no timer.
    for f in (root/'analysis').glob('*/SCORING_RECEIPT.json'):
        d=load(f);add('scoring_plot:'+f.parent.name,'scoring_plotting','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=0,RHS=0),reason='Scorer/figure CPU and wall were not recorded where null; not zero',scope='Offline scoring/plot process only; source online costs excluded')
    for f in (root/'analysis').glob('*/PLOT_RECEIPT.json'):
        d=load(f);add('measured_plot:'+f.parent.name,'scoring_plotting','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=d.get('new_Maxwell_LU',0),RHS=0),scope='Measured plot-only interval')
    for f in (root/'analysis').rglob('AGGREGATION_RECEIPT.json'):
        if any('fixture' in part.lower() or 'fake' in part.lower() for part in f.parts):continue
        d=load(f);add('table_aggregation:'+str(f.relative_to(root)),'scoring_plotting','Apple_local',f,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),counter=dict(LU=0,RHS=0),reason='Table/bootstrap/plot aggregation timer not measured; original online/source costs are not new work')
    supplemental_receipts(root,add)
    site=root/'SITE_BUILD.json'
    if site.is_file():
        d=load(site);add('static_site_build','release_site','Apple_local',site,cpu=None,wall=d.get('wall_seconds'),counter=dict(LU=0,RHS=0),reason='Builder recorded wall only; historic builder reruns were not separately timed')
    summary={}
    for row in rows:
        key=row['host']+' / '+row['phase'];group=summary.setdefault(key,dict(host=row['host'],phase=row['phase'],known_CPU_seconds=0.,unknown_CPU_records=0,records=0))
        group['records']+=1
        if row['CPU_known_sum_include']:
            if row['CPU_seconds'] is None:group['unknown_CPU_records']+=1
            else:group['known_CPU_seconds']+=row['CPU_seconds']
    physics={}
    for row in rows:
        if not row['physical_known_sum_include']:continue
        key=row['host']+' / '+str(row['system_dimension']);group=physics.setdefault(key,dict(host=row['host'],system_dimension=row['system_dimension'],known_LU=0,known_complex_RHS=0,unknown_physics_records=0))
        if row['physical_LU'] is not None:group['known_LU']+=row['physical_LU']
        if row['physical_complex_RHS'] is not None:group['known_complex_RHS']+=row['physical_complex_RHS']
        if row['physical_LU'] is None or row['physical_complex_RHS'] is None:group['unknown_physics_records']+=1
    record=dict(scope='Known lower-bound metadata inventory, not full final study bill or wall duration',created_at_utc=datetime.now(timezone.utc).isoformat(),entries=rows,CPU_by_host_phase=list(summary.values()),physical_counts_by_host_dimension=list(physics.values()),notes=notes,inventory_CPU_seconds=time.process_time()-started,inventory_wall_seconds=time.perf_counter()-wall0,truth_or_material_arrays_read=False,new_numerical_work=0)
    (out/'COST_INVENTORY.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    with (out/'COST_INVENTORY.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in rows)
    (out/'COST_SUMMARY.json').write_text(json.dumps(dict(CPU_by_host_phase=list(summary.values()),physical_counts_by_host_dimension=list(physics.values()),completion='PARTIAL_METADATA_SNAPSHOT; unknowns retained'),indent=2)+'\n')
    md=['# R13 已知费用清单（元数据快照）','','这是按host/阶段分列的已知消耗下界，不是最终全费用或并行历时。当前无终态回执的队列保留为不完整；缺失/中断/评分计时为空，不补零。','','|Host / 阶段|已知CPU秒|未知CPU记录|','|---|---:|---:|']
    md += [f"|{r['host']} / {r['phase']}|{r['known_CPU_seconds']:.6f}|{r['unknown_CPU_records']}|" for r in summary.values()]
    md += ['','外部child CPU只从权威队列计一次；内部计时、分析COST_TOTALS与恢复子队列为引用。队列wall不与child wall求和。不同矩阵维数的LU/RHS分列，projector/current/Gram动作不能当作物理LU RHS。','','## 明确未知或重叠','']
    md += [f"- `{r['cost_id']}`：{r['missing_or_overlap_reason']}" for r in rows if r['missing_or_overlap_reason'] and (r['CPU_seconds'] is None or not r['CPU_known_sum_include'])]
    (out/'COST_SUMMARY.md').write_text('\n'.join(md)+'\n');print(json.dumps(dict(entries=len(rows),output=str(out),unknown_CPU_records=sum(g['unknown_CPU_records'] for g in summary.values()),new_numerical_work=0)))
if __name__=='__main__':main()
