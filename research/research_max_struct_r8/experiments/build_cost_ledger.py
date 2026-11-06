"""Read-only receipt accounting. Never imports a solver, reads arrays/GT, or scores.

Run from any directory with Python standard library only. Writes under
 delegated/cost_package by default. Registry patterns below are intentional;
 HISTORY, PROGRESS, AUXILIARY_COST and queue child totals are NOT extra runs.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import argparse
import json
import os
import csv

ROOT = Path(__file__).resolve().parents[1]


def number(value):
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def read(path):
    if 'sealed_truth' in path.parts:
        raise ValueError('Private truth read prohibited')
    return json.loads(path.read_text())


def build(root):
    rows, contexts, failures, gaps = [], [], [], []
    used = set()
    remote = set(read(root/'audit/HARDWARE_ASSIGNMENT.json').get('remote_objects', []))
    relocation=root/'audit/HARDWARE_REALLOCATION.json'
    if relocation.exists():remote-=set(read(relocation).get('objects',[]))
    hardware_expected={}
    postcheck=root/'audit/HARDWARE_REALLOCATION_POSTCHECK.json'
    if postcheck.exists():
        pc=read(postcheck);hardware_expected=dict(primary=pc.get('primary_actual_solver_tasks'),anchor=pc.get('anchor_actual_solver_tasks'),note='Final dispatch allocation, not measured completed physics. Relocated queue acknowledgments never count as solves.')

    def add(path, category, family, data=None, cpu=None, wall=None, note='', additive=True, suffix=''):
        path = Path(path)
        key = (str(path.relative_to(root)), suffix)
        if key in used:
            raise ValueError('Duplicate accounting row: '+str(key))
        used.add(key)
        d = read(path) if data is None else data
        cpu = number(d.get('cpu_seconds', d.get('CPU_seconds', d.get('cpu_s')))) if cpu is None else number(cpu)
        wall = number(d.get('wall_seconds', d.get('wall_s'))) if wall is None else number(wall)
        case = d.get('case_id', '')
        host = d.get('execution_platform',{}).get('host') or (('XINAN' if case in remote else 'Mac') if (family.startswith('warm') or family=='history_warm') else d.get('host','SOURCE_RECEIPT_SCOPE'))
        host_raw=host
        platform=d.get('execution_platform',{}).get('system','').lower();exe=str(d.get('executable',''))
        if platform=='windows' or exe.startswith('D:') or 'xinan' in str(host).lower():host='XINAN'
        elif platform=='darwin' or exe.startswith('/Volumes/') or 'mac-mini' in str(host).lower() or host=='Mac':host='Mac'
        model = d.get('physical_model_cost', d.get('physical_cost', d.get('counts', d.get('operation_counts', {k:d[k] for k in ['LU','complex_RHS','physical_LU','physical_complex_RHS'] if k in d}))))
        if not isinstance(model,dict): model={}
        row = dict(source=key[0], record=suffix or 'outer_receipt', category=category,
                   family=family, status=d.get('status','RECORDED'), hardware=host,hardware_original_label=host_raw,
                   cpu_seconds=cpu, task_wall_seconds=wall, additive=additive,
                   task_id=d.get('task_id'), matrix_dimension=model.get('full_factorization_dimension'),
                   operation_counts=model, monetary_cost=None, energy=None, note=note,
                   nonadditive_phase_breakdown={k:d[k] for k in ['setup_wall_seconds','setup_cpu_seconds',
                       'setup_shared','auxiliary_wall_seconds','auxiliary_cpu_seconds',
                       'auxiliary_execution_receipt','corrector','physical_cost_training',
                       'physical_cost_prediction','physical_cost_optimizer','physical_cost_validation','optimization_cpu_seconds','model_load_cpu_seconds','model_load_wall_seconds'] if k in d})
        rows.append(row)
        if additive and cpu is None:
            gaps.append(dict(source=key[0], gap='CPU not metered/recorded; not zero'))
        return row

    def matches(pattern):
        return sorted(root.glob(pattern))

    # Full generator receipt includes feature/preprocessing and per-object work.
    # Never sum its per_object children in addition to the outer receipt.
    for p in matches('audit/generation*.json'):
        r=read(p);add(p,'data_generation','warm_data',r,
                     note='Outer generator: includes features; per_object rows overlap and are not added.')
    p=root/'audit/validation_downstream_generation.json'
    if p.exists():
        add(p,'data_generation','warm_validation_observations',
            note='Separate prepare_warm_validation invocation regenerates eight public observations from existing validation materials; outer receipt includes feature setup. New materials=0.')
    for p in matches('RESULTS/training/*/RECEIPT.json'):
        add(p,'training','neural_training',note='Includes in-epoch validation/model selection; training history cumulative timers not added.')
    for p in matches('data_public/initializers/*/RECEIPT.json'):
        add(p,'inference','frozen_initializer',note='One actual frozen initializer invocation; per-case inference/model load nested. Public feature cost belongs to generator.')
    for p in matches('audit/*FEATURE_INFERENCE*.json'):
        r=read(p);add(p,'validation','deployment_measurement',r,
            note='Separate measured deployment/counterfactual benchmark; not substituted for actual original initializer cost.')

    # Frozen task union deduplicates overlapping local/remote/all manifests.
    taskmap={}
    for p in matches('configs/warm_*manifest.json'):
        r=read(p)
        for t in (r if isinstance(r,list) else r.get('tasks',[])):
            ident=t['task_id'];entry=taskmap.setdefault(ident,dict(task=t,manifests=[]))
            entry['manifests'].append(str(p.relative_to(root)))
    expected={k:v for k,v in taskmap.items() if '/warm_primary/' in v['task'].get('output_dir','')}
    completed=[];missing=[];running=[]
    for family in ['warm_primary','warm_development','warm_anchor_development','warm_anchor_development_extended','warm_anchor_confirmation','history_warm']:
        for p in matches('RESULTS/'+family+'/*/RECEIPT.json'):
            r=read(p);st=r.get('status','')
            if st=='SEALED_PUBLIC_PREDICTION':
                online=add(p,'online_inversion',family,r,note='Outer process minus separately metered setup only. Remaining imports/load/auxiliary/corrector/validation/serialization included; nested timers NOT added.')
                sc=number(r.get('setup_cpu_seconds'));sw=number(r.get('setup_wall_seconds'))
                if sc is not None and sw is not None and online['cpu_seconds'] is not None and online['task_wall_seconds'] is not None:
                    if sc>online['cpu_seconds'] or sw>online['task_wall_seconds']:raise ValueError('Invalid nested setup scope '+str(p))
                    online['outer_receipt_cpu_seconds']=online['cpu_seconds'];online['outer_receipt_wall_seconds']=online['task_wall_seconds']
                    online['cpu_seconds']-=sc;online['task_wall_seconds']-=sw
                    setup=add(p,'setup',family,r,cpu=sc,wall=sw,suffix='embedded_setup',note='Nested setup partition subtracted from online row; not additional process.')
                    setup['operation_counts']={};setup['nonadditive_phase_breakdown']={}
                if family=='warm_primary':completed.append(r['task_id'])
            elif st.startswith('FAIL'):
                add(p,'failed',family,r);failures.append(dict(source=str(p.relative_to(root)),status=st))
            else:
                running.append(r.get('task_id',str(p.parent.name)))
                contexts.append(dict(source=str(p.relative_to(root)),status=st,note='Nonterminal/control acknowledgment has no final additive solver cost; RELOCATED_BEFORE_EXECUTION is dispatch, never physics.'))
    for task in sorted(expected):
        if task not in completed:missing.append(task)
    for p in matches('audit/queue*/RECEIPT.json'):
        r=read(p);contexts.append(dict(source=str(p.relative_to(root)),status=r.get('status'),
            started_utc=r.get('started_utc'),batch_wall_seconds=r.get('batch_wall_seconds'),
            elapsed_wall_seconds=r.get('elapsed_wall_seconds'),
            note='Queue elapsed is context only; child process_wall overlaps run_case and is not added.'))

    # Real failed subprocess attempts are separate outer invocations even when
    # later repaired tasks reuse the same task_id in final RESULTS. Queue wrapper
    # and orchestration explanations do not add these children again.
    for p in matches('audit/failed_missing_threadpoolctl/*/RECEIPT.json'):
        if p.parent.name=='queue':
            contexts.append(dict(source=str(p.relative_to(root)),receipt=read(p),note='Failed batch wrapper overlaps child failed invocations; nonadditive.'));continue
        add(p,'failed','warm_remote_dependency_failure',note='Actual failed run_case process before physics dependency import; separate attempt source path, later successful same task_id not duplicate compute.')
        failures.append(dict(source=str(p.relative_to(root)),status=read(p).get('status')))

    # Incremental fetch history is authoritative for each independently paid
    # remote compression invocation. Latest FETCH.json mirrors are context only.
    fetch_sources=matches('audit/*FETCH_HISTORY.jsonl')
    for p in fetch_sources:
        for i,line in enumerate(p.read_text().splitlines()):
            if not line.strip():continue
            data=json.loads(line);suffix='invocation/'+str(i+1)
            remote_wall=number(data.get('remote_pack_wall_seconds'));control=number(data.get('control_wall_seconds'))
            row=add(p,'publication_preparation','remote_incremental_pack',data,cpu=data.get('remote_pack_cpu_seconds'),wall=remote_wall,
                suffix=suffix+'/remote_pack',note='One actual remote archive creation per fetch history entry. Repeated data still consumes paid compression work; final release packaging not run.')
            row['hardware']='XINAN';row['operation_counts']={}
            remainder=control-remote_wall if control is not None and remote_wall is not None and control>=remote_wall else None
            controlrow=add(p,'publication_preparation','fetch_control_transport',{},cpu=None,wall=remainder,suffix=suffix+'/control_remainder',
                note='Local/SSH/SCP/extract control wall minus nested remote-pack wall only; local CPU NOT_METERED. Queue elapsed and latest FETCH mirror not added.')
            controlrow['hardware']='Mac';controlrow['control_total_wall_seconds']=control
    for p in matches('audit/*FETCH.json'):
        contexts.append(dict(source=str(p.relative_to(root)),receipt=read(p),note='Latest fetch mirror/queue progress; authoritative independent costs only in FETCH_HISTORY lines. Not another invocation.'))
    for p in [root/'audit/HARDWARE_REALLOCATION.json',root/'audit/HARDWARE_REALLOCATION_POSTCHECK.json',root/'audit/XINAN_LAUNCH.json']:
        if not p.exists():continue
        data=read(p)
        if p.name=='HARDWARE_REALLOCATION.json':
            r=add(p,'setup','hardware_relocation_control',data,wall=data.get('control_wall_seconds'),note='Control-only acknowledgment, CPU not metered, NO physics. Actual counterpart solves charged only in their Mac RESULTS receipt.');r['operation_counts']={}
        else:
            contexts.append(dict(source=str(p.relative_to(root)),receipt=data,note='Final assignment/dispatch metadata, not physical solve count or separately metered CPU.'))
            if p.name=='HARDWARE_REALLOCATION_POSTCHECK.json':gaps.append(dict(source=str(p.relative_to(root)),gap='Postcheck control CPU unmetered; expected final Mac405/66 and XINAN345/60 are allocation targets, not actual completion totals.'))

    # Isolated remote anchor checkout/control metadata is not solver compute.
    prepared=root/'audit/ANCHOR_XINAN_PREPARED.json'
    if prepared.exists():
        data=read(prepared)
        add(prepared,'setup','anchor_remote_transfer',data,wall=data.get('transfer_wall_seconds'),
            note='Independently metered transfer wall only; CPU and remaining isolated-checkout preparation unmetered. Not queue child compute.')
    for p in sorted(set(matches('audit/*ANCHOR*XINAN*.json')+matches('audit/*anchor*xinan*.json'))):
        if p==prepared:continue
        contexts.append(dict(source=str(p.relative_to(root)),receipt=read(p),note='Isolated remote checkout/dispatch/control receipt; not independent solver work. Waiting/import/setup scope not inferred.'))

    # Existing paper accountant explicitly subtracts reused auxiliary costs and
    # retains partial failed-process bounds. Import its rows ONCE, not plus raws.
    paper=root/'RESULTS/paper_benchmarks/COST_LEDGER.json'
    if paper.exists():
        r=read(paper)
        for i,d in enumerate(r['rows']):
            kind=d.get('kind','');cat=('failed' if any(x in kind for x in ['failed','interrupted']) else
                'data_generation' if 'forward' in kind else 'validation' if 'offline' in kind else
                'inference' if 'inference' in kind else 'online_inversion')
            row=add(paper,cat,'paper_benchmarks',d,cpu=d.get('CPU_seconds'),wall=d.get('wall_seconds'),
                    suffix='rows/'+str(i),note='Import existing non-reuse executed-cost row; raw source '+str(d.get('path'))+'; '+str(d.get('CPU_status','')))
            row['original_row']=d
            if 'lowerbound' in str(d.get('CPU_status','')).lower():row['cost_qualification']='LOWER_BOUND'
        gaps.append(dict(source=str(paper.relative_to(root)),gap=r.get('warning','Original paper accounting exclusions remain.')))
        shared_seen=set()
        for p in matches('RESULTS/paper_benchmarks/*/*/RECEIPT.json'):
            d=read(p);setup=d.get('setup_shared')
            if not setup:continue
            key=(str(p.parent.parent),json.dumps(setup,sort_keys=True))
            if key in shared_seen:continue
            shared_seen.add(key)
            add(p,'setup','paper_shared_32grid',setup,suffix='setup_shared_once',
                note='One case assembly before per-arm timer in paper_benchmark.py; duplicate timing copied into other arm receipts is not added.')
    for rel,family,category in [
        ('RESULTS/warm_development_analysis/RECEIPT.json','warm_validation_analysis','validation'),
        ('RESULTS/warm_anchor_development_analysis/RECEIPT.json','warm_anchor_development_analysis','validation'),
        ('data_public/history_warm/PREPROCESS_RECEIPT.json','history_warm_preprocess','inference'),
        ('RESULTS/history_warm/SCORING_RECEIPT.json','history_warm_scoring','scoring'),
        ('delegated/warm_anchor_analysis/PILOT_RECEIPT.json','warm_anchor_pilot_plotting','validation'),
        ('delegated/warm_anchor_analysis/DEVELOPMENT_RECEIPT.json','warm_anchor_development_analysis_extended','validation'),
        ('delegated/warm_anchor_analysis/CONFIRMATION_RECEIPT.json','warm_anchor_confirmation_analysis','scoring'),
        ('delegated/warm_mechanism_analysis/RECEIPT.json','warm_initializer_survival_analysis','scoring'),
        ('delegated/warm_analysis/CONFIRMATION_EXECUTION.json','warm_confirmation_tables_figures','scoring'),
        ('RESULTS/neural_survival_visualization/RECEIPT.json','neural_initializer_endpoint_visualization','scoring'),
        ('delegated/warm_direct_budget_analysis/COST.json','direct_physics_budget_extension','scoring'),
        ('audit/DOCUMENT_LINK_REVIEW.json','final_document_link_QA','validation')]:
        p=root/rel
        if p.exists():add(p,category,family,note='Independent explicitly registered outer preprocessing/analysis/scoring invocation; nested records and reconstructed solver ledgers not added.')
        else:gaps.append(dict(source=rel,gap='Expected outer receipt not available yet; no cost imputed.'))
    for p in matches('RESULTS/projection_transport_*/RECEIPT.json'):
        add(p,'validation','public_projection_transport',note='Paid initial states reused; public-only state-fixed diagnostic, no new full-wave solve or GT.')
    for rel in ['RESULTS/history_warm/ALL_PREDICTIONS_SEALED.json','delegated/fresnel/MULTIFREQUENCY_BATCH_RECEIPT.json','delegated/fresnel/FINALIZATION_RECEIPT.json']:
        p=root/rel
        if p.exists():contexts.append(dict(source=rel,receipt=read(p),note='Parent batch/process/report wrapper; may include children, retries, predecessor wait and reuse. Not additive with registered outer children.'))

    # Exactly the published successful-arm table, not process COMPLETE totals.
    # The first PR process retained five arms before failed serialization; its
    # resume COMPLETE process excludes them. No subtraction from unequal scopes.
    vp=root/'RESULTS/vector3d/ONLINE_COST_LEDGER.json'
    if vp.exists():
        r=read(vp)
        for i,d in enumerate(r['successful_arms']):
            row=add(vp,'online_inversion','vector3d',d,suffix='successful_arms/'+str(i),
                note='One successful arm; dimensions not interchangeable with scalar TM LU; process totals excluded.')
            row['original_row']=d
        gaps.append(dict(source=str(vp.relative_to(root)),gap=r.get('scope')))
    for p in matches('RESULTS/vector3d/**/SIMULATION.json')+matches('RESULTS/vector3d/**/COMPLETE.json')+matches('RESULTS/vector3d/**/PR_COMPLETE.json'):
        r=read(p);contexts.append(dict(source=str(p.relative_to(root)),receipt=r,
            note='Shared setup/simulation/process context; reused sealed simulation is not new compute; do not add process total to contained arms.'))
    p=root/'RESULTS/vector3d/original_size/COSTS.json'
    if p.exists():add(p,'online_inversion','vector3d_original_size',note='Whole truncated attempt; includes setup/forward/online paid work and cannot be split exactly from this receipt.')
    p=root/'RESULTS/vector3d/refinement8_uncached_attempt/ENGINEERING_INTERRUPTION.json'
    if p.exists():
        r=read(p);token=r['ps_process_measurement'].split()[1]
        parts=token.split(':');cpu=sum(float(v)*60**i for i,v in enumerate(reversed(parts)))
        row=add(p,'failed','vector3d_engineering',r,cpu=cpu,
            note='OS process CPU LOWER BOUND includes completed uncached LSIE 227.39 s; its arm receipt must not be added again. Wall bound not parsed from ambiguous OS columns.')
        row['cost_qualification']='LOWER_BOUND'
        failures.append(dict(source=str(p.relative_to(root)),receipt=r))
    p=root/'RESULTS/vector3d/RETAINED_FAILURE.json'
    if p.exists():add(p,'failed','vector3d_serialization',note='Failed first A+ CPU NOT_MEASURED; initial array overwritten, snapshot NOT_RETAINED.');failures.append(dict(source=str(p.relative_to(root)),receipt=read(p)))

    # Fresnel outer inversion CPU includes auxiliary/current/corrector and final
    # withheld prediction. Never add AUXILIARY_RECEIPT or PROGRESS counters.
    for p in matches('RESULTS/fresnel/*/RECEIPT.json'):
        r=read(p)
        family='fresnel_multifrequency' if p.parent.name.startswith('FoamTwinDielTM_multifrequency_') else 'fresnel_air0_development' if 'air0' in p.parent.name.lower() else 'fresnel'
        row=add(p,'online_inversion',family,r,note='Whole recorded inversion; setup/current/corrector/physical validation nested. Multi-frequency keeps optimizer and validation counts per frequency; anchor reused, no new charge for old run. No pixel-GT scoring.')
        if family in ['fresnel_multifrequency','fresnel_air0_development']:
            row['operation_counts_by_frequency']=dict(optimizer=r.get('physical_cost_optimizer',{}),validation=r.get('physical_cost_validation',{}))
            row['matrix_dimension']=3600
    for p in matches('RESULTS/fresnel/*/FAILED.json'):
        if (p.parent/'RECEIPT.json').exists():
            contexts.append(dict(source=str(p.relative_to(root)),receipt=read(p),note='Same-directory outer final receipt exists; manual review overlap, not extra additive attempt.'));continue
        family='fresnel_multifrequency' if p.parent.name.startswith('FoamTwinDielTM_multifrequency_') else 'fresnel_engineering'
        add(p,'failed',family,note='Actual failed outer invocation; nested paid incumbent/history not extra compute.')
        failures.append(dict(source=str(p.relative_to(root)),receipt=read(p)))
    for p in matches('RESULTS/fresnel/*/INTERRUPTED*.json'):
        r=read(p)
        if (p.parent/'FAILED.json').exists():
            contexts.append(dict(source=str(p.relative_to(root)),receipt=r,note='Interrupted OS bound overlaps same-directory metered FAILED outer receipt; not added.'));continue
        cpu=r.get('observed_CPU_seconds_lower_bound',r.get('CPU_seconds_lower_bound',r.get('cpu_seconds')))
        wall=r.get('observed_wall_seconds_lower_bound',r.get('observed_elapsed_seconds_lower_bound',r.get('elapsed_wall_seconds_lower_bound')))
        row=add(p,'failed','fresnel_engineering',r,cpu=cpu,wall=wall,
                note='Interrupted process bound; auxiliary/progress included and not separately added.')
        row['cost_qualification']='LOWER_BOUND_OR_UNKNOWN'
        failures.append(dict(source=str(p.relative_to(root)),receipt=r))
    # Known post-solver artifact QA has an independent timer; wrappers with
    # predecessor wait remain contextual and are never solver CPU.
    for rel in ['delegated/fresnel/AIR0_FINALIZATION_RECEIPT.json','delegated/fresnel/AIR0_CORRECTED_FINALIZATION_RECEIPT.json']:
        qa=root/rel
        if qa.exists():
            d=read(qa);add(qa,'validation','fresnel_artifact_qa',d,cpu=d.get('finalization_CPU_seconds'),note='Standalone artifact QA/report/plot CPU; completed physical runs are references, not repeated solver work.')
    packaging=root/'delegated/cost_package/PACKAGING_EXECUTION_RECEIPT.json'
    if packaging.exists():
        gaps.append(dict(source=str(packaging.relative_to(root)),gap='Packaging timer excludes imports/public plan+gates/preflight/post-receipt index refresh; later independent52-archive verification CPU not separately metered. Recorded compact+ZIP CPU is a paid-phase measurement, not complete orchestration CPU.'))
        add(packaging,'publication_preparation','final_release_packaging',note='Measured compact serialization+ZIP construction/member checking invocation, only after root executes. Timer excludes initial imports/gates/plan/preflight and post-receipt index refresh. Contains no new physics; member counters are not extra CPU.')
    else:
        contexts.append(dict(source='delegated/cost_package/PACKAGING_EXECUTION_RECEIPT.json',status='NOT_RUN',note='Actual release packaging not executed yet. Refresh final accountant after execution; do not invent packaging cost.'))
    gaps.append(dict(source='delegated/warm_direct_budget_analysis/EXECUTION_NOTES.md',gap='First syntax/startup failure CPU NOT_METERED; no input data or physics accessed. Successful COST.json is separate; later read-only verification CPU not included/not metered.'))
    contexts.append(dict(source='audit/remote_imports/**',note='Late fetched remote audit copies are provenance mirrors, not new actual execution. Never recursively aggregate their duplicate receipts/FETCH histories.'))
    for pat in ['RESULTS/fresnel/*VALIDATION*.json','RESULTS/fresnel/*DIAGNOSTIC*.json','RESULTS/fresnel/*AUDIT*.json','RESULTS/fresnel/*TRANSFER*.json']:
        for p in matches(pat):contexts.append(dict(source=str(p.relative_to(root)),receipt=read(p),note='Qualification context; additive timer scope needs review.'))
    for pat in ['audit/*FAILURE*.json','audit/failed*.json','audit/*DEPENDENCY_REPAIR*.json']:
        for p in matches(pat):
            failures.append(dict(source=str(p.relative_to(root)),receipt=read(p)))
            gaps.append(dict(source=str(p.relative_to(root)),gap='Orchestration/setup/repair overhead is not independently metered; no invented zero.'))

    # Score receipts may only be consumed after the root-level prediction seal.
    # Merely rerunning this accountant never causes a scorer or truth read.
    seal=root/'audit/WARM_ALL_PUBLIC_SEAL.json'
    scoring=[]
    if seal.exists():
        for pat in ['audit/*SCORING*RECEIPT*.json','RESULTS/warm_analysis/RECEIPT.json','RESULTS/warm_scoring/RECEIPT.json']:
            for p in matches(pat):add(p,'scoring','warm_offline_scoring');scoring.append(str(p.relative_to(root)))
    gaps.extend([
        dict(source='warm_primary',gap=f'{len(completed)}/{len(expected)} manifest tasks sealed locally; remote pending/untransferred indistinguishable without local final receipts.'),
        dict(source='global',gap='No exact full experiment grand total: shared setup, report/plot/scoring/publication and failed first serialization are partly unmetered.'),
        dict(source='global',gap='Money/energy not recorded. No currency or zero-cost claim.'),
        dict(source='setup',gap='Standalone public/operator setup not universally metered; warm embedded setup partition is subtracted from online; paper_benchmark shared setup added once per case/timer; vector shared process scopes remain contextual.'),
        dict(source='scoring',gap='Primary scoring NOT_RUN/NOT_IMPORTED unless prediction seal and explicit scoring receipt exist; never read sealed_truth.')])
    cohort_status={}
    for rel,family in [('configs/warm_anchor_development_manifest.json','warm_anchor_development'),('configs/warm_anchor_development_extended_manifest.json','warm_anchor_development_extended'),('configs/warm_anchor_confirmation_manifest.json','warm_anchor_confirmation'),('experiments/history_warm_manifest.json','history_warm')]:
        p=root/rel
        if p.exists():
            data=read(p);tasks=data.get('tasks',[]);sealed=[r['task_id'] for r in rows if r['family']==family and r['category']=='online_inversion']
            cohort_status[family]=dict(manifest=rel,expected=len(tasks),sealed_receipts=len(sealed),missing=[t['task_id'] for t in tasks if t['task_id'] not in sealed],scope=('SEPARATE_FROZEN126_CONFIRMATION; not original750; all complete raw states retained' if family=='warm_anchor_confirmation' else 'DEVELOPMENT_ONLY; not part of750primary; old baseline costs not re-added'))
    sums=defaultdict(lambda:dict(rows=0,known_cpu_seconds=0.0,cpu_unknown_rows=0,sum_task_wall_seconds=0.0,task_wall_unknown_rows=0))
    for r in rows:
        if not r['additive']:continue
        s=sums[(r['category'],r['family'],r['hardware'])];s['rows']+=1
        if r['cpu_seconds'] is None:s['cpu_unknown_rows']+=1
        else:s['known_cpu_seconds']+=r['cpu_seconds']
        if r['task_wall_seconds'] is None:s['task_wall_unknown_rows']+=1
        else:s['sum_task_wall_seconds']+=r['task_wall_seconds']
    if not missing:
        gaps=[g for g in gaps if g['source']!='warm_primary']
        contexts.append(dict(source='warm_primary',status='ALL750_RECEIPTS_AVAILABLE',note='Final original cohort receipt completeness; allocation and physics costs stay separate.'))
    if scoring:
        gaps=[g for g in gaps if g['source']!='scoring']
        contexts.append(dict(source='scoring',status='PRIMARY_SCORING_RECEIPT_IMPORTED',receipts=scoring,note='Only independent scoring outer CPU added; scores and reused solver totals are not new compute.'))
    return dict(schema='R8_NONOVERLAPPING_RECEIPT_LEDGER_V1',status='IN_PROGRESS',
        generated_utc=datetime.now(timezone.utc).isoformat(),truth_content_read=False,solver_invocations=0,hash_checks='NONE',
        accounting_rule='Only explicit outer invocation/retained executed-cost rows additive. Nested timers/queue/process totals contextual. Task-wall sums are workload, NEVER project elapsed. LU/RHS remain per physical model/dimension; no scalar/vector equivalence.',
        hardware_assignment_expected=hardware_expected,frozen_primary_expected=len(expected),additional_cohorts=cohort_status,primary_completed_receipts=len(completed),primary_missing_receipts=missing,
        running_receipts=running,scoring_receipts=scoring,rows=rows,
        grouped_known_costs=[dict(category=k[0],family=k[1],hardware=k[2],**v) for k,v in sorted(sums.items())],
        nonadditive_contexts=contexts,failures=failures,gaps=gaps,
        scientific_completion='Parent must judge only after all 750 primary tasks, Fresnel jobs, seals and scoring; receipt availability alone is not scientific acceptance.')


def packaging_inventory(root):
    """stat only; never opens payload, GT, original chats or third-party data."""
    files=[]
    for base,dirs,names in os.walk(root,followlinks=False):
        dirs[:]=sorted(d for d in dirs if d not in ['.git'])
        for name in sorted(names):
            p=Path(base)/name;rel=p.relative_to(root);s=str(rel);parts=[v.lower() for v in rel.parts];lo=s.lower()
            if parts[:2]==['delegated','cost_package']:continue
            size=p.lstat().st_size
            category='git';reason='Own code/config/report/compact result or metadata'
            if p.is_symlink():category='exclude';reason='Symlink target not followed or packaged'
            elif any(v in parts for v in ['__pycache__','.venv','site-packages','remote_vendor','node_modules']) or p.suffix in ['.pyc','.pyo']:
                category='exclude';reason='Local vendor/dependency/cache'
            elif any(v in parts for v in ['data_local','inputs']) or '.local.' in lo or any(v in lo for v in ['original_chat','source_chat','chat_history','chat_export']):
                category='exclude';reason='Original chat, third-party raw/local fulltext/download or local input boundary'
            elif (p.suffix.lower() in ['.pdf'] and 'FIGURES' not in rel.parts) or p.suffix.lower() in ['.epub','.docx','.eml','.pem','.key','.pub','.download'] or (p.suffix.lower()=='.html' and 'delegated' in parts) or name.lower().startswith(('id_rsa','id_ed25519','known_hosts','authorized_keys')):
                category='exclude';reason='Third-party fulltext/download or SSH identity; manual review required for other ownership'
            elif (parts[0]=='audit' and p.suffix.lower() in ['.zip','.gz','.tgz','.tar']) or any(v in lo for v in ['ssh_identity','ssh_config']):
                category='exclude';reason='Temporary transfer archive or SSH connection artifact'
            elif 'sealed_truth' in parts or any('private' in part for part in parts):
                category='hold_private_truth';reason='Metadata-only inventory; release only after root confirms all predictions sealed and offline scoring/publication authorized'
            elif 'data_training' in parts or ('training' in parts and p.suffix.lower() in ['.pt','.pth','.npz','.npy']):
                category='release_training';reason='Own training data/weights; include full training provenance and exact epoch receipts in Git'
            elif any(v in lo for v in ['history','budget','paid_state','paid-state','trajectory','snapshots']) or name=='CORRECTOR.json' or (p.suffix.lower()=='.npz' and name.lower().startswith(('stage_','corrector_','auxiliary_candidates'))) or size>=10*1024*1024:
                category='release_evidence';reason='Large own histories/budget/state/evidence; Git retains compact final predictions/configs/receipts'
            if category=='git' and size>90*1024*1024:
                category='release_evidence';reason='Git blob size boundary'
            files.append(dict(path=s,bytes=size,extension=p.suffix.lower() or '(none)',category=category,reason=reason))
    grouped=defaultdict(lambda:dict(files=0,bytes=0))
    for f in files:grouped[f['category']]['files']+=1;grouped[f['category']]['bytes']+=f['bytes']
    packages=[]
    # Deterministic finite packs <=600 MiB by uncompressed size, no actual archive.
    cap=600*1024*1024
    for cat in ['release_evidence','release_training']:
        current=[];nbytes=0;idx=1
        for f in sorted((f for f in files if f['category']==cat),key=lambda x:x['path']):
            if current and nbytes+f['bytes']>cap:
                packages.append(dict(name=f'r8-{cat}-part{idx:03d}.zip',uncompressed_bytes=nbytes,members=current));idx+=1;current=[];nbytes=0
            current.append(f['path']);nbytes+=f['bytes']
        if current:packages.append(dict(name=f'r8-{cat}-part{idx:03d}.zip',uncompressed_bytes=nbytes,members=current))
    return dict(status='PROPOSAL_ONLY_NOT_PACKAGED_NOT_UPLOADED',generated_utc=datetime.now(timezone.utc).isoformat(),
        scan='Path/type/stat only; private/fulltext/raw-data contents never opened; cost_package output omitted to avoid recursion.',
        preserves='R1-R7 untouched. No gitignore changes. No hash checks.',groups=dict(grouped),files=files,proposed_packages=packages)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--output',type=Path);ap.add_argument('--finalized',action='store_true',help='Root confirms all authorized work complete; requires both public seal/scoring gates; scientific acceptance remains root-owned');a=ap.parse_args()
    root=a.root.resolve();out=a.output or root/'delegated/cost_package';out.mkdir(parents=True,exist_ok=True)
    ledger=build(root);inv=packaging_inventory(root)
    if a.finalized:
        import importlib.util
        spec=importlib.util.spec_from_file_location('r8_package_gate',root/'delegated/cost_package/package_release.py');gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate);gate.scored_truth_gate(root)
        ledger['status']='FINAL_RECORDED_COSTS_WITH_EXPLICIT_GAPS';ledger['scientific_completion']='Root finalized authorized work; cost gaps retained, no exact full-cost or scientific-success claim.'
    ledger_name='COST_LEDGER_FINAL.json' if a.finalized else 'COST_LEDGER_IN_PROGRESS.json'
    for name,value in [(ledger_name,ledger),('PUBLICATION_INVENTORY.json',inv)]:
        (out/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    with (out/'PUBLICATION_INVENTORY.tsv').open('w') as fh:
        writer=csv.DictWriter(fh,fieldnames=['path','bytes','extension','category','reason'],delimiter='\t');writer.writeheader();writer.writerows(inv['files'])
    text=['# R8成本与发布机械检查（进行中）','',f"本地最终主预测回执 {ledger['primary_completed_receipts']}/{ledger['frozen_primary_expected']}。仅扫描回执与文件元数据；未读GT，未评分，未运行反演。",'',
          '任务wall求和是工作量，不能称项目历时。不同硬件、物理模型、矩阵维度分组，不给统一LU等价费用。已知CPU和包含单独标记的失败进程下界；缺口不当零。','',
          '|类别|分支|硬件|行数|已知CPU秒|未知CPU行|任务wall和秒|','|---|---|---|---:|---:|---:|---:|']
    for r in ledger['grouped_known_costs']:text.append(f"|{r['category']}|{r['family']}|{r['hardware']}|{r['rows']}|{r['known_cpu_seconds']:.6f}|{r['cpu_unknown_rows']}|{r['sum_task_wall_seconds']:.6f}|")
    text+=['','## 额外开发队列','']+[f"- {name}: {v['sealed_receipts']}/{v['expected']}外层回执；{v['scope']}" for name,v in ledger.get('additional_cohorts',{}).items()]
    text+=['','## 明确缺口','']+[f"- {r['source']}: {r['gap']}" for r in ledger['gaps']]
    text+=['','## 失败清单','']+[f"- {r['source']}" for r in ledger['failures']]
    text+=['','## 发布提案','', '|去向|文件数|字节|','|---|---:|---:|']+[f"|{k}|{v['files']}|{v['bytes']}|" for k,v in sorted(inv['groups'].items())]
    text+=['','按R7约600MiB非压缩上限设计分包，具体成员见PUBLICATION_INVENTORY.json。尚未打包/上传，不改gitignore。sealed_truth先HOLD，第三方原测量、原Chat、本地全文下载、临时传输包、SSH身份和本地Python依赖全部排除。Git保留完整实际代码、配置、报告、最终数组、成本回执和失败证据；大history/budget/training数据独立release分包。','',
           '封存后重跑同一命令刷新。输出始终保留IN_PROGRESS标签，完整性和最终科学判断由主线程确认；不凭文件存在自动宣布完成。']
    (out/'COST_GAPS_AND_PACKAGING_ZH.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(dict(status=ledger['status'],primary_completed=ledger['primary_completed_receipts'],expected=ledger['frozen_primary_expected'],ledger_rows=len(ledger['rows']),packaging=inv['groups'],proposed_packages=len(inv['proposed_packages']))))

if __name__=='__main__':main()
