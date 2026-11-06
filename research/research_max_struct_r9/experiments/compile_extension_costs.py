"""Read-only R9 campaign cost union; metadata only, no GT/arrays/solves.

Canonical producer ledgers replace their fallback rows, never supplement them.
CPU is additive only across disjoint processes/work intervals. GPU events,
nested setup, per-arm subtimers and reused continuations are not second charges.
"""
import time
BEGIN_CPU,BEGIN_WALL=time.process_time(),time.perf_counter()
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,platform

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
TARGET=ROOT/'audit/EXTENSION_COSTS.json';MARKDOWN=ROOT/'COST_LEDGER_EXTENSION.md'
MAC='Mac_AppleSilicon';WIN='XINAN_Windows_i9';UNKNOWN='UNKNOWN_PLATFORM'

def utc():return datetime.now(timezone.utc).isoformat()
def host(value):
    text=str(value).lower()
    if any(x in text for x in ('mac','darwin','apple')):return MAC
    if any(x in text for x in ('win','xinan')):return WIN
    return UNKNOWN
def number(value):return value if isinstance(value,(int,float)) and not isinstance(value,bool) else None
def first_number(record,keys):
    for key in keys:
        value=number(record.get(key))
        if value is not None:return value
    return None
def relative(p):return str(Path(p).relative_to(REPO))
def duration(value):
    try:return sum(float(v)*60**i for i,v in enumerate(reversed(str(value).split(':'))))
    except (TypeError,ValueError):return None

class Compiler:
    def __init__(self):self.sources={};self.gaps=[]
    def read(self,path,scope,required=False):
        p=Path(path)
        if not p.is_absolute():p=ROOT/p
        if 'sealed_truth' in p.parts or p.suffix not in ('.json','.jsonl'):raise ValueError('Metadata-only reader refuses this path')
        key=relative(p)
        if not p.exists():
            if required:self.gaps.append(dict(source=key,CPU_seconds=None,scope=scope,reason='Required producer metadata not yet available'))
            return None
        try:
            stat=p.stat();text=p.read_text();data=[json.loads(x) for x in text.splitlines() if x.strip()] if p.suffix=='.jsonl' else json.loads(text)
            self.sources[key]=dict(path=key,bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,scope=scope)
            return data
        except (OSError,ValueError) as exc:
            self.gaps.append(dict(source=key,CPU_seconds=None,scope=scope,reason='Metadata snapshot unavailable/being written',error=str(exc)))
            return None
    def component(self,name,authority,state='COMPLETE_RECORDED_WORK'):
        return dict(id=name,authority=authority,state=state,charges=[],counter_views=[],nonadditive_views=[],gaps=[])
    def charge(self,component,identity,source,platform_name,cpu,scope,wall=None,bound='known_lower_bound',**extra):
        component['charges'].append(dict(id=identity,source=source,platform=platform_name,CPU_seconds=number(cpu),
            CPU_semantics=bound,scope=scope,wall_seconds=number(wall),wall_semantics='Retained source interval only; not added as project makespan',**extra))

    def phase1(self):
        name='audit/ALL_COSTS.json';d=self.read(name,'Canonical PhaseI union; all nested task/component costs already included',True)
        c=self.component('phaseI','canonical',state='IN_PROGRESS' if d is None else 'COMPLETE_RECORDED_WORK')
        if d is None:return c
        launchers=[]
        for queue in d.get('queue_makespans_and_nonoverlap',{}).get('records',[]):
            qpath=REPO/queue['source'];q=self.read(qpath,'PhaseI host provenance only; queue/child costs already in canonical total')
            for row in (q or {}).get('jobs',[]):
                if row.get('command'):launchers.append(row['command'][0])
        established=bool(launchers) and all(str(x).startswith(('/Volumes/','/opt/homebrew/')) for x in launchers)
        self.charge(c,'phaseI-canonical',relative(ROOT/name),MAC if established else UNKNOWN,d['actual_measured_lower_bound']['known_CPU_seconds'],
            '113 Mac inverse receipts + diagnostics/generation/analysis/measured failures/compiler invocations, once; historical R8 excluded',
            platform_evidence=dict(kind='Local PhaseI launcher paths in retained queue metadata',observed_launches=len(launchers),all_Mac_volume_or_Homebrew_paths=established))
        if not established:c['gaps'].append(dict(CPU_seconds=None,reason='PhaseI host not fully established by retained launchers; declared Mac cohort requires review'))
        c['counter_views'].append(dict(source=relative(ROOT/name),dimension='scalar_TM_complex_field_576',
            counters=d['actual_measured_lower_bound'].get('physical',{}),scope=d.get('counter_semantics')))
        c['nonadditive_views']=[dict(source=relative(ROOT/name),excluded='groups/inverse_receipts/compiler_invocations/standaloneweak/portfolio recipes',
            reason='Already in canonical total or attributed standalone method recipes, not additional executed work')]
        c['gaps']+=d.get('unmeasured_costs',[]);c['source_status']=d.get('status');return c

    def small3d(self):
        name='RESULTS/vector_extension/PHASE2_COMPLETE.json';d=self.read(name,'Canonical small3D29 plus migration partial/validation/generation/report union',True)
        c=self.component('small_vector3D','canonical','IN_PROGRESS' if d is None else 'COMPLETE_RECORDED_WORK')
        if d is None:return c
        for h,cpu in d['retained_total_CPU_lower_bounds_by_hardware'].items():
            self.charge(c,'small3D-'+h,relative(ROOT/name),host(h),cpu,
                'Producer deduplicated platform lower bound; includes whole workers and own measured ancillary costs/migration attempts')
        for h,v in d['scientific_workers_by_hardware'].items():
            c['counter_views'].append(dict(source=relative(ROOT/name),platform=host(h),dimension='vector_DDA_complex_field_1536',
                counters={k:v[k] for k in ('physical_LU','physical_source_RHS','physical_adjoint_RHS')},scope='Completed scientific workers only'))
        c['counter_views'].append(dict(source=relative(ROOT/name),dimension='mixed_1536_3000_3408_NOT_POOLED',
            scope=d.get('validation_generation_scope'),counters=d['Mac_incremental_ancillary'].get('aborted_available_physics_lowerbounds',{})))
        c['nonadditive_views'].append(dict(source=relative(ROOT/name),excluded='scientific_workers_by_hardware, Mac_incremental_ancillary, XINAN_ancillary, rows, nested generation/Mie receipts',
            reason='All already included in retained platform total; do not add nested setup/source costs again'))
        c['gaps']=[dict(CPU_seconds=None,reason=x) for x in d.get('missing_or_bounded_costs',[])];return c

    def mie_igt(self):
        name='RESULTS/vector_extension/cell_integrated_mie16/COST.json';d=self.read(name,'Independent one-LU IGT physical diagnostic, outside old29 ledger',True)
        c=self.component('small_sphere_IGT','canonical','IN_PROGRESS' if d is None else 'COMPLETE_RECORDED_WORK')
        if d is None:return c
        self.charge(c,'mieIGT-process',relative(ROOT/name),WIN,d.get('remote_wholeprocess_CPU_seconds'),'Whole remote diagnostic process, checks/setup/oneLU/exports',d.get('remote_script_wall_seconds'))
        self.charge(c,'mieIGT-report',relative(ROOT/name),MAC,d.get('local_postprocess_CPU_seconds'),'Local metadata/figure/report only',d.get('postprocess_wall_seconds'))
        c['counter_views'].append(dict(source=relative(ROOT/name),dimension='vector_pulse_collocation_complex_field_3408',
            counters=dict(physical_LU=d.get('new_forward_LU'),source_RHS=d.get('new_source_RHS')),scope='Single new sphere solve; old CM+RR was reused, not rerun'))
        c['nonadditive_views'].append(dict(source=relative(ROOT/name),excluded='measured_new_CPU_seconds and nested kernel timers',reason='Sum/subsets of the two platform charges'))
        c['gaps'].append(dict(CPU_seconds=d.get('transport_helper_CPU_seconds'),reason=d.get('transport_scope')));return c

    def classic(self):
        canonical='RESULTS/classic_extension_xinan/PHASE2_COST_LEDGER.json'
        d=self.read(canonical,'Canonical classic union replaces every classic fallback row')
        if d is not None:
            c=self.component('classic2D','canonical');c['source_status']=d.get('status')
            views=d.get('by_platform',{});accepted=0
            for h,v in views.items():
                if not isinstance(v,dict):continue
                cpu=first_number(v,('known_CPU_seconds','actual_new_known_CPU_seconds'))
                self.charge(c,'classic-canonical-'+h,relative(ROOT/canonical),host(h),cpu,'Canonical whole classic platform union; categories are nested, not added again')
                accepted+=cpu is not None
            if not accepted:
                c['state']='IN_PROGRESS';c['gaps'].append(dict(CPU_seconds=None,reason='Canonical appeared without explicit per-platform totals; do not guess from nested timers'))
            total=number(d.get('actual_new_known_CPU_seconds'))
            if total is not None and accepted and abs(sum(q['CPU_seconds'] for q in c['charges'] if q['CPU_seconds'] is not None)-total)>1e-6:
                c['state']='IN_PROGRESS';c['gaps'].append(dict(CPU_seconds=None,reason='Canonical per-platform totals do not match declared deduplicated total; review schema'))
            if any(t in str(d.get('status','')).upper() for t in ('RUNNING','IN_PROGRESS','PENDING')):c['state']='IN_PROGRESS'
            c['counter_views'].append(dict(source=relative(ROOT/canonical),dimension='scalar_TM_complex_field_4096',
                counters=d.get('all_inverse_physical_counts',d.get('physical_counts',{})),scope='25 inverse executions; benchmark is separate producer view'))
            if d.get('CUDA_callback_benchmark_separate') is not None:
                c['counter_views'].append(dict(source=relative(ROOT/canonical),dimension='scalar_TM_complex_field_4096_benchmark_separate',
                    counters=d['CUDA_callback_benchmark_separate'],scope='Producer benchmark view only; its CPU is already in Windows platform total, GPUevent is not hostCPU or totalGPUelapsed'))
            c['nonadditive_views'].append(dict(source=relative(ROOT/canonical),excluded='actual_new_known_CPU_seconds plus all original classic task/partial/benchmark/transport receipts',
                reason='Canonical Mac+Windows platform totals already cover all of them; GPU event is a nonadditive device diagnostic'))
            c['gaps']+=d.get('unmeasured_costs',d.get('unknowns',[]));return c

        c=self.component('classic2D','fallback_latest_task_and_process_snapshots','IN_PROGRESS')
        partial='delegated/classic_extension/MAC_PARTIAL_EXECUTION.json';d=self.read(partial,'Mac9 partial aggregate includes all7wait4+2internal lower bounds')
        if d:
            self.charge(c,'classic-Mac9',relative(ROOT/partial),MAC,d.get('known_child_CPU_seconds_sum'),'Whole nine Mac child executions, partial platform cohort; no auxiliary/setup CPU readdition')
            self.charge(c,'classic-Mac-dispatcher',relative(ROOT/partial),MAC,duration(d.get('queue_ps_CPU_time_at_stop')),'Dispatcher ps observed lower bound only; excludes child CPU')
            c['counter_views'].append(dict(source=relative(ROOT/partial),dimension='scalar_TM_complex_field_4096',counters=d.get('physical_counts',{}),scope='Mac9 only; nested auxiliary timers excluded'))
            c['gaps'].append(dict(source=relative(ROOT/partial),CPU_seconds=d.get('queue_CPU_exact'),reason=d.get('unmeasured')))
        prep='RESULTS/classic_extension/PREPARATION_RECEIPT.json';d=self.read(prep,'Own classic manifest/preparation process; not solver setup')
        if d:self.charge(c,'classic-Mac-prepare',relative(ROOT/prep),MAC,first_number(d,('cpu_seconds','CPU_seconds')),'Recorded preparation only',d.get('wall_seconds'))
        # Completed GetProcessTimes rows override their own earlier active snapshots.
        completed={};active={};snapshot='delegated/classic_extension/REMOTE_STATUS_RECEIPT.json';snapshots=self.read(snapshot,'Classic live status/helper history; cumulative task CPU snapshots are not independent charges')
        for poll in snapshots if isinstance(snapshots,list) else ([snapshots] if snapshots else []):
            for row in poll.get('last',[]):completed[row['task_id']]=row
            for row in poll.get('active_progress',[]):
                task=row['task_id'];value=number(row.get('CPU_seconds'))
                if value is not None and (task not in active or value>active[task]['CPU_seconds']):active[task]=row
        process='RESULTS/classic_extension_xinan/queue/TASK_PROCESS_RECEIPTS.json';harvest=self.read(process,'Complete Windows whole-process receipts when fetched; override matching live snapshots')
        for row in harvest if isinstance(harvest,list) else (harvest or {}).get('tasks',[]):completed[row['task_id']]=row
        for task,row in completed.items():
            self.charge(c,'classic-Windows-'+task,relative(ROOT/process) if harvest else relative(ROOT/snapshot),WIN,row.get('child_process_CPU_seconds'),
                'Exact completed Windows child GetProcessTimes; supersedes all snapshots of this same task',row.get('process_wall_seconds'),bound='measured_process_CPU')
        for task,row in active.items():
            if task not in completed:self.charge(c,'classic-Windows-active-'+task,relative(ROOT/snapshot),WIN,row.get('CPU_seconds'),
                'Latest/highest cumulative active-task snapshot only; incomplete tail not measured')
        benchmark='RESULTS/classic_extension_cuda_callback/RECEIPT.json';d=self.read(benchmark,'Single completed whole callback benchmark; CPU/GPU/setup subtimers are nested')
        if d:
            self.charge(c,'classic-callback-benchmark',relative(ROOT/benchmark),WIN,d.get('full_total_internal_process_CPU_seconds'),
                'Whole benchmark process measured interior, startup/output tail unharvested; not another inverse',d.get('full_total_wall_seconds'))
            c['counter_views'].append(dict(source=relative(ROOT/benchmark),dimension='scalar_TM_complex_field_4096',
                counters=dict(benchmark_setup_LU=d.get('total_setup_LU'),GPU_callback_LU=d.get('GPU_callback_LU'),GPU_source_RHS=d.get('GPU_total_source_RHS'),GPU_adjoint_RHS=d.get('GPU_total_adjoint_RHS'),CPU_callback_RHS=40),
                GPU_callback_event_milliseconds=d.get('GPU_callback',{}).get('CUDA_event_milliseconds'),scope='Benchmark separate from inverse counts; GPU event not CPU or total GPU elapsed'))
        # These fields explicitly measure independent transport/controller/helper
        # processes. Never recursively sum arbitrary CPU fields from task receipts.
        paths=sorted((ROOT/'delegated/classic_extension').glob('REMOTE_*_RECEIPT.json'))
        cuda=ROOT/'delegated/classic_extension/CUDA_TRANSPORT_RECEIPT.json'
        if cuda.exists():paths.append(cuda)
        for path in paths:
            d=self.read(path,'Classic transport helper processes only; scientific remote child compute excluded')
            for i,row in enumerate(d if isinstance(d,list) else ([d] if d else [])):
                for field,h in [('local_controller_CPU_seconds',MAC),('local_transport_child_CPU_seconds',MAC),('remote_fetch_helper_CPU_seconds',WIN),('remote_helper_CPU_seconds',WIN)]:
                    if field in row:self.charge(c,f'classic-helper-{path.stem}-{i}-{field}',relative(path),h,row[field],
                        'Own transport/helper CPU only; each invocation once; no GPU event or solver CPU repeated',row.get('remote_helper_wall_seconds') if h==WIN else None)
        failures='delegated/classic_extension/MIGRATION_ENGINEERING_FAILURES.json';d=self.read(failures,'Zero-physics failed engineering attempts, CPU unknown retained')
        if d:c['gaps']+=d.get('attempts',[])
        c['nonadditive_views'].append(dict(source=relative(ROOT/snapshot),excluded='Earlier active snapshots, repeated completed task rows, aux/setup timers, GPUevent; raw Mac task receipts',
            reason='Latest per-task authority and whole process aggregates used once; shared/preparation ownership deferred to final canonical ledger'))
        c['gaps'].append(dict(CPU_seconds=None,reason='Windows16/FIFO/analysis tails and any not-yet-fetched tasks incomplete; canonical classic ledger awaited',completed_windows_tasks=len(completed),active_tasks_observed=len(set(active)-set(completed))))
        return c

    def interface_mie(self):
        name='RESULTS/vector_extension/interface_mie16/COST.json';d=self.read(name,'Independent two known-interface Mie forwards, outside frozen29/original9 canonical',True)
        c=self.component('known_interface_Mie16','canonical_two_processes','IN_PROGRESS' if d is None else 'COMPLETE_RECORDED_WORK')
        if d is None:return c
        for row in d.get('processes',[]):
            self.charge(c,'interface-Mie-'+row['arm'],relative(ROOT/name),host(row.get('platform')),row.get('CPU_seconds'),
                'One separate known-interface forward process; exactly1LU/1RHS; worker setup is included',row.get('wall_seconds'),raw_artifact=row.get('source'))
        self.charge(c,'interface-Mie-parent',relative(ROOT/name),WIN,d.get('parent_CPU_seconds'),d.get('parent_scope'),d.get('parent_wall_seconds'))
        self.charge(c,'interface-Mie-report',relative(ROOT/name),MAC,d.get('local_report_CPU_seconds'),'Local offline report/figure only',d.get('local_report_wall_seconds'))
        c['counter_views'].append(dict(source=relative(ROOT/name),dimension='vector_IGT_complex_field_3408_real_SPD_material',
            counters=dict(physical_LU=d.get('new_LU'),source_RHS=d.get('new_source_RHS')),scope='OFFLINE knownsphereinterface diagnostic, not inverse recovery'))
        c['nonadditive_views'].append(dict(source=relative(ROOT/name),excluded='known_incremental_CPU_seconds and per-worker operator setup/LUT timers',reason='Sum/subsets of disjoint process and reporting charges'))
        c['gaps'].append(dict(CPU_seconds=d.get('unmetered_transport_helper_CPU_seconds'),reason='Small transfer/packing helpers not CPU-metered; not zero'))
        return c

    def original3d(self):
        base='RESULTS/original_vector3d/';c=self.component('original_electrical_size3D','producer_CM_union_plus_disjoint_IGT_processes')
        selected=None;d=None
        for name in ('PHASE2_COST_LEDGER.json','ACCOUNTING_FINAL.json','ACCOUNTING_CURRENT.json','ACCOUNTING_BEFORE_PHYSICS_DIAG.json'):
            candidate=self.read(base+name,'Original3D canonical deduplicated process/work ledger')
            if candidate is not None:selected=base+name;d=candidate;break
        if d is None:c['state']='IN_PROGRESS';c['gaps'].append(dict(CPU_seconds=None,reason='Original3D canonical ledger missing'));return c
        covered=set();rows=d.get('rows',[])
        # Current producer contract: canonical rows are disjoint process totals.
        for i,row in enumerate(rows):
            artifact=row.get('artifact',row.get('source','unknown'));covered.add(artifact)
            if row.get('additive') is False or row.get('included_in'):
                c['nonadditive_views'].append(dict(source=relative(ROOT/selected),excluded=artifact,
                    reason='Producer explicitly marks nested/nonadditive; included_in='+str(row.get('included_in')),scope=row.get('scope'),
                    recorded_CPU_seconds_NOT_ADDED=row.get('CPU_seconds'),recorded_wall_seconds_NOT_ADDED=row.get('wall_seconds'),recorded_counts_NOT_POOLED=row.get('counts')))
                continue
            hp=host(row.get('platform'))
            self.charge(c,'original-CM-'+artifact,relative(ROOT/selected),hp,row.get('CPU_seconds'),row.get('scope','Producer disjoint process row'),row.get('wall_seconds'),
                raw_artifact=base+artifact,operation=row.get('operation'),recorded_status=row.get('status'),terminal_ledger_status=row.get('terminal_ledger_status'))
            c['counter_views'].append(dict(source=relative(ROOT/selected),artifact=artifact,dimension=row.get('vector_dimension','producer_row_scope_no_cross_dimension_LU_equivalence'),
                counters=row.get('counts',dict(completed_primal_RHS=row.get('completed_primal_RHS'))),scope=row.get('operation',row.get('scope'))))
        if not rows:
            c['state']='IN_PROGRESS';c['gaps'].append(dict(CPU_seconds=None,source=relative(ROOT/selected),reason='Canonical schema has no explicit disjoint platform rows; refusing to split a mixed-platform total'))
        c['gaps']+=[dict(CPU_seconds=None,reason=x) for x in d.get('unknowns',[])]
        c['nonadditive_views'].append(dict(source=relative(ROOT/selected),excluded='recorded_CPU_seconds_deduplicated, benchmark subtimers, shared/wiring inside616CPUpartial, layout mainbody, resumed old52source costs',reason=d.get('exclusions')))
        if selected.endswith('ACCOUNTING_CURRENT.json'):
            c['state']='IN_PROGRESS'
            c['gaps'].append(dict(CPU_seconds=None,reason='Producer CURRENT snapshot; shared nine-endpoint resident process still running, its unknown total is not zero. ACCOUNTING_FINAL awaited.'))
        # Sole new authorized physical diagnostic: each PROCESS_COST includes
        # both grids. The per-grid preparation/solve files are nested views.
        if selected.endswith('ACCOUNTING_BEFORE_PHYSICS_DIAG.json'):
            c['state']='IN_PROGRESS'
            c['gaps'].append(dict(CPU_seconds=None,reason='Final original3D ACCOUNTING_FINAL.json awaited: old CM costs plus IGT two-process diagnostic and shared nine-endpoint resident cohort (crossgrid5 +matchedgrid4). Do not treat old5 path as final authority.'))
            live=self.read(base+'STATUS.json','Original3D current writer/scientific state, no additional compute charge')
            if live is not None:c['source_execution_status']=live.get('status')
            for mode in ('prepare','solve'):
                name=base+'igt_source59/'+mode+'_PROCESS_COST.json';extra=self.read(name,'Whole original3D IGT '+mode+' process, grids30/60 included once')
                if extra is None:
                    c['state']='IN_PROGRESS';c['gaps'].append(dict(source=relative(ROOT/name),CPU_seconds=None,reason='Authorized original3D IGT process not yet terminal/fetched'));continue
                cpu=first_number(extra,('CPU_seconds','wholeprocess_CPU_seconds','process_CPU_seconds','full_total_internal_process_CPU_seconds'))
                self.charge(c,'original-IGT-'+mode,relative(ROOT/name),host(extra.get('platform','XINAN')),cpu,
                    'Whole IGT '+mode+' process; per-grid30/60 costs and GPUevents are nested',first_number(extra,('wall_seconds','wholeprocess_wall_seconds','process_wall_seconds')))
                if any(t in str(extra.get('status','')).upper() for t in ('RUNNING','IN_PROGRESS')):c['state']='IN_PROGRESS'
                c['counter_views'].append(dict(source=relative(ROOT/name),dimension='vector_30cubed_81000_and_60cubed_648000_NOT_POOLED',counters=extra.get('counts',{}),scope='Retain per-grid meaning; no LU speed metric across grids'))
            # A future fresh all-arm process is additional only when not already
            # named in the producer canonical union. Its own arm subtimers and
            # shared setup are never added separately.
            parent=base+'invert_CUDA_COST.json';extra=self.read(parent,'Future whole fresh original3D inverse cohort, includes all arms/shared setup')
            if extra is not None and 'invert_CUDA_COST.json' not in covered:
                self.charge(c,'original-fresh-inverse-cohort',relative(ROOT/parent),host(extra.get('platform','XINAN')),extra.get('CPU_seconds'),
                    'Whole fresh inverse process; SHARED_CUDA_PRECOMPUTE_COST and invert_arm_COST nested',extra.get('wall_seconds'))
                if any(t in str(extra.get('status','')).upper() for t in ('RUNNING','IN_PROGRESS')):c['state']='IN_PROGRESS'
        return c

    def residual_transport(self):
        name='RESULTS/original_residual_transport/COST.json'
        d=self.read(name,'Metadata-only interpretation of already paid endpoint norms',True)
        c=self.component('original_residual_transport','single_metadata_process','IN_PROGRESS' if d is None else 'COMPLETE_RECORDED_WORK')
        if d is not None:
            self.charge(c,'residual-transport-analysis',relative(ROOT/name),host(d.get('platform')),d.get('CPU_seconds'),d.get('scope'),d.get('wall_seconds'))
        return c

    def physical_refinement(self):
        name='RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json'
        d=self.read(name,'Independent three-start physical refinement after frozen nine-arm study',True)
        c=self.component('original3D_physical_refinement','canonical_disjoint_process_rows','IN_PROGRESS' if d is None else 'COMPLETE_RECORDED_WORK')
        if d is None:return c
        c['source_status']=d.get('status')
        for i,row in enumerate(d.get('rows',[])):
            artifact=row.get('artifact',row.get('source',str(i)))
            if row.get('additive') is False or row.get('included_in'):
                c['nonadditive_views'].append(dict(source=relative(ROOT/name),excluded=artifact,
                    reason='Producer nested/attribution view; included_in='+str(row.get('included_in'))))
                continue
            self.charge(c,'physical-refinement-process-'+str(i),relative(ROOT/name),host(row.get('platform')),row.get('CPU_seconds'),
                row.get('scope'),row.get('wall_seconds'),raw_artifact=artifact,operation=row.get('operation'))
            c['counter_views'].append(dict(source=relative(ROOT/name),artifact=artifact,dimension='vector_CM30cubed_complex_field_81000',
                counters=row.get('counts',{}),scope='Physical value/gradient/control; no dimension pooling with older LU counts'))
        if not d.get('rows') or any(t in str(d.get('status','')).upper() for t in ('RUNNING','IN_PROGRESS','PENDING')):c['state']='IN_PROGRESS'
        c['gaps']+=d.get('unknowns',[])
        return c

    def publication(self):
        name='audit/INTERIM_PUBLICATION.json';d=self.read(name,'Publication provenance only, no CPU metering')
        c=self.component('publication_metadata','unmetered')
        if d is not None:
            self.charge(c,'interim-publication',relative(ROOT/name),UNKNOWN,None,'Interim/final commit, review, readback and site helpers unmetered; measured archive/upload excluded; not zero')
            c['source_status']=d.get('status');c['gaps'].append(dict(CPU_seconds=None,reason='Publication/web/network/LLM/energy/money not separately metered'))
        packaging='audit/PACKAGING_RECEIPT.json';packed=self.read(packaging,'One final local archive packaging process; scientific writers already stopped')
        if packed is not None:
            self.charge(c,'final-evidence-packaging',relative(ROOT/packaging),MAC,packed.get('CPU_seconds'),
                'Measured final lossless archive packaging process; no scientific solve; manifest/export tail not separately measured',packed.get('wall_seconds'))
        upload_name='audit/RELEASE_UPLOAD_RECEIPT.json';uploaded=self.read(upload_name,'Measured release upload controller and gh child, excluding archive creation')
        if uploaded is not None:
            self.charge(c,'final-release-upload',relative(ROOT/upload_name),MAC,uploaded.get('CPU_seconds'),
                uploaded.get('scope','Upload controller and child; no scientific solve'),uploaded.get('wall_seconds'))
        return c

def aggregate(components):
    groups={};ids=set()
    for component in components:
        for charge in component['charges']:
            if charge['id'] in ids:raise ValueError('Duplicate charge identity '+charge['id'])
            ids.add(charge['id']);h=charge['platform'];g=groups.setdefault(h,dict(known_CPU_lower_bound_seconds=0.,charge_records=0,measured_charge_records=0,unmeasured_charge_ids=[]))
            g['charge_records']+=1
            if charge['CPU_seconds'] is None:g['unmeasured_charge_ids'].append(charge['id'])
            else:g['known_CPU_lower_bound_seconds']+=charge['CPU_seconds'];g['measured_charge_records']+=1
    for g in groups.values():
        if g['measured_charge_records']==0:g['known_CPU_lower_bound_seconds']=None
    return groups

def markdown(result):
    lines=['# R9 extension cost union','',f'Status: **{result["status"]}**. Generated {result["created_utc"]}. Read-only metadata snapshot; no GT/large-array access or new scientific solves. Existing COST_LEDGER.md and raw receipts remain unchanged.','',
        'CPU sums are retained consumption lower bounds by host platform. They are not cross-platform speedup estimates. Unknown CPU is null, never inferred as zero. Task-wall sums, overlapping queues, paused processes and GPU event times do not become a project makespan.','',
        '|Platform|Known CPU lower bound (s)|Charge records|Unmetered charges|','|---|---:|---:|---:|']
    for h,g in result['by_platform'].items():
        value='null' if g['known_CPU_lower_bound_seconds'] is None else f'{g["known_CPU_lower_bound_seconds"]:.6f}'
        lines.append(f'|{h}|{value}|{g["charge_records"]}|{len(g["unmeasured_charge_ids"])}|')
    lines+=['','|Component|Authority|State|Known CPU lower bound (s)|','|---|---|---|---:|']
    for c in result['components']:
        measured=[x['CPU_seconds'] for x in c['charges'] if x['CPU_seconds'] is not None]
        value=f'{sum(measured):.6f}' if measured else 'null'
        lines.append(f'|{c["id"]}|{c["authority"]}|{c["state"]}|{value}|')
    lines+=['','## Disjoint charge scopes','']
    for c in result['components']:
        lines += [f'### {c["id"]}','']
        for q in c['charges']:
            cpu='null (unmeasured)' if q['CPU_seconds'] is None else f'{q["CPU_seconds"]:.6f}s'
            lines.append(f'- `{q["source"]}` — {q["platform"]}, CPU {cpu}; {q["scope"]}.')
        for v in c['nonadditive_views']:lines.append(f'- **Not added again:** {v.get("excluded")}; {v.get("reason")}.')
        for gap in c['gaps']:lines.append(f'- **Gap/limit:** {json.dumps(gap,ensure_ascii=False)}')
        lines.append('')
    lines+=['## Counts and completion','',
        'Physical counters remain in dimension-tagged source views in audit/EXTENSION_COSTS.json. Scalar576, scalar4096, vector1536/3000/3408 and vector30³/60³ are never collapsed into one LU-equivalent speed metric. FFT calls, current Cholesky, source/adjoint RHS, zero-air RHS and partial GMRES counters retain producer semantics.','',
        'Canonical classic and original ledgers replace their fallback rows. A completed task overrides all of its earlier active CPU snapshots. Small3D29 platform totals already include ancillary/migration costs; IGT sphere control is additional. Original partial616CPU already contains its shared precompute/wiring; the whole layout receipt covers its inner probe. Parent CPU and child CPU are added only when explicitly disjoint.','',
        'Rebuild: `python experiments/compile_extension_costs.py` from this study directory (or use its repository-relative path). When all writers have ended, root can use `--final`; it still reports IN_PROGRESS if required canonical/diagnostic work is missing. The flag does not convert unknown costs into exact totals.','']
    return '\n'.join(lines)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--final',action='store_true');args=parser.parse_args()
    previous=[]
    if TARGET.exists():
        try:previous=json.loads(TARGET.read_text()).get('compiler_invocations',[])
        except (OSError,ValueError):pass
    compiler=Compiler();components=[compiler.phase1(),compiler.small3d(),compiler.mie_igt(),compiler.interface_mie(),compiler.classic(),compiler.original3d(),compiler.residual_transport(),compiler.physical_refinement(),compiler.publication()]
    cpu=time.process_time()-BEGIN_CPU;invocation=dict(created_utc=utc(),CPU_seconds=cpu,wall_seconds=time.perf_counter()-BEGIN_WALL,
        scope='This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve')
    history=previous+[invocation];own=compiler.component('extension_cost_compiler','own_invocation_history')
    own_host=MAC if platform.system()=='Darwin' and platform.machine()=='arm64' else WIN if platform.system()=='Windows' else UNKNOWN
    for i,v in enumerate(history):compiler.charge(own,f'extension-compiler-{i}',relative(Path(__file__)),v.get('platform',own_host),v['CPU_seconds'],v['scope'],v['wall_seconds'])
    invocation['platform']=own_host;components.append(own);groups=aggregate(components)
    running=[c['id'] for c in components if c['state']=='IN_PROGRESS']
    result=dict(schema_version=1,status='IN_PROGRESS' if running else 'FINAL_KNOWN_COST_LOWER_BOUND' if args.final else 'SNAPSHOT_KNOWN_COST_LOWER_BOUND',
        created_utc=utc(),final_requested=args.final,running_or_missing_components=running,by_platform=groups,
        mixed_platform_known_CPU_consumption_lower_bound_seconds=sum(g['known_CPU_lower_bound_seconds'] for g in groups.values() if g['known_CPU_lower_bound_seconds'] is not None),
        mixed_platform_total_scope='Consumption bookkeeping only, never a speedup/performance metric',components=components,
        metadata_sources=list(compiler.sources.values()),metadata_read_gaps=compiler.gaps,compiler_invocations=history,
        wall_aggregation='NOT_SUMMED; raw scope intervals retained; no invented makespan',GPU_aggregation='GPUevent views are not added to hostCPU or labeled totalGPUelapsed',
        monetary_cost=None,energy_cost=None,truth_or_large_array_access='NONE',new_scientific_solves=0,hash_checks='NONE',
        counter_policy='No crossdimension LU/RHS/FFT equivalent-speed aggregation',raw_receipts_modified=False,existing_COST_LEDGER_modified=False)
    TARGET.parent.mkdir(parents=True,exist_ok=True);temp=TARGET.with_suffix('.json.tmp');temp.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n');temp.replace(TARGET)
    temporary_markdown=MARKDOWN.with_suffix('.md.tmp');temporary_markdown.write_text(markdown(result));temporary_markdown.replace(MARKDOWN)
    print(json.dumps(dict(status=result['status'],by_platform=groups,running=running,output=relative(TARGET))))
    return 0

if __name__=='__main__':raise SystemExit(main())
