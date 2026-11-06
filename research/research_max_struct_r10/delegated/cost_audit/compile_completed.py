"""Whitelist-only historical/DEV JSON cost compilation; no arrays or fresh work."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import json,time
O=Path(__file__).resolve().parent;ROOT=O.parents[1]
def read(p):return json.loads(p.read_text())
def stamp(p):
    s=p.stat();return dict(path=str(p.relative_to(ROOT)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def main():
    cpu,wall=time.process_time(),time.perf_counter();rows=[];sources=[];unknown=[]
    def add(ident,group,p,CPU=None,wall_time=None,physical_LU=None,physical_RHS=None,physical_dimension=None,
            shifted_LU=None,shifted_RHS=None,shifted_dimension=None,small=None,status=None,truth=None,nested=None,scope='',extra=None):
        sources.append(stamp(p));small=small or {}
        r=dict(id=ident,group=group,status=status,source=str(p.relative_to(ROOT)),recorded_CPU_seconds=CPU,recorded_process_wall_seconds=wall_time,
            physical_LU=physical_LU,physical_complex_RHS=physical_RHS,physical_LU_dimension=physical_dimension,
            shifted_auxiliary_LU=shifted_LU,shifted_auxiliary_complex_RHS=shifted_RHS,shifted_LU_dimension=shifted_dimension,
            denominator_small_LU=small.get('denominator_small_LU'),denominator_small_LU_dimension_sum=small.get('denominator_small_LU_dimension_sum'),
            denominator_small_triangular_complex_RHS=small.get('denominator_small_triangular_complex_RHS'),truth_access=truth,
            cost_aggregation_role='DISJOINT_RECORDED_TOP_LEVEL_WINDOW',scope=scope,nested_views=nested or {},extra=extra or {})
        rows.append(r)
        if CPU is None:unknown.append(dict(id=ident,field='CPU_seconds',reason='Not recorded; not zero'))
        if wall_time is None:unknown.append(dict(id=ident,field='wall_seconds',reason='Not recorded; not zero'))
        return r
    # Historical69-run449-sparse-state energy audit, all failed/interrupted attempts.
    for suffix in ('development','confirmation','vector','failed_attempt1_dtype','failed_attempt3_confirmation_shape'):
        p=ROOT/'delegated/energy_audit'/suffix/'EXECUTION_RECEIPT.json';d=read(p)
        add('energy_'+suffix,'historical_energy',p,d.get('cpu_seconds'),d.get('wall_seconds'),d.get('new_physical_LU'),d.get('new_complex_RHS'),1536 if suffix=='vector' else 576,
            status=d['status'],truth=d['truth_access'],nested=dict(spectrum_and_joint_space_timing_row_count=len(d.get('reference_linear_algebra_cost',[])),
                spectrum_and_joint_space_known_wall_seconds=sum(v.get('wall_seconds',0) for v in d.get('reference_linear_algebra_cost',[])),raw_field='reference_linear_algebra_cost in source receipt'),
            scope='Offline geometry/energy script measured window, includes its setup/forward/decomposition/serialization; Python import/lifetime overhead not independently measured.',
            extra=dict(historical_runs=d.get('runs'),saved_states=d.get('saved_states'),failure_count=d.get('failure_count')))
    p=ROOT/'delegated/energy_audit/failed_attempt2_vector_serialization_COST.json';d=read(p)
    add('energy_interrupted_vector','historical_energy',p,d.get('observed_CPU_seconds'),d.get('observed_elapsed_seconds'),
        status=d['status'],truth='EXPOSED_OFFLINE_DIAGNOSTIC',scope='Own interrupted diagnostic process; ps-observed CPU/elapsed, no recoverable factorization counters.',extra=d)
    unknown.append(dict(id='energy_interrupted_vector',field='physical_LU/RHS',reason='Stopped process counters unavailable; known physical totals are lower bounds'))
    # Completed scoring/figure invocations are separate windows, not nested inside inverse receipts.
    for suffix,cpu_field,wall_field in [('RECEIPT.json','scoring_CPU_seconds','scoring_wall_seconds'),('FIGURE_RECEIPT.json','CPU_seconds','wall_seconds')]:
        p=ROOT/'delegated/denominator_score'/suffix;d=read(p)
        add('historical_dev_'+('score' if suffix=='RECEIPT.json' else'figures'),'historical_postprocessing',p,d.get(cpu_field),d.get(wall_field),d.get('new_offline_physical_LU'),d.get('new_offline_complex_RHS'),576,
            status=d['status'],truth='EXPOSED_OFFLINE_DIAGNOSTIC',scope='Six-candidate/all-reference saved-array scoring or sparse historical energy/material figure production; no inverse rerun.')
    p=ROOT/'RESULTS/denominator_dev/SELF_RETENTION_SCORE.json';d=read(p)
    add('self_retention_offline_score','historical_postprocessing',p,d.get('CPU_seconds'),d.get('wall_seconds'),0 if d.get('new_physics_calls')==0 else None,0 if d.get('new_physics_calls')==0 else None,576,
        status=d['status'],truth=d['truth_access'],scope='Separate historical self-retention saved-material scoring; physical zero follows explicit new_physics_calls=0.')
    #32 successful fixed-current updates plus16 paid updates in retained failed attempt.
    for suffix in ('denominator_fixed_state','denominator_fixed_state_failed_missing_operator_01'):
        p=ROOT/'RESULTS'/suffix/'RECEIPT.json';d=read(p);c=d['physical_cost']
        add(suffix,'fixed_state_updates',p,d.get('cpu_seconds'),d.get('wall_seconds'),c.get('LU'),c.get('actual_complex_RHS'),576,small=d.get('auxiliary_cost'),
            status=d['status'],truth=d['truth_access'],scope='Entire recorded fixed-current counterfactual update invocation; basis setup/readout timers and all counters are nested.',
            nested=dict(basis_setup=d.get('basis_setup'),auxiliary_cost=d.get('auxiliary_cost'),physical_cost=c),extra=dict(observations=d.get('observations'),failure=d.get('traceback')))
    # Exactly8 completed full historical DEV candidates, including2 self-retained.
    paths=sorted((ROOT/'RESULTS/denominator_dev').glob('*/RECEIPT.json'))
    if len(paths)!=8:raise RuntimeError('Expected exactly8 historical DEV candidate receipts; do not silently expand scope')
    for p in paths:
        d=read(p)
        if d['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError('Requested completed DEV candidate not sealed')
        c=d['physical_model_cost'];a=d.get('auxiliary_execution_receipt') or {};ops=a.get('operation_counts') or {};setup_path=p.parent/'DENOMINATOR_SETUP.json';setup=read(setup_path);sources.append(stamp(setup_path))
        add(p.parent.name,'complete_DEV_candidates',p,d.get('cpu_seconds'),d.get('wall_seconds'),c.get('LU'),c.get('actual_complex_RHS'),576,small=ops,
            status=d['status'],truth=d.get('truth_access'),scope=d.get('cost_scope','Inherited recorded runner window, not external full-child lifetime'),
            nested=dict(setup=setup,auxiliary_CPU_seconds=d.get('auxiliary_cpu_seconds'),auxiliary_wall_seconds=d.get('auxiliary_wall_seconds'),
                auxiliary_execution_CPU_seconds=a.get('cpu_seconds'),auxiliary_execution_wall_seconds=a.get('wall_seconds'),
                corrector=d.get('corrector'),physical_operation_CPU=d.get('physical_operation_cpu'),physical_operation_wall=d.get('physical_operation_wall')),
            extra=dict(execution_platform=d.get('execution_platform'),device=d.get('device'),physical_source_RHS=c.get('source_RHS'),physical_adjoint_RHS=c.get('adjoint_RHS')))
    # Retained interrupted onehop optimizer controls: wrong transported maxiter60,
    # not valid LBFGS8 evidence and never scored. Corrected running controls excluded.
    p=ROOT/'RESULTS/onehop_solver_control/FAILED_ATTEMPTS.json'
    if p.exists():
        attempts=read(p)
        for i,d in enumerate(attempts):
            c=d.get('physical_cost') or {}
            add('onehop_control_failed_'+str(i+1),'onehop_optimizer_controls_failed',p,d.get('CPU_seconds'),d.get('wall_seconds'),c.get('LU'),c.get('actual_complex_RHS'),576,
                status=d['status'],truth=d.get('truth_access'),scope='Explicit retained interrupted own process; actual inner maxiter60 due dispatch mismatch, not an8-step control. No material scoring.',
                extra=dict(source_receipt=d.get('path'),reason=d.get('reason'),physical_source_RHS=c.get('source_RHS'),physical_adjoint_RHS=c.get('adjoint_RHS')))
    for p in sorted((ROOT/'RESULTS/onehop_solver_control').glob('*_corrected_dispatch/RECEIPT.json')):
        d=read(p);schedule_path=p.parent/'NUMERICAL_SCHEDULE.json';schedule=read(schedule_path);sources.append(stamp(schedule_path))
        if d['status']!='SEALED_PUBLIC_PREDICTION' or schedule['material_inner_maxiter']!=8:raise RuntimeError('Formal onehop control not sealed with actual8 schedule')
        c=d['physical_model_cost'];a=d.get('auxiliary_execution_receipt') or {}
        add(p.parent.name,'onehop_optimizer_controls_completed',p,d.get('cpu_seconds'),d.get('wall_seconds'),c.get('LU'),c.get('actual_complex_RHS'),576,small=a.get('operation_counts'),
            status=d['status'],truth=d.get('truth_access'),scope='Corrected dispatch, sealed historical onehop+LBFGS8 optimizer-only control; not a new metric or fresh confirmation.',
            nested=dict(auxiliary_CPU_seconds=d.get('auxiliary_cpu_seconds'),auxiliary_wall_seconds=d.get('auxiliary_wall_seconds'),auxiliary_execution=a.get('cpu_seconds'),corrector=d.get('corrector')),
            extra=dict(material_inner_maxiter=schedule['material_inner_maxiter'],execution_platform=d.get('execution_platform'),physical_source_RHS=c.get('source_RHS'),physical_adjoint_RHS=c.get('adjoint_RHS')))
    p=ROOT/'RESULTS/onehop_solver_control/SCORE.json'
    if p.exists():
        d=read(p)
        add('onehop_optimizer_control_offline_score','historical_postprocessing',p,d.get('CPU_seconds'),d.get('wall_seconds'),0 if d.get('new_physics_calls')==0 else None,0 if d.get('new_physics_calls')==0 else None,576,
            status=d['status'],truth=d['truth_access'],scope='Separate existing postseal historical control scoring; reused recorded metrics, no scorer rerun by cost audit.')
    # One no-GT shifted-denominator24-state diagnostic; physical calls explicitly0.
    p=ROOT/'delegated/denominator_spectrum/RECEIPT.json';d=read(p);s=d['shifted_auxiliary_cost']
    add('denominator_spectrum24','shifted_spectrum',p,d.get('CPU_seconds'),d.get('wall_seconds'),d['physical_LU'],d['physical_complex_RHS'],576,
        s.get('shifted_auxiliary_576_LU'),s.get('shifted_auxiliary_complex_RHS'),576,d.get('lowrank_cost'),status=d['status'],truth=d['truth_access'],
        nested=dict(onehop_cost=d.get('onehop_cost')),scope=d['scope'],extra=dict(completed_states=d['completed_states']))
    for filename in ('FIGURE_RECEIPT.json','FAILED_START_ATTEMPT.json'):
        p=ROOT/'delegated/denominator_spectrum'/filename;d=read(p)
        add('spectrum_'+filename.removesuffix('.json').lower(),'shifted_spectrum',p,d.get('CPU_seconds'),d.get('wall_seconds'),
            0 if filename=='FAILED_START_ATTEMPT.json' else d.get('new_LU'),0 if filename=='FAILED_START_ATTEMPT.json' else d.get('new_RHS'),576,
            shifted_LU=0 if filename=='FAILED_START_ATTEMPT.json' else None,status=d['status'],truth=d.get('GT_access'),scope='Plotting recorded CSV, or parse failure before imports/computation; no new physics.')
    p=ROOT/'audit/self_retention_wiring.json';d=read(p)
    add('self_retention_wiring','self_retention_wiring',p,d.get('CPU_seconds'),d.get('wall_seconds'),small=d.get('cost'),status=d['status'],truth=d['truth_access'],scope=d['purpose'],
        extra=dict(physical_counter_scope='No physical LU/RHS fields recorded; leave null, do not invent measurements'))
    # Public small-H spectral qualification and its separate earlier exploration.
    # This is not an inverse, an Atlas identity rerun, or a physical forward.
    p=ROOT/'delegated/denominator_dissipation/RECEIPT.json'
    if p.exists():
        d=read(p)
        add('denominator_dissipation_public_diagnostic','denominator_dissipation',p,d.get('CPU_seconds'),d.get('wall_seconds'),d.get('physical_LU'),d.get('physical_RHS'),576,
            status=d['status'],truth=d['truth_access'],scope='Recorded public19/28-dimensional H spectral boundary diagnostic; zero physical LU/RHS; not a floating-matrix certificate.',
            extra=dict(bases=[r['basis'] for r in d.get('rows',[])],initial_exploration_charged_separately=True))
        add('denominator_dissipation_initial_exploration','denominator_dissipation',p,d.get('initial_exploration_CPU_seconds'),d.get('initial_exploration_wall_seconds'),d.get('physical_LU'),d.get('physical_RHS'),576,
            status='INITIAL_EXPLORATION_RECORDED',truth=d['truth_access'],scope='Separate initial public-operator exploration recorded in source; CPU additive to later diagnostic, wall unmeasured.')
    # Literature artifacts carry no CPU/wall billing receipts; preserve failed429 separately.
    for p in sorted((ROOT/'audit/literature').glob('*.json')):
        d=read(p);iserror='error' in d
        rr=add('literature_'+p.stem,'literature',p,CPU=d.get('CPU_seconds',d.get('cpu_seconds')),wall_time=d.get('wall_seconds'),status='FAILED_RETRIEVAL' if iserror else d.get('status','RETRIEVAL_ARTIFACT_RECORDED'),
            truth='NONE',scope='External scholarly retrieval artifact, not a process CPU or provider billing receipt; artifact count is not an API-call count.',
            extra=dict(error=d.get('error'),query_run_records=len(d.get('query_runs',[])),operation_errors=d.get('operation_errors'),provider_money_cost=None))
        rr['cost_aggregation_role']='UNMEASURED_EXTERNAL_RETRIEVAL_ARTIFACT';rr['physical_counter_scope']='NOT_APPLICABLE_TO_RETRIEVAL; no numerical physics counted'
    groups={}
    for group in sorted({r['group'] for r in rows}):
        rr=[r for r in rows if r['group']==group]
        groups[group]=dict(records=len(rr),recorded_CPU_seconds_known_sum=sum(r['recorded_CPU_seconds'] for r in rr if r['recorded_CPU_seconds'] is not None) if any(r['recorded_CPU_seconds'] is not None for r in rr) else None,
            recorded_process_wall_seconds_known_sum=sum(r['recorded_process_wall_seconds'] for r in rr if r['recorded_process_wall_seconds'] is not None) if any(r['recorded_process_wall_seconds'] is not None for r in rr) else None,
            CPU_unknown_records=[r['id'] for r in rr if r['recorded_CPU_seconds'] is None],wall_unknown_records=[r['id'] for r in rr if r['recorded_process_wall_seconds'] is None])
    physical=defaultdict(lambda:dict(known_LU=0,known_complex_RHS=0));shifted=defaultdict(lambda:dict(known_LU=0,known_complex_RHS=0));small=defaultdict(lambda:dict(known_LU=0,known_small_complex_RHS=0))
    for r in rows:
        if r['physical_LU_dimension'] is not None:
            v=physical[str(r['physical_LU_dimension'])];v['known_LU']+=r['physical_LU'] or 0;v['known_complex_RHS']+=r['physical_complex_RHS'] or 0
        if r['shifted_LU_dimension'] is not None:
            v=shifted[str(r['shifted_LU_dimension'])];v['known_LU']+=r['shifted_auxiliary_LU'] or 0;v['known_complex_RHS']+=r['shifted_auxiliary_complex_RHS'] or 0
        if r['denominator_small_LU']:
            rank=r['denominator_small_LU_dimension_sum']/r['denominator_small_LU']
            key=str(int(rank)) if rank.is_integer() else'MIXED_DIMENSIONS'
            v=small[key];v['known_LU']+=r['denominator_small_LU'];v['known_small_complex_RHS']+=r['denominator_small_triangular_complex_RHS'] or 0
    result=dict(status='COMPLETED_HISTORICAL_DEV_COST_AUDIT_ONLY',created_utc=datetime.now(timezone.utc).isoformat(),
        included_scope=['historical69-run energy audit and all retained attempts','separate historical scoring/figures','32 completed fixed-state readouts plus16 failed-attempt paid readouts','8 completed DEV candidate inversions including2 self-retained','2 interrupted onehop optimizer-control attempts with incorrect maxiter60','2 sealed corrected-dispatch onehop+LBFGS8 controls and existing offline control scoring','24 no-GT shifted auxiliary diagnostics and its plot/parse failure','self-retention wiring','public denominator dissipation diagnostic plus separate initial exploration','recorded literature retrieval artifacts including429, dissipation primary sources and Twofold abstract followup'],
        excluded_scope=['fresh2D generation/168 running inversions/transport','fresh3D generation/inversions/transport','fresh scoring not authorized or executed','old R6/R9 reference inverse costs: reused artifacts, not newly consumed here','unmeasured human/model reasoning and import/controller overhead outside recorded windows'],
        rows=rows,group_summaries=groups,known_recorded_CPU_seconds_sum=sum(r['recorded_CPU_seconds'] or 0 for r in rows),
        known_recorded_process_wall_seconds_sum=sum(r['recorded_process_wall_seconds'] or 0 for r in rows),
        wall_scope='Sum of per-invocation recorded wall windows, not end-to-end calendar makespan; concurrent windows may overlap',
        physical_counts_by_matrix_dimension=dict(physical),shifted_auxiliary_counts_by_matrix_dimension=dict(shifted),denominator_small_counts_by_dimension=dict(small),
        unknowns=unknown,totals_complete=False,unknown_policy='Unknown CPU/wall/counters remain null. Known sums are lower bounds, not zero-filled complete expenditure.',
        no_double_count_policy='Only disjoint top-level recorded invocation windows contribute CPU/wall totals. Nested setup/aux/corrector/operation timings and existing aggregate ledgers are evidence views only. Counter classes/dimensions stay separate.',
        sources=list({s['path']:s for s in sources}.values()),bookkeeping_CPU_seconds=time.process_time()-cpu,bookkeeping_wall_seconds=time.perf_counter()-wall,
        bookkeeping_new_physical_calls=0,fresh_truth_read=False,fresh_scorer_executed=False)
    previous=O/'COMPLETED_COSTS.json';builds=O/'BUILD_COSTS.jsonl'
    if previous.exists() and not builds.exists():
        old=read(previous)
        with builds.open('a') as f:f.write(json.dumps(dict(created_utc=old['created_utc'],CPU_seconds=old.get('bookkeeping_CPU_seconds'),wall_seconds=old.get('bookkeeping_wall_seconds'),new_physical_calls=0))+'\n')
    previous.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with builds.open('a') as f:f.write(json.dumps(dict(created_utc=result['created_utc'],CPU_seconds=result['bookkeeping_CPU_seconds'],wall_seconds=result['bookkeeping_wall_seconds'],new_physical_calls=0))+'\n')
    print(json.dumps({k:result[k] for k in ('known_recorded_CPU_seconds_sum','known_recorded_process_wall_seconds_sum','physical_counts_by_matrix_dimension','shifted_auxiliary_counts_by_matrix_dimension','denominator_small_counts_by_dimension','group_summaries')},indent=2))
if __name__=='__main__':main()
