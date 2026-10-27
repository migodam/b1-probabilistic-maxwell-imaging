"""Reconcile saved evidence only. No physics, hidden data, or new gates."""
from pathlib import Path
import csv,json
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1]
def read(p): return json.loads((ROOT/p).read_text())
def write(p,d): (ROOT/p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def main():
    g0=read('results/G0_DECISION.json');g1=read('results/G1_FACTUAL_SUMMARY.json')
    events=[json.loads(l) for l in (ROOT/'PHYSICS_LEDGER.jsonl').read_text().splitlines() if l]
    costs=defaultdict(lambda:dict(rhs=0,full_FE=0,process_cpu_seconds=0.,summed_job_wall_seconds=0.,jobs=0,failed_jobs=0))
    dense=defaultdict(int)
    for e in events:
        c=costs[e['phase']]
        if e['kind']=='reservation': c['rhs']+=e.get('rhs',0);c['full_FE']+=e.get('factorizations',0)
        if e.get('label')=='job_completion':
            c['process_cpu_seconds']+=e.get('process_cpu_seconds',0)
            c['summed_job_wall_seconds']+=e.get('wall_seconds',0);c['jobs']+=1
            c['failed_jobs']+=int(e.get('exit_code',0)!=0)
        if e.get('kind_override')=='dense_core_factorization':dense[e.get('dense_dimension',0)]+=e.get('dense_factorizations',0)
    for phase in ['G1','G2']: costs[phase]
    historical=read('results/g1_replay_001/RUN_RECEIPT.json')['historical_bank_paid_setup']
    cost=dict(schema='BGGA_FINAL_ACCOUNTING_V1',phase_totals=dict(costs),
        actual_new_physical_paid_requests=dict(rhs=sum(c['rhs'] for c in costs.values()),full_FE=sum(c['full_FE'] for c in costs.values())),
        measured_wrapper_process_cpu_seconds=sum(c['process_cpu_seconds'] for c in costs.values()),
        summed_job_wall_seconds=sum(c['summed_job_wall_seconds'] for c in costs.values()),
        peak_measured_job_RSS_bytes=max(e.get('children_peak_RSS',0) for e in events),
        dense_factorizations_by_dimension=dict(dense), historical_bank_paid_setup_not_new=historical,
        actual_new_G1_G2_physics=dict(full_FE=0,rhs=0),
        no_avoided_physics_demonstrated=True,cold_deployment_speedup='NOT_ESTABLISHED',
        wrapper_authority='job_completion CPU once per job including failed attempts; kernels and nested events not added again',
        reservations='charged before physical calls, including failed/rejected attempts; informational observed_* counters not charged twice',
        accounting_gap='Early source/PDF extraction, static worker reads, a synthetic helper check, shell inspection and publication/network preparation are not fully timed. They are not zero; no physical solves were omitted. No total deployment speed claim.',
        current_summary_job='not yet completed at this snapshot; its final receipt stays in the append-only ledger',
        disk_bytes=sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file()),
        memory='job child RSS high-water on macOS, not sum of arrays or independent-method deployment RSS')
    write('COST_SUMMARY.json',cost)
    rows=[]
    for r in read('results/g0_qualification/RESULTS.json'):
        rows.append(dict(stage='G0_qualification',case=r['case'],method='fixed_vs_fitted',status='PASS' if r['numerical_qualification'] else 'FAIL',whitened_error=r['fixed_vs_interface_fitted_white'],scope='new discrete full-wave; empirical independent mesh'))
    for r in read('results/g0_operators/RESULTS.json'):
        rows.append(dict(stage='G0_operator',case=r['case'],method=r['method'],status=r['strict_algebra_status'],
            whitened_error=r['full_response_difference_whitened_l2'],kernel_cpu_seconds=r['solve_timing']['measured_kernel_cpu_seconds'],scope='query kernel in common process, not standalone cold deployment'))
    with (ROOT/'results/G0_STRESS_TABLE.csv').open() as f:
        for r in csv.DictReader(f): rows.append(dict(stage='G0_stress',**r,status='NUMERICAL_MECHANISM_FIXTURE',scope='known support;12 fixtures/9 physical parents; empirical rank uncertainty retained'))
    extra={}
    with (ROOT/'results/g1_offline_extra/METHOD_EXTRA_SCORES.csv').open() as f:
        er=list(csv.DictReader(f))
        for r in er:
            extra[(r.get('parent_index',r.get('parent')),r.get('repeat_index',r.get('repeat')),r.get('method_id'))]=r
    with (ROOT/'research/delegated/g1_figures/G1_CASE_TABLE.csv').open() as f:
        for r in csv.DictReader(f):
            e=extra.get((r.get('parent_index'),r.get('repeat_index'),r.get('method_id')), {})
            offline={'offline_'+k:v for k,v in e.items() if k not in ['parent_index','repeat_index','method_id']}
            rows.append(dict(stage='G1_replay',**r,**offline,scope='REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE'))
    for method in ['weighted_adjoint_FWI_TR','multistart_FWI','support_birth_death_exchange','ordinary_Schur_cache','ordinary_goal_ROM','BGGA_reopen','locked_30_parent_confirmation']:
        rows.append(dict(stage='G2',method=method,status='NOT_RUN_NO_NET_SAVING_MECHANISM',scope='new cold deployment/unknown-support campaign not started'))
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (ROOT/'FULL_TABLE.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    decision=dict(schema='BGGA_FINAL_SCIENTIFIC_DECISION_V1',
        overall='NO_GO_CURRENT_FULL_PORT_AND_DEFAULT_INTERVAL_ACCELERATION',
        EXACT_DISCRETE_ASSIMILATION='QUALIFIED_ON_TESTED_COMMON_SPACE',
        TARGET_INVARIANT_PROPAGATION_CACHE='QUALIFIED_EMPIRICALLY_ON_12_TESTED_CONFIGURATIONS_NOT_UNIFORM_PRIOR',
        ANGLE_ONLY_SAFETY='REJECTED_AS_GENERAL_RULE; rank-unresolved rows retained',
        TASK_CONDITIONAL_PROPAGATION_IDENTIFIABILITY='EXPERIMENTALLY_SUPPORTED_NARROW_TX0_RX8_CONTRACT_ONLY',
        FINITE_PRIOR_POSTERIOR='CONDITIONAL_CACHED_NUMERICAL_MODEL_CORRECT_FOR_M0_M4_M5_M6; PHYSICAL_POSTERIOR_NOT_CERTIFIED',
        INTERVAL_RAW_F64='FAIL_POSITIVE_BOUND_VIOLATIONS_RETAINED',
        INTERVAL_NORMALIZATION_REPAIR='CONDITIONAL_FIXED_LOG_POTENTIAL_OUTWARD_NORMALIZATION; 168_CHECKPOINTS_NO_OBSERVED_VIOLATION',
        JOINT_DEPENDENCE='PRESERVED_BY_M0_M4_M5_M6; SINGLE_WALL_M1_M3_FAIL_11_OF_24; M2_FAIL_24_OF_24',
        A0_UNKNOWN_SUPPORT='NOT_RUN_NEW_BGGA_A0; INHERITED_WSLAM_A0_ONLINE_NOT_REACHED_NOT_REPAIRED',
        COLD_COST='NO_GO_CURRENT_IMPLEMENTATION; NO_FRESH_INDEPENDENT_COLD_DEPLOYMENT_ADVANTAGE',
        CONTINUUM_CERTIFICATE='NOT_ESTABLISHED',
        G2='NOT_RUN_SCIENTIFIC_STOP',
        no_new_locked_data=True,new_NN=False,new_3D=False,
        G1_method_gates=[dict(method=r['method_id'],joint_TV_mean=r['joint_TV_parent_mean'],
            joint_probability_and_risk_passes_out_of_24=r['all_probability_and_risk_gate_PASS_cases_out_of_24'],
            queries_each=2950) for r in g1['method_facts']],
        G0_qualification_max=g0['qualification_whitened_max'],
        physical_paid_requests=cost['actual_new_physical_paid_requests'],
        next_step='Retain task-conditional propagation identifiability. A new cheap likelihood/propagation envelope must beat ordinary Schur/goal-ROM before any new BGGA speed campaign.')
    write('GATE_DECISION.json',decision)
    # Metadata only; the frozen numerical values and prior are never changed.
    c=read('CONTRACT.json');c['status']='FROZEN_DESIGN_EXECUTED_G0_G1_REPLAY_G2_STOPPED';c['resource_grant_record']='AUTHORIZATION.json';write('CONTRACT.json',c)
    gates=read('GATES.json');gates['status']='Frozen before science; numerical values unchanged; results in GATE_DECISION.json';write('GATES.json',gates)
    print(json.dumps(dict(cost=cost['phase_totals'],rows=len(rows),decision=decision['overall'])))
if __name__=='__main__':main()
