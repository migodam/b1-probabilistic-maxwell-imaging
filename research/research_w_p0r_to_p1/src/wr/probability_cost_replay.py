"""Observation-only cost replay with equal public response-cache rights.

Import is inert: no Book, Family, PDE, observations, or output construction.
main reads completed P1A quadrature rules, y/parent_ids, saved reference and
cost-only ledger fields. The parent launcher charges CPU. No new Maxwell.
--dataset P1A (default) preserves the original route; P1B tests the unchanged
P1A response banks on new y, against B's independently checked continuous_ref.

Both arms cache unique (branch,float.hex(u),float.hex(c)) response nodes once,
before observations are loaded. Neither caches an observation's likelihood or
posterior across passes. Direct integrates its paid order3 CDF rules; the
18-node interpolant integrates full9 and directly truncated CDF9 intervals.
Both call posterior_math for ONE shared gamma, with identical convergence.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from .posterior_math import log_marginal_likelihood
from .prob_reference import interpolate, posterior, cdf_coherence

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/P1A/cost_comparison'
METHODS = ('direct_order3', 'response_interpolation18')
NUISANCE_ORDERS = (11, 21, 41, 65)
NUISANCE_TOLERANCE = 2e-4
PROBABILITY_TOLERANCE = .02
N_OBSERVATIONS = 32
WARM_ORDER_SEED = 20261024


def node_key(branch, u, contrast):
    """Exact physical-cache identity; never round nearby quadrature nodes."""
    return f'{int(branch)}|{float(u).hex()}|{float(contrast).hex()}'


def _json(path):
    return json.loads(Path(path).read_text())


def _write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    temp.replace(path)


def _csv(path, rows):
    if not rows:
        raise ValueError('Refuse empty CSV')
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def timed(callable_, *args):
    wall, cpu = time.perf_counter(), time.process_time()
    value = callable_(*args)
    return value, {'wall_seconds': time.perf_counter()-wall,
                   'cpu_seconds': time.process_time()-cpu}


def _identities(rule):
    return [node_key(b, u, c) for b, u, c in
            zip(rule['branches'], rule['u'], rule['contrast'])]


def load_rules(root=ROOT):
    """Read exactly full3/5/9 and the five CDF rules at each order, no plans."""
    directory = Path(root) / 'results/P1A/rules'
    completion = _json(directory/'COMPLETE.json')
    names = [f'full_{n}' for n in (3, 5, 9)] + [
        f'cdf_{i}_{n}' for i in range(5) for n in (3, 5, 9)]
    if not set(names).issubset(completion['complete_rules']):
        raise ValueError('Completed paid rules are missing required order3/5/9 inputs')
    rules, cuts = {}, []
    for name in names:
        with np.load(directory/f'{name}.npz', allow_pickle=False) as z:
            rule = {k: z[k].copy() for k in
                    ('F', 'branches', 'u', 'contrast', 'x', 'weights')}
        n = len(rule['weights'])
        if rule['F'].shape != (n, 5, 4, 12) or not np.iscomplexobj(rule['F']):
            raise ValueError('Wrong complex response shape: '+name)
        if any(rule[k].shape != (n,) for k in ('branches','u','contrast','x','weights')):
            raise ValueError('Wrong material-coordinate shape: '+name)
        if any(not np.isfinite(a).all() for a in rule.values()):
            raise ValueError('Nonfinite rule: '+name)
        if np.any(rule['weights'] <= 0) or not np.isin(rule['branches'], [0,1]).all():
            raise ValueError('Invalid positive prior masses/branch: '+name)
        if len(set(_identities(rule))) != n:
            raise ValueError('Duplicate nodes within a rule: '+name)
        if name.startswith('full_') and abs(float(rule['weights'].sum())-1) > 1e-12:
            raise ValueError('Full original prior mass must sum to one')
        rules[name] = rule
    centers = [float(np.mean(rules['full_9']['x'][rules['full_9']['branches']==b] -
                             rules['full_9']['u'][rules['full_9']['branches']==b])) for b in (0,1)]
    for i in range(5):
        cut = float(_json(directory/f'cdf_{i}_9.json')['cut'])
        cuts.append(cut)
        original_mass = .5 * sum(float(np.clip((cut-center+.01)/.02,0,1)) for center in centers)
        for order in (3,5,9):
            r = rules[f'cdf_{i}_{order}']
            if abs(float(r['weights'].sum())-original_mass)>2e-12:
                raise ValueError('Truncated rule lost original prior volume')
    if any(b<=a for a,b in zip(cuts,cuts[1:])):
        raise ValueError('CDF checkpoints must be strictly increasing')
    basis = rules['full_3']
    if len(basis['weights']) != 18:
        raise ValueError('Interpolation training rule must have exactly 18 nodes')
    for b in (0,1):
        ids = basis['branches']==b
        u,c = basis['u'][ids],basis['contrast'][ids]
        if len(u)!=9 or len(np.unique(u))!=3 or len(np.unique(c))!=3:
            raise ValueError('Branch-specific training grid must be 3x3')
        if not np.array_equal(np.lexsort((c,u)),np.arange(9)):
            raise ValueError('Saved interpolation cube ordering is not (u,c)')
    return rules, cuts, completion


def unique_node_registry(rules, names):
    """Union nodes with exact floathex identity and equal duplicate F values."""
    registry = {}
    for name in names:
        r = rules[name]
        for i,k in enumerate(_identities(r)):
            if k in registry:
                if not np.array_equal(registry[k]['F'],r['F'][i]):
                    raise ValueError('Same exact public node has different full response')
                registry[k]['rules'].append(name)
            else:
                registry[k] = dict(branch=int(r['branches'][i]),u=float(r['u'][i]),
                    contrast=float(r['contrast'][i]),F=r['F'][i],rules=[name])
    return registry


def build_response_cache(method, rules):
    """Pure arrays, observation-independent; same node/F reuse rights for both arms."""
    if method not in METHODS:
        raise ValueError('Unknown method')
    order = 3 if method=='direct_order3' else 9
    names = [f'full_{order}']+[f'cdf_{i}_{order}' for i in range(5)]
    registry = unique_node_registry(rules,names)
    keys = list(registry)
    position = {k:i for i,k in enumerate(keys)}
    coords = {k:np.array([r[k] for r in registry.values()]) for k in ('branch','u','contrast')}
    if method=='direct_order3':
        responses = np.stack([r['F'] for r in registry.values()])
        physical_training_keys = keys
        response_error = 0.
    else:
        responses = interpolate(rules['full_3'],coords['branch'],coords['u'],coords['contrast'])
        physical_training_keys = _identities(rules['full_3'])
        paid = np.stack([r['F'] for r in registry.values()])
        response_error = float(np.max(np.linalg.norm((responses-paid).reshape(len(keys),-1),axis=1)/.01))
    if responses.shape!=(len(keys),5,4,12) or not np.isfinite(responses).all():
        raise ValueError('Invalid prebuilt response array')
    responses.setflags(write=False)
    maps = {name:np.array([position[k] for k in _identities(rules[name])],int) for name in names}
    return dict(method=method,order=order,names=names,F=responses,maps=maps,
        unique_evaluation_node_keys=keys,unique_physical_training_keys=physical_training_keys,
        response_error_white_on_evaluated_paid_nodes=response_error,
        memory_bytes=int(responses.nbytes+sum(a.nbytes for a in maps.values())),
        interpolation_calls=1 if method==METHODS[1] else 0,
        interpolation_material_rows=len(keys) if method==METHODS[1] else 0)


def marginal_shared_gamma(f,y):
    """Shared gain over all 240 complex channels; union-node quadrature convergence."""
    previous = None
    history = []
    calls = 0
    node_evaluations = 0
    for order in NUISANCE_ORDERS:
        result = log_marginal_likelihood(f,y,n=order,sigma=.01)
        values = np.asarray(result['logZ'],float)
        calls += 1
        node_evaluations += len(f)*order**2
        error = None if previous is None else float(np.max(np.abs(values-previous)))
        history.append({'order':order,'max_node_logZ_change':error})
        if error is not None and error<=NUISANCE_TOLERANCE:
            break
        previous = values
    return values,dict(converged=bool(error is not None and error<=NUISANCE_TOLERANCE),
        history=history,selected_order=order,log_marginal_calls=calls,
        material_statistic_row_evaluations=len(f)*calls,
        complex_channel_statistic_evaluations=len(f)*calls*240,
        material_gamma_likelihood_evaluations=node_evaluations,
        unique_response_rows=len(f),physical_evaluations=0)


def infer_one(cache,rules,cuts,y,reference):
    ll,nuisance = marginal_shared_gamma(cache['F'],y)
    order=cache['order']
    name=f'full_{order}'
    full=posterior(rules[name],ll[cache['maps'][name]])
    cdf=[]
    for i in range(5):
        name=f'cdf_{i}_{order}'
        r=rules[name]
        numerator=float(logsumexp(np.log(r['weights'])+ll[cache['maps'][name]]))
        cdf.append(float(np.exp(numerator-full['logZ'])))
    refcdf=[float(reference['CDF'][str(cut)]['order9']) for cut in cuts]
    branch_error=abs(full['branch1']-float(reference['branch1_posterior_mass']))
    cdf_errors=np.abs(np.asarray(cdf)-refcdf)
    coherent=cdf_coherence(cdf)
    reference_ok=(reference['reference_status']=='QUALIFIED_EMPIRICAL_QUADRATURE'
        and reference['reference_CDF_coherence'] is True
        and float(reference['quadrature_event_error'])<=PROBABILITY_TOLERANCE
        and float(reference['quadrature_CDF_max_error'])<=PROBABILITY_TOLERANCE
        and cdf_coherence(refcdf))
    reasons=[]
    if not reference_ok:reasons.append('REFERENCE_UNQUALIFIED')
    if not nuisance['converged']:reasons.append('NUISANCE_CONVERGENCE_FAILED')
    if not coherent:reasons.append('CDF_OUT_OF_BOUNDS_OR_NONMONOTONE')
    if branch_error>PROBABILITY_TOLERANCE:reasons.append('BRANCH_EVENT_ERROR_GT_0.02')
    if float(np.max(cdf_errors))>PROBABILITY_TOLERANCE:reasons.append('CDF_ERROR_GT_0.02')
    prefix='direct_order3' if cache['method']==METHODS[0] else 'compressed'
    saved_event=reference.get(prefix+'_branch_error')
    saved_cdf=reference.get(prefix+'_CDF_max_error')
    saved_replay_difference={
        'branch_error_difference':None if saved_event is None else abs(branch_error-float(saved_event)),
        'CDF_max_error_difference':None if saved_cdf is None else abs(float(np.max(cdf_errors))-float(saved_cdf))}
    return dict(observation=int(reference['observation']),parent=int(reference['parent']),
        status='QUALIFIED_ON_FINITE_REPLAY' if not reasons else 'UNQUALIFIED',
        failure_reasons=reasons,reference_qualified=reference_ok,
        branch1_mass=full['branch1'],branch_event_error=branch_error,
        x_mean=full['x_mean'],x_sd=full['x_sd'],contrast_mean=full['contrast_mean'],
        x_mean_error=abs(full['x_mean']-float(reference['x_mean'])),
        x_SD_error=abs(full['x_sd']-float(reference['x_sd'])),
        logZ=full['logZ'],logZ_absolute_error=abs(full['logZ']-float(reference['log_evidence'])),
        CDF={str(cut):dict(value=value,reference_order9=rv,error=float(err))
             for cut,value,rv,err in zip(cuts,cdf,refcdf,cdf_errors)},
        CDF_max_error=float(np.max(cdf_errors)),CDF_coherent=coherent,
        saved_reference_method_error_consistency=saved_replay_difference,nuisance=nuisance)


def replay(cache,rules,cuts,ys,references,observation_order):
    """Always recompute likelihood and posterior for all y, even in warm passes."""
    if sorted(map(int,observation_order))!=list(range(len(ys))):
        raise ValueError('Every complete pass must visit each observation exactly once')
    rows=[]
    for i in observation_order:
        answer,timing=timed(infer_one,cache,rules,cuts,ys[int(i)],references[int(i)])
        answer['timing']=timing
        rows.append(answer)
    return rows


def summarize(rows):
    counts={k:sum(r['nuisance'][k] for r in rows) for k in
        ('log_marginal_calls','material_statistic_row_evaluations',
         'complex_channel_statistic_evaluations','material_gamma_likelihood_evaluations')}
    return dict(observations=len(rows),qualified_count=sum(not r['failure_reasons'] for r in rows),
        unqualified=[dict(observation=r['observation'],parent=r['parent'],reasons=r['failure_reasons'])
                     for r in rows if r['failure_reasons']],
        max_branch_event_error=max(r['branch_event_error'] for r in rows),
        max_CDF_error=max(r['CDF_max_error'] for r in rows),
        max_x_mean_error=max(r['x_mean_error'] for r in rows),
        max_x_SD_error=max(r['x_SD_error'] for r in rows),
        max_logZ_absolute_error=max(r['logZ_absolute_error'] for r in rows),
        all_CDF_coherent=all(r['CDF_coherent'] for r in rows),
        all_nuisance_converged=all(r['nuisance']['converged'] for r in rows),
        mathematical_evaluation_counts=counts,new_Maxwell_rhs=0,new_Maxwell_LU=0)


def read_cost_ledger(root=ROOT,stage='P1A'):
    """Read-only accounting snapshot, no Book; never retain payload/coordinates/GT."""
    path=Path(root)/'PHYSICS_LEDGER.jsonl'
    total=dict(rhs=0,factorizations=0,cpu_seconds=0.)
    jobs=[]
    scopes={}
    with path.open() as stream:
        for line in stream:
            if not line.strip():continue
            row=json.loads(line)
            if row.get('stage')!=stage:continue
            if row.get('event')=='reserve':
                total['rhs']+=int(row.get('rhs',0))
                total['factorizations']+=int(row.get('factorizations',0))
                if row.get('scope') is not None:
                    bucket=scopes.setdefault(str(row['scope']),dict(rhs=0,factorizations=0))
                    bucket['rhs']+=int(row.get('rhs',0))
                    bucket['factorizations']+=int(row.get('factorizations',0))
            total['cpu_seconds']+=float(row.get('charged_cpu_seconds',0))
            if row.get('event')=='job_cost':
                jobs.append({k:row.get(k) for k in ('utc','label','charged_cpu_seconds')})
    return dict(source=f'PHYSICS_LEDGER.jsonl: {stage} cost fields only',
        **{f'{stage}_totals_at_read':total},job_CPU_charges=jobs,scope_reservations=scopes,
        includes_current_replay_launcher_final_charge=False,
        note='Snapshot precedes current run_job completion; launcher charges CPU after main returns. No ledger mutation here.')


def p1b_additional_setup(root=ROOT):
    """Report original four-arm bank cost without loading its 650 F vectors."""
    complete=_json(Path(root)/'results/P1B/full_material_rule/COMPLETE.json')
    if (complete.get('status')!='FULL_MATERIAL_RULE_COMPLETE' or complete.get('nodes')!=650
            or complete.get('material_rule')!=[2,65,5]):
        raise ValueError('Missing/mismatched original P1B 650-node bank cost record')
    ledger=read_cost_ledger(root,'P1B')
    scope=complete['scope']
    scope_cost=ledger['scope_reservations'].get(scope)
    if scope_cost is None:
        raise ValueError('No actual ledger scope cost for the P1B full bank; cannot call it free')
    if int(scope_cost['rhs'])!=int(complete['scope_used_rhs']):
        raise ValueError('P1B full-bank ledger and completion RHS costs disagree')
    return dict(source='results/P1B/full_material_rule/COMPLETE.json and cost-only P1B ledger',
        material_nodes=650,material_rule=[2,65,5],hypothetical_all_unpaid_rhs=13000,
        hypothetical_all_unpaid_LU=3250,actual_additional_bank_scope_cost=scope_cost,
        bank_scope=scope,stage_totals_snapshot=ledger,
        required_for_this_direct3_or_interpolation18_replay=False,
        F_vectors_loaded_for_this_replay=False,
        attribution='Additional setup for original P1B four-arm full-node experiment, not a required direct3/interpolation18 bank cost. Retained in total study/four-arm cost; never treated as free or as an avoided solve.',
        overlap_note='Stage totals include this scope. Do not add scope to stage total again; P1A shared-cache reuse can make actual B charges smaller than all-unpaid 650-vector arithmetic.')


def adapt_p1b_reference(row):
    """Fail closed on absent B reference; preserve explicit unqualified status."""
    c=row.get('continuous_ref')
    if not isinstance(c,dict):
        raise ValueError('P1B continuous_ref missing; replay cannot qualify without reference')
    required=('qualified','status','branch_mass','x_mean','x_sd','logZ','CDF',
              'CDF_coherent','event_order_error','CDF_order_max_error')
    missing=[k for k in required if k not in c]
    if missing:
        raise ValueError('Incomplete P1B continuous_ref, no qualification inferred: '+','.join(missing))
    masses=np.asarray(c['branch_mass'],float)
    metrics=np.asarray([c['x_mean'],c['x_sd'],c['logZ'],c['event_order_error'],c['CDF_order_max_error']],float)
    if (masses.shape!=(2,) or not np.isfinite(masses).all() or np.any(masses<0)
            or abs(float(masses.sum())-1)>1e-10 or not np.isfinite(metrics).all()
            or float(c['x_sd'])<0 or np.any(metrics[-2:]<0)):
        raise ValueError('Invalid numeric P1B continuous reference')
    return dict(observation=int(row['observation']),parent=int(row['parent']),
        reference_status=c['status'] if c['qualified'] is True else 'REFERENCE_UNQUALIFIED',
        branch1_posterior_mass=float(masses[1]),x_mean=float(c['x_mean']),x_sd=float(c['x_sd']),
        contrast_mean=c.get('contrast_mean'),log_evidence=float(c['logZ']),CDF=c['CDF'],
        reference_CDF_coherence=c['CDF_coherent'],quadrature_event_error=float(c['event_order_error']),
        quadrature_CDF_max_error=float(c['CDF_order_max_error']),
        reference_reported_qualified=c['qualified'],reference_reported_status=c['status'],
        source='results/P1B/inference/RESULTS.json continuous_ref only',
        reference_contrast_scored=False)


def load_observations_references(dataset,root=ROOT):
    if dataset not in ('P1A','P1B'):
        raise ValueError('dataset must be P1A or P1B')
    base=Path(root)/'results'/dataset
    with np.load(base/'data/OBSERVATIONS.npz',allow_pickle=False) as z:
        ys=z['y'].copy()
        parent_ids=z['parent_ids'].copy()
    if dataset=='P1A':
        references=_json(base/'reference/RESULTS.json')
    else:
        # The container includes other methods; only identities + continuous_ref
        # are passed downstream. No fitted method output or truth enters replay.
        references=[adapt_p1b_reference(row) for row in _json(base/'inference/RESULTS.json')]
    if ys.shape!=(N_OBSERVATIONS,5,4,12) or not np.iscomplexobj(ys) or not np.isfinite(ys).all():
        raise ValueError(f'Expected exactly32 finite complex {dataset} observations')
    if parent_ids.shape!=(N_OBSERVATIONS,) or len(references)!=N_OBSERVATIONS:
        raise ValueError('Reference/observation identity mismatch')
    for i,row in enumerate(references):
        if int(row['observation'])!=i or int(row['parent'])!=int(parent_ids[i]):
            raise ValueError('Reference observation/parent ordering differs from y')
    return ys,parent_ids,references


def physical_costs(rules,completion,caches,root=ROOT):
    teacher=unique_node_registry(rules,list(rules))
    reference_names=[f'full_{n}' for n in (5,9)]+[f'cdf_{i}_{n}' for i in range(5) for n in (5,9)]
    reference_only=unique_node_registry(rules,reference_names)
    expected=len(teacher)
    if expected!=678 or int(completion['unique_paid_points'])!=expected:
        raise ValueError('Public setup does not match the completed 678-node bank')
    base=dict(public_teacher_and_all_rules_unique_nodes=expected,
        public_teacher_equivalent_rhs=20*expected,public_teacher_equivalent_LU=5*expected,
        reference_orders5_9_unique_nodes=len(reference_only),
        reference_orders5_9_equivalent_rhs=20*len(reference_only),
        reference_orders5_9_equivalent_LU=5*len(reference_only),
        actual_rules_completion_P1A_totals=completion['totals']['P1A'],
        ledger_snapshot=read_cost_ledger(root),
        physical_vector_contract='5 frequencies x (one LU + four source RHS): 5 LU / 20 RHS per complete material vector',
        setup_treatment='All 678 public teacher/validation rules are already paid common setup, never free. Reference-only and method node sets overlap; their costs must not be added as disjoint charges.',
        current_replay_incremental_Maxwell_rhs=0,current_replay_incremental_Maxwell_LU=0)
    actual=completion['totals']['P1A']
    base['P1A_nonbank_residual_at_completion']={
        'rhs':int(actual['rhs'])-20*expected,'factorizations':int(actual['factorizations'])-5*expected,
        'attribution':'Stage total minus 678-node public bank; no hidden endpoint/GT coordinates read or inferred. This residual is separate setup/data-generation accounting.'}
    methods={}
    for method,cache in caches.items():
        keys=cache['unique_physical_training_keys']
        if not set(keys).issubset(teacher):raise ValueError('Method requires a nonpublic physical node')
        methods[method]=dict(minimum_unique_physical_nodes_for_defined_method=len(keys),
            node_keys=keys,hypothetical_unpaid_bank_rhs=20*len(keys),hypothetical_unpaid_bank_LU=5*len(keys),
            actual_access_incremental_rhs=0,actual_access_incremental_LU=0,
            prospective_deployment_qualification='NOT_GUARANTEED',
            qualification_setup='The common 678-node teacher/reference bank was paid to establish the measured qualification on these observations.',
            warning='Minimal deployment-bank arithmetic is conditional on this public family and is not a certified future observation/uniform-likelihood guarantee.')
    base['methods']=methods
    base['comparison_verdict']='NO_AUTOMATIC_GO; array-runtime differences are not Full-solve savings'
    return base,teacher


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',choices=('P1A','P1B'),default='P1A',
        help='Use unchanged P1A response rules on P1A or new P1B observations')
    return parser.parse_args(argv)


def main(argv=None):
    dataset=parse_args(argv).dataset
    out=ROOT/f'results/{dataset}/cost_comparison'
    (rules,cuts,completion),input_timing=timed(load_rules,ROOT)
    # Construct both caches before reading y. No fitting or selection on observations.
    caches={}
    build_times={}
    for method in METHODS:
        caches[method],build_times[method]=timed(build_response_cache,method,rules)
    if [len(caches[m]['unique_physical_training_keys']) for m in METHODS]!=[54,18]:
        raise ValueError('The frozen54-node direct3 /18-node interpolation methods changed')
    (ys,parent_ids,references),observation_IO_timing=timed(load_observations_references,dataset,ROOT)
    for row in references:
        if set(row['CDF'])!=set(map(str,cuts)):
            raise ValueError('Saved reference CDF checkpoints differ from paid rules')
    costs,teacher=physical_costs(rules,completion,caches,ROOT)
    if dataset=='P1B':
        costs['P1B_original_four_arm_additional_setup']=p1b_additional_setup(ROOT)
    out.mkdir(parents=True,exist_ok=True)
    reference_source=('results/P1A/reference/RESULTS.json' if dataset=='P1A'
        else 'results/P1B/inference/RESULTS.json continuous_ref only')
    contract=dict(identity=('P1A_SHARED_F_CACHE_ARRAY_COST_REPLAY_V1' if dataset=='P1A' else
        'P1B_SHARED_P1A_F_CACHE_ARRAY_COST_REPLAY_V1'),dataset=dataset,
        target='CONTROLLED_W_FIXED_SUPPORT_JOINT_V1_CONTINUOUS_APPROXIMATED_BY_POSITIVE_QUADRATURE',
        reference=reference_source+'; full5/full9 and directly truncated CDF5/CDF9',
        observations=f'{dataset} OBSERVATIONS.npz: y and parent_ids only; all32, fixed before comparison',
        observation_scope=('Original32 P1A replay' if dataset=='P1A' else
            'New32 P1B observations relative to P1A setup; observation-only generalization check of unchanged methods; no B adaptation or model/threshold/rule selection'),
        response_banks='Only completed P1A full3/5/9 and cdf_i_3/5/9; B650 responses are not loaded',
        fixed_minimum_method_nodes=dict(direct_order3=54,response_interpolation18=18),
        method_and_threshold_changes_on_P1B=False,reference_missing_policy='FAIL_CLOSED; never infer qualification',
        reference_contrast_scored=False,
        thresholds=dict(branch_event=.02,CDF=.02,nuisance_max_node_logZ_change=NUISANCE_TOLERANCE),
        gain='ONE shared amplitude U[.99,1.01] and phase U[-.01,.01] across all frequencies/sources/receivers, marginalized by posterior_math; no profiling',
        CDF='Integrate original-density truncated intervals, divide by same method fullZ; never atomic node-count CDF or independently renormalized truncation',
        nuisance_orders=list(NUISANCE_ORDERS),
        cache_rights='Both reuse exact unique F arrays and rule-node index maps. Neither retains per-y likelihood or posterior across passes. Shared teacher/reference is already paid.',
        cache_construction='Once per method, prior-rule coordinates only, before loading y',
        cold_definition='First in-process cache construction plus first full32 inference; not a fresh process or cold OS disk cache. Shared rule loading reported separately.',
        warm_definition='Two complete recomputations of all32 y from persistent F cache; no result memoization',
        method_order_by_pass={'cold':list(METHODS),'warm1':list(reversed(METHODS)),'warm2':list(METHODS)},
        warm2_order_seed=WARM_ORDER_SEED,rule_input_load_timing=input_timing,
        observation_reference_IO_timing=observation_IO_timing,
        reference_setup_free=False,uniform_likelihood_bound=None,
        P1B_original_four_arm_extra_bank_cost=('Separately retained in PHYSICAL_SETUP_COSTS.json, not required by either replay method' if dataset=='P1B' else 'NOT_LOADED_FOR_P1A'),
        original_minimum_feature_semantics='UNRESOLVED; no qualification of the original broad wall prior',
        model_scope='Only the public controlled fixed-support two-wall family; no general unknown-wall online claim',
        new_Maxwell_rhs=0,new_Maxwell_LU=0,GT_access=False,plan_access=False,
        decision='NOT_EVALUATED_NO_AUTOMATIC_GO')
    _write(out/'CONTRACT.json',contract)
    _write(out/'PHYSICAL_SETUP_COSTS.json',costs)
    _write(out/'ARRAY_CACHE_BUILD.json',{m:dict(**build_times[m],
        unique_response_rows=len(caches[m]['F']),response_cache_bytes=caches[m]['memory_bytes'],
        interpolation_calls=caches[m]['interpolation_calls'],
        interpolation_material_rows=caches[m]['interpolation_material_rows'],
        response_error_white_on_evaluated_paid_nodes=caches[m]['response_error_white_on_evaluated_paid_nodes'],
        error_scope='Finite evaluated full9/CDF9 nodes only, not a uniform prior-support or likelihood error bound') for m in METHODS})
    permutation=np.random.default_rng(WARM_ORDER_SEED).permutation(N_OBSERVATIONS).tolist()
    passes=[('cold',list(range(N_OBSERVATIONS)),METHODS),
            ('warm1',list(range(N_OBSERVATIONS)),tuple(reversed(METHODS))),
            ('warm2',permutation,METHODS)]
    results={}
    flat=[]
    baseline={}
    for pass_name,observation_order,method_order in passes:
        results[pass_name]={}
        for method in method_order:
            rows,timing=timed(replay,caches[method],rules,cuts,ys,references,observation_order)
            by_id={r['observation']:r for r in rows}
            if method in baseline:
                for i,r in by_id.items():
                    if not np.allclose([r['branch1_mass'],r['x_mean'],r['x_sd'],r['logZ']],
                        [baseline[method][i][k] for k in ('branch1_mass','x_mean','x_sd','logZ')],rtol=0,atol=1e-10):
                        raise ValueError('Warm replay changed deterministic posterior arithmetic')
            else:baseline[method]=by_id
            summary=summarize(rows)
            total=dict(timing)
            if pass_name=='cold':
                total={k:timing[k]+build_times[method][k] for k in timing}
            record=dict(method=method,pass_name=pass_name,observation_order=observation_order,
                inference_timing=timing,total_timing_including_first_build=total,
                summary=summary,rows=rows)
            results[pass_name][method]=record
            _write(out/f'{pass_name}_{method}.json',record)
            for r in rows:
                flat.append(dict(pass_name=pass_name,method=method,observation=r['observation'],parent=r['parent'],
                    status=r['status'],failure_reasons=';'.join(r['failure_reasons']),
                    branch1_mass=r['branch1_mass'],branch_event_error=r['branch_event_error'],CDF_max_error=r['CDF_max_error'],
                    CDF_coherent=r['CDF_coherent'],x_mean=r['x_mean'],x_sd=r['x_sd'],x_mean_error=r['x_mean_error'],
                    logZ=r['logZ'],logZ_absolute_error=r['logZ_absolute_error'],
                    wall_seconds=r['timing']['wall_seconds'],cpu_seconds=r['timing']['cpu_seconds'],
                    nuisance_order=r['nuisance']['selected_order'],nuisance_converged=r['nuisance']['converged'],
                    log_marginal_calls=r['nuisance']['log_marginal_calls'],
                    material_gamma_likelihood_evaluations=r['nuisance']['material_gamma_likelihood_evaluations']))
            _csv(out/'PER_OBSERVATION.csv',flat)
            print(json.dumps(dict(pass_name=pass_name,method=method,**timing,
                qualified=summary['qualified_count'],max_CDF_error=summary['max_CDF_error'])),flush=True)
    comparison=[]
    for pass_name,_,_ in passes:
        direct=results[pass_name][METHODS[0]]
        comp=results[pass_name][METHODS[1]]
        dt=direct['total_timing_including_first_build']
        ct=comp['total_timing_including_first_build']
        comparison.append(dict(pass_name=pass_name,direct_wall_seconds=dt['wall_seconds'],
            interpolation_wall_seconds=ct['wall_seconds'],direct_cpu_seconds=dt['cpu_seconds'],
            interpolation_cpu_seconds=ct['cpu_seconds'],
            direct_over_interpolation_wall_ratio=dt['wall_seconds']/ct['wall_seconds'],
            direct_over_interpolation_cpu_ratio=dt['cpu_seconds']/ct['cpu_seconds'],
            direct_qualified_count=direct['summary']['qualified_count'],
            interpolation_qualified_count=comp['summary']['qualified_count'],
            meaning='Array computation only; not avoided physical solves or an automatic GO'))
    _csv(out/'PASS_COMPARISON.csv',comparison)
    _csv(out/'PHYSICAL_NODES.csv',[dict(key=k,branch=r['branch'],u=r['u'],contrast=r['contrast'],
        rule_membership=';'.join(r['rules']),
        direct_required=k in set(caches[METHODS[0]]['unique_physical_training_keys']),
        interpolation_required=k in set(caches[METHODS[1]]['unique_physical_training_keys'])) for k,r in teacher.items()])
    _write(out/'RESULTS.json',dict(contract=contract,passes=comparison,
        method_accuracy={m:results['cold'][m]['summary'] for m in METHODS},
        physical_costs='PHYSICAL_SETUP_COSTS.json',array_cache_build='ARRAY_CACHE_BUILD.json',
        empirical_scope=f'Measured finite set of32 {dataset} observations with frozen P1A response banks and checked quadrature nodes, not continuous-exact posterior, uniform likelihood bound or proof of future deployment savings.',
        shared_reference_and_teacher_free=False,uniform_likelihood_bound=None,
        new_Maxwell_rhs=0,new_Maxwell_LU=0,scientific_GO='NOT_EVALUATED'))
    _write(out/'COMPLETE.json',dict(status='ARRAY_REPLAY_COMPLETE_NOT_A_CAMPAIGN_GO',dataset=dataset,
        observations=32,passes=3,methods=list(METHODS),new_Maxwell_rhs=0,new_Maxwell_LU=0,
        CPU_charging=f'Parent run_job {dataset}; this module does not instantiate or write Book'))


if __name__=='__main__':
    main()
