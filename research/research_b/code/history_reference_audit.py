"""Fixed first-four-parent joint-history audit, zero additional Maxwell.

The original coefficient prior and one proper persistent truncated gain are
integrated once over each whole selected observation history. This is a finite
proxy-target numerical audit, not a full-Maxwell posterior certification.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import os
import resource
import sys
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from scipy.special import logsumexp

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(REPO))
from src.b1_joint.material import RMAX, geometry, event_ids, from_geometry
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.likelihood import log_likelihood
from active_experiment import Simulator, FREQS, IDS
from e0_candidate_controls import profile
from sensing import normalized, resample_indices, conditional_log_increment, field_l1
from cost import Cost, ROOT, dump

PARENTS = (0,1,2,3)
NOISE = (0,1)
POLICIES = ('fixed_3.2','expected_material_risk')
ORDERS = (13,17)
WIDTHS = np.array([*RMAX,.6,.12])
CHUNK = 8192


def _peak_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == 'darwin' else value*1024)


def _memory_check():
    if _peak_bytes() > 2*1024**3:
        raise MemoryError('Fixed joint-history audit exceeded approximately 2 GiB peak RSS')


def _history(simulator, source_rows, sigma, calibration_sd):
    parent, noise = source_rows[0]['parent'], source_rows[0]['noise']
    blocks = [dict(action='repeat_1.6',frequency=1.6,
                   y=simulator.initial(parent,noise),sigma=np.full(32,sigma[FREQS.index(1.6)]))]
    for row in source_rows[1:]:
        aid = row['selected_action']
        if aid == 'reference_calibration':
            frequency, sd = None, np.array([calibration_sd])
        elif aid in IDS:
            frequency = FREQS[IDS.index(aid)]
            sd = np.full(32,sigma[IDS.index(aid)])
        else:
            raise ValueError('This fixed audit supports selected frequency/calibration histories only')
        blocks.append(dict(action=aid,frequency=frequency,
            y=simulator.measure(parent,noise,row['step']-1,aid),sigma=sd))
    if [block['action'] for block in blocks] != source_rows[-1]['history']:
        raise ValueError('Selected action sequence differs from saved final history')
    return blocks


def _joint_fields(fields, blocks, start=None, end=None):
    count = len(next(iter(fields.values()))) if start is None else end-start
    out = []
    for block in blocks:
        if block['frequency'] is None:
            out.append(np.ones((count,1)))
        else:
            value = fields[block['frequency']]
            out.append(value if start is None else value[start:end])
    return np.concatenate(out,axis=1)


def _joint_ll(fields, blocks, cost):
    """One marginal over the entire history, chunked only for memory."""
    n = len(next(iter(fields.values())))
    output = np.empty(n)
    y = np.concatenate([block['y'] for block in blocks])
    sigma = np.concatenate([block['sigma'] for block in blocks])
    for start in range(0,n,CHUNK):
        end = min(n,start+CHUNK)
        cost.check_time()
        output[start:end] = log_likelihood(_joint_fields(fields,blocks,start,end),y,sigma)
    return output


def _summarize(x, physical, events, weights, log_evidence, selected):
    mean = weights@physical
    variance = weights@((physical-mean)**2)
    risks = {name:float(field_l1(latent,x)[0]@weights) for name,latent in selected.items()}
    _memory_check()
    return dict(event_mass=np.bincount(events,weights=weights,minlength=12),
        mean=mean, SD=np.sqrt(np.maximum(variance,0.)),
        log_evidence=float(log_evidence),selected_candidate_expected_L1=risks,
        effective_weights=float(1/np.sum(weights*weights)))


def _prepare_bank(response, order, needed_frequencies, cost):
    cost.check_time()
    reference = Reference(response)
    source = reference.bank(order)
    # Keep exactly the returned x/prior/physical/event arrays. Only discard
    # redundant cell copies and unused stored frequency arrays from this local
    # audit instance; the Reference algorithm/source and disk assets are intact.
    bank = {key:source[key] for key in ('x','physical','prior','events')}
    fields = {1.6:source['F1']}
    if 2.4 in needed_frequencies:
        fields[2.4] = source['F2']
    reference.cache.clear()
    del source,reference
    gc.collect()
    for frequency in sorted(needed_frequencies):
        if frequency not in fields:
            cost.check_time()
            fields[frequency] = response.predict(bank['x'],[frequency])
        _memory_check()
    return bank,fields


def _finite_replay(response, bank, initial_mass, initial_logz, blocks,
                   candidate_x, candidate_fields, source_rows):
    parent, noise = source_rows[0]['parent'],source_rows[0]['noise']
    ids = resample_indices(initial_mass,2048,26101530+parent*17+noise)
    x = bank['x'][ids]
    frequencies = sorted({block['frequency'] for block in blocks if block['frequency'] is not None})
    fields = {frequency:response.predict(x,[frequency]) for frequency in frequencies}
    first = blocks[0]
    histF, y, sd = fields[1.6].copy(),first['y'].copy(),first['sigma'].copy()
    weights = np.full(2048,1/2048)
    decisions = bank['x'][resample_indices(initial_mass,96,26101531)].copy()
    for step,block in enumerate(blocks[1:],start=1):
        future = np.ones((len(x),1)) if block['frequency'] is None else fields[block['frequency']]
        increment = conditional_log_increment(histF,y,sd,future,block['y'],block['sigma'])[0]
        weights = normalized(np.log(np.maximum(weights,1e-300))+increment)
        histF = np.column_stack((histF,future))
        y,sd = np.r_[y,block['y']],np.r_[sd,block['sigma']]
        if source_rows[0]['audit_stage'] == 'E2':
            decisions = np.concatenate((decisions,x[resample_indices(weights,64,26101580+step)]))
    candidate_residuals,_ = profile(_joint_fields(candidate_fields,blocks),y,sd)
    ci = int(np.argmin(candidate_residuals))
    reference_losses = field_l1(decisions,x)@weights
    ri = int(np.argmin(reference_losses))
    selected = dict(mode1_profile_selected=candidate_x[ci],mode2_reference_selected=decisions[ri])
    joint = log_likelihood(histF,y,sd)
    old = log_likelihood(fields[1.6],first['y'],first['sigma'])
    once_weights = normalized(joint-old)
    approximate_logz = initial_logz+logsumexp(joint-old)-np.log(len(x))
    physical = np.column_stack((geometry(x),x[:,2:]))
    summary = _summarize(x,physical,event_ids(x),weights,approximate_logz,selected)
    last = source_rows[-1]
    return dict(summary=summary, selected_latent=selected,
        sample_count=2048, retained_reference_catalogue=len(decisions), candidate_index=ci,
        saved_ESS=last['posterior_ESS'], ESS_minus_saved=float(summary['effective_weights']-last['posterior_ESS']),
        replay_reference_risk_minus_saved=float(summary['selected_candidate_expected_L1']['mode2_reference_selected']-last['finite_reference_risk']),
        profile_chi2_minus_saved=float(candidate_residuals[ci]-last['candidate_profile_chi2']),
        maximum_sequential_vs_joint_once_weight_difference=float(np.max(abs(weights-once_weights))),
        approximation='original initial weighted-rule stratified resampling, then persistent-gain conditional updates; no continuous proposal-density or iid claim',
        log_evidence_scope='initial13 normalization times finite conditional likelihood-ratio mean; approximate')


def run(cost):
    results = ROOT/'results'
    final_path = results/'HISTORY_REFERENCE_AUDIT.json'
    if final_path.exists():
        raise ValueError('Existing history audit preserved; no implicit repeat')
    response = Response(results/'E1_RESPONSE_BANK.npz')
    simulator = Simulator(results/'E1_FRESH_SIMULATOR.npz')
    anchors = np.load(results/'E1_NOISE_ANCHORS.npz',allow_pickle=False)
    sigma = anchors['sigma']
    calibration_sd = sigma[FREQS.index(1.6)]/np.linalg.norm(anchors['F'][FREQS.index(1.6)])
    labels,unique = [],{}
    for stage in ('E1','E2'):
        rows = json.loads((results/f'{stage}_POLICY_RESULTS.json').read_text())
        for parent in PARENTS:
            for noise in NOISE:
                for policy in POLICIES:
                    source = sorted([dict(row,audit_stage=stage) for row in rows if row['parent']==parent
                        and row['noise']==noise and row['policy']==policy],key=lambda row:row['step'])
                    if not source or source[0]['step']!=0:
                        raise ValueError('Fixed audit input episode unavailable')
                    blocks = _history(simulator,source,sigma,calibration_sd)
                    key = f'{stage}/parent{parent}/noise{noise}/'+','.join(block['action'] for block in blocks)
                    unique.setdefault(key,dict(blocks=blocks,source_rows=source,parent=parent,noise=noise,stage=stage))
                    labels.append(dict(stage=stage,parent=parent,noise=noise,policy=policy,key=key,
                        final_step=source[-1]['step'],saved_ESS=source[-1]['posterior_ESS']))
    needed = {block['frequency'] for row in unique.values() for block in row['blocks'] if block['frequency'] is not None}
    raw = np.load(results/'E1_RESPONSE_BANK.npz',allow_pickle=False)
    mesh = np.stack([v.ravel() for v in np.meshgrid(*[raw[f'axis_{j}'] for j in range(4)],indexing='ij')],axis=1)
    candidate_x = np.array([from_geometry(*row) for row in mesh])
    candidate_fields = {frequency:raw[f'fields_{frequency:g}'].reshape(-1,32) for frequency in needed}
    summaries,finite,bank_info = {},{},{}
    for order in ORDERS:
        bank,fields = _prepare_bank(response,order,needed,cost)
        logprior = np.log(bank['prior'])
        bank_info[order] = dict(nodes=len(bank['x']),prior_total_mass=float(bank['prior'].sum()),
            prior_event_mass=np.bincount(bank['events'],weights=bank['prior'],minlength=12),
            frequencies=sorted(fields),measure='original coefficient-uniform mixed pushforward with zero-geometry atoms; actual quadrature prior weights',
            nodes_are_not_an_equal_weight_material_prior=True)
        initial_cache = {}
        for key,row in unique.items():
            cost.check_time()
            blocks = row['blocks']
            if order == 13:
                ikey = (row['parent'],row['noise'])
                if ikey not in initial_cache:
                    initial_ll = _joint_ll(fields,blocks[:1],cost)
                    initial_logz = float(logsumexp(initial_ll+logprior))
                    initial_cache[ikey] = (np.exp(initial_ll+logprior-initial_logz),initial_logz)
                initial_mass,initial_logz = initial_cache[ikey]
                finite[key] = _finite_replay(response,bank,initial_mass,initial_logz,blocks,
                    candidate_x,candidate_fields,row['source_rows'])
            ll = _joint_ll(fields,blocks,cost)
            logz = float(logsumexp(ll+logprior))
            mass = np.exp(ll+logprior-logz)
            selected = finite[key]['selected_latent']
            summaries[(key,order)] = _summarize(bank['x'],bank['physical'],bank['events'],mass,logz,selected)
            print(json.dumps(dict(order=order,history=key,peak_RSS_bytes=_peak_bytes(),new_Maxwell=0)),flush=True)
        dump(results/'HISTORY_REFERENCE_AUDIT_PARTIAL.json',dict(completed_order=order,
            unique_histories=len(unique),bank_info=bank_info,
            scope='partial numerical audit; no changes to original algorithm or statistical outcomes'))
        del bank,fields,logprior,initial_cache,ll,mass
        gc.collect()
    tolerances = json.loads((REPO/'research/joint_p2/configs/protocol.json').read_text())['reference_tolerances']
    evaluated,flat = {},[]
    for key in unique:
        low,high = summaries[(key,13)],summaries[(key,17)]
        sample = finite[key]['summary']
        refinement = dict(event_mass_max=float(np.max(abs(low['event_mass']-high['event_mass']))),
            mean_fraction_prior_width=float(np.max(abs(low['mean']-high['mean'])/WIDTHS)),
            SD_relative=float(np.max(abs(low['SD']-high['SD'])/np.maximum(high['SD'],WIDTHS*1e-6))),
            log_evidence_abs=float(abs(low['log_evidence']-high['log_evidence'])),
            selected_risk_abs={name:float(abs(low['selected_candidate_expected_L1'][name]-high['selected_candidate_expected_L1'][name]))
                for name in high['selected_candidate_expected_L1']})
        stable = (refinement['event_mass_max']<=tolerances['event_mass'] and
            refinement['mean_fraction_prior_width']<=tolerances['mean_quantile_fraction_prior_width'] and
            refinement['SD_relative']<=tolerances['SD_relative'])
        error = dict(event_mass_max=float(np.max(abs(sample['event_mass']-high['event_mass']))),
            mean_fraction_prior_width=float(np.max(abs(sample['mean']-high['mean'])/WIDTHS)),
            mean_abs=np.abs(sample['mean']-high['mean']),
            selected_risk_abs={name:float(abs(sample['selected_candidate_expected_L1'][name]-high['selected_candidate_expected_L1'][name]))
                for name in high['selected_candidate_expected_L1']})
        evaluated[key] = dict(status='CONSISTENT_13_17_ON_CHECKED_SUMMARIES' if stable else 'UNRESOLVED',
            order13=low,order17=high,finite2048=finite[key],refinement13_17=refinement,
            finite2048_minus_order17=error,selected_observation_history=unique[key]['blocks'])
    for label in labels:
        item = evaluated[label['key']]
        a,b = item['refinement13_17'],item['finite2048_minus_order17']
        flat.append(dict(**label,refinement_status=item['status'],event_refinement13_17=a['event_mass_max'],
            mean_fraction_refinement13_17=a['mean_fraction_prior_width'],SD_relative_refinement13_17=a['SD_relative'],
            log_evidence_refinement13_17=a['log_evidence_abs'],
            finite_event_mass_error_vs17=b['event_mass_max'],finite_mean_fraction_error_vs17=b['mean_fraction_prior_width'],
            candidate_risk13=item['order13']['selected_candidate_expected_L1']['mode1_profile_selected'],
            candidate_risk17=item['order17']['selected_candidate_expected_L1']['mode1_profile_selected'],
            finite_candidate_risk=item['finite2048']['summary']['selected_candidate_expected_L1']['mode1_profile_selected'],
            finite_candidate_risk_error_vs17=b['selected_risk_abs']['mode1_profile_selected'],
            finite_reference_selected_risk_error_vs17=b['selected_risk_abs']['mode2_reference_selected'],
            finite_ESS=item['finite2048']['summary']['effective_weights']))
    maximum = lambda path:max(float(path(row)) for row in evaluated.values())
    report = dict(schema='research-b-fixed-joint-history-audit-v1',
        fixed_parents=PARENTS,noise_repetitions=NOISE,policies=POLICIES,orders=ORDERS,
        policy_rows=len(labels),unique_histories=len(unique),independent_parent_count=4,
        physical_mean_coordinates=['central_radius','annular_halfwidth','epsilon_real','epsilon_imag'],
        target='same original coefficient prior, continuous saved-response proxy, one proper persistent truncated N(1,.15^2) gain on[.7,1.3]',
        likelihood='joint old plus actually selected observations once; no multiplication of separate gain-marginal blocks',
        future_truth_read=False,unselected_observations_used=False,new_Maxwell=0,GPU=0,
        bank_info=bank_info,refinement_tolerances_from_prior_P2_protocol=tolerances,
        risk_refinement_tolerance='not newly thresholded; actual absolute differences reported',
        unresolved_unique_histories=[key for key,row in evaluated.items() if row['status']=='UNRESOLVED'],
        maximum_finite_event_mass_error_vs17=maximum(lambda row:row['finite2048_minus_order17']['event_mass_max']),
        maximum_finite_mean_fraction_error_vs17=maximum(lambda row:row['finite2048_minus_order17']['mean_fraction_prior_width']),
        maximum_finite_selected_candidate_risk_error_vs17=maximum(lambda row:row['finite2048_minus_order17']['selected_risk_abs']['mode1_profile_selected']),
        maximum_finite_reference_selected_risk_error_vs17=maximum(lambda row:row['finite2048_minus_order17']['selected_risk_abs']['mode2_reference_selected']),
        maximum_event_refinement13_17=maximum(lambda row:row['refinement13_17']['event_mass_max']),
        maximum_mean_fraction_refinement13_17=maximum(lambda row:row['refinement13_17']['mean_fraction_prior_width']),
        maximum_SD_relative_refinement13_17=maximum(lambda row:row['refinement13_17']['SD_relative']),
        maximum_candidate_risk_refinement13_17=maximum(lambda row:row['refinement13_17']['selected_risk_abs']['mode1_profile_selected']),
        maximum_replay_ESS_difference=maximum(lambda row:abs(row['finite2048']['ESS_minus_saved'])),
        maximum_replay_reference_risk_difference=maximum(lambda row:abs(row['finite2048']['replay_reference_risk_minus_saved'])),
        maximum_profile_chi2_replay_difference=maximum(lambda row:abs(row['finite2048']['profile_chi2_minus_saved'])),
        peak_RSS_bytes=_peak_bytes(),rows=flat,unique_history_details=evaluated,
        scope='first4 parents only, two noises; duplicate policy histories are not independent evidence. No CDF, order21, full64-case audit, full-Maxwell posterior calibration or continuum certificate.',
        algorithm_and_primary_statistics_unchanged=True,
        interpretation_limit='Observed selected-MAE and candidate-profile decisions are not replaced. Reference probability/risk interpretation must include these empirical numerical discrepancies; final scientific judgment remains with parent session.')
    dump(final_path,report)
    with (results/'HISTORY_REFERENCE_AUDIT.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    short = '\n'.join([
        '# 固定 joint-history 数值审计', '',
        f'4个固定parent、各2noise；32政策标签行，{len(unique)}个唯一history。原proper coefficient prior + persistent bounded gain；只用已选观测。', '',
        f'- finite2048 对17阶最大event质量差：{report["maximum_finite_event_mass_error_vs17"]:.9g}',
        f'- finite2048 对17阶最大mean/prior-width差：{report["maximum_finite_mean_fraction_error_vs17"]:.9g}',
        f'- finite2048 对17阶最大已选candidate expected-L1-risk差：{report["maximum_finite_selected_candidate_risk_error_vs17"]:.9g}',
        f'- 13/17最大event差：{report["maximum_event_refinement13_17"]:.9g}',
        f'- 13/17最大mean/prior-width差：{report["maximum_mean_fraction_refinement13_17"]:.9g}',
        f'- 13/17最大SD相对差：{report["maximum_SD_relative_refinement13_17"]:.9g}',
        f'- UNRESOLVED唯一history：{len(report["unresolved_unique_histories"])} / {len(unique)}',
        f'- 峰值RSS：{report["peak_RSS_bytes"]/1024**2:.1f} MiB；新增Maxwell0、GPU0。', '',
        '这仅检查固定子集、同一连续proxy target的数值近似。未增加21阶/CDF/全案例campaign；17阶本身仍只由13/17比较支持。',
        '原算法、Selected-MAE统计和候选决策保持原样；参考event mass、expected risk与校准主张必须同时披露上述误差及未稳项。主会话保留最终科研裁决。', ''
    ])
    (results/'HISTORY_REFERENCE_AUDIT_SHORT.md').write_text(short)
    print(json.dumps({key:report[key] for key in report if key.startswith('maximum_')}),flush=True)
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--job',default='e2-joint-history-reference-audit-v1')
    arguments=parser.parse_args()
    with Cost().job(arguments.job,'E2') as cost:
        run(cost)
