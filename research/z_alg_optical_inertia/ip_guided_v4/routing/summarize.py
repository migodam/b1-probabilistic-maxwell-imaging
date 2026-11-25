#!/usr/bin/env python3
"""Compile paired performance, frozen L2 receipts and plotting inputs."""
import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from routing_run import HERE, STUDY, LOSS, RING, SRAP, read, write


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path,data):
    with path.open('w') as f:
        writer = csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in data for k in row)))
        writer.writeheader()
        writer.writerows(data)


def bounds():
    output = defaultdict(dict)

    def add(ident,cap,d,r,profile,source,t=None):
        cap,d,r = float(cap),float(d),float(r)
        previous = output[ident].get(cap)
        if previous is None or d>previous['distance_lower']:
            output[ident][cap] = dict(distance_lower=d,radius=r,distance_over_radius=d/r,
                exceeds_radius=d>r,cap_upper=cap,profile=profile,receipt=source,
                archived_receipt_wall_seconds=t)

    src = STUDY/'submission/curve/noise_curve.csv'
    for row in rows(src):
        d = max([float(row[k]) for k in ('strict_global_distance_lower','strict_local_distance_lower',
                                        'fixed_0_distance_lower','fixed_2_distance_lower') if row[k]] or [0.])
        add(row['id'],row['cap'],d,row['radius'],'largest applicable original one-sided full-wave profile',
            str(src.relative_to(STUDY)))
    src = SRAP/'cases.csv'
    for row in rows(src):
        add(row['id'],row['cap'],row['joint_distance_lower'],row['padded_radius'],row['profile'],
            str(src.relative_to(STUDY)))
    src = RING/'paired_results.csv'
    for row in rows(src):
        if row['arm'].endswith('_joint') and row['certified_lower']:
            add(row['case']+'__24to24-12tx-48rx',row['cap'],row['certified_lower'],row['radius'],row['arm'],
                str(src.relative_to(STUDY)))
    for i in range(2):
        src = LOSS/f'finite_{i}_state_verified.json'
        for row in read(src)['rows']:
            for arm,data in row['state'].items():
                if data and 'distance_lower' in data:
                    add(row['id'],row['cap'],data['distance_lower'],row['rho'],arm,
                        str(src.relative_to(STUDY)),data.get('wall_seconds'))
        src = LOSS/f'finite_{i}_verified.json'
        for group in read(src)['projected']:
            for row in group['rows']:
                add(row['id'],group['cap'],row['dual_lower'],row['rho'],'verified receiver ellipsoid dual',
                    str(src.relative_to(STUDY)),row['wall_seconds'])
        src = LOSS/f'finite_{i}_directional.json'
        if not src.exists():
            continue
        for group in read(src)['profiles']:
            for row in group['rows']:
                for name,result in row['methods'].items():
                    if 'distance_lower' in result:
                        add(row['id'],group['cap'],result['distance_lower'],row['rho'],name,
                            str(src.relative_to(STUDY)),result.get('wall_seconds'))
    return output


def stat_low(x):
    x = [z for z in x if z is not None]
    if not x:
        return None
    return dict(minimum=float(np.min(x)),median=float(np.median(x)),maximum=float(np.max(x)))


def report_rows(source,cie_id,b):
    f = HERE/'online'/(source+'.json')
    a = read(f)
    gt = read(HERE/'offline'/(source+'.json'))
    ci = HERE/cie_id/(source+'.json')
    radas = a['radius']
    field, l1 = a['field_norm'], a['born_certificate']
    subset = list(b[source].values())
    highest = max((p['cap_upper'] for p in subset if p['exceeds_radius']),default=None)
    for p in subset:
        if abs(p['radius']/radas-1)>2.e-14:
            raise ValueError(f'Inherited L2 radius mismatch {source} {p["radius"]} {radas}')
    base = dict(id=source,object_id=a['object_id'],cohort=a['cohort'],noise_pct=a['noise_pct'],
                physical_grid=a['forward_grid'],inversion_grid=a['inversion_grid'],
                width_m=a['width_m'],frequency_GHz=a['frequency_Hz']/1e9 if a['frequency_Hz'] else None,
                ka=(2*np.pi*a['frequency_Hz']/3e8*a['width_m']/np.sqrt(2)) if a['frequency_Hz'] else None,
                tangent=a['tangent'],n_tx=a['n_tx'],n_rx=a['n_rx'],radius=radas,field_norm=field,
                gt_peak=gt['ground_truth_peak'],
                true_full_wave_training_fit_over_radius=gt['clean_train_fit_to_noisy_over_radius'],
                holdout=gt['holdout_kind'],l1_lower=l1['lower'],
                l1_lower_relative=None if l1['lower'] is None else l1['lower']/field,
                l1_lower_over_radius=None if l1['lower'] is None else l1['lower']/radas,
                l1_excludes=l1['excludes'],l1_online_wall_seconds=l1['wall_seconds'],
                l1_online_cpu_seconds=l1['cpu_seconds'],l2_highest_excluded_cap=highest,
                l2_attempt_status='strict profiles checked, no class decision' if subset and highest is None else
                    'strictly excluded listed cap' if highest is not None else 'no matching stored L2 proof',
                l2_cap05_distance_over_radius=b[source].get(.5,{}).get('distance_over_radius'),
                l2_cap075_distance_over_radius=b[source].get(.75,{}).get('distance_over_radius'),
                l2_cap1_distance_over_radius=b[source].get(1.,{}).get('distance_over_radius'),
                l2_cap5_distance_over_radius=b[source].get(5.,{}).get('distance_over_radius'),
                linear_born_fit_residual_over_radius=a['born_linear_fit']['residual']/radas,
                linear_born_fit_wall_seconds=a['born_linear_fit']['wall_seconds'],
                linear_born_fit_cpu_seconds=a['born_linear_fit']['cpu_seconds'],
                route_selected=a['route']['selected'],route_online_wall_seconds=a['route']['charged_wall_seconds'],
                route_online_cpu_seconds=a['route']['charged_cpu_seconds'],source=a['source'])
    if a['cohort']=='r14_ring7':
        with np.load(RING/'public_inputs.npz',allow_pickle=False) as common:
            base['ka'] = float(common['k'])*a['width_m']/np.sqrt(2)
            base['wavenumber_from_original_units'] = float(common['k'])
    methods = {}
    for k,v in a['arms'].items():
        result = v.get('audit',{})
        calls = result.get('counts',{}).get('native_full_wave_calls',0)
        final = result.get('residual')
        t = {name:v[name] for name in ('wall_seconds','cpu_seconds')}
        if k in ('born_only','born_then_full'):
            t = {name:val+a['born_linear_fit'][name] for name,val in t.items()}
        methods[k] = dict(cost=t,full_calls=calls+(v.get('fit') or {}).get('counts',{}).get('native_full_wave_calls',0),
                          residual=final,material_update_grid=a['inversion_grid'],failure=v.get('failure'))
    if ci.exists():
        for v in read(ci)['arms']:
            k = v['arm']
            result = v.get('audits',[{}])[-1].get('native',{})
            t = {name:v[name] for name in ('wall_seconds','cpu_seconds')}
            if k=='cie_from_born':
                t = {name:val+a['born_linear_fit'][name] for name,val in t.items()}
            methods[k] = dict(cost=t,full_calls=result.get('counts',{}).get('native_full_wave_calls',0),
                              residual=result.get('residual'),material_update_grid=96,
                              failure=result.get('failure',v.get('failure')),pr_steps=v.get('step'))
    pairs = []
    for k,v in methods.items():
        r = v['residual']
        pairs.append(dict(**base,arm=k,arm_inversion_grid=v['material_update_grid'],pr_steps=v.get('pr_steps'),
                          online_wall_seconds=v['cost']['wall_seconds'],online_cpu_seconds=v['cost']['cpu_seconds'],
                          full_wave_calls=v['full_calls'],true_residual_over_radius=None if r is None else r/radas,
                          fits=r<=radas if r is not None else False,solver_failure=v['failure'],
                          **gt['arms'].get(k,{})))
    res = {v['arm']:v for v in pairs}
    chosen = res[base['route_selected']]
    simple = res['full_zero']
    warm = res['born_then_full']
    end = dict(**base,route_fit=chosen['fits'],route_true_residual_over_radius=chosen['true_residual_over_radius'],
               route_material_error=chosen.get('full_relative_l2'),route_internal_error=chosen.get('inside_relative_l2'),
               route_background_rms=chosen.get('background_rms'),route_held_error=chosen.get('held_error'),
               route_full_calls=chosen['full_wave_calls'],
               zero_true_residual_over_radius=simple['true_residual_over_radius'],
               warm_true_residual_over_radius=warm['true_residual_over_radius'],
               zero_material_error=simple.get('full_relative_l2'),warm_material_error=warm.get('full_relative_l2'),
               zero_wall_seconds=simple['online_wall_seconds'],warm_wall_seconds=warm['online_wall_seconds'],
               zero_cpu_seconds=simple['online_cpu_seconds'],warm_cpu_seconds=warm['online_cpu_seconds'],
               cie_zero_true_residual_over_radius=res.get('cie_from_zero',{}).get('true_residual_over_radius'),
               cie_born_true_residual_over_radius=res.get('cie_from_born',{}).get('true_residual_over_radius'),
               cie_zero_wall_seconds=res.get('cie_from_zero',{}).get('online_wall_seconds'),
               cie_born_wall_seconds=res.get('cie_from_born',{}).get('online_wall_seconds'),
               cie_zero_cpu_seconds=res.get('cie_from_zero',{}).get('online_cpu_seconds'),
               cie_born_cpu_seconds=res.get('cie_from_born',{}).get('online_cpu_seconds'))
    return end,pairs,subset


def group_summary(data):
    return dict(views=len(data),independent=len({r['object_id'] for r in data}),
                l1_exclusions=sum(r['l1_excludes'] for r in data),
                l1_lower_relative=stat_low([r['l1_lower_relative'] for r in data]),
                l1_median_instance_wall_seconds=stat_low([r['l1_online_wall_seconds'] for r in data]),
                routes=dict(Counter(r['route_selected'] for r in data)),
                required_l2=dict(Counter(str(r['l2_highest_excluded_cap']) for r in data)),
                route_success=sum(r['route_fit'] for r in data),
                warm_full_residual_wins=sum(1 for r in data if r['zero_true_residual_over_radius'] is not None and
                                           r['warm_true_residual_over_radius'] is not None and
                                           r['warm_true_residual_over_radius']<r['zero_true_residual_over_radius']),
                warm_material_wins=sum(1 for r in data if r['zero_material_error'] is not None and
                                           r['warm_material_error'] is not None and
                                           r['warm_material_error']<r['zero_material_error']),
                sum_wall_seconds={f:float(sum(r[f] or 0. for r in data)) for f in
                                ('zero_wall_seconds','warm_wall_seconds','route_online_wall_seconds',
                                 'cie_zero_wall_seconds','cie_born_wall_seconds')},
                sum_cpu_seconds={f:float(sum(r[f] or 0. for r in data)) for f in
                                ('zero_cpu_seconds','warm_cpu_seconds','route_online_cpu_seconds',
                                 'cie_zero_cpu_seconds','cie_born_cpu_seconds')},
                medians_by_method={f:stat_low([r[f] for r in data]) for f in
                                 ('route_true_residual_over_radius','route_material_error','route_internal_error',
                                  'route_background_rms','route_held_error','zero_true_residual_over_radius',
                                  'warm_true_residual_over_radius','cie_zero_true_residual_over_radius',
                                  'cie_born_true_residual_over_radius')})


def lab_input():
    records = []
    for x in sorted((HERE/'laboratory').glob('*.json')):
        if x.stem.endswith('__seal'):
            continue
        a = read(x)
        p = a['previous_fixed_born_certificate']
        metrics = a['port_fit_metrics']
        levels = [v for v in a['previous_conditional_peak_margins']
                  if v['material_family']=='passive_sector' or 'fresnel' in a['id']]
        captures = [v['cap'] for v in levels if (v.get('remaining_proxy_relative_margin') or
                                                v.get('remaining_proxy_margin') or 0.)>0.]
        records.append(dict(id=a['id'],object_id=a['object_id'],frequency_GHz=a['frequency_GHz'],
                            train_channels=metrics['train']['channels'],train_error=metrics['train']['relative_rms'],
                            held_tx_channels=metrics['held_source']['channels'],
                            held_tx_error=metrics['held_source']['relative_rms'],
                            held_rx_channels=metrics['held_receiver']['channels'],
                            held_rx_error=metrics['held_receiver']['relative_rms'],
                            geometry_verified_on_previous_selected_ports=p['geometry_valid'],
                            historical_born_lower_relative=p['lower_relative_to_selected_field'],
                            historical_exceeds_both_error_proxies=p['exceeds_both_empirical_scales'],
                            highest_historical_conditional_cap=max(captures,default=None),
                            strict_loaded_error_radius=None,material_error=None))
    return records


def cie_input():
    records = []
    for file in sorted((HERE/'cie_diagnostic').glob('*.json')):
        summary = read(file)
        scored = HERE/'offline'/file.name
        offline = read(scored)['arms'] if scored.exists() else {}
        for arm,info in summary['methods'].items():
            for audit in info['audits']:
                k = audit['step']
                d = audit['native'].get('residual')
                raw = audit.get('native_mode_contents',{})
                capped = audit.get('clipped_to_cap_one',{})
                result = offline.get(f'{arm}_{k}',{})
                records.append(dict(id=summary['id'],object_id=summary['object_id'],arm=arm,step=k,
                    n_rx=summary['n_rx'],n_tx=summary['n_tx'],rank=info['rank'],beta=info['beta'],
                    required_raw_current_norm_cap_one=summary['required_current']['negative_current_norm_lower'],
                    native_residual_over_radius=None if d is None else d/summary['radius'],
                    native_peak=audit['native'].get('peak'),native_failure=audit['native'].get('failure'),
                    material_l2=result.get('full_relative_l2'), held_error=result.get('held_error'),
                    raw_negative_current_norm=raw.get('raw_negative_mode_projection_norm'),
                    weighted_centered_negative_coordinate_energy=raw.get('weighted_centered_negative_coordinate_energy'),
                    cap_one_residual_over_radius=capped['audit']['residual']/summary['radius'] if capped else None,
                    cap_one_raw_current_norm=capped.get('raw_negative_mode_projection_norm'),
                    total_arm_wall_seconds=info['wall_seconds'],total_arm_cpu_seconds=info['cpu_seconds']))
    return records


def early_exit_input():
    records = []
    for dirname in ('cie_early_exit','cie_early_exit_narrow'):
        for file in sorted((HERE/dirname).glob('*.json')):
            a = read(file)
            strict = a['strict_box']
            r = a['radius']
            cont = a['arms']['continue_to_16']['audits'][-1]['native']['residual']/r
            rest = a['arms']['born_restart']['audits'][-1]['native']['residual']/r
            records.append(dict(id=a['id'],object_id=a['object_id'],box_halfwidth=a['material_box_halfwidth'],
                point_residual_over_radius=a['one_point_relative_discrepancy'],
                strict_box_distance_over_radius=strict.get('distance_over_radius_lower'),
                neumann_distance_over_radius=strict.get('center',{}).get('neumann_residual_over_radius_lower'),
                l4_excludes=strict['strictly_excludes'],discrepancy_action=a['simple_discrepancy']['selected'],
                l4_action=a['l4_rule']['selected'],continued_cie_final_residual_over_radius=cont,
                restarted_cie_final_residual_over_radius=rest,
                number_of_continued_cells_outside_box=a['unrestricted_continuation_cells_outside_box'],
                continuation_linf_distance=a['unrestricted_continuation_linf_from_checkpoint'],
                l4_cost_seconds=strict['total_new_cost']['wall_seconds'],
                l4_full_path_wall_seconds=a['l4_rule']['charged']['wall_seconds'],
                discrepancy_full_path_wall_seconds=a['simple_discrepancy']['charged']['wall_seconds']))
    return records


def main():
    b = bounds()
    data, arms, proofs = [], [], []
    cohorts = defaultdict(list)
    seals = sorted((HERE/'online').glob('*__seal.json'))
    for seal in seals:
        for obs in read(seal)['ids']:
            case,bs,ls = report_rows(obs,'cie_principal_common_box',b)
            data.append(case)
            cohorts[case['cohort']].append(case)
            arms.extend(bs)
            proofs.extend(dict(id=obs,**ln) for ln in ls)
    main = [r for r in data if r['cohort'] in ('optfresh','srap22')]
    assert len(main)==28 and len({r['object_id'] for r in main})==28
    write_csv(HERE/'principal_28.csv',main)
    write_csv(HERE/'all_routing_cases.csv',data)
    write_csv(HERE/'paired_arms.csv',arms)
    write_csv(HERE/'historical_l2.csv',proofs)
    lab = lab_input()
    write_csv(HERE/'measured_port_cases.csv',lab)
    write_csv(HERE/'cie_diagnostics.csv',cie_input())
    write_csv(HERE/'early_exit_diagnostics.csv',early_exit_input())
    answer = dict(primary=group_summary(main),cohorts={k:group_summary(rows) for k,rows in cohorts.items()},
                  measured_views=len(lab),measured_physical_objects=len({r['object_id'] for r in lab}),
                  post_seal_only=['ground_truth','offset synthetic ports','actual held source/acquired receiver fields'],
                  warm_gate_wall_seconds_once_by_suite={q.stem:read(q)['setup'] for q in seals},
                  l2_source='frozen historical strict receipts; no new L2 time claim')
    write(HERE/'summary.json',answer)
    print(answer['primary'])


if __name__=='__main__': main()
