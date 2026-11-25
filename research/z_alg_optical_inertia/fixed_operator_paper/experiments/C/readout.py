"""Read saved predictive supports first, then update the joint full-wave test."""
import csv
from collections import Counter,defaultdict
import json
from time import time

from flint import arb
import numpy as np
from threadpoolctl import threadpool_limits

from design import form_border, joint_tail, real_trace, source_pair
from model import (BATCHES,COMPONENT,CONTRACTS,HERE,OLD,POLICIES,WEIGHTS,ball,clocks,elapsed,energy,
                   geometry,hash_arrays,hash_file,high,init_arithmetic,low,noise_radius,radius_squared,read_json,
                   require_p0,star,write_json)


def read_array(path,expected):
    with np.load(path,allow_pickle=False) as raw:
        y=raw['Y']
    if hash_arrays(y)!=expected:
        raise ArithmeticError(f'Frozen acquisition content changed: {path}')
    return y


def apply_saved_support(action,archive,ynew):
    start=clocks()
    proof=action['verified_support']
    if proof is None or not proof['support_query_index']:
        return dict(status='NO_SAVED_SUPPORT_AVAILABLE',excluded=[],**elapsed(start))
    aw=archive[proof['support_query_index']]
    if hash_arrays(aw)!=proof['support_query_hash']:
        raise ArithmeticError('A frozen support direction changed.')
    w=ball(aw)
    projection=real_trace(star(w)*ball(ynew))
    columns=[low(projection),high(projection)]
    result={}
    for label, info in proof['supports'].items():
        lo,hi=info['noisy_support_interval']
        gap=max(low(arb(lo)-arb(columns[1])),low(arb(columns[0])-arb(hi)))
        result[label]=dict(status='EXCLUDED' if gap>0 else 'UNRESOLVED',
                           support_gap_lower=gap,declared_noisy_interval=[lo,hi])
    return dict(status='SAVED_SUPPORTS_TESTED_BEFORE_JOINT_UPDATE',measured_projection_interval=columns,
        by_class=result,excluded=[cl for cl,r in result.items() if r['status']=='EXCLUDED'],**elapsed(start))


def joint_update(borders,z,labels,r2):
    start=clocks()
    result={}
    r=r2.sqrt()
    for label in labels:
        checked=[]
        for b in borders[label]:
            lo,hi,error=joint_tail(b,z)
            checked.append(dict(weight=b.weight,negative_index=b.negative,
                distance_lower=max(0.,low(lo.sqrt())),distance_upper=high(hi.sqrt()),
                margin_lower=low(lo.sqrt()-r),excluded=(lo>r2),
                source_eigensystem_defect_upper=error))
        best=max(checked,key=lambda x:x['distance_lower'],default=None)
        excludes=[x for x in checked if x['excluded']]
        result[label]=dict(status='EXCLUDED' if excludes else 'UNRESOLVED',
                           best=best,checked=checked)
    return dict(by_class=result,excluded=[cl for cl,r in result.items() if r['status']=='EXCLUDED'],
                **elapsed(start))


def prefix_reports(record,rows):
    reasons=Counter()
    counts=Counter()
    for plan in record['source_rows']:
        for label,evidence in plan['old_classes'].items():
            counts['candidate_class_contract_instances']+=1
            if evidence['chosen_envelope']:
                counts['candidate_class_has_bounded_envelope']+=1
            for proof in evidence['profiles']:
                base=proof.get('base','state_unresolved')
                counts[f'candidate_{base}_attempted']+=1
                if proof.get('saturated'): counts[f'candidate_{base}_saturated']+=1
                if proof.get('receiver_positive'): counts[f'candidate_{base}_receiver_pd']+=1
                reasons[f'candidate_{proof["stage"]}']+=1
        for proof in plan['indefinite_stress_profiles']:
            base=proof.get('base','state_unresolved')
            counts[f'stress_{base}_attempted']+=1
            if proof.get('saturated'): counts[f'stress_{base}_saturated']+=1
            if proof.get('receiver_positive'): counts[f'stress_{base}_receiver_pd']+=1
            reasons[f'stress_{proof["stage"]}']+=1
        for data in plan['batches'].values():
            if data['status']=='STOPPED_OLD_SINGLETON':
                counts['old_singleton_budgets']+=1
                continue
            counts['all_menu_actions_across_conditions']+=data['finite_menu_size']
            counts['menu_actions_outward_checked']+=data['attempted_interval_actions']
            counts['menu_actions_strictly_separating']+=data['strict_separating_actions']
    counters={}
    for policy in POLICIES:
        for contract in CONTRACTS:
            for size in BATCHES:
                matching=[r for r in rows if (r['policy'],r['energy_contract'],r['batch_size'])==(policy,contract,size)]
                counters[f'{policy}|{contract}|{size}']=dict(
                    tested_target_count=len(matching),
                    wrong_class_excluded=sum(r['correct_single_class'] for r in matching),
                    removed_by_saved_support=sum(r['support_removed_wrong_class'] for r in matching),
                    removed_by_saved_support_exclusively=sum(
                        r['support_removed_wrong_class'] and not r['joint_removed_wrong_class']
                        for r in matching),
                    removed_by_joint=sum(r['joint_removed_wrong_class'] for r in matching),
                    strict_predicted_edge=sum(r['planned_strict_edge'] for r in matching),
                    mean_survivors=float(np.mean([r['surviving_classes'] for r in matching])) if matching else None)
    best_costs={}
    for contract in CONTRACTS:
        for policy in POLICIES:
            targets = {r['target_id'] for r in rows if (r['policy'],r['energy_contract'])==(policy,contract)}
            best_costs[f'{policy}|{contract}']={target:min((r['batch_size'] for r in rows if
                   (r['target_id'],r['policy'],r['energy_contract'])==(target,policy,contract)
                   and r['correct_single_class']),default=None) for target in sorted(targets)}
    return dict(prerequisite_counts=counts,prerequisite_stage_counts=reasons,by_policy_contract_batch=counters,
                minimum_successful_tested_batch=best_costs)


def main():
    start=clocks()
    gate=require_p0()
    plan=read_json(HERE/'proposal.json')
    acq=read_json(HERE/'acquired/acquisition.json')
    old_meta=read_json(HERE/'old/old_acquisition.json')
    registry=read_json(HERE/'registration.json')
    if (not (plan['wall_time_unix']<acq['acquisition_end_unix'])
            or acq['frozen_proposal_hash']!=hash_file(HERE/'proposal.json')):
        raise ArithmeticError('Future acquisition did not follow proposal.')
    rows=[]
    raw_verify={}
    shared={}
    with np.load(HERE/'public_model.npz',allow_pickle=False) as pub, \
            np.load(HERE/'proposal_arrays.npz',allow_pickle=False) as archives:
        bank={}
        for design in plan['source_rows']:
            target,contract=design['scene'],design['energy_contract']
            info=registry['scenes'][target]
            label_list=info['classes']
            named=f'{target}__{contract}'
            old=read_array(HERE/'old'/f'{named}.npz',design['old_data_hash'])
            if target not in bank:
                t=clocks()
                g=geometry(info['layout'])
                l,u=pub[target+'__lo'],pub[target+'__hi']
                bank[target]={label:[b for w in WEIGHTS if (b:=form_border(g,lo,hi,w)[0]) is not None]
                                for label,lo,hi in zip(label_list,l,u,strict=True)}
                shared[target]=dict(certified_borders=sum(len(x) for x in bank[target].values()),**elapsed(t))
            for key,batch in design['batches'].items():
                size=int(key)
                eo,en=energy(contract,size)
                for policy,action in batch['actions'].items():
                    s=clocks()
                    aid=f'{named}__m{size}__{policy}'
                    truth=acq['policy_rows'][aid]
                    physical=acq['reading_files'][truth['acquisition_key']]
                    future=read_array(HERE/'acquired'/physical['filename'],physical['data_hash'])
                    support=apply_saved_support(action,archives,future)
                    # This is deliberately after the stored predictive test.
                    z=source_pair((*OLD,*action['source_indices']),arb(1),np.column_stack((old,future)))
                    # Replace each coefficient with its exact physical amplitude,
                    # without touching any of the old observed data.
                    for i,j in enumerate(OLD):
                        z[j,i]=eo.sqrt()
                    for i,j in enumerate(action['source_indices'],len(OLD)):
                        z[j,i]=en.sqrt()
                    r2=radius_squared(len(OLD)+size,noise_radius(registry,target))
                    updated=joint_update(bank[target],z,label_list,r2)
                    expected=truth['witness_class']
                    with np.load(COMPONENT/'.sealed'/f'{target}.npz',allow_pickle=False) as f:
                        x=f['material']
                    owner=label_list.index(expected)
                    if hash_arrays(x)!=truth['material_hash'] or not np.all(
                            (pub[target+'__lo'][owner]<=x)&(x<=pub[target+'__hi'][owner])):
                        raise ArithmeticError('Offline reference identity/box verification failed.')
                    err2=arb(old_meta['observations'][named]['exact_forward_distance_upper'])**2 + \
                         arb(physical['exact_forward_distance_upper'])**2
                    if not err2<r2:
                        raise ArithmeticError('Paired original forward witness failed joint radius.')
                    removed_old=set(label_list)-set(design['coverage_before_update'])
                    removed_support=set(support['excluded'])
                    removed_joint=set(updated['excluded'])
                    removed=removed_old|removed_support|removed_joint
                    if expected in removed:
                        write_json(HERE/'false_exclusion_FAILURE.json',dict(id=aid,support=support,joint=updated))
                        raise ArithmeticError('Strict false exclusion of true material, stopped experiment.')
                    survivors=[label for label in label_list if label not in removed]
                    wrong=label_list[1-owner]
                    strict_edges=bool(action['verified_support'] and action['verified_support']['strict_edge'])
                    if strict_edges and len(removed_support)==0:
                        raise ArithmeticError('Saved strict disjoint-support edge did not execute.')
                    result=dict(id=aid,target_id=target,pair_id=target,condition_id=f'{named}__m{size}',
                       material_hash=truth['material_hash'],layout=design['layout'],energy_contract=contract,
                       batch_size=size,actual_new_measurements=size,actual_old_measurements=len(OLD),
                       policy=policy,old_data_hash=truth['old_data_hash'],acquired_data_hash=physical['data_hash'],
                       new_source_indices=action['source_indices'],new_source_degrees=action['source_angles_degrees'],
                       new_total_energy=action['total_new_energy'],
                       total_energy=action['old_total_energy']+action['total_new_energy'],
                       r_joint=high(r2.sqrt()),true_witness_distance_upper=high(err2.sqrt()),
                       planned_strict_edge=strict_edges,found_expected_witness=True,
                       support_removed_wrong_class=(wrong in removed_support),
                       joint_removed_wrong_class=(wrong in removed_joint),
                       correct_single_class=(survivors==[expected]),surviving_classes=len(survivors),
                       predicted_cover=list(design['coverage_before_update']),survivor_names=survivors,
                       class_state={label:('EXCLUDED' if label in removed else
                                           'WITNESS_FOUND' if label==expected else 'UNRESOLVED') for label in label_list},
                       **elapsed(s))
                    raw_verify[aid]=dict(old_support=support,joint_update=updated,**result)
                    rows.append(result)
    result=dict(status='PROSPECTIVE_POLICY_REPLAY_COMPLETE',gate=gate,
                comparisons_are_conditional_on_registered_two_box_prior=True,
                distinct_targets=len({r['target_id'] for r in rows}),
                policy_paths=len(rows),true_witness_retained=sum(r['found_expected_witness'] for r in rows),
                interval_bits=128,shared_joint_borders=shared,
                unique_measurement_batches=acq['unique_physical_evaluations'],
                **prefix_reports(plan,rows),**elapsed(start))
    write_json(HERE/'results.json',result)
    write_json(HERE/'interval_readback.json',raw_verify)
    with (HERE/'outcomes.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    init_arithmetic()
    with threadpool_limits(1):
        main()
