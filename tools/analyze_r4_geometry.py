#!/usr/bin/env python3
"""Read completed R4 geometry arrays/rows only; no physics or GT imports."""
from pathlib import Path
from collections import Counter,defaultdict
import json,numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r4_conditional_transport'
RESULT=OUT/'results/III/dynamic_geometry_first'


def stats(values):
    a=np.asarray(values,float)
    return dict(n=len(a),minimum=float(a.min()),q25=float(np.quantile(a,.25)),median=float(np.median(a)),q75=float(np.quantile(a,.75)),maximum=float(a.max())) if len(a) else dict(n=0)


def main():
    receipt=json.loads((OUT/'III/jobs/dynamic_geometry_first/CLOSED.json').read_text())
    result=json.loads((RESULT/'RESULTS.json').read_text())
    if receipt['status']!='COMPLETE' or result['status']!='COMPLETE':raise ValueError('completed physical job required')
    configs=[json.loads(p.read_text()) for p in RESULT.glob('R4DEV_*/d*_r*/SUMMARY.json')]
    skips=[];fds=[];convergence=[];group=defaultdict(list);allarms=[]
    for path in RESULT.glob('R4DEV_*/d*_r*/SKIPS.jsonl'):
        for line in path.read_text().splitlines():skips.append(dict(configuration=str(path.parent.relative_to(RESULT)),**json.loads(line)))
    for path in RESULT.glob('R4DEV_*/d*_r*/BRACKET_ROWS.jsonl'):
        for line in path.read_text().splitlines():fds.append(json.loads(line))
    for path in RESULT.glob('R4DEV_*/d*_r*/FD_CONVERGENCE.json'):
        for item in json.loads(path.read_text()):convergence.append(dict(configuration=str(path.parent.relative_to(RESULT)),**item))
    for c in configs:
        for row in c.get('arm_rows',[]):
            item=dict(parent=c['parent'],dimension=c['dimension'],rank=c['rank'],**row)
            group[(c['dimension'],c['rank'],row['arm'],row['radius'])].append(item);allarms.append(item)
    bracket_groups=[]
    for d in (64,66):
      for rank in (2,4,8):
       for h in (.001,.002,.004):
        rows=[z for z in fds if (z['dimension'],z['rank'],z['h'])==(d,rank,h)]
        bracket_groups.append(dict(dimension=d,rank=rank,h=h,count=len(rows),
            material_amplitude=stats([z['material_amplitude'] for z in rows]),
            normal_material_amplitude=stats([z['material_normal_amplitude'] for z in rows]),
            normal_fraction=stats([z['material_normal_fraction'] for z in rows]),
            latent_norm=stats([z['latent_norm'] for z in rows]),
            any_seam_count=sum(z['any_seam_crossing'] for z in rows),rank_change_count=sum(z['any_rank_change'] for z in rows),
            near_crossing_count=sum(z['any_near_crossing'] for z in rows),
            projector_rotation=stats([z['max_projector_rotation'] for z in rows])))
    arm_groups=[]
    for (d,rank,arm,radius),rows in sorted(group.items()):
        arm_groups.append(dict(dimension=d,rank=rank,arm=arm,radius=radius,completed=len(rows),planned=32,
            actual_RMS=stats([z['actual_RMS'] for z in rows]),radius_ratio=stats([z['actual_RMS']/radius for z in rows]),
            matched_within_1pct=sum(abs(z['actual_RMS']/radius-1)<=.01 for z in rows),
            Phi=stats([z['Phi'] for z in rows]),DeltaPhi=stats([z['DeltaPhi'] for z in rows]),
            prior_delta=stats([z['prior_delta'] for z in rows]),normal_RMS=stats([z['actual_endpoint_normal']['actual_normal_RMS'] for z in rows]),
            normal_fraction=stats([z['actual_endpoint_normal']['actual_normal_fraction'] for z in rows]),
            marginal_RHS=stats([z['marginal_cost']['RHS'] for z in rows]),wall_seconds=stats([z['wall_seconds'] for z in rows]),
            endpoint_seam_count=sum(z['seam_crossed_cells']>0 for z in rows),
            sampled_ball_violation_count=sum(z.get('sampled_outside_reference_ball',0)>0 for z in rows)))
    index={(z['parent'],z['dimension'],z['rank'],z['arm'],z['radius']):z for z in allarms}
    pairs=[]
    for z in allarms:
        if z['arm'] not in ('circuit_forward','circuit_reverse'):continue
        direct=index.get((z['parent'],z['dimension'],z['rank'],'direct_bracket_'+z['arm'].split('_')[-1],z['radius']))
        if direct:
            pairs.append(dict(parent=z['parent'],dimension=z['dimension'],rank=z['rank'],radius=z['radius'],circuit_arm=z['arm'],
                circuit_actual_radius_ratio=z['actual_RMS']/z['radius'],direct_actual_radius_ratio=direct['actual_RMS']/z['radius'],
                marginal_RHS_ratio=z['marginal_cost']['RHS']/direct['marginal_cost']['RHS'],
                cold_RHS_ratio_including_all_three_FD_h=(480+z['marginal_cost']['RHS'])/(480+direct['marginal_cost']['RHS']),
                circuit_DeltaPhi=z['DeltaPhi'],direct_DeltaPhi=direct['DeltaPhi'],
                circuit_prior_delta=z['prior_delta'],direct_prior_delta=direct['prior_delta'],
                circuit_actual_normal_RMS=z['actual_endpoint_normal']['actual_normal_RMS'],direct_actual_normal_RMS=direct['actual_endpoint_normal']['actual_normal_RMS'],
                circuit_max_sampled_DeltaPhi=z['max_sampled_Phi']-(z['Phi']-z['DeltaPhi']),
                same_actual_material_precision_claim=False,continuous_path_certificate=False))
    reason_counts=Counter(z['reason'] for z in skips)
    circuit=[z for z in allarms if z['arm'].startswith('circuit')]
    matched=[z for z in allarms if not z['arm'].startswith('circuit')]
    summary=dict(status='DESCRIPTIVE_COMPLETED_GEOMETRY_PARENT_SCIENTIFIC_REVIEW_REQUIRED',
        independent_new_development_parents=32,configurations=len(configs),planned_configurations=192,
        configuration_status_counts=dict(Counter(z['status'] for z in configs)),
        saved_Full_endpoints=result['endpoint_rows'],saved_bracket_diagnostics=len(fds),
        physical_new_costs=result['new_costs'],controlled_wall_seconds=receipt['controlled_wall_seconds'],profile=receipt['profile'],
        maximum_parent_arm_states=max(z['states'] for z in result['arm_costs'].values()),
        skip_count=len(skips),skip_reasons=dict(reason_counts),
        expected_66_strong_not_run=288,
        completed_circuits=len(circuit),circuits_actual_radius_ratio=stats([z['actual_RMS']/z['radius'] for z in circuit]),
        circuits_within_1pct_radius=sum(abs(z['actual_RMS']/z['radius']-1)<=.01 for z in circuit),
        circuits_with_sampled_reference_ball_violation=sum(z['sampled_outside_reference_ball']>0 for z in circuit),
        completed_direct_or_two_segment_controls=len(matched),
        FD_relative_difference_small_middle=stats([z['relative_difference'] for z in convergence if z['h1']==.001]),
        FD_relative_difference_middle_large=stats([z['relative_difference'] for z in convergence if z['h1']==.002]),
        all_FD_any_seam=sum(z['any_seam_crossing'] for z in fds),all_FD_rank_change=sum(z['any_rank_change'] for z in fds),all_FD_near_crossing=sum(z['any_near_crossing'] for z in fds),
        paired_circuit_direct_count=len(pairs),marginal_RHS_ratio=stats([z['marginal_RHS_ratio'] for z in pairs]),
        cold_RHS_ratio=stats([z['cold_RHS_ratio_including_all_three_FD_h'] for z in pairs]),
        bracket_groups=bracket_groups,arm_groups=arm_groups,
        no_GT=True,continuous_path_certificate=False,posterior_MCSE='NOT_RUN',full_Lie_rank='NOT_PROVEN',
        scientific_scope='one frozen pair/config; finite Euler circuits are not exact flows, reversible maps or posterior kernels; repeated configurations do not increase independent n')
    for filename,data in [('OBSERVED_METRICS.json',summary),('PAIRED_CIRCUIT_DIRECT.json',pairs),('ALL_SKIPS.json',skips)]:
        (RESULT/filename).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('bracket_groups','arm_groups')},indent=2))

if __name__=='__main__':main()
