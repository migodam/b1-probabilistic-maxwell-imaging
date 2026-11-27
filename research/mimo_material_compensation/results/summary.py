"""Collect the distinct claims and their source receipts without merging noise contracts."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent.parent

def read(path):
    return json.loads((ROOT/path).read_text())

def outputs():
    paths=read('paths/summary.json')
    topology=read('topology/fixed_pair/summary.json')
    allocation=read('environment/band_allocation.json')
    post=read('environment/postmeasurement/REPORT_RECEIPT.json')
    six=read('controls/same_data_path.json')
    one_pulse={}
    for side in (32,96):
        base='environment/fluence_modes'
        classes=read(f'{base}/classes_{side}.json')
        proposal=read(f'{base}/n{side}/proposal.json')
        verify=read(f'{base}/n{side}/independent_check.json')
        keys={'class_1','center_1'}
        one_pulse[str(side)]=dict(
            one_shot_directions={a['key']:a['fullclass_lower_over_2r'] for a in classes['rows'] if a['key'] in keys},
            added_energy=proposal['available_energy_upper'],new_reads=proposal['complex_reads'],
            all_sealed_labels_matched=verify['matched_labels'],whole_class_allowance=proposal['beta'],
            sources=f'{base}/n{side}/proposal.json')
    counter=read('environment/fluence_modes/reattack/dual_verification.json')
    return dict(
        physical_scope='Nominal discrete 2-D TM, fixed passive coupled ports. Real nondispersive materials in [0,4.5]. No new chamber readings.',
        source_pulses='Equal duration, available source power integrated over a common duration. Absolute deterministic Frobenius balls unless specifically stated.',
        exact_algebraic_path={k:six[k] for k in ('all_intervals_verified','precision_bits','number_intervals','physical_ROI_endpoints','verified_width_if_all')},
        fixed_witness_competition=dict(independent_reference_objects=len(topology['cases']),undirected_cuts=topology['cut_count'],leave_one_material_out=topology['heldout_medians']),
        actual_old_ball_paths=[{k:r[k] for k in ('case','method','scope','area_target_m2','segments','mean_task_lower','mean_region_lower','residual_upper_over_r')}
                               for r in paths['paths'] if r['scope'].startswith('original noisy')],
        specification_energy=read('environment/decision_limits.json')['records'],
        prospective_one_complex_number=one_pulse,
        same_native96_one_pulse_wide_box=dict(
            pixels=9216,area_m2=counter['area_m2'],core_uniform=[0,2],frequency_MHz=[400,460],
            actual_data_ids={k:dict(inner_diameter=v['diameter_lower'],
                                    accepted={m:w['residual_upper_over_r'] for m,w in v['materials'].items() if w['both_balls']})
                             for k,v in counter['samples'].items()},
            exterior_budgets={k:{s:r[s] for s in ('excluded_exterior_Linf_through','verified_witness_Linf_upper','area_L1_witness_upper_m2','area_L2_squared_witness_upper_m2')}
                              for k,r in read('environment/fluence_modes/reattack/background_budget.json')['rows'].items()},
            native_replay='environment/fluence_modes/reattack/dual_verification.json',
            independent_loaded_24_port_replay='environment/fluence_modes/reattack/loaded_verification.json'),
        mixed_frequency_common_history=allocation['records'],
        recorded_two_ball_refits=[dict(id=k,core=v['core'],old_upper=v['independent_audit']['400']['strict_residual_over_r_upper'],
                                      new_upper=v['independent_audit']['460']['strict_residual_over_r_upper'],both_balls=v['both_balls']) for k,v in post['arms'].items()],
        replay= dict(source='replication/replay_receipts/replay_receipt.json',covered=['whole noisy path','two full specification classes','all 6435 fixed cuts on seven targets','fresh postmeasurement material fits']),
    )

if __name__=="__main__":
    result=outputs()
    (ROOT/'results'/'main.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
