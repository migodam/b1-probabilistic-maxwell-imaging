"""Read reported v3 moment and low-peak physical controls without regenerating data."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]


def load(path):
    return json.loads((STUDY/path).read_text())


def main():
    ring_source = 'validation/max_struct/offline_scalar_scores.json'
    ring = load(ring_source)['rows']['r14_new_0018']
    assert ring['truth_peak'] < 5
    assert ring['full_wave']['full_wave_distance_hi'] <= ring['public_CN99_radius']
    moments_source = 'ip_submission_v2/energy_sdp/verified_alt025.json'
    moments = load(moments_source)
    physical_source = 'ip_submission_v2/energy_bound/receipts/physical_current.json'
    physical = load(physical_source)
    physical_energy = max(row['coherent_energy_lower'] for row in physical['cases'].values())
    matched = {}
    for name, record in moments['observations'].items():
        proof = record['proofs']['alt025.npz']
        assert proof['fit_scale'] == 1
        assert min(proof['weight_eigenvalue_floors']) > 0
        assert max(proof['sector_eigenvalue_ceilings']) < 0
        assert proof['uniform_witness_current_energy_upper'] < physical_energy
        matched[name] = dict(cap=moments['material_upper'], noise_radius=record['inflated_radius'],
                             minimum_covariance_eigenvalue_lower=min(proof['weight_eigenvalue_floors']),
                             maximum_local_sector_eigenvalue_upper=max(proof['sector_eigenvalue_ceilings']),
                             coherent_current_energy_upper=proof['uniform_witness_current_energy_upper'],
                             mean_receiver_defect_interval_upper=proof['interval_mean_receiver_residual_upper'],
                             moment_squared_distance=0)
    settled = dict(scope='Historical public receipt crosscheck; exact toy proof is independent.',
                   small_ring_feasible_control=dict(source=ring_source, observation='r14_new_0018', cap=5,
                                                    peak=ring['truth_peak'], observation_radius=ring['public_CN99_radius'],
                                                    full_wave_distance_upper=ring['full_wave']['full_wave_distance_hi'],
                                                    margin_lower=ring['coverage_margin_lower']),
                   four_region_material_cap=moments['material_upper'],
                   four_region_moments_with_energy_cap=dict(source=moments_source, **matched),
                   source_for_admissible_physical_current_energy_lower=physical_source,
                   admissible_physical_current_energy_lower=physical_energy,
                   original_joint_lower_distance_source='ip_submission_v3/experiments.tex',
                   original_joint_lower_distances=[0.16069878, 0.11060787])
    (HERE/'historical_controls.json').write_text(json.dumps(settled, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
