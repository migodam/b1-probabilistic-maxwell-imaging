"""Verify parent metadata corrections from current saved derived records."""
from pathlib import Path
import json
import time
import sys
wall_start = time.perf_counter()
cpu_start = time.process_time()
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
ledger = (ROOT / 'results/PHYSICS_LEDGER.jsonl').read_text()
receipt = json.loads((ROOT / 'results/cli_gain_regression_corrected/RUN_RECEIPT.json').read_text())
target = json.loads((ROOT / 'results/cases/refined_gain_marginalized/probabilities.json').read_text())
calibration = json.loads((ROOT / 'results/CALIBRATION_RESULTS.json').read_text())
primary = json.loads((ROOT / 'results/PRIMARY_SUMMARY.json').read_text())
initial_reanalysis = json.loads((OUT / 'RECOMPUTATION_NATIVE.json').read_text())
assert receipt['measurement_target']['gain_sd'] == .15
assert receipt['new_maxwell_evaluations'] == 0
assert target['likelihood']['gain_sd'] == .15
assert np.array_equal(receipt['result']['topology_probabilities'], target['result']['topology_probabilities'])
comparison = {}
with np.load(ROOT / 'results/cases/refined_gain_marginalized/posterior.npz', allow_pickle=False) as target_npz:
    with np.load(ROOT / 'results/cli_gain_regression_corrected/posterior.npz', allow_pickle=False) as cli_npz:
        for key in ['data', 'noise_real_sd', 'topology_prior', 'gain_sd', 'posterior_weights', 'theta', 'H', 'forward']:
            error = float(np.max(np.abs(target_npz[key] - cli_npz[key])))
            assert error < 1e-12, (key, error)
            comparison[key] = error
        assert float(target_npz['gain_sd']) == float(cli_npz['gain_sd']) == .15
for name in ['gain_marginalized', 'refined_gain_marginalized']:
    with np.load(ROOT / f'results/cases/{name}/posterior.npz', allow_pickle=False) as file:
        assert float(file['gain_sd']) == .15
expected_material_coverage = [.9277777777777778, .8680555555555556]
assert np.array_equal(calibration['marginal_material_90_coverage'], expected_material_coverage)
assert np.array_equal(primary['calibration']['marginal_material_90_coverage'], expected_material_coverage)
assert 'conditional_material_90_coverage' not in calibration
assert 'conditional_material_90_coverage' not in primary['calibration']
topk = target['top_k']
assert topk['outside_adopted_prior_support_mass'] == 0.
assert topk['unmodeled_geometry_mass'] is None
assert topk['exact_image_probability'] == 0.
assert abs(sum(item['posterior_mass'] for item in topk['items']) + topk['remaining_modeled_mass'] - 1) < 1e-12
assert (ROOT / 'results/PHYSICS_LEDGER.jsonl').read_text() == ledger
result = dict(status='PARENT_CORRECTIONS_VERIFIED_FROM_SAVED_RECORDS',
              gain_replay_probabilities=receipt['result']['topology_probabilities'], npz_field_max_errors=comparison,
              gain_default_initial_probability_change=initial_reanalysis['cases']['gain_marginalized']['npz_loader_default_max_probability_change'],
              initial_gain_issue='NPZ loader gain_sd=0 default plus missing gain metadata silently changed likelihood target',
              final_gain_metadata=.15, corrected_cli_new_physics=0,
              material_coverage_label='marginal_material_90_coverage', material_coverage_values_unchanged=expected_material_coverage,
              topk_support_semantics='outside adopted prior is zero by construction; broader physical model mass remains unquantified',
              new_physics_evaluations=0, physics_before=384, physics_after=384, imported_B1_modules=False,
              shared_source_config_results_modified_by_reviewer=False,
              executable=sys.executable, wall_seconds=time.perf_counter() - wall_start,
              process_cpu_seconds=time.process_time() - cpu_start)
(OUT / 'CORRECTION_VERIFICATION_NATIVE.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
print(json.dumps(dict(status=result['status'], new_physics=0, gain_probability_error=0, wall_seconds=result['wall_seconds'], cpu_seconds=result['process_cpu_seconds'])))
