"""Compare frozen inference methods on exactly reconstructed saved calibration.

No physical solver or B1 module is imported. Existing paid responses are reused.
The importance proposal and its 48 samples remain frozen from competitive y0.
"""
from pathlib import Path
import json
import time
import sys
from collections import Counter
wall_start = time.perf_counter()
cpu_start = time.process_time()
import numpy as np
import scipy
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
READS = []


def read_json(path):
    READS.append(str(path.relative_to(ROOT)))
    return json.loads(path.read_text())


def read_npz(path):
    READS.append(str(path.relative_to(ROOT)))
    with np.load(path, allow_pickle=False) as file:
        return {key: file[key].copy() for key in file.files}


def dump(path, value):
    def convert(item):
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, np.generic):
            return item.item()
        raise TypeError(type(item).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False) + '\n')


ledger_path = ROOT / 'results/PHYSICS_LEDGER.jsonl'
ledger_before = ledger_path.read_text()
READS.append('results/PHYSICS_LEDGER.jsonl')
physics_before = sum(json.loads(line)['event'] == 'BEGIN' for line in ledger_before.splitlines())
assert physics_before == 384
config = read_json(ROOT / 'configs/protocol.json')
raw = read_json(ROOT / 'results/CALIBRATION_RAW.json')
manifest = read_json(ROOT / 'results/CALIBRATION_GENERATION_MANIFEST.json')
samples = read_json(ROOT / 'results/importance/SAMPLES.json')
importance_receipt = read_json(ROOT / 'results/importance/RESULT.json')
calibration_receipt = read_json(ROOT / 'results/CALIBRATION_RESULTS.json')
cache_index = read_json(ROOT / 'results/forward_cache/index.json')
cache_entries = [(json.loads(key), filename) for key, filename in cache_index.items()]
banks = {order: read_npz(ROOT / f'results/banks/order{order}.npz') for order in [3, 5, 7]}


def paid_response(h, theta):
    matches = [filename for parameters, filename in cache_entries
               if parameters['H'] == h and parameters['er'] == theta[0] and parameters['ei'] == theta[1]
               and parameters['implementation'] == 'matching' and parameters['nmax'] == config['nmax'] and parameters['k'] == config['wavenumber']]
    assert len(matches) == 1, (h, theta, len(matches))
    return read_npz(ROOT / 'results/forward_cache' / matches[0])['data']


# Replay the exact published seed, not a new generative experiment.
lo, hi = banks[7]['bounds'][:, 0], banks[7]['bounds'][:, 1]
rng = np.random.default_rng(config['seed'] + 400)
parent_theta = np.asarray([lo + rng.uniform(size=2) * (hi - lo) for _ in range(24)])
assert np.array_equal(parent_theta, np.asarray([parent['theta'] for parent in manifest['parents']]))
assert np.array_equal(parent_theta, np.asarray([parent['theta'] for parent in raw]))
sigma = .65 * np.sqrt(np.mean(banks[7]['F'][banks[7]['H'] == 0]**2))
data, truth, parent_ids, replicates, recorded_probability = [], [], [], [], []
for parent, theta in zip(raw, parent_theta):
    response = paid_response(parent['H'], theta)
    for rep in range(30):
        data.append(response + sigma * rng.normal(size=len(response)))
        truth.append(parent['H'])
        parent_ids.append(parent['id'])
        replicates.append(rep)
        recorded_probability.append(parent['replicates'][rep]['probabilities'])
data = np.asarray(data)
truth = np.asarray(truth)
parent_ids = np.asarray(parent_ids)
replicates = np.asarray(replicates)
recorded_probability = np.asarray(recorded_probability)


def likelihood_matrix(forward):
    # Saved real vector responses; common fixed sigma, no gain nuisance here.
    residual = data[:, None, :] - forward[None, :, :]
    q = np.einsum('nsm,nsm->ns', residual, residual) / sigma**2
    logl = -.5 * (q + forward.shape[1] * np.log(2 * np.pi * sigma**2))
    return logl, q


def integrate_bank(bank):
    logl, q = likelihood_matrix(bank['F'])
    log_evidence = np.stack([logsumexp(logl[:, bank['H'] == h] + np.log(bank['weights'][bank['H'] == h]), axis=1) for h in range(3)], axis=1)
    log_joint = log_evidence + np.log(np.repeat(1 / 3, 3))
    logp = log_joint - logsumexp(log_joint, axis=1)[:, None]
    return np.exp(logp), logp, logl, q


method_probabilities = {}
method_log_probabilities = {}
grid_likelihoods = {}
for order in [3, 5, 7]:
    p, logp, logl, q = integrate_bank(banks[order])
    name = f'gauss{order}'
    method_probabilities[name] = p
    method_log_probabilities[name] = logp
    grid_likelihoods[order] = (logl, q)
identity_error = float(np.max(np.abs(method_probabilities['gauss7'] - recorded_probability)))
assert identity_error < 1e-12

# Fixed proposal bank: each branch estimator is mean L/q_z, unnormalized.
importance_H = np.asarray([sample['H'] for sample in samples])
importance_theta = np.asarray([sample['theta'] for sample in samples])
importance_q = np.asarray([sample['proposal_density_z'] for sample in samples])
assert np.all(importance_q >= .5)
assert Counter(importance_H) == {0: 16, 1: 16, 2: 16}
importance_F = np.asarray([paid_response(sample['H'], sample['theta']) for sample in samples])
is_logl, _ = likelihood_matrix(importance_F)
is_evidence = np.stack([logsumexp(is_logl[:, importance_H == h] - np.log(importance_q[importance_H == h]), axis=1) - np.log(16) for h in range(3)], axis=1)
is_log_joint = is_evidence + np.log(1 / 3)
is_logp = is_log_joint - logsumexp(is_log_joint, axis=1)[:, None]
method_probabilities['frozen_importance48'] = np.exp(is_logp)
method_log_probabilities['frozen_importance48'] = is_logp

logl7, q7 = grid_likelihoods[7]
minimum_q = np.stack([q7[:, banks[7]['H'] == h].min(axis=1) for h in range(3)], axis=1)
residual_logits = -.5 * minimum_q + np.log(1 / 3)
residual_logp = residual_logits - logsumexp(residual_logits, axis=1)[:, None]
method_probabilities['grid7_min_residual_softmax_HEURISTIC'] = np.exp(residual_logp)
method_log_probabilities['grid7_min_residual_softmax_HEURISTIC'] = residual_logp
map_H = banks[7]['H'][np.argmax(logl7 + np.log(1 / 3), axis=1)]
hard_p = np.eye(3)[map_H]
method_probabilities['grid7_hard_MAP_ONE_HOT_DECISION'] = hard_p

reference = method_probabilities['gauss7']
reference_log_bf = method_log_probabilities['gauss7'] - method_log_probabilities['gauss7'][:, :1]
correct_one_hot = np.eye(3)[truth]
method_metrics = {}
per_datum = {}
for name, p in method_probabilities.items():
    brier = np.sum((p - correct_one_hot)**2, axis=1)
    error = np.max(np.abs(p - reference), axis=1)
    total_variation = .5 * np.sum(np.abs(p - reference), axis=1)
    accuracy = np.argmax(p, axis=1) == truth
    passes = error <= config['probability_ref_tolerance']
    if name in method_log_probabilities:
        logscore = -method_log_probabilities[name][np.arange(720), truth]
        infinite = np.zeros(720, bool)
    else:
        infinite = ~accuracy
        # Mean log score is +infinity if any incorrect one-hot prediction occurs.
        # Null represents infinite entries in strict JSON, with explicit flags.
        logscore = np.zeros(720)
    metrics = np.stack([brier, logscore, error, passes.astype(float), accuracy.astype(float)], axis=1)
    method_metrics[name] = metrics
    per_datum[name] = dict(probabilities=p, Brier=brier, log_score=[None if flag else float(value) for flag, value in zip(infinite, logscore)],
                           log_score_is_positive_infinity=infinite, max_probability_error_vs_gauss7=error,
                           total_variation_vs_gauss7=total_variation, within_probability_tolerance=passes,
                           MAP_classification_correct=accuracy)
    if name in ['gauss3', 'gauss5', 'gauss7', 'frozen_importance48']:
        method_log_bf = method_log_probabilities[name] - method_log_probabilities[name][:, :1]
        bf_error = np.max(np.abs(method_log_bf - reference_log_bf), axis=1)
        per_datum[name]['max_log_BF_difference_vs_gauss7'] = bf_error
        per_datum[name]['within_both_reference_tolerances'] = passes & (bf_error <= config['log_bayes_factor_ref_tolerance'])

bootstrap_rng = np.random.default_rng(config['seed'] + 401)
resamples = bootstrap_rng.integers(0, 8, size=(2000, 3, 8))
parent_reference = method_metrics['gauss7'].reshape(24, 30, 5).mean(axis=1)
summary = {}
metric_names = ['Brier_multiclass', 'log_score', 'mean_max_probability_error_vs_gauss7', 'fraction_within_probability_tolerance', 'MAP_accuracy']
costs = dict(gauss3=27, gauss5=75, gauss7=147, frozen_importance48=75,
             grid7_min_residual_softmax_HEURISTIC=147, grid7_hard_MAP_ONE_HOT_DECISION=147)
for name, metrics in method_metrics.items():
    parent_means = metrics.reshape(24, 30, 5).mean(axis=1)
    grouped = parent_means.reshape(3, 8, 5)
    boot = np.asarray([np.mean([grouped[h, indices[h]].mean(axis=0) for h in range(3)], axis=0) for indices in resamples])
    paired = parent_means - parent_reference
    paired_grouped = paired.reshape(3, 8, 5)
    paired_boot = np.asarray([np.mean([paired_grouped[h, indices[h]].mean(axis=0) for h in range(3)], axis=0) for indices in resamples])
    estimates = {}
    for index, metric_name in enumerate(metric_names):
        if name == 'grid7_hard_MAP_ONE_HOT_DECISION' and metric_name == 'log_score' and np.any(per_datum[name]['log_score_is_positive_infinity']):
            estimates[metric_name] = dict(estimate=None, status='POSITIVE_INFINITY_ON_ANY_WRONG_ONE_HOT_PREDICTION',
                                          infinite_datums=int(np.sum(per_datum[name]['log_score_is_positive_infinity'])),
                                          parent_stratified_bootstrap95=None, difference_from_gauss7=None,
                                          probability_floor_applied=False)
        else:
            estimates[metric_name] = dict(estimate=float(parent_means[:, index].mean()),
                                          parent_stratified_bootstrap95=np.quantile(boot[:, index], [.025, .975]),
                                          paired_difference_from_gauss7=float(paired[:, index].mean()),
                                          paired_difference_parent_stratified_bootstrap95=np.quantile(paired_boot[:, index], [.025, .975]))
    error = metrics[:, 2]
    worst = int(np.argmax(error))
    summary[name] = dict(cold_nominal_material_frequency_states=costs[name], actual_new_physics_for_this_reanalysis=0,
                          metrics=estimates, empirical_max_probability_error=float(error.max()),
                          empirical_probability_error_quantiles_50_90_95_99=np.quantile(error, [.5, .9, .95, .99]),
                          per_topology_mean_probability_error=[float(error[truth == h].mean()) for h in range(3)],
                          worst_datum=dict(parent_id=int(parent_ids[worst]), replicate=int(replicates[worst]), truth_H=int(truth[worst]),
                                           probabilities=method_probabilities[name][worst], gauss7_reference=reference[worst]),
                          posterior_status='ONE_HOT_POINT_DECISION_NOT_POSTERIOR' if name == 'grid7_hard_MAP_ONE_HOT_DECISION'
                          else 'MINIMUM_RESIDUAL_HEURISTIC_NOT_POSTERIOR' if 'HEURISTIC' in name else 'CONTINUOUS_PRIOR_EVENT_MASS_NUMERICAL_ESTIMATE')
    if 'max_log_BF_difference_vs_gauss7' in per_datum[name]:
        bf_error = per_datum[name]['max_log_BF_difference_vs_gauss7']
        both_passes = per_datum[name]['within_both_reference_tolerances']
        parent_bf = bf_error.reshape(3, 8, 30).mean(axis=2)
        parent_pass = both_passes.reshape(3, 8, 30).mean(axis=2)
        boot_bf = [np.mean([parent_bf[h, indices[h]].mean() for h in range(3)]) for indices in resamples]
        boot_pass = [np.mean([parent_pass[h, indices[h]].mean() for h in range(3)]) for indices in resamples]
        summary[name]['reference_gate_vs_gauss7'] = dict(max_log_BF_difference=float(bf_error.max()),
                                                        mean_log_BF_difference=float(bf_error.mean()),
                                                        mean_log_BF_difference_parent_bootstrap95=np.quantile(boot_bf, [.025, .975]),
                                                        fraction_within_both_reference_tolerances=float(both_passes.mean()),
                                                        both_reference_tolerance_fraction_parent_bootstrap95=np.quantile(boot_pass, [.025, .975]),
                                                        probability_tolerance=config['probability_ref_tolerance'], log_BF_tolerance=config['log_bayes_factor_ref_tolerance'],
                                                        certificate=False)

parent_is = method_metrics['frozen_importance48'].reshape(24, 30, 5).mean(axis=1)
parent_g5 = method_metrics['gauss5'].reshape(24, 30, 5).mean(axis=1)
paired_is_minus_g5 = (parent_is - parent_g5).reshape(3, 8, 5)
paired_is_g5_boot = np.asarray([np.mean([paired_is_minus_g5[h, indices[h]].mean(axis=0) for h in range(3)], axis=0) for indices in resamples])
equal_cold_cost_comparison = dict(methods=['frozen_importance48', 'gauss5'], nominal_cold_states_each=75,
                                  sign='frozen_importance48 minus gauss5; positive Brier/log-score differences favor gauss5 on this panel',
                                  frozen_importance_not_adapted_per_datum=True,
                                  metrics={metric: dict(paired_mean_difference=float(paired_is_minus_g5[:, :, index].mean()),
                                                        parent_stratified_bootstrap95=np.quantile(paired_is_g5_boot[:, index], [.025, .975]))
                                           for index, metric in enumerate(metric_names)})

assert abs(summary['gauss7']['metrics']['Brier_multiclass']['estimate'] - calibration_receipt['Brier_multiclass']) < 1e-12
assert abs(summary['gauss7']['metrics']['log_score']['estimate'] - calibration_receipt['log_score']) < 1e-12
assert ledger_path.read_text() == ledger_before
np.savez_compressed(OUT / 'DERIVED_CALIBRATION_OBSERVATIONS_NATIVE.npz',
                    data=data, truth_H=truth, parent_id=parent_ids, replicate=replicates,
                    parent_theta=parent_theta, noise_real_sd=sigma, topology_prior=np.repeat(1 / 3, 3),
                    gauss7_probabilities=reference)
scope = dict(target='Adopted 3-template continuous uniform material prior, 32 real Gaussian components, no gain nuisance',
              calibration_design='24 previously generated continuous material parents, 8 per H, 30 previously recorded noise realizations each',
              reconstruction='Published seed replay; exact parent manifest and all 720 Gauss7 probability records verified before comparison',
              original_probability_identity_max_error=identity_error,
              reference='Original paid order7 full-Mie Gauss quadrature; exposed development reference, not certified truth',
              frozen_importance_proposal='Single competitive y0 proposal fitted only to its order3 posterior; samples and q_z frozen for all 720 other data',
              importance_adapted_per_datum=False, importance_redrawn_for_calibration=False,
              cold_cost_semantics='Hypothetical standalone method setup states from existing paid banks, not additional charges. Gauss5 3*25=75; frozen IS27+48=75; Gauss7 reference 3*49=147. Common finite solver audit and reference validation costs are separate.',
              gauss5_actual_initial_pipeline_new_calls=72, gauss5_initial_reused_center_states=3,
              bootstrap='2,000 parent-stratified resamples, 8 parents per true H, paired comparisons preserve shared data',
              method_superiority_scope='Within this already exposed prior/noise calibration panel only; no general sampler superiority or real-world calibration claim',
              hard_map_log_score='Exact one-hot decision: +infinity on wrong predictions. Strict JSON uses null plus explicit positive-infinity flags; no floor.',
              no_new_experiment=True, no_uniform_error_certificate=True, locked_or_historical_arrays_read=False,
              imported_B1_modules=False, neural_training=False, installs=False)
result = dict(schema='B1-derived-frozen-method-comparison-native-v1', scope=scope,
              datum_metadata=dict(parent_id=parent_ids, replicate=replicates, truth_H=truth, noise_real_sd=float(sigma)),
              summaries=summary, equal_nominal_cold_cost_paired_comparison=equal_cold_cost_comparison,
              per_datum=per_datum, new_physics_evaluations=0,
              physics_before=384, physics_after=384, independent_material_clusters=24,
              executable=sys.executable, numpy_version=np.__version__, scipy_version=scipy.__version__,
              read_paths=sorted(set(READS)), wall_seconds=time.perf_counter() - wall_start,
              process_cpu_seconds=time.process_time() - cpu_start)
dump(OUT / 'DERIVED_METHOD_COMPARISON_NATIVE.json', result)
receipt = {key: result[key] for key in ['schema', 'new_physics_evaluations', 'physics_before', 'physics_after', 'independent_material_clusters', 'executable', 'numpy_version', 'scipy_version', 'read_paths', 'wall_seconds', 'process_cpu_seconds']}
receipt['scope'] = scope
receipt['outputs'] = ['DERIVED_METHOD_COMPARISON_NATIVE.json', 'DERIVED_CALIBRATION_OBSERVATIONS_NATIVE.npz']
run_receipts_path = OUT / 'DERIVED_METHOD_COMPARISON_RUNS_NATIVE.json'
if run_receipts_path.exists():
    run_receipts = json.loads(run_receipts_path.read_text())
else:
    previous_receipt_path = OUT / 'DERIVED_METHOD_COMPARISON_RECEIPT_NATIVE.json'
    run_receipts = []
    if previous_receipt_path.exists():
        previous = json.loads(previous_receipt_path.read_text())
        run_receipts.append(dict(run=1, wall_seconds=previous['wall_seconds'], process_cpu_seconds=previous['process_cpu_seconds'], new_physics_evaluations=0,
                                 purpose='Initial probability, Brier, log-score, and clustered uncertainty comparison'))
run_receipts.append(dict(run=len(run_receipts) + 1, wall_seconds=receipt['wall_seconds'], process_cpu_seconds=receipt['process_cpu_seconds'], new_physics_evaluations=0,
                         purpose='Additional logBF reference gates and direct equal-cold-cost paired comparison'))
dump(run_receipts_path, run_receipts)
receipt['total_comparison_wall_seconds'] = sum(item['wall_seconds'] for item in run_receipts)
receipt['total_comparison_process_cpu_seconds'] = sum(item['process_cpu_seconds'] for item in run_receipts)
receipt['comparison_runs'] = len(run_receipts)
dump(OUT / 'DERIVED_METHOD_COMPARISON_RECEIPT_NATIVE.json', receipt)
print(json.dumps(dict(new_physics=0, physics_budget='384/384', wall_seconds=result['wall_seconds'], cpu_seconds=result['process_cpu_seconds'],
                      identity_error=identity_error,
                      methods={name: dict(Brier=record['metrics']['Brier_multiclass']['estimate'], log_score=record['metrics']['log_score']['estimate'],
                                          mean_probability_error=record['metrics']['mean_max_probability_error_vs_gauss7']['estimate'],
                                          max_probability_error=record['empirical_max_probability_error'],
                                          pass_fraction=record['metrics']['fraction_within_probability_tolerance']['estimate']) for name, record in summary.items()})))
