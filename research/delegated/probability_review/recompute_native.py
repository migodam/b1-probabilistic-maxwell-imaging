"""Read-only reanalysis of already paid B1 arrays; imports no B1 solver code."""
from pathlib import Path
import json
import sys
import time
import resource
from collections import Counter

wall_start = time.perf_counter()
cpu_start = time.process_time()
import numpy as np
import scipy
from scipy.special import logsumexp
from scipy.stats import beta, chi2
from numpy.polynomial.legendre import leggauss, legvander2d

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
READS = []


def read_json(path):
    READS.append(str(path.relative_to(ROOT)))
    return json.loads(path.read_text())


def read_npz(path):
    READS.append(str(path.relative_to(ROOT)))
    with np.load(path, allow_pickle=False) as arr:
        return {key: arr[key].copy() for key in arr.files}


def dump(path, value):
    def convert(item):
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, np.generic):
            return item.item()
        raise TypeError(type(item).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False) + '\n')


def likelihood(forward, data, sigma, gain_sd=0.):
    # Independently expand Sherman-Morrison and matrix determinant lemma.
    residual = data - forward
    field_norm2 = np.einsum('ij,ij->i', forward, forward)
    residual_dot = np.einsum('ij,ij->i', residual, forward)
    residual_norm2 = np.einsum('ij,ij->i', residual, residual)
    correction = gain_sd**2 * residual_dot**2 / (sigma**2 + gain_sd**2 * field_norm2)
    q = (residual_norm2 - correction) / sigma**2
    logdet = forward.shape[1] * np.log(sigma**2) + np.log1p(gain_sd**2 * field_norm2 / sigma**2)
    return -.5 * (q + logdet + forward.shape[1] * np.log(2 * np.pi)), q


def quantiles(nodes, posterior_weights, prior_weights, bounds):
    values = np.unique(nodes)
    mass = np.asarray([posterior_weights[nodes == value].sum() for value in values])
    mass = mass / mass.sum()
    prior_mass = np.asarray([prior_weights[nodes == value].sum() for value in values])
    prior_mass = prior_mass / prior_mass.sum()
    edges = bounds[0] + (bounds[1] - bounds[0]) * np.r_[0., prior_mass.cumsum()]
    edges[-1] = bounds[1]
    return np.interp([.05, .5, .95], np.r_[0., mass.cumsum()], edges)


def posterior(bank, data, sigma, topology_prior, gain_sd=0.):
    logl, q = likelihood(bank['F'], data, sigma, gain_sd)
    H, prior_weights = bank['H'], bank['weights']
    evidence = np.asarray([logsumexp(logl[H == h] + np.log(prior_weights[H == h])) for h in range(3)])
    total = logsumexp(np.log(topology_prior) + evidence)
    p = np.exp(np.log(topology_prior) + evidence - total)
    W = np.exp(np.log(np.asarray(topology_prior)[H]) + logl + np.log(prior_weights) - total)
    qs, means, covariances, representatives = [], [], [], []
    for h in range(3):
        ids = np.flatnonzero(H == h)
        conditional = np.exp(logl[ids] + np.log(prior_weights[ids]) - evidence[h])
        mean = conditional @ bank['theta'][ids]
        centered = bank['theta'][ids] - mean
        covariance = (centered.T * conditional) @ centered
        means.append(mean)
        covariances.append(covariance)
        qs.append([quantiles(bank['theta'][ids, j], conditional, prior_weights[ids], bank['bounds'][j]) for j in range(2)])
        representatives.append(int(ids[np.argmax(logl[ids])]))
    return dict(p=p, logz=evidence, logbf=evidence - evidence[0], W=W,
                conditional_quantiles=np.asarray(qs), means=np.asarray(means),
                covariances=np.asarray(covariances), representatives=representatives,
                ppc_tail=float(W @ chi2.sf(q, bank['F'].shape[1])))


def polynomial_bank(coefficients, order, config, bounds):
    x, w = leggauss(order)
    unit = (x + 1) / 2
    theta = np.asarray([(a, b) for a in unit for b in unit])
    theta = bounds[:, 0] + theta * (bounds[:, 1] - bounds[:, 0])
    weights = np.asarray([a * b / 4 for a in w for b in w])
    degree = int(np.sqrt(coefficients.shape[1])) - 1
    standardized = 2 * (theta - bounds[:, 0]) / (bounds[:, 1] - bounds[:, 0]) - 1
    V = legvander2d(standardized[:, 0], standardized[:, 1], [degree, degree]).reshape(len(theta), -1)
    return dict(H=np.repeat(np.arange(3), len(theta)), theta=np.tile(theta, (3, 1)),
                weights=np.tile(weights, 3), F=np.concatenate([V @ coefficients[h] for h in range(3)]), bounds=bounds)


def reference_difference(coarse, fine, config):
    pd = float(np.max(np.abs(coarse['p'] - fine['p'])))
    bd = float(np.max(np.abs(coarse['logbf'] - fine['logbf'])))
    return dict(max_probability_difference=pd, max_log_BF_difference=bd,
                numerical_pass=pd <= config['probability_ref_tolerance'] and bd <= config['log_bayes_factor_ref_tolerance'])


ledger_path = ROOT / 'results/PHYSICS_LEDGER.jsonl'
ledger_before = ledger_path.read_text()
READS.append('results/PHYSICS_LEDGER.jsonl')
rows = [json.loads(line) for line in ledger_before.splitlines()]
begins = [row for row in rows if row['event'] == 'BEGIN']
ends = [row for row in rows if row['event'] == 'END']
config = read_json(ROOT / 'configs/protocol.json')
primary = read_json(ROOT / 'results/PRIMARY_SUMMARY.json')
refined = read_json(ROOT / 'results/refinement/SUMMARY.json')
checks = read_json(ROOT / 'results/refinement/FULL_FORWARD_CHECKS.json')
contract = read_json(ROOT / 'results/BUDGET_CONTRACT.json')
resources = read_json(ROOT / 'results/RESOURCE_SUMMARY.json')
job_rows = [json.loads(line) for line in (ROOT / 'results/JOB_LEDGER.jsonl').read_text().splitlines()]
READS.append('results/JOB_LEDGER.jsonl')
banks = {order: read_npz(ROOT / f'results/banks/order{order}.npz') for order in [3, 5, 7]}
polynomials = read_npz(ROOT / 'results/refinement/FORWARD_POLYNOMIAL.npz')
poly517 = polynomial_bank(polynomials['coeff5'], 17, config, polynomials['bounds'])
poly717 = polynomial_bank(polynomials['coeff7'], 17, config, polynomials['bounds'])
poly725 = polynomial_bank(polynomials['coeff7'], 25, config, polynomials['bounds'])

accounting = dict(new_evaluations=len(begins), cap=contract['cap'], remaining=contract['cap'] - len(begins),
                  contiguous_unique_ids=[row['id'] for row in begins] == list(range(1, len(begins) + 1)),
                  all_finished_once=Counter(row['id'] for row in ends) == Counter(row['id'] for row in begins),
                  failures=sum(row['status'] == 'FAIL' for row in ends),
                  by_kind=dict(sorted(Counter(row['kind'] for row in begins).items())),
                  by_event=dict(sorted(Counter(row['event'] for row in rows).items())),
                  physics_wall_seconds=sum(row['wall_seconds'] for row in ends),
                  job_physics=sum(row['physics_evaluations'] for row in job_rows),
                  job_wall_seconds=sum(row['wall_seconds'] for row in job_rows),
                  job_cpu_seconds=sum(row['cpu_seconds'] for row in job_rows),
                  ledger_high_water=contract['high_water_reservations'], historical_declared_only=contract['historical_2d_calls'],
                  historical_arrays_read=False)
assert accounting['new_evaluations'] == accounting['cap'] == accounting['ledger_high_water'] == 384
assert accounting['contiguous_unique_ids'] and accounting['all_finished_once'] and accounting['failures'] == 0
assert accounting['by_kind'] == resources['by_kind'] and accounting['job_physics'] == 384

normalizations = {}
for order, bank in banks.items():
    normalizations[str(order)] = dict(prior_weight_sums=[float(bank['weights'][bank['H'] == h].sum()) for h in range(3)],
                                   prior_means=[(bank['weights'][bank['H'] == h] @ bank['theta'][bank['H'] == h]).tolist() for h in range(3)])

case_results = {}
dense_covariance_errors = []
for name, record in primary['cases'].items():
    saved = read_npz(ROOT / f'results/cases/{name}/posterior.npz')
    gain_sd = .15 if name == 'gain_marginalized' else 0.
    data, sigma, prior = saved['data'], float(saved['noise_real_sd']), saved['topology_prior']
    r5, r7 = [posterior(banks[order], data, sigma, prior, gain_sd) for order in [5, 7]]
    r517, r717, r725 = [posterior(bank, data, sigma, prior, gain_sd) for bank in [poly517, poly717, poly725]]
    saved_result = record['result']
    saved_refined = refined['cases'][name]['result']
    remaining = float(np.sort(r7['p'])[0])
    case_results[name] = dict(original_probabilities=r7['p'], refined_probabilities=r725['p'],
                             original_saved_probability_max_error=float(np.max(np.abs(r7['p'] - saved_result['topology_probabilities']))),
                             original_saved_log_evidence_max_error=float(np.max(np.abs(r7['logz'] - saved_result['log_evidence']))),
                             original_saved_posterior_weight_max_error=float(np.max(np.abs(r7['W'] - saved['posterior_weights']))),
                             refined_saved_probability_max_error=float(np.max(np.abs(r725['p'] - saved_refined['topology_probabilities']))),
                             original_5_vs_7=reference_difference(r5, r7, config),
                             poly_integration_17_vs_25=reference_difference(r717, r725, config),
                             poly_model_5_vs_7=reference_difference(r517, r717, config),
                             top2_remaining_original=remaining, top2_remaining_refined=float(np.sort(r725['p'])[0]),
                             original_conditional_quantiles=r7['conditional_quantiles'], refined_conditional_quantiles=r725['conditional_quantiles'],
                             max_conditional_quantile_difference_17_vs_25=np.max(np.abs(r717['conditional_quantiles'] - r725['conditional_quantiles']), axis=(0, 2)),
                             posterior_weight_sum=float(r7['W'].sum()), original_ppc_tail=r7['ppc_tail'],
                             refined_ppc_tail=r725['ppc_tail'], gain_sd_used=gain_sd,
                             saved_npz_contains_gain_sd='gain_sd' in saved)
    assert case_results[name]['original_saved_probability_max_error'] < 1e-12
    assert case_results[name]['refined_saved_probability_max_error'] < 1e-12
    if gain_sd:
        for index in [0, 36, 81, 146]:
            f = banks[7]['F'][index]
            covariance = sigma**2 * np.eye(len(f)) + gain_sd**2 * np.outer(f, f)
            residual = data - f
            # Stable dense covariance reference, independently factorized.
            cholesky = np.linalg.cholesky(covariance)
            whitened = np.linalg.solve(cholesky, residual)
            dense = -.5 * (whitened @ whitened + 2 * np.log(np.diag(cholesky)).sum() + len(f) * np.log(2 * np.pi))
            fast = likelihood(f[None], data, sigma, gain_sd)[0][0]
            assert np.isfinite(dense) and np.isfinite(fast)
            dense_covariance_errors.append(float(abs(dense - fast)))
        loader_default = posterior(poly725, data, sigma, prior, 0.)
        case_results[name]['npz_loader_default_gain_sd'] = 0.
        case_results[name]['npz_loader_default_refined_probabilities'] = loader_default['p']
        case_results[name]['npz_loader_default_max_probability_change'] = float(np.max(np.abs(loader_default['p'] - r725['p'])))

interpolation = dict(check_count=len(checks['checks']), new_paid_calls=checks['actual_new_calls'],
                     max_component_error=float(max(np.max(np.abs(np.asarray(item['F']) - item['F_interp7'])) for item in checks['full_outputs'])),
                     max_relative_field_error=float(max(np.linalg.norm(np.asarray(item['F']) - item['F_interp7']) / np.linalg.norm(item['F']) for item in checks['full_outputs'])),
                     max_independent_treams_reported_error=max(item['relative_field_error'] for item in checks['independent_treams']),
                     independent_treams_outputs_saved=False, uniform_error_bound_available=False)

importance_saved = read_json(ROOT / 'results/importance/RESULT.json')
samples = read_json(ROOT / 'results/importance/SAMPLES.json')
cache_index = read_json(ROOT / 'results/forward_cache/index.json')
cache_entries = [(json.loads(key), filename) for key, filename in cache_index.items()]


def paid_response(h, theta):
    matches = [(parameters, filename) for parameters, filename in cache_entries
               if parameters['H'] == h and parameters['er'] == theta[0] and parameters['ei'] == theta[1]
               and parameters['implementation'] == 'matching' and parameters['nmax'] == config['nmax'] and parameters['k'] == config['wavenumber']]
    assert len(matches) == 1, (h, theta, len(matches))
    return read_npz(ROOT / 'results/forward_cache' / matches[0][1])['data']


competitive = read_npz(ROOT / 'results/cases/competitive/posterior.npz')
importance_checks = []
importance_logz = []
for h in range(3):
    diagnostics = importance_saved['diagnostics'][h]
    beta_a, beta_b = np.asarray(diagnostics['beta_a']), np.asarray(diagnostics['beta_b'])
    records = [record for record in samples if record['H'] == h]
    independent_weights = []
    density_errors, logweight_errors = [], []
    for record in records:
        theta = np.asarray(record['theta'])
        z = (theta - banks[7]['bounds'][:, 0]) / np.diff(banks[7]['bounds'], axis=1)[:, 0]
        density = .5 + .5 * np.prod(beta.pdf(z, beta_a, beta_b))
        full = paid_response(h, theta)
        logl = likelihood(full[None], competitive['data'], float(competitive['noise_real_sd']))[0][0]
        logweight = float(logl - np.log(density))
        density_errors.append(abs(density - record['proposal_density_z']))
        logweight_errors.append(abs(logweight - record['logweight']))
        independent_weights.append(logweight)
    independent_weights = np.asarray(independent_weights)
    logz = logsumexp(independent_weights) - np.log(len(independent_weights))
    importance_logz.append(float(logz))
    normalized = np.exp(independent_weights - logsumexp(independent_weights))
    importance_checks.append(dict(H=h, points=len(records), max_density_error=float(max(density_errors)),
                                  max_logweight_error=float(max(logweight_errors)), ESS=float(1 / (normalized @ normalized))))
importance_logz = np.asarray(importance_logz)
importance_p = np.exp(np.log(competitive['topology_prior']) + importance_logz - logsumexp(np.log(competitive['topology_prior']) + importance_logz))
importance = dict(probabilities=importance_p, per_branch=importance_checks,
                  max_saved_probability_error=float(np.max(np.abs(importance_p - importance_saved['topology_probabilities']))),
                  max_probability_error_vs_original_order7=float(np.max(np.abs(importance_p - case_results['competitive']['original_probabilities']))),
                  max_probability_error_vs_refined=float(np.max(np.abs(importance_p - case_results['competitive']['refined_probabilities']))),
                  probability_tolerance=config['probability_ref_tolerance'], meets_declared_tolerance=False,
                  density_semantics='q_z; uniform theta prior is 1/volume and q_theta=q_z/volume, so ratio is 1/q_z',
                  total_paid_setup_plus_samples=27 + 48, finite_ratio_bias=True, bootstrap_missed_mode_certification=False)
assert importance['max_saved_probability_error'] < 1e-12

linear = read_npz(ROOT / 'results/LINEAR_MODELS.npz')
mechanism_saved = read_json(ROOT / 'results/MODEL_CORRECTION_MECHANISM.json')
mechanism = {}
for name, record in primary['cases'].items():
    if name not in mechanism_saved:
        continue
    obs = read_npz(ROOT / f'results/cases/{name}/posterior.npz')
    sigma, data = float(obs['noise_real_sd']), obs['data']
    logcheap = likelihood(linear['F_dipole_linear'], data, sigma)[0]
    logfull = likelihood(banks[7]['F'], data, sigma)[0]
    delta = logcheap - logfull
    ranks, identities, remainder = [], [], []
    predicted = []
    for h in range(3):
        ids = banks[7]['H'] == h
        s = banks[7]['theta'][ids] - linear['theta0']
        Jfull, Jcheap = linear['J_full'][h] / sigma, linear['J_dipole'][h] / sigma
        residual = (linear['centers'][h] - data) / sigma
        U, singular, Vh = np.linalg.svd(Jfull - Jcheap, full_matrices=False)
        ranks.append(int(np.linalg.matrix_rank(Jfull - Jcheap)))
        Y, N = U * singular, Vh
        coordinates = s @ N.T
        signature = np.sum((residual @ Y + s @ Jcheap.T @ Y) * coordinates, axis=1) + .5 * np.sum((coordinates @ Y.T)**2, axis=1)
        direct = .5 * np.sum((residual + s @ Jfull.T)**2, axis=1) - .5 * np.sum((residual + s @ Jcheap.T)**2, axis=1)
        identities.append(float(np.max(np.abs(signature - direct))))
        remainder.extend((delta[ids] - signature).tolist())
        cheap_evidence = logsumexp(np.log(banks[7]['weights'][ids]) + logcheap[ids])
        log_q = np.log(banks[7]['weights'][ids]) + logcheap[ids] - cheap_evidence
        predicted.append(cheap_evidence + logsumexp(log_q - delta[ids]))
    direct_evidence = posterior(banks[7], data, sigma, obs['topology_prior'])['logz']
    mechanism[name] = dict(affine_ranks=ranks, affine_identity_max_error=max(identities),
                          nonlinear_remainder_abs_max=float(np.max(np.abs(remainder))),
                          exact_reweight_evidence_max_error=float(np.max(np.abs(predicted - direct_evidence))),
                          native_A17_exchange=False, paid_full_target_used_for_reweight=True)

manifest = read_json(ROOT / 'results/CALIBRATION_GENERATION_MANIFEST.json')
raw = read_json(ROOT / 'results/CALIBRATION_RAW.json')
calibration_saved = read_json(ROOT / 'results/CALIBRATION_RESULTS.json')
rng = np.random.default_rng(config['seed'] + 400)
lo, hi = banks[7]['bounds'][:, 0], banks[7]['bounds'][:, 1]
reconstructed_theta = np.asarray([lo + rng.uniform(size=2) * (hi - lo) for _ in range(24)])
manifest_theta = np.asarray([item['theta'] for item in manifest['parents']])
assert np.array_equal(reconstructed_theta, manifest_theta)
assert np.array_equal(reconstructed_theta, np.asarray([item['theta'] for item in raw]))
calibration_sigma = .65 * np.sqrt(np.mean(banks[7]['F'][banks[7]['H'] == 0]**2))
parent_metrics = []
probability_identity_errors = []
saved_metric_errors = []
all_probabilities = []
all_truth = []
interval_records = []
for item, theta in zip(raw, reconstructed_theta):
    h = item['H']
    response = paid_response(h, theta)
    metric_rows = []
    for rep in range(30):
        data = response + calibration_sigma * rng.normal(size=len(response))
        result = posterior(banks[7], data, calibration_sigma, np.repeat(1 / 3, 3))
        old = item['replicates'][rep]
        probability_identity_errors.append(float(np.max(np.abs(result['p'] - old['probabilities']))))
        ordering = np.argsort(-result['p'])
        credible = ordering[:np.searchsorted(np.cumsum(result['p'][ordering]), .9) + 1]
        marginal = np.asarray([quantiles(banks[7]['theta'][:, j], result['W'], banks[7]['weights'], banks[7]['bounds'][j]) for j in range(2)])
        conditional = result['conditional_quantiles'][h]
        marginal_covered = (marginal[:, 0] <= theta) & (theta <= marginal[:, 2])
        conditional_covered = (conditional[:, 0] <= theta) & (theta <= conditional[:, 2])
        assert list(marginal_covered) == old['material_90_contains_truth']
        brier = float(np.sum((result['p'] - np.eye(3)[h])**2))
        logscore = float(-np.log(result['p'][h]))
        saved_metric_errors.extend([abs(brier - old['brier']), abs(logscore - old['log_score'])])
        top2_covered = float(h in ordering[:2])
        top2_mass = float(result['p'][ordering[:2]].sum())
        metric_rows.append([brier, logscore, float(h in credible), top2_covered, top2_mass,
                            top2_covered - top2_mass, *marginal_covered.astype(float),
                            *conditional_covered.astype(float), float(result['p'][credible].sum())])
        all_probabilities.append(result['p'])
        all_truth.append(h)
        interval_records.append(dict(parent_id=item['id'], rep=rep, truth_H=h, marginal_quantiles=marginal,
                                     truth_H_conditional_quantiles=conditional, marginal_contains_truth=marginal_covered,
                                     truth_H_conditional_contains_truth=conditional_covered))
    parent_metrics.append(np.mean(metric_rows, axis=0))
assert max(probability_identity_errors) < 1e-12
assert max(saved_metric_errors) < 1e-12
parent_metrics = np.asarray(parent_metrics)
resample_rng = np.random.default_rng(config['seed'] + 401)
resamples = resample_rng.integers(0, 8, size=(2000, 3, 8))
grouped = parent_metrics.reshape(3, 8, -1)
bootstrap = np.asarray([np.mean([grouped[h, indices[h]].mean(axis=0) for h in range(3)], axis=0) for indices in resamples])
metric_names = ['brier', 'log_score', 'topology_credible90_coverage', 'top2_coverage', 'top2_reported_mass',
                'top2_coverage_minus_mass', 'marginal_real90_coverage', 'marginal_imag90_coverage',
                'truth_H_conditional_real90_coverage', 'truth_H_conditional_imag90_coverage', 'topology_credible90_mean_reported_mass']
cluster_metrics = {name: dict(estimate=float(parent_metrics[:, index].mean()),
                            parent_stratified_bootstrap95=np.quantile(bootstrap[:, index], [.025, .975]).tolist())
                   for index, name in enumerate(metric_names)}
assert np.max(np.abs(cluster_metrics['brier']['parent_stratified_bootstrap95'] - np.asarray(calibration_saved['Brier_cluster_bootstrap_95']))) < 1e-12
assert np.max(np.abs(cluster_metrics['log_score']['parent_stratified_bootstrap95'] - np.asarray(calibration_saved['log_score_cluster_bootstrap_95']))) < 1e-12
all_probabilities = np.asarray(all_probabilities)
all_truth = np.asarray(all_truth)
bin_intervals = []
for bin_record in calibration_saved['calibration_bins']:
    h, low, high = bin_record['H'], bin_record['low'], bin_record['high']
    selected = (all_probabilities[:, h] >= low) & (all_probabilities[:, h] < (high if high < 1 else 1 + 1e-14))
    selected = selected.reshape(24, 30)
    parent_count = selected.sum(axis=1).reshape(3, 8)
    parent_success = (selected * (all_truth.reshape(24, 30) == h)).sum(axis=1).reshape(3, 8)
    parent_p = (selected * all_probabilities[:, h].reshape(24, 30)).sum(axis=1).reshape(3, 8)
    gaps = []
    frequencies = []
    for indices in resamples:
        count = sum(parent_count[topology, indices[topology]].sum() for topology in range(3))
        if count:
            frequency = sum(parent_success[topology, indices[topology]].sum() for topology in range(3)) / count
            probability = sum(parent_p[topology, indices[topology]].sum() for topology in range(3)) / count
            frequencies.append(frequency)
            gaps.append(frequency - probability)
    bin_intervals.append(dict(H=h, low=low, high=high, replicates=int(selected.sum()),
                              material_parents_with_selected_replicate=int(np.sum(selected.any(axis=1))),
                              observed_frequency_cluster95=np.quantile(frequencies, [.025, .975]).tolist(),
                              observed_minus_reported_cluster95=np.quantile(gaps, [.025, .975]).tolist(),
                              effective_nonempty_bootstraps=len(frequencies)))
calibration = dict(reconstruction='exact deterministic seed replay of previously recorded simulation; saved paid parent arrays only',
                   manifest_theta_max_error=float(np.max(np.abs(reconstructed_theta - manifest_theta))),
                   saved_probability_max_error=max(probability_identity_errors), saved_score_max_error=max(saved_metric_errors),
                   independent_material_parents=24, stratification='8 per true topology', noise_replicates=720,
                   noise_real_sd=float(calibration_sigma), bootstrap_resamples=2000, bootstrap_seed=config['seed'] + 401,
                   interval_method='same approximate prior-mass-bin CDF used by original code; no continuous CDF error certificate',
                   metrics=cluster_metrics, calibration_bins_with_cluster_uncertainty=bin_intervals,
                   initial_field_label='conditional_material_90_coverage', correct_initial_semantics='marginal across H',
                   true_H_conditional_intervals_scope='oracle conditioning on known generative H; descriptive model check only',
                   excludes_solver_representation_bias=True, public_real_world_calibration=False)
dump(OUT / 'CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json', dict(calibration=calibration, reconstructed_interval_records=interval_records))

assert ledger_path.read_text() == ledger_before, 'Shared physics ledger changed during reanalysis'
receipt = dict(executable=sys.executable, numpy_version=np.__version__, scipy_version=scipy.__version__,
               imported_B1_modules=False, new_physics_evaluations=0, physics_before=384, physics_after=384,
               new_models_or_installs=False, neural_training=False, historical_or_locked_arrays_read=False,
               shared_source_config_results_modified=False, wall_seconds=time.perf_counter() - wall_start,
               process_cpu_seconds=time.process_time() - cpu_start,
               ru_user_cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime,
               ru_system_cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_stime,
               read_paths=sorted(set(READS)), output_scope=str(OUT))
run_receipts_path = OUT / 'RECOMPUTATION_RUNS_NATIVE.json'
if run_receipts_path.exists():
    run_receipts = json.loads(run_receipts_path.read_text())
else:
    initial_path = OUT / 'RECOMPUTATION_NATIVE.json'
    run_receipts = []
    if initial_path.exists():
        initial_receipt = json.loads(initial_path.read_text())['receipt']
        run_receipts.append(dict(run=1, wall_seconds=initial_receipt['wall_seconds'], process_cpu_seconds=initial_receipt['process_cpu_seconds'],
                                 new_physics_evaluations=0, note='Initial numerical replay; dense covariance slogdet emitted RuntimeWarnings. Stable Cholesky follow-up added finite assertions.'))
run_receipts.append(dict(run=len(run_receipts) + 1, wall_seconds=receipt['wall_seconds'], process_cpu_seconds=receipt['process_cpu_seconds'],
                         new_physics_evaluations=0, note='Stable Cholesky covariance cross-check and saved NPZ gain-default reproducibility check; no physical solver imported.'))
dump(run_receipts_path, run_receipts)
receipt['total_reanalysis_wall_seconds'] = sum(run['wall_seconds'] for run in run_receipts)
receipt['total_reanalysis_process_cpu_seconds'] = sum(run['process_cpu_seconds'] for run in run_receipts)
receipt['reanalysis_runs'] = len(run_receipts)
result = dict(receipt=receipt, accounting=accounting, uniform_prior_normalization=normalizations,
              dense_gain_covariance_max_log_likelihood_error=max(dense_covariance_errors),
              cases=case_results, interpolation=interpolation, importance=importance,
              model_correction=mechanism, calibration=calibration)
dump(OUT / 'RECOMPUTATION_NATIVE.json', result)
print(json.dumps(dict(new_physics=0, wall_seconds=receipt['wall_seconds'], cpu_seconds=receipt['process_cpu_seconds'],
                      calibration_probability_identity_max_error=calibration['saved_probability_max_error'],
                      marginal_real90=cluster_metrics['marginal_real90_coverage'], marginal_imag90=cluster_metrics['marginal_imag90_coverage'],
                      true_H_conditional_real90=cluster_metrics['truth_H_conditional_real90_coverage'],
                      true_H_conditional_imag90=cluster_metrics['truth_H_conditional_imag90_coverage'],
                      top2_gap=cluster_metrics['top2_coverage_minus_mass'])))
