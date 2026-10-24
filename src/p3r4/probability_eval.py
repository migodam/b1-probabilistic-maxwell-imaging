"""Pure full66 generator diagnostics. No Maxwell backend, solver, or fitting.

Latent/material marginal rank diagnostics are distinct from the paid,
data-dependent residual SBC and Gamma256 tests, which are NOT evaluated here.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri

from .evaluation import checkpoint_identity, write_new_json

SAMPLE_COUNT = 128
FIXED_SAMPLE_INDEX = 0
BOOTSTRAP_SEED = 20261025099
CI_LEVEL = 0.9833333333333333
COVERAGE_LEVELS = (.5, .8, .9, .95)
PHASE_THRESHOLDS = (.25, .5, .75)


def safe_id(value):
    return isinstance(value, str) and bool(value) and all(c.isalnum() or c in '_-' for c in value)


def source_major_pack(observation):
    a = np.asarray(observation)
    if a.shape != (32, 8) or not np.iscomplexobj(a) or not np.isfinite(a).all():
        raise ValueError('one finite complex observation[32,8] required')
    return np.concatenate((a.real.T, a.imag.T), axis=1).reshape(1, 512)


def read_online_observation(path, parent, *, need_base=False, target_id=None):
    """Materialize only allowlisted observation/identity/y-derived base arrays.

    Never enumerate or access q, x, chi, truth, noise, prediction, or labels.
    Supports a single parent file or a packed NPZ with parent metadata.
    """
    path = Path(path)
    if 'offline_only' in path.resolve().parts:
        raise ValueError('offline_only is prohibited during online inference')
    with np.load(path, allow_pickle=False) as z:
        observation = z['observation']
        if observation.shape == (32, 8):
            index = None
        else:
            ids = z['parent'].astype(str)
            rows = np.flatnonzero(ids == parent)
            if len(rows) != 1 or observation.ndim != 3:
                raise ValueError('packed observation must have one matching parent')
            index = int(rows[0])
            observation = observation[index]
        # Only scalar acquisition metadata and the declared base are allowed.
        if 'target_id' in z:
            tid = np.asarray(z['target_id'])
            tid = str(tid.item()) if tid.ndim == 0 else str(tid[index])
            if target_id is not None and tid != target_id:
                raise ValueError('online target_id mismatch')
        base = {}
        if need_base:
            for key, shape in (('base_mean', (66,)), ('base_chol', (66, 66))):
                value = z[key]
                value = value if index is None else value[index]
                if value.shape != shape or np.iscomplexobj(value) or not np.isfinite(value).all():
                    raise ValueError('invalid y-only conditional Gaussian base')
                base[key] = value[None]
    return source_major_pack(observation), base


def decoder_from_geometry(config, points):
    # Import the existing *pure* original prior; never import p3r4.core.Target.
    try:
        from p1_material.prior import MaterialPrior
    except ModuleNotFoundError:
        from src.p1_material.prior import MaterialPrior
    prior = MaterialPrior(config, points)
    return dict(config=config, points=np.asarray(points).tolist(), A=prior.A.tolist(),
                roi=np.flatnonzero(prior.mask).tolist())


def read_online_anchor_context(path, parent, *, j_summary='none', coordinate_preconditioner='none'):
    """Explicit context_x key only; ambiguous x may be a supervised GT label."""
    path = Path(path)
    if 'offline_only' in path.resolve().parts:
        raise ValueError('offline_only context forbidden online')
    with np.load(path, allow_pickle=False) as z:
        value = z['context_x']
        index = None
        if value.shape != (66,):
            rows = np.flatnonzero(z['parent'].astype(str) == parent)
            if len(rows) != 1:
                raise ValueError('one matching packed anchor parent required')
            index = int(rows[0])
        def read(key):
            a = z[key]
            return (a if index is None else a[index])[None]
        result = dict(context_x=read('context_x'))
        for key in ('base_mean', 'base_chol'):
            if key in z:
                result[key] = read(key)
        if j_summary != 'none' or coordinate_preconditioner != 'none':
            result['context_J'] = read('context_J' if 'context_J' in z else 'J')
    # The predictor checks all dimensions, finiteness and Gaussian base rank.
    return result


def decode_gaussianized(x, decoder):
    """Vectorized original C2 decoder, with no prior truncation/clipping."""
    try:
        from p1_material.prior import saturation
    except ModuleNotFoundError:
        from src.p1_material.prior import saturation
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 66 or not np.isfinite(x).all():
        raise ValueError('finite Gaussianized x[S,66] required')
    cfg = decoder['config']
    field = cfg['field_mean']+x[:, :64]@np.asarray(decoder['A']).T
    fraction = saturation(field/cfg['phase_transition_halfwidth'])
    mask = np.zeros(576, bool)
    mask[np.asarray(decoder['roi'])] = True
    fraction *= mask[None]
    u = ndtr(x[:, 64:])
    re, im = cfg['inclusion_increment_real_interval'], cfg['inclusion_increment_imag_interval']
    delta = re[0]+(re[1]-re[0])*u[:, 0]+1j*(im[0]+(im[1]-im[0])*u[:, 1])
    chi = complex(*cfg['host_chi'])+delta[:, None]*fraction
    return chi, fraction, int(np.sum((u == 0) | (u == 1)))


def fixed_projections(decoder, seed=20261025098):
    """16 real material random projections +4 low-frequency spatial ones."""
    roi = np.asarray(decoder['roi'])
    rng = np.random.default_rng(seed)
    random = rng.normal(size=(16, 1152))
    random /= np.linalg.norm(random, axis=1, keepdims=True)
    points = np.asarray(decoder['points'])
    spatial = np.zeros((4, 1152))
    spatial[0, roi] = 1
    spatial[1, 576+roi] = 1
    spatial[2, roi] = points[roi, 0]-points[roi, 0].mean()
    spatial[3, 576+roi] = points[roi, 1]-points[roi, 1].mean()
    spatial /= np.linalg.norm(spatial, axis=1, keepdims=True)
    return dict(seed=seed, matrix=np.concatenate((random, spatial)).tolist(),
                names=[f'random_real_material_{i:02d}' for i in range(16)]+
                      ['roi_real_constant', 'roi_imag_constant', 'roi_real_x', 'roi_imag_y'],
                realification='concat(Re chi, Im chi)/24; unit Euclidean projection vectors',
                spatial_count=4, random_count=16)


def real_material(chi):
    chi = np.asarray(chi)
    return np.concatenate((chi.real, chi.imag), axis=-1)/24


def energy_score(samples, truth):
    """Unbiased distinct-pair empirical energy score under full-grid RMS."""
    samples, truth = np.asarray(samples), np.asarray(truth)
    first = np.linalg.norm(samples-truth, axis=1).mean()
    pair = 0.
    for i in range(1, len(samples)):
        pair += np.linalg.norm(samples[:i]-samples[i], axis=1).sum()
    return float(first-pair/(len(samples)*(len(samples)-1)))


def scalar_crps(samples, truth):
    """Coordinate-wise energy/CRPS using the same distinct-pair convention."""
    samples, truth = np.asarray(samples), np.asarray(truth)
    first = np.abs(samples-truth).mean(0)
    pair = np.zeros_like(first)
    for i in range(1, len(samples)):
        pair += np.abs(samples[:i]-samples[i]).sum(0)
    return first-pair/(len(samples)*(len(samples)-1))


def rank_diagnostics(samples, truth, seed):
    """Randomized inclusive ranks handle exact plateau atoms and discrete ties.

    Central rank coverage is distinct from ordinary empirical interval coverage.
    These are marginal diagnostics, never simultaneous whole-image coverage.
    """
    samples, truth = np.asarray(samples), np.asarray(truth)
    lower = (samples < truth).sum(0)
    ties = (samples == truth).sum(0)
    u = np.random.default_rng(seed).uniform(size=truth.shape)
    rank = (lower+u*(ties+1))/(len(samples)+1)
    result = dict(rank=rank, tie_count=ties)
    for level in COVERAGE_LEVELS:
        alpha = (1-level)/2
        lo, hi = np.quantile(samples, [alpha, 1-alpha], axis=0)
        key = str(int(round(level*100)))
        result['randomized_rank_coverage_'+key] = ((rank >= alpha) & (rank <= 1-alpha)).astype(float)
        result['empirical_interval_coverage_'+key] = ((truth >= lo) & (truth <= hi)).astype(float)
        result['empirical_interval_width_'+key] = hi-lo
    return result


def parent_scores(x_samples, chi_samples, fraction_samples, x_truth, chi_truth,
                  fraction_truth, projections, roi, seed):
    x_samples, chi_samples = np.asarray(x_samples), np.asarray(chi_samples)
    k = len(x_samples)
    if x_samples.shape != (128, 66) or chi_samples.shape != (128, 576):
        raise ValueError('frozen128 full66 draws and full576 decoded worlds required')
    if not np.isfinite(x_samples).all() or not np.isfinite(chi_samples).all():
        raise ValueError('nonfinite online output')
    vectors = real_material(chi_samples)
    truth_vector = real_material(chi_truth)
    matrix = np.asarray(projections['matrix'])
    projected, projected_truth = vectors@matrix.T, truth_vector@matrix.T
    mean = chi_samples.mean(0)
    result = dict(raw_sample_mean_full_mse=float(np.mean(abs(mean-chi_truth)**2)),
                  raw_sample_mean_full_mae=float(np.mean(abs(mean-chi_truth))),
                  material_energy_score=energy_score(vectors, truth_vector),
                  projection_crps=scalar_crps(projected, projected_truth))
    material_roi = np.r_[np.asarray(roi), 576+np.asarray(roi)]
    for name, sample, truth in (('latent', x_samples, x_truth),
                               ('material_real', vectors, truth_vector),
                               ('material_real_roi', vectors[:, material_roi], truth_vector[material_roi]),
                               ('projection', projected, projected_truth)):
        rank = rank_diagnostics(sample, truth, seed)
        for key, value in rank.items():
            result[name+'__'+key] = value
    threshold = np.asarray(PHASE_THRESHOLDS)
    events = fraction_samples[:, roi, None] >= threshold
    observed = fraction_truth[roi, None] >= threshold
    count = events.sum(0)
    probability = count/k
    # This finite predictive action is expressly separate from raw sample law.
    jeffreys = (count+.5)/(k+1)
    result['phase_event_probability'] = probability
    result['phase_truth_event'] = observed.astype(np.uint8)
    result['phase_brier_per_threshold'] = ((probability-observed)**2).mean(0)
    result['phase_log_score_jeffreys_per_threshold'] = -(observed*np.log(jeffreys)+(1-observed)*np.log1p(-jeffreys)).mean(0)
    result['phase_raw_logscore_zero_support_count'] = (((probability == 0) & observed)|((probability == 1) & ~observed)).sum(0)
    result['phase_raw_logscore_per_threshold'] = np.asarray([
        np.inf if result['phase_raw_logscore_zero_support_count'][j] else
        -np.log(np.where(observed[:, j], probability[:, j], 1-probability[:, j])).mean()
        for j in range(3)])
    # Fixed visible phase probability bins; counts are descriptive pixel counts,
    # not independent sample sizes (parent is the inference unit).
    bins = np.minimum((probability*10).astype(int), 9)
    result['phase_reliability_count'] = np.stack([(bins == b).sum(0) for b in range(10)])
    result['phase_reliability_foreground'] = np.stack([((bins == b)*observed).sum(0) for b in range(10)])
    result['phase_reliability_probability_sum'] = np.stack([((bins == b)*probability).sum(0) for b in range(10)])
    return result


def read_truth_after_gate(path):
    """Call only after the CLI has checked freeze + separate-PID online gate."""
    with np.load(path, allow_pickle=False) as z:
        q, chi = z['q'], z['chi']
    if q.shape != (66,) or chi.shape != (576,) or not np.iscomplexobj(chi) or not np.isfinite(q).all() or not np.isfinite(chi).all():
        raise ValueError('invalid original target truth')
    if np.any(q[64:] <= 0) or np.any(q[64:] >= 1):
        raise ValueError('GT uniform property reached probit boundary; no clipping allowed')
    return np.r_[q[:64], ndtri(q[64:])], chi


def bootstrap_mean(values, draws=4000, seed=BOOTSTRAP_SEED, level=CI_LEVEL):
    values = np.asarray(values, float)
    if values.ndim != 1 or not np.isfinite(values).all() or len(values) < 2:
        return dict(status='NOT_EVALUATED_INCOMPLETE_OR_TOO_FEW_PARENTS', planned_parents=len(values),
                    finite_parents=int(np.isfinite(values).sum()), mean_all_parents=None)
    rng = np.random.default_rng(seed)
    means = []
    for start in range(0, draws, 100):
        idx = rng.integers(0, len(values), (min(100, draws-start), len(values)))
        means.extend(values[idx].mean(1))
    alpha = (1-level)/2
    return dict(status='EVALUATED', n_parent=len(values), mean_all_parents=float(values.mean()),
                CI=np.quantile(means, [alpha, 1-alpha]).tolist(), CI_level=level,
                bootstrap_draws=draws, bootstrap_seed=seed, unit='independent parent')


def check_online_gate(frozen, gate, online_directory):
    if gate.get('status') != 'ONLINE_COMPLETE' or gate.get('truth_read') is not False or gate.get('freeze') != frozen:
        raise ValueError('missing/mismatched observation-only online gate')
    if gate['pid'] == os.getpid():
        raise ValueError('GT scoring requires a separate PID')
    if gate['completed_time_ns'] <= frozen['freeze_time_ns']:
        raise ValueError('online completion must follow method freeze')
    # Check every sealed prediction before any truth is opened. No checksums.
    for record in gate['artifacts']:
        if checkpoint_identity(Path(online_directory)/record['relative_path']) != record['identity']:
            raise ValueError('sealed online artifact changed')
    for method in frozen['methods']:
        if checkpoint_identity(method['checkpoint']) != method['checkpoint_identity']:
            raise ValueError('checkpoint changed since method freeze')
