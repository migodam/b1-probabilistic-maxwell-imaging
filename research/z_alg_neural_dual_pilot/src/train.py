"""Train separate observation-aware and featureless dual multiplier proposals."""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter
import json
import os
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

sys.dont_write_bytecode = True
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
for env in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[env] = '1'

import numpy as np
import torch
from threadpoolctl import threadpool_limits

from chart import MU_SCALE, pullback
from coarse import Coarse
from network import DualProposal, observation_features

STUDY = Path(__file__).resolve().parents[1]
TRAIN_DATA = STUDY / 'data/synthetic_observations.npz'
SIDES = (4, 8)
CAPS = (0.50, 0.75, 1.00)
VALIDATION_NOISE = (0.01, 0.02)
SEED = 202611231991
EPOCHS = 4
BATCH = 12
LEARNING_RATE = 0.004
ANCHOR_PENALTY = 0.001
MAX_GRADIENT_COORD_NORM = 0.35
CUTOFF_GUARD_RELATIVE = 1e-8
EXPECTED_KEYS = {'Y', 'radius', 'split', 'object_id', 'noise_fraction'}

OUTPUTS = {'neural': 'neural_dual', 'no_feature': 'no_feature'}


def write_json(path, data):
    temp = path.with_suffix('.json.tmp')
    temp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def training_inputs():
    """Accept only the offline generator's label-free train/validation schema."""
    with np.load(TRAIN_DATA, allow_pickle=False) as data:
        if set(data.files) != EXPECTED_KEYS:
            raise ValueError('Unexpected fields in the training package.')
        y, radius = data['Y'], data['radius']
        split, object_ids = data['split'], data['object_id']
        levels = data['noise_fraction']
    if set(split) != {'train', 'validation'} or len(set(object_ids)) != len(split):
        raise ValueError('Expected separate and independent train/validation objects.')
    if any(not name.startswith(f'dualpilot_{kind}_') for name, kind in zip(object_ids, split, strict=True)):
        raise ValueError('Only material drawn by the new independent generator is accepted.')
    if y.shape != (len(split), len(levels), 20, 10) or radius.shape != y.shape[:2]:
        raise ValueError('Training inputs have incompatible geometry or condition axes.')
    return y, radius, split, object_ids, levels


class ScoreLedger:
    """Charge each numerical factorization and its included operations."""
    def __init__(self):
        self.operations = Counter()
        self.status = Counter()
        self.inertia = Counter()
        self.wall_s = 0.
        self.cpu_s = 0.
        self.calls = 0
        self.gradient_skips = Counter()
        self.clipped = 0

    def add(self, row):
        self.calls += 1
        self.operations.update(row['operations'])
        self.status[str(row['status'])] += 1
        self.inertia[str(row['nu'])] += 1
        self.wall_s += row['wall_s']
        self.cpu_s += row['cpu_s']

    def summary(self):
        return dict(calls=self.calls, status=dict(self.status), inertia_counts=dict(self.inertia),
                    operations=dict(self.operations), score_wall_s=self.wall_s, score_cpu_s=self.cpu_s,
                    gradient_skips=dict(self.gradient_skips), gradient_clips=self.clipped)


def minibatches(rows, rng):
    shuffled = rng.permutation(rows)
    return (shuffled[start:start+BATCH] for start in range(0, len(shuffled), BATCH))


def gradient(y, raw_mu, row, ledger):
    norm2 = float(np.vdot(y, y).real)
    if row['nu'] is None:
        ledger.gradient_skips['unavailable_stationary_branch'] += 1
        return None
    if row['nu'] >= y.shape[1]:
        if row.get('rescue_K_eigenvalue') is None:
            ledger.gradient_skips['unavailable_exhaustion_guide'] += 1
            return None
        if row['rescue_K_gap'] <= CUTOFF_GUARD_RELATIVE * MU_SCALE:
            ledger.gradient_skips['near_state_rescue_tie'] += 1
            return None
    elif row['eigenvalues'] is None:
        ledger.gradient_skips['unavailable_source_spectrum'] += 1
        return None
    gap = row['source_cutoff_gap']
    if gap is not None and gap <= CUTOFF_GUARD_RELATIVE * norm2:
        ledger.gradient_skips['near_source_cutoff_tie'] += 1
        return None
    a, b = row['grad_mu']/norm2, row['grad_lam']/norm2
    scaled_norm = np.linalg.norm(np.concatenate(pullback(raw_mu, a, b)))
    if not np.isfinite(scaled_norm):
        ledger.gradient_skips['nonfinite_gradient'] += 1
        return None
    scale = min(1., MAX_GRADIENT_COORD_NORM / max(float(scaled_norm), 1e-30))
    ledger.clipped += int(scale < 1)
    return a * scale, b * scale


def run_train_epoch(model, optimizer, co, ys, features, train, epoch, n_noise):
    score_wall, score_cpu = perf_counter(), process_time()
    model.train()
    ledger = ScoreLedger()
    losses = []
    smooth_losses = []
    rng = np.random.default_rng(SEED + epoch)
    noise_idx = rng.integers(0, n_noise, size=len(ys))
    sample = np.array([(i, noise_idx[i], cap) for i in train for cap in range(len(CAPS))])
    for batch in minibatches(sample, rng):
        optimizer.zero_grad(set_to_none=True)
        encoded = model.encoder(torch.from_numpy(features[tuple(batch.T)]))
        loss = 0.
        for side in SIDES:
            mu, lam, raw = model.from_encoded(encoded, side)
            mu_values, lam_values = mu.detach().numpy(), lam.detach().numpy()
            raw_mu = raw[:, 0].detach().numpy()
            a, b = np.zeros_like(mu_values), np.zeros_like(lam_values)
            for i, (obj, noise, cap) in enumerate(batch):
                y = ys[obj, noise]
                row = co.eval(y, CAPS[cap], mu_values[i], lam_values[i], smooth=True)
                ledger.add(row)
                g = gradient(y, raw_mu[i], row, ledger)
                if g is not None:
                    a[i], b[i] = g
                norm2 = float(np.vdot(y, y).real)
                losses.append(float(row['d2_exact']) / norm2)
                smooth_losses.append(float(row['d2']) / norm2)
            # Values are computed by the numerical solver. Use its local derivative
            # as a detached cotangent for the differentiable multiplier chart.
            loss -= ((mu * torch.from_numpy(a)).sum() + (lam * torch.from_numpy(b)).sum()) / len(batch)
            loss += ANCHOR_PENALTY * raw.square().mean()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
    return dict(**ledger.summary(), mean_hard_d2_over_observation_energy=float(np.mean(losses)),
                mean_smooth_guide_over_observation_energy=float(np.mean(smooth_losses)),
                wall_s=perf_counter()-score_wall, cpu_s=process_time()-score_cpu)


@torch.inference_mode()
def validate(model, co, ys, radii, features, validation, noise_indexes, object_ids):
    started, cpu = perf_counter(), process_time()
    model.eval()
    ledger = ScoreLedger()
    rows = []
    for obj in validation:
        values = []
        for noise in noise_indexes:
            y, r = ys[obj, noise], radii[obj, noise]
            y2 = float(np.vdot(y, y).real)
            batch = torch.from_numpy(features[obj, noise])
            encoded = model.encoder(batch)
            for side in SIDES:
                mu, lam, _ = model.from_encoded(encoded, side)
                for cap_ix, cap in enumerate(CAPS):
                    row = co.eval(y, cap, mu[cap_ix].numpy(), lam[cap_ix].numpy(), smooth=False)
                    ledger.add(row)
                    values.append(dict(side=side, cap=cap, noise_index=int(noise),
                                       hard_d2=float(row['d2_exact']),
                                       energy_ratio=float(row['d2_exact'])/y2,
                                       d_over_noise_radius=float(np.sqrt(row['d2_exact'])/r),
                                       nu=row['nu']))
        rows.append(dict(object_id=str(object_ids[obj]), values=values))
    stats = {}
    for side in SIDES:
        stats[str(side)] = {}
        for cap in CAPS:
            group = [v for obj in rows for v in obj['values'] if v['side'] == side and v['cap'] == cap]
            stats[str(side)][str(cap)] = dict(
                mean_hard_d2_over_observation_energy=float(np.mean([v['energy_ratio'] for v in group])),
                mean_hard_d_over_radius=float(np.mean([v['d_over_noise_radius'] for v in group])),
                count_numerical_margin_positive=sum(v['d_over_noise_radius'] > 1. for v in group))
    return dict(**ledger.summary(),
                selector=float(np.mean([np.mean([v['energy_ratio'] for v in obj['values']]) for obj in rows])),
                coarse_rows_by_side_and_cap=stats, objects=rows,
                wall_s=perf_counter()-started, cpu_s=process_time()-cpu)


def cache_inputs(ys, radii, use_observation):
    cached = []
    for group, rr in zip(ys, radii, strict=True):
        features = [observation_features(y, cap, r, use_observation=use_observation)
                    for y, r in zip(group, rr, strict=True) for cap in CAPS]
        cached.append(np.stack(features).reshape(len(rr), len(CAPS), *features[0].shape))
    return np.stack(cached)


def training_recipe(model):
    return dict(data_file=str(TRAIN_DATA.relative_to(STUDY)), caps=CAPS, block_sides=SIDES,
                epochs=EPOCHS, batch_size=BATCH, seed=SEED, optimizer='Adam',
                learning_rate=LEARNING_RATE, anchor_penalty=ANCHOR_PENALTY,
                sample_noise='One input per training object per epoch, uniformly from five paired levels.',
                objective='Mean smoothed retained-source d² divided by observation energy, with chart anchor penalty. '
                          'Source exhaustion uses the shared geometry-only signed state-spectrum guide.',
                gradient_description='Analytic stable-branch gradient. Hard inertia. Source cutoff near-ties skipped, '
                                     'finite scaled per-observation gradient norm clipped.',
                source_cutoff_gap_relative=CUTOFF_GUARD_RELATIVE,
                scaled_per_observation_gradient_clip=MAX_GRADIENT_COORD_NORM,
                validation_selector='Largest object-mean hard d²/||Y||² across validation noise, caps, and '
                                    'both heads. Ties retain the earliest epoch, including initialization.',
                validation_noise_levels=VALIDATION_NOISE, feature_mode='measured' if model.use_observation else 'geometry_and_cap_only',
                chart=model.checkpoint()['chart'])


def train(mode, output):
    started, cpu = perf_counter(), process_time()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.manual_seed(SEED)
    model = DualProposal(use_observation=mode == 'neural')
    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
    ys, rr, split, ids, levels = training_inputs()
    train_ids = np.flatnonzero(split == 'train')
    validation_ids = np.flatnonzero(split == 'validation')
    if len(train_ids) != 60 or len(validation_ids) != 12:
        raise ValueError('The frozen recipe requires 60 training and 12 validation objects.')
    noise = [bisect_left(levels, p) for p in VALIDATION_NOISE]
    if tuple(levels[noise]) != VALIDATION_NOISE:
        raise ValueError('The specified validation levels must exist in the offline data.')
    input_wall, input_cpu = perf_counter(), process_time()
    co = Coarse()
    with threadpool_limits(limits=1, user_api='blas'):
        feats = cache_inputs(ys, rr, model.use_observation)
    precompute = dict(wall_s=perf_counter()-input_wall, cpu_s=process_time()-input_cpu,
                      feature_evaluations=int(np.prod(feats.shape[:3])), geometry_setup=co.setup_cost)
    history = []
    best_value = -1.
    best_trained_value = -1.
    output.mkdir(parents=True, exist_ok=True)
    base = OUTPUTS[mode]
    weights = output / f'{base}_weights.pt'
    pending = output / f'.{base}_selected.pending.pt'
    result_path = output / f'{base}_training.json'
    best_trained = output / f'{base}_best_trained_weights.pt'
    last_trained = output / f'{base}_last_trained_weights.pt'
    recipe = training_recipe(model)
    for epoch in range(EPOCHS + 1):
        t = None if epoch == 0 else run_train_epoch(model, optimizer, co, ys, feats, train_ids, epoch, len(levels))
        val = validate(model, co, ys, rr, feats, validation_ids, noise, ids)
        history.append(dict(epoch=epoch, train=t, validation=val))
        checkpoint = dict(**model.checkpoint(), selected_epoch=epoch, validation_selector=val['selector'])
        if epoch:
            if val['selector'] > best_trained_value:
                best_trained_value = val['selector']
                best_trained_epoch = epoch
                torch.save(checkpoint, best_trained)
            if epoch == EPOCHS:
                torch.save(checkpoint, last_trained)
        if val['selector'] > best_value:
            best_value = val['selector']
            best_epoch = epoch
            torch.save(checkpoint, pending)
        print(f'{mode} epoch={epoch} validation={val["selector"]:.7g} selected={best_epoch} '
              f'wall={perf_counter()-started:.2f}s', flush=True)
    operations = Counter()
    slices = [r[key] for r in history for key in ('train', 'validation') if r[key] is not None]
    for phase in slices:
        operations.update(phase['operations'])
    pending.replace(weights)
    summary = dict(
        freeze='Weights and validation rule finalized before any holdout request; no holdout assets opened.',
        recipe=recipe, selected_epoch=best_epoch, selected_validator=best_value,
        param_count=sum(p.numel() for p in model.parameters()), dataset_train=list(ids[train_ids]),
        dataset_validation=list(ids[validation_ids]), feature_preparation=precompute, history=history,
        total_solve_calls=sum(r['calls'] for r in slices), total_operations=dict(operations),
        host=platform.platform(), torch_version=torch.__version__, numpy_version=np.__version__,
        wall_s=perf_counter()-started, cpu_s=process_time()-cpu, checkpoint=weights.name,
        best_trained_checkpoint=best_trained.name, best_trained_epoch=best_trained_epoch,
        best_trained_validator=best_trained_value, last_trained_checkpoint=last_trained.name)
    write_json(result_path, summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=tuple(OUTPUTS), required=True)
    parser.add_argument('--out', type=Path, default=STUDY / 'results')
    args = parser.parse_args()
    train(args.mode, args.out)
