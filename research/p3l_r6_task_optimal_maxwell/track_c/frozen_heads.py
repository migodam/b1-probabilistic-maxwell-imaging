"""Observation/anchor-only Frozen Head Robustness and current-action callbacks.

This module never reads independent-test truths or historical confirmation
scores. The historical frozen checkpoint loader is imported lazily and in place.
The returned predictions can be scored by a separate offline root evaluator.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import json
import sys
import time
from pathlib import Path
import numpy as np

TRACK = Path(__file__).resolve().parent
REPO = TRACK.parents[2]
CONFIRM_SOURCE = REPO.parent/'b1-p3l-r5-b-confirm/research/p3l_r5_b_independent_confirmation/src'
METHODS = ('physics', 'strongCNN', 'deterministic_matched', 'distilled_w128_n32768',
           'original_K32_frozen')
SEEDS = (2026102601, 2026102602, 2026102603)
ONLINE_KEYS = ('observation', 'residual', 'q0', 'chi0', 'J')


def load_action(method, seed=SEEDS[0]):
    if method not in METHODS or int(seed) not in SEEDS:
        raise ValueError('frozen original method/seed required')
    if str(CONFIRM_SOURCE) not in sys.path:
        sys.path.insert(0, str(CONFIRM_SOURCE))
    import torch
    torch.set_num_threads(1)
    from model_actions import Action
    return Action(method, int(seed), device='cpu')


def packed_measurement(y):
    """Exact Original Target.pack: source-major real32 then imag32 per source."""
    y = np.asarray(y)
    if y.shape != (32, 8) or not np.iscomplexobj(y) or not np.isfinite(y).all():
        raise ValueError('finite complex32x8 observation required')
    return np.concatenate((y.real.T, y.imag.T), axis=1).ravel()


def predict_contexts(contexts, action, parents=None, batch_size=64):
    """No truth/metadata is passed into the original frozen action."""
    if not contexts:
        raise ValueError('nonempty online contexts required')
    if parents is None:
        parents = [str(c['parent']) for c in contexts]
    if len(parents) != len(contexts):
        raise ValueError('parent/context length mismatch')
    online = {key: np.stack([np.asarray(c[key]) for c in contexts]) for key in ONLINE_KEYS}
    y = np.stack([packed_measurement(c['observation']) for c in contexts])
    return action.predict(y, online, parents=list(parents), batch_size=batch_size)


def decision_predictor(action):
    """Callback requires a J already paid by the optimizer; never calls Target."""
    def predict(y, point, parent):
        if point.J is None:
            raise RuntimeError('current action unavailable without cached J; no extra solve')
        context = dict(observation=y, residual=y-point.state.data, q0=point.q,
                       chi0=point.state.chi, J=point.J, parent=np.array(str(parent)))
        return predict_contexts([context], action, [str(parent)], batch_size=1)[0]
    return predict


def evaluate_frozen_heads(context_paths, destination, methods=METHODS, seeds=SEEDS):
    """Save all requested original head predictions and measured inference costs.

    Loading checkpoint, prediction wall/CPU, and physics-zero status are separate
    fields. This is Frozen Head Robustness; a swapped anchor never constitutes a
    matched-head training result.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if (destination/'FROZEN_HEADS.json').exists():
        raise FileExistsError('existing frozen-head evaluation; no silent overwrite')
    contexts = []
    for path in context_paths:
        with np.load(path, allow_pickle=False) as z:
            contexts.append({key: z[key].copy() for key in (*ONLINE_KEYS, 'parent')})
    records = []
    for method in methods:
        for seed in ((SEEDS[0],) if method == 'physics' else seeds):
            start, cpu = time.perf_counter(), time.process_time()
            action = load_action(method, seed)
            load_wall, load_cpu = time.perf_counter()-start, time.process_time()-cpu
            start, cpu = time.perf_counter(), time.process_time()
            prediction = predict_contexts(contexts, action)
            prediction_wall, prediction_cpu = time.perf_counter()-start, time.process_time()-cpu
            path = destination/f'{method}_s{seed}.npz'
            if path.exists():
                raise FileExistsError(path)
            np.savez_compressed(path, prediction=prediction,
                                parent=np.array([str(c['parent']) for c in contexts]))
            records.append(dict(method=method, seed=int(seed), parents=len(contexts),
                                output=str(path), checkpoint_load_wall_seconds=load_wall,
                                checkpoint_load_CPU_seconds=load_cpu,
                                prediction_wall_seconds=prediction_wall,
                                prediction_CPU_seconds=prediction_cpu,
                                new_physics_RHS=0, new_LU=0, GT_read=False,
                                classification='FROZEN_HEAD_ROBUSTNESS'))
            del action
    report = dict(status='COMPLETE', records=records, GT_read=False,
                  original_R5_frozen_checkpoint_loader=str(CONFIRM_SOURCE/'model_actions.py'),
                  matched_lowcost_heads_status='NOT_RUN_REQUIRES_LOW_ANCHOR_TRAINING_CONTEXTS_AND_ACCOUNTED_TRAINING_JOB')
    (destination/'FROZEN_HEADS.json').write_text(json.dumps(report, indent=2)+'\n')
    return report
