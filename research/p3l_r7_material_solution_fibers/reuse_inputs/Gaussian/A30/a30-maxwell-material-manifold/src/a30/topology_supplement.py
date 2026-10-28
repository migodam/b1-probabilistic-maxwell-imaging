#!/usr/bin/env python3
"""Offline digital-topology supplement for a COMPLETE saved-pool A29.1 run.
No physical backend, training, candidate generation or locked data access.
Thresholds depend on truth ONLY for offline diagnostics. Counts are raster
quantities, not certificates of continuous object topology.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import label, generate_binary_structure


def topology(mask: np.ndarray) -> tuple[int, int]:
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError('A two-dimensional mask is required')
    _, components = label(mask, generate_binary_structure(2, 1))
    background_labels, total = label(~mask, generate_binary_structure(2, 2))
    border = np.concatenate((background_labels[0], background_labels[-1],
                             background_labels[:, 0], background_labels[:, -1]))
    exterior = set(int(x) for x in border if x != 0)
    holes = int(total) - len(exterior)
    return int(components), holes


def boundary(mask: np.ndarray) -> np.ndarray:
    p = np.pad(mask, 1, constant_values=False)
    interior = p[:-2, 1:-1] & p[2:, 1:-1] & p[1:-1, :-2] & p[1:-1, 2:]
    return mask & ~interior


def simplex_mean(w: np.ndarray, chi: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    if w.shape != (len(chi),) or not np.isfinite(w).all():
        raise ValueError('Saved weights and pool are inconsistent')
    if np.min(w) < -1e-10 or abs(float(w.sum())-1.) > 1e-8:
        raise ValueError('Invalid saved simplex weights; do not silently repair')
    return w @ chi


def run(root: Path, results: Path) -> None:
    target_path = root/'results/a29/OFFLINE_TARGETS.npz'
    weights_path = results/'A29_1_WEIGHTS.json'
    cases_path = results/'A29_1_CASES.csv'
    for p in (target_path, weights_path, cases_path):
        if not p.is_file():
            raise FileNotFoundError('Full A29.1 inputs/results missing: '+str(p))
    with np.load(target_path, allow_pickle=False) as z:
        targets = {k: z[k].copy() for k in z.files}
    ids = [str(x) for x in targets['case_id']]
    parents = np.asarray(targets['parent_id'], dtype=int)
    if len(ids) != 32 or len(set(ids)) != 32 or set(parents) != set(range(250065,250081)):
        raise ValueError('Exactly 32 exposed cases / 16 original parents required')
    with cases_path.open(newline='') as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 32 or {r['case_id'] for r in rows} != set(ids):
        raise ValueError('The preceding A29.1 run is incomplete')
    weight_rows = json.loads(weights_path.read_text())
    if len(weight_rows) != 32 or {r['case_id'] for r in weight_rows} != set(ids):
        raise ValueError('Saved oracle weights are incomplete')
    weights = {r['case_id']: r for r in weight_rows}
    source_path = Path(__file__).with_name('a29_1_convex_diagnostic.py')
    spec = importlib.util.spec_from_file_location('original_a29_1', source_path)
    if spec is None or spec.loader is None:
        raise ImportError('Cannot load the supplied diagnostic script')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    out = []
    for i, case in enumerate(ids):
        with np.load(root/'results/a29/pools'/f'{case}.npz', allow_pickle=False) as z:
            pool = {k: z[k].copy() for k in z.files}
        chi = np.asarray(pool['chi'], dtype=np.complex128)
        small = np.asarray(pool['small_indices'], dtype=int)
        truth = np.asarray(targets['truth'][i], dtype=np.complex128).ravel()
        bg = np.asarray(targets['background'][i], dtype=np.complex128).ravel()
        den = float(np.linalg.norm(truth-bg))
        if den <= 0 or not np.isfinite(den) or truth.size != 576:
            raise ValueError('Expected a finite nonzero 24x24 contrast')
        saved = weights[case]
        if not np.array_equal(small, np.asarray(saved['small_indices'], int)):
            raise ValueError('Small-pool index mismatch')
        hard = int(np.argmin(np.linalg.norm(chi-truth, axis=1)))
        hard_small = int(small[np.argmin(np.linalg.norm(chi[small]-truth, axis=1))])
        lin_full = original.actual_linear_decision(pool, np.arange(len(chi)), truth, den)
        lin_small = original.actual_linear_decision(pool, small, truth, den)
        methods = {
            'N3': chi[0], 'full_hard_oracle': chi[hard], 'small_hard_oracle': chi[hard_small],
            'full_convex_oracle': simplex_mean(saved['full_oracle'], chi),
            'small_convex_oracle': simplex_mean(saved['small_oracle'], chi[small]),
            'full_L_MAP': chi[lin_full['hard_index']], 'small_L_MAP': chi[lin_small['hard_index']],
            'full_L_mean': lin_full['mean'], 'small_L_mean': lin_small['mean']}
        expected = np.abs(truth-bg).reshape(24,24)
        for fraction in (.15, .25, .35):
            threshold = fraction*float(expected.max())
            tm = expected > threshold
            tc, th = topology(tm)
            te = boundary(tm)
            for name, image in methods.items():
                image = np.asarray(image).ravel()
                pm = np.abs(image-bg).reshape(24,24) > threshold
                pc, ph = topology(pm)
                pe = boundary(pm)
                union = int(np.count_nonzero(tm | pm))
                ec = int(te.sum()+pe.sum())
                out.append({'case_id':case, 'parent_id':int(parents[i]),
                    'background_id':str(targets['background_id'][i]),
                    'method':name, 'threshold_fraction':fraction,
                    'delta_nrmse':float(np.linalg.norm(image-truth)/den),
                    'support_iou':float(np.count_nonzero(tm & pm)/union) if union else 1.,
                    'edge_f1':float(2*np.count_nonzero(te & pe)/ec) if ec else 1.,
                    'pred_components':pc, 'pred_holes':ph, 'pred_euler':pc-ph,
                    'truth_components':tc, 'truth_holes':th, 'truth_euler':tc-th,
                    'count_match_only':pc==tc and ph==th,
                    'passive':bool(np.all(image.imag >= -1e-10))})
    path = results/'A29_1_TOPOLOGY_SUPPLEMENT.csv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0])); writer.writeheader(); writer.writerows(out)
    metadata = {'conditions':32, 'parents':16, 'new_Maxwell_states':0,
        'foreground_connectivity':4, 'background_connectivity':8,
        'thresholds_truth_dependent_offline_only':[.15,.25,.35],
        'primary_continuity_threshold':.25,
        'interpretation':'Raster counts only. Matching counts does not imply matching topology location or continuous geometry.',
        'L_estimator':'Data-dependent bank weighted estimates, not calibrated posterior'}
    (results/'A29_1_TOPOLOGY_METADATA.json').write_text(json.dumps(metadata,indent=2))
    print(path)


def self_test() -> None:
    x = np.zeros((5,5), bool); x[1:4,1:4] = True; x[2,2] = False
    assert topology(x) == (1,1)
    assert topology(np.zeros((5,5),bool)) == (0,0)
    assert topology(np.ones((5,5),bool)) == (1,0)
    x = np.eye(2,dtype=bool)
    assert topology(x) == (2,0)
    print('PASS: four synthetic digital topology cases; no physical solver')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',type=Path)
    parser.add_argument('--results',type=Path)
    parser.add_argument('--self-test',action='store_true')
    a=parser.parse_args()
    if a.self_test:
        self_test()
    else:
        if a.root is None or a.results is None:
            parser.error('--root and --results are required')
        run(a.root.resolve(),a.results.resolve())
