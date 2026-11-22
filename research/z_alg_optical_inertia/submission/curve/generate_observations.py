"""Offline-only paired noise construction and true-material ball audit."""

import json
from pathlib import Path
from time import perf_counter, process_time
import sys

import numpy as np
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY))
from core.certificate import ball, frobenius_squared, safe_float

NOISE_PCT = (0.5, 1, 2, 5, 10)
GUARD = arb(1) / 10**6


def generate():
    started, cpu = perf_counter(), process_time()
    ctx.prec = 128
    public = STUDY / 'synthetic/public'
    evaluation = STUDY / 'synthetic/evaluation'
    old = json.loads((STUDY / 'results/theory/fine_feasibility.json').read_text())
    bounds = {r['id']: r for r in old['observations']}
    out = HERE / 'observations'
    out.mkdir(exist_ok=True)
    rows, audit = [], []
    for idx in range(1, 7):
        object_id = f'optfresh_{idx:02d}'
        stem = f'{object_id}__pw360_10x20'
        old_id = stem + '__n02'
        with np.load(public / f'{old_id}.npz', allow_pickle=False) as inp:
            archived = inp['Y']
        with np.load(evaluation / f'{old_id}.npz', allow_pickle=False) as inp:
            clean = inp['Y_fine_clean']
        clean_error = arb(bounds[old_id]['true_forward_to_saved_numerical_clean_upper'])
        direction = archived - clean
        exact_clean, exact_archive = ball(clean), ball(archived)
        for pct in NOISE_PCT:
            y = archived if pct == 2 else clean + (pct / 2.) * direction
            exact_y = ball(y)
            d = frobenius_squared(exact_y - exact_clean).sqrt()
            norm_y = frobenius_squared(exact_y).sqrt()
            sample_radius = (bounds[old_id]['sample_radius_binary64'] if pct == 2
                             else safe_float(d, upper=True))
            radius = (bounds[old_id]['inflated_radius'] if pct == 2
                      else safe_float(arb(sample_radius) + GUARD * norm_y, upper=True))
            true_distance = (d + clean_error).upper()
            construction_error = frobenius_squared(
                exact_y - exact_clean - (exact_archive - exact_clean) * (arb(pct) / 2)
            ).sqrt()
            if not arb(radius) > true_distance:
                raise ArithmeticError(f'Unverified padded ball for {object_id}, {pct}%')
            case_id = f'{stem}__pair{pct:g}pct'
            file = f'{case_id}.npz'
            np.savez_compressed(out / file, Y=y)
            rows.append(dict(id=case_id, file=file, object_id=object_id,
                             geometry='PW360_10x20', noise_pct=pct,
                             sample_radius=sample_radius, padded_radius=radius))
            record = dict(id=case_id, forward_enclosure_source=old_id,
                          exact_archived_2pct_center=bool(np.array_equal(y, archived)),
                          nominal_direction_affine_error_upper=safe_float(construction_error, upper=True),
                          clean_forward_error_upper=float(clean_error),
                          true_distance_upper=safe_float(true_distance, upper=True),
                          admissibility_margin_lower=safe_float(arb(radius) - true_distance, upper=False))
            if pct == 10:
                independent = stem + '__n10'
                with np.load(public / f'{independent}.npz', allow_pickle=False) as inp:
                    record['archived_independent_10pct_id'] = independent
                    record['center_distance_to_independent_10pct_lower'] = safe_float(
                        frobenius_squared(exact_y - ball(inp['Y'])).sqrt(), upper=False)
                if not record['center_distance_to_independent_10pct_lower'] > 0:
                    raise ArithmeticError('Paired and independent 10% centers must differ.')
            audit.append(record)
    (out / 'manifest.json').write_text(json.dumps(dict(observations=rows), indent=2) + '\n')
    (HERE / 'observation_audit.json').write_text(json.dumps(dict(
        scope='Saved binary64 96² operator, B and S; offline reference used solely for paired input construction and forward-ball audit',
        provenance='synthetic/evaluation saved clean field plus rescaled increment from synthetic/public archived full-ring 2%',
        original_forward_certificate='results/theory/fine_feasibility.json',
        common_relative_guard=1e-6, observations=audit,
        wall_seconds=perf_counter()-started, cpu_seconds=process_time()-cpu), indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    generate()
