"""Independent NumPy full-VIE check of frozen materials; never reads truth.

This code deliberately does not import any solver, reduced kernel, physics
Model, retained LU, or scoring code. NumPy solve forms a fresh full matrix for
each returned candidate. Operator construction is shared with generation.
"""
import argparse
import datetime
import json
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def frozen_files(manifest_path, disclosed_replay=False):
    manifest = json.loads(manifest_path.read_text())
    for item in manifest['files']:
        path = ROOT / item['snapshot']
        stat = path.stat()
        unchanged = stat.st_size == item['bytes'] if disclosed_replay else (stat.st_size, stat.st_mtime_ns) == (item['bytes'], item['mtime_ns'])
        if not unchanged:
            raise RuntimeError('Frozen snapshot identity changed: ' + str(path))
        yield path, json.loads(path.read_text())


def validate(path, result, split, mode='all'):
    case = result.get('case_id')
    if case is None:
        case = next(x.stem for x in (ROOT/'public_data'/split).glob('*.npz') if x.stem in path.name)
    with np.load(ROOT / 'public_data' / split / (case + '.npz'), allow_pickle=False) as f:
        G, E0, S, y, box = (f[x] for x in ('G', 'E0', 'S', 'y', 'box'))
    n, ns = E0.shape
    indices = range(len(result['candidates'])) if mode == 'all' else sorted(set([result['top1_index']] + [i for i,c in enumerate(result['candidates']) if c['relative_residual'] <= 1e-6]))
    started = time.perf_counter()
    rows, rhs, factors = [], 0, 0
    for index in indices:
        c = result['candidates'][index]
        m = np.asarray(c['m'], dtype=float)
        row = dict(index=index, box_admissible=bool(m.shape == (n,) and np.all(np.isfinite(m)) and np.all(m >= box[0]-1e-8) and np.all(m <= box[1]+1e-8)))
        try:
            A = np.eye(n, dtype=complex) - m[:, None] * G
            b = m[:, None] * E0
            rhs += ns
            factors += 1
            current = np.linalg.solve(A, b)
            predicted = S @ current
            sigma = np.linalg.svd(A, compute_uv=False)
            relative = float(np.linalg.norm(predicted-y) / max(np.linalg.norm(y), 1e-300))
            equation = float(np.linalg.norm(A@current-b) / max(np.linalg.norm(b), 1e-300))
            field = E0 + G@current
            constitutive = float(np.linalg.norm(current-m[:,None]*field) / max(np.linalg.norm(current), 1e-300))
            disagreement = abs(relative-c['relative_residual'])
            row.update(relative_residual=relative, state_equation_relative=equation,
                       constitutive_relative=constitutive, resolvent_sigma_min=float(sigma[-1]),
                       resolvent_relative_gap=float(sigma[-1]/sigma[0]),
                       stored_residual_abs_disagreement=disagreement,
                       data_accepted=bool(relative <= 1e-6),
                       validation_pass=bool(row['box_admissible'] and sigma[-1]/sigma[0] > 1e-12 and equation < 1e-10 and constitutive < 1e-10 and disagreement < 1e-8))
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as error:
            row.update(validation_pass=False, error=type(error).__name__+': '+str(error))
        rows.append(row)
    return dict(file=str(path.relative_to(ROOT)), case_id=case, method=result['method'],
                candidate_count=len(result['candidates']), validated_count=len(rows),
                validation_mode=mode, rows=rows, truth_read=False,
                cost=dict(full_complex_rhs=rhs, full_factorizations=factors, full_dimension=n,
                          independent_numpy_solve_calls=factors, wall_seconds=time.perf_counter()-started),
                all_validated_pass=all(x['validation_pass'] for x in rows))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--split', choices=['dev','holdout','confirm'], required=True)
    p.add_argument('--freeze-manifest', type=Path, required=True)
    p.add_argument('--mode', choices=['all','critical'], default='all')
    p.add_argument('--output', type=Path, required=True)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    rows=[validate(path,result,a.split,a.mode) for path,result in frozen_files(a.freeze_manifest)]
    payload=dict(validated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 freeze_manifest=str(a.freeze_manifest), rows=rows, truth_read=False,
                 full_complex_rhs=sum(r['cost']['full_complex_rhs'] for r in rows),
                 full_factorizations=sum(r['cost']['full_factorizations'] for r in rows),
                 wall_seconds=sum(r['cost']['wall_seconds'] for r in rows))
    a.output.write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:payload[k] for k in ('full_complex_rhs','full_factorizations','wall_seconds')}))


if __name__ == '__main__':
    main()
