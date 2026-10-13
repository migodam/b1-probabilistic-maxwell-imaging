"""Independent direct image moments from paid posterior quadrature records."""
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
config = json.loads((ROOT / 'configs/protocol.json').read_text())
records = {}
for folder in sorted((ROOT / 'results/cases').iterdir()):
    if not folder.is_dir():
        continue
    with np.load(folder / 'posterior.npz', allow_pickle=False) as file:
        H, theta, W = file['H'], file['theta'], file['posterior_weights']
    with np.load(folder / 'material_fields.npz', allow_pickle=False) as file:
        saved_mean, saved_sd, saved_void = file['posterior_mean'], file['material_complex_sd'], file['void_probability']
        axis = file['axis']
    X, Y = np.meshgrid(axis, axis)
    radial = np.hypot(X, Y)
    inside = radial <= config['outer_radius']
    complex_material = theta[:, 0] + 1j * theta[:, 1]
    mean = np.ones(radial.shape, complex)
    second_moment = np.ones(radial.shape)
    void_probability = np.zeros(radial.shape)
    for h in range(3):
        void = np.zeros(radial.shape, bool)
        for lower, upper in config['topologies'][h]['void_intervals']:
            void |= (radial >= lower) & (radial <= upper) & inside
        material = inside & ~void
        ids = H == h
        mean += material * (W[ids] @ (complex_material[ids] - 1))
        second_moment += material * (W[ids] @ (np.abs(complex_material[ids])**2 - 1))
        void_probability += float(W[ids].sum()) * void
    variance = second_moment - np.abs(mean)**2
    mean_error = float(np.max(np.abs(mean - saved_mean)))
    variance_error = float(np.max(np.abs(variance - saved_sd**2)))
    void_error = float(np.max(np.abs(void_probability - saved_void)))
    assert mean_error < 1e-12 and variance_error < 1e-12 and void_error < 1e-12
    probability_record = json.loads((folder / 'probabilities.json').read_text())
    p = np.asarray(probability_record['result']['topology_probabilities'])
    displayed = [item['posterior_mass'] for item in probability_record['top_k']['items']]
    remaining = probability_record['top_k']['remaining_modeled_mass']
    assert abs(sum(displayed) + remaining - 1) < 1e-12
    assert abs(remaining - p.min()) < 1e-12
    records[folder.name] = dict(posterior_mean_max_error=mean_error, complex_variance_max_error=variance_error,
                                void_probability_max_error=void_error, top2_remaining_mass=remaining,
                                scope='2D central section of fixed 3D concentric geometry; center is not inferred')

assert (ROOT / 'results/PHYSICS_LEDGER.jsonl').read_text() == ledger
result = dict(new_physics_evaluations=0, physics_before=384, physics_after=384, imported_B1_modules=False,
              executable=sys.executable, cases=records, wall_seconds=time.perf_counter() - wall_start,
              process_cpu_seconds=time.process_time() - cpu_start)
(OUT / 'FIELD_AUDIT_NATIVE.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
print(json.dumps(dict(cases=len(records), new_physics=0, wall_seconds=result['wall_seconds'], cpu_seconds=result['process_cpu_seconds'],
                      max_field_mean_error=max(item['posterior_mean_max_error'] for item in records.values()),
                      max_field_variance_error=max(item['complex_variance_max_error'] for item in records.values()))))
