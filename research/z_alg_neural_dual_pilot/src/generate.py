"""Create label-free synthetic training observations on Paper 1's 96² operator."""
from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from scipy.fft import fft2, ifft2
from threadpoolctl import threadpool_limits

STUDY = Path(__file__).resolve().parents[1]
RESEARCH = STUDY.parent
GEOMETRY = RESEARCH / 'z_alg_optical_inertia/synthetic/public/fine_geometry.npz'
sys.path.insert(0, str(RESEARCH / 'srap_cie_array_study_20261121'))
from legacy_la import fine_solve

FINE_SIDE = 96
TRAIN_COUNT = 60
VALIDATION_COUNT = 12
SHAPE_SEEDS = {'train': 202611231371, 'validation': 202611231372}
NOISE_SEEDS = {'train': 202611231393, 'validation': 202611231394}
NOISE_LEVELS = np.array([0.005, 0.01, 0.02, 0.05, 0.10])
GENERATING_CAP = 4.0
SOLVE_RTOL = 1e-10

_GEOMETRY = None


def atomic_json(path, obj):
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def open_geometry():
    global _GEOMETRY
    with np.load(GEOMETRY, allow_pickle=False) as g:
        spectrum, s, b = g['spectrum'], g['PW360_10x20__S'], g['PW360_10x20__B']
    diagonal = ifft2(spectrum)[0, 0]
    _GEOMETRY = spectrum, s, b, diagonal


def random_material(seed):
    """Continuous random shapes, sampled at 96² cell centers in a 2 m square."""
    rng = np.random.default_rng(seed)
    coords = 2 * (np.arange(FINE_SIDE) + 0.5) / FINE_SIDE - 1
    x, y = np.meshgrid(coords, coords, indexing='ij')
    x, y = np.stack((x, y)) - rng.uniform(-0.27, 0.27, size=(2, 1, 1))
    theta = rng.uniform(-np.pi, np.pi)
    c, s = np.cos(theta), np.sin(theta)
    x, y = c*x + s*y, -s*x + c*y
    x /= rng.uniform(0.37, 0.68)
    y /= rng.uniform(0.31, 0.67)
    family = int(rng.integers(4))
    radial = np.hypot(x, y)
    phi = np.arctan2(y, x)
    heights = rng.uniform(0.80, 1.0)
    if family == 0:
        boundary = 0.84 + rng.uniform(0.06, 0.16)*np.sin(int(rng.integers(3, 8))*phi + theta)
        boundary += rng.uniform(0.03, 0.16)*np.cos(2*phi - theta)
        fraction = (radial < boundary) * (heights + (1-heights)*np.maximum(0, 1-radial))
    elif family == 1:
        ring = radial < 0.91 + 0.10*np.cos(3*phi + theta)
        hx, hy = rng.uniform(-0.22, 0.22, 2)
        hole = ((x-hx)/rng.uniform(0.25, 0.54))**2 + ((y-hy)/rng.uniform(0.22, 0.49))**2 < 1
        fraction = (ring & ~hole) * (heights + (1-heights)*np.cos(y)**2)
    elif family == 2:
        fraction = np.zeros(x.shape)
        for _ in range(int(rng.integers(2, 6))):
            mx, my = rng.uniform(-0.61, 0.61, 2)
            lobe = ((x-mx)/rng.uniform(0.26, 0.66))**2 + ((y-my)/rng.uniform(0.23, 0.65))**2 < 1
            fraction = np.maximum(fraction, rng.uniform(0.56, 1.0)*lobe)
    else:
        p = rng.uniform(2.6, 5.0)
        body = np.abs(x)**p + np.abs(y)**p < 0.92
        z0, z1 = rng.uniform(-0.44, 0.44, 2)
        slit = (np.abs(y-z0-0.22*x) < rng.uniform(0.06, 0.19)) & (x < z1)
        slit |= (np.abs(x-z1+0.16*y) < rng.uniform(0.075, 0.21)) & (y > z0)
        fraction = (body & ~slit)*(heights + (1-heights)*np.cos(1.8*x + 1.1*y)**2)
    peak = rng.uniform(3.0, GENERATING_CAP)
    return (peak * (fraction / fraction.max())).ravel(), family


def bound_norm(values):
    norm = np.linalg.norm(values)
    return np.nextafter(norm * (1 + 8 * values.size * np.finfo(float).eps), np.inf)


def generate_object(split, index):
    s0, c0 = perf_counter(), process_time()
    spectrum, sensor, incident, diagonal = _GEOMETRY
    green_calls = 0

    def apply(currents):
        nonlocal green_calls
        green_calls += currents.shape[1]
        grid = currents.reshape(FINE_SIDE, FINE_SIDE, -1)
        ft = fft2(grid, s=spectrum.shape, axes=(0, 1))
        projected = ifft2(ft*spectrum[:, :, None], axes=(0, 1))[:FINE_SIDE, :FINE_SIDE]
        return projected.reshape(currents.shape)

    seed = np.random.SeedSequence([SHAPE_SEEDS[split], index])
    contrast, family = random_material(seed)
    with threadpool_limits(limits=1):
        clean, stats = fine_solve(contrast, incident, sensor, apply, diagonal, rtol=SOLVE_RTOL)
    rng = np.random.default_rng(np.random.SeedSequence([NOISE_SEEDS[split], index]))
    direction = rng.standard_normal(clean.shape) + 1j*rng.standard_normal(clean.shape)
    direction *= np.linalg.norm(clean) / np.linalg.norm(direction)
    observations = clean[None, :, :] + NOISE_LEVELS[:, None, None] * direction[None, :, :]
    radii = np.array([bound_norm(y - clean) + 1e-6*bound_norm(y) for y in observations])
    object_id = f'dualpilot_{split}_{index:03d}'
    receipt = dict(object_id=object_id, split=split, shape_seed=[SHAPE_SEEDS[split], index],
                   noise_seed=[NOISE_SEEDS[split], index], generator_family=family, **stats,
                   fft_green_rhs=green_calls, nominal_noise_levels=NOISE_LEVELS.tolist(),
                   wall_s=perf_counter()-s0, cpu_s=process_time()-c0)
    return observations, radii, receipt


def generate(out, workers):
    t0, cpu0 = perf_counter(), process_time()
    tasks = [('train', idx) for idx in range(TRAIN_COUNT)]
    tasks.extend(('validation', idx) for idx in range(VALIDATION_COUNT))
    out.mkdir(parents=True, exist_ok=True)
    records = {}
    with ProcessPoolExecutor(max_workers=workers, initializer=open_geometry) as executor:
        futures = {executor.submit(generate_object, *task): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            try:
                records[task] = future.result()
            except Exception as error:
                atomic_json(out / 'generation_failure.json',
                            dict(failed_task=task, error=str(error),
                                 completed=[value[2] for value in records.values()],
                                 wall_s=perf_counter()-t0))
                raise
            n = len(records)
            if n % 6 == 0:
                print(f'Generated {n}/{len(tasks)} native objects in {perf_counter()-t0:.1f}s', flush=True)
    ys, radii, receipts = zip(*(records[task] for task in tasks), strict=True)
    write_start = perf_counter()
    temp = out / 'synthetic_observations.tmp.npz'
    np.savez_compressed(temp, Y=np.stack(ys), radius=np.stack(radii), noise_fraction=NOISE_LEVELS,
                        split=np.array([r['split'] for r in receipts]),
                        object_id=np.array([r['object_id'] for r in receipts]))
    temp.replace(out / 'synthetic_observations.npz')
    receipt = dict(
        contract='2D lossless scalar-TM fixed Green; sparse full ring PW360_10x20, 400 MHz',
        exact_geometry_file=str(GEOMETRY.relative_to(RESEARCH.parent)),
        historical_forward_entrypoint='research/srap_cie_array_study_20261121/legacy_la.py:fine_solve',
        geometry={'fine_side': FINE_SIDE, 'coarse_side_for_training_surrogate': 32, 'frequency_Hz': 400e6,
                  'width_m': 2., 'transmitters': 10, 'receivers': 20},
        independent_objects={'train': TRAIN_COUNT, 'validation': VALIDATION_COUNT},
        generating_peak_range=[3., GENERATING_CAP], train_file_keys=['Y', 'radius', 'noise_fraction', 'split', 'object_id'],
        gt_handling='Generating contrast and clean fields are transient in worker memory and absent from training artifacts.',
        noise='Object-specific complex white direction shared across five noise levels, with 1e-6 relative numerical allowance.',
        allowance_scope='Numerical training data only; no interval forward enclosure is claimed.',
        worker_count=workers, host=platform.platform(), fine_solver_rtol=SOLVE_RTOL, rows=list(receipts),
        all_fine_solve_rhs=sum(r['forward_rhs'] for r in receipts),
        all_fine_gmres_iterations=sum(r['fine_gmres_iterations'] for r in receipts),
        all_fine_fft_green_rhs=sum(r['fft_green_rhs'] for r in receipts),
        all_worker_cpu_s=sum(r['cpu_s'] for r in receipts),
        export_wall_s=perf_counter()-write_start, wall_s=perf_counter()-t0, controller_cpu_s=process_time()-cpu0)
    atomic_json(out / 'generation_receipt.json', receipt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=STUDY / 'data')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    generate(args.out, args.workers)
