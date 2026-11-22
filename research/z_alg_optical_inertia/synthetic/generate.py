#!/usr/bin/env python3
"""Generate six fresh SRAP-physics objects and hold observation truth apart."""
import os

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import argparse
from pathlib import Path
import platform
import sys
import time

import numpy as np
from scipy.linalg import solve

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
sys.path.insert(0, str(WORKSPACE / 'research/srap_cie_array_study_20261121'))
from data2d import (FREQUENCY, SIDE, WIDTH, SyntheticOperators, area_average, ellipse,
                    history, rotated_xy, write_json, write_npz)
from research.research_max_struct_r8.src.paper_physics import acquisition

FINE = 96
CAP = 4.0
SHAPE_SEED = 20261122091
NOISE_SEED = 20261122092
ARRAYS = ('PW360_10x20', 'PW180_10x20')
FRACTIONS = (0.10, 0.02)
FAMILIES = ('forked_block', 'quadrupole_lobes', 'terraced_chevrons',
            'crescent_and_plate', 'staggered_channels', 'graded_pinwheel')


def material(centers, index):
    rng = np.random.default_rng(np.random.SeedSequence([SHAPE_SEED, index]))
    angle = rng.uniform(-0.32, 0.32)
    shift = rng.uniform(-0.10, 0.10, 2)
    scale = rng.uniform(0.88, 1.08, 2)
    contrast = rng.uniform(3.25, CAP)
    x, y = (rotated_xy(centers, angle, shift) / scale).T
    if index == 0:
        block = (np.abs(x) <= 0.53) & (np.abs(y) <= 0.48)
        slots = (y > -0.12) & ((np.abs(x - 0.24) < 0.075) | (np.abs(x + 0.24) < 0.075))
        fraction = (block & ~slots) * (0.86 + 0.14 * (y >= 0))
    elif index == 1:
        fraction = np.maximum.reduce([
            ellipse(x, y, -0.36, -0.27, 0.30, 0.24),
            0.82 * ellipse(x, y, 0.26, -0.32, 0.22, 0.29),
            0.91 * ellipse(x, y, -0.24, 0.33, 0.26, 0.29),
            0.78 * ellipse(x, y, 0.37, 0.24, 0.32, 0.19),
            0.68 * ellipse(x, y, 0.04, 0.01, 0.23, 0.22),
        ])
    elif index == 2:
        bottom = history.stroke(x, y, [[-0.60, -0.23], [-0.25, -0.41],
                                       [0.13, -0.28], [0.52, -0.10]], 0.16)
        upper = history.stroke(x, y, [[-0.43, 0.32], [-0.09, 0.11],
                                      [0.26, 0.31], [0.53, 0.43]], 0.14)
        fraction = np.maximum(bottom, 0.82 * upper)
    elif index == 3:
        crescent = ellipse(x, y, -0.12, -0.03, 0.6, 0.49) & ~ellipse(x, y, 0.13, 0.12, 0.47, 0.38)
        plate = (np.abs(x - 0.45 + 0.35 * y) < 0.14) & (np.abs(y - 0.03) < 0.43)
        fraction = np.maximum(crescent, 0.76 * plate)
    elif index == 4:
        slab = (np.abs(x) < 0.59) & (np.abs(y + 0.02) < 0.40)
        cut = ((np.abs(y - 0.17) < 0.065) & (x < 0.28)) | ((np.abs(y + 0.18) < 0.072) & (x > -0.26))
        fraction = (slab & ~cut) * (0.85 + 0.15 * np.cos(2.5 * x + 1.1 * y) ** 2)
    else:
        radial = np.hypot(x, y)
        azimuth = np.arctan2(y, x)
        blades = (radial < 0.66) & ((radial < 0.18) | (np.cos(4 * azimuth - 5 * radial) > 0.2))
        fraction = blades * (0.72 + 0.28 * (1 - np.minimum(radial / 0.66, 1)))
    return contrast * fraction


def norm_up(x, axis=None):
    """Pad the measured binary64 array norm before storing a public radius."""
    values = np.asarray(x)
    count = values.size if axis is None else values.shape[axis]
    return np.nextafter(np.linalg.norm(values, axis=axis) * (1 + 8 * count * np.finfo(float).eps), np.inf)


def generate(out):
    wall0, cpu0 = time.perf_counter(), time.process_time()
    public, evaluation = out / 'public', out / 'evaluation'
    for folder in (public, evaluation):
        folder.mkdir(parents=True, exist_ok=True)
    ops = SyntheticOperators(FINE, ARRAYS)
    G = acquisition(SIDE, FREQUENCY)['G']
    write_npz(public / 'coarse_green.npz', G=G, centers=ops.coarse_centers)
    manifest = dict(contract='synthetic/FREEZE.md', coarse_side=SIDE, fine_side=FINE,
                    width_m=WIDTH, frequency_hz=FREQUENCY, upper_contrast=CAP,
                    green_file='coarse_green.npz', geometries={}, observations=[])
    receipt = dict(host=platform.platform(), srapphysics_setup=ops.cost, objects={})
    for array in ops.arrays:
        rx, tx = ops.rx_ids[array.name], ops.tx_ids[array.name]
        file = array.name.lower() + '.npz'
        write_npz(public / file, S=ops.coarse_s[rx], B=ops.coarse_e[:, tx],
                  rx_m=array.receiver_points, tx_plane_wave_deg=array.tx[:, 0])
        manifest['geometries'][array.name] = dict(file=file, green_file='coarse_green.npz',
                                                  rx_count=len(rx), tx_count=len(tx))
    for index, family in enumerate(FAMILIES):
        started, cpu = time.perf_counter(), time.process_time()
        object_id = f'optfresh_{index+1:02}'
        m_fine = material(ops.fine_centers, index)
        m = area_average(m_fine, FINE, SIDE)
        if not 0 <= m_fine.min() <= m_fine.max() <= CAP:
            raise ValueError('material outside frozen contrast range')
        clean, fine_stats = history.fine_solve(m_fine, ops.fine_e, ops.fine_s,
                                               ops.green.apply, ops.green.diagonal, rtol=1e-10)
        active = np.flatnonzero(m)
        block = G[np.ix_(active, active)]
        A = np.eye(len(active)) - m[active, None] * block
        rhs = m[active, None] * ops.coarse_e[active]
        current = solve(A, rhs, check_finite=False)
        proxy = ops.coarse_s[:, active] @ current
        born = ops.coarse_s[:, active] @ rhs
        receipt['objects'][object_id] = dict(fine_stats, family=family,
            fine_max_contrast=float(m_fine.max()), fine_nonzero_fraction=float(np.mean(m_fine > 0)),
            coarse_proxy_vs_coarse_born=float(np.linalg.norm(proxy-born)/np.linalg.norm(proxy)),
            coarse_proxy_state_residual=float(np.linalg.norm(A @ current-rhs)/np.linalg.norm(rhs)),
            coarse_proxy_rhs=int(rhs.shape[1]), coarse_proxy_active_cells=len(active),
            wall_s=time.perf_counter()-started, cpu_s=time.process_time()-cpu)
        write_npz(evaluation / (object_id + '.npz'), m=m, m_fine=m_fine, family=np.array(family))
        for ai, array in enumerate(ops.arrays):
            ids = np.ix_(ops.rx_ids[array.name], ops.tx_ids[array.name])
            clean_view, proxy_view = clean[ids], proxy[ids]
            for ni, fraction in enumerate(FRACTIONS):
                obs_id = f'{object_id}__{array.name.lower()}__n{int(fraction*100):02}'
                file = obs_id + '.npz'
                seed = np.random.SeedSequence([NOISE_SEED, index, ai, ni])
                y = clean_view + history.complex_noise(clean_view, fraction, seed)
                noise = y - clean_view
                stats = dict(id=obs_id, object_id=object_id, geometry=array.name, file=file,
                             noise_fraction=fraction, upper_contrast=CAP,
                             radius_noise_fro=float(norm_up(noise)),
                             radii_noise_by_source=norm_up(noise, axis=0).tolist())
                write_npz(public / file, Y=y)
                write_npz(evaluation / file, Y_fine_clean=clean_view, Y_coarse_proxy=proxy_view,
                          discrepancy_fro=norm_up(clean_view-proxy_view),
                          discrepancy_by_source=norm_up(clean_view-proxy_view, axis=0),
                          proxy_error_fro=norm_up(y-proxy_view), proxy_error_by_source=norm_up(y-proxy_view, axis=0),
                          proxy_relative_discrepancy=np.linalg.norm(clean_view-proxy_view)/np.linalg.norm(clean_view))
                manifest['observations'].append(stats)
        print(f'{object_id}: {family}, {receipt["objects"][object_id]["wall_s"]:.2f} s', flush=True)
    write_json(public / 'manifest.json', manifest)
    receipt['wall_s'] = time.perf_counter() - wall0
    receipt['cpu_s'] = time.process_time() - cpu0
    write_json(evaluation / 'generation_receipt.json', receipt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE)
    generate(parser.parse_args().out)
