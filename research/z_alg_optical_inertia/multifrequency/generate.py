#!/usr/bin/env python3
"""Add an independently simulated 500 MHz half-ring on the same six objects."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import solve

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SYNTH = HERE.parent / 'synthetic'
sys.path.insert(0, str(SYNTH))
from generate import CAP, FAMILIES, FINE, FRACTIONS, SIDE, area_average, material, norm_up
import data2d as physics
from research.research_max_struct_r8.src.paper_physics import acquisition

FREQUENCY = 500e6
ARRAY = 'PW180_10x20'
NOISE_SEED = 20261122093


def main():
    t, c = time.perf_counter(), time.process_time()
    public, eval_dir = HERE / 'data/public', HERE / 'data/evaluation'
    for path in (public, eval_dir):
        path.mkdir(parents=True, exist_ok=True)
    physics.FREQUENCY = FREQUENCY
    ops = physics.SyntheticOperators(FINE, (ARRAY,))
    G = acquisition(SIDE, FREQUENCY)['G']
    B, S, Bf, Sf = ops.coarse_e, ops.coarse_s, ops.fine_e, ops.fine_s
    tx, rx = ops.tx_ids[ARRAY], ops.rx_ids[ARRAY]
    ids = np.ix_(rx, tx)
    np.savez_compressed(public / 'pw180_500MHz.npz', G=G, B=B[:, tx], S=S[rx],
                        rx_m=ops.arrays[0].receiver_points, tx_plane_wave_deg=ops.arrays[0].tx[:, 0])
    records, costs = [], {}
    for i, family in enumerate(FAMILIES):
        start, cpu = time.perf_counter(), time.process_time()
        object_id = f'optfresh_{i+1:02}'
        fine = material(ops.fine_centers, i)
        chi = area_average(fine, FINE, SIDE)
        clean, stats = physics.history.fine_solve(fine, Bf, Sf, ops.green.apply, ops.green.diagonal, rtol=1e-10)
        active = np.flatnonzero(chi)
        chi = chi[active]
        lhs = np.eye(len(active)) - chi[:, None]*G[np.ix_(active,active)]
        rhs = chi[:, None]*B[active]
        sol = solve(lhs, rhs, check_finite=False)
        proxy = (S[:, active] @ sol)[ids]
        clean = clean[ids]
        costs[object_id] = dict(stats, family=family, proxy_state_relative_residual=float(np.linalg.norm(lhs@sol-rhs)/np.linalg.norm(rhs)),
                               proxy_relative_discrepancy=float(np.linalg.norm(proxy-clean)/np.linalg.norm(clean)),
                               wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu)
        for j, fraction in enumerate(FRACTIONS):
            observation = f'{object_id}__pw180_10x20__n{int(100*fraction):02}'
            seed = np.random.SeedSequence([NOISE_SEED, i, j])
            Y = clean + physics.history.complex_noise(clean, fraction, seed)
            file = observation + '.npz'
            physics.write_npz(public / file, Y=Y)
            physics.write_npz(eval_dir / file, Y_fine_clean=clean, Y_coarse_proxy=proxy,
                              discrepancy_fro=norm_up(clean-proxy), proxy_error_fro=norm_up(Y-proxy))
            records.append(dict(id=observation, file=file, object_id=object_id, radius_noise_fro=float(norm_up(Y-clean)),
                                noise_fraction=fraction, upper_contrast=CAP))
    physics.write_json(public / 'manifest.json', dict(frequency_hz=FREQUENCY, paired_frequency_hz=400e6,
            side=SIDE, fine_side=FINE, geometry='pw180_500MHz.npz', observations=records))
    physics.write_json(eval_dir / 'cost.json', dict(fine_grid_setup=ops.cost, objects=costs,
                       wall_seconds=time.perf_counter()-t, cpu_seconds=time.process_time()-c))
    print('500 MHz:', len(records), 'observations', flush=True)


if __name__ == '__main__':
    main()
