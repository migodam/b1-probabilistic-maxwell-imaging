#!/usr/bin/env python3
"""Score the matched stored physical truth after every online profile is complete."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys
from types import SimpleNamespace

from scan import ARCHIVE, HERE, STUDY, clock, load_public, read, timing, write

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

sys.path.insert(0,str(STUDY/'submission/diagnostic'))
from audit import dense_audit, hi, lo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,default=ARCHIVE/'sealed_truth/independent24',
                        help='Authorized local directory with original sealed material references.')
    args = parser.parse_args()
    if not args.reference.is_dir():
        raise SystemExit('Original R14 sealed reference directory not available. Use the bundled offline scalar scores.')

    frozen = read(HERE/'frozen.json')
    online = read(HERE/'verified.json')
    if len(online['receipts']) != online['jobs_planned']:
        raise SystemExit('The online run is still unfinished.')
    if read(HERE/'numerical_grid.json')['cases'] != frozen['cases']:
        raise ValueError('Online case list differs from the fixed selection.')
    write(HERE/'online_seal.json',
          dict(cases=frozen['cases'], stages=['numerical_grid.json','verified.json'],
               online_profiles=len(online['receipts']), utc=datetime.now(timezone.utc).isoformat(),
               status='ONLINE_PROFILE_RANKING_AND_INTERVAL_RESULTS_SEALED_BEFORE_SCORING',
               tuning_after_seal=False))
    started = clock()
    G, S, B, Y, radii = load_public(HERE/'public_inputs.npz', frozen)
    geometry = SimpleNamespace(G=G,S=S,E0=B)
    ctx.prec,ctx.threads=128,1
    rows = {}
    for case, y, radius in zip(frozen['cases'],Y,radii,strict=True):
        at = clock()
        with np.load(args.reference/f'{case}.npz',allow_pickle=False) as f:
            m = f['material_forward']
            m_coarse = f['material_coarse']
            clean = f['Y_clean']
            noise = f['noise']
            if int(f['forward_side'])!=24 or not np.array_equal(m,m_coarse):
                raise ValueError(f'Forward and inverse grids are not matched in sealed reference {case}.')
        physical = dense_audit(geometry,m,y,clean)
        physical['operator_to_saved_clean_upper'] = physical.pop('online_forward_error_hi')
        lo_res,hi_res = (physical[k] for k in ('full_wave_distance_lo','full_wave_distance_hi'))
        peak = max(arb(float(x)) for x in m)
        minimum = min(arb(float(x)) for x in m)
        delta = arb(float(radius))
        # No clipping or truth substitution is used to improve physical fit.
        rows[case] = dict(status='inside' if hi_res<=radius else 'outside' if lo_res>radius else 'unresolved',
                          public_CN99_radius=float(radius), truth_peak=float(peak),
                          peak_interval=[lo(peak),hi(peak)], truth_minimum=float(minimum),
                          generating_model_matched=True, full_wave=physical,
                          realized_sampled_noise_fro=float(np.linalg.norm(noise)),
                          realized_sample_over_radius=float(np.linalg.norm(noise)/radius),
                          true_residual_over_radius=[lo(arb(lo_res)/delta),hi(arb(hi_res)/delta)],
                          coverage_margin_lower=lo(delta-arb(hi_res)), **timing(at))
        print(case, rows[case]['status'], 'peak',float(peak), 'full-wave',lo_res,hi_res,
              'delta',float(radius),flush=True)
    write(HERE/'offline_scalar_scores.json',
          dict(scope='Authorized original R14 synthetic reference; scalar scores only, no sealed material or fields exported',
               precision_bits=ctx.prec, Arb_threads=ctx.threads,
               source_rhs=sum(r['full_wave']['source_rhs'] for r in rows.values()),
               receiver_rhs=sum(r['full_wave']['receiver_rhs'] for r in rows.values()),
               inverse_rhs=sum(r['full_wave']['inverse_rhs'] for r in rows.values()), rows=rows, **timing(started)))


if __name__=='__main__':
    with threadpool_limits(1):
        main()
