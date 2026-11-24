"""Single post-seal existence test for the preassigned half-ring observation."""

import argparse
from copy import copy
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from run import DEFAULT_ROOT, DEFAULT_SEAL, HERE, STUDY, NativeForward, elapsed, nearest_lift, timers, write_json
from research.z_alg_xudong_cie.holdout.data import (
    NATIVE_SIDE, get_native_case, public_geometry, require_freeze, write_npz,
)

sys.path.insert(0, str(STUDY/'minpeak_native'))
from verify import ExactFourierGeometry, ExactGeometry, ball, hi, lo, norm, verified_ports
from gn_recovery import gauss_newton
from search import Model

OBJECT_ID = 'b3m_0001'
OBSERVATION = f'{OBJECT_ID}__pw180_10x20__n10'
CAP = 20.
STEPS = 6


def adjudicate(exact, fft, m, y, r):
    start = timers()
    try:
        f, e, work = verified_ports(exact, fft, m)
        distance = norm(f-ball(y))
        lower = max(0., lo(distance.lower()-e))
        upper = hi(distance.upper()+e)
        status = 'inside' if upper <= r else ('outside' if lower > r else 'unresolved')
        return dict(status=status, residual_lower=lower, residual_upper=upper,
                    admission_margin_lower=lo(arb(r)-arb(upper)), inverse_and_state=work, **elapsed(start))
    except (ArithmeticError, ValueError) as exc:
        return dict(status='unresolved', error=str(exc), **elapsed(start))


def verify_saved(root, output, m0, m, p, y, r):
    verify_start = timers()
    ctx.prec, ctx.threads = 128, 1
    exact = ExactFourierGeometry(public_geometry(root, NATIVE_SIDE))
    exact.S, exact.B = p['S'], p['E']
    fft = ExactGeometry(p['S'], p['E'], p['spectrum'], [])
    proof = dict(observation_id=OBSERVATION, interval_bits=ctx.prec, padded_radius=r,
                 original_projection_archive='projections.npz', original_projection_key=f'{OBJECT_ID}__m32',
                 corrected_archive='witness.npz', corrected_key='m32',
                 original_geometry_preparation=elapsed(verify_start))
    proof['projection'] = adjudicate(exact, fft, m0, y, r)

    # Reuse the older interval validator's internal cap of 4 without changing
    # the physical equation. m'(B'+G'J) = m(B+GJ) for the exact binary scale.
    scale_start = timers()
    scale = next(c for c in (1., 2., 4., 8.) if m.max() <= c*4.)
    scaled_m, scaled_b, scaled_spec = m/scale, p['E']*scale, p['spectrum']*scale
    if not (np.array_equal(scale*scaled_m, m) and np.array_equal(scaled_b/scale, p['E'])
            and np.array_equal(scaled_spec/scale, p['spectrum'])):
        raise ArithmeticError('Binary state rescaling must round-trip exactly.')
    scaled_geometry, scaled_fft = exact, fft
    if scale != 1.:
        scaled_geometry = copy(exact)
        scaled_geometry.B = scaled_b
        scaled_geometry.exact = [scale*v for v in exact.exact]
        scaled_geometry.centers = scale*exact.centers
        scaled_geometry.error_max = exact.error_max*scale
        scaled_geometry.kernel_center_max = exact.kernel_center_max*scale
        scaled_fft = ExactGeometry(p['S'], scaled_b, scaled_spec, [])
    proof['power_of_two_state_rescaling'] = dict(scale=scale, exact_round_trip=True, **elapsed(scale_start))
    proof['corrected'] = adjudicate(scaled_geometry, scaled_fft, scaled_m, y, r)
    write_json(output/'strict_witness.json', {**proof, **elapsed(verify_start)})
    print('STRICT', [(label, proof[label]['status'], proof[label].get('residual_upper'))
                     for label in ('projection', 'corrected')], flush=True)

    clipped_start = timers()
    clipped = np.minimum(m, 4.)
    native = NativeForward(p['spectrum'], NATIVE_SIDE)
    ports, forward = native.forward(clipped, p['E'], p['S'])
    distance = float(np.linalg.norm(ports-y))
    write_npz(output/'witness_cap4.npz', m32=clipped.reshape(32,3,32,3)[:,0,:,0].ravel(),
              m96=clipped, numeric_ports=ports, observation=np.array(OBSERVATION))
    clipped_receipt = dict(candidate_archive='witness_cap4.npz', peak=float(clipped.max()),
                           clipped_native_cells=int(np.count_nonzero(m!=clipped)),
                           numeric_residual=distance, numeric_forward=forward)
    if distance < r:
        clipped_receipt['strict'] = adjudicate(exact, fft, clipped, y, r)
    else:
        clipped_receipt['strict'] = dict(status='not_attempted_for_numeric_outside')
    proof['clipped_cap4'] = {**clipped_receipt, **elapsed(clipped_start)}
    write_json(output/'strict_witness.json', {**proof, **elapsed(verify_start)})
    print('CLIPPED', distance, clipped_receipt['strict']['status'], flush=True)


def run(root, output, seal, verify_only):
    start = timers()
    stamp = datetime.now(timezone.utc).isoformat()
    seal = require_freeze(seal)
    _, p = get_native_case(OBSERVATION, root)
    y, r = p['Y'], float(p['delta'])
    with np.load(output/'projections.npz', allow_pickle=False) as z:
        m0 = z[f'{OBJECT_ID}__m96']
        numeric_projection = float(np.linalg.norm(z[f'{OBJECT_ID}__PW180_10x20_prediction']-y))
    if verify_only:
        with np.load(output/'witness.npz', allow_pickle=False) as z:
            if not (str(z['observation']) == OBSERVATION and np.array_equal(y, z['Y'])
                    and r == float(z['padded_radius'])):
                raise ArithmeticError('The archived witness uses a different measurement or radius.')
            m = nearest_lift(z['m32'])
            if not (np.array_equal(m, z['m96']) and np.all((0. <= m) & (m <= CAP))):
                raise ArithmeticError('The archived witness is outside the required lift or box.')
        verify_saved(root, output, m0, m, p, y, r)
        return

    model = Model(p['spectrum'], p['S'], p['E'], y)
    candidate, history = gauss_newton(model, m0, 32, np.arange(32**2), CAP, STEPS, r)
    m = candidate['m']
    m32 = m.reshape(32, 3, 32, 3)[:, 0, :, 0].ravel()
    if not (np.array_equal(m, nearest_lift(m32)) and np.all((0. <= m32) & (m32 <= CAP))):
        raise ArithmeticError('The existence witness must lie in the frozen box and lift.')
    state_relative = float(np.linalg.norm(candidate['j']-m[:, None]*(p['E']+model.green(candidate['j']))) /
                           np.linalg.norm(m[:, None]*p['E']))
    write_npz(output/'witness.npz', m32=m32, m96=m, numeric_ports=candidate['ports'],
              observation=np.array(OBSERVATION), Y=y, padded_radius=np.array(r))
    payload = dict(scope='Offline existence only; same upper material box and frozen native observation',
                   online_seal=str(seal), observation_id=OBSERVATION, started_utc=stamp,
                   assignment='Smallest anonymous ID, half-ring, 10%, fixed before this correction',
                   archive='witness.npz', cap=CAP, maximum_GN_steps=STEPS, numerical_threads=1,
                   numeric_projection_residual=numeric_projection, padded_radius=r,
                   numeric_solution_residual=float(candidate['residual']),
                   numeric_state_relative_residual=state_relative, correction=history,
                   native_algorithm_counts=dict(model.count), floating_stage=elapsed(start))
    write_json(output/'witness_search.json', payload)
    print('FLOAT', OBSERVATION, 'initial', numeric_projection, 'final', candidate['residual'],
          'radius', r, 'work', dict(model.count), flush=True)

    verify_saved(root, output, m0, m, p, y, r)
    write_json(output/'witness_search.json', {**payload, 'finished_utc': datetime.now(timezone.utc).isoformat(),
                                            'strict_observations': 'strict_witness.json', **elapsed(start)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT)
    parser.add_argument('--output', type=Path, default=HERE)
    parser.add_argument('--online-seal', type=Path, default=DEFAULT_SEAL)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.root, args.output, args.online_seal, args.verify_only)
