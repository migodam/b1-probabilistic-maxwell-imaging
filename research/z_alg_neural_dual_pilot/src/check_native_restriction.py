"""Compare independent native grid restrictions and full interval scores on training inputs."""
from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import LinAlgWarning, eigvalsh, solve
from threadpoolctl import threadpool_limits
import warnings

from native import Native, hermitian
from native_check import training_profiles
from projected_upper import NativeRestriction, PAPER

PILOT = Path(__file__).resolve().parents[1]
TILING = 16
CHECKED_PROFILES = (0, 1, 3, 4, 5)


def maximum_error(a, b):
    return float(np.max(np.abs(a - b)))


def construct_reference():
    """Build a literal pixel indicator basis and apply the original FFT operator."""
    with np.load(PAPER/'synthetic/public/fine_geometry.npz', allow_pickle=False) as f:
        symbol, sensor, incident = (f[k] for k in ('spectrum', 'PW360_10x20__S', 'PW360_10x20__B'))
    n = int(np.sqrt(incident.shape[0]))
    ry, rx = np.indices((n, n))
    q = n//TILING
    tile_id = ((ry//q)*TILING + rx//q).ravel()
    k = TILING**2
    v = (tile_id[:, None] == np.arange(k)[None, :]) / q
    gv = np.empty(v.shape, dtype=np.complex128)
    for pos in range(0, k, 16):
        sub = v[:, pos:pos+16].reshape(n, n, -1)
        padded = fft2(sub, s=symbol.shape, axes=(0, 1))
        gv[:, pos:pos+16] = ifft2(symbol[:, :, None]*padded, axes=(0, 1))[:n, :n].reshape(n*n, -1)
    return dict(v=v, gv=gv, s=sensor@v, b=v.T@incident, ry=ry.ravel(), rx=rx.ravel(),
                n=n, trial_dim=k, native_fft_green_rhs=k)


def restricted_distance(k, l, y):
    nu = int(np.count_nonzero(eigvalsh(k, check_finite=False) < 0))
    with warnings.catch_warnings():
        warnings.simplefilter('error', LinAlgWarning)
        x = solve(k, l, assume_a='her', check_finite=False)
    eig = eigvalsh(hermitian(y.conj().T@y-l.conj().T@x), check_finite=False)[::-1]
    retained = np.maximum(eig[nu:], 0.)
    return float(np.sqrt(retained.sum())), float(retained.max(initial=0.)), nu


def validate_one(reference, restriction, native, case, y, radius, u, mu, lam):
    n, v, gv = reference['n'], reference['v'], reference['gv']
    w = int(np.sqrt(len(mu)))
    span = n//w
    material_id = (reference['ry']//span)*w + reference['rx']//span
    z_pixels = (u*mu[material_id] + 1j*lam[material_id])/2
    right = v.T @ (z_pixels[:, None]*gv)
    k0 = hermitian(reference['s'].conj().T@reference['s'] + v.T@(mu[material_id, None]*v)
                   - right - right.conj().T)
    # Each tile indicator has one constant multiplier. Multiplication here is
    # on the native support before applying the independently assembled basis.
    with np.load(PAPER/'synthetic/public/fine_geometry.npz', allow_pickle=False) as data:
        incident = data['PW360_10x20__B']
    l0 = reference['s'].conj().T @ y + v.T@(z_pixels[:, None]*incident)

    ty, tx = np.indices((TILING, TILING))
    side_tiles = TILING//w
    idx = ((ty//side_tiles)*w + tx//side_tiles).ravel()
    z = (u*mu[idx] + 1j*lam[idx])/2
    zg = z[:, None]*restriction.mid_g
    k1 = hermitian(restriction.ss + np.diag(mu[idx]) - zg - zg.conj().T)
    l1 = restriction.mid_s.conj().T@y + z[:, None]*restriction.mid_b
    orient = dict(K_max_entry_error=maximum_error(k0,k1), L_max_entry_error=maximum_error(l0,l1))
    if max(orient.values()) > 5e-12:
        raise ArithmeticError('Projected pencil has inconsistent indexing or orientation.')

    numerical_d, top, restricted_nu = restricted_distance(k0, l0, y)
    diagnostic_radius = max(radius, 1.05*np.sqrt(y.shape[1]*top)+1e-10)
    proof = native.score(y,u,mu,lam,radius=radius)
    nominal_upper = restriction.check(y,u,mu,lam,radius)
    diagnostic_upper = restriction.check(y,u,mu,lam,diagnostic_radius)
    if proof['upper'] is None or diagnostic_upper['upper'] is None:
        raise ArithmeticError(f'Expected an overlapping pair of resolved checks for {case}.')
    if proof['lower'] > numerical_d + 2e-9 or proof['lower'] > diagnostic_upper['upper']:
        raise ArithmeticError(f'Restricted value is below the rigorous full lower bound for {case}.')
    if nominal_upper['upper'] is not None and proof['lower'] > nominal_upper['upper']:
        raise ArithmeticError(f'Nominal upper bound conflicts with the rigorous full lower bound for {case}.')
    return dict(case=case, cap=u, block_side=w, orientation=orient,
                actual_training_radius=radius, diagnostic_radius=diagnostic_radius,
                numerical_restricted_distance=numerical_d,
                numerical_restricted_top_retained_source_eigenvalue=top,
                numerical_restricted_state_inertia=restricted_nu,
                native_full=proof, restricted_at_actual_radius=nominal_upper,
                restricted_at_diagnostic_radius=diagnostic_upper)


def main():
    wall,cpu=perf_counter(),process_time()
    rows=[]
    reference_start, reference_cpu=perf_counter(),process_time()
    with threadpool_limits(limits=1, user_api='blas'):
        reference=construct_reference()
    ref_setup=dict(wall_s=perf_counter()-reference_start,cpu_s=process_time()-reference_cpu,
                   native_fft_green_rhs=reference['native_fft_green_rhs'])
    restriction=NativeRestriction(TILING)
    expected_g = reference['v'].T @ reference['gv']
    geometry_error=dict(G_max_entry_error=maximum_error(expected_g,restriction.mid_g),
                        S_max_entry_error=maximum_error(reference['s'],restriction.mid_s),
                        B_max_entry_error=maximum_error(reference['b'],restriction.mid_b),
                        column_basis_orthogonality_error=maximum_error(reference['v'].T @ reference['v'],np.eye(TILING**2)))
    if max(geometry_error.values()) > 5e-12:
        raise ArithmeticError('Projected geometry does not match direct native FFT actions.')
    native=Native()
    for i,(case,y,r,u,mu,lam) in enumerate(training_profiles()):
        if i not in CHECKED_PROFILES:
            continue
        with threadpool_limits(limits=1, user_api='blas'):
            row=validate_one(reference,restriction,native,case,y,r,u,mu,lam)
        rows.append(row)
        print(case,u,row['native_full']['lower'],row['numerical_restricted_distance'],
              row['restricted_at_actual_radius']['status'],row['restricted_at_diagnostic_radius']['upper'],flush=True)

    result=dict(scope='dualpilot_train_000, training only, fixed 2% noise. Profiles in native 96² units.',
                note='The actual radius is checked first. A second diagnostic radius tests the method even when '
                     'the native distance is physically above the actual radius or the trial space is too loose.',
                tiling=TILING, independent_reference_setup=ref_setup, restricted_interval_setup=restriction.setup,
                geometry_entry_errors=geometry_error, full_interval_calls=len(rows), restricted_radius_calls=2*len(rows),
                independent_dense_trial_K_solves=len(rows), independent_projected_rows=rows,
                native_cache=native.cache_summary(), wall_s=perf_counter()-wall,cpu_s=process_time()-cpu)
    (PILOT/'results/native_restriction_training_check.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
