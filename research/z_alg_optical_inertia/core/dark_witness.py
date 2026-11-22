"""Arb exact-nullspace witness at the full-ring fresh material comparison cap."""
import argparse
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh, qr
from threadpoolctl import threadpool_limits

from certificate import ball, frobenius_squared, safe_float, star
from inertia import adjoint, clock, hermitian, timing
from scan import RESULTS, load_cases, read_json, save_json


def negative_witness(p, u, ratio, precision=256, return_proposal=False):
    started = clock()
    Q, _ = qr(adjoint(p.S), mode='full', check_finite=False)
    N = Q[:, p.r:]
    C = np.eye(p.n) - u * p.Gh + ratio * p.Ga
    small = hermitian(adjoint(N) @ (C @ N))
    value, vec = eigh(small, subset_by_index=(0, 0), check_finite=False)
    ctx.prec = precision
    proposal = N @ vec
    G, S, approximate = map(ball, (p.G, p.S, proposal))
    St = star(S)
    exact_null = approximate - St * (S * St).solve(S * approximate)
    norm_squared = frobenius_squared(exact_null)
    xGx = (star(exact_null) * (G * exact_null))[0, 0]
    result = norm_squared - arb(u) * xGx.real + arb(ratio) * xGx.imag
    upper = safe_float(result, upper=True)
    receipt = {'status': 'verified' if upper < 0. else 'inconclusive',
            'cap': u, 'lambda_over_mu': ratio, 'precision_bits': precision,
            'dark_form_upper': upper, 'dark_form_numeric': float(value[0]),
            'dark_norm_squared_lower': safe_float(norm_squared, upper=False),
            'dark_radiation_energy_upper': safe_float(xGx.imag, upper=True),
            'dark_radiation_energy_lower': safe_float(xGx.imag, upper=False),
            'dark_G_H_quadratic_lower': safe_float(xGx.real, upper=False),
            'receiver_Gram_order': p.r, **timing(started)}
    return (receipt, proposal) if return_proposal else receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULTS / 'fresh_dark_witness.json')
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api='blas'):
        index = 'optfresh_01__pw360_10x20__n02'
        p, *_ = next(iter(load_cases('synthetic').values()))[index]
        receipts = [negative_witness(p, .75, t) for t in (0., 10., 100., 1000., 10_000.)]
        result = {'geometry_representative': index, 'same_G_and_S_as': [
            'optfresh_05__pw360_10x20__n02', 's2d_007__pw360_10x20'],
            'certificates': receipts, 'scope': 'chosen ratios as binary64; for each ratio, every mu>0'}
        save_json(args.output, result)
        print(result, flush=True)


if __name__ == '__main__':
    main()
