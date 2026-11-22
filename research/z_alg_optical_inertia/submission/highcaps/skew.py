#!/usr/bin/env python3
"""Numerically check selected more-skewed quadrants at the cap 1."""
import json
from time import perf_counter, process_time

import numpy as np
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh
from threadpoolctl import threadpool_limits

from screen import APERTURE, FINE, HERE, N, blocks, green

WEIGHTS = (
    (1.5, .5, .5, .5), (1.5, .5, 1.5, .5), (.45, 1.55, .68, 1.32),
    (1., .3, 1., .3), (1., 1., .1, .1), (1., .1, .1, .1),
)


def main():
    start = perf_counter(), process_time()
    with np.load(FINE, allow_pickle=False) as inputs:
        S, spectrum = inputs[APERTURE+'__S'], inputs['spectrum']
    input_weights = np.empty(N)
    output = []
    for vector in WEIGHTS:
        for region, value in zip(blocks(), vector, strict=True):
            input_weights[region] = value
        for t in (0., .4, 1.):
            mu = .001

            def apply(x):
                weighted_g = input_weights * green(x, spectrum).ravel()
                adjoint_g = green(input_weights*x, spectrum.conj()).ravel()
                optical = green(x, spectrum.imag).ravel()
                return S.conj().T @ (S @ x) / mu + input_weights*x - .5*(weighted_g + adjoint_g) + t*optical

            matrix = LinearOperator((N, N), dtype=complex, matvec=apply)
            v0 = np.random.default_rng(545).standard_normal(N)
            record = dict(w=vector, lambda_ratio_all_blocks=t, mu=mu, cap=1.)
            try:
                spectrum_low = eigsh(matrix, which='SA', k=12, ncv=32, tol=1.e-8, maxiter=350,
                                    v0=v0, return_eigenvectors=False)
            except ArpackNoConvergence as err:
                spectrum_low = err.eigenvalues
                record['reason'] = 'ARPACK iteration cap'
            record.update(lowest_values=np.sort(spectrum_low).tolist(),
                          negatives_returned=int(np.count_nonzero(spectrum_low < -1.e-8)))
            output.append(record)
            print(vector, t, 'numerical negatives returned', record['negatives_returned'], flush=True)
    info = dict(method='numerical 12 lowest eigenvalues of K/mu on native grid',
                w_quadrant_order='row-major', evidence='numerical check',
                cases=output, wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])
    (HERE/'skew.json').write_text(json.dumps(info, indent=2)+'\n')


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
