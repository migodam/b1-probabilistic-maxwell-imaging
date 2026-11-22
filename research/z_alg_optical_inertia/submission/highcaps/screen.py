#!/usr/bin/env python3
"""Find shared 10-dimensional hidden-current subspaces for local high-cap tests."""
import itertools
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import svd
from scipy.sparse.linalg import LinearOperator, eigsh
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
FINE = STUDY / 'synthetic/public/fine_geometry.npz'
SIDE = 96
N = SIDE**2
APERTURE = 'PW360_10x20'
BITS = np.array(list(itertools.product((-1, 1), repeat=4)))
VERTEX_SIGNS = np.array(list(itertools.product(BITS, BITS)))
BASE_TS = (-4., -2., -1., 0., 1., 2., 4., 8., 16.)
TARGET_BOXES = ((.125, .125, -4., 16.), (.5, .06, 0., 0.))


def hermitian(a):
    return (a + a.swapaxes(-1, -2).conj()) / 2.


def green(v, spec):
    values = np.asarray(v).reshape(SIDE, SIDE, -1)
    return ifft2(fft2(values, s=spec.shape, axes=(0, 1)) * spec[..., None],
                 axes=(0, 1))[:SIDE, :SIDE].reshape(N, -1)


def blocks():
    grid = np.arange(N).reshape(SIDE, SIDE)
    split = SIDE//2
    return tuple(grid[y:y+split, x:x+split].ravel() for y, x in itertools.product((0, split), repeat=2))


def forms(v, gv):
    a = np.stack([v[i].conj().T @ v[i] for i in blocks()])
    c = np.stack([v[i].conj().T @ gv[i] for i in blocks()])
    return hermitian(a), hermitian(c), hermitian(c/1j)


def propose(spec, receive, t):
    _, _, basis = svd(receive, full_matrices=False)
    basis = basis.conj().T

    def project(z):
        return z - basis @ (basis.conj().T @ z)

    def action(v):
        pv = project(v)
        return project(pv + green(pv, -spec.real + t*spec.imag).ravel())

    matrix = LinearOperator((N, N), matvec=action, dtype=complex)
    initial = project(np.random.default_rng(613).normal(size=N))
    _, vectors = eigsh(matrix, which='SA', k=10, ncv=35, tol=1.e-10, v0=initial)
    v = np.linalg.qr(np.column_stack([project(c) for c in vectors.T]))[0]
    return v, forms(v, green(v, spec))


def max_eigenvalue(matrices, real_width, optical_width, t_lo, t_hi):
    a, h, o = matrices
    real_choices = 1 + real_width*VERTEX_SIGNS[:, 0]
    optical_jitter = optical_width*VERTEX_SIGNS[:, 1]
    mx = float('-inf')
    for u, t in itertools.product((1., 1.5), (t_lo, t_hi)):
        corners = np.einsum('vb,bij->vij', real_choices, a-u*h)
        corners += np.einsum('vb,bij->vij', t+optical_jitter, o)
        mx = max(mx, float(np.linalg.eigvalsh(corners)[:, -1].max()))
    return mx


def screen():
    start = perf_counter(), process_time()
    with np.load(FINE, allow_pickle=False) as data:
        spec, receive = data['spectrum'], data[APERTURE+'__S']
    proposals = [propose(spec, receive, t) for t in BASE_TS]
    configs = [
        (w, e, a, b) for w, e in (
            (.02, .02), (.06, .06), (.125, .125), (.25, .25), (.5, .5), (.75, .75), (.875, .875),
            (.5, .06), (.875, .06), (.98, .02)
        ) for a, b in ((-1., 1.), (-2., 2.), (-4., 4.), (0., 8.), (0., 16.), (-4., 16.),
                       (0., 4.), (0., 2.), (.25, .75), (0., 0.))
    ]
    hidden = {}
    for cfg in configs:
        scores = [max_eigenvalue(m, *cfg) for v, m in proposals]
        idx = int(np.argmin(scores))
        hidden[' '.join(map(str, cfg))] = dict(args=cfg, score=scores[idx], center=BASE_TS[idx])
    centers = {hidden[' '.join(map(str, cfg))]['center'] for cfg in TARGET_BOXES}
    np.savez_compressed(HERE/'subspace_proposals.npz', **{
        f't_{t:g}': v for t, (v, _) in zip(BASE_TS, proposals, strict=True) if t in centers
    })
    info = dict(aperture=APERTURE, state_cells=N, centers=BASE_TS, candidates=hidden,
                wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])
    (HERE/'screen.json').write_text(json.dumps(info, indent=2)+'\n')
    for cfg in ((.25,.25,-4.,16.), (.5,.5,-4.,16.), (.875,.875,-4.,16.),
                (.98,.02,-4.,16.), (.25,.25,-1.,1.), (.5,.5,-1.,1.), (.875,.875,-1.,1.)):
        print(cfg, hidden[' '.join(map(str, cfg))], flush=True)


if __name__ == '__main__':
    with threadpool_limits(1):
        screen()
