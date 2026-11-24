"""Explore measured row-space minima, and trace-efficient projected rays."""

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, norm, svd
from scipy.sparse.linalg import eigsh, LinearOperator
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "submission/single"))
from local_moment import CASES, CAP, N, NumericGeometry, blocks, h, load


def approx_forms(g, x, b=None):
    e = g.green(x)
    if b is not None:
        e += b
    q, t = [], []
    for ix in g.indices:
        xp = x[ix]
        xe = xp.conj().T @ e[ix]
        q.append(xp.conj().T @ xp - CAP * h(xe))
        t.append((xe - xe.conj().T) / (2j))
    return np.array(q), np.array(t)


def eig_bnds(a):
    ev = np.linalg.eigvalsh(a)
    return ev[..., [0,-1]]


def eigen_directions(g, operator, k, which):
    def action(x):
        return g.project_off(operator(g.project_off(x, g.v)), g.v)
    ev, x = eigsh(LinearOperator((N,N), matvec=action, dtype=complex), k=k,
                  tol=1e-8, which=which, v0=np.random.default_rng(1241).normal(size=N))
    return ev, x


def reactive_rays(g):
    def real_green(x):
        return g.green(x, g.symbol.real)
    return eigen_directions(g, real_green, 18, "LA")


def radiation_rays(g):
    def imag_green(x):
        return g.green(x, g.symbol.imag)
    return eigen_directions(g, imag_green, 10, "LA")


def fit_rows(g, y):
    return g.v @ ((g.u.conj().T @ y) / g.sigma[:, None])


def run():
    s, b, symbol, data, _ = load()
    g = NumericGeometry(s,b,symbol,blocks(2))
    print('S singular values', g.sigma)
    bp = g.project_off(b, g.v)
    print('Bker singulars', svd(bp, compute_uv=False))
    evr, raysr = reactive_rays(g)
    evg, raysg = radiation_rays(g)
    print("reactive", evr, 'rays q', approx_forms(g,raysr)[0].diagonal(axis1=1,axis2=2).real)
    print("radiative", evg, 'rays q', approx_forms(g,raysg)[0].diagonal(axis1=1,axis2=2).real)
    print("old rays total optical", approx_forms(g, np.load(HERE.parents[1]/'submission/single/local_moment.npz')['rays'])[1].sum(axis=0).diagonal().real)
    for key,y in data.items():
        print('case',key,'||Y||_2',norm(y,2))
        xp = fit_rows(g, y)
        for x in [xp, np.load(HERE.parents[1]/'submission/single/local_moment.npz')[key]]:
            q,t=approx_forms(g,x,b)
            print('||x||²',norm(x,2)**2,'q t',eig_bnds(q),eig_bnds(t),'joint q t',eig_bnds(q.sum(axis=0)),eig_bnds(t.sum(axis=0)))
        print('visibility of optical gradient', svd(g.project_off(g.green(xp, g.symbol.imag)-.5j*b,g.v),compute_uv=False))


if __name__ == '__main__':
    with threadpool_limits(1):
        run()
