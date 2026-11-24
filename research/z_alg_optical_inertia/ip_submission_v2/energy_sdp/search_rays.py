"""Form a shared low-energy uniform witness from four projected dark currents."""

import numpy as np
from scipy.linalg import eigh, norm, svd
from scipy.optimize import differential_evolution, minimize
from threadpoolctl import threadpool_limits

from probe import (approx_forms, blocks, fit_rows, load, NumericGeometry, reactive_rays,
                   radiation_rays, h, HERE, N)


def local_optical_gradient(g, v):
    ev = g.green(v)
    gradients = []
    for ix in g.indices:
        pv = np.zeros(N, complex)
        pev = np.zeros(N, complex)
        pv[ix], pev[ix] = v[ix], ev[ix]
        gradients.append((pev - g.green(pv, g.symbol.conj())) / (2j))
    return g.project_off(np.column_stack(gradients), g.v)


def max_ev(x):
    return np.linalg.eigvalsh(x)[..., -1]


def make_ray_span(g):
    e, v = reactive_rays(g)
    k = np.argmax(e)
    ref = v[:, k]
    if norm(ref.real) >= norm(ref.imag):
        ref = ref.real / norm(ref.real)
    else:
        ref = ref.imag / norm(ref.imag)
    ref = g.project_off(ref, g.v)
    ref /= norm(ref)
    opt_grad = local_optical_gradient(g, ref)
    frame, _, _ = svd(opt_grad - opt_grad.mean(axis=1, keepdims=True), full_matrices=False)
    rad_e, rad_v = radiation_rays(g)
    radiative = rad_v[:, np.argsort(rad_e)[-4:]]
    # Include only the strongest reactive mode and response gradients in the shared search.
    ray_basis, _ = np.linalg.qr(g.project_off(np.column_stack((ref, frame[:, :3], radiative)),g.v))
    q, t = approx_forms(g, ray_basis)
    return ray_basis, q, t


def rays(q, t, theta):
    tetrahedron = np.array([[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])
    rotation = np.array([[1, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1]]).T
    coeff = np.vstack((np.ones(4), tetrahedron * theta[:3,None],
                       np.ones(4) * theta[3], 1j * rotation * theta[4:,None]))
    coeff /= np.linalg.norm(coeff, axis=0)[None, :]
    values = [np.einsum('ij,bik,kj->bj',coeff.conj(), a,coeff).real for a in [q,t]]
    return coeff, *values


def mean_span(g, y):
    x0 = fit_rows(g, y)
    residual_gradient = g.project_off(g.green(x0, g.symbol.imag) - .5j * g.b, g.v)
    row_e = g.green(x0)
    block_gradient = []
    sector_gradient = []
    for ix in g.indices:
        part_x, part_e = np.zeros_like(x0), np.zeros_like(x0)
        part_x[ix], part_e[ix] = x0[ix], (row_e + g.b)[ix]
        response = g.green(part_x, g.symbol.conj())
        block_gradient.append((part_e - response) / (2j))
        sector_gradient.append(part_x - .5*.75*(part_e + response))
    _, reaction = reactive_rays(g)
    opt = g.project_off(np.column_stack([residual_gradient] + block_gradient + sector_gradient + [reaction]),g.v)
    basis, values, _ = svd(opt, full_matrices=False)
    basis = basis[:, values > 1e-10*values[0]]
    eb = g.green(basis)
    q0, t0 = approx_forms(g, x0, g.b)
    q, t = approx_forms(g, basis)
    lq = []
    lt = []
    for ix in g.indices:
        v, x, e1, e2 = basis[ix], x0[ix], eb[ix], (row_e+g.b)[ix]
        lq.append(v.conj().T @ x - .5*.75 * (v.conj().T @ e2 + e1.conj().T @ x))
        lt.append((v.conj().T @ e2 - e1.conj().T @ x) / (2j))
    lq,lt = np.array(lq),np.array(lt)
    return x0, basis, q0, t0, q, t, lq, lt


def mean(span, ridge, omega):
    x0, basis, q0, t0, q, t, lq, lt = span
    h0 = q.sum(axis=0) + omega*t.sum(axis=0)
    if ridge + np.linalg.eigvalsh(h0)[0] < 1e-5:
        return None
    z = -np.linalg.solve(h0 + ridge*np.eye(h0.shape[0]), lq.sum(axis=0)+omega*lt.sum(axis=0))
    forms = []
    for f0,a,l in [(q0,q,lq),(t0,t,lt)]:
        cross = z.conj().T @ l
        forms.append(f0 + z.conj().T @ a @ z + cross + cross.conj().swapaxes(-1,-2))
    moment = z.conj().T @ z + x0.conj().T @ x0
    return z, *forms, moment


def objective(params, span, rq, rt):
    coeff, qr, tr = rays(rq,rt,params[:7])
    built = mean(span, np.exp(params[7]), np.exp(params[8]))
    if built is None:
        return 1e8
    z, qx, tx, moment = built
    try:
        weights = np.linalg.solve(tr, -tx.reshape(4,-1)).reshape(tx.shape)
    except np.linalg.LinAlgError:
        return 1e8
    q = qx + np.einsum('ij,jkl->ikl',qr,weights)
    bad = max(0, 1e-3 + max_ev(-weights).max()) + 20*max(0, .01 + max_ev(q).max())
    energy = max_ev(moment + weights.sum(axis=0))
    return (max(0,energy) + 1e4*bad) / 1e3


def run():
    s,b,symbol,data,_ = load()
    g = NumericGeometry(s,b,symbol,blocks(2))
    basis, rq, rt = make_ray_span(g)
    proposals = {}
    for case,y in data.items():
        span = mean_span(g,y)
        print(case, 'X original norm', norm(span[0],2)**2, flush=True)
        bounds = [(-.35,.35)]*3+[(-.13,.13)]*4 + [(np.log(.06),np.log(70)),(np.log(.005),np.log(50))]
        result = differential_evolution(objective, bounds, args=(span,rq,rt),
                                        maxiter=220, popsize=15, rng=1821,tol=1e-6,
                                        callback=lambda x,cv: (print('candidate', objective(x,span,rq,rt), flush=True) if cv > .002 else None))
        params = result.x
        coeff,qr,tr = rays(rq,rt,params[:7])
        z,qx,tx,moment = mean(span, np.exp(params[7]),np.exp(params[8]))
        weights = np.linalg.solve(tr, -tx.reshape(4,-1)).reshape(tx.shape)
        q = qx + np.einsum('ij,jkl->ikl',qr,weights)
        print('num cost',g.green_columns,'status',result.message,'params',params,'score',result.fun)
        print('ray norms', np.linalg.norm(coeff,axis=0),'q',qr,'tau',tr,'radiation',tr.sum(axis=0))
        print('weights floor',-max_ev(-weights),'q max',max_ev(q),'energy', max_ev(moment+weights.sum(axis=0)),flush=True)
        proposals[case]=span[0]+span[1] @ z
        proposals[case+'__rays']=basis@coeff
    np.savez_compressed(HERE/'ray_search_proposals.npz',**proposals)


if __name__=='__main__':
    with threadpool_limits(1):
        run()
