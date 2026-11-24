"""Find a finite-energy uniform witness by optimizing the dark part of its mean."""

import argparse
import json
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import norm
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from search_rays import (approx_forms, blocks, load, NumericGeometry, make_ray_span,
                         local_optical_gradient, rays, mean_span, mean, HERE)
from probe import eigen_directions


def real_coordinates(z):
    return np.stack([z.real, z.imag], axis=0).ravel()


def spectral_penalty(a, floor):
    d, v = np.linalg.eigh(a)
    w = np.maximum(0, (d-floor))
    return np.sum(w*w)/2, (v*w[...,None,:]) @ v.conj().swapaxes(-1,-2)


class Objective:
    def __init__(self, span, ray_q, ray_t, anchor_q=None, anchor_t=None):
        self.span = span
        self.ray_q = ray_q
        self.t_inverse = np.linalg.inv(ray_t)
        self.gram0 = span[0].conj().T @ span[0]
        self.shape = span[-1].shape[1:]
        size = self.shape[1]
        self.q_shift = 0 if anchor_q is None else anchor_q[:,None,None]*np.eye(size)
        self.t_shift = 0 if anchor_t is None else anchor_t[:,None,None]*np.eye(size)
        self.current_anchor = 0

    def forms(self, z):
        x0, basis, q0, t0, q, t, lq, lt = self.span
        zstar = z.conj().T
        out = []
        for f0, a, lin in [(q0,q,lq),(t0,t,lt)]:
            cross = zstar @ lin
            out.append(f0 + zstar @ a @ z + cross + cross.conj().swapaxes(-1,-2))
        qx, tx = out
        qx += self.current_anchor*self.q_shift
        tx += self.current_anchor*self.t_shift
        weights = -np.einsum('ib,bjk->ijk',self.t_inverse,tx)
        sector = qx + np.einsum('bi,ijk->bjk',self.ray_q,weights)
        moment = self.gram0 + zstar @ z + weights.sum(axis=0) + self.current_anchor*np.eye(z.shape[1])
        return weights, sector, moment

    def evaluate(self, p, cap, penalty):
        zre, zim = p.reshape((2, *self.shape))
        z = zre + 1j*zim
        weights, sector, moment = self.forms(z)
        q, t, lq, lt = self.span[4:]
        fw, gw = spectral_penalty(-weights, -.25)
        fq, gq = spectral_penalty(sector, -.1)
        fe, ge = spectral_penalty(moment, cap)
        gw = -gw + 20*np.einsum('bi,bjk->ijk',self.ray_q,gq) + ge[None,:,:]
        gt = -np.einsum('ib,ijk->bjk',self.t_inverse,gw)
        df = ((q@z+lq) @ (20*gq)).sum(axis=0) + ((t@z+lt) @ gt).sum(axis=0) + z @ ge
        return (fw+20*fq+fe)*penalty, real_coordinates(2*df)*penalty


def maxev(a):
    return np.linalg.eigvalsh(a)[..., -1]


def load_rays(g, artifact, case, regularization=None, spread=.08):
    if artifact is not None:
        v = artifact[case+'__rays']
        hidden = g.project_off(v, g.v)
        hidden /= np.linalg.norm(hidden, axis=0)[None,:]
        out = [a.diagonal(axis1=1,axis2=2).real for a in approx_forms(g,hidden)]
        return hidden, *out
    if regularization is not None:
        symbol = .75*g.symbol.real - regularization*g.symbol.imag
        ev, v = eigen_directions(g, lambda x:g.green(x,symbol), 2, 'LA')
        ref = v[:,np.argmax(ev)]
        opt_grad = local_optical_gradient(g, ref)
        frames, _, _ = np.linalg.svd(opt_grad-opt_grad.mean(axis=1, keepdims=True),full_matrices=False)
        tetrahedron = np.array([[1, 1, -1, -1], [1,-1,1,-1],[1,-1,-1,1]])
        z = g.project_off(ref[:,None] + spread*frames[:,:3]@tetrahedron,g.v)
        z /= np.linalg.norm(z,axis=0)[None,:]
        out = [a.diagonal(axis1=1,axis2=2).real for a in approx_forms(g,z)]
        return z, *out
    basis, rq, rt = make_ray_span(g)
    theta = np.array([.08027864, -.08363065, .10121346, .01370608,
                      .01830999, -.00771984, -.00735677])
    coeff, qr, tr = rays(rq,rt,theta)
    return basis @ coeff, qr, tr


def make_anchor(g, radiative_penalty):
    if radiative_penalty is None:
        source = HERE.parents[1]/'submission/single/local_moment.npz'
        with np.load(source) as prior:
            stable = prior['rays'][:, :1]
    else:
        symbol = .75*g.symbol.real - radiative_penalty*g.symbol.imag
        v, q = eigen_directions(g, lambda x:g.green(x,symbol), 2, 'LA')
        mode = q[:,np.argmax(v)].reshape(96,96)
        stable = np.column_stack([np.rot90(mode,j).ravel() for j in range(4)])
    stable = g.project_off(stable, g.v)
    stable /= np.linalg.norm(stable, axis=0)[None,:]
    sq, st = [a.diagonal(axis1=1,axis2=2).real.mean(axis=1) for a in approx_forms(g,stable)]
    return stable, sq, st


def run(args):
    wall,cpu = perf_counter(),process_time()
    s,b,symbol,data,_ = load()
    g = NumericGeometry(s,b,symbol,blocks(2))
    artifact = np.load(args.rays) if args.rays else None
    shared_start,shared_cpu = perf_counter(),process_time()
    shared_green0 = g.green_columns
    anchors = make_anchor(g, args.anchor_penalty) if args.anchors else None
    shared = {'numeric_fft_columns':g.green_columns-shared_green0,
              'seconds_wall':perf_counter()-shared_start, 'seconds_cpu':process_time()-shared_cpu}
    proposals, costs = {},{}
    for case,observed in data.items():
        y = args.data_scale*observed
        cw, cc, gc = perf_counter(),process_time(),g.green_columns
        dark,qr,tr = load_rays(g, artifact, case,args.ray_penalty,args.ray_spread)
        span = mean_span(g,y)
        if args.anchors:
            stable,sq,st = anchors
            o = Objective(span,qr,tr,sq,st)
        else:
            stable = None
            o = Objective(span,qr,tr)
        z,*_ = mean(span,2.0,2.4)
        result = None
        accepted = None
        attempts = []
        plan = [(cap,0.) for cap in args.caps] if not args.anchors else [
            (anchor+cap,anchor) for anchor in args.anchors for cap in args.caps]
        accepted_score = np.inf
        for cap,anchor in plan:
            o.current_anchor = anchor
            result = minimize(o.evaluate, real_coordinates(z), jac=True, method='L-BFGS-B',
                              args=(cap, 1.e-3), options={'maxiter': args.maxiter,
                                                         'gtol': 1e-8, 'ftol': 1e-12})
            zre, zim = result.x.reshape((2, *o.shape))
            z = zre + 1j*zim
            w,q,m = o.forms(z)
            energies, neg_weights, sectors = maxev(m), maxev(-w), maxev(q)
            check = {'fit_target': cap, 'anchor_energy':anchor, 'form_energy': float(energies),
                     'q_upper': float(np.max(sectors)), 'weight_floor': -float(np.max(neg_weights)),
                     'status': str(result.message), 'iterations': int(result.nit),
                     'evaluations': int(result.nfev)}
            print(case, check, flush=True)
            attempts.append(check)
            if sectors.max() < -0.08 and neg_weights.max() < -.08 and energies < accepted_score:
                accepted = z.copy()
                accepted_score = energies
                accepted_anchor = anchor
        if accepted is None:
            raise ArithmeticError(f'No feasible proposal for {case}.')
        proposals[case]=span[0]+span[1] @ accepted
        proposals[case+'__rays']=dark
        proposals[case+'__fit_scale']=args.data_scale
        if stable is not None:
            proposals[case+'__auxiliary_ray'] = stable
            proposals[case+'__auxiliary_weight'] = accepted_anchor
        costs[case] = {'attempts': attempts, 'numeric_fft_columns': g.green_columns-gc,
                       'seconds_wall': perf_counter()-cw, 'seconds_cpu': process_time()-cc}
    np.savez_compressed(args.output.with_suffix('.npz'), **proposals)
    args.output.write_text(json.dumps({'cases': costs, 'fit_scale':args.data_scale,
                                       'anchor_penalty':args.anchor_penalty,'shared_preparation':shared,
                                       'seconds_wall': perf_counter()-wall,
                                       'seconds_cpu': process_time()-cpu}, indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=type(HERE), default=HERE/'optimized_mean.json')
    parser.add_argument('--rays', type=type(HERE))
    parser.add_argument('--ray-penalty', type=float)
    parser.add_argument('--ray-spread', type=float,default=.08)
    parser.add_argument('--anchors', type=float, nargs='+')
    parser.add_argument('--anchor-penalty', type=float)
    parser.add_argument('--data-scale', type=float, default=1.)
    parser.add_argument('--maxiter', type=int, default=900)
    parser.add_argument('--caps', type=float, nargs='+', default=[50000, 22000, 10000, 5000])
    with threadpool_limits(1):
        run(parser.parse_args())
