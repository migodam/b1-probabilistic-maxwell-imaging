"""Numerical global and matched local retention-gradient proposals from public Y."""

import argparse
from collections import Counter
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, eigvalsh, qr
from scipy.linalg.lapack import get_lapack_funcs
from threadpoolctl import threadpool_limits

from curve_inputs import CAPS, HERE, full_ring_inputs

# Both scalar and block controls use these same multivariate step sizes.
STEPS = (.02, .06, .15, .30)
# Include the previous fine ratio grid and a broad range in both directions.
RATIOS = (-10000., -100., -10., -3., -1., -.5, -.25, -.1, 0., .05, .1, .2,
          .3, .4, .55, .75, 1., 1.5, 2., 3., 5., 10., 30., 100., 10000.)
MU = np.unique(np.r_[np.geomspace(1e-6, .08, 27),
                     .00001, .0001, .00035, .001, .0025, .005])


def h(a):
    return (a + a.conj().T) / 2


def times(wall, cpu):
    return dict(wall_seconds=perf_counter()-wall, cpu_seconds=process_time()-cpu)


def save_json(path, obj):
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


class Projected:
    """Use a small numerical search space. Proof never assumes this space is complete."""

    def __init__(self, spec, s, b, q):
        n = int(np.sqrt(b.shape[0]))
        p = len(spec)
        indices = np.arange(n)
        modes = np.arange(-q, q+1)
        wave = np.exp(2j*np.pi*indices[:, None]*modes/p)
        wave = (wave[:, None, :, None]*wave[None, :, None, :]).reshape(n*n, -1)
        wave, _ = qr(wave, mode='economic', check_finite=False)

        def green(x):
            x = x.reshape(n, n, -1)
            fourier = fft2(x, s=spec.shape, axes=(0, 1))
            return ifft2(spec[:, :, None]*fourier, axes=(0, 1))[:n, :n].reshape(n*n, -1)

        ports = np.column_stack((s.conj().T, b))
        # Add aperture columns and their Green images beyond the low Fourier square.
        more = np.column_stack((ports, green(ports)))
        more -= wave @ (wave.conj().T@more)
        extra, sv, _ = np.linalg.svd(more, full_matrices=False)
        v = np.column_stack((wave, extra[:, sv > 1.e-9*sv[0]]))
        v, _ = qr(v, mode='economic', check_finite=False)
        gv = green(v)
        self.basis = v
        self.S = s @ v
        self.S2 = h(self.S.conj().T @ self.S)
        self.dimension = v.shape[1]
        self.full_dimension = n*n
        groups = ((2*indices[:, None]//n)*2 + 2*indices[None, :]//n).ravel()
        self.blocks = []
        for j in range(4):
            active = groups == j
            proj = v[active].conj().T
            pg = proj @ gv[active]
            self.blocks.append((h(proj@v[active]), h(pg), h(pg/1j), proj@b[active]))
        self.Gh = sum(row[1] for row in self.blocks)
        self.Ga = sum(row[2] for row in self.blocks)
        self.B = sum(row[3] for row in self.blocks)
        self.s = b.shape[1]
        self.I = np.eye(self.dimension)
        self.cost = Counter(projected_green_columns=v.shape[1]+ports.shape[1])

    def scalar(self, u, t, data, mu_grid=MU):
        """Shared A factor and sensor update across all mu and observations."""
        d, v = eigh(h(self.I-u*self.Gh+t*self.Ga), check_finite=False, driver='evd')
        self.cost['projected_scalar_state_eigendecompositions'] += 1
        if np.min(np.abs(d)) <= 1.e-12:
            return {}
        sv = (self.S@v).conj().T
        zv = v.conj().T@self.B * (.5*(u+1j*t))
        yy = np.concatenate([v.y for v in data], axis=1)
        dy = sv@yy
        yg = [h(r.y.conj().T@r.y) for r in data]
        svd = sv / d[:, None]
        dr, vr = eigh(h(sv.conj().T@svd), check_finite=False)
        self.cost['scalar_receiver_eigendecompositions'] += 1
        vp = vr.conj().T @ sv.conj().T
        drd, drz = dy/d[:, None], zv/d[:, None]
        vpd, vpz = vp@drd, vp@drz
        out = {}
        for mu in mu_grid:
            ev = 1 + dr/mu
            if np.min(np.abs(ev)) <= 1.e-10:
                continue
            nu = int(np.count_nonzero(d < 0) - np.count_nonzero(ev < 0))
            if nu < 0:
                continue
            self.cost['projected_scalar_profiles'] += 1
            if nu >= self.s:
                out[float(mu)] = (nu, np.zeros(len(data)), np.zeros(len(data)))
                continue
            h0 = (dy+mu*np.tile(zv, (1, len(data)))).reshape(self.dimension, len(data), self.s)
            xx = ((drd+mu*np.tile(drz, (1, len(data)))) -
                  (svd@vr)@((vpd+mu*np.tile(vpz, (1, len(data))))/ev[:, None])/mu) / mu
            xx = xx.reshape(self.dimension, len(data), self.s)
            hh = np.array(yg)-np.einsum('kji,kjl->jil', h0.conj(), xx)
            vals = np.linalg.eigvalsh((hh+hh.conj().transpose(0, 2, 1))/2)[:, ::-1]
            self.cost['projected_source_eigendecompositions'] += len(data)
            out[float(mu)] = (nu, np.sqrt(np.maximum(vals[:, nu:], 0.).sum(axis=1)), vals[:,nu])
        return out

    def local(self, u, w, lam, y, *, need_gradient=False, with_sources=False):
        a = self.S2.copy()
        rhs = self.S.conj().T@y
        for weight, optical, (p, re, im, b) in zip(w, lam, self.blocks, strict=True):
            a += weight*(p-u*re)+optical*im
            rhs += .5*(weight*u+1j*optical)*b
        factor, lwork, solve = get_lapack_funcs(('hetrf', 'hetrf_lwork', 'hetrs'), (a,))
        work, _ = lwork(len(a), lower=True)
        lu, pivots, info = factor(a, lower=True, lwork=int(work.real))
        self.cost['projected_local_state_ldl'] += 1
        if info:
            return dict(nu=None, d=0., status='LDL unresolved')
        diagonal = lu.diagonal().real.copy()
        i = 0
        while i < len(a):
            if pivots[i] < 0:
                middle = (diagonal[i]+diagonal[i+1])/2
                width = np.hypot((diagonal[i]-diagonal[i+1])/2, abs(lu[i+1, i]))
                diagonal[i], diagonal[i+1] = middle-width, middle+width
                i += 1
            i += 1
        nu = int(sum(diagonal < 0))
        if nu >= self.s:
            return dict(nu=nu, d=0., status='illumination budget exhausted')
        x, info = solve(lu, pivots, rhs, lower=True)
        if info or np.linalg.norm(a@x-rhs) > 1.e-8*np.linalg.norm(rhs):
            return dict(nu=nu, d=0., status='stationary solve unresolved')
        vals, v = eigh(h(y.conj().T@y-rhs.conj().T@x), check_finite=False)
        vals, v = vals[::-1], v[:, ::-1]
        self.cost['projected_source_eigendecompositions'] += 1
        d = float(np.sqrt(np.maximum(vals[nu:], 0.).sum()))
        row = dict(nu=nu, d=d, first_retained_eigenvalue=float(vals[nu]), status='evaluated')
        if with_sources:
            row['top_source_modes'] = v
        if need_gradient:
            active = (np.arange(self.s) >= nu) & (vals > 0.)
            positive_count = int(active.sum())
            if not positive_count:
                active[nu] = True
            v = v[:, active]
            xv = x @ v
            dm, dl = [], []
            for p, re, im, b in self.blocks:
                vxb = np.vdot(xv, b@v)
                dm.append((np.vdot(xv, (p-u*re)@xv).real-u*vxb.real))
                dl.append((np.vdot(xv, im@xv).real+vxb.imag))
            row.update(gradient_mu=dm, gradient_lambda=dl, retained_positive_count=positive_count,
                       gradient_objective='positive_tail_sum' if positive_count else 'first_signed_retained_eigenvalue')
        return row


def profile(mu, t, d, nu, front):
    return dict(mu=mu, t=t, d=float(d), nu=nu, first_retained_eigenvalue=float(front))


def best(proposals, *, nz=False):
    can = [r for r in proposals if r['nu'] is not None and (not nz or r['nu']>0)]
    return max(can, key=lambda r: (r['d'],r['first_retained_eigenvalue']), default=None)


def local_arms(p, u, base, obs):
    """Fixed inertia across all reported gradient improvements."""
    if base is None:
        return dict(status='no admissible global seed', trials=[])
    mu, t = base['mu'], base['t']
    q = p.local(u, np.repeat(mu, 4), np.repeat(mu*t, 4), obs.y, need_gradient=True)
    if 'gradient_mu' not in q:
        return dict(status='unable to evaluate a retained source mode', baseline=q, trials=[])
    gm = np.array(q['gradient_mu'])
    gl = np.array(q['gradient_lambda'])
    trials = []
    for step in STEPS:
        # Use the same absolute multiplier displacement limit in both controls.
        for arm in ('real_scalar', 'real_blocks', 'joint_scalar', 'joint_blocks', 'optical_scalar', 'optical_blocks'):
            block = arm.endswith('blocks')
            dm = np.sign(gm if block else np.full(4, gm.sum()))
            dl = np.sign(gl if block else np.full(4, gl.sum()))
            w = mu*(1+step*dm) if not arm.startswith('optical') else np.repeat(mu, 4)
            el = mu*(t+step*max(1, abs(t))*dl) if not arm.startswith('real') else np.repeat(mu*t, 4)
            r = p.local(u, w, el, obs.y)
            trials.append(dict(arm=arm, step=step, weights=w.tolist(), optical=el.tolist(),
                               same_inertia=(r['nu']==q['nu']), **r))
    return dict(status='evaluated', baseline=dict(mu=mu, t=t, **q), trials=trials)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=7)
    parser.add_argument('--cap', type=float, action='append')
    parser.add_argument('--out', type=Path, default=HERE/'projected_search.json')
    args = parser.parse_args()
    start, cpu = perf_counter(), process_time()
    spec, s, b, data = full_ring_inputs()
    p = Projected(spec, s, b, args.q)
    rows = []
    setup = times(start, cpu)
    by_cap = {}
    for u in args.cap or CAPS:
        twall, tcpu = perf_counter(), process_time()
        before = p.cost.copy()
        proposals = [[] for i in data]
        for t in RATIOS:
            for mu, (nu, d, first) in p.scalar(u, t, data).items():
                for i, d0 in enumerate(d):
                    proposals[i].append(profile(mu, t, d0, nu, first[i]))
        for obs, options in zip(data, proposals, strict=True):
            admissible = [p0 for p0 in options if p0['nu'] < b.shape[1]]
            result = {arm: best([r for r in admissible if accept(r)]) for arm, accept in (
                ('global_real', lambda p0: p0['t']==0.),
                ('global_optical', lambda p0: p0['t']!=0.),
                ('global_indefinite', lambda p0: p0['nu']>0))}
            family = dict()
            for method, initial in (
                ('from_real', result['global_real']),
                ('from_best', best([r for r in result.values() if r is not None])),
            ):
                family[method] = local_arms(p, u, initial, obs)
            if u >= .75 and result['global_indefinite'] != result['global_optical']:
                family['from_indefinite'] = local_arms(p, u, result['global_indefinite'], obs)
            rank = sorted(admissible, key=lambda r:r['d'], reverse=True)[:12]
            nus = Counter(z['nu'] for z in options)
            best_surviving = min(options, key=lambda r: (r['nu'], abs(r['t']), abs(r['mu']-.00035)))
            rows.append(dict(id=obs.id, cap=u, radius=obs.padded_radius, object_id=obs.object_id,
                             noise_pct=obs.noise_pct, arms=result, local=family,
                             best_global_profiles=rank, minimum_inertia_profile=best_surviving,
                             surveyed_inertia_counts=dict(nus), n_profiles=len(options)))
        by_cap[str(u)] = dict(incremental_operations=dict(p.cost-before), **times(twall, tcpu))
        print(u, [(r['noise_pct'],r['object_id'],round(max(
            (p['d'] for p in r['arms'].values() if p is not None), default=0.)/r['radius'],3))
                  for r in rows[-len(data):]], by_cap[str(u)], flush=True)
        save_json(args.out, dict(scope='Numerical proposals from projected native FFT; full-space verification required',
                                 q=args.q, search_dimension=p.dimension, native_dimension=p.full_dimension,
                                 object_count=6, public_center_count=len(data), ratios=RATIOS,
                                 mu_grid=MU.tolist(), local_step_factors=STEPS,
                                 shared_setup=setup, by_cap=by_cap, rows=rows,
                                 **times(start,cpu)))


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
