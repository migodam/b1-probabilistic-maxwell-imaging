#!/usr/bin/env python3
"""Run the two 2026-11-20 Maxwell comparisons. See VERDICT_REPORT_ZH.md.

`all` generates observations, writes every prediction, then scores them.
`prepare`, `solve` and `score` allow staging on separate hosts. Solvers open only
training observations. Array-valued predictions and results use plain NPZ.
"""
from __future__ import annotations

import os
for flag in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[flag] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import csv
import importlib.util
import json
from pathlib import Path
import platform
import signal
import sys
import time

sys.dont_write_bytecode = True
import numpy as np
from scipy import linalg as la
from scipy.fft import ifftn
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, gmres

OUT = Path(__file__).resolve().parent
WORKSPACE = Path(os.environ.get('MAXWELL_KILLTEST_WORKSPACE', str(OUT.parents[1])))
WORKTREES = WORKSPACE / 'research/delegated/max_struct_r1'
R8_REPO = WORKTREES / 'r8_worktree'
if not (R8_REPO / 'research/research_max_struct_r8/src/paper_physics.py').is_file():
    R8_REPO = WORKSPACE
sys.path.insert(0, str(R8_REPO))
R14 = WORKTREES / 'r14_worktree/research/research_max_struct_r14'
if not (R14 / 'src/dualmetric_core.py').is_file():
    R14 = WORKSPACE / 'research/research_max_struct_r14'
A24 = Path(os.environ.get('MAXWELL_KILLTEST_A24', str(WORKSPACE / 'Gaussian/A26/a26-om-feedback-prior-completion')))
A17 = WORKSPACE / 'Gaussian/A17'
FRESNEL = WORKSPACE / 'research/delegated/a2_literature/data/2005_iop_21_6_S09'
RUNTIME = OUT / '.runs'
C0 = 299792458.0
A_ARMS = ('paper_nie', 'cie_dual', 'block_cie_01_sw', 'block_cie_11_sw', 'dr_standalone')
B_ARMS = ('receiver_svd', 'projected_internal', 'a17_proj_fast', 'full_gradient_matched_time')
A_BUDGETS = (15., 30.)
LIGHT_SEEDS = (12., 24.)
A_AUDIT_PERIOD = 1.
B_BUDGET = 15.
A17_PRIOR = 1e-5
A17_SEED = 171030

# Fixed before model fitting / scoring. The new geometries use seeded rotations
# and translations and are distinct from the R8/R13/R14 cohorts.
A_CASES = (
    ('austria_eps4_400MHz', '1A', 4., 400., 150, 'full'),
    ('austria_eps5p5_400MHz', '1A', 5.5, 400., 150, 'full'),
    ('austria_eps7_400MHz', '1A', 7., 400., 150, 'full'),
    ('austria_eps3p5_600MHz', '1A', 3.5, 600., 150, 'full'),
    ('austria_eps3p5_400MHz_limited180', '1B', 3.5, 400., 150, 'half'),
    ('austria_eps4_400MHz_sparse8tx', '1B', 4., 400., 150, 'sparse'),
    ('new_double_rings', '1B', 5.5, 400., 96, 'full'),
    ('new_eccentric_shell', '1B', 6.5, 400., 96, 'full'),
    ('new_coupled_triplet', '1B', 8., 400., 96, 'full'),
    ('new_kinked_slot', '1B', 4.5, 400., 96, 'full'),
)
F_CASES = (('FoamTwinDielTM', 8, 13), ('FoamTwinDielTM', 8, 6),
           ('FoamDielIntTM', 6, 6), ('FoamDielIntTM', 8, 6),
           ('FoamDielExtTM', 6, 6), ('FoamDielExtTM', 8, 6))
B_NAMES = ('austria_like', 'two_holes', 'eccentric_shell', 'digit5', 'digit8',
           'gapped_islands', 'thin_bridge', 'two_phase_cavity')


def module_at(name, path):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def specs():
    cases = [dict(case=name, suite=suite, eps=eps, frequency=mhz*1e6,
                  fine=fine, side=32, width=2., aperture=aperture, index=i)
             for i, (name, suite, eps, mhz, fine, aperture) in enumerate(A_CASES)]
    cases += [dict(case=f'{name}_{freq}GHz_{ns}tx', suite='1C', object=name,
                   frequency=freq*1e9, ns=ns, side=60, width=.2, fine=0)
              for name, freq, ns in F_CASES]
    cases += [dict(case='a17_'+name, suite='2A', index=i, side=24, fine=48,
                   width=1., frequency=2*C0) for i, name in enumerate(B_NAMES)]
    return cases


def grid(side, width, xfirst=False):
    h = width/side
    a = (np.arange(side)+.5)*h-width/2
    yy, xx = np.meshgrid(a, a, indexing='ij')
    return np.c_[xx.ravel(), yy.ravel()] if not xfirst else np.c_[yy.ravel(), xx.ravel()]


def circle(x, y, cx, cy, radius):
    return (x-cx)**2+(y-cy)**2 <= radius**2


def stroke(x, y, vertices, radius):
    p = np.c_[x, y]
    distance = np.full(len(p), np.inf)
    for a, b in zip(vertices[:-1], vertices[1:]):
        a, b = np.asarray(a), np.asarray(b)
        v = b-a
        h = np.clip((p-a)@v/(v@v), 0, 1)
        distance = np.minimum(distance, np.sum((p-a-h[:, None]*v)**2, axis=1))
    return distance <= radius**2


def material(case, coords):
    """Called only by observation generation and after prediction sealing."""
    from research.research_max_struct_r8.src.paper_physics import austria
    i = case['index']
    if case['suite'] in ('1A', '1B') and i < 6:
        return austria(coords, case['eps'])
    campaign_b = case['suite'] == '2A'
    rng = np.random.default_rng((20261120100 if campaign_b else 20261120000)+i)
    offset = rng.uniform(-.045, .045, 2)
    angle = rng.uniform(-.22, .22)
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    x, y = ((coords-offset) @ rot).T
    if not campaign_b:
        if i == 6:
            mask = ((circle(x,y,-.39,.05,.34) & ~circle(x,y,-.39,.05,.20)) |
                    (circle(x,y,.35,-.17,.30) & ~circle(x,y,.35,-.17,.17)))
            contrast = (case['eps']-1)*np.where(x<0, 1., .88)
        elif i == 7:
            mask = circle(x,y,.02,-.06,.54) & ~circle(x,y,.16,.03,.30)
            core = circle(x,y,.12,.01,.13)
            return (case['eps']-1)*(mask+.72*core)
        elif i == 8:
            mask = (circle(x,y,-.425,-.12,.23) | circle(x,y,.065,-.12,.23) |
                    circle(x,y,-.18,.315,.23))
            contrast = (case['eps']-1)*np.where(y>0, 1., .80)
        else:
            mask = (np.abs(x)<.65) & (np.abs(y)<.43)
            mask &= ~stroke(x, y, np.array([[-.70,-.22],[-.22,-.22],[-.09,.12],[.45,.20],[.70,.15]]), .065)
            contrast = (case['eps']-1)*np.where(y>0, 1., .83)
        return contrast*mask
    if i == 0:
        mask = ((circle(x,y,0,-.10,.29) & ~circle(x,y,0,-.10,.145)) |
                circle(x,y,.155,.255,.087) | circle(x,y,-.155,.255,.087))
    elif i == 1:
        mask = (np.abs(x)<.30) & (np.abs(y)<.225)
        mask &= ~(circle(x,y,-.142,0.,.092) | circle(x,y,.142,.051,.101))
    elif i == 2:
        mask = circle(x,y,0,0,.31) & ~circle(x,y,.11,.030,.207)
    elif i == 3:
        mask = stroke(x,y,np.array([[.22,.32],[-.20,.32],[-.23,.04],[.105,.04],
                                   [.23,-.07],[.21,-.25],[.06,-.32],[-.22,-.30]]),.063)
    elif i == 4:
        top = (x/.205)**2+((y-.161)/.194)**2
        bot = (x/.236)**2+((y+.162)/.205)**2
        mask = ((top<=1)&(top>=.36)) | ((bot<=1)&(bot>=.36))
    elif i == 5:
        mask = circle(x,y,-.21,-.13,.132) | circle(x,y,.08,-.13,.135) | circle(x,y,-.07,.111,.130)
    elif i == 6:
        mask = (circle(x,y,-.24,-.14,.138) | circle(x,y,.23,.11,.134) |
                stroke(x,y,np.array([[-.24,-.14],[.01,-.055],[.23,.11]]),.041))
    else:
        mask = (x/.33)**2+(y/.24)**2 < 1
        mask &= ~circle(x,y,.10,-.016,.117)
    c = (3.55 + .35*np.sin(6*x+1.2*i)*np.cos(4*y))
    loss = .14 + .22*(.5+.5*np.cos(5*x-7*y+i))
    if i in (1,5,7):
        c = np.where(x+y<.02, c, .52*c)
        loss = np.where(x+y<.02, loss, .6*loss)
    return mask*(c+1j*loss)


def fine_solve(chi, E, S, apply, diagonal, *, rtol=2e-8):
    """R8 FFT-embedded VIE, restricted to exactly nonzero generating cells."""
    active = np.flatnonzero(chi != 0)
    local = chi[active]
    n = len(chi)
    def mv(z):
        extended = np.zeros(n, complex)
        extended[active] = z
        return z-local*apply(extended[:, None])[:, 0][active]
    operator = LinearOperator((len(active), len(active)), matvec=mv, dtype=complex)
    jacobi = LinearOperator(operator.shape, matvec=lambda z: z/(1-local*diagonal), dtype=complex)
    J = np.empty((len(active), E.shape[1]), complex)
    max_residual, counts = 0., []
    for s, inc in enumerate(E[active].T):
        rhs = local*inc
        history = []
        result, info = gmres(operator,rhs,M=jacobi,rtol=rtol,atol=0.,restart=120,
                             maxiter=12,callback=history.append,callback_type='pr_norm')
        residual = np.linalg.norm(mv(result)-rhs)/np.linalg.norm(rhs)
        if info or residual > 3*rtol:
            raise RuntimeError(f'fine forward source {s}: GMRES {info}, residual {residual:.3g}')
        J[:, s] = result
        max_residual = max(max_residual, float(residual))
        counts.append(len(history))
    return S[:, active]@J, dict(fine_active_cells=len(active), forward_rhs=E.shape[1],
                                fine_gmres_iterations=sum(counts), fine_state_residual=max_residual)


def complex_noise(clean, fraction, seed):
    rng = np.random.default_rng(seed)
    direction = rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape)
    return direction*(fraction*np.linalg.norm(clean)/np.linalg.norm(direction))


def prepare_one(case):
    from research.research_max_struct_r8.src.paper_physics import acquisition, EqualAreaGreen
    from research.research_max_struct_r8.src.fresnel_adapter import load_file, operators, incident_fit
    path = public_path(case)
    if path.exists():
        return case['case']
    start = time.perf_counter()
    name, suite, n = case['case'], case['suite'], case['side']
    fine = case['fine']
    case = dict(case)
    if suite == '1C':
        directory = Path(os.environ.get('FRESNEL_2005_DIR', str(FRESNEL)))
        _, Y, I = load_file(directory/(case['object']+'.exp'), case['frequency']/1e9)
        ids = np.flatnonzero(np.any(np.isfinite(Y),axis=0))
        Y, I = Y[:, ids], I[:, ids]
        ns = len(ids)
        tx_angle = ids*2*np.pi/ns
        rx_angle = np.arange(Y.shape[0])*2*np.pi/Y.shape[0]
        sources = 1.67*np.c_[np.cos(tx_angle),np.sin(tx_angle)]
        receivers = 1.67*np.c_[np.cos(rx_angle),np.sin(rx_angle)]
        centers = grid(n, case['width'])
        k = 2*np.pi*case['frequency']/C0
        p = operators(centers,sources,receivers,case['width']/n,k,dense=False)
        E, calibration, _ = incident_fit(I,sources,receivers,centers,k)
        if ns == 18 and case['ns'] == 13:
            held = np.array([0,4,8,12,16]); train = np.setdiff1d(np.arange(ns),held)
        elif ns == 18:
            train = np.array([0,3,6,9,12,15]); held = np.setdiff1d(np.arange(ns),train)
        else:
            held = np.array([0,4]); train = np.setdiff1d(np.arange(ns),held)
        rx_ids = np.stack([np.flatnonzero(np.isfinite(Y[:,s])) for s in train])
        Y_train = np.stack([Y[rx_ids[j],s] for j,s in enumerate(train)],axis=1)
        case.update(train_sources=(train+1).tolist(), held_sources=(held+1).tolist(),
                    train_channels=int(Y_train.size), held_channels=int(np.isfinite(Y[:,held]).sum()),
                    calibration_median_validation=float(np.median([r['validation_relative'] for r in calibration])),
                    calibration_max_validation=max(r['validation_relative'] for r in calibration),
                    data_ball_ratio=.15, noise_relative=None, source_specific=True)
        np.savez(RUNTIME/'heldout'/f'{name}.npz', E=E[:,held], S=p['S'], Y=Y[:,held])
        arrays = dict(E=E[:,train],S=p['S'],rx_ids=rx_ids,Y=Y_train,centers=centers,
                      k=k,radius=p['radius'],delta=.15*np.linalg.norm(Y_train))
    elif suite != '2A':
        f = acquisition(fine,case['frequency'],dense=False)
        g = acquisition(n,case['frequency'],dense=False)
        tx = np.arange(20); rx = np.arange(40)
        if case['aperture'] == 'half':
            tx, rx = tx[:10], rx[:20]
        elif case['aperture'] == 'sparse':
            tx = np.floor(np.arange(8)*20/8).astype(int)
            rx = np.floor(np.arange(24)*40/24).astype(int)
        E_high, S_high = f['E0'][:,tx], f['S'][rx]
        G_high = EqualAreaGreen(f['centers'],(fine,fine),f['radius'],f['k'])
        clean, meta = fine_solve(material(case,f['centers']),E_high,S_high,G_high.apply,G_high.diagonal)
        noise = complex_noise(clean,.10,202611201000+case['index'])
        case.update(meta,noise_relative=.10,noise_seed=202611201000+case['index'],
                    train_channels=int(clean.size),tx=tx.tolist(),rx=rx.tolist(),
                    forward_FFT=G_high.cost['Green_FFT_2D_transforms']+G_high.cost['Green_IFFT_2D_transforms']+1)
        arrays = dict(E=g['E0'][:,tx],S=g['S'][rx],Y=clean+noise,centers=g['centers'],
                      k=g['k'],radius=g['radius'],delta=np.linalg.norm(noise))
    else:
        a2 = module_at('killtest_a2', A24/'vendor/a2/physics.py')
        def vie(side, aperture, nr, integrated):
            geometry = a2.Geometry(n=side,side=1.,n_tx=8,n_rx=nr,radius_tx=.75,
                                   radius_rx=.8,aperture=aperture)
            return a2.VIE(geometry,frequency_hz=case['frequency'],cell_integrated=integrated)
        # Fine square-cell integration, coarse point-collocation/nonzero disk self term.
        # One set of fine currents generates all receiver channels.
        f = vie(fine,'full',32,True)
        ang = np.linspace(-np.pi/3,np.pi/3,10)
        limited_rx = .8*np.c_[np.cos(ang),np.sin(ang)]
        s2_high = f.k**2*f._cell_integral(limited_rx,f.points)
        chi = material(case,f.points)
        clean, meta = fine_solve(chi,f.E,np.vstack((f.S,s2_high)),f.D.matmat,f.D.kernel[fine-1,fine-1])
        g = vie(n,'full',32,False)
        s2 = g.k**2*g.h**2*a2.point_green(limited_rx,g.points,g.k)
        case.update(meta,noise_relative=.02,tx=list(range(8)),limited_arc_deg=120.,
                    noise_seeds=[202611202000+case['index'],202611202100+case['index']],
                    forward_integration='fine GL4x4 square-cell with square self; coarse center and disk self')
        arrays = dict(E=g.E,S=g.S,S_limited=s2,G=g.D.block(np.arange(n*n),np.arange(n*n)),
                      Y=clean[:32]+complex_noise(clean[:32],.02,case['noise_seeds'][0]),
                      Y_limited=clean[32:]+complex_noise(clean[32:],.02,case['noise_seeds'][1]),
                      centers=g.points,k=g.k,volume=g.h**2)
    case['data_wall'] = time.perf_counter()-start
    arrays['case_json'] = np.array(json.dumps(case,ensure_ascii=False))
    np.savez(path,**arrays)
    print(f'prepared {name} in {case["data_wall"]:.2f}s',flush=True)
    return name


def public_path(case):
    return RUNTIME/'public'/f'{case["case"]}.npz'


def run_jobs(fn, jobs, workers):
    if workers == 1:
        return [fn(*job) for job in jobs]
    import multiprocessing as mp
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
        pending = [pool.submit(fn,*job) for job in jobs]
        return [p.result() for p in as_completed(pending)]


class TimeLimit(BaseException):
    """A timed operation may overrun; its completed result is then ineligible."""


class Clock:
    def __init__(self, limit):
        self.start = time.perf_counter()
        self.cpu_start = time.process_time()
        self.limit = limit
        self.counts = Counter()
        self.cuts = {}

    def watch(self,*times):
        self.cuts.update({t: None for t in times})

    def capture(self):
        now = self.wall()
        for t, row in self.cuts.items():
            if row is None and now>=t:
                self.cuts[t] = (now,self.cpu(),dict(self.counts))
        return now

    def wall(self):
        return time.perf_counter()-self.start

    def cpu(self):
        return time.process_time()-self.cpu_start

    def check(self):
        if self.capture() >= self.limit:
            raise TimeLimit()

    @contextmanager
    def span(self, name, **counts):
        self.check()
        self.counts.update(counts)
        yield


def array_file(path):
    with np.load(path,allow_pickle=False) as z:
        return {k:z[k] for k in z.files}


def stage_table(case):
    if case['suite'] == '1C':
        modes,beta,rank = (6,10,16),(6.,2.,.5),5
    elif case['eps'] < 4:
        modes,beta,rank = (5,7,10),(6.,2.,.5),15
    else:
        modes,beta,rank = (5,8,10),(8.,2.,.5),15
    return modes,beta,rank


class Ports:
    """Common receiver operator with the gather/scatter adjoint of actual masks."""
    def __init__(self, S, mask):
        self.S = S
        self.SH = S.conj().T.copy()
        self.mask = mask

    def receive(self, J):
        y = self.S@J
        return y if self.mask is None else np.take_along_axis(y,self.mask.T,axis=0)

    def adjoint(self, y):
        if self.mask is None:
            return self.SH@y
        r = np.zeros((len(self.S),y.shape[1]),complex)
        np.put_along_axis(r,self.mask.T,y,axis=0)
        return self.SH@r

    def indices(self, ns):
        if self.mask is None:
            return [np.arange(len(self.S))]*ns
        return list(self.mask)


def dense_from_fft(green):
    """Materialize the already defined finite Toeplitz kernel by lag lookup."""
    n0,n1 = green.shape
    axis = np.arange(n0*n1)
    a,b = np.divmod(axis,n1)
    kernel = ifftn(green.spectrum)
    green.cost['Green_IFFT_2D_transforms'] += 1
    return kernel[(a[:,None]-a)%green.pad[0],(b[:,None]-b)%green.pad[1]]


class PaperModel:
    """Full R8 VIE state and analytic adjoint, with R8 source-specific port masks."""
    def __init__(self, public, case, clock):
        from research.research_max_struct_r8.src.paper_physics import EqualAreaGreen
        from research.research_z_alg_r4_diversity.src.physics import Regularizer
        self.clock = clock
        self.shape = (case['side'],)*2
        self.E = public['E']
        self.n,self.ns = self.E.shape
        self.ports = Ports(public['S'],public.get('rx_ids'))
        self.green = EqualAreaGreen(public['centers'],self.shape,float(public['radius']),
                                    float(public['k']),clock.counts)
        self.G = dense_from_fft(self.green)
        self.out = None
        self.reg = Regularizer(self.shape,smooth=1e-4,tv=1e-4,epsilon=.01)
        clock.check()

    def state(self,m):
        self.clock.check()
        if self.out is not None and np.array_equal(self.out['m'],m):
            return self.out
        if not np.any(m):
            j = np.zeros_like(self.E)
            factor = None
        else:
            matrix = np.asfortranarray(-m[:,None]*self.G)
            matrix.flat[::self.n+1] += 1.
            self.clock.counts['LU'] += 1
            factor = la.lu_factor(matrix,overwrite_a=True,check_finite=False)
            self.clock.counts['complex_RHS'] += self.ns
            j = la.lu_solve(factor,m[:,None]*self.E,check_finite=False)
        q = self.E+self.green.apply(j) if factor is not None else self.E
        self.out = dict(m=m.copy(),j=j,Q=q,Y=self.ports.receive(j),lu=factor)
        return self.out

    def objective(self,m,Y):
        s = self.state(m)
        r = (s['Y']-Y)/np.linalg.norm(Y)
        prior = float(self.reg.evaluate(m)[0])
        return .5*float(np.vdot(r,r).real)+prior, float(np.linalg.norm(r)), s['Y']

    def data_value_gradient(self,m,Y):
        self.clock.check()
        s = self.state(m)
        scale2 = np.vdot(Y,Y).real
        r = s['Y']-Y
        a = self.ports.adjoint(r)
        if s['lu'] is not None:
            self.clock.counts['complex_RHS'] += self.ns
            a = la.lu_solve(s['lu'],a,trans=2,check_finite=False)
        prior, grad, _ = self.reg.evaluate(m)
        return .5*float(np.vdot(r,r).real)/scale2+prior, np.sum((a.conj()*s['Q']).real,axis=1)/scale2+grad


class FourierCurrents:
    """R8 Fourier APCS and its adjoint with one fixed DPCS basis per port group."""
    def __init__(self,shape,modes,groups,cost):
        from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
        self.groups = [(ids,NativeFourierAPCS(shape,modes,v,cost)) for ids,v in groups]
        self.frequencies = self.groups[0][1].frequencies
        self.k = len(self.frequencies)
        self.ns = sum(len(ids) for ids,_ in groups)
        self.n = int(np.prod(shape))

    def apply(self,coefficients):
        result = np.empty((self.n,self.ns),complex)
        for ids,op in self.groups:
            result[:,ids] = op.apply(coefficients[:,ids])
        return result

    def adjoint(self,J):
        result = np.empty((self.k,self.ns),complex)
        for ids,op in self.groups:
            result[:,ids] = op.adjoint(J[:,ids])
        return result


def dpcs(model,Y,rank):
    indices = model.ports.indices(model.ns)
    groups = {}
    for j,ids in enumerate(indices):
        groups.setdefault(tuple(ids),[]).append(j)
    JD = np.empty_like(model.E)
    bases = []
    S = model.ports.S
    # All measured ports gather rows of the same receiver kernel. One Gram
    # supplies the prescribed leading DPCS modes for every real source mask.
    gram = S@model.ports.SH
    for ids,columns in groups.items():
        model.clock.check()
        ids = np.asarray(ids)
        P = S[ids]
        H = gram[np.ix_(ids,ids)]
        k = min(rank,len(ids))
        values,u = la.eigh(H,subset_by_index=(len(ids)-k,len(ids)-1),check_finite=False)
        model.clock.counts['receiver_SVD'] += 1
        good = values>np.finfo(float).eps*max(P.shape)*max(values[-1],0.)
        sigma = np.sqrt(values[good][::-1]); u = u[:,good][:,::-1]
        v = (P.conj().T@u)/sigma
        JD[:,columns] = v@((u.conj().T@Y[:,columns])/sigma[:,None])
        bases.append((columns,v))
    return JD,bases


def cie_steps(model,Y,case,method,m):
    """Common sourcewise R8 PR; literal NIE/00/01/11 material/current operators."""
    from research.research_max_struct_r2.src.cie_global import _baseline_readout
    from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption
    from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
    cost = model.clock.counts
    modes, betas, k = stage_table(case)
    JD,groups = dpcs(model,Y,k)
    QD = model.E+model.green.apply(JD)
    nie = method=='paper_nie'
    block_current = method=='block_cie_11_sw'
    block_material = method in ('block_cie_01_sw','block_cie_11_sw')
    blocks = OneHopAbsorption(model.G,model.shape,cost) if block_material else None
    wd = 1/np.maximum(np.sum(abs(Y)**2,axis=0),1e-300)
    ws = 1/np.maximum(np.sum(abs(model.E)**2,axis=0),1e-300)
    Os,Ls = np.diag(ws),np.diag(np.sqrt(ws))
    box = np.zeros(model.n),np.full(model.n,20.)
    reg = lambda z:model.reg.evaluate(z)[:2]
    cfg = dict(aplus_regularizer=model.reg,aplus_sweeps=2,aplus_arc_candidates=True,
               aplus_inherited_gtol=True,material_inner_maxiter=60,
               material_inner_gtol=1e-8,material_inner_ftol=1e-12,
               cache_material_value_gradient=True)
    coeff = {}
    count = 0
    for stage,(mf,bm) in enumerate(zip(modes,betas)):
        model.clock.check()
        bj = bm if nie or block_current else betas[0]
        U = FourierCurrents(model.shape,mf,groups,cost)
        frequencies = [tuple(f) for f in U.frequencies]
        A = np.array([coeff.get(f,np.zeros(model.ns,complex)) for f in frequencies])
        fluctuation = U.apply(A)
        J,Q = JD+fluctuation,QD+model.green.apply(fluctuation)
        old_g = old_direction = None
        for step in range(500):
            model.clock.check()
            if block_current:
                _,d,a = blocks.parameters(m,bj)
                p = (bj/d)[:,None]
                p_apply = lambda z: p*z+a[:,None]*blocks.off_action(p*z)
                p_adj = lambda z: p.conj()*(z+blocks.off_action(a.conj()[:,None]*z,adjoint=True))
            else:
                p = (bj/(1+bj*m))[:,None]
                p_apply = lambda z:p*z
                p_adj = lambda z:p*z
            D = model.ports.receive(J)-Y
            R = p_apply(J-m[:,None]*Q)
            v = p_adj(R*ws)
            gradient = U.adjoint(model.ports.adjoint(D*wd)+v-model.green.adjoint(m[:,None]*v))
            direction = -gradient
            if old_g is not None:
                pr = np.maximum(0.,np.sum((gradient.conj()*(gradient-old_g)).real,axis=0)/
                                np.maximum(np.sum(abs(old_g)**2,axis=0),1e-300))
                direction += old_direction*pr
                bad = np.sum((gradient.conj()*direction).real,axis=0)>=0
                direction[:,bad] = -gradient[:,bad]
            dJ = U.apply(direction)
            dQ = model.green.apply(dJ)
            dd = model.ports.receive(dJ)
            dr = p_apply(dJ-m[:,None]*dQ)
            denom = wd*np.sum(abs(dd)**2,axis=0)+ws*np.sum(abs(dr)**2,axis=0)
            alpha = -np.sum((gradient.conj()*direction).real,axis=0)/np.maximum(denom,1e-300)
            change = direction*alpha
            delta2 = float(np.mean(np.sum(abs(change)**2,axis=0)/np.maximum(np.sum(abs(A)**2,axis=0),1e-300)))
            A += change
            J += dJ*alpha
            Q += dQ*alpha
            if block_material:
                m,_ = mobius_material_readout(blocks,J,Q,m,bm,Os,Ls,1.,box,reg,cfg,cost)
            elif nie:
                B = Q+bm*J
                z = bm*np.sum((B.conj()*J).real*ws,axis=1)/np.maximum(np.sum(abs(B)**2*ws,axis=1),1e-300)
                z = np.clip(z,0.,bm*20/(1+bm*20))
                m = z/(bm*(1-z))
            else:
                m,_ = _baseline_readout(J,Q,m,bm,ws,1.,box,reg,True,cfg,cost)
            old_g,old_direction = gradient,direction
            count += 1
            cost['PR_CG_steps'] += 1
            end = step==499 or (step>=1 and delta2<(1e-4 if stage==0 else 1e-6))
            model.clock.check()
            yield m,count,stage,end
            if end:
                break
        coeff = {f:a.copy() for f,a in zip(frequencies,A)}


class ScaledFFT:
    def __init__(self,green,scale,adjoint=False):
        self.green,self.scale,self.is_adjoint = green,scale,adjoint

    def __matmul__(self,v):
        return self.scale*self.green.apply(v,adjoint=self.is_adjoint)

    def conj(self):
        return self

    @property
    def T(self):
        return ScaledFFT(self.green,self.scale,not self.is_adjoint)


def make_projector(model,Y,delta):
    """Exact R14 graph-ball projection, shared receiver-row solves for masks."""
    r15 = R14.parent/'research_max_struct_r15/runtime/field_current.py'
    r12 = R14.parent/'research_max_struct_r12/src/field_current.py'
    path = r15 if r15.is_file() else r12
    mod = module_at('killtest_field_current',path)
    P = object.__new__(mod.DenseFieldCurrentProjectors)
    P.E,P.G,P.cells,P.sources,P.components = model.E,model.G,model.n,model.ns,1
    P.scale,P.delta,P.grouped_ports = 6.,float(delta),True
    P.lower,P.upper = np.zeros(model.n),np.full(model.n,20.)
    P.costs = model.clock.counts
    B = np.asfortranarray(P.scale*model.G)
    H = la.blas.zherk(1.,B,trans=2,lower=True)
    H.flat[::model.n+1] += 1.
    model.clock.counts['H_cholesky'] += 1
    P.H = la.cho_factor(H,lower=True,overwrite_a=True,check_finite=False)
    P.B = ScaledFFT(model.green,P.scale)
    model.clock.check()
    T = P.scale*model.ports.S
    indices = model.ports.indices(model.ns)
    active = np.unique(np.concatenate(indices))
    V = la.cho_solve(P.H,T[active].conj().T,check_finite=False)
    model.clock.counts['H_setup_rhs'] += len(active)
    remap = np.zeros(len(T),int)
    remap[active] = np.arange(len(active))
    P.ports,P.targets,P.port_groups = [],list(Y.T),[]
    groups = {}
    for j,ids in enumerate(indices):
        groups.setdefault(tuple(ids),[]).append(j)
    floor2 = 0.
    for rows,cols in groups.items():
        model.clock.check()
        i = np.asarray(rows)
        z = V[:,remap[i]]
        gram = T[i]@z
        values,U = la.eigh((gram+gram.conj().T)*.5,check_finite=False)
        model.clock.counts['port_eigh'] += 1
        tol = 32*np.finfo(float).eps*max(len(i),model.n)*max(np.max(abs(values)),np.finfo(float).tiny)
        pos = values>tol
        p = mod._Port(T[i],z,np.where(pos,values,0.),U,pos)
        P.ports.append(p)
        P.port_groups.append((p,np.array(cols)))
        floor2 += float(np.linalg.norm(U[:,~pos].conj().T@Y[:,cols])**2)
    P.minimum_reachable_residual = float(np.sqrt(floor2))
    P._energy_scale = max(float(np.vdot(Y,Y).real),P.delta**2,np.finfo(float).tiny)
    P._ball_tolerance = 128*np.finfo(float).eps*P._energy_scale
    model.clock.counts['DR_setup_wall'] = model.clock.wall()
    model.clock.counts['DR_ball_floor_rel'] = P.minimum_reachable_residual/np.linalg.norm(Y)
    if P.delta**2<floor2-P._ball_tolerance:
        raise mod.UnreachableDataBall(P.delta,P.minimum_reachable_residual)
    P.last_A_info = {}
    model.clock.check()
    return P


def dr_steps(model,Y,delta,m):
    P = make_projector(model,Y,delta)
    s = model.state(m)
    x,q = s['j']/P.scale,s['Q'].copy()
    for i in range(3000):
        model.clock.check()
        cx,cq,cm = P.material(x,q)
        ax,aq = P.affine_ball(2*cx-x,2*cq-q)
        x,q = x+ax-cx,q+aq-cq
        model.clock.check()
        yield cm,i+1,0,i==2999


def count_pack(counts):
    fft = sum(counts.get(k,0) for k in ('Green_kernel_FFT_setup','Green_FFT_2D_transforms',
                'Green_IFFT_2D_transforms','coefficient_FFT_2D_transforms','coefficient_IFFT_2D_transforms'))
    return dict(LU=int(counts.get('LU',counts.get('full_LU_factorizations',0))),
                FFT=int(fft), complex_RHS=int(counts.get('complex_RHS',counts.get('linear_solve_rhs',0))),
                H_cholesky=int(counts.get('H_cholesky',0)),H_rhs=int(counts.get('H_setup_rhs',0)+counts.get('H_iteration_rhs',0)),
                reduced_cholesky=int(counts.get('small_cholesky',0)),
                receiver_SVD=int(counts.get('receiver_SVD',0)),
                DR_setup_wall=counts.get('DR_setup_wall'), DR_ball_floor_rel=counts.get('DR_ball_floor_rel'))


def snapshot(clock,m,pred,value,res,step,stage=0,wall_offset=0.,cpu_offset=0.,base_counts=None):
    costs = clock.counts if base_counts is None else Counter(base_counts)+clock.counts
    item = dict(m=m.copy(),pred=pred.copy(),phi=float(value),residual=float(res),step=step,stage=stage,
                available_wall=wall_offset+clock.wall(),available_CPU=cpu_offset+clock.cpu(),counts=dict(costs))
    return item


def best_before(history,limit):
    ready = [p for p in history if p['available_wall']<=limit and np.isfinite(p['phi'])]
    return min(ready,key=lambda p:p['phi']) if ready else None


def run_light(model,Y,seed,used_wall,used_cpu,base_counts,budget):
    clock = Clock(budget-used_wall)
    model.clock = clock
    model.green.cost = clock.counts
    model.out = None
    points = [seed]
    iterations = 0
    status = 'CORRECTOR_DONE'
    def callback(x):
        nonlocal iterations
        clock.check()
        iterations += 1
        v,r,p = model.objective(x,Y)
        points.append(snapshot(clock,x,p,v,r,iterations,wall_offset=used_wall,cpu_offset=used_cpu,base_counts=base_counts))
        clock.check()
    try:
        result = minimize(lambda m:model.data_value_gradient(m,Y),seed['m'],jac=True,
                         method='L-BFGS-B',bounds=[(0.,20.)]*model.n,callback=callback,
                         options=dict(maxiter=20,maxfun=80,maxls=20,ftol=1e-12,gtol=1e-8))
        status = str(result.message)
    except TimeLimit:
        status='WALL_LIMIT'
    except Exception as exc:
        status=f'{type(exc).__name__}: {exc}'
    elapsed = used_wall+clock.wall()
    result = best_before(points,budget)
    return result,dict(wall=elapsed,CPU=used_cpu+clock.cpu(),status=status,iterations=iterations,
                       overrun=max(0.,elapsed-budget),cost=Counter(base_counts)+clock.counts),points


def solve_a_arm(public,case,arm):
    clock = Clock(30.)
    clock.watch(12.,15.,24.)
    history = []
    m = np.zeros(case['side']**2)
    last_step,last_stage = 0,-1
    model = None
    status = 'SCHEDULE_COMPLETE'
    Y = public['Y']
    try:
        model = PaperModel(public,case,clock)
        v,r,p = model.objective(m,Y)
        history.append(snapshot(clock,m,p,v,r,0,-1))
        # Same audit time grid for all arms. Audits consume proposal budget.
        due = time.perf_counter()
        sequence = dr_steps(model,Y,public['delta'],m) if arm=='dr_standalone' else cie_steps(model,Y,case,arm,m)
        for candidate,last_step,last_stage,end in sequence:
            clock.check()
            if time.perf_counter()>=due or end:
                v,r,p = model.objective(candidate,Y)
                history.append(snapshot(clock,candidate,p,v,r,last_step,last_stage))
                due = time.perf_counter()+A_AUDIT_PERIOD
            clock.check()
    except TimeLimit:
        status = 'WALL_LIMIT'
    except Exception as exc:
        status = f'FAILED: {type(exc).__name__}: {exc}'
    clock.capture()
    raw_elapsed,raw_cpu,raw_counts = clock.wall(),clock.cpu(),dict(clock.counts)
    # Snapshot all prefix costs before taking either corrector branch.
    # If the proposal stops early, charge only its actual completed work.
    ended = (raw_elapsed,raw_cpu,raw_counts)
    spent = {lim:clock.cuts.get(lim) or ended for lim in (12.,15.,24.,30.)}
    rows = []
    for T,target in zip(A_BUDGETS,LIGHT_SEEDS):
        raw = best_before(history,T)
        seed = best_before(history,target)
        used,cpu_spent,raw_cost = spent[T]
        rows.append((T,'raw',raw,dict(wall=used,CPU=cpu_spent,
                        status=status if raw_elapsed<=T or T==30. else 'PAID_PREFIX',iterations=last_step if T==30. else None,
                        overrun=max(0.,used-T),cost=raw_cost)))
        if model is not None and seed is not None:
            pref_wall,pref_cpu,pref_cost = spent[target]
            best,rec,_ = run_light(model,Y,seed,pref_wall,pref_cpu,pref_cost,T)
            rec.update(seed_wall=seed['available_wall'],seed_step=seed['step'],seed_phi=seed['phi'],
                       proposal_wall=pref_wall)
        else:
            best,rec = seed,dict(wall=used,CPU=cpu_spent,status='PROPOSAL_FAILED_OR_NO_SEED',iterations=0,
                                 overrun=max(0.,used-T),cost=raw_cost)
        rows.append((T,'light',best,rec))
    return rows,history,dict(wall=raw_elapsed,CPU=raw_cpu,status=status,iterations=last_step,stage=last_stage,cost=raw_counts)


def record_row(arm,T,view,item,work,**extra):
    row = dict(method=arm,budget=T,view=view,wall_s=work['wall'],cpu_s=work['CPU'],
               status=work['status'],iterations=work['iterations'],overrun_s=work['overrun'],
               seed_wall=work.get('seed_wall'),seed_step=work.get('seed_step'),proposal_wall=work.get('proposal_wall'),
               **count_pack(work['cost']),**extra)
    if item is not None:
        row.update(phi=item['phi'],relative_data_residual=item['residual'],
                   candidate_wall=item['available_wall'],candidate_step=item['step'],candidate_stage=item['stage'])
    return row


def write_predictions(path, arrays, rows):
    arrays['rows_json'] = np.array(json.dumps(rows,ensure_ascii=False,default=lambda v:v.item()))
    temp = path.with_suffix('.tmp.npz')
    np.savez_compressed(temp,**arrays)
    temp.replace(path)


def solve_a_case(case):
    p = RUNTIME/'pred'/f'{case["case"]}.npz'
    if p.exists():
        return case['case']
    from research.research_max_struct_r8.src import paper_physics
    from research.research_max_struct_r2.src import cie_global
    from research.research_max_struct_r4.src import cie_onehop
    from research.research_max_struct_r7.src import aplus_readout
    from research.research_max_struct_r6.src import fft_tsom_ops
    path = R14.parent/'research_max_struct_r15/runtime/field_current.py'
    if path.is_file():
        module_at('killtest_field_current',path)
    else:
        module_at('killtest_field_current',R14.parent/'research_max_struct_r12/src/field_current.py')
    public = array_file(public_path(case))
    case = json.loads(str(public['case_json']))
    arrays,rows = {},[]
    for arm in A_ARMS:
        result,history,report = solve_a_arm(public,case,arm)
        if history:
            arrays[f'{arm}_raw_audit_material'] = np.stack([v['m'] for v in history])
            arrays[f'{arm}_raw_audit'] = np.array([[v[k] for k in ('available_wall','available_CPU','phi','residual','step','stage')] for v in history])
        for T,view,value,work in result:
            key = f'{arm}_{int(T)}_{view}'
            row = record_row(arm,T,view,value,work,case=case['case'],suite=case['suite'],key=key,
                             proposal_status=report['status'],proposal_iterations=report['iterations'],
                             native_threads=1,host=platform.node())
            if value is not None:
                arrays[key] = value['m']
                arrays[key+'_Y'] = value['pred']
            rows.append(row)
        print(f'{case["case"]} {arm}: {report["status"]}, {report["iterations"]} proposal steps, '
              f'{report["wall"]:.2f}s',flush=True)
    write_predictions(p,arrays,rows)
    return case['case']


def digits_model(public,aperture,clock):
    sys.path.insert(0,str(A24/'src'))
    digits = A17/'DIGITS_2D_20261030/src/physics.py'
    if digits.is_file():
        cls = module_at('killtest_digits_physics',digits).DigitsPhysics
    else:
        # The custom geometry below uses only the inherited A24 state/solve
        # methods, so the public A24 checkout is sufficient for reproduction.
        from a24.physics import Physics2D as cls
    model = object.__new__(cls)
    model.book,model.config,model.grid = clock,{'backward_residual_tolerance':1e-10},24
    model.n = model.N = len(public['centers'])
    model.points,model.G = public['centers'],public['G']
    model.volume = float(public['volume'])
    model.sqrt_volume = np.sqrt(model.volume)
    model.k = float(public['k'])
    model.components = 1
    model.incident = model.sqrt_volume*public['E']
    suffix = '_limited' if aperture=='limited' else ''
    model.S = public['S'+suffix]/model.sqrt_volume
    model.q,model.s = model.S.shape[0],model.incident.shape[1]
    return model,public['Y'+suffix]


def orth_extend(U,bank,k):
    """Deterministic MGS in column order, preserving the informative columns."""
    for v in bank.T:
        if U.shape[1]>=k:
            break
        u = v.copy()
        for _ in range(2):
            u -= U@(U.conj().T@u)
        norm = np.linalg.norm(u)
        if norm>1e-10*np.linalg.norm(v):
            U = np.column_stack((U,u/norm))
    return U


def receiver_basis(S,points,k,clock):
    clock.check()
    _,s,vh = la.svd(S,full_matrices=False,check_finite=False)
    clock.counts['receiver_SVD'] += 1
    n_visible = int(np.count_nonzero(s>max(S.shape)*np.finfo(float).eps*s[0]))
    R = vh[:n_visible].conj().T
    fq = sorted(((i,j) for i in range(-4,5) for j in range(-4,5)),key=lambda ij:(ij[0]**2+ij[1]**2,ij))
    fourier = np.exp(2j*np.pi*(points@np.asarray(fq).T))
    # Ten independent received rows cannot supply twelve or sixteen nonzero
    # right singular directions. Smooth deterministic null completion is common.
    U = orth_extend(R[:,:k],fourier,k)
    return U,R,fourier,n_visible,float(s[-1]/s[0])


def internal_basis(G,R,k,clock):
    U = R[:,:min(k//2,R.shape[1])]
    project = lambda x:x-U@(U.conj().T@x)
    rng = np.random.default_rng(A17_SEED)
    test = rng.standard_normal((len(G),k+12))+1j*rng.standard_normal((len(G),k+12))
    Q,_ = la.qr(project(test),mode='economic',check_finite=False)
    for _ in range(2):
        clock.check()
        Q,_ = la.qr(project(G.conj().T@(G@Q)),mode='economic',check_finite=False)
        clock.counts['internal_G_rhs'] += 2*Q.shape[1]
    _,_,vh = la.svd(G@Q,full_matrices=False,check_finite=False)
    clock.counts['internal_G_rhs'] += Q.shape[1]
    return orth_extend(U,Q@vh.conj().T,k)


def local_transfer(L,S,U):
    UH = U.conj().T
    small = UH@(L@U)
    return (S@U)@la.solve(small,UH,assume_a='gen',check_finite=False)


def complex_endpoint(T,field,r,lam,ell,scale,clock):
    # Complex equivalent of the original packed-real damped endpoint. Same
    # mass coordinates and penalty on both real and imaginary components.
    J = (T[:,None,:]*field.T[None,:,:]/scale).reshape(r.size,-1)
    base = -ell/lam
    H = J@J.conj().T
    H.flat[::len(H)+1] += lam
    clock.check()
    clock.counts['small_cholesky'] += 1
    return base-J.conj().T@la.cho_solve(la.cho_factor(H,lower=True,overwrite_a=True,check_finite=False),
                                       r.ravel()+J@base,check_finite=False)


def clip_mass(x,h):
    return np.clip(x.real,-.5*h,4.*h)+1j*np.clip(x.imag,0.,.5*h)


def norm2(x):
    return float(np.vdot(x,x).real)


def quadratic(r_action,d,lam,ell):
    return .5*norm2(r_action)+.5*lam*norm2(d)+float(np.vdot(ell,d).real)


def top_exchange(model,state,U,R,fourier,x,delta,T,field,r,lam,ell,scale,clock):
    """Screen one-column swaps by the exact Galerkin Schur term and one projected
    quadratic line step. No full tangent or endpoint solve for rejected swaps.
    """
    clock.check()
    G,S,L = model.G,model.S,state.L
    cut = int(np.argmin(np.sum(abs(S@U)**2,axis=0)))
    V = np.delete(U,cut,axis=1)
    VH,LV = V.conj().T,L@V
    B = VH@LV
    base_factor = la.lu_factor(B,check_finite=False)
    C = la.lu_solve(base_factor,VH,check_finite=False)
    SV = S@V
    Tv = SV@C
    rng = np.random.default_rng(A17_SEED)
    k,n = U.shape[1],model.n
    bank = np.column_stack([R[:,k:k+2],G@R[:,4:6],fourier[:,1:3],
                            state.current[:,:2],S.conj().T@r[:,:2],
                            rng.standard_normal((n,2))+1j*rng.standard_normal((n,2))])
    bank -= V@(VH@bank)
    lengths = np.linalg.norm(bank,axis=0)
    nonzero = lengths>1e-11
    bank = bank[:,nonzero]/lengths[nonzero]
    LC = L@bank
    right = bank.conj().T-(bank.conj().T@LV)@C
    column = bank-V@(C@LC)
    schur = np.sum(bank.conj()*(LC-LV@(C@LC)),axis=0)
    observed = S@column
    existing = r+T@(delta[:,None]*field)/scale
    base_value = quadratic(existing,delta,lam,ell)
    best = None
    for i in range(bank.shape[1]):
        if abs(schur[i])<1e-12:
            continue
        clock.check()
        # Low-rank port transfer: T_new = T_V + O_i B_i, exactly at fixed m.
        op,ip = observed[:,i]/schur[i], right[i]
        Ti = Tv+op[:,None]*ip[None,:]
        rc = r+Ti@(delta[:,None]*field)/scale
        g = np.sum(field.conj()*(Ti.conj().T@rc),axis=1)/scale+lam*delta+ell
        diagonal = (np.sum(abs(Ti)**2,axis=0)*np.sum(abs(field)**2,axis=1)/scale**2)+lam
        d = clip_mass(x+delta-g/diagonal,model.sqrt_volume)-(x+delta)
        Jd = Ti@(d[:,None]*field)/scale
        slope = float(np.vdot(g,d).real)
        alpha = float(np.clip(-slope/(norm2(Jd)+lam*norm2(d)+1e-300),0.,1.))
        estimate = base_value-quadratic(rc+alpha*Jd,delta+alpha*d,lam,ell)
        clock.counts['screened_swaps'] += 1
        if best is None or estimate>best[0]:
            best = (float(estimate),np.column_stack((V,bank[:,i])),Ti)
    return best


def run_b_arm(public,case,aperture,k,arm,wall_limit=B_BUDGET):
    clock = Clock(wall_limit)
    model,Y = digits_model(public,aperture,clock)
    scale = np.linalg.norm(Y)
    chi0 = np.full(model.n,.1+.01j)
    x0 = model.sqrt_volume*chi0
    history,trace = [],[]
    status,events,accepts = 'MAX_UPDATES',0,0
    n_visible,s_ratio = None,None
    init_wall = None
    def keep(state,value,r,step):
        history.append(snapshot(clock,state.chi,state.data,value,np.linalg.norm(r),step))
    try:
        state = model.full_state(chi0)
        x = x0.copy()
        r = (state.data-Y)/scale
        value = .5*norm2(r)
        keep(state,value,r,0)
        if arm!='full_gradient_matched_time':
            U,R,fourier,n_visible,s_ratio = receiver_basis(model.S,model.points,k,clock)
            if arm=='projected_internal':
                U = internal_basis(model.G,R,k,clock)
        init_wall = clock.wall()
        max_updates = 300 if arm=='full_gradient_matched_time' else 18
        for it in range(max_updates):
            clock.check()
            if np.linalg.norm(r)<=.022:
                status = 'DISCREPANCY'
                break
            field = state.field/model.sqrt_volume
            lam = A17_PRIOR+.01*.3**(it//3)
            ell = A17_PRIOR*(x-x0)
            full = arm=='full_gradient_matched_time'
            if full:
                T = model.solve_factor(state.factor,model.S.conj().T,adjoint=True,role='full_J_receiver_adjoint').conj().T
                g = np.sum(field.conj()*(T.conj().T@r),axis=1)/scale+ell
            else:
                T = local_transfer(state.L,model.S,U)
                load = model.solve_factor(state.factor,model.S.conj().T@r,adjoint=True,role='material_adjoint')
                g = np.sum(field.conj()*load,axis=1)/scale+ell
            step = complex_endpoint(T,field,r,lam,ell,scale,clock)
            delta = clip_mass(x+step,model.sqrt_volume)-x
            row = dict(step=it+1,prediction_remainder=None,projection_fraction=float(np.linalg.norm(step-delta)/max(np.linalg.norm(step),1e-300)),
                       prior=float(.5*A17_PRIOR*norm2(x-x0)),projection_active_cells=int(np.count_nonzero(abs(step-delta)>1e-12)),
                       accepted_swap=False,verified=False,relative_data_residual=float(np.linalg.norm(r)))
            exact_action = None
            if arm=='a17_proj_fast' and (it+1)%4==0:
                events += 1
                row['exchange_event'] = events
                begin = clock.wall()
                top = top_exchange(model,state,U,R,fourier,x,delta,T,field,r,lam,ell,scale,clock)
                if top is not None:
                    score,V,Tnew = top
                    trial = complex_endpoint(Tnew,field,r,lam,ell,scale,clock)
                    candidate = clip_mass(x+trial,model.sqrt_volume)-x
                    directions = np.c_[delta[:,None]*field,candidate[:,None]*field]
                    linear = model.S@model.solve_factor(state.factor,directions,role='projected_exchange')/scale
                    db,dc = np.split(linear,2,axis=1)
                    gain = quadratic(r+db,delta,lam,ell)-quadratic(r+dc,candidate,lam,ell)
                    reduced_b = T@(delta[:,None]*field)/scale
                    reduced_c = Tnew@(candidate[:,None]*field)/scale
                    approx_unprojected = quadratic(r+T@(step[:,None]*field)/scale,step,lam,ell)-quadratic(r+Tnew@(trial[:,None]*field)/scale,trial,lam,ell)
                    approx_projected = quadratic(r+reduced_b,delta,lam,ell)-quadratic(r+reduced_c,candidate,lam,ell)
                    row.update(verified=True,prescreen_score=score,Q_proj=float(gain),
                               reduced_unprojected_gain=approx_unprojected,reduced_projected_gain=approx_projected,
                               base_projection_fraction=row['projection_fraction'],
                               swap_projection_fraction=float(np.linalg.norm(trial-candidate)/max(np.linalg.norm(trial),1e-300)),
                               base_transfer_relative_error=float(np.linalg.norm(reduced_b-db)/max(np.linalg.norm(db),1e-300)),
                               swap_transfer_relative_error=float(np.linalg.norm(reduced_c-dc)/max(np.linalg.norm(dc),1e-300)))
                    exact_action = db
                    if gain>1e-12*max(1.,value):
                        U,step,delta,exact_action = V,trial,candidate,dc
                        row['accepted_swap'] = True
                        accepts += 1
                row['exchange_wall_s'] = clock.wall()-begin
            row['projection_fraction'] = float(np.linalg.norm(step-delta)/max(np.linalg.norm(step),1e-300))
            slope = float(np.vdot(g,delta).real)
            if slope>=0. or np.linalg.norm(delta)<1e-14:
                delta = clip_mass(x-g/max(1.,np.linalg.norm(g)),model.sqrt_volume)-x
                slope = float(np.vdot(g,delta).real)
                exact_action = None
                row['fallback'] = True
            row['slope'] = slope
            row['halvings'] = None
            if slope>=0. or np.linalg.norm(delta)<1e-14:
                status = 'PROJECTED_STATIONARY'
                trace.append(row)
                break
            for half in range(17):
                clock.check()
                alpha = 2.**(-half)
                xtrial = x+alpha*delta
                new = model.full_state(xtrial/model.sqrt_volume)
                rnew = (new.data-Y)/scale
                vnew = .5*norm2(rnew)+.5*A17_PRIOR*norm2(xtrial-x0)
                clock.check()
                if half==0:
                    row['unit_trial_phi_change'] = vnew-value
                if vnew<=value+1e-4*alpha*slope:
                    if exact_action is not None:
                        remainder = rnew-r-alpha*exact_action
                        row['prediction_remainder'] = float(np.linalg.norm(remainder)/max(alpha*np.linalg.norm(exact_action),1e-300))
                    row.update(halvings=half,accepted_phi_change=vnew-value,
                               lower_active_cells=int(np.sum(new.chi.real<=-.5+1e-8)),
                               upper_active_cells=int(np.sum(new.chi.real>=4.-1e-8)),
                               loss_lower_active_cells=int(np.sum(new.chi.imag<=1e-8)))
                    state,x,r,value = new,xtrial,rnew,vnew
                    keep(state,value,r,it+1)
                    break
            trace.append(row)
            if row['halvings'] is None:
                status = 'LINE_SEARCH_REJECTED'
                break
            if it>=8 and alpha*np.linalg.norm(delta)<=1e-6*max(1.,np.linalg.norm(x)):
                status = 'SMALL_STEP'
                break
    except TimeLimit:
        status = 'WALL_LIMIT'
    except Exception as exc:
        status = f'FAILED: {type(exc).__name__}: {exc}'
    elapsed = clock.wall()
    result = best_before(history,wall_limit)
    work = dict(wall=elapsed,CPU=clock.cpu(),status=status,iterations=len(history)-1,
                overrun=max(0.,elapsed-wall_limit),cost=dict(clock.counts),
                init_wall=init_wall,events=events,accepts=accepts,n_visible=n_visible,
                tail_ratio=s_ratio)
    return result,history,trace,work


def collect_b_rows(case,arrays,rows,aperture,k,arm,run,cut,full_limit=None):
    value,history,trace,work = run
    root = f'{aperture}_k{k}_{arm}'
    if history:
        arrays[root+'_trace_m'] = np.stack([v['m'] for v in history])
        arrays[root+'_trace'] = np.array([[v[x] for x in ('available_wall','available_CPU','phi','residual','step')] for v in history])
    arrays[root+'_diagnostics_json'] = np.array(json.dumps(trace))
    for view,pick in [('native',value),('matched',best_before(history,cut))]:
        key = f'{root}_{view}'
        if view=='native' or work['wall']<=cut:
            use = work
        else:
            # Scope of counter and CPU here is the last full evaluation before
            # the cut. The native row retains all attempted work and overruns.
            past = [v for v in history if v['available_wall']<=cut]
            last = past[-1] if past else None
            use = dict(wall=cut,CPU=last['available_CPU'] if last else None,status='RECEIVER_WALL_PREFIX',
                       iterations=last['step'] if last else None,overrun=0.,cost=last['counts'] if last else {})
        row=record_row(arm, (B_BUDGET if full_limit is None else full_limit) if view=='native' else cut,view,pick,use,
                       case=case['case'],suite='2A',aperture=aperture,rank=k,key=key,
                       native_threads=1,host=platform.node(),proposal_status=work['status'],
                       setup_wall=work['init_wall'],receiver_nonnull_rank=work['n_visible'],
                       receiver_tail_ratio=work['tail_ratio'],exchange_events=work['events'],
                       accepted_exchanges=work['accepts'],
                       receiver_completion=max(0,k-(work['n_visible'] or 0)) if arm in ('receiver_svd','a17_proj_fast') else None)
        if pick is not None:
            arrays[key] = pick['m']
            arrays[key+'_Y'] = pick['pred']
        rows.append(row)


def solve_b_case(case,replace_step_capped_full=False):
    p = RUNTIME/'pred'/f'{case["case"]}.npz'
    if p.exists() and not replace_step_capped_full:
        return case['case']
    public = array_file(public_path(case))
    arrays,rows = {},[]
    full = 'full_gradient_matched_time'
    if replace_step_capped_full:
        old = array_file(p)
        rows = [r for r in json.loads(str(old.pop('rows_json'))) if r['method']!=full]
        arrays = {k:v for k,v in old.items() if full not in k}
    for aperture in ('full','limited'):
        runs = {}
        if not replace_step_capped_full:
            for k in (12,16):
                for arm in B_ARMS[:3]:
                    runs[k,arm] = run_b_arm(public,case,aperture,k,arm)
            cuts = {k:min(B_BUDGET,runs[k,'receiver_svd'][3]['wall']) for k in (12,16)}
            cap = min(B_BUDGET,max(r[3]['wall'] for r in runs.values()))
        else:
            native = [r for r in rows if r['view']=='native' and r['aperture']==aperture]
            cuts = {k:min(B_BUDGET,next(r['wall_s'] for r in native if r['method']=='receiver_svd' and r['rank']==k)) for k in (12,16)}
            cap = min(B_BUDGET,max(r['wall_s'] for r in native))
        exact = run_b_arm(public,case,aperture,12,full,wall_limit=cap)
        for k,cut in cuts.items():
            if not replace_step_capped_full:
                for arm in B_ARMS[:3]:
                    collect_b_rows(case,arrays,rows,aperture,k,arm,runs[k,arm],cut)
            collect_b_rows(case,arrays,rows,aperture,k,full,exact,cut,full_limit=cap)
        print(f'{case["case"]} {aperture}: full gradient {exact[3]["iterations"]} steps, {exact[3]["wall"]:.2f}s, receiver limits {cuts}',flush=True)
    write_predictions(p,arrays,rows)
    return case['case']


def ensure_dirs():
    for name in ('public','heldout','pred'):
        (RUNTIME/name).mkdir(parents=True,exist_ok=True)


def solve_one(case):
    return solve_b_case(case) if case['suite']=='2A' else solve_a_case(case)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',nargs='?',default='all',choices=('all','prepare','solve','score'))
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--campaign',choices=('A','B','both'),default='both')
    ap.add_argument('--case',help='Execute only this case during prepare/solve')
    args = ap.parse_args()
    ensure_dirs()
    case_list = [c for c in specs() if (not args.case or c['case']==args.case) and
                 (args.campaign=='both' or (c['suite']=='2A')==(args.campaign=='B'))]
    if args.mode in ('all','prepare'):
        run_jobs(prepare_one,[(c,) for c in case_list],args.workers)
    if args.mode in ('all','solve'):
        run_jobs(solve_one,[(c,) for c in case_list],args.workers)
    if args.mode in ('all','score'):
        score_all(args.workers)


def coarse_truth(case,centers):
    # Area averages of the specified indicators, sampled on an 8x8 subcell rule.
    h = case['width']/case['side']
    u = ((np.arange(8)+.5)/8.-.5)*h
    a,b = np.meshgrid(u,u)
    offsets = np.c_[a.ravel(),b.ravel()]
    z = (centers[:,None,:]+offsets).reshape(-1,2)
    return material(case,z).reshape(len(centers),-1).mean(axis=1)


def otsu_support(a):
    positive = np.maximum(a.real,0.)
    if np.ptp(positive)==0:
        return positive>0
    w,edges = np.histogram(positive,bins=256)
    centers = (edges[1:]+edges[:-1])*.5
    left = np.cumsum(w)
    right = left[-1]-left
    mass = np.cumsum(w*centers)
    m = np.divide(mass,left,out=np.zeros_like(mass),where=left>0)
    n = np.divide(mass[-1]-mass,right,out=np.zeros_like(mass),where=right>0)
    threshold = centers[np.argmax(left*right*(m-n)**2)]
    return positive>threshold


def iou(a,b):
    union = np.sum(a|b)
    return float(np.sum(a&b)/union) if union else 1.


def errors(m,truth):
    diff = m-truth
    level = .5*np.max(truth.real)
    mask = truth.real>level
    return dict(relativeL2=float(np.linalg.norm(diff)/np.linalg.norm(truth)),
                amplitudeNRMSE=float(np.sqrt(norm2(diff)/len(diff))/np.max(np.abs(truth))),
                realRelativeL2=float(np.linalg.norm(diff.real)/np.linalg.norm(truth.real)),
                lossRelativeL2=float(np.linalg.norm(diff.imag)/np.linalg.norm(truth.imag)) if np.iscomplexobj(truth) else None,
                OtsuIoU=iou(otsu_support(m),mask),half_contrast_IoU=iou(m.real>level,mask))


def score_case(case):
    s = case['case']; name = s+'.npz'
    data = array_file(public_path(case))
    case = json.loads(str(data['case_json']))
    a = array_file(RUNTIME/'pred'/name)
    rows = json.loads(str(a.pop('rows_json')))
    expected = 32 if case['suite']=='2A' else 20
    if len(rows)!=expected:
        raise ValueError(f'{s}: {len(rows)}/{expected} terminal rows')
    a['case_json'],a['centers'] = data['case_json'],data['centers']
    refs = {}
    if case['suite']=='1C':
        held = array_file(RUNTIME/'heldout'/name)
        valid = np.isfinite(held['Y'])
        hm = PaperModel(dict(data,E=held['E'],S=held['S'],rx_ids=None),case,Clock(np.inf))
        a['heldout_mask'] = valid
        a['train_rx_ids'] = data['rx_ids']
        cache = []
    else:
        reference = coarse_truth(case,data['centers'])
        a['truth_coarse'] = reference
        a['Y_train'] = data['Y']
        if case['suite']=='2A':
            a['Y_train_limited'] = data['Y_limited']
            a['unit_step_diagnostics_json'] = np.array(json.dumps(diagnose_b_unit_steps(data,a)))
            for aperture in ('full','limited'):
                model,y = digits_model(data,aperture,Clock(np.inf))
                state = model.full_state(reference)
                refs[aperture] = np.linalg.norm(state.data-y)/np.linalg.norm(y)
        else:
            model = PaperModel(data,case,Clock(np.inf))
            refs['A'] = model.objective(reference,data['Y'])[1]
    for r in rows:
        r.update(fine_side=case['fine'],inverse_side=case['side'],noise_relative=case['noise_relative'],
                 frequency_Hz=case['frequency'],data_preparation_wall_s=case['data_wall'],
                 fine_state_residual=case.get('fine_state_residual'))
        if r['key'] not in a:
            r['score_status'] = 'NO_ELIGIBLE_PREDICTION'
            continue
        m = a[r['key']]
        r.update(score_status='SCORED',contrast_min=float(m.real.min()),contrast_max=float(m.real.max()))
        if case['suite']=='1C':
            before = hm.clock.wall()
            reuse = next((v for z,v in cache if np.array_equal(z,m)),None)
            pred = reuse if reuse is not None else hm.state(m)['Y'].copy()
            if reuse is None:
                cache.append((m,pred))
            a[r['key']+'_held_Y'] = pred
            r.update(heldout_relative_complex_error=float(np.linalg.norm((pred-held['Y'])[valid])/np.linalg.norm(held['Y'][valid])),
                     train_channels=case['train_channels'],held_channels=case['held_channels'],
                     heldout_scoring_wall_s=hm.clock.wall()-before,
                     empty_calibration_median=case['calibration_median_validation'],
                     empty_calibration_max=case['calibration_max_validation'],
                     data_ball_ratio=case['data_ball_ratio'])
        else:
            r.update(errors(m,reference),reference_data_residual=float(refs[r.get('aperture','A')]))
            if case['suite']!='2A':
                r['data_ball_ratio'] = float(data['delta']/np.linalg.norm(data['Y']))
    return case,rows,a


def diagnose_b_unit_steps(public,bundle):
    """Reconstruct the full trial ray from saved endpoints and line-search step size."""
    clock = Clock(np.inf)
    model,Y = digits_model(public,'limited',clock)
    scale = np.linalg.norm(Y)
    x0 = model.sqrt_volume*np.full(model.n,.1+.01j)
    report = {}
    for k in (12,16):
        key = f'limited_k{k}_a17_proj_fast'
        trace = bundle[key+'_trace']
        material = dict(zip(trace[:,4].astype(int),bundle[key+'_trace_m']))
        events = json.loads(str(bundle[key+'_diagnostics_json']))
        checks = []
        for e in events:
            i = e['step']
            if not e['accepted_swap'] or e.get('fallback') or i not in material:
                continue
            alpha = 2.**(-e['halvings'])
            a,b = material[i-1],material[i]
            unit = clip_mass(model.sqrt_volume*(a+(b-a)/alpha),model.sqrt_volume)/model.sqrt_volume
            base = model.full_state(a)
            exact_tangent = model.S@model.solve_factor(
                base.factor,(unit-a)[:,None]*base.field,role='unit_trial_diagnostic')/scale
            trial = model.full_state(unit)
            r = (base.data-Y)/scale
            exact = (trial.data-base.data)/scale
            prior_change = .5*A17_PRIOR*(norm2(model.sqrt_volume*unit-x0)-norm2(model.sqrt_volume*a-x0))
            linear_gain = .5*(norm2(r+exact_tangent)-norm2(r))+prior_change
            actual_gain = .5*(norm2(r+exact)-norm2(r))+prior_change
            checks.append(dict(step=i,alpha=alpha,tangent_model_phi_change=linear_gain,
                               actual_phi_change=actual_gain,
                               recorded_phi_difference=actual_gain-e['unit_trial_phi_change'],
                               nonlinear_data_remainder=float(np.linalg.norm(exact-exact_tangent)/max(np.linalg.norm(exact_tangent),1e-300))))
        report[str(k)] = checks
    report['cost'] = dict(wall_s=clock.wall(),cpu_s=clock.cpu(),**count_pack(clock.counts))
    return report


def numerical_summary(rows,bundles):
    def lookup(case,method,view,budget=None,aperture=None,rank=None):
        return next(r for r in rows if r['case']==case and r['method']==method and r['view']==view and
                    (budget is None or r['budget']==budget) and r.get('aperture')==aperture and r.get('rank')==rank)
    def val(*args,metric='relativeL2',**kwargs):
        return lookup(*args,**kwargs).get(metric,float('inf'))
    syn = [c['case'] for c in specs() if c['suite'] in ('1A','1B')]
    fre = [c['case'] for c in specs() if c['suite']=='1C']
    ids_b = [c['case'] for c in specs() if c['suite']=='2A']
    comparisons_a = {}
    for arm in A_ARMS:
        output = {}
        for T in A_BUDGETS:
            for view in ('raw','light'):
                x = np.array([val(c,arm,view,T) for c in syn]); ref = np.array([min(val(c,b,view,T) for b in A_ARMS[:2]) for c in syn])
                rr = [lookup(c,arm,view,T) for c in syn]
                output[f'{int(T)}_{view}'] = dict(median_l2=float(np.median(x)),wins=int(np.sum(x<ref)),
                       median_paired_gain=float(np.median(1-x/ref)),gains_by_case=dict(zip(syn,(1-x/ref).tolist())),
                       median_data_residual=float(np.median([r.get('relative_data_residual',float('inf')) for r in rr])),
                       median_half_iou=float(np.median([r.get('half_contrast_IoU',0.) for r in rr])),
                       median_wall=float(np.median([r['wall_s'] for r in rr])),
                       median_lu=float(np.median([r['LU'] for r in rr])))
        rr,rl = output['30_raw'],output['30_light']
        both = [c for c in syn if rr['gains_by_case'][c]>0 and rl['gains_by_case'][c]>0]
        matches = {v:{c:val(c,arm,v,30.,metric='heldout_relative_complex_error') for c in fre} for v in ('raw','light')}
        frost = {v:{c:min(val(c,b,v,30.,metric='heldout_relative_complex_error') for b in A_ARMS[:2]) for c in fre} for v in ('raw','light')}
        output['both_wins30'] = both
        output['synth_qualifies'] = len(both)>=6 and min(rr['median_paired_gain'],rl['median_paired_gain'])>=.20
        output['fresnel30'] = matches
        # The user's full/sparse requirement applies to the paired Twin fields.
        twin = fre[:2]
        output['fresnel_qualifies'] = all(matches[v][c]<=frost[v][c] for c in twin for v in matches) and all(matches[v][twin[1]]<frost[v][twin[1]] for v in matches)
        comparisons_a[arm] = output
    comparisons_b,diag = {},{}
    for aper in ('full','limited'):
        for k in (12,16):
            key = f'{aper}_{k}'
            d = {}
            for view in ('native','matched'):
                nrmse = {arm:np.array([val(c,arm,view,aperture=aper,rank=k,metric='amplitudeNRMSE') for c in ids_b]) for arm in B_ARMS}
                cand = nrmse['a17_proj_fast']
                best = np.minimum(nrmse['receiver_svd'],nrmse['projected_internal'])
                full = nrmse[B_ARMS[-1]]
                d[view] = dict(median={arm:float(np.median(v)) for arm,v in nrmse.items()},
                               both_wins=int(np.sum(cand<best)), median_lowrank_gain=float(np.median(1-cand/best)),
                               receiver_wins=int(np.sum(cand<nrmse['receiver_svd'])),
                               median_receiver_gain=float(np.median(1-cand/nrmse['receiver_svd'])),
                               full_wins=int(np.sum(full<cand)),median_full_gain=float(np.median(1-full/cand)),
                               paired_nrmse={arm:arr.tolist() for arm,arr in nrmse.items()})
            sums = {arm:sum(lookup(c,arm,'native',aperture=aper,rank=k)['wall_s'] for c in ids_b) for arm in B_ARMS}
            d['wall_totals'] = sums
            d['wall_ratio'] = sums['a17_proj_fast']/sums['receiver_svd']
            n,m = d['native'],d['matched']
            d['clear_full_advantage'] = m['full_wins']>=6 and m['median_full_gain']>=.15
            d['qualifies'] = (aper=='limited' and n['both_wins']>=6 and n['median_lowrank_gain']>=.15 and
                              d['wall_ratio']<=1.5 and not d['clear_full_advantage'])
            comparisons_b[key] = d
            events = []
            for c in ids_b:
                root = f'{aper}_k{k}_a17_proj_fast'
                tr = json.loads(str(bundles[c][root+'_diagnostics_json']))
                matches = [dict(r,case=c) for r in tr if 'exchange_event' in r]
                events.extend(matches)
            verified = [e for e in events if e['verified']]
            accepted = [e for e in verified if e['accepted_swap']]
            def med(items,field):
                vs = [r[field] for r in items if r.get(field) is not None]
                return float(np.median(vs)) if vs else None
            diag[key] = dict(events=len(events),verified=len(verified),positive_q=len(accepted),
                 projection_sign_reversals=sum(e['reduced_unprojected_gain']>0>=e['reduced_projected_gain'] for e in verified),
                 projected_to_exact_sign_reversals=sum(e['reduced_projected_gain']>0>=e['Q_proj'] for e in verified),
                 projection_median=med(verified,'swap_projection_fraction'),
                 base_transfer_error_median=med(verified,'base_transfer_relative_error'),
                 swap_transfer_error_median=med(verified,'swap_transfer_relative_error'),
                 nonlinear_remainder_median=med(accepted,'prediction_remainder'),
                 positive_q_unit_uphill=sum(e.get('unit_trial_phi_change',0)>0 for e in accepted),
                 positive_q_backtracked=sum(e.get('halvings',0) not in (0,None) for e in accepted),
                 exchange_wall_s=sum(e.get('exchange_wall_s',0) for e in events),
                 q_max=max((e['Q_proj'] for e in verified),default=None),detail=events)
    return dict(A=comparisons_a,B=comparisons_b,B_events=diag,
                B_unit_step_diagnostics={c:json.loads(str(bundles[c]['unit_step_diagnostics_json'])) for c in ids_b},
                MAX_STRUCT_IMAGING_VERDICT='GO' if any(comparisons_a[arm]['synth_qualifies'] and comparisons_a[arm]['fresnel_qualifies'] for arm in A_ARMS[2:]) else 'KILL',
                A17_IMAGING_VERDICT='GO' if any(x['qualifies'] for x in comparisons_b.values()) else 'KILL')


def score_all(workers):
    cases = specs()
    missing = [c['case'] for c in cases if not (RUNTIME/'pred'/f'{c["case"]}.npz').is_file()]
    if missing:
        raise RuntimeError('Finish all predictions before offline material or held-out scoring: '+', '.join(missing))
    print(f'All {len(cases)} prediction files exist. Starting offline scoring.',flush=True)
    results = run_jobs(score_case,[(case,) for case in cases],min(workers,4))
    bundles = {case['case']:a for case,rows,a in results}
    ordered = {case['case']:rows for case,rows,a in results}
    rows = [r for spec in cases for r in ordered[spec['case']]]
    verdict = numerical_summary(rows,bundles)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/'metrics.csv').open('w',newline='',encoding='utf-8') as f:
        csv_writer = csv.DictWriter(f,fieldnames=fields)
        csv_writer.writeheader()
        csv_writer.writerows(rows)
    arrays = {case+'__'+name:v for case,bundle in bundles.items() for name,v in bundle.items()}
    arrays['verdict_summary_json'] = np.array(json.dumps(verdict,ensure_ascii=False))
    arrays['keywords'] = np.array('per-case arrays use case__variable; *_raw_audit columns wall,CPU,Phi,train_residual,step,stage; *_trace columns wall,CPU,Phi,train_residual,step')
    np.savez_compressed(OUT/'predictions.npz',**arrays)
    print(f'Scored {len(rows)} retained output rows. {verdict["MAX_STRUCT_IMAGING_VERDICT"]=} {verdict["A17_IMAGING_VERDICT"]=}',flush=True)
    make_figures(rows,bundles)


def make_figures(rows,bundles):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib.colors import PowerNorm, Normalize
    from matplotlib.lines import Line2D

    figures = OUT/'figures'
    figures.mkdir(exist_ok=True)
    labels = dict(zip(A_ARMS,('NIE','Dual-CIE','Block 01 SW','Block 11 SW','DR')))
    labels.update(zip(B_ARMS,('Rx-SVD','Projected G','A17 proj fast','Full gradient')))
    colors = dict(zip(A_ARMS,('#646b76','#15729d','#b45528','#714c93','#2b8571')))
    colors.update(zip(B_ARMS,('#15729d','#714c93','#b45528','#2b8571')))
    index = {(r['case'],r['method'],r['view'],r.get('rank'),r.get('aperture'),r['budget']):r for r in rows}
    bindex = {(r['case'],r['method'],r['view'],r['rank'],r['aperture']):r for r in rows if r['suite']=='2A'}
    specs_all = specs()
    synthetic = [c for c in specs_all if c['suite'] in ('1A','1B')]
    fresnel = [c for c in specs_all if c['suite']=='1C']
    bcases = [c for c in specs_all if c['suite']=='2A']

    def case_label(case):
        s = case['case']
        if case['suite']=='1C':
            return s.replace('Foam','').replace('TM_','\n').replace('_',' ')
        return s.removeprefix('new_').removeprefix('a17_').replace('_',' ').replace('MHz ','MHz\n')

    def image_map(ax,val,side,norm,xfirst=False):
        ax.set_xticks([])
        ax.set_yticks([])
        if val is None:
            ax.set_facecolor('#ebedf0')
            ax.text(.5,.5,'Unavailable',ha='center',va='center',fontsize=7,transform=ax.transAxes)
            return None
        raster = val.real.reshape(side,side)
        return ax.imshow(raster.T if xfirst else raster,origin='lower',
                         interpolation='nearest',norm=norm,cmap='viridis')

    def save(fig,name):
        fig.savefig(figures/(name+'.png'),dpi=160,facecolor='white')
        plt.close(fig)

    def finish_montage(fig,axes,example,scale,title,subtitle):
        fig.suptitle(title+'\n'+subtitle,fontsize=10.5)
        if example is not None:
            fig.colorbar(example,ax=axes.ravel().tolist(),label=scale,fraction=.018,pad=.019)

    for T in A_BUDGETS:
        for view in ('raw','light'):
            path = f'{int(T)}_{view}'
            fig,ax = plt.subplots(len(synthetic),6,figsize=(9.9,15.2),layout='constrained')
            norm = PowerNorm(.65,0,20)
            for i,case in enumerate(synthetic):
                a = bundles[case['case']]
                ex = image_map(ax[i,0],a['truth_coarse'],32,norm)
                ax[i,0].set_ylabel(case_label(case),rotation=0,ha='right',va='center',fontsize=7)
                for j,arm in enumerate(A_ARMS,1):
                    r = index[case['case'],arm,view,None,None,T]
                    image_map(ax[i,j],a.get(r['key']),32,norm)
                    if r['score_status']=='SCORED':
                        ax[i,j].set_xlabel(f"L2 {r['relativeL2']:.2f}  /  {r['relative_data_residual']:.2f}",fontsize=7)
            for x,name in zip(ax[0],['Reference']+[labels[a] for a in A_ARMS]):
                x.set_title(name,fontsize=9)
            finish_montage(fig,ax,ex,'Real permittivity contrast',
                           f'MAX-STRUCT / cross-grid TM / {int(T)} s {view}',
                           'Each prediction shows material relative L2 / training data residual.')
            save(fig,'max_struct_'+path)

            fig,ax = plt.subplots(len(fresnel),5,figsize=(8.9,9.8),layout='constrained')
            norm = PowerNorm(.55,0,20)
            ex = None
            for i,case in enumerate(fresnel):
                a = bundles[case['case']]
                ax[i,0].set_ylabel(case_label(case),rotation=0,ha='right',va='center',fontsize=7)
                for j,arm in enumerate(A_ARMS):
                    r = index[case['case'],arm,view,None,None,T]
                    p = image_map(ax[i,j],a.get(r['key']),60,norm)
                    if p is not None:
                        ex = p
                        ax[i,j].set_xlabel(f"{r['relative_data_residual']:.3f} / {r['heldout_relative_complex_error']:.3f}",fontsize=7)
            for x,arm in zip(ax[0],A_ARMS):
                x.set_title(labels[arm],fontsize=9)
            finish_montage(fig,ax,ex,'Real permittivity contrast',
                           f'Institut Fresnel 2005 / 60 x 60 / {int(T)} s {view}',
                           'Labels are train / held-source complex relative error.')
            save(fig,'fresnel_'+path)

    for aperture in ('full','limited'):
        for rank in (12,16):
            fig,ax = plt.subplots(len(bcases),9,figsize=(13.5,12.8),layout='constrained')
            norm = Normalize(-.5,4)
            for i,case in enumerate(bcases):
                a = bundles[case['case']]
                ex = image_map(ax[i,0],a['truth_coarse'],24,norm,xfirst=True)
                ax[i,0].set_ylabel(case_label(case),rotation=0,ha='right',va='center',fontsize=7)
                for j,(view,arm) in enumerate(((v,a) for v in ('native','matched') for a in B_ARMS),1):
                    r = bindex[case['case'],arm,view,rank,aperture]
                    image_map(ax[i,j],a.get(r['key']),24,norm,xfirst=True)
                    if r['score_status']=='SCORED':
                        ax[i,j].set_xlabel(f"{r['amplitudeNRMSE']:.3f} / {r['wall_s']:.2f}s",fontsize=7)
            ax[0,0].set_title('Reference',fontsize=9)
            for j,(view,arm) in enumerate(((v,a) for v in ('native','matched') for a in B_ARMS),1):
                ax[0,j].set_title(f'{labels[arm]}\n{view}',fontsize=8)
            finish_montage(fig,ax,ex,'Real permittivity contrast',
                           f'A17 / {aperture} aperture / rank {rank} / cross-grid TM',
                           'Labels are complex-material amplitude NRMSE / charged wall time.')
            save(fig,f'a17_{aperture}_k{rank}')

    fig,ax = plt.subplots(2,2,figsize=(11.0,7.2),layout='constrained')
    times = np.arange(1,31)
    for arm in A_ARMS:
        traces = []
        for case in synthetic:
            a = bundles[case['case']]
            key = f'{arm}_raw_audit'
            reference = a['truth_coarse']
            if key not in a:
                traces.append(np.full(len(times),np.nan))
                continue
            tr = a[key]
            l2 = np.linalg.norm(a[key+'_material']-reference[None,:],axis=1)/np.linalg.norm(reference)
            j = np.minimum.accumulate(tr[:,2])
            steps = np.flatnonzero(np.r_[True,j[1:]<j[:-1]])
            ir = np.searchsorted(tr[steps,0],times,side='right')-1
            values = np.full(len(times),np.nan)
            values[ir>=0] = l2[steps[ir[ir>=0]]]
            traces.append(values)
        ax[0,0].plot(times,np.nanmedian(traces,axis=0),color=colors[arm],label=labels[arm])
        for T in A_BUDGETS:
            items = [index[c['case'],arm,'light',None,None,T] for c in synthetic]
            error = np.median([r.get('relativeL2',np.nan) for r in items])
            ax[0,0].scatter(T,error,marker='D',color=colors[arm],s=23,zorder=4)
        for j,case in enumerate(fresnel[:2]):
            r0 = [index[case['case'],arm,'raw',None,None,T] for T in A_BUDGETS]
            r1 = [index[case['case'],arm,'light',None,None,T] for T in A_BUDGETS]
            ax[0,1].plot(A_BUDGETS,[r.get('heldout_relative_complex_error',np.nan) for r in r0],
                         color=colors[arm],ls='-' if j else '--',marker='o',ms=3)
            ax[0,1].scatter(A_BUDGETS,[r.get('heldout_relative_complex_error',np.nan) for r in r1],
                            color=colors[arm],marker='D',s=20)
    ax[0,0].set(title='A / 10 cross-grid cases',xlabel='Paid wall cap (s)',ylabel='Median relative material L2')
    handles,legend_labels = ax[0,0].get_legend_handles_labels()
    handles += [Line2D([0],[0],linestyle='',marker='D',color='#1e2022',markersize=4)]
    ax[0,0].legend(handles,legend_labels+['Light corrector'],ncol=2,fontsize=7)
    ax[0,1].set(title='A / Twin Fresnel (dashed 13 tx, solid 6 tx)',xlabel='Paid wall cap (s)',
                ylabel='Held-source complex relative error')
    for col,rank in enumerate((12,16)):
        aperture = 'limited'
        frac = np.linspace(.2,1.,17)
        for arm in B_ARMS:
            traces = []
            for case in bcases:
                a = bundles[case['case']]
                key = f'{aperture}_k{rank}_{arm}_trace'
                ref = a['truth_coarse']
                horizon = bindex[case['case'],'receiver_svd','matched',rank,aperture]['budget']
                if key not in a:
                    traces.append(np.full(len(frac),np.nan))
                    continue
                tr = a[key]
                error = np.sqrt(np.mean(np.abs(a[key+'_m']-ref[None,:])**2,axis=1))/np.max(np.abs(ref))
                minima = np.minimum.accumulate(tr[:,2])
                steps = np.flatnonzero(np.r_[True,minima[1:]<minima[:-1]])
                ir = np.searchsorted(tr[steps,0],horizon*frac,side='right')-1
                values = np.full(len(frac),np.nan)
                values[ir>=0] = error[steps[ir[ir>=0]]]
                traces.append(values)
            ax[1,col].plot(frac,np.nanmedian(traces,axis=0),color=colors[arm],label=labels[arm])
        ax[1,col].set(title=f'B / limited-view (rank {rank})',xlabel='Fraction of receiver-SVD wall',
                      ylabel='Median complex-material NRMSE')
    ax[1,1].legend(ncol=2,fontsize=7)
    for a in ax.flat:
        a.grid(color='#abb1b8',alpha=.24,lw=.7)
    fig.suptitle('Error versus paid time with physics-selected prefixes',fontweight='bold')
    save(fig,'time_error')


if __name__ == '__main__':
    main()
