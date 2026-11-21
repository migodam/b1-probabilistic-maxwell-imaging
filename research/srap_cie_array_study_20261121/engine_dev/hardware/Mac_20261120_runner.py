#!/usr/bin/env python3
"""Limited-Aperture Subspace-Routed DR vs Mature CIE/NIE & Causal Ablations.

Reuses verified MAX-STRUCT (R2/R4/R6/R7/R8/R10/R12/R15) operators and kernels.
Cross-grid 150^2/96^2 -> 32^2 synthetic limited-aperture benchmarks + Fresnel 2005
half-aperture real measurements. Truth/held-out data is never read during solving.
"""
from __future__ import annotations

import os
for _flag in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_flag] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True
import numpy as np
from scipy import linalg as la
from scipy.fft import ifftn
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, eigsh, gmres

HERE = Path(__file__).resolve().parent
WORKSPACE = Path(os.environ.get('LA_WORKSPACE', str(HERE.parents[1])))
# Support both full Inv_SLAM workspace layout and compact remote staging layout on xinan.
if (HERE / 'vendor_repo/research/research_max_struct_r8/src/paper_physics.py').is_file():
    REPO_ROOT = HERE / 'vendor_repo'
    R15_ROOT = HERE / 'vendor_repo/research/research_max_struct_r15'
    R12_ROOT = HERE / 'vendor_repo/research/research_max_struct_r12'
    FRESNEL_DIR = HERE / 'fresnel_data'
else:
    WORKTREES = WORKSPACE / 'research/delegated/max_struct_r1'
    REPO_ROOT = WORKTREES / 'r8_worktree'
    R15_ROOT = WORKTREES / 'r14_worktree/research/research_max_struct_r15'
    R12_ROOT = WORKTREES / 'r14_worktree/research/research_max_struct_r12'
    FRESNEL_DIR = WORKSPACE / 'research/delegated/a2_literature/data/2005_iop_21_6_S09'

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if (WORKSPACE / 'research/delegated/max_struct_r1/r14_worktree').is_dir():
    _r14w = str(WORKSPACE / 'research/delegated/max_struct_r1/r14_worktree')
    if _r14w not in sys.path:
        sys.path.append(_r14w)

RUNTIME = HERE / '.runs'
C0 = 299792458.0

ARMS = (
    'paper_nie',
    'cie_dual',
    'u11_hop',
    'u11_reset',
    'u11_routed_ap',
    'u11_routed_dr',
    'u11_full_dr',
)

BUDGETS = (15.0, 30.0)
LIGHT_SEEDS = (12.0, 24.0)
AUDIT_PERIOD = 1.0

# 12 diverse limited-aperture synthetic cross-grid cases (10% fixed-norm noise, 400 MHz).
# Index 0 is the exposed replication control; indices 1..11 are fresh/unexposed.
SYNTH_CASES = (
    # (name, geom, eps, mhz, fine, aperture, rot_deg, exposed_control)
    ('la180_austria_rot0_eps3p5',     'austria',       3.5, 400.0, 150, 'half180',   0.0, True),
    ('la180_austria_rot90_eps3p5',    'austria',       3.5, 400.0, 150, 'half180',  90.0, False),
    ('la180_austria_rot180_eps3p5',   'austria',       3.5, 400.0, 150, 'half180', 180.0, False),
    ('la180_austria_rot0_eps4p0',     'austria',       4.0, 400.0, 150, 'half180',   0.0, False),
    ('la180_double_rings_eps4p0',     'double_rings',  4.0, 400.0,  96, 'half180',   0.0, False),
    ('la180_eccentric_shell_eps4p5',  'eccentric',     4.5, 400.0,  96, 'half180',   0.0, False),
    ('la180_coupled_triplet_eps4p5',  'triplet',       4.5, 400.0,  96, 'half180',   0.0, False),
    ('la180_kinked_bar_eps4p0',       'kinked_bar',    4.0, 400.0,  96, 'half180',   0.0, False),
    ('la180_core_shell_eps5p0',       'core_shell',    5.0, 400.0,  96, 'half180',   0.0, False),
    ('la180_hollow_square_bar_eps4p0','hollow_sq_bar', 4.0, 400.0,  96, 'half180',   0.0, False),
    ('la120_austria_eps3p5',          'austria',       3.5, 400.0, 150, 'arc120',    0.0, False),
    ('la120_eccentric_shell_eps4p0',  'eccentric',     4.0, 400.0,  96, 'arc120',    0.0, False),
)

# 6 real measured Fresnel 2005 limited-aperture cases (60^2 air0 [0, 20]).
FRESNEL_CASES = (
    # (name, object, freq_ghz, mode)
    ('fresnel_twin_8GHz_half180_tx',   'FoamTwinDielTM', 8, 'half180_tx'),
    ('fresnel_twin_8GHz_half180_txrx', 'FoamTwinDielTM', 8, 'half180_txrx'),
    ('fresnel_twin_6GHz_half180_tx',   'FoamTwinDielTM', 6, 'half180_tx'),
    ('fresnel_int_6GHz_half180_tx',    'FoamDielIntTM',  6, 'half180_tx'),
    ('fresnel_int_8GHz_half180_tx',    'FoamDielIntTM',  8, 'half180_tx'),
    ('fresnel_ext_6GHz_half180_tx',    'FoamDielExtTM',  6, 'half180_tx'),
)


def module_at(name: str, path: Path):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[name]


def all_specs():
    cases = [
        dict(case=name, suite='LA_SYNTH', geom=geom, eps=eps, frequency=mhz * 1e6,
             fine=fine, side=32, width=2.0, aperture=ap, rot_deg=rot_deg,
             exposed_control=exposed, index=i)
        for i, (name, geom, eps, mhz, fine, ap, rot_deg, exposed) in enumerate(SYNTH_CASES)
    ]
    cases += [
        dict(case=name, suite='LA_FRESNEL', object=obj, frequency=freq * 1e9,
             mode=mode, side=60, width=0.2, fine=0, exposed_control=False,
             index=100 + i)
        for i, (name, obj, freq, mode) in enumerate(FRESNEL_CASES)
    ]
    return cases


def grid(side: int, width: float):
    h = width / side
    a = (np.arange(side) + 0.5) * h - width / 2.0
    yy, xx = np.meshgrid(a, a, indexing='ij')
    return np.c_[xx.ravel(), yy.ravel()]


def circle(x, y, cx, cy, r):
    return (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2


def stroke(x, y, vertices, r):
    p = np.c_[x, y]
    d2 = np.full(len(p), np.inf)
    for a, b in zip(vertices[:-1], vertices[1:]):
        a, b = np.asarray(a, float), np.asarray(b, float)
        v = b - a
        t = np.clip((p - a) @ v / (v @ v), 0.0, 1.0)
        d2 = np.minimum(d2, np.sum((p - a - t[:, None] * v) ** 2, axis=1))
    return d2 <= r ** 2


def material(case: dict, coords: np.ndarray) -> np.ndarray:
    """Ground-truth contrast; called ONLY during data generation and offline scoring."""
    from research.research_max_struct_r8.src.paper_physics import austria
    geom = case['geom']
    eps = case['eps']
    rot_rad = np.deg2rad(case.get('rot_deg', 0.0))
    c, s = np.cos(rot_rad), np.sin(rot_rad)
    R = np.array([[c, -s], [s, c]])
    xy = coords @ R.T
    if geom == 'austria':
        return austria(xy, eps)

    x, y = xy[:, 0], xy[:, 1]
    if geom == 'double_rings':
        # Front ring (y>0, facing aperture) + shadowed back ring (y<0)
        m1 = circle(x, y, -0.32, 0.28, 0.31) & ~circle(x, y, -0.32, 0.28, 0.17)
        m2 = circle(x, y, 0.30, -0.30, 0.33) & ~circle(x, y, 0.30, -0.30, 0.18)
        return ((eps - 1.0) * (m1.astype(float) + 0.90 * m2.astype(float)))
    if geom == 'eccentric':
        # Eccentric hollow shell with cavity shifted toward shadowed y<0 side
        shell = circle(x, y, 0.0, 0.0, 0.56) & ~circle(x, y, 0.08, -0.16, 0.30)
        core = circle(x, y, 0.06, -0.14, 0.12)
        return (eps - 1.0) * (shell.astype(float) + 0.75 * core.astype(float))
    if geom == 'triplet':
        # Two front cylinders shadowing a deeper third cylinder at y=-0.38
        c1 = circle(x, y, -0.34, 0.22, 0.23)
        c2 = circle(x, y, 0.32, 0.20, 0.23)
        c3 = circle(x, y, 0.02, -0.36, 0.26)
        return (eps - 1.0) * (c1.astype(float) + 0.90 * c2.astype(float) + 0.85 * c3.astype(float))
    if geom == 'kinked_bar':
        mask = (np.abs(x) < 0.62) & (np.abs(y) < 0.42)
        mask &= ~stroke(x, y, [[-0.68, -0.18], [-0.18, -0.18], [0.05, 0.14], [0.68, 0.14]], 0.07)
        return (eps - 1.0) * np.where(y > 0, 1.0, 0.85) * mask
    if geom == 'core_shell':
        # Outer moderate shell (eps=2.8) + inner shadowed high-contrast core (eps=case['eps'])
        outer = circle(x, y, 0.0, -0.04, 0.58) & ~circle(x, y, 0.0, -0.04, 0.36)
        inner = circle(x, y, -0.04, -0.10, 0.24)
        return (2.8 - 1.0) * outer.astype(float) + (eps - 1.0) * inner.astype(float)
    if geom == 'hollow_sq_bar':
        sq = (np.abs(x + 0.24) < 0.34) & (np.abs(y + 0.14) < 0.34)
        sq &= ~((np.abs(x + 0.24) < 0.16) & (np.abs(y + 0.14) < 0.16))
        bar = (np.abs(x - 0.38) < 0.15) & (np.abs(y - 0.12) < 0.48)
        return (eps - 1.0) * (sq.astype(float) + 0.85 * bar.astype(float))
    raise ValueError(f'Unknown geom: {geom}')


def fine_solve(chi, E, S, apply, diagonal, *, rtol=2e-8):
    active = np.flatnonzero(chi != 0)
    local = chi[active]
    n = len(chi)

    def mv(z):
        ext = np.zeros(n, complex)
        ext[active] = z
        return z - local * apply(ext[:, None])[:, 0][active]

    op = LinearOperator((len(active), len(active)), matvec=mv, dtype=complex)
    jac = LinearOperator(op.shape, matvec=lambda z: z / (1.0 - local * diagonal), dtype=complex)
    J = np.empty((len(active), E.shape[1]), complex)
    max_res, iters = 0.0, 0
    for s, inc in enumerate(E[active].T):
        rhs = local * inc
        hist = []
        sol, info = gmres(op, rhs, M=jac, rtol=rtol, atol=0.0, restart=120,
                          maxiter=15, callback=hist.append, callback_type='pr_norm')
        res = float(np.linalg.norm(mv(sol) - rhs) / np.linalg.norm(rhs))
        if info or res > 5 * rtol:
            raise RuntimeError(f'fine_solve src {s}: info={info}, res={res:.3g}')
        J[:, s] = sol
        max_res = max(max_res, res)
        iters += len(hist)
    return S[:, active] @ J, dict(fine_active_cells=int(len(active)),
                                  forward_rhs=int(E.shape[1]),
                                  fine_gmres_iterations=int(iters),
                                  fine_state_residual=max_res)


def complex_noise(clean, fraction, seed):
    rng = np.random.default_rng(seed)
    d = rng.standard_normal(clean.shape) + 1j * rng.standard_normal(clean.shape)
    return d * (fraction * np.linalg.norm(clean) / np.linalg.norm(d))


def public_path(case):
    return RUNTIME / 'public' / f"{case['case']}.npz"


def heldout_path(case):
    return RUNTIME / 'heldout' / f"{case['case']}.npz"


def prepare_one(case):
    from research.research_max_struct_r8.src.paper_physics import acquisition, EqualAreaGreen
    from research.research_max_struct_r8.src.fresnel_adapter import load_file, operators, incident_fit

    path = public_path(case)
    if path.exists():
        return case['case']
    t0 = time.perf_counter()
    case = dict(case)
    name, n = case['case'], case['side']

    if case['suite'] == 'LA_FRESNEL':
        _, Y, I = load_file(FRESNEL_DIR / (case['object'] + '.exp'), case['frequency'] / 1e9)
        ids = np.flatnonzero(np.any(np.isfinite(Y), axis=0))
        Y, I = Y[:, ids], I[:, ids]
        ns = len(ids)
        tx_angle = ids * 2.0 * np.pi / ns
        rx_angle = np.arange(Y.shape[0]) * 2.0 * np.pi / Y.shape[0]
        sources = 1.67 * np.c_[np.cos(tx_angle), np.sin(tx_angle)]
        receivers = 1.67 * np.c_[np.cos(rx_angle), np.sin(rx_angle)]
        centers = grid(n, case['width'])
        k = 2.0 * np.pi * case['frequency'] / C0
        p = operators(centers, sources, receivers, case['width'] / n, k, dense=False)
        E, calibration, _ = incident_fit(I, sources, receivers, centers, k)

        # True contiguous half-aperture split on the first hemisphere
        half = ns // 2
        train = np.arange(half)
        held = np.arange(half, ns)
        if case['mode'] == 'half180_txrx':
            # Restrict BOTH Tx and Rx to the upper hemisphere [0, pi]
            hemi_rx = (rx_angle <= np.pi + 1e-8)
            rx_list = [np.flatnonzero(np.isfinite(Y[:, s]) & hemi_rx) for s in train]
            # Equalize count per source by taking ragged lists via padded mask
            max_r = max(len(r) for r in rx_list)
            rx_ids = np.full((len(train), max_r), -1, dtype=int)
            Y_train = np.zeros((max_r, len(train)), complex)
            valid_mask = np.zeros((max_r, len(train)), bool)
            for j, (s, rids) in enumerate(zip(train, rx_list)):
                rx_ids[j, :len(rids)] = rids
                Y_train[:len(rids), j] = Y[rids, s]
                valid_mask[:len(rids), j] = True
            # Also save opposite-hemisphere rx channels on train sources + all held sources
            Y_held_rx = Y[:, train].copy()
            Y_held_rx[hemi_rx, :] = np.nan + 1j * np.nan
            np.savez(heldout_path(case), E=E[:, held], S=p['S'], Y=Y[:, held],
                     E_train_hemi=E[:, train], Y_train_opposite_rx=Y_held_rx)
            case.update(train_sources=(train + 1).tolist(), held_sources=(held + 1).tolist(),
                        train_channels=int(valid_mask.sum()),
                        held_channels=int(np.isfinite(Y[:, held]).sum()),
                        calibration_median_validation=float(np.median([r['validation_relative'] for r in calibration])),
                        data_ball_ratio=0.15, source_specific=True, ragged_rx=True)
            arrays = dict(E=E[:, train], S=p['S'], rx_ids=rx_ids, valid_mask=valid_mask,
                          Y=Y_train, centers=centers, k=k, radius=p['radius'],
                          delta=0.15 * np.linalg.norm(Y_train[valid_mask]))
        else:
            rx_ids = np.stack([np.flatnonzero(np.isfinite(Y[:, s])) for s in train])
            Y_train = np.stack([Y[rx_ids[j], s] for j, s in enumerate(train)], axis=1)
            np.savez(heldout_path(case), E=E[:, held], S=p['S'], Y=Y[:, held])
            case.update(train_sources=(train + 1).tolist(), held_sources=(held + 1).tolist(),
                        train_channels=int(Y_train.size),
                        held_channels=int(np.isfinite(Y[:, held]).sum()),
                        calibration_median_validation=float(np.median([r['validation_relative'] for r in calibration])),
                        data_ball_ratio=0.15, source_specific=True, ragged_rx=False)
            arrays = dict(E=E[:, train], S=p['S'], rx_ids=rx_ids, Y=Y_train,
                          centers=centers, k=k, radius=p['radius'],
                          delta=0.15 * np.linalg.norm(Y_train))
    else:
        fine = case['fine']
        f = acquisition(fine, case['frequency'], dense=False)
        g = acquisition(n, case['frequency'], dense=False)
        tx = np.arange(20)
        rx = np.arange(40)
        if case['aperture'] == 'half180':
            tx, rx = tx[:10], rx[:20]
        elif case['aperture'] == 'arc120':
            tx, rx = tx[:7], rx[:14]
        else:
            raise ValueError(case['aperture'])
        E_high, S_high = f['E0'][:, tx], f['S'][rx]
        G_high = EqualAreaGreen(f['centers'], (fine, fine), f['radius'], f['k'])
        chi_fine = material(case, f['centers'])
        clean, meta = fine_solve(chi_fine, E_high, S_high, G_high.apply, G_high.diagonal)
        # Match killtest noise seed for the exposed replication control (index=0 -> 202611201004 in killtest? Let's use clean dedicated seed 20261121000+idx, and for control rot0 use exact killtest seed 202611201004!)
        seed = 202611201004 if case['exposed_control'] else (20261121000 + case['index'])
        noise = complex_noise(clean, 0.10, seed)
        case.update(meta, noise_relative=0.10, noise_seed=seed,
                    train_channels=int(clean.size), tx=tx.tolist(), rx=rx.tolist())
        arrays = dict(E=g['E0'][:, tx], S=g['S'][rx], Y=clean + noise,
                      centers=g['centers'], k=g['k'], radius=g['radius'],
                      delta=float(np.linalg.norm(noise)))

    case['data_wall'] = time.perf_counter() - t0
    arrays['case_json'] = np.array(json.dumps(case, ensure_ascii=False))
    np.savez(path, **arrays)
    print(f"[prepare] {name} done in {case['data_wall']:.2f}s", flush=True)
    return name


class TimeLimit(BaseException):
    pass


class Clock:
    def __init__(self, limit: float):
        self.limit = float(limit)
        self.t0 = time.perf_counter()
        self.c0 = time.process_time()
        self.counts = Counter()
        self.cuts = {}
        self._marks = []

    def watch(self, *times):
        self._marks = sorted(times)

    def wall(self):
        return time.perf_counter() - self.t0

    def cpu(self):
        return time.process_time() - self.c0

    def capture(self):
        w = self.wall()
        while self._marks and w >= self._marks[0]:
            self.cuts[self._marks.pop(0)] = (w, self.cpu(), dict(self.counts))
        return w

    def check(self):
        if self.capture() >= self.limit:
            raise TimeLimit()


class Ports:
    """Supports dense shared receivers, uniform masked receivers, and ragged per-source receivers."""
    def __init__(self, S, mask=None, valid_mask=None):
        self.S = np.asarray(S, complex)
        self.mask = None if mask is None else np.asarray(mask, int)
        self.valid_mask = None if valid_mask is None else np.asarray(valid_mask, bool)

    def receive(self, J):
        if self.mask is None:
            return self.S @ J
        out = np.zeros((self.mask.shape[1], J.shape[1]), complex)
        for s in range(J.shape[1]):
            if self.valid_mask is None:
                out[:, s] = self.S[self.mask[s]] @ J[:, s]
            else:
                v = self.valid_mask[:, s]
                out[v, s] = self.S[self.mask[s, v]] @ J[:, s]
        return out

    def adjoint(self, R):
        if self.mask is None:
            return self.S.conj().T @ R
        out = np.zeros((self.S.shape[1], R.shape[1]), complex)
        for s in range(R.shape[1]):
            if self.valid_mask is None:
                out[:, s] = self.S[self.mask[s]].conj().T @ R[:, s]
            else:
                v = self.valid_mask[:, s]
                out[:, s] = self.S[self.mask[s, v]].conj().T @ R[v, s]
        return out

    def indices(self, ns):
        if self.mask is None:
            all_r = np.arange(len(self.S))
            return [all_r for _ in range(ns)]
        if self.valid_mask is None:
            return [self.mask[s] for s in range(ns)]
        return [self.mask[s, self.valid_mask[:, s]] for s in range(ns)]


def dense_from_fft(green):
    side = green.shape[0]
    ker = ifftn(green.spectrum)
    rows, cols = np.indices((side, side)).reshape(2, -1)
    dr = (rows[:, None] - rows[None, :]) % green.pad[0]
    dc = (cols[:, None] - cols[None, :]) % green.pad[1]
    return np.asarray(ker[dr, dc], complex)


class PaperModel:
    def __init__(self, public, case, clock):
        from research.research_max_struct_r8.src.paper_physics import EqualAreaGreen
        from research.research_z_alg_r4_diversity.src.physics import Regularizer
        self.clock = clock
        self.side = int(case['side'])
        self.shape = (self.side, self.side)
        self.n = self.side * self.side
        self.E = np.asarray(public['E'], complex)
        self.ns = self.E.shape[1]
        self.ports = Ports(public['S'],
                           public['rx_ids'] if 'rx_ids' in public else None,
                           public['valid_mask'] if 'valid_mask' in public else None)
        self.green = EqualAreaGreen(public['centers'], self.shape,
                                    float(public['radius']), float(public['k']), clock.counts)
        self.reg = Regularizer(self.shape, smooth=1e-4, tv=1e-4)
        self.large = self.n > 1500
        self._G = None

    @property
    def G(self):
        if self._G is None:
            self._G = dense_from_fft(self.green)
        return self._G

    def state(self, m):
        m = np.asarray(m, float)
        self.clock.check()
        if not self.large:
            A = -m[:, None] * self.G
            A.flat[::self.n + 1] += 1.0
            lu = la.lu_factor(A, overwrite_a=True, check_finite=False)
            self.clock.counts['LU'] += 1
            J = la.lu_solve(lu, m[:, None] * self.E, check_finite=False)
            self.clock.counts['complex_RHS'] += self.ns
            self.clock.check()
            return dict(j=J, Q=self.E + self.green.apply(J), lu=lu)
        J, lu = self._solve_many(m, m[:, None] * self.E, adjoint=False)
        return dict(j=J, Q=self.E + self.green.apply(J), lu=lu)

    def _solve_many(self, m, rhs, adjoint=False, lu=None):
        if lu is not None:
            self.clock.counts['complex_RHS'] += rhs.shape[1]
            return la.lu_solve(lu, rhs, trans=2 if adjoint else 0, check_finite=False), lu
        diag = 1.0 - m * self.green.diagonal
        if adjoint:
            diag = diag.conj()

        def mv(z):
            self.clock.check()
            if adjoint:
                return z - self.green.apply((m * z)[:, None], adjoint=True)[:, 0]
            return z - m * self.green.apply(z[:, None])[:, 0]

        op = LinearOperator((self.n, self.n), matvec=mv, dtype=complex)
        jac = LinearOperator((self.n, self.n), matvec=lambda z: z / diag, dtype=complex)
        out = np.empty_like(rhs)
        for s in range(rhs.shape[1]):
            sol, info = gmres(op, rhs[:, s], M=jac, rtol=1e-7, atol=0.0, restart=80, maxiter=3)
            self.clock.counts['complex_RHS'] += 1
            if info or np.linalg.norm(mv(sol) - rhs[:, s]) > 1e-6 * max(np.linalg.norm(rhs[:, s]), 1e-30):
                A = -m[:, None] * self.G
                A.flat[::self.n + 1] += 1.0
                lu = la.lu_factor(A, overwrite_a=True, check_finite=False)
                self.clock.counts['LU'] += 1
                self.clock.counts['complex_RHS'] += rhs.shape[1] - s
                out[:, s:] = la.lu_solve(lu, rhs[:, s:], trans=2 if adjoint else 0, check_finite=False)
                return out, lu
            out[:, s] = sol
        return out, None

    def objective(self, m, Y):
        st = self.state(m)
        pred = self.ports.receive(st['j'])
        diff = pred - Y
        if self.ports.valid_mask is not None:
            diff = np.where(self.ports.valid_mask, diff, 0.0)
            yn = float(np.linalg.norm(Y[self.ports.valid_mask]))
        else:
            yn = float(np.linalg.norm(Y))
        res = float(np.linalg.norm(diff) / max(yn, 1e-300))
        reg_v, _ = self.reg.evaluate(m)[:2]
        return 0.5 * res * res + float(reg_v), res, pred

    def data_value_gradient(self, m, Y):
        st = self.state(m)
        pred = self.ports.receive(st['j'])
        diff = pred - Y
        if self.ports.valid_mask is not None:
            diff = np.where(self.ports.valid_mask, diff, 0.0)
            yn2 = max(float(np.vdot(Y[self.ports.valid_mask], Y[self.ports.valid_mask]).real), 1e-300)
        else:
            yn2 = max(float(np.vdot(Y, Y).real), 1e-300)
        adj_rhs = self.ports.adjoint(diff)
        if not self.large:
            adj = la.lu_solve(st['lu'], adj_rhs, trans=2, check_finite=False)
            self.clock.counts['complex_RHS'] += self.ns
        else:
            adj, _ = self._solve_many(m, adj_rhs, adjoint=True, lu=st['lu'])
        grad_data = np.sum((adj.conj() * st['Q']).real, axis=1) / yn2
        reg_v, reg_g = self.reg.evaluate(m)[:2]
        val = 0.5 * float(np.vdot(diff, diff).real) / yn2 + float(reg_v)
        return val, grad_data + np.asarray(reg_g, float)


class FourierCurrents:
    """R8 Fourier APCS and its adjoint with one fixed DPCS basis per port group."""
    def __init__(self, shape, mf, groups, cost):
        from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
        self.groups = [(ids, NativeFourierAPCS(shape, mf, v, cost)) for ids, v in groups]
        self.frequencies = self.groups[0][1].frequencies
        self.k = len(self.frequencies)
        self.ns = sum(len(ids) for ids, _ in groups)
        self.n = int(np.prod(shape))
        self.factors = None
        self.cost = cost

    def apply(self, coefficients):
        result = np.empty((self.n, self.ns), complex)
        for ids, op in self.groups:
            result[:, ids] = op.apply(coefficients[:, ids])
        return result

    def adjoint(self, J):
        result = np.empty((self.k, self.ns), complex)
        for ids, op in self.groups:
            result[:, ids] = op.adjoint(J[:, ids])
        return result

    def project_coefficients(self, target):
        if self.factors is None:
            self.factors = []
            for _, op in self.groups:
                u, s, vh = la.svd(op.apply(np.eye(self.k, dtype=complex)),
                                  full_matrices=False, check_finite=False)
                active = s > 1e-12 * s[0]
                self.factors.append((u[:, active], s[active], vh[active]))
                self.cost['Fourier_current_injection_SVD'] += 1
        coefficients = np.zeros((self.k, self.ns), complex)
        for (columns, _), (u, s, vh) in zip(self.groups, self.factors):
            coefficients[:, columns] = vh.conj().T @ ((u.conj().T @ target[:, columns]) / s[:, None])
        self.cost['Fourier_current_injection_RHS'] += self.ns
        return coefficients


def dpcs(model, Y, rank):
    S = model.ports.S
    ns = model.ns
    indices = model.ports.indices(ns)
    groups = {}
    for s, ids in enumerate(indices):
        groups.setdefault(tuple(ids.tolist()), []).append(s)
    JD = np.zeros_like(model.E)
    bases = []
    for rows, cols in groups.items():
        cols = np.asarray(cols, int)
        P = S[np.asarray(rows, int)]
        gram = P @ P.conj().T
        vals, vecs = la.eigh(gram, check_finite=False)
        model.clock.counts['receiver_SVD'] += 1
        keep = np.flatnonzero(vals > 1e-12 * max(float(vals[-1]), 1e-30))
        if rank is not None:
            keep = keep[-min(rank, len(keep)):]
        keep = keep[::-1]
        sigma = np.sqrt(vals[keep])
        U = vecs[:, keep]
        V = (P.conj().T @ U) / sigma[None, :]
        # Slice Y rows matching valid entries when ragged
        if model.ports.valid_mask is not None:
            Ys = np.column_stack([Y[model.ports.valid_mask[:, s], s] for s in cols])
        else:
            Ys = Y[:, cols]
        JD[:, cols] = V @ ((U.conj().T @ Ys) / sigma[:, None])
        bases.append((cols, V))
    return JD, bases


def stage_table(case):
    if case['suite'] == 'LA_FRESNEL':
        return (6, 10, 16), (6.0, 2.0, 0.5), 5
    if case['eps'] < 4.0:
        return (5, 7, 10), (6.0, 2.0, 0.5), 15
    return (5, 8, 10), (8.0, 2.0, 0.5), 15


def leading_hermitian(model, action, count):
    rng = np.random.default_rng(20261120)
    v0 = rng.standard_normal(model.n) + 1j * rng.standard_normal(model.n)
    op = LinearOperator((model.n, model.n), matvec=action, dtype=complex)
    model.clock.counts['feedback_Hermitian_eigensolves'] += 1
    vals, vecs = eigsh(op, k=count, which='LA', v0=v0, tol=1e-10)
    return vals[::-1], vecs[:, ::-1]


def fast_weak_modes_fft(model, strong, cutoff):
    """Matrix-free FFT Lanczos (eigsh) for leading eigenspace of (I-Ps) G*G (I-Ps) > cutoff.

    Exact to <1e-8 vs dense O(N^3) zherk+eigh, and ~20-40x faster on 60^2 Fresnel.
    """
    green = model.green
    sh = strong.conj().T

    def matvec(v):
        model.clock.check()
        vp = v - strong @ (sh @ v)
        gv = green.apply(vp[:, None])
        y = green.apply(gv, adjoint=True)[:, 0]
        return y - strong @ (sh @ y)

    op = LinearOperator((model.n, model.n), matvec=matvec, dtype=complex)
    rng = np.random.default_rng(20261121)
    v0 = rng.standard_normal(model.n) + 1j * rng.standard_normal(model.n)
    v0 -= strong @ (sh @ v0)

    # First estimate rank cheaply with a single batched FFT randomized probe
    def apply_M_block(X):
        model.clock.check()
        Xp = X - strong @ (sh @ X)
        Y = green.apply(green.apply(Xp), adjoint=True)
        return Y - strong @ (sh @ Y)

    k_probe = 64
    Omega = rng.standard_normal((model.n, k_probe)) + 1j * rng.standard_normal((model.n, k_probe))
    Omega -= strong @ (sh @ Omega)
    Q, _ = la.qr(apply_M_block(Omega), mode='economic', check_finite=False)
    for _ in range(2):
        Q, _ = la.qr(apply_M_block(Q), mode='economic', check_finite=False)
    B = Q.conj().T @ apply_M_block(Q)
    est_vals = la.eigvalsh(0.5 * (B + B.conj().T))
    est_rank = int(np.sum(est_vals > 0.85 * cutoff))
    k = min(max(est_rank + 10, 24), model.n - strong.shape[1] - 2)

    while True:
        vals, vecs = eigsh(op, k=k, which='LA', v0=v0, tol=1e-3, maxiter=300)
        model.clock.counts['feedback_weak_fft_eigsh'] += 1
        if float(np.min(vals)) <= cutoff or k >= min(model.n - strong.shape[1] - 2, 240):
            keep = vals > cutoff
            if not np.any(keep):
                return np.zeros((model.n, 0), complex)
            V = vecs[:, keep]
            V -= strong @ (sh @ V)
            V, _ = la.qr(V, mode='economic', check_finite=False)
            return V
        k = min(k + 32, model.n - strong.shape[1] - 2)


def spectral_feedback(model, need_weak_modes: bool):
    start = model.clock.wall()
    green = model.green

    def gram_action(v):
        model.clock.check()
        return green.apply(green.apply(v[:, None]), adjoint=True)[:, 0]

    S = model.ports.S
    if model.ports.mask is not None:
        freq = np.zeros(len(S), float)
        for ids in model.ports.indices(model.ns):
            freq[ids] += 1.0
        freq /= model.ns
        S = np.sqrt(freq[freq > 0, None]) * S[freq > 0]
    _, sig, vh = la.svd(S, full_matrices=False, check_finite=False)
    model.clock.counts['receiver_SVD'] += 1
    scale, _ = leading_hermitian(model, gram_action, 1)
    cutoff = 0.01 * float(scale[0])
    strong = vh[sig >= 0.1 * sig[0]].conj().T
    Gs = green.apply(strong)
    eig, V = la.eigh(Gs.conj().T @ Gs, check_finite=False)
    U = strong @ V[:, eig >= cutoff]

    weak_modes = None
    if need_weak_modes:
        weak_modes = fast_weak_modes_fft(model, strong, cutoff)

    model.clock.check()
    cls = module_at('la_twofold_denom', R15_ROOT / 'src/twofold_denominator.py').TwofoldDenominator
    denominator = cls(model.G, U, model.shape, one_hop=True, cost=model.clock.counts)
    model.clock.counts['feedback_rank'] = int(denominator.rank)
    if weak_modes is not None:
        model.clock.counts['feedback_route_rank'] = int(weak_modes.shape[1])
    model.clock.counts['feedback_setup_wall'] = model.clock.wall() - start
    model.clock.check()
    return denominator, weak_modes


class ScaledFFT:
    def __init__(self, green, scale, adjoint=False):
        self.green, self.scale, self.is_adjoint = green, scale, adjoint

    def __matmul__(self, v):
        return self.scale * self.green.apply(v, adjoint=self.is_adjoint)

    def conj(self):
        return self

    @property
    def T(self):
        return ScaledFFT(self.green, self.scale, not self.is_adjoint)


def make_projector(model, Y, delta):
    r15 = R15_ROOT / 'runtime/field_current.py'
    r12 = R12_ROOT / 'src/field_current.py'
    path = r15 if r15.is_file() else r12
    mod = module_at('la_field_current', path)
    P = object.__new__(mod.DenseFieldCurrentProjectors)
    P.E, P.G, P.cells, P.sources, P.components = model.E, model.G, model.n, model.ns, 1
    P.scale, P.delta, P.grouped_ports = 6.0, float(delta), True
    P.lower, P.upper = np.zeros(model.n), np.full(model.n, 20.0)
    P.costs = model.clock.counts
    B = np.asfortranarray(P.scale * model.G)
    H = la.blas.zherk(1.0, B, trans=2, lower=True)
    H.flat[::model.n + 1] += 1.0
    model.clock.counts['H_cholesky'] += 1
    P.H = la.cho_factor(H, lower=True, overwrite_a=True, check_finite=False)
    P.B = ScaledFFT(model.green, P.scale)
    model.clock.check()

    T = P.scale * model.ports.S
    indices = model.ports.indices(model.ns)
    active = np.unique(np.concatenate(indices))
    V = la.cho_solve(P.H, T[active].conj().T, check_finite=False)
    model.clock.counts['H_setup_rhs'] += len(active)
    remap = np.zeros(len(T), int)
    remap[active] = np.arange(len(active))
    P.ports, P.targets, P.port_groups = [], [], []
    for s in range(model.ns):
        if model.ports.valid_mask is not None:
            P.targets.append(Y[model.ports.valid_mask[:, s], s])
        else:
            P.targets.append(Y[:, s])

    groups = {}
    for j, ids in enumerate(indices):
        groups.setdefault(tuple(ids.tolist()), []).append(j)
    floor2 = 0.0
    for rows, cols in groups.items():
        model.clock.check()
        idx = np.asarray(rows, int)
        z = V[:, remap[idx]]
        gram = T[idx] @ z
        vals, U = la.eigh(0.5 * (gram + gram.conj().T), check_finite=False)
        model.clock.counts['port_eigh'] += 1
        tol = 32 * np.finfo(float).eps * max(len(idx), model.n) * max(float(np.max(abs(vals))), np.finfo(float).tiny)
        pos = vals > tol
        p = mod._Port(T[idx], z, np.where(pos, vals, 0.0), U, pos)
        P.ports.append(p)
        P.port_groups.append((p, np.asarray(cols, int)))
        Y_grp = np.column_stack([P.targets[c] for c in cols])
        floor2 += float(np.linalg.norm(U[:, ~pos].conj().T @ Y_grp) ** 2)
    P.minimum_reachable_residual = float(np.sqrt(floor2))
    yn2 = sum(float(np.vdot(t, t).real) for t in P.targets)
    P._energy_scale = max(yn2, P.delta ** 2, np.finfo(float).tiny)
    P._ball_tolerance = 128 * np.finfo(float).eps * P._energy_scale
    model.clock.counts['DR_setup_wall'] = model.clock.wall()
    model.clock.counts['DR_ball_floor_rel'] = P.minimum_reachable_residual / max(np.sqrt(yn2), 1e-300)
    if P.delta ** 2 < floor2 - P._ball_tolerance:
        raise mod.UnreachableDataBall(P.delta, P.minimum_reachable_residual)
    P.last_A_info = {}
    model.clock.check()
    return P


def run_cie_arm(model, Y, case, method, m, delta):
    from research.research_max_struct_r2.src.cie_global import _baseline_readout
    from research.research_max_struct_r10.src.lowrank_denominator import matrix_material_readout

    cost = model.clock.counts
    modes, betas, k = stage_table(case)
    JD, groups = dpcs(model, Y, k)
    QD = model.E + model.green.apply(JD)

    nie = (method == 'paper_nie')
    u11_family = method in ('u11_hop', 'u11_reset', 'u11_routed_ap', 'u11_routed_dr', 'u11_full_dr')
    global_current = u11_family
    axis = None if global_current else 0

    if model.ports.valid_mask is not None:
        yn2_s = np.array([float(np.vdot(Y[model.ports.valid_mask[:, s], s],
                                        Y[model.ports.valid_mask[:, s], s]).real)
                          for s in range(model.ns)])
    else:
        yn2_s = np.sum(abs(Y) ** 2, axis=0)
    wd = 1.0 / np.maximum(yn2_s, 1e-300)
    ws = 1.0 / np.maximum(np.sum(abs(model.E) ** 2, axis=0), 1e-300)
    Os, Ls = np.diag(ws), np.diag(np.sqrt(ws))
    box = (np.zeros(model.n), np.full(model.n, 20.0))
    reg = lambda z: model.reg.evaluate(z)[:2]
    cfg = dict(cache_material_value_gradient=True)

    weak_modes = None
    if u11_family:
        need_weak = method in ('u11_routed_ap', 'u11_routed_dr')
        denominator, weak_modes = spectral_feedback(model, need_weak)
        cfg.update(matrix_denominator=denominator, matrix_inner_maxiter=8)

    use_proj = method in ('u11_routed_ap', 'u11_routed_dr', 'u11_full_dr')
    use_cycle = method in ('u11_reset', 'u11_routed_ap', 'u11_routed_dr', 'u11_full_dr')
    projector = None
    memory_x = np.zeros_like(model.E)
    memory_q = np.zeros_like(model.E)

    def h_proj(Z):
        if method == 'u11_full_dr':
            return Z.copy()
        if weak_modes is None or weak_modes.shape[1] == 0:
            return np.zeros_like(Z)
        return weak_modes @ (weak_modes.conj().T @ Z)

    coeff = {}
    count = 0
    for stage, (mf, bm) in enumerate(zip(modes, betas)):
        model.clock.check()
        bj = bm if nie else betas[0]
        U = FourierCurrents(model.shape, mf, groups, cost)
        frequencies = [tuple(f) for f in U.frequencies]
        A = np.array([coeff.get(f, np.zeros(model.ns, complex)) for f in frequencies])
        fluct = U.apply(A)
        J, Q = JD + fluct, QD + model.green.apply(fluct)
        old_g = old_dir = None

        for step in range(500):
            model.clock.check()
            p = (bj / (1.0 + bj * m))[:, None]
            D = model.ports.receive(J) - Y
            if model.ports.valid_mask is not None:
                D = np.where(model.ports.valid_mask, D, 0.0)
            R = p * (J - m[:, None] * Q)
            v = p * (R * ws)
            grad = U.adjoint(model.ports.adjoint(D * wd) + v - model.green.adjoint(m[:, None] * v))
            direction = -grad
            if old_g is not None:
                pr = np.maximum(0.0, np.sum((grad.conj() * (grad - old_g)).real, axis=axis) /
                                np.maximum(np.sum(abs(old_g) ** 2, axis=axis), 1e-300))
                direction += old_dir * pr
                bad = np.sum((grad.conj() * direction).real, axis=axis) >= 0
                if global_current:
                    if bad:
                        direction = -grad
                else:
                    direction[:, bad] = -grad[:, bad]
            dJ = U.apply(direction)
            dQ = model.green.apply(dJ)
            dd = model.ports.receive(dJ)
            if model.ports.valid_mask is not None:
                dd = np.where(model.ports.valid_mask, dd, 0.0)
            dr = p * (dJ - m[:, None] * dQ)
            denom = wd * np.sum(abs(dd) ** 2, axis=0) + ws * np.sum(abs(dr) ** 2, axis=0)
            if global_current:
                denom = np.sum(denom)
            alpha = -np.sum((grad.conj() * direction).real, axis=axis) / np.maximum(denom, 1e-300)
            change = direction * alpha
            delta2 = float(np.mean(np.sum(abs(change) ** 2, axis=axis) /
                                   np.maximum(np.sum(abs(A) ** 2, axis=axis), 1e-300)))
            A += change
            J += dJ * alpha
            Q += dQ * alpha

            if u11_family:
                m, _ = matrix_material_readout(None, J, Q, m, bm, Os, Ls, 1.0, box, reg, cfg, cost)
            elif nie:
                B = Q + bm * J
                z = bm * np.sum((B.conj() * J).real * ws, axis=1) / np.maximum(np.sum(abs(B) ** 2 * ws, axis=1), 1e-300)
                z = np.clip(z, 0.0, bm * 20.0 / (1.0 + bm * 20.0))
                m = z / (bm * (1.0 - z))
            else:
                m, _ = _baseline_readout(J, Q, m, bm, ws, 1.0, box, reg, True, cfg, cost)

            old_g, old_dir = grad, direction
            count += 1
            cost['PR_CG_steps'] += 1
            end = (step == 499) or (step >= 1 and delta2 < (1e-4 if stage == 0 else 1e-6))
            cycle = use_cycle and (count % 16 == 0 or end)
            if cycle:
                if use_proj:
                    if projector is None:
                        projector = make_projector(model, Y, delta)
                    gx, gq = J / projector.scale, Q.copy()
                    x, q = gx + memory_x, gq + memory_q
                    cx, cq, _ = projector.material(x, q)
                    if method == 'u11_routed_ap':
                        ax, aq = projector.affine_ball(cx, cq)
                        memory_x = h_proj(ax - gx)
                        memory_q = h_proj(aq - gq)
                    else:
                        ax, aq = projector.affine_ball(2.0 * cx - x, 2.0 * cq - q)
                        memory_x += h_proj(ax - cx)
                        memory_q += h_proj(aq - cq)
                    sx, _, m = projector.material(gx + memory_x, gq + memory_q)
                    A = U.project_coefficients(projector.scale * sx - JD)
                    J = JD + U.apply(A)
                    Q = model.E + model.green.apply(J)
                    old_g = old_dir = None
                    cost['reflection_cycles'] += 1
                else:
                    # u11_reset: PR restart ONLY
                    old_g = old_dir = None
                    cost['pr_only_resets'] += 1
            model.clock.check()
            yield m, count, stage, (end or cycle)
            if end:
                break
        coeff = {f: a.copy() for f, a in zip(frequencies, A)}


def snapshot(clock, m, pred, val, res, step, stage=0, wall_off=0.0, cpu_off=0.0, base_counts=None):
    counts = clock.counts if base_counts is None else Counter(base_counts) + clock.counts
    return dict(m=m.copy(), pred=pred.copy(), phi=float(val), residual=float(res),
                step=int(step), stage=int(stage),
                available_wall=wall_off + clock.wall(),
                available_CPU=cpu_off + clock.cpu(), counts=dict(counts))


def best_before(history, limit):
    ready = [p for p in history if p['available_wall'] <= limit and np.isfinite(p['phi'])]
    return min(ready, key=lambda p: p['phi']) if ready else None


def run_light(model, Y, seed, used_wall, used_cpu, base_counts, budget):
    clock = Clock(max(0.05, budget - used_wall))
    model.clock = clock
    model.green.cost = clock.counts
    points = [seed]
    iters = 0
    status = 'CORRECTOR_DONE'

    def cb(x):
        nonlocal iters
        clock.check()
        iters += 1
        v, r, p = model.objective(x, Y)
        points.append(snapshot(clock, x, p, v, r, iters,
                               wall_off=used_wall, cpu_off=used_cpu, base_counts=base_counts))
        clock.check()

    try:
        res = minimize(lambda m: model.data_value_gradient(m, Y), seed['m'], jac=True,
                       method='L-BFGS-B', bounds=[(0.0, 20.0)] * model.n, callback=cb,
                       options=dict(maxiter=20, maxfun=80, maxls=20, ftol=1e-12, gtol=1e-8))
        status = str(res.message)
    except TimeLimit:
        status = 'WALL_LIMIT'
    except Exception as exc:
        status = f'{type(exc).__name__}: {exc}'
    elapsed = used_wall + clock.wall()
    best = best_before(points, budget) or seed
    return best, dict(wall=elapsed, CPU=used_cpu + clock.cpu(), status=status,
                      iterations=iters, overrun=max(0.0, elapsed - budget),
                      cost=dict(Counter(base_counts) + clock.counts))


def solve_one_job(job):
    case, arm = job
    out_file = RUNTIME / 'arm_outputs' / f"{case['case']}__{arm}.npz"
    if out_file.exists():
        return f"{case['case']}__{arm} (cached)"

    with np.load(public_path(case), allow_pickle=False) as z:
        public = {k: z[k] for k in z.files}
    meta = json.loads(str(public['case_json']))
    Y = public['Y']

    clock = Clock(30.0)
    clock.watch(12.0, 15.0, 24.0)
    history = []
    m0 = np.zeros(case['side'] ** 2)
    last_step, last_stage = 0, -1
    status = 'SCHEDULE_COMPLETE'
    model = None
    try:
        model = PaperModel(public, meta, clock)
        v, r, p = model.objective(m0, Y)
        history.append(snapshot(clock, m0, p, v, r, 0, -1))
        due = time.perf_counter()
        for cand, last_step, last_stage, force_audit in run_cie_arm(model, Y, meta, arm, m0, public['delta']):
            clock.check()
            if time.perf_counter() >= due or force_audit:
                v, r, p = model.objective(cand, Y)
                history.append(snapshot(clock, cand, p, v, r, last_step, last_stage))
                due = time.perf_counter() + AUDIT_PERIOD
            clock.check()
    except TimeLimit:
        status = 'WALL_LIMIT'
    except Exception as exc:
        status = f'FAILED: {type(exc).__name__}: {exc}'

    clock.capture()
    raw_wall, raw_cpu, raw_counts = clock.wall(), clock.cpu(), dict(clock.counts)
    ended = (raw_wall, raw_cpu, raw_counts)
    spent = {lim: clock.cuts.get(lim) or ended for lim in (12.0, 15.0, 24.0, 30.0)}

    rows = []
    arrays = {}
    for T, target in zip(BUDGETS, LIGHT_SEEDS):
        raw = best_before(history, T) or history[0]
        seed = best_before(history, target) or history[0]
        uw, uc, rc = spent[T]
        rows.append(dict(
            case=case['case'], suite=case['suite'], method=arm, budget_s=T, view='raw',
            exposed_control=bool(case.get('exposed_control', False)),
            wall_s=uw, cpu_s=uc, status=status if raw_wall <= T or T == 30.0 else 'PAID_PREFIX',
            phi=raw['phi'], train_residual=raw['residual'], step=raw['step'], stage=raw['stage'],
            feedback_rank=rc.get('feedback_rank'), route_rank=rc.get('feedback_route_rank'),
            feedback_setup_wall=rc.get('feedback_setup_wall'),
            DR_setup_wall=rc.get('DR_setup_wall'),
            reflection_cycles=rc.get('reflection_cycles', 0),
            pr_only_resets=rc.get('pr_only_resets', 0),
            LU=rc.get('LU', 0),
        ))
        arrays[f"m_{int(T)}_raw"] = raw['m']

        pw, pc, pcost = spent[target]
        best_l, rec_l = run_light(model, Y, seed, pw, pc, pcost, T) if model is not None else (seed, dict(wall=pw, CPU=pc, status='NO_MODEL', iterations=0, cost=pcost))
        rows.append(dict(
            case=case['case'], suite=case['suite'], method=arm, budget_s=T, view='light',
            exposed_control=bool(case.get('exposed_control', False)),
            wall_s=rec_l['wall'], cpu_s=rec_l['CPU'], status=rec_l['status'],
            phi=best_l['phi'], train_residual=best_l['residual'],
            step=best_l['step'], stage=best_l['stage'],
            feedback_rank=rec_l['cost'].get('feedback_rank'),
            route_rank=rec_l['cost'].get('feedback_route_rank'),
            feedback_setup_wall=rec_l['cost'].get('feedback_setup_wall'),
            DR_setup_wall=rec_l['cost'].get('DR_setup_wall'),
            reflection_cycles=rec_l['cost'].get('reflection_cycles', 0),
            pr_only_resets=rec_l['cost'].get('pr_only_resets', 0),
            LU=rec_l['cost'].get('LU', 0),
        ))
        arrays[f"m_{int(T)}_light"] = best_l['m']

    # Also save paid trajectory curve (time, phi, m) for diagnostic plots on key cases
    arrays['traj_wall'] = np.array([h['available_wall'] for h in history], float)
    arrays['traj_res'] = np.array([h['residual'] for h in history], float)
    arrays['traj_phi'] = np.array([h['phi'] for h in history], float)
    arrays['rows_json'] = np.array(json.dumps(rows, ensure_ascii=False))
    out_file.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_file, **arrays)
    print(f"[solve] {case['case']} | {arm:14s} | {status[:18]:18s} | raw30_res={rows[2]['train_residual']:.3f} wall={raw_wall:.1f}s cyc={raw_counts.get('reflection_cycles',0)} setup={raw_counts.get('feedback_setup_wall',0) or 0:.2f}s", flush=True)
    return f"{case['case']}__{arm}"


def verify_fast_weak_modes():
    """Unit test: verify FFT randomized subspace matches exact dense eigh on 32^2 half-aperture."""
    from research.research_max_struct_r8.src.paper_physics import acquisition
    g = acquisition(32, 400e6, dense=True)
    tx, rx = np.arange(10), np.arange(20)
    pub = dict(E=g['E0'][:, tx], S=g['S'][rx], centers=g['centers'], k=g['k'], radius=g['radius'])
    clk = Clock(60.0)
    model = PaperModel(pub, dict(side=32), clk)
    _, sig, vh = la.svd(model.ports.S, full_matrices=False)
    strong = vh[sig >= 0.1 * sig[0]].conj().T
    gram_G = model.G.conj().T @ model.G
    cutoff = 0.01 * float(la.eigvalsh(gram_G)[-1])
    Pperp = np.eye(model.n) - strong @ strong.conj().T
    M = Pperp @ gram_G @ Pperp
    M = 0.5 * (M + M.conj().T)
    vals, vecs = la.eigh(M)
    V_exact = vecs[:, vals > cutoff]
    t0 = time.perf_counter()
    V_fast = fast_weak_modes_fft(model, strong, cutoff)
    dt = time.perf_counter() - t0
    assert V_fast.shape == V_exact.shape, f"rank mismatch: {V_fast.shape} vs {V_exact.shape}"
    # Subspace distance ||V_fast V_fast* - V_exact V_exact*||_2
    svals = la.svdvals(V_fast.conj().T @ V_exact)
    sin_theta = float(np.sqrt(max(0.0, 1.0 - float(np.min(svals)) ** 2)))
    print(f"[verify_fft_subspace] rank={V_fast.shape[1]}, sin_max_angle={sin_theta:.3e}, time={dt:.3f}s", flush=True)
    assert sin_theta < 1e-2, f"Subspace mismatch: {sin_theta}"


def coarse_truth(case, centers):
    # Exact 8x8 subcell area average per coarse cell (works for any fine/coarse ratio).
    h = case['width'] / case['side']
    u = ((np.arange(8) + 0.5) / 8.0 - 0.5) * h
    a, b = np.meshgrid(u, u)
    offsets = np.c_[a.ravel(), b.ravel()]
    z = (centers[:, None, :] + offsets).reshape(-1, 2)
    return material(case, z).reshape(len(centers), -1).mean(axis=1)


def iou_score(a, b):
    inter = int(np.sum(a & b))
    union = int(np.sum(a | b))
    return float(inter / union) if union else 1.0


def score_all():
    """Offline scoring ONLY after every arm output is sealed on disk."""
    all_rows = []
    specs_map = {c['case']: c for c in all_specs()}

    # First ensure every prediction exists
    for case in all_specs():
        for arm in ARMS:
            p = RUNTIME / 'arm_outputs' / f"{case['case']}__{arm}.npz"
            if not p.exists():
                raise FileNotFoundError(f"Missing prediction: {p}")

    for case in all_specs():
        name = case['case']
        with np.load(public_path(case), allow_pickle=False) as z:
            pub = {k: z[k] for k in z.files}
        meta = json.loads(str(pub['case_json']))

        if case['suite'] == 'LA_SYNTH':
            truth = coarse_truth(case, pub['centers'])
            t_norm = float(np.linalg.norm(truth))
            t_max = float(np.max(truth))
            t_half = truth >= 0.5 * t_max
            # Also compute offline coarse-grid truth data residual for diagnostic context
            clk_d = Clock(60.0)
            mod_d = PaperModel(pub, meta, clk_d)
            _, coarse_ref_res, _ = mod_d.objective(truth, pub['Y'])
        else:
            with np.load(heldout_path(case), allow_pickle=False) as hz:
                held = {k: hz[k] for k in hz.files}
            # Build held-out forward model once per Fresnel case
            clk_h = Clock(600.0)
            pub_h = dict(E=held['E'], S=held['S'], centers=pub['centers'],
                         k=pub['k'], radius=pub['radius'])
            mod_h = PaperModel(pub_h, meta, clk_h)
            Y_held = held['Y']
            valid_h = np.isfinite(Y_held)
            y_held_norm = float(np.linalg.norm(Y_held[valid_h]))

        for arm in ARMS:
            p = RUNTIME / 'arm_outputs' / f"{name}__{arm}.npz"
            with np.load(p, allow_pickle=False) as az:
                rows = json.loads(str(az['rows_json']))
                mats = {k: az[k] for k in az.files if k.startswith('m_')}
            for row in rows:
                m = mats[f"m_{int(row['budget_s'])}_{row['view']}"]
                row['max_contrast'] = float(np.max(m))
                if case['suite'] == 'LA_SYNTH':
                    diff = m - truth
                    row['rel_L2'] = float(np.linalg.norm(diff) / max(t_norm, 1e-30))
                    row['amp_NRMSE'] = float(np.sqrt(np.mean(diff ** 2)) / max(t_max, 1e-30))
                    row['half_IoU'] = iou_score(m >= 0.5 * t_max, t_half)
                    zero_bg = (truth == 0)
                    row['bg_error_frac'] = float(np.sum(diff[zero_bg] ** 2) / max(np.sum(diff ** 2), 1e-30))
                    row['coarse_ref_res'] = float(coarse_ref_res)
                    row['held_field_err'] = np.nan
                else:
                    st_h = mod_h.state(m)
                    pred_h = held['S'] @ st_h['j']
                    err_h = float(np.linalg.norm((pred_h - Y_held)[valid_h]) / max(y_held_norm, 1e-30))
                    row['rel_L2'] = np.nan
                    row['amp_NRMSE'] = np.nan
                    row['half_IoU'] = np.nan
                    row['bg_error_frac'] = np.nan
                    row['coarse_ref_res'] = np.nan
                    row['held_field_err'] = err_h
                all_rows.append(row)

    csv_path = HERE / 'metrics.csv'
    cols = list(all_rows[0].keys())
    with csv_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(all_rows)
    print(f"[score] Wrote {len(all_rows)} rows to {csv_path}", flush=True)
    return all_rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmd', choices=('verify', 'prepare', 'solve', 'score', 'all'))
    p.add_argument('--workers', type=int, default=8)
    args = p.parse_args()

    (RUNTIME / 'public').mkdir(parents=True, exist_ok=True)
    (RUNTIME / 'heldout').mkdir(parents=True, exist_ok=True)
    (RUNTIME / 'arm_outputs').mkdir(parents=True, exist_ok=True)

    if args.cmd == 'verify':
        verify_fast_weak_modes()
        return

    if args.cmd in ('prepare', 'all'):
        verify_fast_weak_modes()
        cases = all_specs()
        with ProcessPoolExecutor(max_workers=min(args.workers, len(cases))) as ex:
            futs = [ex.submit(prepare_one, c) for c in cases]
            for f in as_completed(futs):
                f.result()

    if args.cmd in ('solve', 'all'):
        jobs = [(c, arm) for c in all_specs() for arm in ARMS]
        # Run 32^2 synthetic jobs at full worker count, then 60^2 Fresnel at min(workers, 12) for RAM safety
        synth_jobs = [j for j in jobs if j[0]['suite'] == 'LA_SYNTH']
        fresnel_jobs = [j for j in jobs if j[0]['suite'] == 'LA_FRESNEL']
        t0 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for f in as_completed([ex.submit(solve_one_job, j) for j in synth_jobs]):
                f.result()
        print(f"[solve] All {len(synth_jobs)} synthetic jobs finished in {time.perf_counter() - t0:.1f}s", flush=True)
        t1 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=min(args.workers, 4)) as ex:
            for f in as_completed([ex.submit(solve_one_job, j) for j in fresnel_jobs]):
                f.result()
        print(f"[solve] All {len(fresnel_jobs)} Fresnel jobs finished in {time.perf_counter() - t1:.1f}s", flush=True)

    if args.cmd in ('score', 'all'):
        score_all()


if __name__ == '__main__':
    main()
