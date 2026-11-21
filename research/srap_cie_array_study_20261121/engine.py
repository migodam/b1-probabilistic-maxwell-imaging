"""Sourcewise Dual-CIE plus optional routed alternating field/current projection."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
legacy_path = Path(os.environ.get('SRAP_LEGACY_FILE', HERE / 'legacy_la.py'))
if not legacy_path.is_file():
    legacy_path = HERE.parent / 'limited_aperture_routed_dr_20261120/run_la_study.py'
spec = importlib.util.spec_from_file_location('srap_legacy_la', legacy_path)
legacy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = legacy
spec.loader.exec_module(legacy)

import numpy as np
from scipy import linalg as la
from scipy.fft import fftn
from scipy.optimize import brentq

METHODS = ('dual', 'srap', 'reset', 'full_ap', 'graph', 'srap_cap', 'dual_cap', 'u11')
CYCLED = ('srap', 'reset', 'full_ap', 'graph', 'srap_cap')
ROUTED = ('srap', 'graph', 'srap_cap')


def schedule(case):
    if case.get('suite') == 'LA_FRESNEL':
        return (6, 10, 16), (6., 2., .5), 5
    cfg = case.get('public_schedule')
    if isinstance(cfg, dict):
        return tuple(cfg['modes']), tuple(cfg['betas']), int(cfg['rank'])
    if isinstance(cfg, (list, tuple)) and len(cfg) == 3:
        return tuple(cfg[0]), tuple(cfg[1]), int(cfg[2])
    if cfg in ('legacy_moderate', 'historical_low', 'low', 'legacy_low'):
        return (5, 7, 10), (6., 2., .5), 15
    if cfg not in (None, 'default', 'legacy_high', 'high'):
        raise ValueError(f'Unknown public schedule {cfg!r}')
    return (5, 8, 10), (8., 2., .5), 15


def observed(case_ports, Y, s):
    return Y[:, s] if case_ports.valid_mask is None else Y[case_ports.valid_mask[:, s], s]


def dpcs(model, Y, rank, cap):
    groups = {}
    for j, ids in enumerate(model.ports.indices(model.ns)):
        groups.setdefault(tuple(ids), []).append(j)
    JD = np.zeros_like(model.E)
    bases = []
    for rows, cols in groups.items():
        model.clock.check()
        cols = np.asarray(cols)
        P = model.ports.S[np.asarray(rows)]
        vals, u = la.eigh(P @ P.conj().T, check_finite=False)
        model.clock.counts['receiver_SVD'] += 1
        keep = np.flatnonzero(vals > 1e-12 * max(vals[-1], 1e-30))[-rank:][::-1]
        if cap:
            keep = keep[vals[keep] >= .01 * vals[-1]]
        sigma = np.sqrt(vals[keep])
        V = (P.conj().T @ u[:, keep]) / sigma
        target = np.column_stack([observed(model.ports, Y, j) for j in cols])
        JD[:, cols] = V @ ((u[:, keep].conj().T @ target) / sigma[:, None])
        bases.append((cols, V))
    return JD, bases


def route_basis(model):
    """Legacy weak-receive/strong-propagation eigenspace, without U11 setup."""
    start = model.clock.wall()
    S = model.ports.S
    if model.ports.mask is not None:
        weights = np.zeros(len(S))
        for ids in model.ports.indices(model.ns):
            np.add.at(weights, ids, 1. / model.ns)
        S = np.sqrt(weights[weights > 0, None]) * S[weights > 0]
    # Only singular vectors above .1 of the largest value enter this space.
    vals, u = la.eigh(S @ S.conj().T, check_finite=False)
    keep = vals >= .01 * vals[-1]
    strong = S.conj().T @ (u[:, keep] / np.sqrt(vals[keep]))
    model.clock.counts['receiver_SVD'] += 1

    def gram(v):
        model.clock.check()
        return model.green.adjoint(model.green.apply(v[:, None]))[:, 0]

    scale, _ = legacy.leading_hermitian(model, gram, 1)
    weak = legacy.fast_weak_modes_fft(model, strong, .01 * scale[0])
    model.clock.counts['feedback_route_rank'] = weak.shape[1]
    model.clock.counts['strong_receiver_rank'] = strong.shape[1]
    model.clock.counts['feedback_setup_wall'] = model.clock.wall() - start
    return weak, strong


def column_energy(z):
    return np.einsum('ij,ij->j', z.conj(), z).real


class BallProjector:
    """FFT-PCG solves for H=I+B*B, grouped ports and one global multiplier."""
    scale = 6.0

    def __init__(self, model, Y, delta, recycled=None):
        self.model, self.E = model, model.E
        self.cells, self.sources, self.components = model.n, model.ns, 1
        self.lower, self.upper = np.zeros(model.n), np.full(model.n, 20.)
        self.delta, self.costs = float(delta), model.clock.counts
        self._coarse = None
        self.last_A_info = {}
        self.warm = None
        self.cg_tolerance = 2e-10
        r15 = legacy.R15_ROOT / 'runtime/field_current.py'
        mod = legacy.module_at('srap_field_current', r15 if r15.is_file() else legacy.R12_ROOT / 'src/field_current.py')
        self.material = mod.DenseFieldCurrentProjectors.material.__get__(self)
        self._fields = mod.DenseFieldCurrentProjectors._fields.__get__(self)
        start = model.clock.wall()
        if recycled is not None and recycled.shape[1]:
            # Spectral inverse on the shared strong+weak physics space, identity outside.
            gb = self.B(recycled)
            vals, rot = la.eigh(gb.conj().T @ gb, check_finite=False)
            basis = recycled @ rot
            self._coarse = (basis, (vals / (1. + vals))[:, None] * basis.conj().T)
        grp = {}
        for j, ids in enumerate(model.ports.indices(model.ns)):
            grp.setdefault(tuple(ids), []).append(j)
        self.blocks = []
        active = np.unique(np.concatenate([np.asarray(rows) for rows in grp]))
        T = self.scale * model.ports.S[active]
        # Pivoted QR compresses only receiver RHS, preserving high-precision port maps.
        q, r, piv = la.qr(T.conj().T, mode='economic', pivoting=True, check_finite=False)
        tail2 = np.cumsum(np.sum(abs(r)**2, axis=1)[::-1])[::-1]
        rank = int(np.count_nonzero(tail2 > 1e-24 * tail2[0]))
        transfer = np.empty((rank, len(active)), complex)
        transfer[:, piv] = r[:rank]
        self.costs['port_rhs_compressed_rank'] = rank
        self.costs['port_rhs_rows'] = len(active)
        self.costs['port_compression_relative'] = float(np.sqrt(tail2[rank]/tail2[0])) if rank < len(tail2) else 0.
        z = self.solve(q[:, :rank]) @ transfer
        self.costs['H_setup_rhs'] += rank
        yy = floor2 = 0.
        for rows, cols in grp.items():
            model.clock.check()
            rr = np.searchsorted(active, rows)
            tj, zj = T[rr], z[:, rr]
            K = tj @ zj
            vals, u = la.eigh((K + K.conj().T) * .5, check_finite=False)
            self.costs['port_eigh'] += 1
            tol = 32 * np.finfo(float).eps * max(len(rows), model.n) * np.max(abs(vals))
            if vals[0] < -tol:
                raise ArithmeticError(f'negative port Gram eigenvalue {vals[0]}, tolerance {tol}')
            keep = vals > tol
            target = np.column_stack([observed(model.ports, Y, j) for j in cols])
            floor2 += np.linalg.norm(u[:, ~keep].conj().T @ target)**2
            yy += np.linalg.norm(target)**2
            self.blocks.append((np.array(cols), tj, zj, vals[keep], u[:, keep], target))
        self.floor2, self.yy = float(floor2), float(yy)
        self.bt = 128 * np.finfo(float).eps * max(yy, self.delta**2, 1e-300)
        if floor2 > self.delta**2 + self.bt:
            raise mod.UnreachableDataBall(self.delta, np.sqrt(floor2))
        self.costs['DR_setup_wall'] = model.clock.wall() - start
        self.costs['DR_ball_floor_rel'] = np.sqrt(floor2 / max(yy, 1e-300))

    def B(self, x):
        return self.scale * self.model.green.apply(x)

    def BH(self, q):
        return self.scale * self.model.green.adjoint(q)

    def H(self, x):
        return x + self.BH(self.B(x))

    def precondition(self, r):
        if self._coarse is None:
            return r.copy()
        basis, dual = self._coarse
        return r - basis @ (dual @ r)

    def solve(self, b, start=None):
        # Independent per-column PCG in FFT batches. Check the true H residual.
        x = np.zeros_like(b) if start is None else start.copy()
        if start is not None:
            x[:, ~np.any(b,axis=0)] = 0.
        r = b.copy() if start is None else b - self.H(x)
        z = self.precondition(r)
        p = z.copy()
        rz = np.sum((r.conj()*z).real, axis=0)
        bn = np.maximum(column_energy(b), 1e-300)
        for it in range(801):
            self.model.clock.check()
            active = column_energy(r) > self.cg_tolerance**2 * bn
            if not np.any(active):
                err = np.sqrt(column_energy(b - self.H(x))/bn)
                worst = float(np.max(err, initial=0.))
                self.costs['H_max_true_relative_residual'] = max(worst, self.costs['H_max_true_relative_residual'])
                if worst > 2 * self.cg_tolerance:
                    raise ArithmeticError(f'PCG true residual {worst:g}')
                self.costs['H_CG_batch_iterations'] += it
                return x
            if it == 800:
                break
            hp = self.H(p[:, active])
            alpha = rz[active] / np.sum((p[:, active].conj()*hp).real, axis=0)
            x[:, active] += p[:, active] * alpha
            r[:, active] -= hp * alpha
            za = self.precondition(r[:, active])
            rz_next = np.sum((r[:, active].conj()*za).real, axis=0)
            p[:, active] = za + p[:, active] * (rz_next / rz[active])
            rz[active] = rz_next
            self.costs['H_CG_rhs_iterations'] += int(active.sum())
        raise ArithmeticError('FFT-PCG failed to reach the fixed tolerance')

    def affine_ball(self, X, Q):
        self.costs['A_projections'] += 1
        self.costs['H_iteration_rhs'] += self.sources
        xu = self.solve(X + self.BH(Q - self.E), start=self.warm)
        self.warm = xu.copy()
        parts = []
        initial2 = floor2 = 0.
        for cols, t, v, ev, u, y in self.blocks:
            res = t @ xu[:, cols] - y
            initial2 += np.linalg.norm(res)**2
            z = u.conj().T @ res
            floor2 += np.linalg.norm(res - u @ z)**2
            parts.append((cols, v, ev, u, z, np.sum(abs(z)**2, axis=1)))
        mu = 0.
        if initial2 > self.delta**2:
            if floor2 > self.delta**2 + self.bt:
                raise ArithmeticError('Current projection discrepancy floor exceeds radius')
            if self.delta**2 <= floor2 + self.bt:
                mu = float('inf')
            else:
                def discrepancy(a):
                    self.costs['multiplier_evaluations'] += 1
                    return floor2 + sum(power @ ((1 + a*ev)**-2)
                                        for _, _, ev, _, _, power in parts) - self.delta**2
                hi = 1. / max(ev[-1] for _, _, ev, _, _, _ in parts if len(ev))
                while discrepancy(hi) > 0.:
                    hi *= 4
                    if not np.isfinite(hi):
                        raise ArithmeticError('No finite projection multiplier')
                mu = brentq(discrepancy, 0., hi, xtol=np.nextafter(0.,1.), rtol=8*np.finfo(float).eps)
            for cols, v, ev, u, z, _ in parts:
                xu[:, cols] -= v @ (u @ (z / (ev[:, None] + 1./mu)))
        actual = np.sqrt(sum(np.linalg.norm(t @ xu[:, ids] - y)**2
                             for ids, t, _, _, _, y in self.blocks))
        self.last_A_info = dict(multiplier=mu, affine_data_residual=float(actual),
                                radius=self.delta, unconstrained_residual=float(np.sqrt(initial2)),
                                relative_excess=float(max(0.,actual-self.delta)/max(np.sqrt(self.yy),1e-300)))
        if self.last_A_info['relative_excess'] > 1e-7:
            raise ArithmeticError(f'Projected port discrepancy exceeds radius: {self.last_A_info}')
        return xu, self.E + self.B(xu)


class FourierCurrents(legacy.FourierCurrents):
    def project_coefficients(self, target):
        """Exact small-rank decomposition of U=(I−VV*)F and its pseudoinverse."""
        if self.factors is None:
            self.factors = []
            for _, op in self.groups:
                fv = fftn(op.Vd.reshape(op.shape + (-1,)), axes=(0,1), norm='ortho')[tuple(op.frequencies.T)]
                w, _, _ = la.svd(fv, full_matrices=False, check_finite=False)
                u, s, vh = la.svd(op.apply(w), full_matrices=False, check_finite=False)
                keep = s > 1e-12
                self.factors.append((w, u[:,keep], s[keep], vh[keep]))
                self.cost['Fourier_current_injection_SVD'] += 1
        out = np.empty((self.k, self.ns), complex)
        for (ids, op), (w, u, s, vh) in zip(self.groups, self.factors):
            t = target[:, ids]
            natural = op.adjoint(t)
            out[:, ids] = natural + w @ (vh.conj().T @ ((u.conj().T @ t)/s[:,None]) - w.conj().T @ natural)
        self.cost['Fourier_current_injection_RHS'] += self.ns
        return out


def ratio(a, b):
    return float(np.linalg.norm(a) / max(np.linalg.norm(b), 1e-300))


def run_cie_arm(model, Y, case, method, m, delta):
    from research.research_max_struct_r2.src.cie_global import _baseline_readout
    from research.research_max_struct_r10.src.lowrank_denominator import matrix_material_readout
    if method not in METHODS:
        raise ValueError(method)
    if not hasattr(model, 'trace_route'):
        model.trace_route = []
    trace = model.trace_route
    hook = getattr(model, 'checkpoint_hook', None)
    cost = model.clock.counts
    modes, betas, k = schedule(case)
    start = model.clock.wall()
    JD, groups = dpcs(model, Y, k, method.endswith('_cap'))
    QD = model.E + model.green.apply(JD)
    yn2 = np.array([np.linalg.norm(observed(model.ports, Y, s))**2 for s in range(model.ns)])
    wd = 1. / np.maximum(yn2, 1e-300)
    ws = 1. / np.maximum(np.sum(abs(model.E)**2, axis=0), 1e-300)
    box = np.zeros(model.n), np.full(model.n, 20.)
    reg = lambda z: model.reg.evaluate(z)[:2]
    cfg = dict(cache_material_value_gradient=True)
    axis = None if method == 'u11' else 0
    need_projection = method in (*ROUTED, 'full_ap', 'u11')
    weak = strong = projector = None
    if method == 'u11':
        denominator, weak = legacy.spectral_feedback(model, True)
        cfg.update(matrix_denominator=denominator, matrix_inner_maxiter=8)
        Os, Ls = np.diag(ws), np.diag(np.sqrt(ws))
    elif method in ROUTED:
        weak, strong = route_basis(model)
    trace.append(dict(event='setup', wall=model.clock.wall(), setup_s=model.clock.wall()-start,
                      dpcs_ranks=[v.shape[1] for _, v in groups],
                      weak_rank=0 if weak is None else weak.shape[1],
                      strong_rank=0 if strong is None else strong.shape[1]))
    mx, mq = np.zeros_like(JD), np.zeros_like(QD)
    coeff = {}
    count = 0
    for stage, (mf, bm) in enumerate(zip(modes, betas)):
        model.clock.check()
        bj = betas[0]
        typ = legacy.FourierCurrents if method == 'u11' else FourierCurrents
        U = typ(model.shape, mf, groups, cost)
        frequencies = [tuple(f) for f in U.frequencies]
        A = np.array([coeff.get(f, np.zeros(model.ns, complex)) for f in frequencies])
        fluct = U.apply(A)
        J, Q = JD + fluct, QD + model.green.apply(fluct)
        old_g = old_dir = None
        trace.append(dict(event='stage', stage=stage, step=count, mf=mf, beta_j=bj,
                          beta_m=bm, wall=model.clock.wall()))
        for step in range(500):
            model.clock.check()
            p = (bj / (1. + bj*m))[:,None]
            D = model.ports.receive(J) - Y
            if model.ports.valid_mask is not None:
                D = np.where(model.ports.valid_mask, D, 0.)
            R = p * (J - m[:,None]*Q)
            v = p * (R*ws)
            gradient = U.adjoint(model.ports.adjoint(D*wd) + v - model.green.adjoint(m[:,None]*v))
            direction = -gradient
            if old_g is not None:
                pr = np.maximum(0., np.sum((gradient.conj()*(gradient-old_g)).real,axis=axis) /
                                np.maximum(np.sum(abs(old_g)**2,axis=axis),1e-300))
                direction += old_dir * pr
                bad = np.sum((gradient.conj()*direction).real,axis=axis) >= 0
                if axis is None:
                    if bad:
                        direction = -gradient
                else:
                    direction[:,bad] = -gradient[:,bad]
            dJ = U.apply(direction)
            dQ = model.green.apply(dJ)
            dd = model.ports.receive(dJ)
            if model.ports.valid_mask is not None:
                dd = np.where(model.ports.valid_mask,dd,0.)
            dr = p * (dJ - m[:,None]*dQ)
            den = wd*np.sum(abs(dd)**2,axis=0) + ws*np.sum(abs(dr)**2,axis=0)
            if axis is None:
                den = np.sum(den)
            alpha = -np.sum((gradient.conj()*direction).real,axis=axis)/np.maximum(den,1e-300)
            change = direction*alpha
            delta2 = float(np.mean(np.sum(abs(change)**2,axis=axis) /
                                   np.maximum(np.sum(abs(A)**2,axis=axis),1e-300)))
            A += change
            J += dJ*alpha
            Q += dQ*alpha
            if method == 'u11':
                m, _ = matrix_material_readout(None,J,Q,m,bm,Os,Ls,1.,box,reg,cfg,cost)
            else:
                m, _ = _baseline_readout(J,Q,m,bm,ws,1.,box,reg,True,cfg,cost)
            old_g, old_dir = gradient, direction
            count += 1
            cost['PR_CG_steps'] += 1
            end = step == 499 or (step >= 1 and delta2 < (1e-4 if stage == 0 else 1e-6))
            checkpoint = count % 16 == 0 or end
            cycle = checkpoint and (method in CYCLED or method == 'u11')
            if hook is not None and checkpoint:
                hook('before', count, stage, dict(J=J,Q=Q,m=m,A=A,mx=mx,mq=mq,gradient=gradient,direction=direction))
            if cycle and need_projection:
                start = model.clock.wall()
                if projector is None:
                    if method == 'u11':
                        projector = legacy.make_projector(model,Y,delta)
                    else:
                        recycle = None if weak is None else np.column_stack((strong,weak))
                        projector = BallProjector(model,Y,delta,recycle)
                gx, gq = J/projector.scale, Q
                cx,cq,_ = projector.material(gx+mx,gq+mq)
                ax,aq = projector.affine_ball(cx,cq)
                dx,dq = ax-gx,aq-gq
                if method == 'full_ap':
                    mx,mq = dx,dq
                else:
                    mx = weak @ (weak.conj().T @ dx)
                    mq = weak @ (weak.conj().T @ dq)
                bmx = projector.scale * model.green.apply(mx)
                projected_q = mq
                if method == 'graph':
                    mq = bmx
                sx, _, m = projector.material(gx+mx,gq+mq)
                target = projector.scale*sx-JD
                A = U.project_coefficients(target)
                injected = U.apply(A)
                J = JD+injected
                Q = model.E+model.green.apply(J)
                old_g = old_dir = None
                cost['reflection_cycles'] += 1
                cost['AP_cycles'] += 1
                ports_full = model.ports.receive(projector.scale*dx)
                ports_route = model.ports.receive(projector.scale*mx)
                trace.append(dict(event='AP',stage=stage,step=count,wall=model.clock.wall(),
                                  cycle_s=model.clock.wall()-start,
                                  x_route_fraction=ratio(mx,dx),q_route_fraction=ratio(mq,dq),
                                  commutator_vs_field=ratio(projected_q-bmx,dq),
                                  commutator_vs_routed_field=ratio(projected_q-bmx,projected_q),
                                  graph_error_vs_field=ratio(mq-bmx,dq),
                                  aperture_leak_vs_full=ratio(ports_route,ports_full),
                                  aperture_leak_vs_data=float(np.linalg.norm(ports_route)/max(np.sqrt(yn2.sum()),1e-300)),
                                  reinjection_error=ratio(injected-target,target),
                                  reinjection_change_error=ratio(J-projector.scale*sx,projector.scale*sx-projector.scale*gx),
                                  projector_ball_excess=projector.last_A_info.get('relative_excess'),
                                  setup_s=float(cost.get('DR_setup_wall',0.))))
                if hook is not None:
                    hook('AP',count,stage,dict(dx=dx,dq=dq,component_x=mx,component_q=projected_q,
                         graph_q=bmx,mx=mx,mq=mq,sx=sx,J=J,Q=Q,A=A,m=m))
            elif cycle:
                old_g = old_dir = None
                cost['pr_only_resets'] += 1
            if hook is not None and checkpoint:
                hook('after',count,stage,dict(J=J,Q=Q,m=m,A=A,mx=mx,mq=mq))
            model.clock.check()
            yield m, count, stage, end or cycle
            if end:
                break
        coeff = {f:a.copy() for f,a in zip(frequencies,A)}
