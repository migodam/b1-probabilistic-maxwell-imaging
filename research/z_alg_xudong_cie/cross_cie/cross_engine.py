#!/usr/bin/env python3
"""Fourier-current optical corrections to the mature sourcewise Dual-CIE arm."""

from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

SRAP = Path(__file__).resolve().parents[2] / "srap_cie_array_study_20261121"
sys.path.insert(0, str(SRAP))
import engine as mature

METHODS = ("dual", "reset", "diag", "full_ot", "full", "regional4", "regional16", "pixel")


def hermitian(matrix):
    return (matrix + matrix.conj().T)/2


def optical(matrix):
    return (matrix - matrix.conj().T)/(2j)


def positive_part(matrix):
    value, basis = eigh(matrix, check_finite=False)
    return (basis * np.maximum(value, 0.)) @ basis.conj().T


def opt_test(matrix, diagonal, passive):
    if diagonal:
        entries = np.diag(matrix).real
        return np.diag(np.maximum(entries, 0.) if passive else entries)
    return positive_part(matrix) if passive else matrix


def norm_squared(value):
    return float(np.vdot(value, value).real)


def blocks(shape, partitions):
    indexes = np.arange(np.prod(shape)).reshape(shape)
    return [stratum.ravel() for band in np.array_split(indexes, partitions, axis=0)
            for stratum in np.array_split(band, partitions, axis=1)]


class FrozenLaw:
    """Dimensionless local forms with source weights and checkpoint denominators."""

    def __init__(self, J, E, B, cap, subdomains, diagonal=False, passive=False, sector_weight=0.):
        self.cap = np.asarray(cap).reshape(-1, 1)
        self.weight = 1./np.sqrt(np.maximum(np.sum(abs(B)**2, axis=0), 1.e-300))
        self.diagonal = diagonal
        self.passive = passive
        self.sector_weight = sector_weight
        jw, ew = J*self.weight, E*self.weight
        self.regions = []
        for pix in subdomains:
            optical_scale = max(norm_squared(ew[pix]), 1.e-30)
            sector_scale = max((np.linalg.norm(jw[pix]) + np.linalg.norm(self.cap[pix]*ew[pix]))**2, 1.e-30)
            self.regions.append((pix, optical_scale, sector_scale))

    def gradient_parts(self, J, E):
        """Real Euclidean gradient, separated into direct and Green-adjoint terms."""
        jw, ew = J*self.weight, E*self.weight
        direct = np.zeros_like(J)
        adjoint_rhs = np.zeros_like(J)
        opt_loss = sector_loss = 0.
        for pix, os, ss in self.regions:
            j, e = jw[pix], ew[pix]
            q = opt_test(optical(j.conj().T @ e), self.diagonal, self.passive)
            opt_loss += .5*norm_squared(q/os)
            w = q/(os*os)
            direct[pix] -= 1j*(e @ w)*self.weight
            adjoint_rhs[pix] += 1j*(j @ w)*self.weight
            if self.sector_weight:
                u = self.cap[pix]
                raw_sector = hermitian(j.conj().T @ j) - hermitian(j.conj().T @ (u*e))
                sector = positive_part(raw_sector)
                sector_loss += .5*norm_squared(sector/ss)
                w = self.sector_weight*sector/(ss*ss)
                direct[pix] += ((2*j-u*e) @ w)*self.weight
                adjoint_rhs[pix] -= (u*j @ w)*self.weight
        return direct, adjoint_rhs, opt_loss, sector_loss

    def line(self, J, E, dJ, dE):
        """Quadratic coefficients before projection or norms on the real line."""
        jw, ew, djw, dew = [z*self.weight for z in (J, E, dJ, dE)]
        terms = []
        for pix, os, ss in self.regions:
            j, e, dj, de = (z[pix] for z in (jw, ew, djw, dew))
            jh, djh = j.conj().T, dj.conj().T
            o = [optical(Z)/os for Z in (jh@e, jh@de+djh@e, djh@de)]
            sector = None
            if self.sector_weight:
                u = self.cap[pix]
                sector = [hermitian(Z)/ss for Z in (jh @ (j-u*e),
                          jh @ (dj-u*de)+djh @ (j-u*e), djh @ (dj-u*de))]
            terms.append((o, sector))
        return terms

    def line_value(self, terms, alpha):
        opt_loss = sector_loss = 0.
        for o, sector in terms:
            q = opt_test(o[0]+alpha*(o[1]+alpha*o[2]), self.diagonal, self.passive)
            opt_loss += .5*norm_squared(q)
            if sector is not None:
                r = positive_part(sector[0]+alpha*(sector[1]+alpha*sector[2]))
                sector_loss += .5*norm_squared(r)
        return opt_loss, sector_loss


class PixelLaw:
    """Exact per-pixel loss with a fixed dark-row regularization energy."""

    def __init__(self, B, E, passive):
        self.weight = 1./np.sqrt(np.maximum(np.sum(abs(B)**2, axis=0), 1.e-300))
        self.delta = max(1.e-12*norm_squared(E*self.weight)/len(E), 1.e-30)
        self.passive = passive

    def gradient_parts(self, J, E):
        j, e = J*self.weight, E*self.weight
        a = np.sum(abs(j)**2, axis=1, keepdims=True)
        b = np.sum(abs(e)**2, axis=1, keepdims=True)
        c = np.sum(j.conj()*e, axis=1, keepdims=True)
        safe = b+self.delta
        if self.passive:
            radius = np.sqrt(np.maximum(a*b-c.real*c.real, 0.))
            positive = np.maximum(.5*(c.imag+radius), 0.)
            loss = positive*positive
            slope = positive/safe
            ratio = np.divide(slope,radius,out=np.zeros_like(radius),where=radius>0.)
            direct = -1j*slope*e+ratio*(b*j-c.real*e)
            ge = 1j*slope*j+ratio*(a*e-c.real*j)-2.*loss*e/(safe*safe)
        else:
            power = np.maximum(a*b-(c*c).real, 0.)
            loss = .5*power
            direct = (b*j-c*e)/safe
            ge = (a*e-c.conj()*j)/safe-power*e/(safe*safe)
        return direct*self.weight, ge*self.weight, float(np.sum(loss/safe)), 0.

    def line(self, J, E, dJ, dE):
        j, e, dj, de = [z*self.weight for z in (J, E, dJ, dE)]
        return (
            (np.sum(abs(j)**2,axis=1),2*np.sum((j.conj()*dj).real,axis=1),np.sum(abs(dj)**2,axis=1)),
            (np.sum(abs(e)**2,axis=1),2*np.sum((e.conj()*de).real,axis=1),np.sum(abs(de)**2,axis=1)),
            (np.sum(j.conj()*e,axis=1),np.sum(j.conj()*de+dj.conj()*e,axis=1),np.sum(dj.conj()*de,axis=1)),
        )

    def line_value(self, terms, alpha):
        a,b,c = [s0+alpha*(s1+alpha*s2) for s0,s1,s2 in terms]
        if self.passive:
            radius = np.sqrt(np.maximum(a*b-c.real*c.real, 0.))
            loss = np.maximum(.5*(c.imag+radius), 0.)**2
        else:
            loss = .5*np.maximum(a*b-(c*c).real, 0.)
        return float(np.sum(loss/(np.maximum(b, 0.)+self.delta))), 0.


def masked(model, values):
    mask = model.ports.valid_mask
    return np.where(mask, values, 0.) if mask is not None else values


def correct_current(model, basis, A, J, E, m, Y, cap, p, wd, ws, subdomains, method, passive,
                    strength, sector_strength, trust, count, stage):
    """Take one safeguarded projected law step, keeping the shared material fixed."""
    cost = model.clock.counts
    cost["cross_checkpoint_trials"] += 1
    t0 = model.clock.wall()
    model.clock.check()
    law = (PixelLaw(model.E, E, passive) if method=="pixel" else
           FrozenLaw(J, E, model.E, cap, subdomains, method=="diag", passive, sector_strength))
    direct, adj_rhs, opt0, sec0 = law.gradient_parts(J, E)
    grad = direct + model.green.adjoint(adj_rhs)
    dA = -basis.adjoint(grad)
    projected = norm_squared(dA)
    info = dict(step=count, stage=stage, accepted=False, opt_before=opt0, sector_before=sec0,
                coefficient_gradient_Fro=float(np.sqrt(projected)), line_evaluations=0,
                spatial_regions=len(J) if method=="pixel" else len(subdomains))
    if projected <= np.finfo(float).tiny:
        return False, dict(info, cause="zero_projected_gradient", elapsed=model.clock.wall()-t0)
    dJ = basis.apply(dA)
    dE = model.green.apply(dJ)
    D = masked(model, model.ports.receive(J)-Y)
    dD = masked(model, model.ports.receive(dJ))
    C = J-m[:, None]*E
    dC = dJ-m[:, None]*dE
    source_weight = np.sqrt(ws)
    r, dr = p*C*source_weight, p*dC*source_weight
    a = norm_squared(r)/max(norm_squared(C*source_weight), 1.e-30)
    wD, wdD = D*np.sqrt(wd), dD*np.sqrt(wd)
    base0 = .5*(norm_squared(r)+norm_squared(wD))
    base_derivative = float((np.vdot(wD,wdD)+np.vdot(r,dr)).real)
    law_derivative = float(np.vdot(grad, dJ).real)
    derivative = base_derivative + strength*a*law_derivative
    info.update(cie_state_scale=a, base_before=base0, base_derivative=base_derivative,
                law_derivative=law_derivative, merit_derivative=derivative)
    if derivative >= 0:
        return False, dict(info, cause="no_merit_descent", elapsed=model.clock.wall()-t0)
    base_curvature = norm_squared(wdD)+norm_squared(dr)
    alpha = min(-derivative/max(base_curvature, 1.e-30),
                trust*max(np.linalg.norm(J), 1.e-30)/max(np.linalg.norm(dJ), 1.e-30))
    terms = law.line(J, E, dJ, dE)
    law0 = opt0+sector_strength*sec0
    merit0 = base0+strength*a*law0
    accepted = False
    for i in range(6):
        model.clock.check()
        opt1, sec1 = law.line_value(terms, alpha)
        base1 = base0 + alpha*(base_derivative+.5*alpha*base_curvature)
        merit1 = base1+strength*a*(opt1+sector_strength*sec1)
        info["line_evaluations"] += 1
        cost["cross_merit_evaluations"] += 1
        if merit1 < merit0 + .01*alpha*derivative and merit1 < merit0 - 1.e-15*max(1., merit0):
            accepted = True
            break
        alpha *= .5
    if not accepted:
        return False, dict(info, cause="no_admissible_line_step", elapsed=model.clock.wall()-t0)
    coefficient_sizes = np.sum(abs(A)**2, axis=0)
    coefficient_move = float(np.mean(np.sum(abs(alpha*dA)**2,axis=0)/np.maximum(coefficient_sizes,1.e-300)))
    A += alpha*dA
    J += alpha*dJ
    E += alpha*dE
    cost["cross_corrections_accepted"] += 1
    return True, dict(info, accepted=True, alpha=float(alpha), optical_after=opt1, sector_after=sec1,
                      base_after=base1, merit_before=merit0, merit_after=merit1,
                      relative_current_step=float(alpha*np.linalg.norm(dJ)/max(np.linalg.norm(J), 1.e-30)),
                      relative_coefficient_step_sq=coefficient_move,
                      elapsed=model.clock.wall()-t0)


def iterate(model, Y, case, method, m0, upper=20., *, strength=16., sector_strength=1.,
            passive=False, trust=.08, cadence=16, stationary_stages=False,
            unconditional_pr_reset=False):
    """Yield the mature tuple; optional projected corrections share the original wall timer."""
    from research.research_max_struct_r2.src.cie_global import _baseline_readout
    if method not in METHODS:
        raise ValueError(f"Unknown Dual-CIE arm {method!r}.")
    m = np.array(m0, float, copy=True).reshape(model.n)
    u = np.broadcast_to(np.asarray(upper, float), (model.n,)).copy()
    if not np.all(np.isfinite(u)) or np.any(u<0.) or not np.all(np.isfinite(m)):
        raise ValueError("The initial material and the upper box must be finite and real.")
    if strength < 0 or sector_strength < 0 or not 0 < trust <= 1 or cadence < 1:
        raise ValueError("Invalid optical step configuration.")
    np.clip(m,0.,u,out=m)
    box = np.zeros(model.n), u
    modes, betas, rank = mature.schedule(case)
    model.trace_route = []
    model.trace_cross_ot = []
    hook = getattr(model, "checkpoint_hook", None)
    cost = model.clock.counts
    started = model.clock.wall()
    JD, groups = mature.dpcs(model,Y,rank,False)
    QD = model.E+model.green.apply(JD)
    yn2 = np.array([np.linalg.norm(mature.observed(model.ports,Y,s))**2 for s in range(model.ns)])
    wd = 1./np.maximum(yn2,1.e-300)
    ws = 1./np.maximum(np.sum(abs(model.E)**2,axis=0),1.e-300)
    reg = lambda z: model.reg.evaluate(z)[:2]
    cfg = dict(cache_material_value_gradient=True)
    model.trace_route.append(dict(event="setup", wall=model.clock.wall(),
                                  setup_s=model.clock.wall()-started,
                                  dpcs_ranks=[v.shape[1] for _,v in groups], weak_rank=0,strong_rank=0))
    subdomains = blocks(model.shape, {"regional4":2,"regional16":4}.get(method,1))
    sector = sector_strength if method in ("full","regional4","regional16") else 0.
    coeff = {}
    count = 0
    bj = betas[0]
    for stage,(mf,bm) in enumerate(zip(modes,betas)):
        model.clock.check()
        U = mature.FourierCurrents(model.shape,mf,groups,cost)
        frequencies = [tuple(f) for f in U.frequencies]
        A = np.array([coeff.get(f,np.zeros(model.ns,complex)) for f in frequencies])
        fluct = U.apply(A)
        J,Q = JD+fluct,QD+model.green.apply(fluct)
        old_g = old_dir = None
        model.trace_route.append(dict(event="stage",stage=stage,step=count,mf=mf,beta_j=bj,beta_m=bm,
                                      wall=model.clock.wall()))
        for step in range(500):
            model.clock.check()
            p = (bj/(1.+bj*m))[:,None]
            D = masked(model,model.ports.receive(J)-Y)
            R = p*(J-m[:,None]*Q)
            v = p*(R*ws)
            gradient = U.adjoint(model.ports.adjoint(D*wd)+v-model.green.adjoint(m[:,None]*v))
            direction = -gradient
            if old_g is not None:
                pr = np.maximum(0.,np.sum((gradient.conj()*(gradient-old_g)).real,axis=0) /
                                 np.maximum(np.sum(abs(old_g)**2,axis=0),1.e-300))
                direction += old_dir*pr
                bad = np.sum((gradient.conj()*direction).real,axis=0)>=0
                direction[:,bad] = -gradient[:,bad]
            dJ = U.apply(direction)
            dQ = model.green.apply(dJ)
            dd = masked(model,model.ports.receive(dJ))
            dr = p*(dJ-m[:,None]*dQ)
            den = wd*np.sum(abs(dd)**2,axis=0)+ws*np.sum(abs(dr)**2,axis=0)
            alpha = -np.sum((gradient.conj()*direction).real,axis=0)/np.maximum(den,1.e-300)
            change = direction*alpha
            delta2 = float(np.mean(np.sum(abs(change)**2,axis=0) /
                                   np.maximum(np.sum(abs(A)**2,axis=0),1.e-300)))
            A += change
            J += dJ*alpha
            Q += dQ*alpha
            m,_ = _baseline_readout(J,Q,m,bm,ws,1.,box,reg,True,cfg,cost)
            old_g,old_dir = gradient,direction
            count += 1
            cost["PR_CG_steps"] += 1
            tolerance = 1.e-4 if stage==0 else 1.e-6
            end = step==499 or (step>=1 and delta2 < tolerance)
            checkpoint = count%cadence==0 or end
            if hook is not None and checkpoint:
                hook("before",count,stage,dict(J=J,Q=Q,m=m,A=A,gradient=gradient,direction=direction,
                                              mx=np.zeros_like(J),mq=np.zeros_like(Q)))
            accepted = False
            early_end_trial = stationary_stages and end and step<499
            if method not in ("dual","reset") and strength>0 and (count%cadence==0 or early_end_trial):
                p_new = (bj/(1.+bj*m))[:,None]
                accepted,detail = correct_current(model,U,A,J,Q,m,Y,u,p_new,wd,ws,subdomains,
                                                  method,passive,strength,sector,trust,count,stage)
                model.trace_cross_ot.append(detail)
                if accepted:
                    if early_end_trial and detail["relative_coefficient_step_sq"]>=tolerance:
                        end = False
                        cost["cross_stage_end_deferred"] += 1
            if accepted or (count%cadence==0 and
                            (method=="reset" or (method!="dual" and unconditional_pr_reset))):
                old_g = old_dir = None
                cost["cross_PR_resets"] += 1
            if hook is not None and checkpoint:
                hook("after",count,stage,dict(J=J,Q=Q,m=m,A=A,mx=np.zeros_like(J),mq=np.zeros_like(Q)))
            model.clock.check()
            yield m,count,stage,end or accepted
            if end:
                break
        coeff = {f:a.copy() for f,a in zip(frequencies,A)}
