"""Shared trace Gram compression and interval certificates for physical rows."""
from dataclasses import dataclass
from time import perf_counter

from common import REGION, herm, star
from bounds import (acb_mat, arb, assemble, ball, block, bt, frob, gersh, low, midpoint,
                    propose_frame, rows, sector, sector_ball, select, shift, strict_pd,
                    tail_proposal, threshold_options)
import numpy as np
from scipy.linalg import qr, solve


class TraceBank:
    def __init__(self, acq, b, mu):
        start = perf_counter()
        n = acq.n
        q = len(REGION)
        self.mu, self.b = mu, b
        self.ell, self.upp = np.zeros(n), np.full(n, b)
        self.upp[REGION] = 0
        self.weight = np.full(n, mu)
        self.weight[REGION] = 0
        self.opt = np.zeros(n)
        self.T = np.vstack((np.eye(n)[REGION], acq.G[REGION]))
        self.Tb = ball(self.T)
        I = np.eye(n)
        self.K = herm(star(acq.S) @ acq.S +
                      sector(I, acq.G, self.ell, self.upp, self.weight, self.opt))
        self.Kb = bt(acq.Sb) * acq.Sb + sector_ball(ball(I), acq.Gb, self.ell, self.upp,
                                                 self.weight, self.opt)
        self.n = n
        self.d = 2*q
        self.verification = {}
        self.available = False
        self.numeric_ready = False
        try:
            _, _, order = qr(self.T, pivoting=True)
            pivot, free = order[:self.d], order[self.d:]
            Tp, Tf = self.T[:, pivot], self.T[:, free]
            self.U = np.zeros((n,self.d),complex)
            self.W = np.zeros((n,n-self.d),complex)
            self.U[pivot] = solve(Tp,np.eye(self.d),check_finite=False)
            self.W[pivot] = -self.U[pivot] @ Tf
            self.W[free] = np.eye(n-self.d)
            self.Cplus = herm(star(self.W) @ self.K @ self.W)
            self.TK = self.U - self.W @ solve(self.Cplus,star(self.W) @ self.K @ self.U,
                                               assume_a="her",check_finite=False)
            self.M = herm(star(self.TK) @ self.K @ self.TK)
            self.numeric_ready = True
            inverse = ball(Tp).solve(ball(np.eye(self.d)))
            dependency = -(inverse*ball(Tf))
            Ub, Wb = acb_mat(n,self.d), acb_mat(n,n-self.d)
            for k, index in enumerate(pivot):
                for j in range(self.d):
                    Ub[int(index),j] = inverse[k,j]
                for j in range(n-self.d):
                    Wb[int(index),j] = dependency[k,j]
            for j, index in enumerate(free):
                Wb[int(index),j] = 1
            self.Wb = Wb
            self.Cplusb = bt(self.Wb)*self.Kb*self.Wb
            kpos, krec = strict_pd(self.Cplusb)
            self.verification.update(Cplus=krec, trace_pivot_pixels=pivot.tolist(),
                                     exact_trace_nullspace=True, trace_rank=self.d,
                                     trace_coordinate_inverse_verified=True)
            if not kpos:
                self.verification["error"] = "public trace-null complement positive gap not established"
                return
            TKb = Ub-self.Wb*self.Cplusb.solve(bt(self.Wb)*self.Kb*Ub)
            self.TKb = TKb
            self.Mb = bt(TKb)*self.Kb*TKb
            self.available = True
        except (ZeroDivisionError, ValueError, np.linalg.LinAlgError) as exc:
            self.verification["error"] = str(exc)
        finally:
            self.preprocess = perf_counter()-start

    def base(self, acq, verified):
        """Apply the cached exact trace-null Schur map to a new observation block."""
        if verified:
            if not self.available:
                raise ArithmeticError(self.verification["error"])
            base = acq.receiver_ball + sector_ball(acq.Jb, acq.Eb, self.ell, self.upp,
                                                    self.weight, self.opt)
            n = self.n
            ab = select(base, range(n), range(n, n+acq.s))
            cc = select(base, range(n, n+acq.s), range(n, n+acq.s))
            null_coupling = bt(self.Wb)*ab
            return block(self.Mb, bt(self.TKb)*ab,
                         cc-bt(null_coupling)*self.Cplusb.solve(null_coupling))
        base = acq.receiver + sector(acq.J, acq.E, self.ell, self.upp, self.weight, self.opt)
        ab = base[:self.n,self.n:]
        cc = base[self.n:,self.n:]
        null_coupling = star(self.W) @ ab
        diag = cc-star(null_coupling) @ solve(self.Cplus,null_coupling,assume_a="her",check_finite=False)
        return herm(block(self.M,star(self.TK) @ ab,diag))


def local_update(acq, ell, upp, mu, local, optical, selected=(), strict=False):
    q = len(REGION)
    w = np.full(q, mu*local)
    o = np.full(q, mu*optical)
    if selected:
        w[list(selected)] = 0.
        o[list(selected)] = 0.
    j = np.hstack((np.eye(q), np.zeros((q, q+acq.s))))
    e = np.hstack((np.zeros((q, q)), np.eye(q), acq.B[REGION]))
    if strict:
        return sector_ball(ball(j), ball(e), ell[REGION], upp[REGION], w, o)
    return sector(j, e, ell[REGION], upp[REGION], w, o)


def chart(acq, selected, material):
    q = len(REGION)
    d = 2*q
    scalar = dict(zip(selected, material))
    free = [i for i in range(d) if i not in scalar]
    matrix = np.zeros((d+acq.s, len(free)+acq.s), dtype=complex)
    matrix[free, np.arange(len(free))] = 1.
    matrix[d:, len(free):] = np.eye(acq.s)
    for i in selected:
        matrix[i, free.index(q+i)] = scalar[i]
        matrix[i, len(free):] = scalar[i] * acq.B[REGION[i]]
    return matrix


def exact_chart(acq, selected, material):
    q = len(REGION)
    d = 2*q
    free = [i for i in range(d) if i not in selected]
    rows_e = {}
    for i in selected:
        e = acb_mat(1, len(free)+acq.s)
        e[0, free.index(q+i)] = 1
        for j in range(acq.s):
            e[0, len(free)+j] = acq.Bb[int(REGION[i]), j]
        rows_e[i] = e
    mapping = ball(chart(acq, selected, [0.]*len(selected)))
    for i, value in zip(selected, material, strict=True):
        for j in range(mapping.ncols()):
            mapping[i, j] = value * rows_e[i][0, j]
    return mapping, rows_e


def verify_chart(acq, trial, exact, selected, ell, upp, terms, radius, save_frame):
    # Direct coordinates make the material pivot exactly the identity for all cells.
    mid = [(arb(float(a))+arb(float(b)))/2 for a, b in zip(ell, upp, strict=True)]
    dev = [(arb(float(b))-arb(float(a)))/2 for a, b in zip(ell, upp, strict=True)]
    P, rows_e = exact_chart(acq, selected, mid)
    pf = chart(acq, selected, [float(v) for v in mid])
    states = P.ncols() - acq.s
    reduced = herm(star(pf) @ trial @ pf)
    curvature = select(exact, selected, selected)
    quadratic_psd, curvature_record = strict_pd(curvature)
    terms_record, distance2 = [], arb(0)
    proposal_failure = None
    for t, j, copies in terms:
        frame, edge = propose_frame(shift(reduced, states, t), states+j)
        if frame is None:
            return None, "frame_proposal_not_positive"
        F = ball(frame)
        M = P * F
        At = shift(exact, 2*len(REGION), t)
        AM = At * M
        H = bt(M) * AM
        d_rows = [rows_e[i]*F for i in selected]
        linear = []
        for k, i in enumerate(selected):
            lpart = bt(d_rows[k]) * rows(AM, [i])
            linear.append(lpart + bt(lpart))
        dim = F.ncols()
        # Second-order is PSD whenever the selected-current principal block is PD.
        # Otherwise bound the entire second-order form by ||A_II||_F ||Delta M||_F².
        quadratic = arb(0)
        if not quadratic_psd:
            quadratic = (frob(curvature) *
                         sum((abs(d)**2 * frob(e)**2 for d, e in zip(dev, d_rows, strict=True)),
                             arb(0))).upper()
        perturbation = []
        for a in range(dim):
            perturbation.append(sum((d * sum((abs(Q[a, b]) for b in range(dim)), arb(0))
                                     for d, Q in zip(dev, linear, strict=True)), arb(0)).upper())
        upper = max(perturbation)
        gap = (gersh(H) - upper - quadratic).lower()
        if not gap > 0:
            proposal_failure = ("uniform_positive_gap_nonpositive"
                                if gersh(H) > 0 else "midpoint_positive_gap_nonpositive")
            return None, proposal_failure
        terms_record.append(dict(frame=save_frame(frame), threshold=t, j=j, term_multiplicity=copies,
                                 positive_count=states+j, floated_positive_edge=edge,
                                 stretched_positive_gap_lower=low(gap), interval_linear_perturbation_upper=float(upper),
                                 interval_quadratic_perturbation_upper=float(quadratic),
                                 frame_rank_proved_by_positive_pullback=True))
        distance2 += arb(t)*copies
    if not distance2 > arb(radius)**2:
        return None, "below_noise"
    return dict(terms=terms_record, states=states, proof_lower=low(distance2.sqrt()),
                strict_distance_over_radius_lower=low(distance2.sqrt()/arb(radius)),
                pivot="identity", pivot_condition_number=1.0,
                principal_block_quadratic_positive=quadratic_psd,
                curvature_check=curvature_record), None


def interval_search(acq, numeric_base, strict_base, ell, upp, selected, params, save_frame,
                    deadline, max_depth=4, max_leaves=16):
    mu, local, optical = params
    trial = numeric_base + local_update(acq, ell, upp, mu, local, optical, selected)
    exact = strict_base + local_update(acq, ell, upp, mu, local, optical, selected, strict=True)
    states = trial.shape[0]-acq.s-len(selected)
    l = ell[REGION[list(selected)]]
    u = upp[REGION[list(selected)]]
    pending = [(l, u, 0, "root")]
    leaves = []
    all_excluded = True
    times = dict(proposal=0., proof=0.)
    while pending:
        if perf_counter() >= deadline:
            all_excluded = False
            leaves.extend(dict(box=[a.tolist(), b.tolist()], depth=depth, path=path,
                               reason="wall_budget") for a, b, depth, path in pending)
            break
        a, b, depth, path = pending.pop()
        started = perf_counter()
        pf = chart(acq, selected, .5*(a+b))
        candidate = tail_proposal(star(pf) @ trial @ pf, states, acq.radius)
        times["proposal"] += perf_counter()-started
        proof, reason = None, "midpoint_tail_at_or_below_noise"
        for terms in threshold_options(candidate, acq.radius):
            started = perf_counter()
            proof, reason = verify_chart(acq, trial, exact, selected, a, b, terms,
                                         acq.radius, save_frame)
            times["proof"] += perf_counter()-started
            if proof is not None:
                break
        if proof is not None:
            leaves.append(dict(box=[a.tolist(), b.tolist()], depth=depth, path=path, verified=proof))
            continue
        if reason != "midpoint_tail_at_or_below_noise" and depth < max_depth and len(leaves)+len(pending)+2 <= max_leaves:
            # Width scaled by the center-dependent magnitude of the shared constitutive row.
            index = int(np.argmax((b-a) / np.maximum(1., .5*(a+b))))
            if a[index] < b[index]:
                middle = .5*(a[index]+b[index])
                uleft, lright = b.copy(), a.copy()
                uleft[index], lright[index] = middle, middle
                pending.append((lright, b, depth+1, path+"R"))
                pending.append((a, uleft, depth+1, path+"L"))
                continue
        all_excluded = False
        leaves.append(dict(box=[a.tolist(), b.tolist()], depth=depth, path=path, reason=reason,
                           strength_numeric=candidate["strength"]))
    return dict(excluded=all_excluded, local_rows=REGION[list(selected)].tolist(),
                trace_rows_within_joint_bank=list(selected), leaves=leaves,
                analytical_cover="closed midpoint bisection covers successive parent intervals",
                searched_weights=list(params), elapsed_components=times)
