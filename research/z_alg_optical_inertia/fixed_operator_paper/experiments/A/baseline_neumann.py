"""Receiver-direction Neumann bounds with weighted and spectral remainder checks."""
from time import perf_counter

from common import herm, star
from bounds import acb_mat, arb, ball, bt, bdiag, frob, high, low, strict_pd
import numpy as np
from scipy.linalg import solve, svd, svdvals


def abs_matrix(matrix):
    return acb_mat(matrix.nrows(), matrix.ncols(), [abs(a) for a in matrix.entries()])


def weighted_step(matrix, scaling):
    ratios = matrix * scaling
    return max((ratios[i, 0].real / scaling[i, 0].real).upper()
               for i in range(matrix.nrows()))


def total_direction(rows, values):
    return sum(((a * b).real for a, b in zip(rows.entries(), values.entries())), arb(0))


def direction_arrays(W, residual, port_map, e0, h):
    left = bt(W)*port_map
    weights = abs_matrix(bdiag(h, left.transpose()))
    constant = sum(((w.conjugate()*x).real for w, x in zip(W.entries(), residual.entries())),
                   arb(0))
    linear = sum((a*abs(sum(((left[j, i]*e0[i, j]).real for j in range(e0.ncols())), arb(0)))
                  for i, a in enumerate(h)), arb(0))
    return weights, constant, linear


def strict_forward_terms(acq, ell, upp):
    pairs = [(arb(float(a)), arb(float(b))) for a, b in zip(ell, upp)]
    m = [(a+b)/2 for a, b in pairs]
    h = [(b-a)/2 for a, b in pairs]
    A = ball(np.eye(acq.n)) - bdiag(m, acq.Gb)
    J0 = A.solve(bdiag(m, acq.Bb))
    E0 = acq.Bb + acq.Gb*J0
    residual = acq.Yb - acq.Sb*J0
    # Right actions avoid any unverified floating inverse.
    Gmap = A.transpose().solve(acq.Gb.transpose()).transpose()
    Smap = A.transpose().solve(acq.Sb.transpose()).transpose()
    TH = bdiag(h, Gmap.transpose()).transpose()
    return E0, residual, TH, Smap, h


def strict_direction(E0, residual, TH, Smap, h, W, kind, proposal):
    weights, constant, linear = direction_arrays(W, residual, Smap, E0, h)
    M = abs_matrix(TH)
    q0 = abs_matrix(E0)
    if kind == "weighted":
        scaling = ball(proposal[:, None])
        kappa = weighted_step(M, scaling)
        if not kappa < 1:
            return None
        norm_record = {"kind": kind, "kappa_upper": high(kappa),
                       "weights": proposal.tolist(), "proof": "diagonal Collatz weighted row bound"}
    else:
        kappa = arb(proposal)
        pd, certificate = strict_pd(acb_mat(TH.nrows(), TH.ncols(),
                                            [kappa*kappa if i == j else 0
                                             for i in range(TH.nrows()) for j in range(TH.ncols())])
                                     - bt(TH)*TH)
        if not pd:
            return None
        norm_record = {"kind": kind, "kappa_upper": proposal, "proof": certificate}
    qk = q0
    current = linear
    best = None
    for power in range(1, 9):
        qk = M*qk
        if power in (1, 2, 4, 8):
            if kind == "weighted":
                total = arb(0)
                for j in range(qk.ncols()):
                    factor = max((qk[i, j].real / scaling[i, 0].real).upper()
                                 for i in range(qk.nrows())) / (1-kappa)
                    total += factor*sum(((weights[i,j]*scaling[i,0]).real
                                         for i in range(qk.nrows())), arb(0))
            else:
                total = frob(weights)*kappa**power*frob(E0)/(1-kappa)
            bound = (constant-current-total)/frob(W)
            if best is None or bound.lower() > best[0].lower():
                best = (bound, norm_record | {"order": power, "receiver_inner_product_lower": low(constant),
                                             "directional_variation_upper": high(current+total),
                                             "geometric_tail_upper": high(total)})
        current += total_direction(weights, qk)
    return best


def neumann(acq, ell, upp, save_frame, deadline):
    started = perf_counter()
    m, h = (ell+upp)/2, (upp-ell)/2
    n = acq.n
    A = np.eye(n)-m[:, None]*acq.G
    try:
        T = solve(A.T, acq.G.T, check_finite=False).T * h
        J0 = solve(A, m[:, None]*acq.B, check_finite=False)
        residual = acq.Y-acq.S @ J0
        sigma = float(svdvals(T, check_finite=False)[0])
        M = np.abs(T)
        spectral = float(max(abs(np.linalg.eigvals(T))))
        rho_abs = float(max(abs(np.linalg.eigvals(M))))
        alternatives = []
        if sigma < 1.:
            alternatives.append(("spectral", sigma+.01*(1-sigma)))
        if rho_abs < 1.:
            d = solve(np.eye(n)-M, np.ones(n), check_finite=False)
            if np.all(d > 0.):
                d /= max(d)
                alternatives.append(("weighted", d))
        reasons = []
        contraction = dict(stored_midpoint_spectral_norm=sigma,
                           stored_midpoint_weighted_absolute_spectral_radius=rho_abs,
                           positive_endpoint_perturbation_spectral_radius_float=spectral)
        if not alternatives:
            return dict(excluded=False, **contraction, proposal_seconds=perf_counter()-started,
                        strict_seconds=0., reason="neither tested uniform contraction proposed")
        Ud, sd, Vd = svd(residual, full_matrices=False, check_finite=False)
        directions = [residual, sd[0]*np.outer(Ud[:, 0], Vd[0])]
        receiver_basis, _, _ = svd(acq.S, full_matrices=True, check_finite=False)
        weakest = receiver_basis[:, -4:]
        directions.append(weakest @ (star(weakest) @ residual))
        proposal_time = perf_counter()-started
        started = perf_counter()
        exact = strict_forward_terms(acq, ell, upp)
        best = None
        for W in directions:
            if np.linalg.norm(W) == 0. or perf_counter() >= deadline:
                continue
            for kind, proposal in alternatives:
                checked = strict_direction(*exact, ball(W), kind, proposal)
                if checked is None:
                    reasons.append(kind + "_rounding_enclosure_failed")
                    continue
                bound, record = checked
                if best is None or bound.lower() > best[0].lower():
                    best = bound, W, record
        if best is None:
            return dict(excluded=False, **contraction, proposal_seconds=proposal_time,
                        strict_seconds=perf_counter()-started, reasons=reasons)
        bound, W, record = best
        excluded = bound.lower() > arb(acq.radius)
        return dict(excluded=bool(excluded), **contraction, proposal_seconds=proposal_time,
                    strict_seconds=perf_counter()-started, **record, reasons=reasons,
                    strict_distance_over_radius_lower=low(bound/arb(acq.radius)),
                    receiver_direction=save_frame(W) if excluded else None)
    except (ZeroDivisionError, ValueError, np.linalg.LinAlgError) as exc:
        return dict(excluded=False, reason=str(exc), proposal_seconds=perf_counter()-started,
                    strict_seconds=0.)


def replay_neumann(acq, ell, upp, rec, frame):
    parameter = np.array(rec["weights"]) if rec["kind"] == "weighted" else rec["kappa_upper"]
    checked = strict_direction(*strict_forward_terms(acq, ell, upp), ball(frame),
                               rec["kind"], parameter)
    return checked is not None and bool(checked[0].lower() > arb(acq.radius))
