"""Spatial PSD and local/matrix QCQP dual lower bounds within a finite budget."""
from time import perf_counter

from common import REGION, herm, star
from bounds import (acb_mat, arb, ball, bt, high, low, propose_frame, sector,
                    sector_ball, select, validate_frame, weight_vector)
import numpy as np
from scipy.linalg import eigvalsh, solve, hadamard


def build_duals(acq, ell, upp, numeric, deadline):
    """Sample nonnegative spatial diagonal and signed real optical multipliers."""
    trials = {}
    modes = (np.eye(acq.s), hadamard(acq.s, dtype=np.int8))
    for key, border, nu in numeric:
        if nu != 0:
            continue
        w, o = weight_vector(acq.n, *key)
        trial = summarize(acq, border, w, o, modes)
        if trial is not None:
            trials[tuple(w)+tuple(o)] = trial
    # Full native diagonal spatial operators are legal. Under the search cap,
    # optimize coordinates independently in the local region for each class.
    initial = list(trials.values())
    for candidates in (initial, [trial for trial in initial if np.all(trial["optical"] == 0)]):
        if not candidates:
            continue
        best = max(candidates, key=lambda result: result["trace"])
        initial_w = best["weight"]
        initial_o = best["optical"]
        for index in REGION:
            for factor in (.2, 5.):
                if perf_counter() >= deadline:
                    return modes, list(trials.values())
                w = initial_w.copy()
                w[index] *= factor
                key = tuple(w)+tuple(initial_o)
                if key in trials:
                    continue
                border = acq.receiver + sector(acq.J, acq.E, ell, upp, w, initial_o)
                trial = summarize(acq, border, w, initial_o, modes)
                if trial is not None:
                    trials[key] = trial
            if np.any(initial_o != 0):
                o = initial_o.copy()
                o[index] *= -1
                key = tuple(initial_w)+tuple(o)
                if key in trials:
                    continue
                border = acq.receiver+sector(acq.J, acq.E, ell, upp, initial_w, o)
                trial = summarize(acq, border, initial_w, o, modes)
                if trial is not None:
                    trials[key] = trial
    return modes, list(trials.values())


def summarize(acq, border, weight, optical, modes):
    n = acq.n
    K = border[:n, :n]
    vals = eigvalsh(K, check_finite=False)
    if vals[0] <= 1.e-12 * max(1., vals[-1]):
        return None
    ab = border[:n, n:]
    H = herm(border[n:, n:]-star(ab) @ solve(K, ab, assume_a="pos", check_finite=False))
    direct = np.real(np.diag(H))
    rotated = np.real(np.diag(star(modes[1]) @ H @ modes[1]))/acq.s
    return dict(weight=weight, optical=optical, border=border,
                trace=float(np.trace(H).real), mode=(direct, rotated))


def strict_schur(acq, ell, upp, weight, optical, frame):
    if np.any(weight < 0) or frame.shape != (acq.n, acq.n):
        return None, None
    n = acq.n
    exact = acq.receiver_ball + sector_ball(acq.Jb, acq.Eb, ell, upp, weight, optical)
    K = select(exact, range(n), range(n))
    gap, krec = validate_frame(K, frame)
    if not gap > 0:
        return None, krec
    ab = select(exact, range(n), range(n, n+acq.s))
    cc = select(exact, range(n, n+acq.s), range(n, n+acq.s))
    H = cc-bt(ab)*K.solve(ab)
    return H, krec


def potential_plans(modes, trials):
    if not trials:
        return []
    plans = []
    variants = (("spatial_psd", [(i, row) for i, row in enumerate(trials) if np.all(row["optical"] == 0)]),
                ("local_qcqp", list(enumerate(trials))))
    for label, possibilities in variants:
        for i, row in sorted(possibilities, key=lambda item: item[1]["trace"], reverse=True)[:2]:
            plans.append((max(0., row["trace"]), label, [(i, None)]))
    for typ in range(len(modes)):
        assembled = [trial["mode"][typ] for trial in trials]
        best = np.argmax(assembled, axis=0)
        values = np.max(assembled, axis=0)
        group = [(int(i), direction) for direction, (i, v) in enumerate(zip(best, values)) if v > 0.]
        plans.append((float(np.sum(np.maximum(values,0.))), ("source_diagonal_qcqp",
                                                             "matrix_qcqp_hadamard")[typ],
                      group))
    return sorted(plans, reverse=True, key=lambda row: row[0])


def psd_qcqp(acq, ell, upp, numeric, save_frame, deadline):
    started = perf_counter()
    modes, trials = build_duals(acq, ell, upp, numeric, deadline)
    plans = potential_plans(modes, trials)
    proposal = perf_counter()-started
    results, seen = [], {}
    started = perf_counter()
    proof = None
    for value, variant, group in plans:
        if value <= acq.radius**2:
            results.append(dict(variant=variant, numeric_bound_squared=value, reason="at_or_below_noise"))
            continue
        if perf_counter() >= deadline:
            results.append(dict(variant=variant, numeric_bound_squared=value, reason="wall_budget"))
            continue
        try:
            lower = arb(0)
            pieces = []
            typ = int(variant == "matrix_qcqp_hadamard")
            for i, q in group:
                row = trials[i]
                if i not in seen:
                    frame, edge = propose_frame(row["border"][:acq.n, :acq.n], acq.n)
                    if frame is None:
                        raise ArithmeticError("positive K frame not proposed")
                    H, krec = strict_schur(acq, ell, upp, row["weight"], row["optical"], frame)
                    if H is None:
                        raise ArithmeticError("positive K frame failed Arb check")
                    seen[i] = H, krec, frame, edge
                H, krec, frame, edge = seen[i]
                if q is None:
                    term = sum((H[l,l].real for l in range(acq.s)), arb(0))
                    vec = None
                else:
                    vec = ball(modes[typ][:, q:q+1])
                    term = (bt(vec)*H*vec)[0, 0].real/(acq.s if typ else 1)
                lower += max(arb(0), term.lower())
                pieces.append((i, q, term, krec, frame, edge))
            exclusion = lower > arb(acq.radius)**2
            results.append(dict(variant=variant, numeric_bound_squared=value,
                                strict_bound_squared_lower=low(lower), excluded=bool(exclusion)))
            if exclusion and proof is None:
                legs = []
                for i, q, term, krec, frame, edge in pieces:
                    row = trials[i]
                    legs.append(dict(weights=row["weight"].tolist(), optical=row["optical"].tolist(),
                                     source_mode=q, source_basis=("I" if typ == 0 else "Hadamard"),
                                     normalized_term_lower=low(term), frame=save_frame(frame),
                                     positive_inertia=krec, floating_K_gap=edge))
                proof = dict(variant=variant, dual_terms=legs,
                             strict_distance_over_radius_lower=low(lower.sqrt()/arb(acq.radius)))
        except (ZeroDivisionError, ValueError, ArithmeticError) as exc:
            results.append(dict(variant=variant, numeric_bound_squared=value, reason=str(exc)))
    return dict(excluded=proof is not None, proof=proof, attempted=results,
                positive_definite_trials=len(trials), proposal_seconds=proposal,
                strict_seconds=perf_counter()-started)


def replay_qcqp(acq, ell, upp, record, frames):
    lower = arb(0)
    identity, H = np.eye(acq.s), hadamard(acq.s, dtype=np.int8)
    for row in record["proof"]["dual_terms"]:
        w, o = np.array(row["weights"]), np.array(row["optical"])
        schur, _ = strict_schur(acq, ell, upp, w, o, frames[row["frame"]])
        if schur is None:
            return False
        if row["source_mode"] is None:
            term = sum((schur[l, l].real for l in range(acq.s)), arb(0))
        else:
            q = row["source_mode"]
            mat, norm = (H, acq.s) if row["source_basis"] == "Hadamard" else (identity, 1)
            v = ball(mat[:, q:q+1])
            term = (bt(v)*schur*v)[0,0].real/norm
        lower += max(arb(0), term.lower())
    return bool(lower > arb(acq.radius)**2)
