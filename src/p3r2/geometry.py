"""GT-blind, nonlinear candidate trajectories; these are not posterior samples.

Importing this module starts no solver or job. Exact decoder constraints retain
task-scale material progress after a strong-space corrector. Only the caller's
Target performs physical work; every attempted endpoint is charged there.
"""
from pathlib import Path
import json
import time

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize

from .core import save, write

ARMS = ("frozen_fixed", "dynamic_fixed", "dynamic_curvature", "dynamic_bandit", "ordinary_random")
FEATURES = ("bias", "radius", "step_fraction", "direction_fraction", "sign",
    "log_generalized_cost", "predicted_Phi_change", "prior_change",
    "latent_step", "saturation_fraction", "seam_fraction",
    "previous_Maxwell_error", "previous_projector_rotation", "material_curve_error")


def material_frame(J, D, weak_dimension=4):
    """Weak directions per unit material RMS, including original prior cost.

    M=[Re D;Im D]/sqrt(n); T whitens its identifiable Gaussian64 tangent.
    Eigensystem T' (J'J+I) T ranks the original likelihood/prior quadratic
    cost per unit material displacement. No measure determinant is inserted.
    """
    J, D = np.asarray(J, float)[:, :64], np.asarray(D)[:, :64]
    M = np.r_[D.real, D.imag]/np.sqrt(len(D))
    _, s, Vh = np.linalg.svd(M, full_matrices=False)
    keep = s > max(float(s[0])*1e-8, 1e-14)
    if np.count_nonzero(keep) < weak_dimension:
        raise ValueError("insufficient identifiable material tangent dimension")
    T = Vh[keep].T/s[keep]
    H = (J@T).T@(J@T)+T.T@T
    lam, U = eigh((H+H.T)/2)
    directions = T@U[:, :weak_dimension]
    # Reproducible signs; comparisons use projectors rather than basis signs.
    for k in range(weak_dimension):
        i = int(np.argmax(abs(directions[:, k])))
        directions[:, k] *= 1 if directions[i, k] >= 0 else -1
    Q = np.linalg.qr(directions, mode="complete")[0]
    return dict(directions=directions, Vw=Q[:, :weak_dimension],
        Vs=Q[:, weak_dimension:], projector=Q[:, :weak_dimension]@Q[:, :weak_dimension].T,
        eigenvalues=np.maximum(lam, 0), material_singular_values=s,
        material_rank=int(np.count_nonzero(keep)), weak_dimension=weak_dimension)


def _factor(point):
    if point.T is None:
        raise ValueError("receiver transfer must be paid/cached before material action")
    return np.einsum("qp,ps->qsp", point.T, point.state.field).reshape(256, 576)


def _white_material(target, C, dchi):
    return target.white((C@dchi).reshape(32, 8))


def _strata(prior, q):
    _, aux = prior.decode(q)
    f = aux["field"][prior.mask]/prior.temperature
    return np.where(f <= -1, -1, np.where(f >= 1, 1, 0)), f


def _progress(prior, q, origin, velocity):
    chi, _ = prior.decode(q)
    d = chi-origin.state.chi
    return float(np.linalg.norm(d)/np.sqrt(len(d))), float(np.vdot(velocity, d).real/len(d))


def _valid(prior, q, origin, velocity, radius):
    if not prior.legal(q) or not np.array_equal(q[64:], origin.q[64:]):
        return False
    rms, progress = _progress(prior, q, origin, velocity)
    return bool(radius-1e-10 <= rms <= 1.02*radius and progress >= .5*radius-1e-9
        and np.linalg.norm(q[:64]-origin.q[:64]) <= 2+1e-9)


def radius_seed(prior, origin, vector, radius):
    v = np.r_[np.asarray(vector, float), 0., 0.]
    velocity = prior.derivative(origin.q)@v
    scale = np.linalg.norm(velocity)/np.sqrt(len(velocity))
    if scale <= 1e-14:
        return None
    v, velocity = v/scale, velocity/scale
    # Frozen Gaussian displacement cap is two prior SD in aggregate per step.
    cap = 2/max(np.linalg.norm(v[:64]), 1e-30)
    hi = min(radius, cap)
    def fun(t):
        return prior.material_rms(prior.decode(origin.q+t*v)[0], origin.state.chi)-radius
    for _ in range(18):
        if fun(hi) >= 0:
            root = brentq(fun, 0, hi, xtol=1e-13)
            # Tiny outward bias prevents numerical radius undershoot.
            q = origin.q+min(cap, root*(1+1e-9))*v
            return (q, velocity, v) if _valid(prior, q, origin, velocity, radius) else None
        if hi >= cap:
            break
        hi = min(cap, 2*hi)
    return None


def strong_correct(target, origin, predictor, model, frame, velocity, radius, y,
                   maxiter=35):
    """Optimize only Vs coordinates while retaining actual RMS and progress."""
    P, Vs = target.prior, frame["Vs"]
    C = _factor(model)
    residual = target.white(model.state.data-y)
    n = len(origin.state.chi)
    scale = max(1., .5*float(residual@residual))
    def q_from(b):
        q = predictor.q.copy()
        q[:64] += Vs@b
        return q
    cache = {}
    def terms(b):
        if "b" not in cache or not np.array_equal(cache["b"], b):
            q = q_from(b)
            chi, _ = P.decode(q)
            D = P.derivative(q)[:, :64]@Vs
            r = residual+_white_material(target, C, chi-model.state.chi)
            A = np.column_stack([_white_material(target, C, D[:, k]) for k in range(D.shape[1])])
            cache.update(b=b.copy(), q=q, chi=chi, D=D, r=r, A=A)
        return cache
    def objective(b):
        z = terms(b)
        return float((.5*z["r"]@z["r"]+.5*z["q"][:64]@z["q"][:64])/scale)
    def gradient(b):
        z = terms(b)
        return (z["A"].T@z["r"]+Vs.T@z["q"][:64])/scale
    def sphere(b):
        d = terms(b)["chi"]-origin.state.chi
        return float(np.vdot(d, d).real/(n*radius*radius)-1)
    def jsphere(b):
        z = terms(b)
        return 2*np.real(z["D"].conj().T@(z["chi"]-origin.state.chi))/(n*radius*radius)
    def progress(b):
        return float(np.vdot(velocity, terms(b)["chi"]-origin.state.chi).real/(n*radius)-.5)
    def jprogress(b):
        return np.real(terms(b)["D"].conj().T@velocity)/(n*radius)
    def trust(b):
        d = terms(b)["q"][:64]-origin.q[:64]
        return float(1-d@d/4)
    def jtrust(b):
        return -.5*Vs.T@(terms(b)["q"][:64]-origin.q[:64])
    fit = minimize(objective, np.zeros(Vs.shape[1]), jac=gradient,
        method="SLSQP", constraints=[dict(type="eq", fun=sphere, jac=jsphere),
            dict(type="ineq", fun=progress, jac=jprogress), dict(type="ineq", fun=trust, jac=jtrust)],
        options=dict(maxiter=int(maxiter), ftol=1e-9))
    q = q_from(fit.x)
    valid = _valid(P, q, origin, velocity, radius)
    return q, dict(optimizer_success=bool(fit.success), iterations=int(fit.nit),
        message=str(fit.message), radius_and_progress_retained=valid,
        weak_coordinate_change=float(np.linalg.norm(frame["Vw"].T@(q[:64]-predictor.q[:64]))),
        actual_RMS=_progress(P, q, origin, velocity)[0],
        signed_progress=_progress(P, q, origin, velocity)[1])


def feature_vector(target, origin, seed, frame, action, radius, step, steps,
                   C, y, previous_error=0., previous_rotation=0.):
    q, velocity, _ = seed
    chi, _ = target.prior.decode(q)
    dchi = chi-origin.state.chi
    linear = _white_material(target, C, dchi)
    residual = target.white(origin.state.data-y)
    strata, f = _strata(target.prior, origin.q)
    tangent = target.prior.derivative(origin.q)@(q-origin.q)
    k, sign = action//2, (1 if action % 2 == 0 else -1)
    values = np.array([1., radius, step/max(1, steps-1), k/max(1, frame["weak_dimension"]-1), sign,
        np.log1p(frame["eigenvalues"][k]),
        (.5*np.dot(residual+linear, residual+linear)-.5*residual@residual)/max(1., .5*residual@residual),
        .5*(q[:64]@q[:64]-origin.q[:64]@origin.q[:64]),
        np.linalg.norm(q[:64]-origin.q[:64]), np.mean(strata != 0), np.mean(abs(abs(f)-1) < .1),
        previous_error, previous_rotation, np.linalg.norm(dchi-tangent)/max(np.linalg.norm(dchi), 1e-30)], float)
    if not np.isfinite(values).all():
        raise ValueError("nonfinite online feature")
    return values


def load_policy(path):
    obj = json.loads(Path(path).read_text())
    if obj.get("feature_names") != list(FEATURES) or obj.get("kind") != "ridge_contextual_bandit":
        raise ValueError("policy schema mismatch")
    if obj.get("training_split") != "train" or not obj.get("frozen_before_evaluation"):
        raise ValueError("policy not frozen train-only")
    # This exported inference object contains no truth arrays or rewards.
    allowed = {"kind", "feature_names", "mean", "scale", "coefficients",
        "ridge", "training_split", "frozen_before_evaluation", "fit_rows", "fit_parents",
        "reward_definition", "no_posterior_claim", "schema_version"}
    if set(obj)-allowed:
        raise ValueError("unexpected inference policy fields")
    return obj


def policy_score(policy, features):
    x = (np.asarray(features)-np.asarray(policy["mean"]))/np.asarray(policy["scale"])
    return x@np.asarray(policy["coefficients"])


def _persist(target, p, directory):
    # Full Target state/LU plus receiver transfer, and a separate explicit J.
    J = target.jacobian(p)
    path = Path(directory)/"cache"/f"{p.key}.npz"
    target.persist_cache(p, path)
    save(Path(directory)/"receiver_J"/f"{p.key}.npz", q=p.q, J=J,
        material_derivative=target.prior.derivative(p.q), receiver_T=p.T)
    return str(path)


def run_trajectory(target, parent, anchor, *, arm, directory, steps=8,
                   radii=(.02, .04), weak_dimension=4, policy=None, seed=2026102203,
                   max_corrector_iterations=35):
    if arm not in ARMS or not 1 <= int(steps) <= 32:
        raise ValueError("unsupported trajectory arm or step count")
    if any(r not in (.02, .04) for r in radii):
        raise ValueError("frozen task radii are .02/.04")
    directory = Path(directory)
    if (directory/"TRACE.jsonl").exists():
        raise FileExistsError("trajectory trace already exists")
    directory.mkdir(parents=True, exist_ok=True)
    y = target.observation(parent)
    _persist(target, anchor, directory)
    frozen = material_frame(target.jacobian(anchor), target.prior.derivative(anchor.q), weak_dimension)
    save(directory/"BASE_FRAME.npz", **frozen)
    current = anchor
    previous_error = previous_rotation = 0.
    rows = []
    initial_strata, _ = _strata(target.prior, anchor.q)
    start = time.perf_counter()
    initial_counts = dict(target.book.counts)
    for step in range(int(steps)):
        radius = float(radii[step % len(radii)])
        dynamic = material_frame(target.jacobian(current), target.prior.derivative(current.q), weak_dimension)
        frame = frozen if arm == "frozen_fixed" else dynamic
        if arm == "ordinary_random":
            rng = np.random.default_rng(np.random.SeedSequence([int(seed), step]))
            directions = rng.normal(size=(64, weak_dimension))
            Q = np.linalg.qr(directions, mode="complete")[0]
            frame = dict(dynamic, directions=directions, Vw=Q[:, :weak_dimension],
                Vs=Q[:, weak_dimension:], projector=Q[:, :weak_dimension]@Q[:, :weak_dimension].T)
        C = _factor(anchor if arm == "frozen_fixed" else current)
        candidates = []
        for action in range(2*weak_dimension):
            direction = frame["directions"][:, action//2]*(1 if action % 2 == 0 else -1)
            proposal = radius_seed(target.prior, current, direction, radius)
            if proposal is not None:
                features = feature_vector(target, current, proposal, frame, action, radius, step, steps,
                    C, y, previous_error, previous_rotation)
                candidates.append(dict(action=action, seed=proposal, features=features))
        if not candidates:
            rows.append(dict(step=step, status="NO_FEASIBLE_TASK_SCALE_DIRECTION", target_RMS=radius))
            break
        if arm == "dynamic_bandit" and policy is not None:
            scores = [float(policy_score(policy, c["features"])) for c in candidates]
            chosen = int(np.argmax(scores)); rule = "frozen_train_contextual_reward_model"
        elif arm == "dynamic_curvature":
            scores = [float(-c["features"][6]-.05*c["features"][7]-c["features"][-1]) for c in candidates]
            chosen = int(np.argmax(scores)); rule = "linear_data_prior_and_decoder_curvature"
        else:
            # Development behavior coverage is deterministic and independent of GT.
            chosen = (int(seed)+step) % len(candidates) if arm == "dynamic_bandit" else step % len(candidates)
            scores = None; rule = "deterministic_exploration" if arm == "dynamic_bandit" else "fixed_direction_cycle"
        selected = candidates[chosen]
        q, velocity, _ = selected["seed"]
        before = dict(target.book.counts)
        # Persist the decision before the first physical action.
        decision = dict(parent=parent, arm=arm, step=step, target_RMS=radius,
            action=int(selected["action"]), features=selected["features"], feature_names=FEATURES,
            candidates=[dict(action=c["action"], features=c["features"]) for c in candidates],
            decision_scores=scores, decision_rule=rule, GT_used_online=False)
        write(directory/f"step_{step:03d}_DECISION.json", decision)
        predictor = target.evaluate(q, purpose=f"{parent}_{arm}_step{step}_predictor")
        predictor_cache = _persist(target, predictor, directory)
        model = anchor if arm == "frozen_fixed" else predictor
        cq, correction = strong_correct(target, current, predictor, model, frame, velocity, radius, y,
            maxiter=max_corrector_iterations)
        endpoints = [(predictor, "predictor", predictor_cache)]
        if correction["radius_and_progress_retained"] and np.linalg.norm(cq-predictor.q) > 1e-9:
            corrected = target.evaluate(cq, purpose=f"{parent}_{arm}_step{step}_strong_corrector")
            endpoints.append((corrected, "strong_corrector", _persist(target, corrected, directory)))
        accepted, role, cache_path = min(endpoints, key=lambda item: target.phi(item[0], y)+.5*item[0].q[:64]@item[0].q[:64])
        delta = accepted.state.chi-current.state.chi
        actual_response = target.white(accepted.state.data-current.state.data)
        model_response = _white_material(target, C, delta)
        previous_error = float(np.linalg.norm(actual_response-model_response)/max(1., np.linalg.norm(actual_response)))
        new_frame = material_frame(target.jacobian(accepted), target.prior.derivative(accepted.q), weak_dimension)
        previous_rotation = float(np.linalg.norm(new_frame["projector"]-dynamic["projector"], ord=2))
        strata, f = _strata(target.prior, accepted.q)
        full_cost = {k:float(target.book.counts.get(k, 0)-before.get(k, 0)) for k in ("full_forward_calls", "linear_solve_rhs")}
        measured = []
        for point, endpoint_role, endpoint_cache in endpoints:
            measured.append(dict(endpoint=point.key, role=endpoint_role, cache_path=endpoint_cache,
                Phi=target.phi(point, y), objective=target.phi(point, y)+.5*float(point.q[:64]@point.q[:64]),
                actual_RMS=target.prior.material_rms(point.state.chi, current.state.chi),
                signed_progress=_progress(target.prior, point.q, current, velocity)[1],
                backward_residual=float(point.state.residual)))
        row = dict(**decision, status="MEASURED_NONLINEAR_CANDIDATE_STEP", endpoint=accepted.key,
            endpoint_role=role, cache_path=cache_path, endpoints=measured, corrector=correction,
            origin_endpoint=current.key, initial_anchor_endpoint=anchor.key,
            Phi=target.phi(accepted, y), delta_Phi_from_anchor=target.phi(accepted, y)-target.phi(anchor, y),
            actual_step_RMS=target.prior.material_rms(accepted.state.chi, current.state.chi),
            anchor_displacement_RMS=target.prior.material_rms(accepted.state.chi, anchor.state.chi),
            signed_progress=_progress(target.prior, accepted.q, current, velocity)[1],
            prior_penalty=.5*float(accepted.q[:64]@accepted.q[:64]),
            Gaussian_step_norm=float(np.linalg.norm(accepted.q[:64]-current.q[:64])),
            saturation_strata_changed_from_anchor=int(np.sum(strata != initial_strata)),
            saturation_fraction=float(np.mean(strata != 0)),
            minimum_seam_distance=float(np.min(abs(abs(f)-1))),
            projector_rotation_step=previous_rotation,
            projector_rotation_from_anchor=float(np.linalg.norm(new_frame["projector"]-frozen["projector"], ord=2)),
            frozen_dynamic_projector_difference=float(np.linalg.norm(dynamic["projector"]-frozen["projector"], ord=2)),
            Maxwell_secant_relative_error=previous_error, cost=full_cost,
            wall_seconds=time.perf_counter()-start, GT_used_online=False,
            posterior_sample=False, continuous_path_verified=False)
        save(directory/f"step_{step:03d}_FRAME.npz", **new_frame)
        with (directory/"TRACE.jsonl").open("a") as stream:
            from .core import plain
            stream.write(json.dumps(plain(row), allow_nan=False)+"\n")
        rows.append(row)
        current = accepted
        write(directory/"RESULTS.json", dict(parent=parent, arm=arm, rows=rows,
            completed_steps=len(rows), status="IN_PROGRESS", GT_used_online=False, posterior_samples=False))
    result = dict(parent=parent, arm=arm, rows=rows, completed_steps=sum("endpoint" in r for r in rows),
        requested_steps=int(steps), status="COMPLETE" if len(rows) == steps and "endpoint" in rows[-1] else "STOPPED_NO_FEASIBLE_DIRECTION",
        wall_seconds=time.perf_counter()-start,
        cost={k:float(target.book.counts.get(k, 0)-initial_counts.get(k, 0)) for k in ("full_forward_calls", "linear_solve_rhs")},
        GT_used_online=False, posterior_samples=False, exact_fiber_claim=False,
        continuous_path_verified=False, candidate_probabilities=None)
    write(directory/"RESULTS.json", result)
    return result
