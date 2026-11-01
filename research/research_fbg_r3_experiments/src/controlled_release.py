"""N1 bounded auxiliary-current release, finite seeds and physical readback.

No labels or hidden endpoints enter this module. The physical corrector may be
provided as a callback, or applied uniformly by the experiment orchestrator.
"""
from __future__ import annotations

import numpy as np
from scipy import linalg

try:
    from .current_normal import current_action, normal_diagnostic
    from .joint_current_model import JointCurrentModel, ReleaseCost, ReleaseBudgetExceeded, regularizer, whiten_data, whiten_physical
    from .profile_schur import profile_schur, release_candidate, least_squares, solve_symmetric, box_quadratic_step
    from .current_mode_selection import build_common_pool, evaluate_modes
except ImportError:
    from current_normal import current_action, normal_diagnostic
    from joint_current_model import JointCurrentModel, ReleaseCost, ReleaseBudgetExceeded, regularizer, whiten_data, whiten_physical
    from profile_schur import profile_schur, release_candidate, least_squares, solve_symmetric, box_quadratic_step
    from current_mode_selection import build_common_pool, evaluate_modes


def finite_normal_seed(profile, action, model, m, state, core, config=None, cost=None, joint_current=None):
    """A common signed tau, held in ORIGINAL physical current coordinates.

    Fixed-radius mode records an unattainable requested scale. With the frozen
    bounded-feasibility repair, the requested RMS is a declared maximum: every
    trial QP/Green action is charged, and achieved scale is explicitly reported.
    No claim of global maximality or monotonic feasibility is made.
    """
    cfg = config or {}
    joint_anchor = state["j"] if joint_current is None else np.asarray(joint_current, complex)
    rel = action["release"]
    D, v, R = rel["D"], rel["v"], rel["R"]
    if not len(v):
        return {"status": "DUPLICATE_CURRENT_RANGE", "branches": []}
    H, g, p = profile["H"], profile["g"], profile["p"]
    w = v + D @ p
    rule = cfg.get("seed_rule", "task")
    if rule in ("random", "fourier"):
        # Fourier uses the selected Fourier group and a fixed real coordinate;
        # random uses the same group's QR coordinates and a frozen random draw.
        if rule == "random":
            w = np.random.default_rng(int(cfg.get("random_seed", 20261101)) + 71).normal(size=len(v))
        else:
            w = np.zeros(len(v)); w[0] = 1.0
    elif np.linalg.norm(w) <= 1e-12:
        K = D @ solve_symmetric(H, D.T)
        w = linalg.eigh(0.5 * (K + K.T), check_finite=False)[1][:, -1]
    if np.linalg.norm(w) <= 1e-14:
        return {"status": "ZERO_SEED_DIRECTION", "branches": []}
    w = w / np.linalg.norm(w)
    response = -solve_symmetric(H, D.T @ w)
    response_rms = float(np.linalg.norm(response) / np.sqrt(len(m)))
    delta = float(cfg.get("seed_material_rms", cfg.get("material_scale", .2)))
    if response_rms <= 1e-14:
        return {"status": "ZERO_FIRST_ORDER_MATERIAL_RESPONSE", "branches": [], "response_rms": response_rms}
    tau_requested = delta / response_rms
    C = action["group"]["modes"]
    ref_current = max(float(cfg.get("seed_current_scale", np.linalg.norm(joint_anchor))), 1e-14)
    ref_defect = max(float(cfg.get("seed_defect_scale", ref_current)), 1e-14)
    trial_count = 0
    search_enabled = bool(cfg.get("bounded_seed_scale", True))

    def build(sign, tau, is_search=False):
        nonlocal trial_count
        trial_count += 1
        u = sign * float(tau) * w
        coeff_C = least_squares(R, u)
        h = box_quadratic_step(H, g + D.T @ u, m, cfg.get("box", (.02, 20.0)))
        coeff_U = -least_squares(profile["A"], profile["b"] + profile["B"] @ h + rel["A_C"] @ coeff_C)
        V = current_action(core, coeff_U) + current_action(C, coeff_C)
        mseed, jseed = m + h, joint_anchor + V
        actual_scale = float(np.linalg.norm(h - p) / np.sqrt(len(m)))
        defect = jseed - mseed[:, None] * (model.E + model.G @ jseed)
        if cost:
            cost.add("green_vector_actions", model.E.shape[1])
            cost.add("material_qp_solves")
            cost.add("seed_feasibility_trials", sign=sign, tau=float(tau), is_search=is_search)
        current_ratio, defect_ratio = float(np.linalg.norm(V) / ref_current), float(np.linalg.norm(defect) / ref_defect)
        normal = normal_diagnostic(model.G, m, state["Q"], V, cfg.get("beta", 0.0), cost=cost)
        reasons = []
        if search_enabled:
            if actual_scale > delta * (1.0 + 1e-8):
                reasons.append("MATERIAL_MAXIMUM_EXCEEDED")
        elif abs(actual_scale - delta) > float(cfg.get("seed_scale_tolerance", cfg.get("scale_tolerance", .2))) * delta:
            reasons.append("MATERIAL_SCALE_UNREACHED_ON_BOX")
        if current_ratio > float(cfg.get("seed_current_cap", 3.0)):
            reasons.append("CURRENT_BUDGET_REJECTED")
        if defect_ratio > float(cfg.get("seed_defect_cap", 3.0)):
            reasons.append("DEFECT_BUDGET_REJECTED")
        if not is_search and (normal["no_op"] or normal["normal_energy"] <= 1e-22):
            reasons.append("NO_NORMAL_ACTIVATION")
        status = reasons[0] if reasons else "READY"
        return {"sign": sign, "tau": float(tau), "u": u, "m": mseed, "j": jseed,
                         "seed_coeff_C": coeff_C, "seed_coeff_U": coeff_U,
                         "requested_material_rms": delta, "actual_material_rms": actual_scale,
                         "current_ratio": current_ratio, "defect_ratio": defect_ratio,
                         "normal": normal, "status": status, "rejection_reasons": reasons}

    branches = [build(sign, tau_requested, is_search=search_enabled) for sign in (1, -1)]
    requested_receipt = [{k: b[k] for k in ("sign", "tau", "actual_material_rms", "current_ratio", "defect_ratio", "rejection_reasons")} for b in branches]
    tau = tau_requested
    scale_status = "REQUESTED_SCALE"
    if search_enabled and any(b["rejection_reasons"] for b in branches):
        base = [build(sign, 0.0, is_search=True) for sign in (1, -1)]
        if any(b["rejection_reasons"] for b in base):
            return {"status": "NO_FEASIBLE_TESTED_SCALE", "branches": [], "tau_requested": tau_requested,
                    "requested_receipt": requested_receipt, "zero_scale_receipt": [{k: b[k] for k in ("sign", "current_ratio", "defect_ratio", "rejection_reasons")} for b in base],
                    "feasibility_evaluations": trial_count, "seed_rule": rule}
        lower, upper, branches = 0.0, tau_requested, base
        # Both signs are tested together; this returns one common absolute tau.
        # Eight finite bisection trials do not certify interval monotonicity.
        for _ in range(min(8, int(cfg.get("seed_scale_trials", 8)))):
            trial_tau = .5 * (lower + upper)
            tested = [build(sign, trial_tau, is_search=True) for sign in (1, -1)]
            if all(not b["rejection_reasons"] for b in tested):
                lower, branches = trial_tau, tested
            else:
                upper = trial_tau
        tau, scale_status = lower, "CAP_LIMITED"
    for branch in branches:
        branch["scale_status"] = scale_status
        branch["tau_requested"] = tau_requested
        if search_enabled and branch["actual_material_rms"] < float(cfg.get("seed_min_material_rms", .02)):
            branch["rejection_reasons"].append("NO_OP_MATERIAL_SCALE")
        if branch["normal"]["no_op"] or branch["normal"]["normal_energy"] <= 1e-22:
            branch["rejection_reasons"].append("NO_NORMAL_ACTIVATION")
        branch["status"] = branch["rejection_reasons"][0] if branch["rejection_reasons"] else "READY"
        if cost:
            cost.add("finite_seed_proposals", sign=branch["sign"], tau=float(tau), status=branch["status"])
    return {"status": "BUILT", "branches": branches, "tau": float(tau), "tau_requested": tau_requested,
            "w": w, "scale_status": scale_status, "requested_receipt": requested_receipt,
            "feasibility_evaluations": trial_count, "response_rms": response_rms, "seed_rule": rule}


def _rollout(model, Y, m0, j0, core, candidate, cfg, cost, history, tag, fixed=False, iterations=0):
    m, j = np.asarray(m0, float).copy(), np.asarray(j0, complex).copy()
    config = dict(cfg)
    lm = float(cfg.get("lm", 1e-4))
    initial_norm = max(float(cfg.get("rollout_current_scale", np.linalg.norm(j0))), 1e-14)
    defect_scale = max(float(cfg.get("rollout_defect_scale", initial_norm)), 1e-14)
    stop = "ITERATION_LIMIT"
    for iteration in range(iterations):
        cost.check()
        # Declared rho stages change the objective; no derivative history is kept.
        schedule = cfg.get("rho_schedule")
        if schedule and not fixed:
            config["rho"] = float(schedule[min(len(schedule)-1, iteration * len(schedule) // max(1, iterations))])
        joint = JointCurrentModel(model, Y, config, cost)
        lin = joint.linearize(m, j, core)
        prof = profile_schur(lin["b"], lin["B"], lin["A"], lin["Lambda"], lin["ell"],
                             lm, m, cfg["box"], cost)
        if fixed or candidate is None or cfg.get("release") == "none":
            h, V = prof["p"], current_action(core, prof["c"])
            gamma, eta = None, None
        else:
            A_C = joint.current_matrix(m, candidate, lin["w"])
            rel = release_candidate(prof, A_C, cfg.get("release", "soft"), cfg.get("eta_target", .1), cfg.get("gamma"), cost)
            h = rel["p"]
            V = current_action(core, rel["c_U"]) + current_action(candidate, rel["c_C"])
            gamma, eta = rel["gamma"], rel["eta"]
        cost.add("joint_iterations")
        step_norm = float(np.linalg.norm(h) + np.linalg.norm(V))
        if step_norm <= cfg.get("joint_step_tolerance", 1e-10):
            stop = "ZERO_JOINT_STEP"
            history.append({"tag": tag, "iteration": iteration, "stage": "hold" if fixed else "release",
                            "status": stop, "step_norm": step_norm, "auxiliary_objective": lin["objective"]})
            break
        accepted = False
        trial, ls = lin, -1
        for ls in range(int(cfg.get("joint_line_search", 4))):
            cost.check()
            step = 0.5 ** ls
            trial_m, trial_j = m + step * h, j + step * V
            if np.linalg.norm(trial_j - j0) / initial_norm > float(cfg.get("rollout_current_cap", 10.0)):
                continue
            trial = joint.residual(trial_m, trial_j)
            if trial["defect_norm"] / defect_scale > float(cfg.get("rollout_defect_cap", 10.0)):
                continue
            # The increment ridge controls the local proposal, while acceptance
            # uses the true fixed-stage nonlinear auxiliary objective.
            if np.isfinite(trial["objective"]) and trial["objective"] <= lin["objective"] - cfg.get("joint_accept_margin", 1e-12):
                m, j, accepted = trial_m, trial_j, True
                lm = max(lm / 2.0, float(cfg.get("lm_min", 1e-8)))
                break
        history.append({"tag": tag, "iteration": iteration, "stage": "hold" if fixed else "release",
                        "status": "ACCEPTED" if accepted else "REJECTED", "line_search_trials": ls + 1,
                        "step_norm": step_norm, "rho": joint.rho, "gamma": gamma, "eta": eta,
                        "auxiliary_objective": trial["objective"] if accepted else lin["objective"],
                        "auxiliary_data_loss": trial["data_loss"] if accepted else lin["data_loss"],
                        "auxiliary_state_loss": trial["state_loss"] if accepted else lin["state_loss"],
                        "defect_norm": trial["defect_norm"] if accepted else lin["defect_norm"]})
        if not accepted:
            lm *= 10.0
            if lm > float(cfg.get("lm_max", 1e6)):
                stop = "AUXILIARY_STEP_REJECTED"
                break
    return {"m": m, "j": j, "stop": stop, "config": config}


def _action_record(action):
    rel, normal = action["release"], action["normal"]
    return {"group": action["group"]["id"], "family": action["group"]["family"],
            "real_current_dof": rel["real_dof"], "residual_rank": rel["residual_rank"],
            "eta0": rel["eta0"], "eta": rel["eta"], "gamma": rel["gamma"] if np.isfinite(rel["gamma"]) else None,
            "proposal_gain": action["proposal_gain"], "proposal_model": action["proposal_model"],
            "physical_gain": action["physical_gain"], "energy_score": action["energy_score"],
            "normal_fraction": normal["normal_fraction"], "normal_energy": normal["normal_energy"],
            "normal_trace": action["signature"]["normal_trace"],
            "normal_reference": action.get("normal_reference", "raw_unweighted_exact_physical_anchor"),
            "joint_anchor_is_physical": action.get("joint_anchor_is_physical", True),
            "material_step_rms": float(np.linalg.norm(rel["p"]) / np.sqrt(len(rel["p"]))),
            "delta_material_step_rms": float(np.linalg.norm(action["delta_p"]) / np.sqrt(len(rel["p"]))),
            "duplicate_current_range": rel["duplicate"]}


def run_release(model, Y, m0, config=None):
    """Return {m, candidates, history, cost}; m is physically read back.

    Main controls: selector (material_step/energy/normal/random), release
    (hard/soft/none), operation (augment/replace/dressing), finite_seed,
    pool_groups, group_atoms, episodes, hold_steps, joint_steps, finite scale,
    primitive budgets, and optional common corrector callback.
    """
    cfg = dict(config or {})
    native_affine = cfg.get("affine_mode", "physical_anchor") == "native_tsom"
    if native_affine:
        cfg.setdefault("operation", "replace")
        cfg.setdefault("chart", "cie")
        cfg.setdefault("beta", .2)
        cfg.setdefault("rho", 1.0)
        cfg.setdefault("state_weight", 1.0)
        cfg.setdefault("state_source_whitener", 1/np.maximum(np.linalg.norm(model.E, axis=0), 1e-150))
        source_data = 1/np.maximum(np.linalg.norm(Y, axis=0), 1e-150)
        cfg.setdefault("data_whitener", model.pack((1+1j)*np.ones_like(Y)*source_data[None, :]))
        cfg.setdefault("physical_data_scale", max(float(np.linalg.norm(Y)), 1e-14))
        cfg.setdefault("finite_seed", False)
        cfg.setdefault("local_gain_min", -np.inf)
    cfg.setdefault("box", cfg.get("material_box", (.02, 20.0)))
    cfg.setdefault("data_scale", max(float(np.linalg.norm(Y)), 1e-14))
    cfg.setdefault("rho", 1.0 / max(float(np.linalg.norm(model.E) ** 2), 1e-14))
    cfg.setdefault("release", "soft")
    cfg.setdefault("selector", "material_step")
    cfg.setdefault("operation", "augment")
    cost = ReleaseCost(model, cfg)
    history, candidates = [], []
    m = np.asarray(m0, float).copy()
    if np.any(m < cfg["box"][0]) or np.any(m > cfg["box"][1]):
        raise ValueError("Common physical anchor lies outside the declared material box")
    incumbent_score = np.inf
    stop = "EPISODE_LIMIT"
    active_fixed_dpcs_basis, active_fixed_dpcs_current = None, None

    def checkpoint(mm, jj, tag, anchor_m, aux_cfg=None):
        nonlocal m, incumbent_score
        aux_m = np.asarray(mm, float).copy()
        aux = JointCurrentModel(model, Y, aux_cfg or cfg, cost).residual(mm, jj)
        cost.check()
        correction = None
        callback = cfg.get("corrector")
        if callable(callback):
            cost.add("common_corrector_calls")
            correction = callback(model, Y, mm, cfg.get("corrector_config", {}))
            mm = np.asarray(correction["m"], float)
        state = cost.state(mm, purpose=tag)
        rr = whiten_physical(model.pack(state["Y"] - Y), cfg)
        reg = regularizer(model, mm, cfg)[0]
        value = float(0.5 * rr @ rr + reg)
        raw_relative = float(np.linalg.norm(state["Y"] - Y) / max(np.linalg.norm(Y), 1e-14))
        duplicate = np.linalg.norm(mm - anchor_m) / np.sqrt(len(mm)) < cfg.get("duplicate_material_rms", 1e-5)
        duplicate = duplicate or any(np.linalg.norm(mm - c["m"]) / np.sqrt(len(mm)) < cfg.get("duplicate_material_rms", 1e-5) for c in candidates)
        threshold = cfg.get("qualification_residual")
        status = "DUPLICATE" if duplicate else ("QUALIFIED" if threshold is not None and raw_relative <= threshold else "PHYSICAL_VALIDATED")
        if threshold is not None and raw_relative > threshold:
            status = "PHYSICAL_QUALIFICATION_FAILED"
        record = {"m": np.asarray(mm).copy(), "tag": tag, "status": status,
                  "physical_score": value, "relative_residual": raw_relative,
                  "auxiliary_objective": aux["objective"], "auxiliary_data_loss": aux["data_loss"],
                  "auxiliary_state_loss": aux["state_loss"], "auxiliary_defect_norm": aux["defect_norm"],
                  "auxiliary_relative_residual": float(np.linalg.norm(model.S @ jj - Y) / max(np.linalg.norm(Y), 1e-14)),
                  "m_auxiliary": aux_m,
                  "j_auxiliary": np.asarray(jj).copy(),
                  "material_rms_from_anchor": float(np.linalg.norm(mm-anchor_m)/np.sqrt(len(mm))),
                  "corrector": correction}
        if active_fixed_dpcs_basis is not None:
            record["fixed_dpcs_change_norm"] = float(np.linalg.norm(active_fixed_dpcs_basis.conj().T @ (jj-active_fixed_dpcs_current)))
            record["affine_mode"] = "native_tsom"
        candidates.append(record)
        cost.add("receiver_vector_actions", model.E.shape[1])
        history.append({k: v for k, v in record.items() if k not in ("m", "m_auxiliary", "j_auxiliary", "corrector")})
        if np.isfinite(value) and value < incumbent_score:
            incumbent_score, m = value, np.asarray(mm).copy()
        return record

    try:
        initial = cost.state(m, purpose="common_exact_anchor")
        rr = whiten_physical(model.pack(initial["Y"]-Y), cfg)
        incumbent_score = float(.5 * rr @ rr + regularizer(model, m, cfg)[0])
        history.append({"tag": "common_exact_anchor", "physical_score": incumbent_score,
                        "relative_residual": float(np.linalg.norm(initial["Y"]-Y)/max(np.linalg.norm(Y),1e-14)),
                        "current_dof_per_complex_spatial_atom": 2 * model.E.shape[1]})
        for episode in range(int(cfg.get("episodes", 1))):
            cost.check()
            anchor_m = m.copy()
            state = cost.state(anchor_m, purpose=f"episode_{episode}_anchor")
            dictionary = cfg.get("common_pool") or build_common_pool(model, Y, anchor_m, state, cfg, cost)
            core, groups = dictionary["core"], dictionary["groups"]
            joint_anchor = dictionary.get("joint_anchor_current", state["j"])
            active_fixed_dpcs_basis = dictionary.get("fixed_dpcs_basis")
            active_fixed_dpcs_current = dictionary.get("fixed_dpcs_current")
            if cfg["operation"] == "dressing":
                rng = np.random.default_rng(int(cfg.get("random_seed", 20261101)) + 1701)
                groups = [dict(group, modes=group["modes"] + np.einsum("npd,dk->npk", core, .25 * rng.normal(size=(core.shape[2], group["real_dof"])))) for group in groups]
            old_extra = dictionary["removable"] if cfg["operation"] == "replace" else None
            scoring = evaluate_modes(model, Y, anchor_m, state, core, groups, cfg, cost, old_extra,
                                     joint_current=joint_anchor if "joint_anchor_current" in dictionary else None)
            history.append({"tag": "joint_anchor", "episode": episode,
                            "affine_mode": cfg.get("affine_mode", "physical_anchor"),
                            "joint_anchor_is_physical": scoring["joint_anchor_is_physical"],
                            "joint_defect_norm": scoring["linearization"]["defect_norm"],
                            "joint_data_loss": scoring["linearization"]["data_loss"],
                            "joint_state_loss": scoring["linearization"]["state_loss"],
                            "receiver_audit": dictionary.get("receiver_audit")})
            for action in scoring["actions"]:
                history.append({"tag": "mode_score", "episode": episode, "operation": cfg["operation"], **_action_record(action)})
            selected = scoring["selected"]
            if selected is None:
                stop = "NO_NONREDUNDANT_CURRENT_CANDIDATE"
                break
            C = selected["group"]["modes"]
            history.append({"tag": "mode_selected", "episode": episode,
                            "capacity_before": core.shape[2] + (old_extra.shape[2] if old_extra is not None else 0),
                            "capacity_after": core.shape[2] + C.shape[2], **_action_record(selected)})
            if cfg.get("execute_local", True) and selected["physical_gain"] is not None and selected["physical_gain"] > cfg.get("local_gain_min", 0.0):
                local = _rollout(model, Y, anchor_m, joint_anchor, core, C, cfg, cost, history,
                                 f"episode_{episode}_local", iterations=int(cfg.get("joint_steps", 9)))
                history.append({"tag": "rollout_complete", "episode": episode, "stop": local["stop"],
                                "declared_joint_steps": int(cfg.get("joint_steps", 9))})
                checkpoint(local["m"], local["j"], f"episode_{episode}_local", anchor_m, local["config"])
            elif not cfg.get("finite_seed", True):
                # Essential exact-root no-seed control: calculate the actual zero
                # local step, retaining it without inserting an arbitrary nullstep.
                history.append({"tag": "local_slot_no_positive_physical_gain", "episode": episode,
                                "physical_gain": selected["physical_gain"], "step_rms": float(np.linalg.norm(selected["release"]["p"])/np.sqrt(len(m)))})
            if cfg.get("finite_seed", True):
                seeds = finite_normal_seed(scoring["profile"], selected, model, anchor_m, state, core, cfg, cost,
                                           joint_current=joint_anchor if "joint_anchor_current" in dictionary else None)
                history.append({"tag": "finite_seed_setup", "episode": episode,
                                "status": seeds["status"], "tau": seeds.get("tau"), "tau_requested": seeds.get("tau_requested"),
                                "scale_status": seeds.get("scale_status"), "requested_receipt": seeds.get("requested_receipt"),
                                "zero_scale_receipt": seeds.get("zero_scale_receipt"),
                                "feasibility_evaluations": seeds.get("feasibility_evaluations"), "seed_rule": cfg.get("seed_rule", "task")})
                for seed in seeds["branches"]:
                    label = f"episode_{episode}_seed_{seed['sign']:+d}"
                    history.append({"tag": label, "status": seed["status"], "tau": seed["tau"],
                                    "requested_material_rms": seed["requested_material_rms"], "actual_material_rms": seed["actual_material_rms"],
                                    "current_ratio": seed["current_ratio"], "defect_ratio": seed["defect_ratio"],
                                    "normal_fraction": seed["normal"]["normal_fraction"], "rejection_reasons": seed["rejection_reasons"],
                                    "scale_status": seed["scale_status"], "tau_requested": seed["tau_requested"]})
                    if seed["status"] != "READY":
                        cost.add("seed_branches_rejected", reason=seed["status"], sign=seed["sign"])
                        continue
                    cost.add("finite_seed_insertions", sign=seed["sign"], tau=float(seed["tau"]))
                    try:
                        held = _rollout(model, Y, seed["m"], seed["j"], core, None, cfg, cost, history, label,
                                        fixed=True, iterations=int(cfg.get("hold_steps", 3)))
                        free = _rollout(model, Y, held["m"], held["j"], core, C, cfg, cost, history, label,
                                        iterations=int(cfg.get("joint_steps", 9)))
                        checkpoint(free["m"], free["j"], label, anchor_m, free["config"])
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                        history.append({"tag": label, "status": "NUMERICAL_FAILURE", "error": str(exc)})
                        cost.add("failed_branches", reason=str(exc), tag=label)
    except ReleaseBudgetExceeded as exc:
        stop = "BUDGET:" + str(exc)
        history.append({"tag": "stopped", "status": stop})
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        stop = "NUMERICAL_FAILURE"
        history.append({"tag": "stopped", "status": stop, "error": str(exc)})
        cost.add("failed_episodes", reason=str(exc))
    return {"m": m, "candidates": candidates, "history": history, "cost": cost.snapshot(),
            "stop": stop, "physical_score": incumbent_score,
            "config": {k: v for k, v in cfg.items() if k not in ("common_pool", "corrector", "regularizer", "physical_J", "data_whitener")}}


run = run_release
