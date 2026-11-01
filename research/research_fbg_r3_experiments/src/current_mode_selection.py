"""Common public current pool and two-sided material-task release selection."""
from __future__ import annotations

import numpy as np
from scipy import linalg

try:
    from .current_normal import spatial_modes, current_action, normal_diagnostic, quotient_normal_signature
    from .joint_current_model import JointCurrentModel, whiten_data, whiten_physical, regularizer, physical_task_config
    from .profile_schur import profile_schur, release_candidate, physical_task_gain
except ImportError:
    from current_normal import spatial_modes, current_action, normal_diagnostic, quotient_normal_signature
    from joint_current_model import JointCurrentModel, whiten_data, whiten_physical, regularizer, physical_task_config
    from profile_schur import profile_schur, release_candidate, physical_task_gain


def _orth_columns(U, against=None, count=None):
    U = np.asarray(U, complex)
    if U.ndim == 1:
        U = U[:, None]
    if against is not None and against.shape[1]:
        U = U - against @ (against.conj().T @ U)
        U = U - against @ (against.conj().T @ U)
    if not U.shape[1]:
        return U
    u, s, _ = linalg.svd(U, full_matrices=False, check_finite=False)
    rank = int(np.count_nonzero(s > 1e-11 * s.max(initial=0)))
    return u[:, :min(rank, count) if count is not None else rank]


def _fourier_atoms(shape):
    coordinates = np.indices(shape).reshape(len(shape), -1).T
    frequencies = np.indices(shape).reshape(len(shape), -1).T
    centered = np.where(frequencies <= np.asarray(shape) / 2, frequencies, frequencies - np.asarray(shape))
    order = np.argsort(np.sum(centered ** 2, axis=1), kind="stable")
    for k in order:
        yield np.exp(2j * np.pi * (coordinates / np.asarray(shape)) @ centered[k]) / np.sqrt(np.prod(shape))


def _build_native_tsom_pool(model, Y, m, state, config, cost=None):
    """Separate 6x6 native affine chart: fixed DPCS plus sixteen APCS atoms.

    The first candidate restores precisely the native mf=2 Fourier span. Other
    candidates replace its last orthonormal atom, keeping fifteen APCS atoms.
    The joint anchor is jD (a=0), and is deliberately NOT a physical current.
    """
    try:
        from .cie_tsom_baseline import fourier_apcs
    except ImportError:
        from cie_tsom_baseline import fourier_apcs
    n, p = model.E.shape
    if tuple(model.shape) != (6, 6) or n != 36:
        raise ValueError("native_tsom is a separate 6x6-only variant; do not run a new 4x4 variant")
    if int(config.get("group_atoms", 1)) != 1 or int(config.get("native_mf", 2)) != 2:
        raise ValueError("native_tsom freezes mf=2 and one replacement spatial atom")
    if config.get("operation", "replace") != "replace":
        raise ValueError("native_tsom requires same-rank operation='replace'")
    left, singular, vh = linalg.svd(model.S, full_matrices=False, check_finite=False)
    cutoff = float(config.get("backbone_receiver_rtol", config.get("dpcs_rtol", 1e-3)))
    rank_d = int(np.count_nonzero(singular > cutoff * singular[0]))
    if rank_d != 16:
        raise ValueError(f"native_tsom frozen 6x6 observable DPCS rank is16, got {rank_d}")
    Vd = vh[:rank_d].conj().T
    jD = Vd @ ((left[:, :rank_d].conj().T @ Y) / singular[:rank_d, None])
    projected_fourier, native_info = fourier_apcs((6, 6), 2, Vd)
    Ua, sa, _ = linalg.svd(projected_fourier, full_matrices=False, check_finite=False)
    native_rank = int(np.count_nonzero(sa > max(projected_fourier.shape) * np.finfo(float).eps * max(float(sa[0]), 1.0)))
    if native_rank != 16:
        raise ValueError(f"native mf=2 APCS span must have rank16, got {native_rank}")
    Ua = Ua[:, :native_rank]
    core, old = Ua[:, :15], Ua[:, 15:16]
    against = np.column_stack((Vd, core))
    complete_native = np.column_stack((Vd, Ua))
    def projected_atom(atom):
        atom = np.asarray(atom, complex).reshape(n, 1)
        scale = np.linalg.norm(atom)
        for _ in range(2):
            atom -= against @ (against.conj().T @ atom)
        length = np.linalg.norm(atom)
        return atom / length if length > 1e-11 * scale else None
    # Only additional Fourier directions outside the native mf2 span can test
    # replacement; its own final direction is already included as candidate0.
    all_fourier = list(_fourier_atoms((6, 6)))
    extra_fourier = [f for f in all_fourier if np.linalg.norm(f-complete_native@(complete_native.conj().T@f)) > 1e-11*np.linalg.norm(f)]
    fourier_queue = iter(extra_fourier)
    rng = np.random.default_rng(int(config.get("random_seed", 20261101)))
    krylov = np.sum(jD, axis=1)
    if np.linalg.norm(krylov) < 1e-14:
        krylov = np.sum(model.E, axis=1)
    group_count = int(config.get("pool_groups", 8))
    families = tuple(config.get("native_pool_families", ("fourier", "krylov", "random")))
    if not families or any(f not in ("fourier", "krylov", "random") for f in families):
        raise ValueError("native_pool_families supports Fourier/Krylov/random APCS proposals")
    groups = []
    for index in range(group_count):
        family = "fourier" if index == 0 else families[(index-1) % len(families)]
        atom, fallback = old.copy() if index == 0 else None, 0
        if index:
            for _ in range(24):
                if family == "fourier":
                    try:
                        raw = next(fourier_queue)
                    except StopIteration:
                        raw = rng.normal(size=n)+1j*rng.normal(size=n); fallback += 1
                elif family == "krylov":
                    krylov = model.G @ krylov
                    if cost:
                        cost.add("green_vector_actions")
                    raw = krylov
                else:
                    raw = rng.normal(size=n)+1j*rng.normal(size=n)
                atom = projected_atom(raw)
                if atom is not None:
                    break
        if atom is None:
            raise ValueError("No nonzero native APCS candidate in the retained quotient")
        groups.append({"id": "fourier_native_old" if index == 0 else f"{family}_native_{index}",
                       "family": family, "spatial": atom, "modes": spatial_modes(atom, p),
                       "real_dof": 2*p, "active_real_dof": 16*2*p, "backbone_real_dof": 15*2*p,
                       "origin": "native_mf2_old_atom" if index == 0 else ("seeded_random_fallback" if fallback else family),
                       "fallback_draws": fallback, "restores_native_mf2_span": index == 0})
    data_residual = model.S @ jD-Y
    data_norm = max(float(np.linalg.norm(Y)), 1e-300)
    S_core = model.S @ core
    def tail(SU):
        u, s, _ = linalg.svd(SU, full_matrices=False, check_finite=False)
        # Native APCS receiver actions are analytically zero when all sixteen
        # receiver modes are DPCS. Compare roundoff with S, never with SU alone.
        tol = 1e-11 * max(float(singular[0]), float(s.max(initial=0.0)))
        rank = int(np.count_nonzero(s > tol))
        return rank, float(np.linalg.norm(data_residual-u[:, :rank]@(u[:, :rank].conj().T@data_residual))/data_norm)
    rank_core_data, tail_core = tail(S_core)
    for group in groups:
        rr, tt = tail(np.column_stack((S_core, model.S @ group["spatial"])))
        group.update(receiver_data_rank=rr, receiver_tail_relative=tt,
                     dpcs_projection_error=float(np.linalg.norm(Vd.conj().T @ group["spatial"])))
    initial_defect = jD-np.asarray(m)[:, None]*(model.E+model.G@jD)
    audit = {"affine_mode": "native_tsom", "joint_anchor_is_physical": False,
             "receiver_mode": "fixed_native_DPCS", "receiver_relative_cutoff": cutoff,
             "fixed_dpcs_spatial_rank": rank_d, "fixed_dpcs_real_variable_dof": 0,
             "native_mf": 2, "native_mf2_apcs_rank": native_rank,
             "backbone_spatial_rank": 15, "backbone_real_dof": 15*2*p,
             "active_apcs_spatial_rank": 16, "active_real_dof": 16*2*p,
             "available_dark_spatial_rank": n-rank_d,
             "core_data_rank": rank_core_data, "core_receiver_tail_relative": tail_core,
             "fixed_dpcs_data_relative_residual": float(np.linalg.norm(data_residual)/data_norm),
             "joint_anchor_defect_norm": float(np.linalg.norm(initial_defect)),
             "physical_anchor_current_distance": float(np.linalg.norm(jD-state["j"])),
             "best_candidate_receiver_tail_relative": min(g["receiver_tail_relative"] for g in groups),
             "candidate_receiver_tails": {g["id"]: g["receiver_tail_relative"] for g in groups},
             "candidate_origins": {g["id"]: g["origin"] for g in groups},
             "normal_reference": "raw_unweighted_exact_physical_anchor_diagnostic_only",
             "tail_metric": "native_affine_unwhitened_Frobenius_divided_by_norm_Y"}
    audit["receiver_tail_rank_threshold_relative_to_S"] = 1e-11
    if cost:
        cost.add("dictionary_builds")
        cost.add("receiver_svd")
        cost.add("native_APCS_svd", 2)
        cost.add("fourier_atoms_generated", len(all_fourier)+native_info["fourier_columns"])
        cost.add("DPCS_projection_columns", native_info["fourier_columns"])
        cost.add("dictionary_real_dof", group_count*2*p)
        cost.add("dictionary_backbone_real_dof", 15*2*p)
        cost.add("green_vector_actions", p)
        cost.add("receiver_vector_actions", p+15+group_count)
        cost.add("receiver_tail_svd", group_count+1)
        cost.add("receiver_tail_audit", **audit)
    return {"core": spatial_modes(core, p), "core_spatial": core, "removable": spatial_modes(old, p),
            "groups": groups, "source_count": p, "group_real_dof": 2*p,
            "joint_anchor_current": jD, "fixed_dpcs_current": jD, "fixed_dpcs_basis": Vd,
            "native_apcs_basis": Ua, "native_fourier_basis": projected_fourier,
            "receiver_singular_values": singular, "receiver_audit": audit}


def build_common_pool(model, Y, m, state, config=None, cost=None):
    """All selectors receive these exact same equal-real-capacity candidates.

    Fourier, receiver, Green-Krylov and seeded random atoms use only public
    operators, observations and the common anchor. Acquisition is charged.
    ``backbone_receiver_atoms='observable'`` retains all right singular vectors
    above ``backbone_receiver_rtol`` (default 1e-3, matching mature CIE).
    Integer receiver counts preserve the original experiment configuration.
    """
    cfg = config or {}
    if cfg.get("affine_mode", "physical_anchor") == "native_tsom":
        return _build_native_tsom_pool(model, Y, m, state, cfg, cost)
    n, p = model.E.shape
    shape = tuple(model.shape) if np.prod(model.shape) == n else (n,)
    _, receiver_singular, vh = linalg.svd(model.S, full_matrices=False, check_finite=False)
    receiver = vh.conj().T
    receiver_setting = cfg.get("backbone_receiver_atoms", 1)
    observable = receiver_setting == "observable"
    receiver_rtol = float(cfg.get("backbone_receiver_rtol", cfg.get("dpcs_rtol", 1e-3)))
    if observable:
        if not 0.0 < receiver_rtol < 1.0:
            raise ValueError("Observable receiver cutoff must lie in (0,1)")
        nr = int(np.count_nonzero(receiver_singular > receiver_rtol * receiver_singular[0])) if len(receiver_singular) else 0
    else:
        nr = min(int(receiver_setting), receiver.shape[1])
    nf = int(cfg.get("backbone_fourier_atoms", 1))
    fourier = list(_fourier_atoms(shape))
    core = _orth_columns(np.column_stack([receiver[:, :nr], np.column_stack(fourier[:nf])]) if nf else receiver[:, :nr])
    def candidate_orth(U, against, count=None):
        if not observable:
            return _orth_columns(U, against, count)
        # Large observable backbones make exact-span proposals common. Reference
        # the original atoms when testing quotient rank, never roundoff alone.
        U = np.asarray(U, complex)
        if U.ndim == 1:
            U = U[:, None]
        reference = float(np.linalg.norm(U))
        for _ in range(2):
            U = U - against @ (against.conj().T @ U)
        u, s, _ = linalg.svd(U, full_matrices=False, check_finite=False)
        rank = int(np.count_nonzero(s > 1e-11 * max(reference, float(s.max(initial=0)))))
        return u[:, :min(rank, count) if count is not None else rank]
    # The removable group is explicit and has the same real capacity as a proposal.
    k = int(cfg.get("group_atoms", 1))
    q = int(cfg.get("pool_groups", 8))
    rng = np.random.default_rng(int(cfg.get("random_seed", 20261101)))
    source = model.S.conj().T @ (state["Y"] - Y)
    if cost:
        cost.add("receiver_transpose_vector_actions", p)
    krylov = np.sum(source, axis=1)
    if np.linalg.norm(krylov) < 1e-14:
        krylov = np.sum(state["j"], axis=1)
    sequences = {"fourier": iter(fourier[nf:]), "receiver": iter(receiver[:, nr:].T),
                 "random": None, "krylov": None}
    families = cfg.get("pool_families", ("fourier", "receiver", "krylov", "random"))
    groups = []
    for group in range(q):
        family = families[group % len(families)]
        collected = np.empty((n, 0), complex)
        fallback_draws = 0
        for attempt in range(max(12, 4 * k)):
            if collected.shape[1] == k:
                break
            if family == "random":
                atom = rng.normal(size=n) + 1j * rng.normal(size=n)
            elif family == "krylov":
                krylov = model.G @ krylov
                if cost:
                    cost.add("green_vector_actions")
                atom = krylov
            else:
                try:
                    atom = next(sequences[family])
                except StopIteration:
                    atom = rng.normal(size=n) + 1j * rng.normal(size=n)
                    fallback_draws += 1
            candidate = candidate_orth(atom, np.column_stack((core, collected)), 1)
            if candidate.shape[1]:
                collected = np.column_stack((collected, candidate))
        if collected.shape[1] != k:
            raise ValueError("Requested current group capacity exceeds available quotient rank")
        groups.append({"id": f"{family}_{group}", "family": family,
                       "spatial": collected, "modes": spatial_modes(collected, p), "real_dof": 2 * k * p,
                       "fallback_draws": fallback_draws,
                       "origin": "seeded_random_fallback" if fallback_draws else family})
    removable_spatial = candidate_orth(np.column_stack(fourier[nf:nf + max(4 * k, k)]), core, k)
    if removable_spatial.shape[1] != k:
        raise ValueError("No full-rank removable current group")
    # Exact finite-space data representability audit, before any material fit.
    # Each source has independent complex coefficients, hence the same complex
    # range S U is projected source by source. This is a raw Frobenius/Y tail.
    residual = np.asarray(state["Y"]) - Y
    data_norm = max(float(np.linalg.norm(Y)), 1e-300)
    S_core = model.S @ core
    def receiver_tail(SU):
        left, singular, _ = linalg.svd(SU, full_matrices=False, check_finite=False)
        threshold = np.finfo(float).eps * max(SU.shape) * float(singular.max(initial=0.0))
        rank = int(np.count_nonzero(singular > threshold))
        Qs = left[:, :rank]
        tail = residual - Qs @ (Qs.conj().T @ residual)
        return rank, float(np.linalg.norm(tail) / data_norm)
    core_data_rank, core_tail = receiver_tail(S_core)
    core_real_dof = 2 * p * core.shape[1]
    for group in groups:
        rank, tail = receiver_tail(np.column_stack((S_core, model.S @ group["spatial"])))
        group.update(receiver_data_rank=rank, receiver_tail_relative=tail,
                     backbone_real_dof=core_real_dof,
                     active_real_dof=core_real_dof + group["real_dof"])
    receiver_audit = {"receiver_mode": "observable" if observable else "fixed_integer",
                      "receiver_relative_cutoff": receiver_rtol if observable else None,
                      "receiver_atoms": nr, "fourier_atoms": nf,
                      "backbone_spatial_rank": core.shape[1], "backbone_real_dof": core_real_dof,
                      "active_real_dof": core_real_dof + 2 * k * p,
                      "core_data_rank": core_data_rank, "core_receiver_tail_relative": core_tail,
                      "best_candidate_receiver_tail_relative": min((g["receiver_tail_relative"] for g in groups), default=core_tail),
                      "candidate_receiver_tails": {g["id"]: g["receiver_tail_relative"] for g in groups},
                      "candidate_origins": {g["id"]: g["origin"] for g in groups},
                      "tail_metric": "unwhitened_Frobenius_divided_by_norm_Y",
                      "tail_numerical_rank_relative_threshold": "eps*max(SU.shape)"}
    if cost:
        cost.add("dictionary_builds")
        cost.add("receiver_svd")
        cost.add("fourier_atoms_generated", n)
        cost.add("dictionary_real_dof", sum(g["real_dof"] for g in groups))
        cost.add("dictionary_backbone_real_dof", core_real_dof)
        cost.add("receiver_vector_actions", core.shape[1] + q * k)
        cost.add("receiver_tail_svd", q + 1)
        cost.add("receiver_tail_audit", **receiver_audit)
    return {"core": spatial_modes(core, p), "core_spatial": core,
            "removable": spatial_modes(removable_spatial, p),
            "groups": groups, "source_count": p, "group_real_dof": 2 * k * p,
            "receiver_singular_values": receiver_singular, "receiver_audit": receiver_audit}


def paid_physical_jvp(model, m, h, state, config, cost=None):
    if cost:
        cost.check(extra_rhs=model.E.shape[1])
        cost.add("physical_jvp_calls")
        cost.add("expected_complex_rhs", model.E.shape[1])
    if hasattr(model, "tangent"):
        tangent = model.tangent(m, h)[0]
    else:
        rhs = np.asarray(h)[:, None] * state["Q"]
        if state.get("lu") is not None:
            tangent = linalg.lu_solve(state["lu"], rhs, check_finite=False)
        else:
            tangent = linalg.solve(np.eye(model.n) - m[:, None] * model.G, rhs, check_finite=False)
            if cost:
                cost.add("expected_lu")
        if hasattr(model, "cost"):
            model.cost["actual_complex_RHS"] += model.E.shape[1]
            model.cost["tangent_RHS"] += model.E.shape[1]
    if cost:
        cost.add("receiver_vector_actions", model.E.shape[1])
    return whiten_physical(model.pack(model.S @ tangent), config)


def evaluate_modes(model, Y, m, state, core, groups, config=None, cost=None, old_extra=None, joint_current=None):
    cfg = dict(config or {})
    joint = JointCurrentModel(model, Y, cfg, cost)
    joint_anchor = state["j"] if joint_current is None else np.asarray(joint_current, complex)
    lin = joint.linearize(m, joint_anchor, core)
    prof = profile_schur(lin["b"], lin["B"], lin["A"], lin["Lambda"], lin["ell"],
                         cfg.get("lm", 1e-4), m, cfg.get("box", (.02, 20.0)), cost)
    baseline = None
    if old_extra is not None:
        baseline = release_candidate(prof, joint.current_matrix(m, old_extra), release="hard", cost=cost)
        old_p = baseline["p"]
        old_action = current_action(core, baseline["c_U"]) + current_action(old_extra, baseline["c_C"])
    else:
        old_p, old_action = prof["p"], current_action(core, prof["c"])
    _, ellF, LambdaF = regularizer(model, m, cfg)
    rF = whiten_physical(model.pack(state["Y"] - Y), cfg)
    Jprovided = cfg.get("physical_J")
    def proposal_jvp(h):
        if Jprovided is not None:
            return np.asarray(Jprovided) @ h
        # Public frozen-anchor Born map is a proposal model, with no error bound.
        if cost:
            cost.add("receiver_vector_actions", model.E.shape[1])
            cost.add("born_proposal_actions")
        return whiten_physical(model.pack(model.S @ (h[:, None] * state["Q"])), cfg)
    Jp0 = proposal_jvp(old_p)
    actions = []
    for group in groups:
        C = group["modes"]
        A_C = joint.current_matrix(m, C, lin["w"])
        rel = release_candidate(prof, A_C, cfg.get("release", "soft"), cfg.get("eta_target", .1), cfg.get("gamma"), cost)
        dp = rel["p"] - old_p
        actual = current_action(core, rel["c_U"]) + current_action(C, rel["c_C"]) - old_action
        normal = normal_diagnostic(model.G, m, state["Q"], actual, cfg.get("beta", 0.0), cost=cost)
        signature = quotient_normal_signature(model.G, m, state["Q"], core, C, cost=cost)
        prescore = physical_task_gain(rF, Jp0, proposal_jvp(dp), old_p, dp, ellF, LambdaF)
        spatial = group.get("spatial")
        energy = float(np.linalg.norm(spatial.conj().T @ state["j"]) ** 2) if spatial is not None else float(np.linalg.norm(C) ** 2)
        actions.append({"group": group, "release": rel, "delta_p": dp, "old_p": old_p,
                        "current_action": actual, "normal": normal, "signature": signature,
                        "joint_anchor_is_physical": joint_current is None,
                        "normal_reference": "raw_unweighted_exact_physical_anchor_diagnostic_only" if joint_current is not None else "raw_unweighted_exact_physical_anchor",
                        "proposal_gain": prescore, "proposal_model": "provided_exact_J" if Jprovided is not None else "Born_no_error_bound",
                        "energy_score": energy, "physical_gain": None})
    selector = cfg.get("selector", "material_step")
    if selector in ("material_step", "task", "schur"):
        ordered = sorted(actions, key=lambda a: a["proposal_gain"], reverse=True)
    elif selector == "energy":
        ordered = sorted(actions, key=lambda a: a["energy_score"], reverse=True)
    elif selector == "normal":
        ordered = sorted(actions, key=lambda a: a["normal"]["normal_energy"], reverse=True)
    elif selector == "random":
        ordered = [actions[i] for i in np.random.default_rng(int(cfg.get("random_seed", 20261101))).permutation(len(actions))]
    elif selector in ("fourier", "krylov", "receiver"):
        # Fixed native-family control, retaining original public pool order.
        # Every pool group above was still acquired, linearized and costed.
        ordered = sorted(actions, key=lambda a: a["group"]["family"] != selector)
    else:
        raise ValueError(f"Unknown current selector: {selector}")
    ordered = [a for a in ordered if not a["release"]["duplicate"]]
    finalists = ordered[:int(cfg.get("finalist_jvps", 2))]
    if finalists:
        # Old-step action is cached once; each finalist pays p more tangent RHS.
        Jp = paid_physical_jvp(model, m, old_p, state, cfg, cost)
        for action in finalists:
            Jdp = paid_physical_jvp(model, m, action["delta_p"], state, cfg, cost)
            action["physical_gain"] = physical_task_gain(rF, Jp, Jdp, old_p, action["delta_p"], ellF, LambdaF)
        if selector in ("material_step", "task", "schur"):
            finalists.sort(key=lambda a: a["physical_gain"], reverse=True)
    selected = finalists[0] if finalists else None
    return {"profile": prof, "baseline_release": baseline, "old_p": old_p,
            "old_current_action": old_action, "actions": actions, "finalists": finalists,
            "selected": selected, "linearization": lin, "joint_anchor_current": joint_anchor,
            "joint_anchor_is_physical": joint_current is None}


select_modes = evaluate_modes
