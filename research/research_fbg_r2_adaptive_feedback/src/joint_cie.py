"""Exploratory classical CSI/CIE joint-current proposal, not a new theory.

For fixed beta in {0, 0.2}, define r=m/(1+beta*m), H=G+beta*I,
V=E+HJ and the state defect C=J-r*V=(I-r*H)J-r*E. The auxiliary
objective alternates exact least-squares current and real rowwise material
updates. It differs from the feedback Gram objective: J is independent of the
physical Maxwell solve inside each auxiliary step. Full physics is solved only
for the initial point and paid checkpoints. The material step omits the common
spatial regularizer and is a proposal; full-physics checkpoints are ranked
with that common regularizer, and callers may apply a physical corrector.
No convergence or basin claim follows.
"""
from collections import Counter
import time

import numpy as np
import torch


def _host_scalar(x):
    return float(x)


def propose(shared, target, m0, wall_s, config=None):
    """Run a bounded, single-stream joint-current proposal and return JSON.

    Parameters
    ----------
    shared : research_z_alg_cuda_pipeline.src.device.SharedOperators
        Resident complex128 G,E0,S, eye, SH and one CUDA device.
    target : array or tensor
        Public complex observations with shape (receivers,sources).
    m0 : array or tensor
        Initial real pixel material, shape (n,). The initial physical current
        uses one charged Maxwell LU and all sources as one RHS block.
    wall_s : float
        Local bounded wall budget, including setup and all host waits.
    config : mapping
        ``beta`` 0 or .2, ``max_iter`` <=100, ``checkpoint_every`` default 10,
        ``material_lower``/``material_upper`` default .02/20. The fixed
        data scale is ||Y||²; the fixed state scale is ||r_initial E||².

    A step solves the mature Cholesky system
      N=S*S/||Y||² + A*A/||r0 E||²,
      N J=S*Y/||Y||² + A*(r E)/||r0 E||²,
    where A=I-rH, followed by r_i=Re<V_i,J_i>/||V_i||² clipped to the
    exact mapped original material box. Each full-physics checkpoint is
    charged separately. ``material_candidates`` contain only fully scored
    Maxwell materials, including the common R4 spatial regularizer. A failed
    Cholesky ends the run; no hidden ridge is added.
    """
    cfg = dict(config or {})
    beta = float(cfg.get("beta", 0.))
    if beta not in (0., 0.2):
        raise ValueError("this exploratory comparison supports beta=0 or beta=0.2")
    wall_s = float(wall_s)
    if wall_s <= 0:
        raise ValueError("wall_s must be positive")
    max_iter = min(100, max(0, int(cfg.get("max_iter", 100))))
    checkpoint_every = int(cfg.get("checkpoint_every", 10))
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")
    lower = float(cfg.get("material_lower", .02))
    upper = float(cfg.get("material_upper", 20.))
    if not (0 <= lower < upper and 1 + beta * lower > 0 and 1 + beta * upper > 0):
        raise ValueError("invalid physical material box")
    r_lower, r_upper = lower / (1 + beta * lower), upper / (1 + beta * upper)
    n, ns = shared.n, shared.ns
    began = time.perf_counter()
    counts = Counter()
    trace = []
    candidates = []
    status = "OK"
    stop_reason = "iteration_limit"
    stream = torch.cuda.Stream(device=shared.device)
    stream.wait_event(shared.ready)
    wait_seconds = 0.

    def pull(values, material=None):
        """One completion event waits only for this private stream."""
        nonlocal wait_seconds
        packed = torch.stack([v.real.reshape(()) for v in values]).to(torch.float64)
        host = torch.empty(len(values), dtype=torch.float64, pin_memory=True)
        host.copy_(packed, non_blocking=True)
        material_host = None
        if material is not None:
            material_host = torch.empty(n, dtype=torch.float64, pin_memory=True)
            material_host.copy_(material, non_blocking=True)
        done = torch.cuda.Event()
        done.record(stream)
        wait_start = time.perf_counter()
        done.synchronize()
        wait_seconds += time.perf_counter() - wait_start
        counts["host_result_returns"] += 1
        return host.numpy().copy(), None if material_host is None else material_host.numpy().copy()

    with torch.cuda.stream(stream), torch.no_grad():
        Y = torch.as_tensor(target, dtype=torch.complex128, device=shared.device)
        if tuple(Y.shape) != (shared.nr, ns):
            raise ValueError("target shape does not match resident operators")
        m = torch.as_tensor(m0, dtype=torch.float64, device=shared.device).reshape(-1).clone()
        if m.numel() != n:
            raise ValueError("m0 length does not match resident operators")
        m.clamp_(lower, upper)
        H = shared.G + beta * shared.eye
        y_scale = Y.abs().square().sum().real
        r = m / (1 + beta * m)
        state_scale = (r[:, None] * shared.E0).abs().square().sum().real
        scales, _ = pull([y_scale, state_scale])
        if scales[0] <= 0 or scales[1] <= 0:
            raise ValueError("fixed data and state scales must be positive")
        SHY = shared.SH @ Y / y_scale
        SHS = shared.SH @ shared.S / y_scale
        counts.update({"shared_gram_build": 1, "shared_target_product": 1,
                       "feedback_operator_build": 1})

        def physical(material, iteration):
            """One paid original Maxwell forward plus the common R4 prior."""
            L = shared.eye - material[:, None] * shared.G
            lu = torch.linalg.lu_factor_ex(L, check_errors=False)
            Jphys = torch.linalg.lu_solve(lu.LU, lu.pivots, material[:, None] * shared.E0)
            residual = shared.S @ Jphys - Y
            physical_data = .5 * residual.abs().square().sum().real / y_scale
            edge = (material[shared.right] - material[shared.left]) / 19.5
            smooth = float(cfg.get("smooth_weight", 1e-4))
            tv = float(cfg.get("tv_weight", 1e-4))
            tv_eps = float(cfg.get("tv_epsilon", .01))
            physical_reg = smooth * edge.square().mean() + tv * (torch.sqrt(edge.square() + tv_eps ** 2) - tv_eps).mean()
            counts.update({"physical_LU": 1, "physical_source_RHS": ns,
                           "physical_forward": 1, "physical_regularizer": 1})
            packed, saved_m = pull([lu.info.to(torch.float64), physical_data, physical_reg], material)
            info, data, reg = int(packed[0]), float(packed[1]), float(packed[2])
            if info != 0 or not np.isfinite(data + reg):
                return None, Jphys
            row = {"iteration": iteration, "m": saved_m.tolist(),
                   "physical_data_objective": data, "physical_regularizer": reg,
                   "physical_ranking_objective": data + reg, "lu_info": info,
                   "ranking_scope": "original full Maxwell data plus common R4 spatial regularizer"}
            candidates.append(row)
            return row, Jphys

        initial, J = physical(m, 0)
        if initial is None:
            return {"status": "FAIL", "stop_reason": "initial_physical_LU_failed",
                    "beta": beta, "m": None, "material_candidates": [],
                    "physical_incumbent": None, "currentsummary": None,
                    "trace": [], "counts": dict(counts), "wall_s": time.perf_counter() - began,
                    "host_wait_s": wait_seconds, "joint_objective_scope": "proposal_only"}
        incumbent = initial
        if time.perf_counter() - began >= wall_s:
            stop_reason = "wall_after_initial_forward"

        for iteration in range(1, max_iter + 1):
            if time.perf_counter() - began >= wall_s:
                stop_reason = "wall_limit"
                break
            A = shared.eye - r[:, None] * H
            rhs_state = r[:, None] * shared.E0
            normal = SHS + A.mH @ A / state_scale
            normal = .5 * (normal + normal.mH)
            rhs = SHY + A.mH @ rhs_state / state_scale
            chol = torch.linalg.cholesky_ex(normal, check_errors=False)
            Jnew = torch.cholesky_solve(rhs, chol.L)
            V = shared.E0 + H @ Jnew
            numerator = (V.conj() * Jnew).sum(dim=1).real
            denominator = V.abs().square().sum(dim=1).real
            safe_denominator = torch.where(denominator > 0, denominator, torch.ones_like(denominator))
            r_new = torch.where(denominator > 0, numerator / safe_denominator, r)
            r_new = r_new.clamp(r_lower, r_upper)
            m_new = r_new / (1 - beta * r_new)
            data_residual = shared.S @ Jnew - Y
            state_residual = Jnew - r_new[:, None] * V
            data_objective = .5 * data_residual.abs().square().sum().real / y_scale
            state_objective = .5 * state_residual.abs().square().sum().real / state_scale
            material_change = torch.linalg.vector_norm(m_new - m)
            counts.update({"joint_cholesky": 1, "joint_source_RHS": ns,
                           "joint_G_products": 1, "joint_state_gram": 1,
                           "joint_objective_evaluations": 1,
                           "joint_material_row_updates": 1})
            is_checkpoint = iteration % checkpoint_every == 0 or iteration == max_iter
            values, _ = pull([chol.info.to(torch.float64), data_objective,
                              state_objective, material_change])
            info = int(values[0])
            if info != 0 or not np.all(np.isfinite(values[1:])):
                status = "FAIL"
                stop_reason = "joint_cholesky_failed" if info else "joint_nonfinite"
                trace.append({"iteration": iteration, "status": stop_reason,
                              "cholesky_info": info})
                break
            J, r, m = Jnew, r_new, m_new
            trace.append({"iteration": iteration, "auxiliary_data": float(values[1]),
                          "auxiliary_state": float(values[2]),
                          "auxiliary_total": float(values[1] + values[2]),
                          "material_change_norm": float(values[3]),
                          "physical_checkpoint": bool(is_checkpoint)})
            if is_checkpoint:
                candidate, _ = physical(m, iteration)
                if candidate is None:
                    status = "FAIL"
                    stop_reason = "checkpoint_physical_LU_failed"
                    break
                if candidate["physical_ranking_objective"] < incumbent["physical_ranking_objective"]:
                    incumbent = candidate

        if status == "OK" and trace and trace[-1]["iteration"] != candidates[-1]["iteration"]:
            # A wall stop can leave a valid auxiliary material unscreened.
            if time.perf_counter() - began < wall_s:
                candidate, _ = physical(m, trace[-1]["iteration"])
                if candidate is not None and candidate["physical_ranking_objective"] < incumbent["physical_ranking_objective"]:
                    incumbent = candidate
            else:
                stop_reason = "wall_limit_unscreened_final_auxiliary"

        final_data_residual = shared.S @ J - Y
        final_state_residual = J - r[:, None] * (shared.E0 + H @ J)
        summary_values, _ = pull([torch.linalg.vector_norm(J),
                                  torch.linalg.vector_norm(final_data_residual),
                                  torch.linalg.vector_norm(final_state_residual),
                                  torch.sqrt(state_scale)])

    return {"status": status, "stop_reason": stop_reason, "beta": beta,
            "m": incumbent["m"], "material_candidates": candidates,
            "physical_incumbent": incumbent,
            "currentsummary": {"final_current_norm": float(summary_values[0]),
                               "final_auxiliary_data_residual_norm": float(summary_values[1]),
                               "final_auxiliary_state_residual_norm": float(summary_values[2]),
                               "fixed_state_scale_norm": float(summary_values[3]),
                               "fixed_data_scale_squared": float(scales[0]),
                               "fixed_state_scale_squared": float(scales[1])},
            "trace": trace, "counts": dict(counts),
            "wall_s": time.perf_counter() - began, "host_wait_s": wait_seconds,
            "joint_objective_scope": "exploratory_classical_CSI_CIE_proposal_only",
            "material_step_regularizer": "OMITTED_proposal_only",
            "physical_selection_scope": "original_full_Maxwell_data_plus_common_R4_regularizer"}
