"""Development-only Q1 CIE mechanism arms over the unchanged CUDA core.

Run ``python -m research.research_fbg_r2_adaptive_feedback.experiments.cie_mechanism_probe MANIFEST``.
The CLI imports ``src.adaptive_solver``, patches only ``solve`` for explicitly
named methods, and delegates every other task unchanged. ``cie_coordinate``
is exactly the original physical data objective in a rational coordinate;
``joint_cie:0`` and ``joint_cie:0.2`` are exploratory classical CSI/CIE
joint-current proposals followed by an original-physics Native TRF corrector.
Neither arm is a feedback-Gram objective or a new CIE theory.
"""
from collections import Counter
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize

PROJECT = Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from src import adaptive_solver as ad
from src.joint_cie import propose as joint_propose


def _coordinate_phase(lane, x0, stage, cfg, wall_s, beta, phase):
    """One bounded original-data L-BFGS-B phase in normalized rho coordinates.

    c is the standard real coarse coefficient. rho=c/(1+beta*c) is applied
    coefficientwise before the existing bilinear basis, so every legal c maps
    to exactly the same legal physical material as the untransformed solver.
    For z=(rho-rho_lo)/(rho_hi-rho_lo), dc/dz=(rho_hi-rho_lo)
    (1+beta*c)^2. ``lane.evaluate`` supplies dF/dc including the common
    regularizer; this phase applies that exact gradient chain.
    """
    lower, upper = float(cfg.get("lower", .02)), float(cfg.get("upper", 20.))
    if beta < 0 or 1 + beta * lower <= 0 or 1 + beta * upper <= 0:
        raise ValueError("coordinate beta is not monotone on the legal box")
    rho_lo, rho_hi = lower / (1 + beta * lower), upper / (1 + beta * upper)
    rho_width = rho_hi - rho_lo
    start = time.perf_counter()
    x0 = np.clip(np.asarray(x0, dtype=np.float64), lower, upper)
    z0 = (x0 / (1 + beta * x0) - rho_lo) / rho_width
    best_coefficients, best_value = x0.copy(), float("inf")
    trace = []

    def callback(z):
        nonlocal best_coefficients, best_value
        if time.perf_counter() - start >= wall_s:
            raise ad.engine.BudgetStop("WALL_SECONDS")
        rho = rho_lo + rho_width * z
        coefficients = rho / (1 - beta * rho)
        value, gradient_c = lane.evaluate(coefficients)
        physical = float(lane.last_physical_objective)
        if physical < best_value:
            best_coefficients, best_value = coefficients.copy(), physical
        trace.append({"evaluation": len(trace) + 1,
                      "seconds": time.perf_counter() - start,
                      "objective": float(value), "physical_objective": physical})
        gradient_z = np.asarray(gradient_c, dtype=np.float64) * (
            rho_width * (1 + beta * coefficients) ** 2)
        return float(value), gradient_z

    fit = None
    status, stop_reason = "COMPLETED_NUMERICAL", None
    try:
        fit = minimize(callback, np.clip(z0, 0., 1.), jac=True, method="L-BFGS-B",
                       bounds=[(0., 1.)] * len(z0),
                       options={"maxiter": int(cfg.get("maxiter", 2000)),
                                "maxfun": int(cfg.get("maxfun", 4000)),
                                "ftol": float(cfg.get("ftol", 1e-12)),
                                "gtol": float(cfg.get("gtol", 1e-8)),
                                "maxls": int(cfg.get("maxls", 25))})
    except ad.engine.BudgetStop as exc:
        status, stop_reason = "BUDGET_STOP", str(exc)
    result = lane.close(best_coefficients)
    if not np.isfinite(best_value):
        best_value = float(result["ranking_objective"])
    result.update({"wrapper_method": "cie_coordinate", "phase": phase,
                   "coordinate_beta": beta,
                   "coordinate_transform": "rho=c/(1+beta*c), coefficientwise before unchanged physical basis",
                   "gradient_chain": "dF/dz=(dF/dc)*(rho_hi-rho_lo)*(1+beta*c)^2",
                   "objective_scope": "original full Maxwell data plus common regularizer",
                   "optimizer_status": None if fit is None else int(fit.status),
                   "optimizer_success": None if fit is None else bool(fit.success),
                   "optimizer_message": None if fit is None else str(fit.message),
                   "optimizer_nfev": None if fit is None else int(fit.nfev),
                   "wrapper_status": status, "wrapper_stop_reason": stop_reason,
                   "optimizer_wall_s": time.perf_counter() - start,
                   "best_physical_objective": best_value,
                   "evaluation_trace": trace, "_coefficients": best_coefficients,
                   "stage": stage, "truth_access": "NONE_PUBLIC_INPUT_ONLY"})
    return result


def solve_coordinate(shared, public, task, context=None):
    """The standard 8/12/24 schedule with two original-data 55/45 phases."""
    cfg = dict(task.get("config", {}))
    beta = float(cfg.get("beta", .2))
    wall = float(cfg.get("wall_seconds", 45.))
    started = time.perf_counter()
    fractions = cfg.get("stage_fractions", [.2, .25, .55])
    stages = cfg.get("stages", [8, 12, 24])
    if list(map(int, stages)) != [8, 12, 24] or len(fractions) != 3:
        raise ValueError("coordinate probe requires the common 8/12/24 stage schedule")
    coefficients, old_stage = None, None
    phases = []
    for stage, fraction in zip(stages, fractions):
        stage = int(stage)
        remaining = wall - (time.perf_counter() - started)
        if remaining <= 0:
            break
        c = {**cfg, "stage": stage, "wall_seconds": min(remaining, wall * fraction)}
        if coefficients is not None:
            x = ad.engine._resample(coefficients, old_stage, stage)
        elif "initial_pixels" in task:
            x = ad.engine._resample(task["initial_pixels"], 24, stage)
        else:
            x = np.full(stage * stage, task.get("initial_material", .1), dtype=np.float64)
        lane = ad.FeedbackLane(shared, public, stage=stage, objective="data", config=c,
                               task_context=context)
        if context is not None:
            context["lanes"].append(lane)
        stage_left = c["wall_seconds"]
        for index, fraction_of_stage in enumerate((.55, .45)):
            remaining = min(wall - (time.perf_counter() - started), stage_left)
            if remaining <= 0:
                break
            phase_wall = remaining * fraction_of_stage if index == 0 else remaining
            previous = Counter(lane.counts)
            pstart = time.perf_counter()
            result = _coordinate_phase(lane, x, stage, c, phase_wall, beta,
                                       "proposal" if index == 0 else "corrector")
            x = result.pop("_coefficients")
            stage_left -= time.perf_counter() - pstart
            result["counts"] = dict(Counter(result["counts"]) - previous)
            phases.append(result)
        coefficients, old_stage = x, stage
    if not phases:
        raise RuntimeError("no coordinate physical phase completed")
    chosen = min(phases, key=lambda row: row["ranking_objective"])
    counts = Counter()
    for row in phases:
        counts.update(row["counts"])
    return {**chosen, "method": task["method"], "counts": dict(counts),
            "phases": phases, "task_wall_s": time.perf_counter() - started,
            "stages_reached": sorted({row["stage"] for row in phases}),
            "coordinate_beta": beta, "feedback_operator": "NONE",
            "online_selection": "minimum original full Maxwell data plus common regularizer over paid candidates",
            "truth_access": "NONE", "budget": cfg}


def solve_joint(shared, public, task, context=None):
    """Classical joint-current proposal, then paid 24-pixel Native TRF."""
    cfg = dict(task.get("config", {}))
    method = task["method"]
    beta = float(method.split(":", 1)[1])
    wall = float(cfg.get("wall_seconds", 45.))
    started = time.perf_counter()
    if "initial_pixels" in task:
        m0 = np.asarray(task["initial_pixels"], dtype=np.float64)
    else:
        m0 = np.full(shared.n, task.get("initial_material", .1), dtype=np.float64)
    if m0.shape != (shared.n,):
        raise ValueError("joint CIE probe requires a full-resolution initial material")
    Y = public.get("Y", public.get("y"))
    joint_cfg = {**cfg, "beta": beta,
                 "material_lower": float(cfg.get("lower", .02)),
                 "material_upper": float(cfg.get("upper", 20.))}
    # Bound the proposal by the task's total physical-LU/full-size-RHS caps.
    # A later Native TRF can consume the remaining allowance; no joint work
    # receives a free source block merely because it uses Cholesky rather than LU.
    proposed_iterations = min(100, max(0, int(joint_cfg.get("max_iter", 100))))
    checkpoint_every = int(joint_cfg.get("checkpoint_every", 10))
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")
    max_lu = float(cfg.get("max_lu", float("inf")))
    max_rhs = float(cfg.get("max_rhs", float("inf")))
    if max_lu <= 1 or max_rhs <= shared.ns:
        raise ad.engine.BudgetStop("joint proposal lacks allowance for its initial paid full-Maxwell forward")
    while proposed_iterations > 0:
        # One possible final unscreened-material check after an early wall stop.
        checkpoints = (proposed_iterations + checkpoint_every - 1) // checkpoint_every + 1
        if 1 + checkpoints < max_lu and shared.ns * (1 + proposed_iterations + checkpoints) < max_rhs:
            break
        proposed_iterations -= 1
    joint_cfg["max_iter"] = proposed_iterations
    proposal = joint_propose(shared, Y, m0, .45 * wall, joint_cfg)
    if proposal.get("physical_incumbent") is None:
        raise RuntimeError("joint-current proposal has no physically screened material")
    spent = time.perf_counter() - started
    remaining = min(.55 * wall, max(0., wall - spent))
    lane = ad.FeedbackLane(shared, public, stage=24, objective="data", config=cfg,
                           task_context=context)
    if context is not None:
        context["lanes"].append(lane)
    proposal_counts = Counter(proposal["counts"])
    # Charge all full-size current solves and paid Maxwell checkpoints to the
    # same physical RHS cap before the corrector starts.
    lane.counts.update({"LU": proposal_counts["physical_LU"],
                        "source_RHS": proposal_counts["physical_source_RHS"] + proposal_counts["joint_source_RHS"],
                        "joint_cholesky": proposal_counts["joint_cholesky"],
                        "joint_source_RHS": proposal_counts["joint_source_RHS"],
                        "joint_physical_LU": proposal_counts["physical_LU"],
                        "joint_physical_source_RHS": proposal_counts["physical_source_RHS"]})
    charged_proposal_counts = Counter(lane.counts)
    corrector = None
    if remaining > 0:
        corrector_cfg = {**cfg, "stage": 24, "wall_seconds": remaining,
                         "initial_coefficients": proposal["m"]}
        corrector = ad.engine.solve_task(lane, {"method": "native_trf", "stage": 24,
                                                "config": corrector_cfg})
        corrector.pop("_coefficients", None)
        corrector["phase"] = "original_physics_native_trf_corrector"
        corrector["counts"] = dict(Counter(corrector["counts"]) - charged_proposal_counts)
    proposal_accounting = Counter(proposal_counts)
    proposal_accounting.update({"LU": proposal_counts["physical_LU"],
                                "source_RHS": proposal_counts["physical_source_RHS"] + proposal_counts["joint_source_RHS"]})
    proposal_row = {"phase": "exploratory_classical_joint_current_proposal",
                    "ranking_objective": proposal["physical_incumbent"]["physical_ranking_objective"],
                    "relative_residual": float(np.sqrt(2 * proposal["physical_incumbent"]["physical_data_objective"])),
                    "m": proposal["m"],
                    "counts": dict(proposal_accounting), "receipt": proposal,
                    "stage": 24}
    phases = [proposal_row] + ([corrector] if corrector is not None else [])
    chosen = min(phases, key=lambda row: row["ranking_objective"])
    total_counts = Counter(proposal_accounting)
    if corrector is not None:
        total_counts.update(corrector["counts"])
    return {**chosen, "method": method, "counts": dict(total_counts),
            "phases": phases, "joint_proposal": proposal,
            "task_wall_s": time.perf_counter() - started,
            "stages_reached": [24], "feedback_operator": "NONE",
            "schedule_scope": "45 percent joint-current proposal then 55 percent 24-pixel Native TRF; differs from standard multiscale schedule",
            "online_selection": "minimum original full Maxwell data plus common regularizer over paid candidates",
            "truth_access": "NONE", "budget": cfg}


def main():
    original = ad.solve

    def patched(shared, public, task, context=None):
        if task["method"] == "cie_coordinate":
            return solve_coordinate(shared, public, task, context)
        if task["method"] in ("joint_cie:0", "joint_cie:0.2"):
            return solve_joint(shared, public, task, context)
        return original(shared, public, task, context)

    ad.solve = patched
    try:
        ad.main()
    finally:
        ad.solve = original


if __name__ == "__main__":
    main()
