"""Material risk actions from recorded trajectories, with pure GT-free inference.

Convex coefficients are estimator coefficients, never posterior probabilities.
No solver, job, truth-file reader, or training action is imported here.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

REGULARIZATION = 1e-4
METHODS = ("learned_convex", "arc_length_trapezoid", "equal_steps", "last_point", "anchor")


def _simplex(weights, nodes):
    w = np.asarray(weights, float)
    if w.shape != (nodes,) or not np.isfinite(w).all() or np.min(w) < -1e-10 or abs(w.sum()-1) > 1e-8:
        raise ValueError("finite nonnegative estimator coefficients summing to one required")
    w = np.maximum(w, 0.)
    return w/w.sum()


def integrate_materials(materials, *, method="learned_convex", weights=None, expected_nodes=9):
    """Pure online action on anchor+ordered endpoints; no GT inputs or access.

    Incomplete, nonfinite endpoint paths return the original paid anchor and a
    flag. An absent/nonfinite anchor is an error because no fallback exists.
    Equal-steps means equal coefficients over all nine recorded nodes, anchor
    included. Arc-length trapezoids use actual complex material RMS lengths.
    """
    if method not in METHODS:
        raise ValueError("unknown trajectory risk action")
    stack = np.asarray(materials, complex)
    if stack.ndim != 2 or stack.shape[0] < 1 or stack.shape[1] != 576 or not np.isfinite(stack[0]).all():
        raise ValueError("finite original anchor of 576 material cells required")
    anchor = stack[0].copy()
    complete = stack.shape == (int(expected_nodes), 576) and np.isfinite(stack).all()
    if not complete:
        return anchor, dict(method=method, fallback_to_anchor=True, fallback_reason="INCOMPLETE_OR_INVALID_PATH",
            expected_nodes=int(expected_nodes), received_nodes=len(stack), estimator_coefficients=[1.]+[0.]*(int(expected_nodes)-1),
            posterior_probabilities=False)
    w = np.zeros(int(expected_nodes))
    fallback, reason = False, None
    if method == "anchor":
        w[0] = 1.
    elif method == "last_point":
        w[-1] = 1.
    elif method == "equal_steps":
        w[:] = 1/len(w)
    elif method == "arc_length_trapezoid":
        lengths = np.linalg.norm(np.diff(stack, axis=0), axis=1)/np.sqrt(stack.shape[1])
        if lengths.sum() <= 1e-14:
            w[0], fallback, reason = 1., True, "ZERO_ARC_LENGTH"
        else:
            w[:-1] += lengths/2
            w[1:] += lengths/2
            w /= lengths.sum()
    elif weights is None:
        w[0], fallback, reason = 1., True, "NO_FROZEN_COEFFICIENTS_FOR_ARM"
    else:
        w = _simplex(weights, int(expected_nodes))
    estimate = w@stack
    return estimate, dict(method=method, fallback_to_anchor=fallback, fallback_reason=reason,
        expected_nodes=int(expected_nodes), received_nodes=len(stack), estimator_coefficients=w.tolist(),
        posterior_probabilities=False)


def fit_convex_weights(material_stacks, truths, *, regularization=REGULARIZATION):
    """Offline train-only full-grid complex MSE plus anchor regularization.

    Objective = mean_parent mean_cell |sum_node w chi - chi_GT|^2
                + 1e-4 ||w-e_anchor||_2^2, with w>=0 and sum(w)=1.
    The material MSE is neither rescaled nor normalized by anchor risk.
    """
    if float(regularization) != REGULARIZATION:
        raise ValueError("regularization frozen at 1e-4 in material-MSE units")
    stacks, truth = np.asarray(material_stacks, complex), np.asarray(truths, complex)
    if stacks.ndim != 3 or stacks.shape[2] != 576 or truth.shape != (len(stacks), 576):
        raise ValueError("train stacks parent x node x 576 and train truth parent x 576 required")
    if not np.isfinite(stacks).all() or not np.isfinite(truth).all():
        raise ValueError("finite training arrays required")
    parents, nodes, cells = stacks.shape
    e0 = np.zeros(nodes); e0[0] = 1.
    differences = stacks-stacks[:, :1]
    residual0 = stacks[:, 0]-truth
    gram = np.real(np.einsum("pnc,pmc->nm", differences.conj(), differences))/(parents*cells)
    linear = np.real(np.einsum("pnc,pc->n", differences.conj(), residual0))/(parents*cells)
    base = float(np.mean(abs(residual0)**2))
    def objective(w):
        return float(base+2*linear@w+w@gram@w+regularization*np.dot(w-e0, w-e0))
    def gradient(w):
        return 2*(linear+gram@w+regularization*(w-e0))
    fit = minimize(objective, e0, jac=gradient, method="SLSQP",
        bounds=[(0., 1.)]*nodes,
        constraints=[dict(type="eq", fun=lambda w:float(w.sum()-1), jac=lambda w:np.ones(nodes))],
        options=dict(ftol=1e-13, maxiter=1000))
    if not fit.success:
        raise RuntimeError(f"convex weight fit failed: {fit.message}")
    w = _simplex(fit.x, nodes)
    estimates = np.einsum("n,pnc->pc", w, stacks)
    mse = float(np.mean(abs(estimates-truth)**2))
    penalty = float(regularization*np.dot(w-e0, w-e0))
    # Check the returned point against the feasible anchor, preserving failure.
    if mse+penalty > base+1e-11:
        raise RuntimeError("fitted action worse than feasible anchor objective")
    active = w > 1e-8
    grad = gradient(w)
    multiplier = float(np.mean(grad[active]))
    residual = max(float(np.max(abs(grad[active]-multiplier))),
        max(0., float(np.max(multiplier-grad[~active]))) if (~active).any() else 0.)
    return w, dict(optimizer_success=True, optimizer_message=str(fit.message), iterations=int(fit.nit),
        train_parents=parents, nodes=nodes, train_full_grid_complex_MSE=mse,
        anchor_train_full_grid_complex_MSE=base, regularization=regularization,
        regularization_penalty=penalty, penalized_objective=mse+penalty,
        simplex_KKT_residual=residual, coefficient_sum=float(w.sum()), minimum_coefficient=float(w.min()),
        risk_action=True, posterior_integration=False)


def load_frozen(path):
    frozen = json.loads(Path(path).read_text())
    allowed = {"schema_version", "kind", "training_split", "train_parents", "nodes", "arms",
        "regularization", "objective", "frozen_before_fresh_evaluation", "historical_train_exposure",
        "coefficient_semantics", "baselines", "limitations", "source_online_job"}
    if set(frozen)-allowed or frozen.get("kind") != "constant_convex_trajectory_risk_action":
        raise ValueError("unexpected frozen risk-action schema")
    if frozen.get("training_split") != "train" or frozen.get("regularization") != REGULARIZATION:
        raise ValueError("train-only frozen regularization provenance required")
    if not frozen.get("frozen_before_fresh_evaluation"):
        raise ValueError("action must be frozen before fresh evaluation")
    for arm, coefficients in frozen["arms"].items():
        _simplex(coefficients, int(frozen["nodes"]))
    return frozen


def apply_frozen(materials, frozen, arm, *, method="learned_convex"):
    """Pure estimator; missing learned arm and incomplete path fall back visibly."""
    weights = frozen["arms"].get(arm)
    return integrate_materials(materials, method=method, weights=weights, expected_nodes=int(frozen["nodes"]))
