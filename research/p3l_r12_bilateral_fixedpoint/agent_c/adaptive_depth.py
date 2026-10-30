"""Prediction-only diagnostics for hierarchical exploration. No physics or GT IO.

This module deliberately accepts arrays and budget summaries rather than a
Problem/Physics object. The shared engine remains the only physical caller.
"""
from dataclasses import dataclass, asdict
import numpy as np


def orthobasis(A, rtol=1e-10):
    A = np.asarray(A)
    if A.ndim != 2:
        raise ValueError("basis must be a matrix")
    if not A.size:
        return A[:, :0]
    u, s, _ = np.linalg.svd(A, full_matrices=False)
    return u[:, s > rtol * max(float(s[0]), 1e-300)]


def outside_span(vector, basis):
    """Relative and absolute Euclidean distance from an explicitly named span."""
    v = np.asarray(vector).reshape(-1)
    q = orthobasis(basis)
    if len(v) != q.shape[0]:
        raise ValueError("span/vector coordinate mismatch")
    remainder = v - q @ (q.conj().T @ v)
    n = float(np.linalg.norm(v))
    return dict(absolute=float(np.linalg.norm(remainder)),
                relative=0.0 if n == 0 else float(np.linalg.norm(remainder) / n),
                basis_rank=q.shape[1], coordinate_metric="Euclidean in supplied coordinates")


def projector_change(old_basis, new_basis):
    """Subspace change is independent of singular-vector signs/rotations."""
    q0, q1 = orthobasis(old_basis), orthobasis(new_basis)
    if q0.shape[0] != q1.shape[0]:
        raise ValueError("projector coordinate mismatch")
    overlap = q0.conj().T @ q1
    novelty2 = max(0.0, q1.shape[1] - float(np.linalg.norm(overlap, 'fro') ** 2))
    fro2 = max(0.0, q0.shape[1] + q1.shape[1] - 2 * float(np.linalg.norm(overlap, 'fro') ** 2))
    return dict(old_rank=q0.shape[1], new_rank=q1.shape[1],
                new_outside_old_Frobenius=float(np.sqrt(novelty2)),
                new_outside_old_fraction=float(np.sqrt(novelty2 / max(q1.shape[1], 1))),
                projector_Frobenius=float(np.sqrt(fro2)),
                principal_cosines=np.linalg.svd(overlap, compute_uv=False).tolist())


def finite_secant(x0, x1, response0, response1, jacobian0, jacobian1=None):
    """All response/J arguments must use the SAME original whitening/packing.

    A large secant defect is local linearization evidence. It is not a
    statistical information gain, a coverage certificate or a derivative bound.
    """
    dx = np.asarray(x1) - np.asarray(x0)
    df = np.asarray(response1).reshape(-1) - np.asarray(response0).reshape(-1)
    j0 = np.asarray(jacobian0)
    if j0.shape != (len(df), len(dx)):
        raise ValueError("secant response/J coordinates do not agree")
    tangent = j0 @ dx
    defect = df - tangent
    out = dict(step_norm=float(np.linalg.norm(dx)),
               response_change_norm=float(np.linalg.norm(df)),
               tangent_change_norm=float(np.linalg.norm(tangent)),
               secant_defect_norm=float(np.linalg.norm(defect)),
               secant_defect_relative=float(np.linalg.norm(defect) / max(np.linalg.norm(df), 1.0)),
               normalization="max(norm(actual whitened response change), 1)")
    if jacobian1 is not None:
        j1 = np.asarray(jacobian1)
        if j1.shape != j0.shape:
            raise ValueError("refreshed J coordinate mismatch")
        out.update(jacobian_relative_change=float(np.linalg.norm(j1-j0, 'fro') / max(np.linalg.norm(j0, 'fro'), 1e-30)),
                   refreshed_secant_defect_norm=float(np.linalg.norm(df - .5*(j0+j1)@dx)))
    return out


@dataclass(frozen=True)
class DepthPolicy:
    # Deterministic defaults frozen before R12 scoring. Engine may impose a
    # smaller max_depth/budget. Any changed defaults require a new DEV version.
    max_depth: int = 3
    min_depth: int = 1
    saturation_patience: int = 2
    distinct_rms: float = 1e-6
    task_mean_increment: float = 1e-4
    task_energy_increment: float = 1e-4
    task_count_increment: int = 1
    secant_refresh_threshold: float = 0.10
    projector_novelty_threshold: float = 0.05

    def decide(self, completed_depth, history, remaining=None, next_atomic_cost=None):
        """History contains ONLINE completed-layer summaries only.

        Required history entries: distinct_added, mean_range_increment,
        energy_range_increment, count_range_increment. Optional secant and
        projector fields affect the explanation, never override a hard cap.
        Missing task evidence is not interpreted as zero/saturation.
        """
        if completed_depth >= self.max_depth:
            return dict(action="STOP", reason="MAX_DEPTH", policy=asdict(self))
        if remaining is not None and next_atomic_cost is not None:
            if any(int(next_atomic_cost.get(k, 0)) > int(remaining.get(k, 0))
                   for k in ("LU", "RHS", "RHS576")):
                return dict(action="STOP", reason="ATOMIC_BUDGET_INSUFFICIENT", policy=asdict(self))
        if completed_depth < self.min_depth:
            return dict(action="DEEPEN", reason="MINIMUM_DECLARED_DEPTH", policy=asdict(self))
        recent = list(history)[-self.saturation_patience:]
        fields = ("distinct_added", "mean_range_increment", "energy_range_increment", "count_range_increment")
        complete = len(recent) == self.saturation_patience and all(all(k in row for k in fields) for row in recent)
        def saturated(r):
            return (int(r["distinct_added"]) == 0
                    and float(r["mean_range_increment"]) <= self.task_mean_increment
                    and float(r["energy_range_increment"]) <= self.task_energy_increment
                    and int(r["count_range_increment"]) < self.task_count_increment)
        if complete and all(saturated(r) for r in recent):
            return dict(action="STOP", reason="TWO_LAYER_ONLINE_SATURATION", policy=asdict(self))
        last = recent[-1] if recent else {}
        refresh = (float(last.get("secant_defect_relative", 0)) > self.secant_refresh_threshold
                   or float(last.get("new_outside_old_fraction", 0)) > self.projector_novelty_threshold)
        return dict(action="DEEPEN", reason="REFRESHED_LOCAL_MODEL_EVIDENCE" if refresh else "UNSATURATED_OR_INCOMPLETE_ONLINE_EVIDENCE",
                    refresh_recommended=refresh, policy=asdict(self))
