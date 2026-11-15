"""Dense fixed-background field/current projections, with no truth interface.

Rows are cell-major: scalar ``(N, sources)`` or isotropic vector
``(3*N, sources)`` with consecutive x/y/z rows belonging to one cell.
The physical backend and the data-plus-prior objective remain caller-owned.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence
from copy import deepcopy
import json
import time

import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.optimize import brentq


class UnreachableDataBall(ValueError):
    """The fixed receiver map cannot reach the requested global data ball."""

    def __init__(self, delta: float, minimum_residual: float):
        self.delta, self.minimum_residual = delta, minimum_residual
        super().__init__(f"delta={delta:.16g} is below reachable residual "
                         f"{minimum_residual:.16g}")


@dataclass
class PhysicalAudit:
    """One physical solve supplied by the existing backend.

    ``objective`` includes precisely the caller's common data and prior terms.
    ``costs`` contains increments for THIS call, e.g. physical_lu/physical_rhs.
    J is unscaled physical current; Q is the physical total field.
    """
    objective: float
    J: np.ndarray
    Q: np.ndarray
    metrics: Mapping[str, Any] = field(default_factory=dict)
    costs: Mapping[str, float] = field(default_factory=dict)


@dataclass
class AuditRecord:
    iteration: int
    objective: float
    incumbent_objective: float
    incumbent_iteration: int
    cpu_seconds: float
    wall_seconds: float
    gap: float
    candidate_material: np.ndarray
    incumbent_material: np.ndarray
    # Complete trajectory state, both shadows, and physical audited fields.
    states: dict[str, np.ndarray]
    metrics: dict[str, Any]
    incumbent_metrics: dict[str, Any]
    costs: dict[str, float]
    physical_costs: dict[str, float]


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(v) for v in value]
    return value


@dataclass
class ProjectionRun:
    method: str
    material: np.ndarray
    objective: float
    incumbent_iteration: int
    iterations: int
    stop_reason: str
    history: list[AuditRecord]
    last_X: np.ndarray
    last_Q: np.ndarray
    costs: dict[str, float]
    physical_costs: dict[str, float]
    setup_cpu_seconds: float
    setup_wall_seconds: float
    cpu_seconds: float
    wall_seconds: float

    def save(self, directory: str | Path) -> None:
        """Save all audits, incumbents and states without pickle or truth."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        arrays = {"material": self.material, "last_X": self.last_X,
                  "last_Q": self.last_Q}
        records = []
        for k, record in enumerate(self.history):
            prefix = f"audit_{k:05d}_"
            arrays[prefix + "candidate_material"] = record.candidate_material
            arrays[prefix + "incumbent_material"] = record.incumbent_material
            arrays.update({prefix + key: value for key, value in record.states.items()})
            records.append({key: value for key, value in vars(record).items()
                            if key not in {"candidate_material", "incumbent_material", "states"}})
        np.savez_compressed(directory / "states.npz", **arrays)
        summary = {key: value for key, value in vars(self).items()
                   if key not in {"material", "last_X", "last_Q", "history"}}
        summary["history"] = records
        (directory / "history.json").write_text(
            json.dumps(_jsonable(summary), ensure_ascii=False, indent=2) + "\n")


@dataclass
class _Port:
    T: np.ndarray
    V: np.ndarray
    eigenvalues: np.ndarray
    U: np.ndarray
    positive: np.ndarray


class DenseFieldCurrentProjectors:
    """Exact dense A/C projectors in the fixed metric ||J/scale||²+||Q||².

    S may be a common matrix or one matrix per source. Y may likewise be a
    receiver-by-source matrix or a list of source vectors. Receiver masks act
    before whitening. A whitener is None, a scalar, a diagonal vector or a
    matrix; supply one item per source in ``whiteners``. These are independent
    source blocks, all constrained by ONE delta and ONE multiplier.
    """

    def __init__(self, G: np.ndarray, E: np.ndarray,
                 S: np.ndarray | Sequence[np.ndarray],
                 Y: np.ndarray | Sequence[np.ndarray], delta: float,
                 *, scale: float = 6., lower: float | np.ndarray = 0.,
                 upper: float | np.ndarray = 12., components: int = 1,
                 receiver_masks: Sequence[np.ndarray | None] | None = None,
                 whiteners: Sequence[np.ndarray | float | None] | None = None,
                 public_metadata: Mapping[str, Any] | None = None,
                 grouped_ports: bool = False):
        wall, cpu = time.perf_counter(), time.process_time()
        self.G = np.array(G, dtype=complex, copy=True)
        self.E = np.array(E, dtype=complex, copy=True)
        if (self.G.ndim != 2 or self.G.shape[0] != self.G.shape[1]
                or self.E.ndim != 2 or self.E.shape[0] != self.G.shape[0]):
            raise ValueError("G must be square and E must have shape (dofs, sources)")
        if components not in (1, 3) or len(self.G) % components:
            raise ValueError("components must be 1 or 3 with cell-major rows")
        if not np.isfinite(scale) or scale <= 0 or not np.isfinite(delta) or delta < 0:
            raise ValueError("scale must be positive and delta nonnegative and finite")
        if not np.all(np.isfinite(self.G)) or not np.all(np.isfinite(self.E)):
            raise ValueError("G/E must be finite")
        self.components, self.cells = components, len(self.G) // components
        self.sources, self.scale, self.delta = self.E.shape[1], float(scale), float(delta)
        self.lower = np.broadcast_to(np.asarray(lower, float), (self.cells,)).copy()
        self.upper = np.broadcast_to(np.asarray(upper, float), (self.cells,)).copy()
        if (np.any(~np.isfinite(self.lower)) or np.any(~np.isfinite(self.upper))
                or np.any(self.lower > self.upper)):
            raise ValueError("material box must be finite and ordered")
        self.public_metadata = dict(public_metadata or {})
        # Optional v2: default false preserves the initially frozen v1 path.
        self.grouped_ports = bool(grouped_ports)
        self.costs = {key: 0. for key in (
            "H_cholesky", "H_setup_rhs", "H_iteration_rhs", "port_eigh",
            "G_action_calls", "G_action_rhs", "G_adjoint_action_calls",
            "G_adjoint_action_rhs", "receiver_action_calls", "receiver_action_rhs",
            "receiver_adjoint_action_calls", "receiver_adjoint_action_rhs",
            "A_projections", "C_projections", "multiplier_evaluations")}
        self.B = self.scale * self.G
        self.H = cho_factor(np.eye(len(self.G)) + self.B.conj().T @ self.B,
                            lower=True, check_finite=False)
        self.costs["H_cholesky"] = 1.
        masks = list(receiver_masks) if receiver_masks is not None else [None] * self.sources
        white = list(whiteners) if whiteners is not None else [None] * self.sources
        if len(masks) != self.sources or len(white) != self.sources:
            raise ValueError("receiver_masks/whiteners must have one item per source")
        if isinstance(S, np.ndarray) and S.ndim == 2:
            source_S = [S] * self.sources
        else:
            source_S = list(S)
        if isinstance(Y, np.ndarray) and Y.ndim == 2:
            if Y.shape[1] != self.sources:
                raise ValueError("Y columns must match source count")
            source_Y = [Y[:, j] for j in range(self.sources)]
        else:
            source_Y = list(Y)
        if len(source_S) != self.sources or len(source_Y) != self.sources:
            raise ValueError("S and Y must have one block per source")
        self.ports, self.targets = [], []
        unique_ports: list[_Port] = []
        for matrix, target, mask, W in zip(source_S, source_Y, masks, white):
            matrix, target = np.asarray(matrix, complex), np.asarray(target, complex)
            if matrix.ndim != 2 or matrix.shape[1] != len(self.G) or target.shape != (len(matrix),):
                raise ValueError("each S/Y block must have shapes (receivers,dofs)/(receivers,)")
            if mask is not None:
                mask = np.asarray(mask)
                if mask.ndim != 1:
                    raise ValueError("receiver masks must be 1-D booleans or indices")
                matrix, target = matrix[mask], target[mask]
            T, yw = self.scale * matrix, target
            if W is not None:
                W = np.asarray(W)
                if W.ndim == 0:
                    T, yw = W * T, W * yw
                elif W.ndim == 1 and W.shape == (len(matrix),):
                    T, yw = W[:, None] * T, W * yw
                elif W.ndim == 2 and W.shape[1] == len(matrix):
                    T, yw = W @ T, W @ yw
                else:
                    raise ValueError("whitener shape incompatible with masked receivers")
            T, yw = np.array(T, complex, copy=True), np.array(yw, complex, copy=True)
            if not np.all(np.isfinite(T)) or not np.all(np.isfinite(yw)):
                raise ValueError("S/Y/whitening must be finite")
            port = next((p for p in unique_ports if np.array_equal(p.T, T)), None)
            if port is None:
                V = cho_solve(self.H, T.conj().T, check_finite=False)
                self.costs["H_setup_rhs"] += len(T)
                K = T @ V
                values, U = eigh((K + K.conj().T) / 2, check_finite=False)
                self.costs["port_eigh"] += 1
                eigen_scale = max(float(np.max(np.abs(values), initial=0.)), np.finfo(float).tiny)
                tolerance = 32 * np.finfo(float).eps * max(1, len(T), len(self.G)) * eigen_scale
                if np.any(values < -tolerance):
                    raise ArithmeticError("port Gram has a materially negative eigenvalue")
                positive = values > tolerance
                values = np.where(positive, values, 0.)
                port = _Port(T, V, values, U, positive)
                unique_ports.append(port)
            self.ports.append(port)
            self.targets.append(yw)
        self.port_groups = [(p, np.asarray([j for j, block in enumerate(self.ports)
                                           if block is p], dtype=int))
                            for p in unique_ports]
        # A nonempty A requires the whitened target to be reachable up to delta.
        floor2 = sum(float(np.linalg.norm(p.U[:, ~p.positive].conj().T @ y) ** 2)
                     for p, y in zip(self.ports, self.targets))
        self.minimum_reachable_residual = float(np.sqrt(floor2))
        self._energy_scale = max(sum(float(np.vdot(y, y).real) for y in self.targets),
                                 self.delta ** 2, np.finfo(float).tiny)
        self._ball_tolerance = 128 * np.finfo(float).eps * self._energy_scale
        if delta ** 2 < floor2 - self._ball_tolerance:
            raise UnreachableDataBall(self.delta, self.minimum_reachable_residual)
        self.setup_wall_seconds = time.perf_counter() - wall
        self.setup_cpu_seconds = time.process_time() - cpu
        self.last_A_info: dict[str, Any] = {}

    def _fields(self, X: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        X, Q = np.asarray(X, complex), np.asarray(Q, complex)
        if X.shape != self.E.shape or Q.shape != self.E.shape:
            raise ValueError(f"X/Q must both have shape {self.E.shape}")
        if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Q)):
            raise ValueError("projection state must be finite")
        return X, Q

    def affine_ball(self, X: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Project onto propagation graph intersected with the global data ball."""
        if self.grouped_ports:
            return self._affine_ball_grouped(X, Q)
        X, Q = self._fields(X, Q)
        self.costs["A_projections"] += 1
        self.costs["G_adjoint_action_calls"] += 1
        self.costs["G_adjoint_action_rhs"] += self.sources
        xu = cho_solve(self.H, X + self.B.conj().T @ (Q - self.E), check_finite=False)
        self.costs["H_iteration_rhs"] += self.sources
        residuals = []
        for j, (p, y) in enumerate(zip(self.ports, self.targets)):
            residuals.append(p.T @ xu[:, j] - y)
        self.costs["receiver_action_calls"] += self.sources
        self.costs["receiver_action_rhs"] += self.sources
        initial2 = sum(float(np.vdot(r, r).real) for r in residuals)
        multiplier, limit = 0., False
        x = xu.copy()
        if initial2 > self.delta ** 2:
            zs = [p.U.conj().T @ r for p, r in zip(self.ports, residuals)]
            floor2 = sum(float(np.vdot(z[~p.positive], z[~p.positive]).real)
                         for p, z in zip(self.ports, zs))
            if self.delta ** 2 < floor2 - self._ball_tolerance:
                raise UnreachableDataBall(self.delta, float(np.sqrt(floor2)))
            limit = self.delta ** 2 <= floor2 + self._ball_tolerance
            if limit:
                multiplier = float("inf")
            else:
                def discrepancy(mu: float) -> float:
                    self.costs["multiplier_evaluations"] += 1
                    value = floor2
                    for p, z in zip(self.ports, zs):
                        v = p.eigenvalues[p.positive]
                        with np.errstate(over="ignore"):
                            attenuated = z[p.positive] / (1 + mu * v)
                        value += float(np.vdot(attenuated, attenuated).real)
                    return value - self.delta ** 2
                maxeig = max(float(np.max(p.eigenvalues, initial=0.)) for p in self.ports)
                hi = 1. / maxeig
                while discrepancy(hi) > 0:
                    hi *= 4.
                    if not np.isfinite(hi):
                        raise ArithmeticError("finite data-ball multiplier exceeded float range")
                multiplier = brentq(discrepancy, 0., hi, xtol=np.nextafter(0., 1.),
                                     rtol=8 * np.finfo(float).eps)
            for j, (p, z) in enumerate(zip(self.ports, zs)):
                coeff = np.zeros_like(z)
                if limit:
                    coeff[p.positive] = z[p.positive] / p.eigenvalues[p.positive]
                else:
                    coeff[p.positive] = z[p.positive] / (1 / multiplier + p.eigenvalues[p.positive])
                x[:, j] -= p.V @ (p.U @ coeff)
                self.costs["receiver_adjoint_action_calls"] += 1
                self.costs["receiver_adjoint_action_rhs"] += 1
        self.costs["G_action_calls"] += 1
        self.costs["G_action_rhs"] += self.sources
        q = self.E + self.B @ x
        self.last_A_info = {"multiplier": multiplier, "pseudoinverse_limit": limit,
                            "unconstrained_residual": float(np.sqrt(initial2)),
                            "minimum_reachable_residual": self.minimum_reachable_residual}
        return x, q

    def _affine_ball_grouped(self, X: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Optional v2: shared-port matrix RHS and cached eigenspace powers.

        Only exact-identical port maps are grouped. All source blocks still
        share one multiplier. Unlike v1, each discrepancy evaluation uses only
        aggregated eigenspace powers and performs no source matrix actions.
        """
        X, Q = self._fields(X, Q)
        self.costs["A_projections"] += 1
        self.costs["G_adjoint_action_calls"] += 1
        self.costs["G_adjoint_action_rhs"] += self.sources
        xu = cho_solve(self.H, X + self.B.conj().T @ (Q - self.E), check_finite=False)
        self.costs["H_iteration_rhs"] += self.sources
        residuals = []
        for p, columns in self.port_groups:
            target = np.column_stack([self.targets[j] for j in columns])
            residuals.append((p, columns, p.T @ xu[:, columns] - target))
            self.costs["receiver_action_calls"] += 1
            self.costs["receiver_action_rhs"] += len(columns)
        initial2 = sum(float(np.vdot(r, r).real) for _, _, r in residuals)
        multiplier, limit = 0., False
        x = xu.copy()
        if initial2 > self.delta ** 2:
            blocks = []
            floor2 = 0.
            for p, columns, r in residuals:
                z = p.U.conj().T @ r
                power = np.sum(abs(z) ** 2, axis=1)
                floor2 += float(np.sum(power[~p.positive]))
                blocks.append((p, columns, z, power[p.positive]))
            if self.delta ** 2 < floor2 - self._ball_tolerance:
                raise UnreachableDataBall(self.delta, float(np.sqrt(floor2)))
            limit = self.delta ** 2 <= floor2 + self._ball_tolerance
            if limit:
                multiplier = float("inf")
            else:
                def discrepancy(mu: float) -> float:
                    self.costs["multiplier_evaluations"] += 1
                    value = floor2
                    for p, _, _, power in blocks:
                        with np.errstate(over="ignore"):
                            attenuation = 1 / (1 + mu * p.eigenvalues[p.positive])
                        value += float(power @ (attenuation * attenuation))
                    return value - self.delta ** 2
                maxeig = max(float(np.max(p.eigenvalues, initial=0.)) for p in self.ports)
                hi = 1. / maxeig
                while discrepancy(hi) > 0:
                    hi *= 4.
                    if not np.isfinite(hi):
                        raise ArithmeticError("finite data-ball multiplier exceeded float range")
                multiplier = brentq(discrepancy, 0., hi, xtol=np.nextafter(0., 1.),
                                     rtol=8 * np.finfo(float).eps)
            for p, columns, z, _ in blocks:
                coeff = np.zeros_like(z)
                denominator = p.eigenvalues[p.positive] + (0. if limit else 1 / multiplier)
                coeff[p.positive] = z[p.positive] / denominator[:, None]
                x[:, columns] -= p.V @ (p.U @ coeff)
                self.costs["receiver_adjoint_action_calls"] += 1
                self.costs["receiver_adjoint_action_rhs"] += len(columns)
        self.costs["G_action_calls"] += 1
        self.costs["G_action_rhs"] += self.sources
        q = self.E + self.B @ x
        self.last_A_info = {"multiplier": multiplier, "pseudoinverse_limit": limit,
                            "unconstrained_residual": float(np.sqrt(initial2)),
                            "minimum_reachable_residual": self.minimum_reachable_residual,
                            "implementation": "grouped_ports_v2"}
        return x, q

    def material(self, X: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Global cell projection onto one shared real scalar material box.

        All source AND vector-component entries contribute to the same cell's
        Gram. Stable quadratic roots and both box endpoints are evaluated.
        Degenerate/tied cases choose the first minimizer deterministically.
        """
        X, Q = self._fields(X, Q)
        self.costs["C_projections"] += 1
        xc = X.reshape(self.cells, self.components * self.sources)
        qc = Q.reshape(self.cells, self.components * self.sources)
        norm = np.maximum(np.max(np.abs(xc), axis=1), np.max(np.abs(qc), axis=1))
        norm = np.where(norm > 0, norm, 1.)
        xn, qn = xc / norm[:, None], qc / norm[:, None]
        A = np.sum(abs(xn) ** 2, axis=1)
        C = np.sum(abs(qn) ** 2, axis=1)
        B = np.real(np.sum(np.conj(xn) * qn, axis=1))
        d = C - A
        stable = -.5 * (d + np.copysign(np.hypot(d, 2 * B), d))
        root1 = np.divide(stable, B, out=np.full_like(B, np.nan), where=B != 0)
        root2 = np.divide(-B, stable, out=np.full_like(B, np.nan), where=stable != 0)
        # B=0 gives the stationary root alpha=0; endpoints resolve flat cases.
        root1 = np.where((B == 0) & (d != 0), 0., root1)
        lo, hi = self.lower / self.scale, self.upper / self.scale
        candidates = np.stack([lo, hi, root1, root2], axis=1)
        valid = np.isfinite(candidates) & (candidates >= lo[:, None]) & (candidates <= hi[:, None])
        candidates = np.where(valid, candidates, 0.)
        h = np.hypot(1., candidates)
        # Evaluate an actual squared distance, avoiding Gram cancellation.
        differences = (xn[:, None, :] / h[:, :, None]
                       - (candidates / h)[:, :, None] * qn[:, None, :])
        values = np.sum(abs(differences) ** 2, axis=2)
        values[~valid] = np.inf
        alpha = candidates[np.arange(self.cells), np.argmin(values, axis=1)]
        h = np.hypot(1., alpha)
        cq = (qc / h[:, None] + (alpha / h)[:, None] * xc) / h[:, None]
        cx = alpha[:, None] * cq
        return cx.reshape(self.E.shape), cq.reshape(self.E.shape), alpha * self.scale


def run_projection(method: str, projector: DenseFieldCurrentProjectors,
                   initial_material: np.ndarray,
                   audit_callback: Callable[[np.ndarray], PhysicalAudit], *,
                   steps: int, audit_every: int = 50, omega: float = 1.,
                   cpu_budget_seconds: float | None = None,
                   wall_budget_seconds: float | None = None,
                   checkpoint_callback: Callable[[AuditRecord], None] | None = None) -> ProjectionRun:
    """AP or relaxed DR with physical objective selection and complete audits.

    Runner clocks start before initial physical audit and exclude the separately
    measured projector setup. Pass remaining proposal budgets after setup.
    Budgets are checked at iteration boundaries, so an in-progress projection
    or physical audit can overrun. No extra physical solve occurs on budget exit.
    The last scheduled/final audit and its incumbent remain available. Calling
    ``checkpoint_callback`` at every audit permits incremental serialization.
    """
    aliases = {"ap": "alternating_projection", "dr": "douglas_rachford"}
    method = aliases.get(method.lower(), method.lower())
    if method not in ("alternating_projection", "douglas_rachford"):
        raise ValueError("method must be AP or DR")
    if steps < 0 or int(steps) != steps or audit_every < 1 or int(audit_every) != audit_every:
        raise ValueError("steps/audit_every must be nonnegative/positive integers")
    if not 0 < omega < 2:
        raise ValueError("DR relaxation omega must satisfy 0 < omega < 2")
    for budget in (cpu_budget_seconds, wall_budget_seconds):
        if budget is not None and (not np.isfinite(budget) or budget < 0):
            raise ValueError("budgets must be nonnegative finite seconds or None")
    m0 = np.asarray(initial_material, float).copy()
    if (m0.shape != (projector.cells,) or np.any(~np.isfinite(m0))
            or np.any(m0 < projector.lower) or np.any(m0 > projector.upper)):
        raise ValueError("initial material must be a finite cell vector within the box")
    wall, cpu = time.perf_counter(), time.process_time()
    costs_before = dict(projector.costs)
    physical_costs: dict[str, float] = {}
    history: list[AuditRecord] = []

    def physical(m: np.ndarray) -> PhysicalAudit:
        audit = audit_callback(m.copy())
        if not isinstance(audit, PhysicalAudit):
            raise TypeError("audit_callback must return PhysicalAudit")
        if not np.isfinite(audit.objective):
            raise ValueError("physical objective must be finite")
        projector._fields(audit.J, audit.Q)
        for key, value in audit.costs.items():
            if not np.isfinite(value) or value < 0:
                raise ValueError("physical callback costs must be finite nonnegative increments")
            physical_costs[key] = physical_costs.get(key, 0.) + float(value)
        # Backends may recycle field buffers. Keep physical incumbents stable
        # across subsequent calls independently of that implementation detail.
        return PhysicalAudit(float(audit.objective), np.array(audit.J, copy=True),
                             np.array(audit.Q, copy=True), deepcopy(dict(audit.metrics)),
                             dict(audit.costs))

    initial = physical(m0)
    X, Q = initial.J.copy() / projector.scale, initial.Q.copy()
    best, best_m, best_iteration = initial, m0.copy(), 0

    def record(iteration: int, m: np.ndarray, audit: PhysicalAudit,
               cx: np.ndarray, cq: np.ndarray, ax: np.ndarray, aq: np.ndarray,
               gap: float, previous_X: np.ndarray, previous_Q: np.ndarray) -> None:
        nonlocal best, best_m, best_iteration
        if audit.objective < best.objective:
            best, best_m, best_iteration = audit, m.copy(), iteration
        states = {"X": X.copy(), "Q": Q.copy(), "previous_X": previous_X.copy(),
                  "previous_Q": previous_Q.copy(), "C_X": cx.copy(), "C_Q": cq.copy(),
                  "A_X": ax.copy(), "A_Q": aq.copy(), "physical_J": audit.J.copy(),
                  "physical_Q": audit.Q.copy(), "incumbent_J": best.J.copy(),
                  "incumbent_Q": best.Q.copy()}
        item = AuditRecord(iteration, float(audit.objective), float(best.objective),
                           best_iteration, time.process_time() - cpu,
                           time.perf_counter() - wall, float(gap), m.copy(), best_m.copy(),
                           states, {**dict(audit.metrics), "A_projection": dict(projector.last_A_info)},
                           dict(best.metrics), dict(projector.costs), dict(physical_costs))
        history.append(item)
        if checkpoint_callback is not None:
            checkpoint_callback(item)

    record(0, m0, initial, X, Q, X, Q, 0., X, Q)
    completed, stop_reason = 0, "steps_completed"
    for iteration in range(1, int(steps) + 1):
        if cpu_budget_seconds is not None and time.process_time() - cpu >= cpu_budget_seconds:
            stop_reason = "cpu_budget"
            break
        if wall_budget_seconds is not None and time.perf_counter() - wall >= wall_budget_seconds:
            stop_reason = "wall_budget"
            break
        old_X, old_Q = X, Q
        cx, cq, cm = projector.material(old_X, old_Q)
        if method == "douglas_rachford":
            ax, aq = projector.affine_ball(2 * cx - old_X, 2 * cq - old_Q)
            X, Q = old_X + omega * (ax - cx), old_Q + omega * (aq - cq)
        else:
            ax, aq = projector.affine_ball(cx, cq)
            X, Q = ax, aq
        completed = iteration
        if iteration % audit_every == 0 or iteration == steps:
            audit = physical(cm)
            gap = np.hypot(np.linalg.norm(ax - cx), np.linalg.norm(aq - cq))
            record(iteration, cm, audit, cx, cq, ax, aq, gap, old_X, old_Q)
    # Per-run fixed-background costs exclude any earlier use of a shared object;
    # setup is reported once, separately, and included in cumulative audit costs.
    run_costs = {key: value - costs_before[key] for key, value in projector.costs.items()}
    return ProjectionRun(method, best_m.copy(), float(best.objective), best_iteration,
                         completed, stop_reason, history, X.copy(), Q.copy(), run_costs,
                         dict(physical_costs), projector.setup_cpu_seconds,
                         projector.setup_wall_seconds, time.process_time() - cpu,
                         time.perf_counter() - wall)
