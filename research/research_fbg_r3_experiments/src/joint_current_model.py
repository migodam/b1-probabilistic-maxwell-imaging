"""Exact nonlinear CSI/CIE residual and its correctly realified derivative."""
from __future__ import annotations

import time
from collections import Counter
import numpy as np

try:
    from .current_normal import real_pack, current_action
except ImportError:
    from current_normal import real_pack, current_action


class ReleaseBudgetExceeded(RuntimeError):
    pass


class ReleaseCost:
    """Measured backend deltas plus separately counted auxiliary operations."""
    def __init__(self, model, config=None):
        self.model, self.config = model, config or {}
        self.started = time.perf_counter()
        self.initial = dict(getattr(model, "cost", {}))
        self.aux = Counter()
        self.events = []

    def add(self, kind, count=1, **detail):
        self.aux[kind] += count
        if detail:
            self.events.append({"kind": kind, "count": count, **detail})

    def check(self, extra_lu=0, extra_rhs=0, extra_state=0):
        current = self.snapshot()
        limits = {"max_wall_seconds": current["wall_seconds"],
                  "max_full_states": self.aux["full_state_calls"] + extra_state,
                  "max_complex_rhs": current.get("actual_complex_RHS", self.aux["expected_complex_rhs"]) + extra_rhs,
                  "max_lu": current.get("LU", self.aux["expected_lu"]) + extra_lu,
                  "max_green_vector_actions": self.aux["green_vector_actions"]}
        for name, value in limits.items():
            if self.config.get(name) is not None and (value > self.config[name] or (name == "max_wall_seconds" and value >= self.config[name])):
                raise ReleaseBudgetExceeded(name)

    def snapshot(self):
        backend = {k: v - self.initial.get(k, 0) for k, v in dict(getattr(self.model, "cost", {})).items()}
        return {**backend, "auxiliary": dict(self.aux), "wall_seconds": time.perf_counter() - self.started,
                "events": self.events.copy(), "backend_counts_available": bool(hasattr(self.model, "cost"))}

    def state(self, m, purpose):
        cached = getattr(self.model, "_current", None)
        reused = cached is not None and np.array_equal(np.asarray(m), cached.get("m"))
        self.check(extra_lu=0 if reused else 1, extra_rhs=0 if reused else self.model.E.shape[1], extra_state=1)
        self.add("full_state_calls", purpose=purpose)
        if not reused:
            self.add("expected_lu")
            self.add("expected_complex_rhs", self.model.E.shape[1])
        return self.model.state(m)


def regularizer(model, m, config=None):
    config = config or {}
    fun = config.get("regularizer", getattr(model, "regularizer", None))
    if callable(fun):
        value, grad, hessian = fun(np.asarray(m, float))
        if hasattr(hessian, "toarray"):
            hessian = hessian.toarray()
        return float(value), np.asarray(grad, float), np.asarray(hessian, float)
    mat = config.get("regularizer_matrix")
    n = len(m)
    H = np.zeros((n, n)) if mat is None else np.asarray(mat, float)
    if H.ndim == 0:
        H = float(H) * np.eye(n)
    if H.ndim == 1:
        H = np.diag(H)
    d = np.asarray(m) - np.asarray(config.get("regularizer_center", np.zeros(n)))
    return float(0.5 * d @ H @ d), H @ d, H


def state_weights(m, config):
    chart = config.get("chart", "csi")
    if chart in ("csi", "identity", None):
        w, dw = np.ones(len(m)), np.zeros(len(m))
    elif chart in ("cie", "cie_normalized"):
        beta = float(config.get("beta", 0.2))
        den = 1.0 + beta * np.asarray(m)
        if np.any(np.abs(den) < config.get("pole_margin", 1e-8)):
            raise ValueError("CIE state weight pole")
        numerator = beta if chart == "cie" else 1.0
        w, dw = numerator / den, -numerator * beta / den ** 2
    else:
        raise ValueError(f"Unknown state chart: {chart}")
    fixed = np.asarray(config.get("state_weight", 1.0))
    return w * fixed, dw * fixed


def whiten_data(x, config):
    W = config.get("data_whitener")
    if callable(W):
        return np.asarray(W(x), float)
    if W is not None:
        W = np.asarray(W)
        return W * x if W.ndim <= 1 else W @ x
    return x / float(config.get("data_scale", 1.0))


def physical_task_config(config):
    """Optional original-physics whitening, separate from the joint objective."""
    out = dict(config)
    if "physical_data_whitener" in config:
        out["data_whitener"] = config["physical_data_whitener"]
    elif "physical_data_scale" in config:
        out.pop("data_whitener", None)
    if "physical_data_scale" in config:
        out["data_scale"] = config["physical_data_scale"]
    return out


def whiten_physical(x, config):
    """Original physical residual metric; explicit physical scale replaces W_y."""
    return whiten_data(x, physical_task_config(config))


class JointCurrentModel:
    def __init__(self, model, Y, config=None, cost=None):
        self.model, self.Y, self.config = model, np.asarray(Y, complex), config or {}
        self.G, self.E, self.S = np.asarray(model.G), np.asarray(model.E), np.asarray(model.S)
        self.n, self.p = self.E.shape
        self.pack = getattr(model, "pack", real_pack)
        self.rho = float(self.config.get("rho", 1.0))
        if self.rho <= 0:
            raise ValueError("rho must be positive")
        self.state_source_whitener = np.broadcast_to(np.asarray(self.config.get("state_source_whitener", 1.0), float), (self.p,)).copy()
        if np.any(~np.isfinite(self.state_source_whitener)) or np.any(self.state_source_whitener <= 0):
            raise ValueError("state_source_whitener must contain finite positive source weights")
        self.cost = cost

    def residual(self, m, j):
        m, j = np.asarray(m, float), np.asarray(j, complex)
        Q = self.E + self.G @ j
        e = j - m[:, None] * Q
        w, dw = state_weights(m, self.config)
        data = whiten_data(self.pack(self.S @ j - self.Y), self.config)
        state = np.sqrt(self.rho) * self.pack(w[:, None] * e * self.state_source_whitener[None, :])
        reg, grad, Hreg = regularizer(self.model, m, self.config)
        if self.cost:
            self.cost.add("joint_objective_evaluations")
            self.cost.add("green_vector_actions", self.p)
            self.cost.add("receiver_vector_actions", self.p)
        return {"b": np.r_[data, state], "data_residual": data, "state_residual": state,
                "Q": Q, "e": e, "w": w, "dw": dw,
                "objective": float(0.5 * (data @ data + state @ state) + reg),
                "data_loss": float(0.5 * data @ data), "state_loss": float(0.5 * state @ state),
                "defect_norm": float(np.linalg.norm(e)), "regularizer": reg, "ell": grad, "Lambda": Hreg}

    def current_matrix(self, m, modes, w=None):
        modes = np.asarray(modes, complex)
        if w is None:
            w = state_weights(m, self.config)[0]
        d = modes.shape[2]
        if d == 0:
            return np.empty((len(whiten_data(self.pack(np.zeros_like(self.Y)), self.config)) + 2 * self.E.size, 0))
        GV = np.einsum("ij,jpd->ipd", self.G, modes, optimize=True)
        LV = modes - np.asarray(m)[:, None, None] * GV
        SV = np.einsum("ij,jpd->ipd", self.S, modes, optimize=True)
        A = np.column_stack([np.r_[whiten_data(self.pack(SV[:, :, k]), self.config),
                                  np.sqrt(self.rho) * self.pack(w[:, None] * LV[:, :, k] * self.state_source_whitener[None, :])] for k in range(d)])
        if self.cost:
            self.cost.add("green_vector_actions", self.p * d)
            self.cost.add("receiver_vector_actions", self.p * d)
        return A

    def linearize(self, m, j, modes):
        out = self.residual(m, j)
        # For each real material component: -P_i Q_i + P'_i e_i.
        injection = (-out["w"][:, None] * out["Q"] + out["dw"][:, None] * out["e"]) * self.state_source_whitener[None, :]
        nd = len(out["data_residual"])
        B = np.zeros((len(out["b"]), self.n))
        for k in range(self.n):
            z = np.zeros_like(j)
            z[k] = injection[k]
            B[nd:, k] = np.sqrt(self.rho) * self.pack(z)
        out["B"] = B
        out["A"] = self.current_matrix(m, modes, out["w"])
        if self.cost:
            self.cost.add("joint_linearizations")
        return out


def joint_linearization(model, Y, m, j, modes, config=None, cost=None):
    return JointCurrentModel(model, Y, config, cost).linearize(m, j, modes)
