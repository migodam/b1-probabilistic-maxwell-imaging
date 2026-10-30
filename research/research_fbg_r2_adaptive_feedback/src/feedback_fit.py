"""Fixed feedback operators for the scalar-TM current formulation.

The physical convention is L=I-diag(m)G, j=L^-1 diag(m)E, Q=E+Gj,
R=Sj-Y. For a feedback B fixed throughout an optimization, V=Q+Bj.
All complex inner products use conjugation. Material perturbations are real.
Only operator construction/safety checks use small matrix SVDs; objective
callbacks use the retained physical LU and a small source Gram Cholesky.
"""
from dataclasses import dataclass
import numpy as np

try:
    import torch
except ImportError:  # CPU-only inspection remains possible.
    torch = None


def _torch(x):
    return torch is not None and isinstance(x, torch.Tensor)


def _as(x, like):
    if _torch(like):
        return torch.as_tensor(x, dtype=like.dtype, device=like.device)
    return np.asarray(x, dtype=like.dtype)


def _eye(n, like):
    return torch.eye(n, dtype=like.dtype, device=like.device) if _torch(like) else np.eye(n, dtype=like.dtype)


def _herm(x):
    return x.mH if _torch(x) else np.swapaxes(x.conj(), -1, -2)


def _transpose(x):
    return x.transpose(-1, -2) if _torch(x) else np.swapaxes(x, -1, -2)


def _solve(a, b):
    return torch.linalg.solve(a, b) if _torch(a) else np.linalg.solve(a, b)


def _real(x):
    return x.real


def _maximum_zero(x):
    return torch.clamp_min(x, 0) if _torch(x) else np.maximum(x, 0)


def _scalar(x):
    return float(x.detach().cpu()) if _torch(x) else float(x)


@dataclass
class FeedbackOperator:
    """B as zero, diagonal, dense, or independent 2x2 spatial blocks.

    ``blocks`` has shape (number_of_blocks, 4, 4), with each index group
    corresponding to a 2x2 square on the row-major material grid. ``data``
    stays on the backend/device supplied by the caller. No physical G is copied.
    """

    kind: str
    data: object = None
    indices: object = None
    n: int = 0
    name: str = ""
    reference_scale: object = 1.
    details: object = None

    @classmethod
    def zero(cls, n):
        return cls("zero", n=n)

    @classmethod
    def diagonal(cls, values, name=""):
        return cls("diagonal", data=values, n=len(values), name=name)

    @classmethod
    def dense(cls, matrix):
        return cls("dense", data=matrix, n=matrix.shape[0])

    @classmethod
    def blocks(cls, matrices, indices, n):
        return cls("blocks", data=matrices, indices=indices, n=n)

    @property
    def diag(self):
        if self.kind != "diagonal":
            raise AttributeError("diag is only defined for diagonal feedback")
        return self.data

    def apply(self, x):
        if self.kind == "zero":
            return x * 0
        data = _as(self.data, x)
        if self.kind == "diagonal":
            if x.ndim == 1:
                return data * x
            return data.reshape(*([1] * (x.ndim - 2)), self.n, 1) * x
        if self.kind == "dense":
            return data @ x
        if self.kind == "blocks":
            out = x * 0
            idx = torch.as_tensor(self.indices, dtype=torch.long, device=x.device) if _torch(x) else self.indices
            if x.ndim == 1:
                out[idx] = (data @ x[idx][..., None])[..., 0]
            elif x.ndim == 2:
                out[idx] = data @ x[idx]
            else:
                out[..., idx, :] = data @ x[..., idx, :]
            return out
        raise ValueError(f"unknown feedback kind: {self.kind}")

    def adjoint(self, x):
        if self.kind == "zero":
            return x * 0
        data = _as(self.data, x)
        if self.kind == "diagonal":
            d = data.conj()
            if x.ndim == 1:
                return d * x
            return d.reshape(*([1] * (x.ndim - 2)), self.n, 1) * x
        if self.kind == "dense":
            return _herm(data) @ x
        if self.kind == "blocks":
            out = x * 0
            idx = torch.as_tensor(self.indices, dtype=torch.long, device=x.device) if _torch(x) else self.indices
            if x.ndim == 1:
                out[idx] = (_herm(data) @ x[idx][..., None])[..., 0]
            elif x.ndim == 2:
                out[idx] = _herm(data) @ x[idx]
            else:
                out[..., idx, :] = _herm(data) @ x[..., idx, :]
            return out
        raise ValueError(f"unknown feedback kind: {self.kind}")

    def field(self, q, m):
        """Apply (I+B diag(m))Q without materializing a dense matrix."""
        return q + self.apply(m[:, None] * q)

    def tangent(self, G, current_tangent):
        """Exact fixed-B V_h=(G+B)j_h; no material finite differences."""
        return G @ current_tangent + self.apply(current_tangent)

    def scaled(self, scale):
        if self.kind == "zero":
            return self
        result = FeedbackOperator(self.kind, self.data * scale, self.indices, self.n, self.name,
                                  self.reference_scale, self.details)
        for attribute in ("fit_lambda", "fit_ridge_fraction"):
            if hasattr(self, attribute):
                setattr(result, attribute, getattr(self, attribute))
        return result

    def safety(self, m, lower=0.02, upper=20.):
        """One-anchor Mf singular/condition receipt; call at fit time only."""
        if self.kind == "zero":
            return {"min_singular": 1., "max_singular": 1., "condition": 1.,
                    "safe": True, "lower": lower, "upper": upper, "method": "identity"}
        if self.kind == "diagonal":
            sv = (_as(self.data, m + 0j) * m + 1).abs() if _torch(m) else np.abs(1 + self.data * m)
            method = "diagonal_exact"
        elif self.kind == "blocks":
            idx = torch.as_tensor(self.indices, dtype=torch.long, device=m.device) if _torch(m) else self.indices
            B = self.data
            eye = torch.eye(4, dtype=B.dtype, device=B.device) if _torch(B) else np.eye(4, dtype=B.dtype)
            Mf = eye + B * m[idx][..., None, :]
            sv = torch.linalg.svdvals(Mf) if _torch(Mf) else np.linalg.svd(Mf, compute_uv=False)
            method = "batched_4x4_exact"
        else:
            B = self.dense_matrix(self.data)
            Mf = _eye(self.n, B) + B * m[None, :]
            sv = torch.linalg.svdvals(Mf) if _torch(Mf) else np.linalg.svd(Mf, compute_uv=False)
            method = "dense_exact_once"
        smallest, largest = _scalar(sv.min()), _scalar(sv.max())
        condition = largest / smallest if smallest else float("inf")
        return {"min_singular": smallest, "max_singular": largest,
                "condition": condition, "safe": bool(smallest >= lower and condition <= upper),
                "lower": lower, "upper": upper, "method": method}

    def box_certify(self, lower_m=0.02, upper_m=20., perturbation_cap=0.8):
        """Return once-scaled B and a certificate over the entire material box.

        For diagonal B, min/max |1+b_i m_i| are exact for every independent
        real m_i in [lower_m,upper_m]. For block/full B, one spectral norm
        bounds ||B diag(m)||₂≤upper_m||B||₂. Scaling to <=0.8 gives
        sigma_min(Mf)>=0.2, sigma_max<=1.8, condition<=9 throughout the box.
        A returned scale is an engineering safety operation, not fitted data.
        """
        if lower_m < 0 or upper_m < lower_m or not 0 < perturbation_cap < 1:
            raise ValueError("invalid material box or perturbation cap")
        if self.kind == "zero":
            return self, {"method": "identity_exact", "scale": 1., "min_singular_bound": 1.,
                          "max_singular_bound": 1., "condition_bound": 1., "whole_box": True}
        if self.kind == "diagonal":
            def bounds(b):
                abs2 = (b.abs() if _torch(b) else np.abs(b)) ** 2
                real = b.real
                if _torch(b):
                    vertex = torch.where(abs2 > 0, -real / torch.clamp_min(abs2, 1e-300),
                                         torch.full_like(real, lower_m))
                    vertex = torch.clamp(vertex, lower_m, upper_m)
                    lo = torch.minimum(torch.minimum((1 + b * lower_m).abs(), (1 + b * upper_m).abs()),
                                       (1 + b * vertex).abs())
                    hi = torch.maximum((1 + b * lower_m).abs(), (1 + b * upper_m).abs())
                else:
                    vertex = np.where(abs2 > 0, -real / np.maximum(abs2, 1e-300), lower_m)
                    vertex = np.clip(vertex, lower_m, upper_m)
                    lo = np.minimum.reduce([np.abs(1 + b * lower_m), np.abs(1 + b * upper_m),
                                            np.abs(1 + b * vertex)])
                    hi = np.maximum(np.abs(1 + b * lower_m), np.abs(1 + b * upper_m))
                smallest, largest = _scalar(lo.min()), _scalar(hi.max())
                condition = largest / smallest if smallest else float("inf")
                return smallest, largest, condition
            threshold = (1 + perturbation_cap) / (1 - perturbation_cap)
            smallest, largest, condition = bounds(self.data)
            scale = 1.
            if smallest < 1 - perturbation_cap or condition > threshold:
                scale = .25
                smallest, largest, condition = bounds(self.data * scale)
            safe = smallest >= 1 - perturbation_cap and condition <= threshold
            result = self.scaled(scale) if safe else FeedbackOperator("zero", n=self.n, name=self.name)
            return result, {"method": "diagonal_interval_exact", "scale": scale if safe else 0.,
                            "fallback_scale_attempted": scale if scale < 1 else None,
                            "rejected_to_zero": not safe,
                            "min_singular_bound": smallest, "max_singular_bound": largest,
                            "condition_bound": condition, "whole_box": True, "safe": bool(safe)}
        if self.kind == "blocks":
            sv = torch.linalg.svdvals(self.data) if _torch(self.data) else np.linalg.svd(self.data, compute_uv=False)
        else:
            sv = torch.linalg.svdvals(self.data) if _torch(self.data) else np.linalg.svd(self.data, compute_uv=False)
        norm = _scalar(sv.max())
        scale = min(1., perturbation_cap / (upper_m * norm)) if norm else 1.
        delta = scale * upper_m * norm
        result = self.scaled(scale)
        return result, {"method": "spectral_norm_box_bound", "scale": scale,
                        "unscaled_operator_norm": norm, "min_singular_bound": 1 - delta,
                        "max_singular_bound": 1 + delta,
                        "condition_bound": (1 + delta) / (1 - delta),
                        "whole_box": True, "safe": True}

    def dense_matrix(self, like):
        out = _eye(self.n, like) * 0
        if self.kind == "zero":
            return out
        if self.kind == "dense":
            return _as(self.data, like)
        if self.kind == "diagonal":
            d = _as(self.data, like)
            return torch.diag(d) if _torch(d) else np.diag(d)
        data = _as(self.data, like)
        for k, ind in enumerate(self.indices):
            for a, i in enumerate(ind):
                for b, j in enumerate(ind):
                    out[i, j] = data[k, a, b]
        return out


def spatial_2x2_indices(side):
    """Disjoint 2x2 groups on an even, row-major square grid."""
    if side % 2:
        raise ValueError("2x2 spatial blocks require an even side length")
    return np.asarray([[r * side + c, r * side + c + 1,
                        (r + 1) * side + c, (r + 1) * side + c + 1]
                       for r in range(0, side, 2) for c in range(0, side, 2)], dtype=np.int64)


def source_normalized_current_tangents(G, E, m, directions, gram_floor=0.):
    """Return exact current tangents from source-whitened fields at one anchor.

    ``directions`` is n x k (real). H=Q*Q is Cholesky factored once and
    Q_white=Q chol(H)^-*; each current tangent is
    L^-1 diag(h)Q_white. All k source blocks share the same physical L.
    The optional floor is an explicit ridge on H and must be reported by caller.
    """
    n = G.shape[0]
    if directions.ndim == 1:
        directions = directions[:, None]
    L = _eye(n, G) - m[:, None] * G
    J = _solve(L, m[:, None] * E)
    Q = E + G @ J
    H = _herm(Q) @ Q
    if gram_floor:
        H = H + gram_floor * _eye(H.shape[0], H)
    C = torch.linalg.cholesky(H) if _torch(H) else np.linalg.cholesky(H)
    Qw = Q @ _solve(_herm(C), _eye(C.shape[0], C))
    blocks = directions.T[:, :, None] * Qw[None, :, :]
    rhs = (blocks.permute(1, 0, 2) if _torch(blocks) else blocks.transpose(1, 0, 2)).reshape(n, -1)
    return _solve(L, rhs), {"field": Q, "white_field": Qw, "source_gram": H}


def _ridge_rows(X, Y, ridge, real=False):
    """Fit A X≈-Y with Frobenius ridge; X has few input rows."""
    if real:
        gram = _real(X @ _herm(X))
        cross = _real(Y @ _herm(X))
    else:
        gram = X @ _herm(X)
        cross = Y @ _herm(X)
    system = gram + ridge * _eye(gram.shape[-1], gram)
    return -_transpose(_solve(_transpose(system), _transpose(cross)))


def fit_tangent_ridge(G, T, structure="diagonal", ridge=1e-8, side=None):
    """Fit fixed B to attenuate exact transformed-field tangents (G+B)T.

    ``T`` is a concatenation of exact source-normalized current tangents.
    The full candidate solves min_B ||(G+B)T||_F²+ridge||B||_F².
    Diagonal and spatial-block variants enforce their respective sparsity.
    ``real_nonnegative`` sets B=diag(k), k>=0. ``real`` permits either
    sign. ``complex`` permits complex diagonal entries. All fits use ridge>0.
    This fits tangent attenuation at supplied anchors; it does not establish
    projector curvature or an optimization-basin guarantee.
    """
    if ridge <= 0:
        raise ValueError("ridge must be positive")
    n = G.shape[0]
    Y = G @ T
    if structure == "zero":
        return FeedbackOperator.zero(n)
    if structure == "full":
        return FeedbackOperator.dense(_ridge_rows(T, Y, ridge))
    if structure in ("real", "real_nonnegative", "complex", "diagonal"):
        numerator = (Y * T.conj()).sum(axis=1)
        denominator = ((T.abs() if _torch(T) else np.abs(T)) ** 2).sum(axis=1) + ridge
        if structure == "real_nonnegative":
            coeff = _maximum_zero(-_real(numerator) / denominator)
        elif structure == "real":
            coeff = -_real(numerator) / denominator
        else:
            coeff = -numerator / denominator
        return FeedbackOperator.diagonal(coeff)
    if structure in ("block", "block_real"):
        if side is None or side * side != n:
            raise ValueError("block fit requires side with side²=n")
        indices = spatial_2x2_indices(side)
        idx = torch.as_tensor(indices, dtype=torch.long, device=T.device) if _torch(T) else indices
        data = _ridge_rows(T[idx], Y[idx], ridge, real=(structure == "block_real"))
        return FeedbackOperator.blocks(data, indices, n)
    raise ValueError(f"unknown fit structure: {structure}")


def fit_feedback(G, Q, tangents, kind="complex", side=None, ridge=0.05):
    """Adapter for cached exact anchor tangents of shape (directions,n,sources).

    Whiten every source block by the anchor H=Q*Q Cholesky before fitting.
    Consequently the fit minimizes source-normalized transformed-field
    tangent energy. ``Q`` and ``tangents`` must come from the same anchor.
    """
    if tangents.ndim != 3 or tangents.shape[1] != G.shape[0] or tangents.shape[2] != Q.shape[1]:
        raise ValueError("expected tangents with shape (directions,n,sources)")
    H = _herm(Q) @ Q
    C = torch.linalg.cholesky(H) if _torch(H) else np.linalg.cholesky(H)
    right = _solve(_herm(C), _eye(C.shape[0], C))
    white = tangents @ right
    T = (white.permute(1, 0, 2) if _torch(white) else white.transpose(1, 0, 2)).reshape(G.shape[0], -1)
    row_power = ((T.abs() if _torch(T) else np.abs(T)) ** 2).sum(axis=1)
    actual_lambda = ridge * _scalar(row_power.mean())
    if actual_lambda <= 0:
        raise ValueError("tangent power is zero; feedback fit is undefined")
    feedback = fit_tangent_ridge(G, T, structure=kind, ridge=actual_lambda, side=side)
    feedback.fit_lambda = actual_lambda
    feedback.fit_ridge_fraction = ridge
    return feedback


def candidate_feedback(G, name, beta=0., fitted=None, rng=None):
    """Create zero, scalar beta I, self=-diag(G), fitted or shuffled control."""
    n = G.shape[0]
    if name == "zero":
        return FeedbackOperator.zero(n)
    if name == "beta":
        return FeedbackOperator.diagonal(_as(beta, G) * (torch.ones(n, dtype=G.dtype, device=G.device) if _torch(G) else np.ones(n, dtype=G.dtype)))
    if name == "self":
        return FeedbackOperator.diagonal(-(torch.diagonal(G) if _torch(G) else np.diag(G)))
    if name == "fitted":
        if fitted is None:
            raise ValueError("fitted candidate requires a FeedbackOperator")
        return fitted
    if name == "shuffled":
        if fitted is None:
            raise ValueError("shuffled control requires a fitted FeedbackOperator")
        # Permute spatial labels on both axes: Frobenius and spectral norms match.
        order = np.random.default_rng(rng).permutation(n)
        B = fitted.dense_matrix(G)
        if _torch(B):
            p = torch.as_tensor(order, device=B.device)
            B = B[p][:, p]
        else:
            B = B[np.ix_(order, order)]
        return FeedbackOperator.dense(B)
    raise ValueError(f"unknown candidate: {name}")


def cap_feedback_at_anchors(feedback, anchors, min_singular=0.1, cond_cap=100., scale_cap=None, iterations=32):
    """Scale a candidate once so Mf=I+B diag(m) passes fixed anchor checks.

    The function checks only supplied anchors and the final scaled candidate;
    it makes no assertion between anchors. ``scale_cap`` bounds ||B diag(m)||₂
    at each anchor. A monotone-safe bisection starts from alpha=0, and checks
    min singular value and condition number of Mf. This setup is not a callback.
    """
    if len(anchors) == 0:
        raise ValueError("at least one anchor is required")
    if min_singular <= 0 or cond_cap < 1:
        raise ValueError("invalid singular/condition thresholds")
    mats = [feedback.dense_matrix(feedback.data if feedback.data is not None else m) * m[None, :]
            for m in anchors]
    def measure(alpha):
        min_sv, max_cond, max_pert = float("inf"), 0., 0.
        for A in mats:
            Mf = _eye(feedback.n, A) + alpha * A
            sv = torch.linalg.svdvals(Mf) if _torch(A) else np.linalg.svd(Mf, compute_uv=False)
            hi, lo = _scalar(sv.max()), _scalar(sv.min())
            min_sv = min(min_sv, lo)
            max_cond = max(max_cond, hi / lo if lo else float("inf"))
            if scale_cap is not None:
                ps = torch.linalg.svdvals(alpha * A) if _torch(A) else np.linalg.svd(alpha * A, compute_uv=False)
                max_pert = max(max_pert, _scalar(ps.max()))
        return min_sv, max_cond, max_pert
    def passes(v):
        return v[0] >= min_singular and v[1] <= cond_cap and (scale_cap is None or v[2] <= scale_cap)
    lo, hi = 0., 1.
    if not passes(measure(hi)):
        for _ in range(iterations):
            mid = (lo + hi) / 2
            if passes(measure(mid)):
                lo = mid
            else:
                hi = mid
    else:
        lo = 1.
    metrics = measure(lo)
    return feedback.scaled(lo), {"scale": lo, "min_singular_at_anchors": metrics[0],
                                 "max_condition_at_anchors": metrics[1],
                                 "max_feedback_perturbation_at_anchors": metrics[2],
                                 "anchor_count": len(anchors)}


def feedback_value_gradient(shared, state, feedback, factor):
    """Torch device objective/analytic real-material gradient for fixed B.

    Reuses ``ObjectiveLane._forward`` state from CUDA ``src/device.py``.
    ``factor`` is the fixed .5*c normalization chosen from a reference state;
    it must never depend on the current material. With V=Q+Bj and H=V*V,
    psi=factor tr(R H^-1 R*), W=R H^-1, D=W*W. The differential is
    dpsi=2factor Re[<W,dR>-<VD,dV>]. The direct B diag(h)Q term in
    dV=B diag(h)Q+(I+B diag(m))G L^-1 diag(h)Q combines exactly with the
    second term to give dV=(G+B)L^-1 diag(h)Q. The single adjoint RHS is
    S*W-(G+B)*VD. Returns device ``(value, gradient_m, gram_info)``.
    """
    if torch is None or not _torch(state["q"]):
        raise TypeError("feedback_value_gradient requires a Torch lane state")
    q, r, m = state["q"], state["r"], state["m"]
    v = feedback.field(q, m)
    h = v.mH @ v
    h = .5 * (h + h.mH)
    chol, info = torch.linalg.cholesky_ex(h, check_errors=False)
    w = torch.cholesky_solve(r.mH, chol).mH
    d = w.mH @ w
    vd = v @ d
    rhs = shared.SH @ w - (shared.GH @ vd + feedback.adjoint(vd))
    adj = torch.linalg.lu_solve(state["lu"], state["pivots"], rhs, adjoint=True)
    grad = 2 * factor * (adj.conj() * q).sum(dim=1).real
    value = factor * (r.conj() * w).sum().real
    return value, grad, info


def feedback_jvp(G, S, state, feedback, direction):
    """Exact real-material JVP (dR,dV) using the same retained Torch LU."""
    q = state["q"]
    if _torch(q):
        dj = torch.linalg.lu_solve(state["lu"], state["pivots"], direction[:, None] * q)
    else:
        dj = _solve(state["L"], direction[:, None] * q)
    return S @ dj, feedback.tangent(G, dj)
