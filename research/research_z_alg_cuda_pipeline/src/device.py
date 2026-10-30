"""Same complex128 VIE and objectives as R4; one host return per evaluation.

The CUDA context and fixed operators are shared. Each optimizer has a private
stream, retained LU, pinned staging buffers and completion event. There is no
device-wide barrier in the evaluation path. Numerical status is returned with
the requested result rather than causing extra host round trips.
"""
from collections import Counter
import time
import numpy as np
import torch


def basis(fine, coarse):
    if fine == coarse:
        return np.eye(fine * fine)
    t = np.linspace(0., coarse - 1., fine)
    lo = np.minimum(t.astype(int), coarse - 2)
    a = t - lo
    b = np.zeros((fine, coarse))
    b[np.arange(fine), lo] = 1. - a
    b[np.arange(fine), lo + 1] = a
    return np.kron(b, b)


class SharedOperators:
    def __init__(self, public, device="cuda:0"):
        self.device = torch.device(device)
        torch.cuda.set_device(self.device)
        self.G = torch.as_tensor(public["G"], dtype=torch.complex128, device=self.device)
        self.E0 = torch.as_tensor(public.get("E0", public.get("U0")), dtype=torch.complex128, device=self.device)
        self.S = torch.as_tensor(public.get("S", public.get("H")), dtype=torch.complex128, device=self.device)
        self.GH, self.SH = self.G.mH.contiguous(), self.S.mH.contiguous()
        self.n, self.ns, self.nr = self.G.shape[0], self.E0.shape[1], self.S.shape[0]
        self.side = int(np.sqrt(self.n))
        self.eye = torch.eye(self.n, dtype=torch.complex128, device=self.device)
        self.bases = {s: torch.as_tensor(basis(self.side, s), dtype=torch.float64, device=self.device)
                      for s in (8, 12, self.side)}
        grid = np.arange(self.n).reshape(self.side, self.side)
        left = np.concatenate((grid[:, :-1].ravel(), grid[:-1, :].ravel()))
        right = np.concatenate((grid[:, 1:].ravel(), grid[1:, :].ravel()))
        self.left = torch.as_tensor(left, device=self.device)
        self.right = torch.as_tensor(right, device=self.device)
        self.edges = len(left)
        # Same R4 .1 public reference; one shared LU charged as setup.
        f = torch.linalg.lu_factor_ex(self.eye - .1 * self.G, check_errors=False)
        q = self.E0 + self.G @ torch.linalg.lu_solve(f.LU, f.pivots, .1 * self.E0)
        self.reference_field_mean = (q.abs().square().sum() / self.ns).real
        self.ready = torch.cuda.Event()
        self.ready.record()
        self.setup_counts = {"LU": 1, "source_RHS": self.ns, "operator_uploads": 1}

    def score_batch(self, materials, targets):
        """Exact nonlinear candidate screening; independent materials share G/E/S.

        Inputs/outputs are device tensors. Caller chooses chunks/streams and the
        point at which a host result is needed. No implicit Python status read.
        """
        a = self.eye - materials[..., :, None] * self.G
        f = torch.linalg.lu_factor_ex(a, check_errors=False)
        j = torch.linalg.lu_solve(f.LU, f.pivots, materials[..., :, None] * self.E0)
        y = self.S @ j
        r = y - targets
        relative = torch.sqrt(r.abs().square().sum((-2, -1)) / targets.abs().square().sum((-2, -1)))
        return relative, f.info


class ObjectiveLane:
    def __init__(self, shared, public, stage=24, objective="data", config=None):
        self.shared, self.stage, self.objective = shared, stage, objective
        self.config = {"smooth_weight": 1e-4, "tv_weight": 1e-4,
                       "tv_epsilon": .01, **(config or {})}
        self.stream = torch.cuda.Stream(device=shared.device)
        self.stream.wait_event(shared.ready)
        self.counts = Counter()
        self.n_coeffs = stage * stage
        self.started = time.perf_counter()
        self.wait_seconds = 0.
        self.last_x = None
        self.state = None
        self.events = []
        with torch.cuda.stream(self.stream):
            self.target = torch.as_tensor(public.get("Y", public.get("y")), dtype=torch.complex128, device=shared.device)
            self.scale2 = self.target.abs().square().sum()
            self.B = shared.bases[stage]
            self.DB = (self.B[shared.right] - self.B[shared.left]) / 19.5
            self.input_host = torch.empty(self.n_coeffs, dtype=torch.float64, pin_memory=True)
            self.input_device = torch.empty(self.n_coeffs, dtype=torch.float64, device=shared.device)
            self.output_host = torch.empty(self.n_coeffs + 5, dtype=torch.float64, pin_memory=True)
            self.done = torch.cuda.Event()

    def _pull(self, value, buffer=None):
        if buffer is None:
            buffer = torch.empty(value.shape, dtype=value.dtype, pin_memory=True)
        buffer.copy_(value, non_blocking=True)
        self.done.record(self.stream)
        t = time.perf_counter()
        self.done.synchronize()  # This result is required by this optimizer only.
        self.wait_seconds += time.perf_counter() - t
        self.counts["host_result_returns"] += 1
        return buffer.numpy().copy()

    def _forward(self, x):
        x = np.asarray(x, dtype=np.float64)
        if self.last_x is not None and np.array_equal(x, self.last_x):
            return self.state
        s = self.shared
        self.input_host.copy_(torch.from_numpy(x))
        self.input_device.copy_(self.input_host, non_blocking=True)
        m = self.input_device if self.stage == s.side else self.B @ self.input_device
        a = s.eye - m[:, None] * s.G
        f = torch.linalg.lu_factor_ex(a, check_errors=False)
        b = m[:, None] * s.E0
        j = torch.linalg.lu_solve(f.LU, f.pivots, b)
        q, y = s.E0 + s.G @ j, s.S @ j
        r = y - self.target
        data = .5 * r.abs().square().sum() / self.scale2
        d = (m[s.right] - m[s.left]) / 19.5
        tvq = torch.sqrt(d * d + self.config["tv_epsilon"] ** 2)
        reg = (self.config["smooth_weight"] * d.square().mean()
               + self.config["tv_weight"] * (tvq - self.config["tv_epsilon"]).mean())
        self.state = dict(m=m, a=a, lu=f.LU, pivots=f.pivots, info=f.info,
                          j=j, q=q, y=y, r=r, data=data, d=d, tvq=tvq, reg=reg)
        self.last_x = x.copy()
        self.counts.update({"LU": 1, "source_RHS": s.ns, "forward": 1})
        return self.state

    def evaluate(self, x):
        s = self.shared
        with torch.cuda.stream(self.stream), torch.no_grad():
            z = self._forward(x)
            q, r = z["q"], z["r"]
            gram_info = torch.zeros((), dtype=torch.int32, device=s.device)
            if self.objective == "dual":
                h = q.mH @ q
                h = .5 * (h + h.mH)
                chol, gram_info = torch.linalg.cholesky_ex(h, check_errors=False)
                weighted = torch.cholesky_solve(r.mH, chol).mH
                factor = .5 * s.reference_field_mean / self.scale2
                value = factor * (r.conj() * weighted).sum().real
                rhs = s.SH @ weighted - s.GH @ (q @ (weighted.mH @ weighted))
                adj = torch.linalg.lu_solve(z["lu"], z["pivots"], rhs, adjoint=True)
                g = 2 * factor * (adj.conj() * q).sum(1).real
                self.counts["Gram"] += 1
            else:
                value = z["data"]
                adj = torch.linalg.lu_solve(z["lu"], z["pivots"], s.SH @ r, adjoint=True)
                g = (adj.conj() * q).sum(1).real / self.scale2
            d = z["d"]
            gd = (2 * self.config["smooth_weight"] * d + self.config["tv_weight"] * d / z["tvq"]) / (s.edges * 19.5)
            greg = torch.zeros_like(z["m"])
            greg.index_add_(0, s.right, gd)
            greg.index_add_(0, s.left, -gd)
            g = g + greg
            if self.stage != s.side:
                g = self.B.T @ g
            packed = torch.cat(((value + z["reg"]).reshape(1), g,
                                z["data"].reshape(1), z["reg"].reshape(1),
                                z["info"].reshape(1).to(torch.float64), gram_info.reshape(1).to(torch.float64)))
            result = self._pull(packed, self.output_host)
            self.counts["adjoint_RHS"] += s.ns
            self.counts["gradient"] += 1
        # These flags arrived with the requested optimizer result; no new sync.
        if result[-1] != 0 or result[-2] != 0:
            raise ArithmeticError(f"LU/Gram factorization failed: {result[-2:]}")
        self.last_physical_objective = float(result[-4] + result[-3])
        return float(result[0]), result[1:1 + self.n_coeffs]

    def residual(self, x):
        s = self.shared
        with torch.cuda.stream(self.stream), torch.no_grad():
            z = self._forward(x)
            data = torch.view_as_real(z["r"].T.contiguous()).reshape(-1) / torch.sqrt(self.scale2)
            a = np.sqrt(2 * self.config["smooth_weight"] / s.edges)
            b = np.sqrt(2 * self.config["tv_weight"] / s.edges)
            return self._pull(torch.cat((data, a * z["d"], b * z["d"] / torch.sqrt(z["tvq"] + self.config["tv_epsilon"]))))

    def jacobian(self, x):
        s = self.shared
        with torch.cuda.stream(self.stream), torch.no_grad():
            z = self._forward(x)
            t = torch.linalg.lu_solve(z["lu"], z["pivots"], s.S.T.conj(), adjoint=True).conj().T
            jc = torch.einsum("ri,is->sri", t, z["q"]).reshape(s.ns * s.nr, s.n)
            j = torch.stack((jc.real, jc.imag), dim=1).reshape(2 * s.ns * s.nr, s.n)
            if self.stage != s.side:
                j = j @ self.B
            j = j / torch.sqrt(self.scale2)
            a = np.sqrt(2 * self.config["smooth_weight"] / s.edges)
            d, q = z["d"], z["tvq"]
            td = np.sqrt(2 * self.config["tv_weight"] / s.edges) * (
                1 / torch.sqrt(q + self.config["tv_epsilon"])
                - d * d / (2 * q * (q + self.config["tv_epsilon"]) ** 1.5))
            out = self._pull(torch.cat((j, a * self.DB, td[:, None] * self.DB)))
            self.counts.update({"FullJ": 1, "receiver_RHS": s.nr})
            return out

    def close(self, x):
        s = self.shared
        with torch.cuda.stream(self.stream), torch.no_grad():
            z = self._forward(x)
            b = z["m"][:, None] * s.E0
            state_residual = torch.linalg.vector_norm(z["a"] @ z["j"] - b) / torch.linalg.vector_norm(b)
            values = torch.cat((z["m"], (z["data"] + z["reg"]).reshape(1),
                                torch.sqrt(2 * z["data"]).reshape(1), state_residual.reshape(1), z["info"].reshape(1).double()))
            result = self._pull(values)
        return {"m": result[:s.n].tolist(), "ranking_objective": float(result[-4]),
                "relative_residual": float(result[-3]), "state_relative_residual": float(result[-2]),
                "lu_info": int(result[-1]), "counts": dict(self.counts),
                "host_result_wait_seconds": self.wait_seconds,
                "lane_wall_seconds": time.perf_counter() - self.started,
                "dtype": "complex128", "physical_dimension": s.n,
                "all_illuminations_per_LU": s.ns}
