"""Exact active-support VIE response and complete occupancy Jacobian.

The support is STRICTLY a > 0. No positive occupancy is truncated.  Zero
occupancy cells still have response derivatives and are included in J_occ.
All observations and Jacobian rows use source-major packing, F.T.reshape(-1).
"""
from __future__ import annotations

import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve


class Solver:
    """Exact elimination of rows with identically zero contrast current."""

    name = "EXACT_SUPPORT"

    def __init__(self, C, O, E, Fh, sigma, ledger=None):
        self.C = np.asarray(C, dtype=np.complex128)
        self.O = np.asarray(O, dtype=np.complex128)
        self.E = np.asarray(E, dtype=np.complex128)
        self.Fh = np.asarray(Fh, dtype=np.complex128)
        self.sigma = float(sigma)
        self.ledger = ledger
        self.n = self.C.shape[0]
        self.nrec, self.nsrc = self.O.shape[0], self.E.shape[1]
        if self.C.shape != (self.n, self.n):
            raise ValueError("C must be square")
        if self.O.shape[1] != self.n or self.E.shape[0] != self.n:
            raise ValueError("incompatible port dimensions")
        if self.Fh.shape != (self.nrec, self.nsrc):
            raise ValueError("Fh shape must match receiver x source")
        if not self.sigma > 0:
            raise ValueError("sigma must be positive")
        self.calls = 0

    def _record(self, diagnostics):
        """Accept a simple callable or an append/event style ledger."""
        if self.ledger is None:
            return
        if callable(self.ledger):
            self.ledger(dict(diagnostics))
        elif hasattr(self.ledger, "event"):
            self.ledger.event("EXACT_SUPPORT_EVALUATE", **diagnostics)
        elif hasattr(self.ledger, "append"):
            self.ledger.append(dict(diagnostics))
        else:
            raise TypeError("ledger must be callable or have event/append")

    def evaluate(self, a, delta, need_jac=False, y=None):
        started = time.perf_counter()
        a = np.asarray(a, dtype=np.float64).reshape(-1)
        delta = complex(delta)
        if a.shape != (self.n,) or not np.all(np.isfinite(a)):
            raise ValueError("a must have one finite real occupancy per cell")
        if np.any(a < 0):
            raise ValueError("exact support requires nonnegative occupancy")
        if not np.isfinite(delta):
            raise ValueError("delta must be finite")
        S = np.flatnonzero(a > 0)
        ns = int(S.size)
        nsolve = 0
        n_rhs = 0
        LU = None
        hs = np.zeros((ns, self.nsrc), dtype=np.complex128)
        ts = time.perf_counter()
        if ns:
            aS = a[S]
            CS = self.C[np.ix_(S, S)]
            M = np.eye(ns, dtype=np.complex128) - delta * aS[:, None] * CS
            LU = lu_factor(M, check_finite=False)
            hs = lu_solve(LU, delta * aS[:, None] * self.E[S], check_finite=False)
            nsolve += 1
            n_rhs += self.nsrc
        solve_seconds = time.perf_counter() - ts
        F = self.Fh + self.O[:, S] @ hs
        result = {"F": F}
        jac_seconds = 0.0
        if need_jac:
            tj = time.perf_counter()
            e = self.E + self.C[:, S] @ hs
            if ns:
                # Z = O_S (I - delta A_S C_SS)^(-1).
                # Ordinary transpose is essential: the derivative is holomorphic.
                Z = lu_solve(LU, self.O[:, S].T, trans=1, check_finite=False).T
                # On the active block, O_S + delta Z A_S C_SS = Z.
                # Only inactive columns require the indirect scattering term.
                # This identity is especially useful near full support.
                T = np.empty_like(self.O)
                T[:, S] = Z
                I = np.flatnonzero(a == 0)
                if I.size:
                    T[:, I] = self.O[:, I] + delta * (Z * a[S][None, :]) @ self.C[np.ix_(S, I)]
                nsolve += 1
                n_rhs += self.nrec
            else:
                T = self.O
            J_occ = np.einsum("rj,js->srj", delta * T, e, optimize=True)
            J_occ = J_occ.reshape(self.nsrc * self.nrec, self.n)
            J_delta = (T @ (a[:, None] * e)).T.reshape(-1)
            result.update(J_occ=J_occ, J_delta=J_delta)
            jac_seconds = time.perf_counter() - tj
        self.calls += 1
        diagnostics = dict(
            method=self.name, call=self.calls, support_size=ns,
            factor_order=ns, LU=int(ns > 0), solve_calls=nsolve,
            RHS=n_rhs, RHS_by_order={str(ns): n_rhs} if ns else {},
            forward_RHS=self.nsrc if ns else 0,
            adjoint_RHS=self.nrec if ns and need_jac else 0,
            J_occ_columns=self.n if need_jac else 0,
            inactive_derivative_columns=int(self.n - ns) if need_jac else 0,
            strict_support=True, truncation_threshold=0.0,
            packing="source-major F.T.reshape(-1)",
            solve_seconds=solve_seconds, jacobian_seconds=jac_seconds,
            wall_seconds=time.perf_counter() - started,
            no_full_order_factorization=bool(ns < self.n),
        )
        result["diagnostics"] = diagnostics
        if y is not None:
            residual = (F - np.asarray(y)).T.reshape(-1)
            result["residual"] = residual
            result["phi"] = float(np.vdot(residual, residual).real / self.sigma**2)
        self._record(diagnostics)
        return result
