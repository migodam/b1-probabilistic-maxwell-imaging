"""Blind public input, fixed dressed Green, and full-wave forward/adjoint."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from time import perf_counter, process_time
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
WORKSPACE = HERE.parents[3]
UPSTREAM = WORKSPACE / "research/delegated/max_struct_r1/r8_worktree"
if not UPSTREAM.is_dir():
    UPSTREAM = STUDY / "vendor_r8"
sys.path.insert(0, str(UPSTREAM))

from research.research_fbg_r3_experiments.src.cie_tsom_baseline import _readout
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
from research.research_z_alg_r4_diversity.src.physics import Regularizer

FREEZE = json.loads((HERE / "FREEZE.json").read_text())
CASE_DIR = STUDY / "ports/public/cases"
OPERATOR_DIR = STUDY / "ports/public/operators"
OUT = HERE / "results"


def now():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_s=perf_counter() - start[0], cpu_s=process_time() - start[1])


def write_json(path, value):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    tmp.replace(path)


def write_arrays(path, **arrays):
    tmp = path.with_suffix(".tmp")
    with tmp.open("wb") as f:
        np.savez_compressed(f, **arrays)
    tmp.replace(path)


def col_energy(z):
    return np.einsum("ij,ij->j", z.conj(), z).real


def energy(z):
    return float(np.vdot(z, z).real)


class Green:
    """Finite disk-pulse Toeplitz FFT plus exact small port-network correction."""

    def __init__(self, public, side, counts):
        self.side, self.counts = side, counts
        self.spectrum = public["spectrum"]
        self.left = public["low_rank_left"]
        self.right = public["low_rank_right"]
        self.left_h = self.left.conj().T.copy()
        self.right_h = self.right.conj().T.copy()
        self.diagonal = public["diagonal_dressed"]
        self.dense = public["G"] if "G" in public and public["G"].size else None
        if self.left.shape != (side**2, self.right.shape[0]) or self.right.shape[1] != side**2:
            raise ValueError("Incompatible low-rank dressed Green")

    def _apply(self, value, adjoint):
        n = self.side
        view = value.reshape(n, n, -1)
        symbol = self.spectrum.conj() if adjoint else self.spectrum
        transform = fft2(view, s=symbol.shape, axes=(0, 1), workers=1)
        filtered = ifft2(symbol[..., None] * transform, axes=(0, 1), workers=1)
        convolution = filtered[:n, :n].reshape(value.shape)
        self.counts["green_adjoint" if adjoint else "green_forward"] += 1
        self.counts["green_fft_2d"] += 2 * value.shape[1]
        return convolution + (self.right_h @ (self.left_h @ value) if adjoint else
                              self.left @ (self.right @ value))

    def apply(self, value):
        return self._apply(value, False)

    def adjoint(self, value):
        return self._apply(value, True)


class Model:
    def __init__(self, public, y, radius, operator_id, side):
        self.side, self.shape, self.n = side, (side, side), side**2
        self.S = public["S"]
        self.B = public["E0"]
        self.SH = self.S.conj().T.copy()
        self.Y, self.radius, self.operator_id = y, radius, operator_id
        self.ns = self.B.shape[1]
        if self.S.shape[0] != y.shape[0] or self.B.shape != (self.n, y.shape[1]):
            raise ValueError("Training data and public operators have incompatible shape")
        self.y_scale = energy(y)
        self.cost = Counter()
        self.green = Green(public, side, self.cost)
        self.reg = Regularizer(self.shape, smooth=1e-4, tv=1e-4, epsilon=.01)
        self.state = None
        self.factor = None
        self.factor_x = None
        self.last_solutions = {}
        self.coarse = None
        self.coarse_x = None
        self.coarse_factor = None

    @classmethod
    def read(cls, case):
        with np.load(CASE_DIR / f"{case}.npz", allow_pickle=False) as src:
            y, radius, operator_id = src["Y"], float(src["radius"]), str(src["operator_id"].item())
        with np.load(OPERATOR_DIR / f"{operator_id}.npz", allow_pickle=False) as op:
            side = int(round(np.sqrt(op["E0"].shape[0])))
            return cls(op, y, radius, operator_id, side)

    def reset(self):
        self.state = None
        self.factor = None
        self.factor_x = None
        self.last_solutions.clear()
        self.coarse = None
        self.coarse_x = None
        self.coarse_factor = None
        self.cost.clear()

    def coarse_preconditioner(self, x, adjoint):
        started = now()
        if self.coarse is None:
            width = 6
            axes = np.r_[np.arange(width), np.arange(self.side-width, self.side)]
            index = tuple(np.array(np.meshgrid(axes, axes, indexing="ij")).reshape(2,-1))
            rank = len(index[0])
            packed = np.zeros(self.shape + (rank,), complex)
            packed[index[0], index[1], np.arange(rank)] = 1.
            fourier_columns = ifft2(packed, axes=(0,1), norm="ortho", workers=1).reshape(self.n,-1)
            self.coarse = index, self.green.apply(fourier_columns), np.eye(rank, dtype=complex)
            self.cost["native_coarse_basis_builds"] += 1
        index, g_coarse, identity = self.coarse
        if self.coarse_x is None or not np.array_equal(self.coarse_x, x):
            restrict = fft2((x[:,None]*g_coarse).reshape(self.shape + (-1,)),
                            axes=(0,1), norm="ortho", workers=1)[index]
            self.coarse_factor = lu_factor(identity-restrict, check_finite=False)
            self.coarse_x = x.copy()
            self.cost["native_coarse_factors"] += 1
        fact = self.coarse_factor

        def action(v):
            fourier = fft2(v.reshape(self.shape), norm="ortho", workers=1)[index]
            lift = lu_solve(fact, fourier, trans=2 if adjoint else 0, check_finite=False)-fourier
            packed = np.zeros(self.shape, complex)
            packed[index] = lift
            self.cost["native_coarse_preconditioner_calls"] += 1
            return v + ifft2(packed, norm="ortho", workers=1).ravel()

        for key,value in elapsed(started).items():
            self.cost["native_coarse_setup_"+key] += value
        return LinearOperator((self.n,self.n), matvec=action, dtype=complex)

    def solve(self, x, b, adjoint=False, rtol=None):
        tol = FREEZE["baseline"]["physical_rtol"] if rtol is None else rtol
        mode = "physical_adjoint" if adjoint else "physical_forward"
        self.cost[mode + "_rhs"] += b.shape[1]
        if self.green.dense is not None:
            if not adjoint or self.factor is None or not np.array_equal(self.factor_x, x):
                a = -x[:, None] * self.green.dense
                a.flat[::self.n + 1] += 1
                self.factor = lu_factor(a, check_finite=False)
                self.factor_x = x.copy()
                self.cost["physical_lu"] += 1
            answer = lu_solve(self.factor, b, trans=2 if adjoint else 0, check_finite=False)
        else:
            answer = self.solve_fft(x, b, adjoint, mode, tol)
        if adjoint:
            defect = answer - self.green.adjoint(x[:, None] * answer) - b
        else:
            defect = answer - x[:, None] * self.green.apply(answer) - b
        relative = np.sqrt(col_energy(defect) / np.maximum(col_energy(b), 1e-300))
        self.cost["max_" + mode + "_relative_defect"] = max(
            self.cost["max_" + mode + "_relative_defect"], float(relative.max()))
        if not np.all(np.isfinite(answer)) or relative.max() > 5 * tol:
            raise ArithmeticError(f"{mode} original state relative defect = {relative.max():.4g}")
        return answer

    def solve_fft(self, x, b, adjoint, mode, tol):
        def apply(v):
            column = v[:, None]
            if adjoint:
                return v - self.green.adjoint(x[:, None] * column)[:, 0]
            return v - x * self.green.apply(column)[:, 0]

        diagonal = 1 - x * (self.green.diagonal.conj() if adjoint else self.green.diagonal)
        fallback = LinearOperator((self.n, self.n), matvec=lambda v: v / diagonal, dtype=complex)
        matrix = LinearOperator((self.n, self.n), matvec=apply, dtype=complex)
        result = np.zeros_like(b)
        if not np.any(b):
            return result
        preconditioner = self.coarse_preconditioner(x, adjoint)
        old = self.last_solutions.get(adjoint, np.zeros_like(b))
        if old.shape != b.shape:
            old = np.zeros_like(b)
        for source in range(b.shape[1]):
            initial = old[:, source]
            if not np.any(b[:, source]):
                continue
            iterations = []

            def callback(value):
                iterations.append(float(value))

            result[:, source], flag = gmres(matrix, b[:, source], M=preconditioner,
                                             x0=initial, rtol=tol * .9, atol=0,
                                             restart=65, maxiter=5, callback=callback,
                                             callback_type="pr_norm")
            if flag:
                self.cost[mode + "_coarse_fallback_rhs"] += 1
                result[:, source], flag = gmres(matrix, b[:, source], M=fallback,
                                                 x0=initial, rtol=tol * .9, atol=0,
                                                 restart=65, maxiter=16, callback=callback,
                                                 callback_type="pr_norm")
            self.cost[mode + "_gmres_iters"] += len(iterations)
            if flag:
                self.cost[mode + "_failed_rhs"] += 1
                raise ArithmeticError(f"{mode} GMRES source {source}: info={flag}")
        self.last_solutions[adjoint] = result.copy()
        return result

    def forward(self, x, fresh=False):
        if not fresh and self.state is not None and np.array_equal(x, self.state["x"]):
            return self.state
        before = now()
        self.cost["physics_evaluations"] += 1
        j = self.solve(x, x[:, None] * self.B)
        field = self.B + self.green.apply(j)
        pred = self.S @ j
        residual = pred - self.Y
        absolute = np.sqrt(energy(residual))
        self.state = {
            "x": x.copy(), "j": j, "field": field, "pred": pred, "residual": residual,
            "abs_residual": absolute, "residual_ratio": absolute / self.radius,
            "relative_residual": absolute / np.sqrt(self.y_scale),
            "state_relative_defect": np.sqrt(energy(j - x[:, None] * field) /
                                             max(energy(x[:, None] * self.B), 1e-300)),
            **elapsed(before)
        }
        return self.state

    def loss(self, x):
        st = self.forward(x)
        p = self.reg.evaluate(x)[0]
        return .5 * energy(st["residual"]) / self.y_scale + p

    def data_grad(self, x):
        st = self.forward(x)
        adjoint = self.solve(x, self.SH @ st["residual"], adjoint=True)
        grad = np.sum((adjoint.conj() * st["field"]).real, axis=1) / self.y_scale
        return .5 * energy(st["residual"]) / self.y_scale, grad

    def audit(self, x, **labels):
        st = self.forward(x)
        return {
            "physical_l2": st["abs_residual"], "relative_residual": st["relative_residual"],
            "residual_over_radius": st["residual_ratio"],
            "relative_state_defect": st["state_relative_defect"],
            "physical_objective": self.loss(x),
            "within_ball_numerical": bool(st["residual_ratio"] < 1),
            **labels
        }
