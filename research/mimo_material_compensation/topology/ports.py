"""Fixed-termination 24-port scalar-TM response and projection utilities."""
from collections import Counter
from itertools import combinations
from pathlib import Path
from time import perf_counter, process_time
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True

import numpy as np
from scipy.linalg import lu_factor, lu_solve

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent.parent / "z_alg_optical_inertia" / "blind_uq_mimo_theory"
sys.path.insert(0, str(PRIOR / "ports"))
from model import DiskGreen, Green, physical_solve, power_report

P = 24
SWITCHABLE = np.arange(16)
A = np.arange(8)
B = np.arange(0, 16, 2)
TX = np.array([(0, *tail) for tail in combinations(range(1, 16), 7)])
RX = np.array([np.setdiff1d(SWITCHABLE, t) for t in TX])
EDGE_I, EDGE_J = np.triu_indices(16, 1)
MEMBERS = np.any(TX[:, :, None] == SWITCHABLE[None, None, :], axis=1)
CUT_EDGES = MEMBERS[:, EDGE_I] ^ MEMBERS[:, EDGE_J]
INCIDENCE = CUT_EDGES.astype(float)
WATT = 1 / 8


def now():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_s=perf_counter()-start[0], cpu_s=process_time()-start[1])


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
    tmp.replace(path)


def save(path, **value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as fp:
        np.savez_compressed(fp, **value)
    tmp.replace(path)


def load(path, keys=None):
    with np.load(path, allow_pickle=False) as z:
        return {key: z[key] for key in (keys if keys is not None else z.files)}


def cut_id(tx):
    tx = np.asarray(tx)
    if 0 not in tx:
        tx = np.setdiff1d(SWITCHABLE, tx)
    return int(np.flatnonzero((TX == np.sort(tx)).all(axis=1))[0])


def edges(z):
    return z[EDGE_I, EDGE_J, ...]


def project(z, tx, rx=None, power=WATT):
    tx = np.asarray(tx)
    if rx is None:
        rx = np.setdiff1d(SWITCHABLE, tx)
    return z[np.asarray(rx)[:, None], tx[None, :], ...] * np.sqrt(power)


def realify(z):
    return np.concatenate((z.real, z.imag), axis=0)


class PortEngine:
    def __init__(self, side=32, restrict_zero_cells=True):
        self.side = side
        self.restrict_zero_cells = restrict_zero_cells
        self.op = load(PRIOR / "ports" / "public" / "operators" / f"n{side}_f400_a08.npz")
        self.n = side**2
        self.green = Green(self.op)
        self.g = self.op.get("G")
        if self.g is not None and not self.g.size:
            self.g = None
        self.z0 = self.op["M_nominal"]
        self.zref = float(self.op["reference_ohms"])
        self.eye = np.eye(P)
        self.z0_inv = np.linalg.solve(self.z0, self.eye)
        self.drive = 2*np.sqrt(self.zref)*self.eye
        self.b = self.op["E0_bare"] @ (self.z0_inv @ self.drive)
        self.c = -1j*np.sqrt(self.zref)*(self.z0_inv @ self.op["S_bare"])
        self.empty = self.eye - 2*self.zref*self.z0_inv
        self.counters = Counter()

    def state(self, x, jacobian=False, columns=None):
        self.counters["state_solves"] += 1
        if columns is None:
            columns = np.arange(P)
        b = self.b[:, columns]
        self.counters["source_rhs"] += len(columns)
        rhs = x[:, None]*b
        active = np.flatnonzero(x) if self.restrict_zero_cells else np.arange(self.n)
        if self.g is None:
            j, cost = physical_solve(x, self.green.apply, b, self.green.diagonal)
            self.counters["gmres_iterations"] += sum(cost["iterations"])
            factor = None
        else:
            # A zero contrast forces its current exactly to zero regardless of
            # the other pixels. The active principal system is algebraically exact.
            j = np.zeros_like(rhs)
            factor = None
            if active.size:
                g = self.g if active.size == self.n else self.g[np.ix_(active,active)]
                a = np.asfortranarray(-x[active,None]*g)
                a.flat[::len(active)+1] += 1
                factor = lu_factor(a, check_finite=False, overwrite_a=True)
                self.counters["dense_LU"] += 1
                self.counters["sum_LU_active_cells"] += len(active)
                j[active] = lu_solve(factor, rhs[active], check_finite=False)
        # The FFT application evaluates the fixed physical equation independently of the dense LU.
        field = b + self.green.apply(j)
        state_defect = float(np.linalg.norm(j-x[:, None]*field)/max(np.linalg.norm(rhs), 1e-300))
        if not np.isfinite(state_defect) or state_defect > 1e-9:
            raise ArithmeticError(f"Fixed Green state defect {state_defect:g}")
        out = dict(x=x, j=j, field=field, factor=factor, columns=columns, delta=self.c@j,
                   state_defect=state_defect, active=active)
        if jacobian:
            out["jacobian"] = self.jacobian(out)
        return out

    def jacobian(self, st):
        self.counters["Jacobian_builds"] += 1
        # Reciprocity: C is a fixed scalar multiple of B transpose for this circuit.
        # A separate transpose solve is used on the native dense grid, so the derivative
        # is of the stored B and C exactly and does not assume their rounded proportionality.
        if st["factor"] is not None:
            self.counters["transpose_rhs"] += self.c.shape[0]
            active = st["active"]
            ca = lu_solve(st["factor"], self.c[:,active].T, trans=1, check_finite=False).T
            if len(active) == self.n:
                c_loaded = ca
            else:
                c_loaded = self.c + (ca*st["x"][active]) @ self.g[active]
                c_loaded[:, active] = ca
        elif self.g is not None:
            c_loaded = self.c
        else:
            if not np.array_equal(st["columns"], np.arange(P)):
                raise ValueError("FFT Jacobian needs all reciprocal fields")
            c_loaded = self.reciprocal_scale() * st["field"].T
        return np.einsum("ri,ij->rji", c_loaded, st["field"], optimize=True)

    def jvp(self, st, real_direction):
        return self.jacobian(st) @ np.asarray(real_direction, dtype=float)

    def adjoint(self, st, complex_dual):
        return np.einsum("rji,rj->i",self.jacobian(st).conj(),complex_dual).real

    def reciprocal_scale(self):
        factor = -0.5j * float(self.op["pulse"])
        if np.linalg.norm(self.c-factor*self.b.T)/np.linalg.norm(self.c) > 1e-12:
            raise ArithmeticError("B and C are inconsistent under reciprocal scaling")
        return factor

    def full_response(self, x, circuit_state=None):
        if circuit_state is not None and circuit_state!="fixed_matched_50ohm":
            raise ValueError("Projection requires the original fixed matched 50 Ω port loads")
        start = now()
        st = self.state(x)
        delta = st["delta"]
        total = self.empty + delta
        current = (self.eye - total)/np.sqrt(self.zref)
        # Loaded impedance independent of column-order conventions, reconstructed from current.
        antenna = np.linalg.solve(current.T, self.drive.T).T - self.zref*self.eye
        power = power_report(current, antenna, self.drive)
        power["dissipated_plus_radiated_W"] = power["antenna_W"]
        power["ohmic_dissipated_W"] = (8*np.sum(abs(current)**2, axis=0)).tolist()
        power["radiated_W"] = (np.array(power["antenna_W"])-power["ohmic_dissipated_W"]).tolist()
        return dict(total=total, empty=self.empty, delta=delta, current=current, antenna=antenna,
                    power=power, state_defect=st["state_defect"], timing=elapsed(start))

    def independent_response(self, x):
        """Bare VIE and loaded circuit, no dressed Green solve."""
        disk = DiskGreen(self.side, 400)
        action = lambda a: disk.kernel.apply(a[:, None])[:, 0] if a.ndim == 1 else disk.kernel.apply(a)
        dense = None
        if self.side == 32:
            if not hasattr(self, "_bare"):
                self._bare = disk.dense()
            dense = self._bare
        j, cost = physical_solve(x, action, self.op["E0_bare"], disk.kernel.diagonal, dense)
        antenna = self.op["Z_AA"] - 1j*(self.op["S_bare"]@j)
        total = self.eye - 2*self.zref*np.linalg.solve(antenna + self.zref*self.eye, self.eye)
        return total-self.empty, total, cost
