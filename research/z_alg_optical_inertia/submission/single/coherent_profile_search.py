"""Generate 4x4 dual candidates using a native local Galerkin basis.

The Galerkin restriction proposes weights and search directions. Run the
full native diagnostic and the interval verifier before making exclusions.
"""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from numpy.polynomial.legendre import legvander
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from coherent_dual import LABELS
from local_moment import CAP, N, SIDE, NumericGeometry, blocks, elapsed, h, load

MU = 3.5e-4
DEGREE = 6


def compress(geometry):
    width = SIDE // 4
    nlocal = DEGREE**2
    positions = (np.arange(width)+.5)*2/width-1
    basis1 = np.linalg.qr(legvander(positions, DEGREE-1))[0]
    basis = np.kron(basis1, basis1)
    nsmall = nlocal * len(geometry.indices)
    projection = np.zeros((N, nsmall))
    for i, ix in enumerate(geometry.indices):
        projection[ix, i*nlocal:(i+1)*nlocal] = basis
    green = np.empty((nsmall, nsmall), complex)
    for start in range(0, nsmall, 16):
        cols = slice(start, start+16)
        applied = geometry.green(projection[:, cols])
        for i, where in enumerate(geometry.indices):
            green[i*nlocal:(i+1)*nlocal, cols] = basis.T @ applied[where]
    return geometry.s@projection, projection.T@geometry.b, green


class CandidateObjective:
    def __init__(self, s, b, g, y, gap):
        self.b, self.g, self.gap = b, g, gap
        self.ss, self.sy, self.obs_gram = s.conj().T@s/MU, s.conj().T@y/MU, y.conj().T@y
        self.cost = dict(matrix_evaluations=0, projected_minimum_eigensolves=0,
                         projected_cholesky_factorizations=0, projected_stationary_columns=0,
                         source_eigensolves=0)
        self.last = None

    @staticmethod
    def local_forms(x, field):
        je = (x.conj()*field).reshape(16, DEGREE**2).sum(axis=1)
        q_real = (abs(x)**2).reshape(16, DEGREE**2).sum(axis=1) - CAP*je.real
        return np.r_[q_real, je.imag]

    def evaluate(self, p):
        if self.last is not None and np.array_equal(self.last, p):
            return
        self.last = p.copy()
        self.cost["matrix_evaluations"] += 1
        weights = np.repeat(p[:16], DEGREE**2)
        z = np.repeat((CAP*p[:16]+1j*p[16:])/2, DEGREE**2)
        k = self.ss-h(2*z[:, None]*self.g)
        k.flat[::len(k)+1] += weights
        eig, v = eigh(k, subset_by_index=(0, 0), check_finite=False)
        self.cost["projected_minimum_eigensolves"] += 1
        v = v[:, 0]
        self.floor = float(eig[0]-self.gap)
        self.floor_jac = self.local_forms(v, self.g@v)
        if eig[0] < .0001:
            self.obj = .03+max(0., .0001-eig[0])
            self.obj_jac = -self.floor_jac
            return
        l = self.sy+z[:, None]*self.b
        x = cho_solve(cho_factor(k, check_finite=False, lower=True), l, check_finite=False)
        self.cost["projected_cholesky_factorizations"] += 1
        self.cost["projected_stationary_columns"] += l.shape[1]
        elimination = h(self.obs_gram-MU*l.conj().T@x)
        last = len(elimination)-1
        de, modes = eigh(elimination, subset_by_index=(last, last), check_finite=False)
        self.cost["source_eigensolves"] += 1
        a = modes[:, 0]
        current = x@a
        self.obj = -float(de[0])
        self.obj_jac = -MU*self.local_forms(current, self.b@a+self.g@current)

    def objective(self, p):
        self.evaluate(p)
        return self.obj

    def gradient(self, p):
        self.evaluate(p)
        return self.obj_jac

    def gap_constraint(self, p):
        self.evaluate(p)
        return self.floor

    def gap_gradient(self, p):
        self.evaluate(p)
        return self.floor_jac


def run(output, profile_output):
    start, cpu = perf_counter(), process_time()
    s, b, symbol, data, _ = load()
    geometry = NumericGeometry(s, b, symbol, blocks(4))
    s, b, g = compress(geometry)
    setup = dict(**elapsed(start, cpu), native_fft_columns=geometry.green_columns)
    rows, profiles = [], []
    initial = np.r_[np.asarray([1.1, 1.2, .68])[LABELS].ravel(), np.full(16, .1)]
    for case, y in reversed(list(data.items())):
        for gap in (.02, .05):
            row_wall, row_cpu = perf_counter(), process_time()
            objective = CandidateObjective(s, b, g, y, gap)
            result = minimize(objective.objective, initial, jac=objective.gradient, method="SLSQP",
                              bounds=[(.03, 4.)]*16 + [(-.7, 1.8)]*16,
                              constraints=[dict(type="ineq", fun=objective.gap_constraint,
                                                jac=objective.gap_gradient)],
                              options=dict(maxiter=100, ftol=2.e-9))
            w, t = np.round(result.x, 3).reshape(2, 4, 4)
            profiles.append(dict(mu=MU, sector=w.tolist(), optical=t.tolist(),
                                 label=case+" "+str(gap)))
            rows.append(dict(case=case, target_Galerkin_eigenfloor=gap, success=bool(result.success),
                             Galerkin_sq_residual_numeric=float(-result.fun),
                             minimum_Galerkin_eigenvalue_numeric=gap+objective.gap_constraint(result.x),
                             optimization_message=result.message, **objective.cost,
                             **elapsed(row_wall, row_cpu)))
    report = dict(partition=[4, 4], native_cells=N, projected_dimension=len(g), setup=setup,
                  full_state_equation_solves=0, rows=rows, **elapsed(start, cpu))
    output.write_text(json.dumps(report, indent=2)+"\n")
    profile_output.write_text(json.dumps(profiles, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output, args.profiles)
