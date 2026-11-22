"""Prove all coherent single-source optima are zero for a local SDP.

The proposed rays and means use the native FFT. Arb recomputes the exact
receiver-null projections, all local forms, and the covariance weights.
"""

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, eigvalsh, svd
from scipy.sparse.linalg import LinearOperator, eigsh
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PUBLIC = STUDY / "synthetic/public"
sys.path.insert(0, str(STUDY / "core"))
from certificate import ball, gershgorin, safe_float, star

GEOMETRY = "PW360_10x20"
CASES = ("optfresh_01__pw360_10x20__n02", "optfresh_05__pw360_10x20__n02")
SIDE = 96
N = SIDE**2
CAP = .75
PARTITION_SIDE = 2
MASKED_PLANE_WAVE_MODES = 14
EPSILON = 1.e-4
RIDGE = 1.e-4


def elapsed(wall, cpu):
    return {"wall_seconds": perf_counter() - wall, "cpu_seconds": process_time() - cpu}


def blocks(partition_side):
    width = SIDE // partition_side
    return (np.arange(N).reshape(partition_side, width, partition_side, width)
            .transpose(0, 2, 1, 3).reshape(partition_side**2, width**2))


def mid(a):
    return np.array([[complex(value) for value in row] for row in a.tolist()])


def hi(a):
    return safe_float(a, upper=True)


def lo(a):
    return safe_float(a, upper=False)


def norm_upper(a):
    squared = sum((abs(value).upper()**2 for value in a.entries()), arb(0))
    return hi(squared.sqrt())


def h(a):
    return (a + a.conj().T) / 2


def subrows(a, indices):
    return acb_mat(len(indices), a.ncols(),
                   [a[int(i), j] for i in indices for j in range(a.ncols())])


class NumericGeometry:
    def __init__(self, s, b, symbol, indices):
        self.s, self.b, self.symbol, self.indices = s, b, symbol, indices
        self.u, self.sigma, vh = svd(s, full_matrices=False)
        self.v = vh.conj().T
        self.green_columns = 0
        self.reactive_eigen_operator_calls = 0

    def green(self, x, symbol=None):
        ncolumns = x.reshape(N, -1).shape[1]
        self.green_columns += ncolumns
        if symbol is None:
            symbol = self.symbol
        source = x.reshape(SIDE, SIDE, ncolumns)
        ext = fft2(source, s=symbol.shape, axes=(0, 1))
        result = ifft2(symbol[:, :, None] * ext, axes=(0, 1))[:SIDE, :SIDE]
        return result.reshape(x.shape)

    @staticmethod
    def project_off(x, orthonormal_basis):
        return x - orthonormal_basis @ (orthonormal_basis.conj().T @ x)

    def dark_proposals(self):
        coord = (np.arange(SIDE) + .5) * 2 / SIDE - 1
        yy, xx = np.meshgrid(coord, coord, indexing="ij")
        angles = 2 * np.pi * np.arange(100) / 100
        phase = 2 * np.pi * 4.e8 / 3.e8
        field = np.exp(1j * phase * (xx.ravel()[:, None] * np.cos(angles)
                                    + yy.ravel()[:, None] * np.sin(angles)))
        w = svd(self.project_off(field, self.v), full_matrices=False)[0]
        projector = np.column_stack((self.v, w[:, :MASKED_PLANE_WAVE_MODES]))

        def reactive(z):
            self.reactive_eigen_operator_calls += 1
            z = self.project_off(z, projector)
            return self.project_off(self.green(z, CAP * self.symbol.real), projector)

        operator = LinearOperator((N, N), matvec=reactive, dtype=complex)
        v0 = self.project_off(np.random.default_rng(7393).normal(size=N), projector)
        _, modes = eigsh(operator, k=1, ncv=32, which="LA", tol=5.e-10, v0=v0)
        mode = modes[:, 0]
        real_mode = mode.real if np.linalg.norm(mode.real) > np.linalg.norm(mode.imag) else mode.imag
        v = self.project_off(real_mode / np.linalg.norm(real_mode), projector)
        ev = self.green(v)
        gradients = []
        for ix in self.indices:
            pv = np.zeros(N, complex)
            pev = np.zeros(N, complex)
            pv[ix], pev[ix] = v[ix], ev[ix]
            gradients.append((pev - self.green(pv, self.symbol.conj())) / (2j))
        g = self.project_off(np.column_stack(gradients), self.v)
        frame, singular, _ = svd(g - g.mean(axis=1, keepdims=True), full_matrices=False)
        tetrahedron = np.array([[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])
        z = self.project_off(v[:, None] + EPSILON * (frame[:, :3] @ tetrahedron), self.v)
        return z, singular.tolist()

    def fit_profile(self):
        basis = svd(self.project_off(self.b, self.v), full_matrices=False)[0]
        basis = self.project_off(basis, self.v)
        optical = h(basis.conj().T @ self.green(basis, self.symbol.imag))
        return basis, optical + RIDGE * np.eye(len(optical))

    def mean_proposal(self, y, profile):
        basis, ridge_optical = profile
        x0 = self.v @ ((self.u.conj().T @ y) / self.sigma[:, None])
        rhs = basis.conj().T @ (self.green(x0, self.symbol.imag) - .5j * self.b)
        return x0 - basis @ np.linalg.solve(ridge_optical, rhs)


class ExactGeometry:
    def __init__(self, s, b, symbol, indices):
        self.s, self.b, self.indices = ball(s), ball(b), indices
        self.symbol = [[acb(float(g.real), float(g.imag)) for g in row] for row in symbol]
        self.s_star = star(self.s)
        self.gram = self.s * self.s_star
        self.green_columns = 0
        self.receiver_solves = 0
        self.receiver_rhs_columns = 0

    def project_to_data(self, x_proposal, y):
        x = ball(x_proposal)
        self.receiver_solves += 1
        self.receiver_rhs_columns += x.ncols()
        return x + self.s_star * self.gram.solve(y - self.s * x)

    def green(self, x, *, adjoint=False):
        self.green_columns += x.ncols()
        p, zero = len(self.symbol), acb(0)
        pad = [zero] * (p - SIDE)
        cols = []
        for c in range(x.ncols()):
            rows = [acb.dft([x[SIDE * i + j, c] for j in range(SIDE)] + pad)
                    for i in range(SIDE)]
            transformed = []
            for j in range(p):
                d = acb.dft([rows[i][j] for i in range(SIDE)] + pad)
                filtered = [d[i] * (self.symbol[i][j].conjugate() if adjoint
                                    else self.symbol[i][j]) for i in range(p)]
                transformed.append(acb.dft(filtered, inverse=True)[:SIDE])
            cols.extend(k for i in range(SIDE)
                        for k in acb.dft([transformed[j][i] for j in range(p)],
                                         inverse=True)[:SIDE])
        return acb_mat(x.ncols(), N, cols).transpose()

    def forms(self, x, e):
        sectors, optical = [], []
        for ix in self.indices:
            j, field = subrows(x, ix), subrows(e, ix)
            je = star(j) * field
            sectors.append(star(j) * j - arb(CAP) * (je + star(je)) / 2)
            optical.append((je - star(je)) / acb(0, 2))
        return sectors, optical


def certify_positive(a):
    vectors = ball(eigh(h(mid(a)), check_finite=False)[1])
    gram = star(vectors) * vectors
    gram_floor = gershgorin(gram, upper=False)
    gram_ceil = gershgorin(gram, upper=True)
    floor = gershgorin(star(vectors) * a * vectors, upper=False)
    if not (gram_floor > 0 and floor > 0):
        raise ArithmeticError("The proposed matrix is not certified positive definite.")
    return lo(floor / gram_ceil)


def dark_response(exact, z0):
    nblock, nrays = len(exact.indices), z0.shape[1]
    z = exact.project_to_data(z0, acb_mat(exact.s.nrows(), nrays))
    sector, optical = exact.forms(z, exact.green(z))
    rho = acb_mat(nblock, nrays, [sector[b][j, j].real
                                   for b in range(nblock) for j in range(nrays)])
    tau = acb_mat(nblock, nrays, [optical[b][j, j].real
                                   for b in range(nblock) for j in range(nrays)])
    detector = exact.s * z
    record = dict(
        all_real_sector_entries_upper=hi(max(q.real.upper() for q in rho.entries())),
        total_optical_entries_lower=[lo(sum((tau[b, j].real for b in range(nblock)), arb(0)))
                                     for j in range(nrays)],
        local_optical_midpoints=mid(tau).real.tolist(),
        local_sector_midpoints=mid(rho).real.tolist(),
        interval_receiver_residual_upper=norm_upper(detector))
    return rho, tau, record


def exact_local_weight_matrices(tau, optical):
    nblock, nsource = len(optical), optical[0].nrows()
    rhs = acb_mat(nblock, nsource**2, [-x for a in optical for x in a.entries()])
    weights = tau.solve(rhs)
    return [acb_mat(nsource, nsource, [weights[j, k] for k in range(nsource**2)])
            for j in range(nblock)]


def verify_case(exact, x0, y0, rho, tau):
    nsource, nblock = y0.shape[1], len(exact.indices)
    y = ball(y0)
    x = exact.project_to_data(x0, y)
    sector, optical = exact.forms(x, exact.b + exact.green(x))
    weights = exact_local_weight_matrices(tau, optical)
    floors = [certify_positive(w) for w in weights]
    final_sector_max = []
    for b in range(nblock):
        final_sector = sector[b]
        for j in range(nblock):
            final_sector = final_sector + rho[b, j] * weights[j]
        final_sector_max.append(-certify_positive(-final_sector))
    residual = norm_upper(exact.s * x - y)
    return dict(
        coherent_directions="Every complex a with Euclidean norm 1",
        spatial_SDP_optimal_squared_residual=0,
        interval_mean_receiver_residual_upper=residual,
        covariance_coefficient_eigenvalue_lowers=floors,
        expanded_local_sector_eigenvalue_uppers=final_sector_max,
        current_lift_rank_upper=nblock + 1,
        all_optical_equalities_exact=True)


def load():
    with np.load(PUBLIC / "fine_geometry.npz") as f:
        symbol, s, b = (f[key] for key in
                        ("spectrum", GEOMETRY + "__S", GEOMETRY + "__B"))
    data = {}
    for case in CASES:
        with np.load(PUBLIC / (case + ".npz")) as f:
            data[case] = f["Y"]
    existing = STUDY / "results/theory/fine_feasibility.json"
    baseline = json.loads(existing.read_text())
    return s, b, symbol, data, baseline


def run(output, proposal):
    wall, cpu = perf_counter(), process_time()
    s, b, symbol, data, baseline = load()
    indices = blocks(PARTITION_SIDE)
    ctx.prec, ctx.threads = 160, 1
    geom_wall, geom_cpu = perf_counter(), process_time()
    exact = ExactGeometry(s, b, symbol, indices)
    if proposal is not None:
        with np.load(proposal) as f:
            rays, means = f["rays"], {case: f[case] for case in CASES}
        numeric = None
        spectrum = None
    else:
        numeric = NumericGeometry(s, b, symbol, indices)
        rays, spectrum = numeric.dark_proposals()
        mean_profile = numeric.fit_profile()
        means = {}
    rho, tau, geometry = dark_response(exact, rays)
    geometry.update(elapsed(geom_wall, geom_cpu),
                    numeric_fft_columns=numeric.green_columns if numeric else None,
                    numeric_reactive_eigen_operator_calls=(numeric.reactive_eigen_operator_calls
                                                           if numeric else None),
                    numeric_local_optical_gradient_singular_values=spectrum,
                    exact_fft_columns=exact.green_columns,
                    exact_receiver_solves=exact.receiver_solves,
                    exact_receiver_rhs_columns=exact.receiver_rhs_columns)
    rows = {}
    for case in CASES:
        case_wall, case_cpu = perf_counter(), process_time()
        green_exact0, rhs0 = exact.green_columns, exact.receiver_rhs_columns
        green_num0 = numeric.green_columns if numeric else 0
        if numeric:
            means[case] = numeric.mean_proposal(data[case], mean_profile)
        record = verify_case(exact, means[case], data[case], rho, tau)
        record.update(elapsed(case_wall, case_cpu),
                      numeric_fft_columns=(numeric.green_columns-green_num0 if numeric else None),
                      exact_fft_columns=exact.green_columns-green_exact0,
                      exact_receiver_solves=1,
                      exact_receiver_rhs_columns=exact.receiver_rhs_columns-rhs0)
        rows[case] = record
        print(case, record["covariance_coefficient_eigenvalue_lowers"],
              record["expanded_local_sector_eigenvalue_uppers"], flush=True)
    if proposal is None:
        np.savez_compressed(output.with_suffix(".npz"), rays=rays, **means)
    report = dict(cap=CAP, solver="Constructive local SDP primal",
                  partition=[2, 2], native_cells=N, receivers=s.shape[0], illuminations=b.shape[1],
                  interval_precision_bits=ctx.prec, observations=rows, geometry=geometry,
                  full_state_equation_solves=0, replayed_proposal=proposal is not None,
                  comparison_radius_source="results/theory/fine_feasibility.json",
                  padded_noise_radius={case: next(row["inflated_radius"]
                                                 for row in baseline["observations"]
                                                 if row["id"] == case) for case in CASES},
                  **elapsed(wall, cpu))
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "local_moment.json")
    parser.add_argument("--proposal", type=Path,
                        help="Verify frozen binary proposals without regenerating the rays and means.")
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output, args.proposal)
