"""Screen local-moment SDP dual profiles on the full 96x96 states.

Positive Ritz values and stationary solves are numerical proposals. The
interval proof is separate. A numerical value below the radius is a property
of the selected profile, with no implication for the optimized SDP.
"""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, cg, eigsh
from threadpoolctl import threadpool_limits

from local_moment import CAP, HERE, N, SIDE, NumericGeometry, blocks, elapsed, h, load

LABELS = np.array([[0, 1, 1, 0], [1, 2, 2, 1], [1, 2, 2, 1], [0, 1, 1, 0]])
SYMMETRIC_PROFILES = (
    (1.20, 1.26, .97, .10), (1.20, 1.26, .94, .10), (1.20, 1.26, .91, .10),
    (1.10, 1.20, .90, .10), (1.05, 1.05, .85, 0.), (1.10, 1.15, .91, 0.),
    (1.10, 1.15, .86, 0.), (1.02, 1.20, .58, 0.), (1.02, 1.20, .58, .10),
    (1.10, 1.20, .60, .30), (1.10, 1.20, .75, .10), (1.00, 1.00, .60, 0.),
    (1.00, 1.00, .60, .10), (1.00, 1.00, .50, .30),
)


def expand_block(values):
    values = np.asarray(values).reshape(4, 4)
    return values.repeat(SIDE // 4, axis=0).repeat(SIDE // 4, axis=1).ravel()


def candidate_grid():
    for mu in (3.5e-4, 5.5e-4, 1.e-4, 5.e-5, 1.e-5, 2.e-4):
        for corners, edges, center, t in SYMMETRIC_PROFILES:
            yield dict(mu=mu, sector=np.array([corners, edges, center])[LABELS].tolist(),
                       optical=np.full((4, 4), t).tolist())


def verify_stationarity(geometry, observations, proposal, radii):
    wall, cpu = perf_counter(), process_time()
    green0 = geometry.green_columns
    mu = proposal["mu"]
    w = expand_block(proposal["sector"])
    t = expand_block(proposal["optical"])
    z = (CAP * w + 1j * t) / 2
    receiver = geometry.s
    receiver_star = receiver.conj().T

    def action(x):
        return (receiver_star @ (receiver @ x) / mu + w*x - z*geometry.green(x)
                - geometry.green(z.conj()*x, geometry.symbol.conj()))

    def precondition(x):
        factors = geometry.sigma**2 / (mu + geometry.sigma**2)
        return x - geometry.v @ (factors * (geometry.v.conj().T @ x))

    a = LinearOperator((N, N), matvec=action, dtype=complex)
    m = LinearOperator((N, N), matvec=precondition, dtype=complex)
    eigenvalue, vectors = eigsh(a, k=1, which="SA", tol=1.e-10, ncv=50,
                               v0=np.random.default_rng(99223).normal(size=N))
    low = float(eigenvalue[0])
    directional = vectors[:, 0]
    defect = float(np.linalg.norm(action(directional)-low*directional))
    result = dict(**proposal, native_Ritz_min=low, native_Ritz_defect=defect,
                  eigensolve_fft_columns=geometry.green_columns - green0, cases={})
    if low > 0:
        for obs, y in observations.items():
            case_wall, case_cpu = perf_counter(), process_time()
            case_green0 = geometry.green_columns
            rhs = receiver_star @ y / mu + z[:, None] * geometry.b
            solutions, solves = [], []
            for column in rhs.T:
                iterations = []
                x, status = cg(a, column, M=m, rtol=3.e-12, atol=0., maxiter=600,
                               callback=lambda _: iterations.append(1))
                solutions.append(x)
                solves.append(dict(iterations=len(iterations), status=int(status)))
            x = np.column_stack(solutions)
            backward = np.linalg.norm(np.column_stack([action(xx) for xx in x.T])-rhs) / np.linalg.norm(rhs)
            row = dict(relative_stationary_residual=float(backward), cg_solves=solves, radius=radii[obs])
            if backward < 1.e-9 and all(solve["status"] == 0 for solve in solves):
                matrix = h(y.conj().T @ y - mu * rhs.conj().T @ x)
                values, modes = eigh(matrix)
                coherent = modes[:, -1]
                row.update(coherent_residual_floor_numeric=float(np.sqrt(max(0., values[-1]))),
                           original_beam_floor_numeric=float(np.sqrt(max(0., matrix.diagonal().real.max()))),
                           top_coherent_real=coherent.real.tolist(),
                           top_coherent_imaginary=coherent.imag.tolist(),
                           projected_data_norm=float(np.linalg.norm(y @ coherent)),
                           joint_trace_floor_numeric=float(np.sqrt(np.maximum(0., values).sum())),
                           source_form_eigenvalues=values.tolist())
            row.update(elapsed(case_wall, case_cpu),
                       numeric_fft_columns=geometry.green_columns - case_green0)
            result["cases"][obs] = row
    result.update(elapsed(wall, cpu), numeric_fft_columns=geometry.green_columns - green0)
    return result


def run(profiles, output):
    wall, cpu = perf_counter(), process_time()
    s, b, symbol, data, feasibility = load()
    geometry = NumericGeometry(s, b, symbol, blocks(4))
    radii = {row["id"]: row["inflated_radius"] for row in feasibility["observations"] if row["id"] in data}
    setup = elapsed(wall, cpu)
    rows = []
    for index, proposal in enumerate(profiles):
        row = verify_stationarity(geometry, data, proposal, radii)
        row["index"] = index
        floors = [result.get("coherent_residual_floor_numeric") for result in row["cases"].values()]
        print(index, row["mu"], row["native_Ritz_min"], floors, flush=True)
        rows.append(row)
        report = dict(cap=CAP, partition=[4, 4], native_cells=N, setup=setup, rows=rows,
                      full_state_equation_solves=0, numeric_fft_columns=geometry.green_columns,
                      **elapsed(wall, cpu))
        output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, help="JSON array of native dual profile proposals.")
    parser.add_argument("--output", type=Path, default=HERE / "coherent_dual.json")
    args = parser.parse_args()
    profiles = json.loads(args.profiles.read_text()) if args.profiles else candidate_grid()
    with threadpool_limits(1):
        run(profiles, args.output)
