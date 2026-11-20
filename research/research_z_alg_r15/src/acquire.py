"""Generate the broadband TM differential-scattering experiment."""

import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres

from ndt import Air, BACKGROUND, BOUNDS, FREQUENCIES_GHZ, PLATE_PERMITTIVITY
from ndt import PLATE_SIDE, RECEIVERS, SOURCES, differential, grid, host


COARSE = 24
FINE = 48
QUADRATURE_SIDE = 12
NOISE_SEED = 15110820
REFERENCE_ID = "dev_hole"
REFERENCE_GHZ = 16.0
RESULTS = Path(__file__).resolve().parents[1] / "results"


def ellipse(x, y, a, b, angle_deg, delta):
    return dict(kind="ellipse", x_mm=x, y_mm=y, a_mm=a, b_mm=b,
                angle_deg=angle_deg, delta_epsilon=delta)


def circle(x, y, radius, delta=-1.25):
    return ellipse(x, y, radius, radius, 0.0, delta)


def capsule(x, y, half_length, half_width, angle_deg, delta):
    return dict(kind="capsule", x_mm=x, y_mm=y, a_mm=half_length, b_mm=half_width,
                angle_deg=angle_deg, delta_epsilon=delta)


def pair(x, y, spacing, angle_deg, radius):
    angle = np.deg2rad(angle_deg)
    dx = spacing * np.cos(angle) / 2
    dy = spacing * np.sin(angle) / 2
    return [circle(x - dx, y - dy, radius), circle(x + dx, y + dy, radius)]


CASES = [
    ("dev_hole", "development", [circle(-2.15, 2.7, 1.65)]),
    ("dev_pair", "development", pair(1.3, -1.85, 4.2, 33, 1.25)),
    ("dev_slot", "development", [capsule(2.35, 1.4, 2.7, 0.66, -25, -1.15)]),
    ("test_hole", "independent", [circle(4.1, -3.2, 1.30)]),
    ("test_pair_close", "independent", pair(-0.8, 2.2, 3.0, -24, 1.28)),
    ("test_pair_mid", "independent", pair(-0.8, 2.2, 4.4, -24, 1.28)),
    ("test_pair_wide", "independent", pair(-0.8, 2.2, 6.2, -24, 1.28)),
    ("test_ellipse", "independent", [ellipse(-4.18, -2.7, 2.65, 0.91, 49, -1.15)]),
    ("test_slot", "independent", [capsule(3.75, 3.45, 3.18, 0.62, 63, -1.25)]),
    ("test_mixed", "independent", [circle(-4.25, -1.85, 1.48),
                                    ellipse(2.95, 3.25, 1.67, 1.00, -31, 1.50)]),
    ("test_multi", "independent", [circle(-5.35, 4.30, 1.15),
                                    circle(3.30, 4.05, 1.48),
                                    ellipse(1.25, -4.70, 2.03, 0.80, 17, -1.25)]),
    ("test_empty", "independent", []),
]


def material(case_shapes, n=FINE, samples_per_side=QUADRATURE_SIDE):
    side_mm = 1e3 * PLATE_SIDE / n
    offsets = ((np.arange(samples_per_side) + 0.5) / samples_per_side - 0.5) * side_mm
    offset_x, offset_y = np.meshgrid(offsets, offsets, indexing="xy")
    xy = 1e3 * grid(n)
    sample_x = xy[:, 0, None] + offset_x.ravel()
    sample_y = xy[:, 1, None] + offset_y.ravel()
    contrast = np.zeros_like(sample_x)
    for primitive in case_shapes:
        angle = np.deg2rad(primitive["angle_deg"])
        x = sample_x - primitive["x_mm"]
        y = sample_y - primitive["y_mm"]
        u = np.cos(angle) * x + np.sin(angle) * y
        v = -np.sin(angle) * x + np.cos(angle) * y
        a, b = primitive["a_mm"], primitive["b_mm"]
        if primitive["kind"] == "ellipse":
            inside = (u / a) ** 2 + (v / b) ** 2 <= 1
        else:
            end = np.maximum(abs(u) - (a - b), 0.0)
            inside = end**2 + v**2 <= b**2
        contrast[inside] = primitive["delta_epsilon"]
    return contrast.mean(axis=1).reshape(n, n)


def fft_field(air, apply, material, initial=None):
    m = np.asarray(material).ravel()
    size = m.size
    operator = LinearOperator((size, size), matvec=lambda e: e - apply(m * e), dtype=complex)
    diagonal = 1 - air.self_term * m
    pre = LinearOperator((size, size), matvec=lambda e: e / diagonal, dtype=complex)
    field = np.empty_like(air.E0)
    residuals = []
    iterations = []
    for tx, e0 in enumerate(air.E0.T):
        history = []
        e, info = gmres(operator, e0, x0=None if initial is None else initial[:, tx],
                        M=pre, rtol=1e-11, atol=0, restart=80, maxiter=30,
                        callback=history.append, callback_type="pr_norm")
        residual = relative(operator @ e, e0)
        if info or residual > 1e-9:
            raise RuntimeError(f"GMRES status {info} and relative residual {residual:.3g}")
        field[:, tx] = e
        residuals.append(residual)
        iterations.append(len(history))
    return field, dict(max_relative_residual=max(residuals), krylov_iterations=sum(iterations))


def dense_absolute(air, green, delta, background=BACKGROUND):
    """Direct solve in total polarization current for a numerical comparison."""
    m = background + np.asarray(delta).reshape(-1, 1)
    matrix = np.eye(m.size, dtype=complex) - m * green
    rhs = m * air.E0
    current = lu_solve(lu_factor(matrix, check_finite=False), rhs, check_finite=False)
    return air.S @ current, relative(matrix @ current, rhs)


def relative(value, reference):
    return float(np.linalg.norm(value - reference) / max(np.linalg.norm(reference), 1e-300))


def main():
    wall_start, cpu_start = perf_counter(), process_time()
    fine_truth = np.array([material(shapes) for _, _, shapes in CASES])
    coarse_truth = fine_truth.reshape(len(CASES), COARSE, 2, COARSE, 2).mean(axis=(2, 4))
    if fine_truth.min() < BOUNDS[0] or fine_truth.max() > BOUNDS[1]:
        raise ValueError("Generated truth exceeds the stated material bounds")
    shape = (len(CASES), len(FREQUENCIES_GHZ), len(RECEIVERS), len(SOURCES))
    clean = np.zeros(shape, dtype=complex)
    clean_same = np.zeros_like(clean)
    background = np.empty(shape[1:], dtype=complex)
    background_same = np.empty_like(background)
    diagnostics = []
    for i, f_ghz in enumerate(FREQUENCIES_GHZ):
        f_wall_start, f_cpu_start = perf_counter(), process_time()
        air = Air(FINE, f_ghz)
        apply = air.convolution()
        e_bg, stats = fft_field(air, apply, np.full(FINE**2, BACKGROUND))
        background[i] = air.S @ (BACKGROUND * e_bg)
        coarse_air = Air(COARSE, f_ghz)
        model = host(COARSE, f_ghz)
        background_same[i] = coarse_air.S @ (BACKGROUND * model.Eh)
        n_iter = stats["krylov_iterations"]
        max_residual = stats["max_relative_residual"]
        for j, delta in enumerate(fine_truth):
            if not np.any(delta):
                continue
            m = BACKGROUND + delta.reshape(-1, 1)
            field, stats = fft_field(air, apply, m, e_bg)
            clean[j, i] = air.S @ (m * field) - background[i]
            clean_same[j, i] = differential(model, coarse_truth[j])
            n_iter += stats["krylov_iterations"]
            max_residual = max(max_residual, stats["max_relative_residual"])

        record = dict(frequency_ghz=float(f_ghz), fine_krylov_iterations=n_iter,
                      fine_tx_rhs=(len(CASES)) * len(SOURCES),
                      max_fine_relative_state_residual=max_residual)
        if f_ghz in (2, 16, 30):
            fine_green = air.dense()
            fine_direct, direct_residual = dense_absolute(air, fine_green, fine_truth[0])
            fine_bg_direct, _ = dense_absolute(air, fine_green, np.zeros(FINE**2))
            coarse_green = coarse_air.dense()
            coarse_direct, _ = dense_absolute(coarse_air, coarse_green, coarse_truth[0])
            coarse_bg_direct, _ = dense_absolute(coarse_air, coarse_green, np.zeros(COARSE**2))
            record.update(
                dense_fft_relative_matvec=relative(apply(e_bg[:, 0]), fine_green @ e_bg[:, 0]),
                fine_fft_direct_relative_differential=relative(clean[0, i], fine_direct - fine_bg_direct),
                fine_direct_relative_state_residual=direct_residual,
                coarse_host_direct_relative_differential=relative(clean_same[0, i], coarse_direct - coarse_bg_direct),
            )
        record["elapsed_seconds"] = perf_counter() - f_wall_start
        record["cpu_seconds"] = process_time() - f_cpu_start
        diagnostics.append(record)
        print(f"{f_ghz:5.2f} GHz: {n_iter:5d} Krylov iterations, "
              f"residual {max_residual:.2e}, {record['elapsed_seconds']:.2f}s", flush=True)

    reference_i = int(np.flatnonzero(FREQUENCIES_GHZ == REFERENCE_GHZ)[0])
    reference_j = next(j for j, (case_id, _, _) in enumerate(CASES) if case_id == REFERENCE_ID)
    reference_amplitude = float(np.sqrt(np.mean(np.abs(clean[reference_j, reference_i]) ** 2)))
    rng = np.random.default_rng(NOISE_SEED)
    noise = (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)
    np.savez_compressed(RESULTS / "acquisition.npz", frequencies_ghz=FREQUENCIES_GHZ,
                        case_ids=[case_id for case_id, _, _ in CASES],
                        splits=[split for _, split, _ in CASES], clean=clean,
                        clean_same=clean_same, noise=noise, truth=coarse_truth,
                        truth_fine=fine_truth, reference_amplitude=reference_amplitude,
                        background=background, background_same=background_same,
                        sources_m=SOURCES, receivers_m=RECEIVERS)
    record = dict(
        host_relative_permittivity=PLATE_PERMITTIVITY, side_mm=1e3 * PLATE_SIDE,
        coarse_cells_per_side=COARSE, fine_cells_per_side=FINE,
        fine_cell_geometry_quadrature_samples=QUADRATURE_SIDE**2,
        source_radius_mm=1e3 * np.linalg.norm(SOURCES[0]),
        receiver_radius_mm=1e3 * np.linalg.norm(RECEIVERS[0]),
        frequency_ghz=FREQUENCIES_GHZ.tolist(), material_deviation_bounds=BOUNDS,
        reference_case_id=REFERENCE_ID, reference_frequency_ghz=REFERENCE_GHZ,
        reference_amplitude=reference_amplitude,
        reference_statistic="Real scalar RMS of the 128 complex fine-grid differential readings",
        unit_power_incident_coefficient="(i/4) sqrt(f_GHz/16) times H_0^(1)(k r)",
        noise_draw="numpy PCG64, (N(0,1)+i*N(0,1))/sqrt(2), fixed per case/frequency/receiver/transmitter",
        noise_seed=NOISE_SEED, case_definitions=[dict(id=i, split=s, primitives=p) for i, s, p in CASES],
        elapsed_seconds=perf_counter() - wall_start, cpu_seconds=process_time() - cpu_start,
        relative_cross_grid_differential=[
            None if not case_shapes else [relative(clean_same[j, i], clean[j, i])
            for i in range(len(FREQUENCIES_GHZ))] for j, (_, _, case_shapes) in enumerate(CASES)
        ], diagnostics=diagnostics,
    )
    (RESULTS / "acquisition.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Reference RMS {reference_amplitude:.7g}, total {record['elapsed_seconds']:.2f}s")


if __name__ == "__main__":
    main()
