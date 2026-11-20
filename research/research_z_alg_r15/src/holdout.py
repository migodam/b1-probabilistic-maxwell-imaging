"""Generate a seeded, paired lossless/Debye holdout on the R15 acquisition."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

from acquire import COARSE, FINE, RESULTS, capsule, circle, ellipse, material, pair
from acquire import dense_absolute, fft_field, relative
from ndt import Air, FREQUENCIES_GHZ, PLATE_PERMITTIVITY, PLATE_SIDE
from ndt import RECEIVERS, SOURCES, debye_epsilon, dielectric_factor, differential, host, host_debye


SHAPE_SEED = 15112026
NOISE_SEED = 15112027
NOISE_FRACTION = 0.02


def cases(seed=SHAPE_SEED):
    rng = np.random.default_rng(seed)
    records = []

    def record(shape, primitives, **extra):
        records.append(dict(id=f"blind_{len(records):02d}", split="independent",
                            shape=shape, primitives=primitives, **extra))

    def position(limit):
        return (float(v) for v in rng.uniform(-limit, limit, 2))

    def draw(lo, hi):
        return float(rng.uniform(lo, hi))

    def angle():
        return draw(0, 360)

    x, y = position(5.8)
    record("single_air_hole", [circle(x, y, draw(1.1, 1.65))])
    for label, center_limit, dist_bounds in (("close", 4.2, (2.65, 3.15)),
                                             ("wide", 3.6, (5.3, 6.6))):
        x, y = position(center_limit)
        spacing = draw(*dist_bounds)
        record("two_air_holes", pair(x, y, spacing, angle(), draw(1.05, 1.27)),
               spacing_mm=spacing, spacing_class=label)

    x, y = position(5.6)
    record("elliptical_air_hole", [ellipse(x, y, draw(2.2, 3.05), draw(0.7, 1.0), angle(), -1.25)])
    x, y = position(5.4)
    record("thin_slot", [capsule(x, y, draw(2.5, 3.65), draw(0.36, 0.52), angle(), -1.25)])

    x, y = position(3.8)
    theta = np.deg2rad(angle())
    offset = 0.5 * draw(4.4, 6.8) * np.array((np.cos(theta), np.sin(theta)))
    opposite = np.array((x, y))
    hole, load = opposite - offset, opposite + offset
    record("mixed_air_hole_and_dielectric_load", [
        circle(float(hole[0]), float(hole[1]), draw(1.05, 1.5)),
        ellipse(float(load[0]), float(load[1]), draw(1.25, 1.85), draw(0.65, 1.05),
                angle(), draw(0.9, 1.5)),
    ])

    x, y = position(2.7)
    azimuth = angle() + 120 * np.arange(3) + rng.uniform(-14, 14, 3)
    directions = np.column_stack((np.cos(np.deg2rad(azimuth)), np.sin(np.deg2rad(azimuth))))
    offsets = rng.uniform(3.65, 5.2, 3)[:, None] * directions
    marks = offsets + np.array((x, y))
    record("three_air_defects", [
        circle(*map(float, marks[0]), draw(0.8, 1.2)),
        circle(*map(float, marks[1]), draw(0.85, 1.3)),
        ellipse(*map(float, marks[2]), draw(1.3, 1.95), draw(0.58, 0.95), angle(), -1.25),
    ])

    x, y = position(5.3)
    azimuth = angle() + np.array([-50., 12., 86.]) + rng.uniform(-8, 8, 3)
    directions = np.column_stack((np.cos(np.deg2rad(azimuth)), np.sin(np.deg2rad(azimuth))))
    lengths = rng.uniform(1.25, 1.95, 3)
    vertices = np.vstack((np.zeros(2), np.cumsum(lengths[:, None] * directions, axis=0)))
    vertices += np.array((x, y)) - vertices.mean(axis=0)
    midpoints = 0.5 * (vertices[:-1] + vertices[1:])
    width = draw(0.44, 0.58)
    record("bent_slot_unseen_geometry", [
        capsule(float(cx), float(cy), float(length) / 2 + width, width, float(theta), -1.25)
        for (cx, cy), length, theta in zip(midpoints, lengths, azimuth)
    ], centerline_vertices_mm=vertices.tolist(), excluded_from_development_shapes=True)
    record("empty", [])
    return records


def forward_batch(f_ghz, truth, truth_fine, epsilon, factor, host_model, verify=False):
    fine_side = truth_fine.shape[1]
    air = Air(fine_side, f_ghz)
    apply = air.convolution()
    m_h = epsilon - 1
    empty, stats = fft_field(air, apply, np.full(fine_side**2, m_h))
    y_h = air.S @ (m_h * empty)
    shape = (len(truth_fine), len(RECEIVERS), len(SOURCES))
    clean = np.zeros(shape, dtype=complex)
    clean_same = np.zeros_like(clean)
    model = host_model(COARSE, f_ghz)
    air_same = Air(COARSE, f_ghz)
    y_h_same = air_same.S @ (m_h * model.Eh)
    max_residual = stats["max_relative_residual"]
    iterations = stats["krylov_iterations"]
    for j, q_ref in enumerate(truth_fine):
        if not np.any(q_ref):
            continue
        m = m_h + factor * q_ref.reshape(-1, 1)
        field, stats = fft_field(air, apply, m, empty)
        clean[j] = air.S @ (m * field) - y_h
        clean_same[j] = differential(model, factor * truth[j])
        max_residual = max(max_residual, stats["max_relative_residual"])
        iterations += stats["krylov_iterations"]
    diag = dict(frequency_ghz=float(f_ghz), fine_tx_rhs=len(truth_fine) * len(SOURCES),
                fine_krylov_iterations=iterations, max_fine_relative_state_residual=max_residual)
    if verify:
        fine_green = air.dense()
        direct, residual = dense_absolute(air, fine_green, factor * truth_fine[0], m_h)
        direct_bg, _ = dense_absolute(air, fine_green, np.zeros(fine_side**2), m_h)
        coarse_green = air_same.dense()
        coarse_direct, _ = dense_absolute(air_same, coarse_green, factor * truth[0], m_h)
        coarse_bg, _ = dense_absolute(air_same, coarse_green, np.zeros(COARSE**2), m_h)
        diag.update(dense_fft_relative_fine_differential=relative(clean[0], direct - direct_bg),
                    dense_host_relative_coarse_differential=relative(clean_same[0], coarse_direct - coarse_bg),
                    dense_fine_relative_state_residual=residual)
    return clean, clean_same, y_h, y_h_same, diag


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, default=RESULTS / "acquisition.npz",
                    help="Load the original dev_hole reference RMS without recalibration")
    args = ap.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    with np.load(args.reference, allow_pickle=False) as data:
        reference = float(data["reference_amplitude"])
    definitions = cases()
    truth_fine = np.array([material(case["primitives"]) for case in definitions])
    truth = truth_fine.reshape(len(definitions), COARSE, 2, COARSE, 2).mean(axis=(2, 4))
    shape = (len(definitions), len(FREQUENCIES_GHZ), len(RECEIVERS), len(SOURCES))
    rng = np.random.default_rng(NOISE_SEED)
    noise = (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)
    common = dict(frequencies_ghz=FREQUENCIES_GHZ, case_ids=[case["id"] for case in definitions],
                  splits=[case["split"] for case in definitions], noise=noise,
                  reference_amplitude=reference, sigma_2pct=NOISE_FRACTION * reference,
                  sources_m=SOURCES, receivers_m=RECEIVERS)
    for name, eps, factors, host_model in (
        ("holdout", np.full_like(FREQUENCIES_GHZ, PLATE_PERMITTIVITY),
         np.ones_like(FREQUENCIES_GHZ), host),
        ("holdout_debye", debye_epsilon(FREQUENCIES_GHZ), dielectric_factor(FREQUENCIES_GHZ), host_debye),
    ):
        wall, cpu = perf_counter(), process_time()
        clean = np.empty(shape, dtype=complex)
        same = np.empty_like(clean)
        backgrounds = np.empty(shape[1:], dtype=complex)
        backgrounds_same = np.empty_like(backgrounds)
        diagnostics = []
        for i, (f, e, factor) in enumerate(zip(FREQUENCIES_GHZ, eps, factors)):
            clean[:, i], same[:, i], backgrounds[i], backgrounds_same[i], rec = forward_batch(
                float(f), truth, truth_fine, e, factor, host_model,
                verify=(name == "holdout_debye" and f in (16., 30.)),
            )
            diagnostics.append(rec)
            print(f"{name:13s} {f:5.2f} GHz  {rec['fine_krylov_iterations']:5d} iterations  "
                  f"residual {rec['max_fine_relative_state_residual']:.2e}", flush=True)
        observations = dict(**common, clean=clean, clean_same=same,
                            background=backgrounds, background_same=backgrounds_same,
                            background_epsilon=eps, dielectric_factor=factors)
        np.savez_compressed(RESULTS / f"{name}_observations.npz", **observations)
        np.savez_compressed(RESULTS / f"{name}.npz", **observations, truth=truth, truth_fine=truth_fine)
        meta = dict(name=name, reference_case_id="dev_hole", reference_frequency_ghz=16.0,
                    reference_amplitude=reference, sigma_2pct=NOISE_FRACTION * reference,
                    loss_law="2.10+0.35/(1-i*2*pi*f*10 ps)" if name == "holdout_debye" else "none",
                    background_epsilon_real_imag=np.column_stack((eps.real, eps.imag)).tolist(),
                    dielectric_factor_real_imag=np.column_stack((factors.real, factors.imag)).tolist(),
                    noise_seed=NOISE_SEED, shape_seed=SHAPE_SEED, case_definitions=definitions,
                    grid_side_mm=1000 * PLATE_SIDE, coarse_cells_per_side=COARSE, fine_cells_per_side=FINE,
                    frequency_ghz=FREQUENCIES_GHZ.tolist(),
                    cross_grid_relative=[None if not case["primitives"] else [
                        relative(same[j, i], clean[j, i]) for i in range(len(eps))]
                        for j, case in enumerate(definitions)],
                    diagnostics=diagnostics, wall_seconds=perf_counter() - wall,
                    cpu_seconds=process_time() - cpu)
        (RESULTS / f"{name}.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
