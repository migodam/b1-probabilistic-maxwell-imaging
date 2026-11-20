"""Compare equal-readout Fresnel TM protocols by locating one plastic cylinder."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lstsq
from scipy.optimize import minimize
from scipy.special import h1vp, hankel1, jv, jvp

from ndt import C0, ring


TARGETS = ("FoamDielIntTM", "FoamDielExtTM")
SOURCE_IDS = np.arange(8)
SOURCES = ring(8, 1.67)
RECEIVERS = ring(360, 1.67)
RECEIVER_SUBSET = RECEIVERS[::5]
TRAIN_ANGLES = np.arange(60, 300, 10)
TEST_ANGLES = TRAIN_ANGLES + 5
TRAIN_RX = ((45 * SOURCE_IDS[:, None] + TRAIN_ANGLES) % 360) // 5
TEST_RX = ((45 * SOURCE_IDS[:, None] + TEST_ANGLES) % 360) // 5
FREQUENCIES = np.arange(2, 11)
FEW_FREQUENCIES = (2, 6, 10)
MORE_FREQUENCIES = (2, 3, 4, 6, 8, 10)
HOLD_FREQUENCIES = (5, 7, 9)
ARMS = {
    "6f4tx_0deg": (MORE_FREQUENCIES, SOURCE_IDS[::2]),
    "6f4tx_45deg": (MORE_FREQUENCIES, SOURCE_IDS[1::2]),
    "3f8tx": (FEW_FREQUENCIES, SOURCE_IDS),
}
FOAM_R = 0.040
PLASTIC_R = 0.0155
FOAM_EPS = 1.45
PLASTIC_EPS = 3.0
SEARCH_MM = {"nested": (0., 24.4), "exterior": (55.7, 86.)}
RESULT_ROOT = Path(__file__).resolve().parents[1]


def load_measurements(path):
    rows = np.loadtxt(path)
    if rows.shape != (9 * 8 * 241, 7) or not np.isfinite(rows).all():
        raise ValueError(f"Invalid Fresnel TM input shape or contents: {path}")
    tx = rows[:, 0].astype(int) - 1
    rx = rows[:, 1].astype(int) - 1
    f = rows[:, 2].astype(int) - FREQUENCIES[0]
    indexes = (f, rx, tx)
    if len(set(zip(*indexes))) != len(rows):
        raise ValueError(f"Duplicate measurement: {path}")
    total = np.full((9, 360, 8), np.nan + 1j * np.nan)
    incident = total.copy()
    # The file is conjugated to e^(-i omega t) / outgoing H_0^(1).
    total[indexes] = rows[:, 3] - 1j * rows[:, 4]
    incident[indexes] = rows[:, 5] - 1j * rows[:, 6]
    rel_rx = (rx - 45 * tx) % 360
    if rel_rx.min() != 60 or rel_rx.max() != 300:
        raise ValueError(f"Unexpected receiver coordinates: {path}")
    return total - incident, incident


def wave(points, center, modes, k, kind, normals=None):
    v = points - center
    r = np.linalg.norm(v, axis=1)
    theta = np.arctan2(v[:, 1], v[:, 0])
    phase = np.exp(1j * theta[:, None] * modes)
    fn, derivative = (hankel1, h1vp) if kind == "H" else (jv, jvp)
    w = fn(modes[None, :], k * r[:, None])
    if normals is None:
        return phase * w
    u = v / r[:, None]
    dr = np.sum(u * normals, axis=1)
    tang = np.column_stack((-u[:, 1], u[:, 0]))
    dt = np.sum(tang * normals, axis=1) / r
    dw = k * derivative(modes[None, :], k * r[:, None]) * dr[:, None]
    dw = dw + 1j * w * modes * dt[:, None]
    return phase * w, phase * dw


def incident_basis(points, source, k, normals=None):
    angle = np.arctan2(-source[1], -source[0])
    modes = np.arange(-2, 3)
    # H_{-m}=(-1)^m H_m. This rephases the basis toward the antenna axis.
    factors = 0.25j * np.exp(-1j * modes * angle) * np.where(modes < 0, (-1.) ** modes, 1.)
    values = wave(points, source, modes, k, "H", normals)
    if normals is None:
        return values * factors
    return values[0] * factors, values[1] * factors


def calibrate(incident, k):
    """Fit only target-absent fields, as in the R8 incident continuation."""
    coefs, validation, opposite_comparator = [], [], []
    for s, source in enumerate(SOURCES):
        opposite = (180 + 45 * s) % 360
        near = (opposite + np.arange(-30, 31)) % 360
        keep = np.arange(len(near)) % 3 != 0
        basis = incident_basis(RECEIVERS[near], source, k)
        scale = np.linalg.norm(basis[keep], axis=0)
        x, _, _, _ = lstsq(basis[keep] / scale, incident[near[keep], s], cond=1e-4)
        coef = x / scale
        coefs.append(coef)
        reference = incident[near[~keep], s]
        cyl = basis[:, 2] / basis[30, 2] * incident[opposite, s]
        validation.append(rel(basis[~keep] @ coef, reference))
        opposite_comparator.append(rel(cyl[~keep], reference))
    return np.array(coefs), validation, opposite_comparator


def incident_at(points, normals, k, coefs):
    items = [incident_basis(points, source, k, normals) for source in SOURCES]
    electric = np.stack([b @ c for (b, _), c in zip(items, coefs)], axis=1)
    normal_derivative = np.stack([d @ c for (_, d), c in zip(items, coefs)], axis=1)
    return np.vstack((electric, normal_derivative / k))


def rel(value, reference):
    return float(np.linalg.norm(value - reference) / np.linalg.norm(reference))


class Cylinders:
    """Solve the coupled scalar TM circular interfaces in a multipole basis."""

    def __init__(self, f, coefs, extra_order=0):
        self.k = 2 * np.pi * f * 1e9 / C0
        self.kh = self.k * np.sqrt(FOAM_EPS)
        self.kp = self.k * np.sqrt(PLASTIC_EPS)
        self.coefs = coefs
        order_h = int(np.ceil(self.kh * FOAM_R)) + 7 + extra_order
        order_p = int(np.ceil(self.kp * PLASTIC_R)) + 7 + extra_order
        self.h_modes = np.arange(-order_h, order_h + 1)
        self.p_modes = np.arange(-order_p, order_p + 1)
        self.nh, self.np = self.h_modes.size, self.p_modes.size
        self.normal_h = ring(2 * self.nh, 1.0)
        self.normal_p = ring(2 * self.np, 1.0)
        self.pts_h = FOAM_R * self.normal_h
        self.pts_p = PLASTIC_R * self.normal_p
        zero = np.zeros(2)
        h = self.surface(self.pts_h, self.normal_h, zero, self.h_modes, self.k, "H")
        j = self.surface(self.pts_h, self.normal_h, zero, self.h_modes, self.kh, "J")
        j_p = self.surface(self.pts_p, self.normal_p, zero, self.p_modes, self.kp, "J")
        h_p_air = self.surface(self.pts_p, self.normal_p, zero, self.p_modes, self.k, "H")
        h_p_foam = self.surface(self.pts_p, self.normal_p, zero, self.p_modes, self.kh, "H")
        self.fixed = h, j, j_p, h_p_air, h_p_foam
        self.h_inc = incident_at(self.pts_h, self.normal_h, self.k, self.coefs)
        self.h_receive = wave(RECEIVER_SUBSET, zero, self.h_modes, self.k, "H")

    def surface(self, pts, normals, center, modes, k, kind):
        b, db = wave(pts, center, modes, k, kind, normals)
        return np.vstack((b, db / self.k))

    def response(self, center_mm, topology, measure_residual=False):
        p_center = np.asarray(center_mm) * 1e-3
        p_surface = p_center + self.pts_p
        zero = np.zeros(2)
        h, j, j_p, h_p_air, h_p_foam = self.fixed
        z_hp = np.zeros((h.shape[0], self.np), dtype=complex)
        z_ph = np.zeros((j_p.shape[0], self.nh), dtype=complex)
        if topology == "nested":
            inner_out = self.surface(self.pts_h, self.normal_h, p_center, self.p_modes, self.kh, "H")
            host_to_inner = self.surface(p_surface, self.normal_p, zero, self.h_modes, self.kh, "J")
            a = np.block([[h, -j, -inner_out, z_hp],
                          [z_ph, host_to_inner, h_p_foam, -j_p]])
            rhs = -np.vstack((self.h_inc, np.zeros((j_p.shape[0], 8), dtype=complex)))
        else:
            plastic_to_host = self.surface(self.pts_h, self.normal_h, p_center, self.p_modes, self.k, "H")
            host_to_plastic = self.surface(p_surface, self.normal_p, zero, self.h_modes, self.k, "H")
            a = np.block([[h, -j, plastic_to_host, z_hp],
                          [host_to_plastic, z_ph, h_p_air, -j_p]])
            rhs = -np.vstack((self.h_inc, incident_at(p_surface, self.normal_p, self.k, self.coefs)))
        scales = np.linalg.norm(a, axis=0)
        state, _, rank, _ = lstsq(a / scales, rhs, cond=1e-11, lapack_driver="gelsy", check_finite=False)
        state /= scales[:, None]
        result = self.h_receive @ state[:self.nh]
        if topology == "exterior":
            p_receive = wave(RECEIVER_SUBSET, p_center, self.p_modes, self.k, "H")
            result += p_receive @ state[2 * self.nh:2 * self.nh + self.np]
        if measure_residual:
            return result, rel(a @ state, rhs), int(rank)
        return result


def center(polar):
    r, angle = polar
    return r * np.array((np.cos(angle), np.sin(angle)))


def candidates():
    for topology, radii, n_angle in [("nested", [0, 8, 16, 23], 8),
                                      ("exterior", [55.7, 64, 75, 86], 24)]:
        for r in radii:
            for theta in (np.array([0.]) if r == 0 else np.arange(n_angle) * 2 * np.pi / n_angle):
                yield topology, (r, theta)


def channels(field, rx, tx=SOURCE_IDS):
    return field[rx[tx], tx[:, None]]


def fit_locations(models, measured):
    """Coarse identical pose dictionary and bounded same-objective refinement."""
    def training_error(f, pred, tx):
        diff = channels(pred - measured[f][::5], TRAIN_RX, tx)
        return float(np.vdot(diff, diff).real)

    totals = {name: sum(float(np.linalg.norm(channels(measured[f][::5], TRAIN_RX, tx))**2)
                        for f in ff) for name, (ff, tx) in ARMS.items()}
    coarse_results = {name: {"nested": [], "exterior": []} for name in ARMS}
    for topology, p in candidates():
        predictions = {f: models[f].response(center(p), topology) for f in MORE_FREQUENCIES}
        for name, (ff, tx) in ARMS.items():
            score = sum(training_error(f, predictions[f], tx) for f in ff) / totals[name]
            coarse_results[name][topology].append((score, p))

    result = {}
    for name, (ff, tx) in ARMS.items():
        branches = {}
        for topology, evidence in coarse_results[name].items():
            def objective(p):
                return sum(training_error(f, models[f].response(center(p), topology), tx)
                           for f in ff) / totals[name]

            coarse_score, p0 = min(evidence, key=lambda item: item[0])
            start = np.array(p0)
            if start[0] == 0:
                # A centered start has no angular derivative; use a finite simplex.
                simplex = np.array([[0, 0], [2., 0], [2., np.pi / 2]])
            else:
                dr = -1.5 if start[0] == SEARCH_MM[topology][1] else 1.5
                simplex = np.array([start, start + [dr, 0], start + [0, np.pi / 30]])
            fit = minimize(objective, start, method="Nelder-Mead", bounds=[SEARCH_MM[topology], (None, None)],
                           options={"maxfev": 90, "xatol": 0.025, "fatol": 1e-7, "initial_simplex": simplex})
            branches[topology] = dict(center_mm=center(fit.x).tolist(),
                                      training_relative=float(np.sqrt(fit.fun)),
                                      dictionary_min_relative=float(np.sqrt(coarse_score)),
                                      refinement_evaluations=int(fit.nfev),
                                      refinement_converged=bool(fit.success))
        topology = min(branches, key=lambda item: branches[item]["training_relative"])
        result[name] = dict(topology=topology, center_mm=branches[topology]["center_mm"],
                            training_relative=branches[topology]["training_relative"],
                            branches=branches, training_frequencies_ghz=list(ff),
                            training_source_degrees=(45 * tx).tolist(),
                            training_complex_readouts=len(ff) * len(tx) * len(TRAIN_ANGLES))
        print(f"  {name}: {topology}, center {result[name]['center_mm']}, "
              f"train relative {result[name]['training_relative']:.4f}", flush=True)
    return result


def prediction_error(pred, data, mask):
    p, y = pred[mask], data[mask]
    error = np.angle(p * y.conj())
    return dict(complex_relative=rel(p, y), magnitude_relative=rel(abs(p), abs(y)),
                target_power_weighted_circular_phase_rms_degrees=float(np.degrees(
                    np.sqrt(np.average(error**2, weights=abs(y)**2)))),
                complex_readouts=int(mask.sum()))


def score(models, fits, measured):
    target_cube = np.stack([measured[f][::5] for f in FREQUENCIES])
    available = np.isfinite(target_cube)
    rx_hold = np.zeros_like(available)
    rx_hold[:, TEST_RX, SOURCE_IDS[:, None]] = True
    withheld_f = np.isin(FREQUENCIES, HOLD_FREQUENCIES)
    common = np.isin(FREQUENCIES, FEW_FREQUENCIES)
    predictions = []
    for name, fit in fits.items():
        fields = [models[f].response(fit["center_mm"], fit["topology"], measure_residual=True)
                  for f in FREQUENCIES]
        cube = np.stack([field[0] for field in fields])
        predictions.append(cube)
        fit["max_relative_interface_residual"] = max(field[1] for field in fields)
        fit["held_receivers_common_training_f"] = prediction_error(cube, target_cube, rx_hold & common[:, None, None])
        fit["held_frequencies"] = prediction_error(cube, target_cube, available & withheld_f[:, None, None])
        fit["held_frequencies_per_f"] = {
            str(f): prediction_error(cube[i], target_cube[i], available[i])
            for i, f in enumerate(FREQUENCIES) if withheld_f[i]
        }
    return np.stack(predictions)


def noise_proxy(all_data, all_inc):
    # Difference of the two empty sweeps includes drift. As a rough scale for
    # total-minus-empty, the ratio is not a calibrated device noise estimate.
    a, b = (all_inc[name] for name in TARGETS)
    rows = {}
    for name in TARGETS:
        rows[name] = {}
        for f in FREQUENCIES:
            uncertainty = channels((a[f] - b[f])[::5], TRAIN_RX)
            y = channels(all_data[name][f][::5], TRAIN_RX)
            rows[name][int(f)] = dict(aperture_scattered_rms=float(np.sqrt(np.mean(abs(y)**2))),
                                      empty_sweep_difference_rms=float(np.sqrt(np.mean(abs(uncertainty)**2))),
                                      empty_sweep_proxy_snr_db=float(20 * np.log10(np.linalg.norm(y) / np.linalg.norm(uncertainty))))
    return rows


def make_figure(result, path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    colors = {"6f4tx_0deg": "#007f88", "6f4tx_45deg": "#3b9ca7", "3f8tx": "#bb5c2c"}
    labels = {"6f4tx_0deg": "6 frequencies × 4 views (0°)",
              "6f4tx_45deg": "6 frequencies × 4 views (45°)",
              "3f8tx": "3 frequencies × 8 views"}
    markers = {"6f4tx_0deg": "o", "6f4tx_45deg": "^", "3f8tx": "s"}
    titles = {"FoamDielIntTM": "Embedded inclusion", "FoamDielExtTM": "Adjacent inclusion"}
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 5.2), width_ratios=(0.82, 1.18))
    for row, (name, case) in enumerate(result["targets"].items()):
        geometry, errors = axes[row]
        g_inset = geometry.inset_axes([0.56, 0.24, 0.42, 0.62])
        geo_ref = next(iter(case["arms"].values()))["schematic_center_mm"]
        offsets = np.array([np.array(a["center_mm"]) - geo_ref for a in case["arms"].values()])
        inset_extent = max(3, np.ceil(1.25 * abs(offsets).max()))
        geometry.add_patch(Circle((0, 0), 40, facecolor="#ebeff4", edgecolor="#8796a7", lw=1.1))
        geometry.add_patch(Circle(geo_ref, 15.5, fill=False, edgecolor="#414957", lw=1.1, ls="--"))
        geometry.scatter(*geo_ref, marker="+", color="#414957", s=85, zorder=5)
        g_inset.scatter(0, 0, marker="+", color="#414957", s=55, label="Nominal center", zorder=4)
        g_inset.axhline(0, lw=0.5, color="0.77")
        g_inset.axvline(0, lw=0.5, color="0.77")
        for offset, (key, arm) in zip((-0.14, 0, 0.14), case["arms"].items()):
            xy = arm["center_mm"]
            delta = np.array(xy) - geo_ref
            g_inset.scatter(*delta, color=colors[key], s=22, marker=markers[key],
                            linewidth=0.8, alpha=0.83, zorder=5)
            metrics = [arm["held_frequencies_per_f"][str(f)]["complex_relative"] for f in HOLD_FREQUENCIES]
            metrics.append(arm["held_receivers_common_training_f"]["complex_relative"])
            errors.scatter(np.arange(len(metrics)) + offset, 100 * np.array(metrics), color=colors[key],
                           marker=markers[key], s=24, label=labels[key])
        geometry.set(xlim=(-83, 87), ylim=(-45, 45), aspect="equal", xlabel="x (mm)", ylabel="y (mm)")
        geometry.set_title(titles[name], loc="left", fontsize=10, fontweight="bold")
        g_inset.set(xlim=(-inset_extent, inset_extent), ylim=(-inset_extent, inset_extent), aspect="equal")
        g_inset.set_xticks([-inset_extent, 0, inset_extent])
        g_inset.set_yticks([-inset_extent, 0, inset_extent])
        g_inset.tick_params(labelsize=6, length=2, pad=1)
        g_inset.set_title("Center offset (mm)", fontsize=7, pad=2)
        errors.set_xticks(range(4), labels=["5 GHz", "7 GHz", "9 GHz", "Receiver\nangles"])
        errors.set(ylabel="Relative complex error (%)", ylim=(0, 35), xlim=(-0.4, 3.4))
        errors.axvline(2.5, color="0.8", lw=0.8, ls="--")
        errors.set_title("Held-out frequencies and receiver angles", loc="left", fontsize=10, fontweight="bold")
        errors.grid(axis="y", color="0.89", lw=0.7)
    handles, legend_labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.08, bottom=0.14, right=0.98, top=0.94, hspace=0.58, wspace=0.31)
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Directory containing the original Institut Fresnel .exp files")
    args = parser.parse_args()
    wall, cpu = perf_counter(), process_time()
    data, inc = {}, {}
    for name in TARGETS:
        scatter, incident = load_measurements(args.data_dir / (name + ".exp"))
        data[name] = dict(zip(FREQUENCIES, scatter))
        inc[name] = dict(zip(FREQUENCIES, incident))

    records, pred = {}, []
    for name in TARGETS:
        print(name, flush=True)
        calibration = {f: calibrate(inc[name][f], 2 * np.pi * f * 1e9 / C0) for f in FREQUENCIES}
        models = {f: Cylinders(f, calibration[f][0]) for f in FREQUENCIES}
        fits = fit_locations(models, data[name])
        pred.append(score(models, fits, data[name]))
        # Check the full-wave truncation at the endpoints of the test band.
        for fit in fits.values():
            fit["modal_refinement_5_9ghz"] = {}
            for f in (5, 9):
                p0 = models[f].response(fit["center_mm"], fit["topology"])
                p1, residual, rank = Cylinders(f, calibration[f][0], 7).response(
                    fit["center_mm"], fit["topology"], measure_residual=True)
                fit["modal_refinement_5_9ghz"][str(f)] = dict(
                    relative_field_change=rel(p0, p1), relative_interface_residual=residual, rank=rank)
        records[name] = dict(arms=fits, incident_validation={
            str(f): dict(order2_relative_errors=calibration[f][1],
                         opposite_monopole_relative_errors=calibration[f][2]) for f in FREQUENCIES})

    # Published engineering schematic centers are used only after all fits.
    nominal = {"FoamDielIntTM": [-5., 0.], "FoamDielExtTM": [-55.5, 0.]}
    for name, record in records.items():
        for fit in record["arms"].values():
            fit["schematic_center_mm"] = nominal[name]
            fit["center_error_to_schematic_mm"] = float(np.linalg.norm(np.array(fit["center_mm"]) - nominal[name]))
    result = dict(targets=records, retained_real_measurement_snr_proxy=noise_proxy(data, inc),
                  train_receiver_relative_deg=TRAIN_ANGLES.tolist(), held_receiver_relative_deg=TEST_ANGLES.tolist(),
                  held_frequencies_ghz=list(HOLD_FREQUENCIES), source_receiver_radius_m=1.67,
                  prior=dict(foam_radius_mm=FOAM_R * 1e3, plastic_radius_mm=PLASTIC_R * 1e3,
                             foam_epsilon=FOAM_EPS, plastic_epsilon=PLASTIC_EPS, search_radius_mm=SEARCH_MM),
                  elapsed_seconds=perf_counter() - wall, cpu_seconds=process_time() - cpu)
    out = RESULT_ROOT / "results"
    (out / "fresnel_measured.json").write_text(json.dumps(result, indent=2) + "\n")
    np.savez_compressed(out / "fresnel_predictions.npz", targets=TARGETS, arms=list(ARMS),
                        frequencies_ghz=FREQUENCIES, receiver_global_deg=np.arange(0, 360, 5),
                        source_global_deg=45 * SOURCE_IDS, prediction=np.stack(pred))
    make_figure(result, RESULT_ROOT / "figures" / "fresnel_equal_readout.png")


if __name__ == "__main__":
    main()
