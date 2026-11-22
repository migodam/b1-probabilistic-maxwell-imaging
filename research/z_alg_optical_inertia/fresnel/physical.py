"""Actual Fresnel readings and calibrated 2D scalar-TM disk collocation."""
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import solve
from scipy.special import hankel1, jv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE.parent / "results/fresnel"
R8 = ROOT / "research/delegated/max_struct_r1/r8_worktree"
if R8.is_dir():
    sys.path.insert(0, str(R8))
from research.research_max_struct_r8.src.fresnel_adapter import incident_fit, load_file

RAW = ROOT / "research/delegated/a2_literature/data/2005_iop_21_6_S09"
NAMES = ("FoamDielIntTM", "FoamDielExtTM", "FoamTwinDielTM")
FREQUENCIES = (4, 6, 8)
SOURCES = dict(zip(NAMES, (8, 8, 18), strict=True))
SIDE = 60
WIDTH = 0.20
DISTANCE = 1.67


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(zip(("wall_seconds", "cpu_seconds"),
                    (x - y for x, y in zip(clock(), start)), strict=True))


def adjoint(array):
    return array.conj().T


def hermitian(array):
    return (array + adjoint(array)) * .5


def unit(degrees):
    angle = np.deg2rad(degrees)
    return np.c_[np.cos(angle), np.sin(angle)]


def grid(side):
    spacing = WIDTH / side
    axis = (np.arange(side) + .5) * spacing - WIDTH / 2
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    return np.c_[xx.ravel(), yy.ravel()]


def disk_kernel(k, spacing):
    radius = spacing / np.sqrt(np.pi)
    amplitude = .5j * np.pi * k * radius * jv(1, k * radius)
    diagonal = .5j * np.pi * k * radius * hankel1(1, k * radius) - 1
    return amplitude, diagonal


def operators(frequency, side=SIDE):
    k = 2 * np.pi * frequency * 1e9 / 299792458.
    spacing = WIDTH / side
    centers = grid(side)
    receivers = DISTANCE * unit(np.arange(360))
    amplitude, diagonal = disk_kernel(k, spacing)

    # Same R8 Richmond kernel, assembled from a translation lookup to avoid
    # re-evaluating special functions for the same displacement.
    offsets = spacing * np.arange(side)
    distances = np.hypot(offsets[:, None], offsets)
    distances[0, 0] = 1.
    kernel = amplitude * hankel1(0, k * distances)
    kernel[0, 0] = diagonal
    iy, ix = np.indices((side, side)).reshape(2, -1)
    green = kernel[np.abs(iy[:, None] - iy), np.abs(ix[:, None] - ix)]
    sensor = amplitude * hankel1(0, k * np.linalg.norm(receivers[:, None] - centers, axis=2))
    return dict(k=k, G=green, S=sensor, receivers=receivers, centers=centers)


def nominal_contrast(name, centers, side):
    # In the axes used by the historical R8 nominal Twin control.  The source
    # paper's approximate Figure 5 offsets are diagnostic geometry.
    h = WIDTH / side
    offsets = (np.arange(8) + .5) / 8 - .5
    material = np.zeros(len(centers))
    for ox in offsets:
        for oy in offsets:
            sample = centers + h * np.array([ox, oy])
            foam = np.hypot(sample[:, 0], sample[:, 1]) < .04
            plastic = np.zeros(len(centers), dtype=bool)
            if name in {NAMES[0], NAMES[2]}:
                plastic |= np.hypot(sample[:, 0] + .005, sample[:, 1]) < .0155
            if name in {NAMES[1], NAMES[2]}:
                plastic |= np.hypot(sample[:, 0] + .060, sample[:, 1]) < .0155
            material += np.where(plastic, 2., np.where(foam, .45, 0.)) / 64
    return material


def nominal_forward(name, k, receiver, incident, side):
    centers = grid(side)
    material = nominal_contrast(name, centers, side)
    support = material > 0
    centers, material = centers[support], material[support]
    amplitude, diagonal = disk_kernel(k, WIDTH / side)
    distances = np.linalg.norm(centers[:, None] - centers, axis=2)
    np.fill_diagonal(distances, 1.)
    system = -material[:, None] * amplitude * hankel1(0, k * distances)
    np.fill_diagonal(system, 1. - material * diagonal)
    rhs = material[:, None] * incident(centers)
    current = solve(system, rhs, assume_a="gen", check_finite=False)
    s = amplitude * hankel1(0, k * np.linalg.norm(receiver[:, None] - centers, axis=2))
    record = dict(side=side, active_cells=int(support.sum()), direct_solves=1,
                  source_rhs=rhs.shape[1],
                  relative_state_residual=float(np.linalg.norm(system @ current - rhs) /
                                                np.linalg.norm(rhs)))
    return s @ current, record


def load_case(name, frequency, geometry, *, model_audit):
    start = clock()
    path = RAW / f"{name}.exp"
    raw, measurement, empty = load_file(path, frequency)
    tx = np.unique(raw[:, 0]).astype(int) - 1
    ns = SOURCES[name]
    if not np.array_equal(tx, np.arange(ns)):
        raise ValueError(f"{name} unexpected raw transmit indices: {tx}")
    measurement, empty = measurement[:, tx], empty[:, tx]
    mask = np.isfinite(measurement)
    source_degrees = np.arange(ns) * (360 // ns)
    source_positions = DISTANCE * unit(source_degrees)
    relative_angle = (np.arange(360)[:, None] - source_degrees) % 360
    if not np.array_equal(mask, (60 <= relative_angle) & (relative_angle <= 300)):
        raise ValueError("Raw gaps disagree with the primary acquisition geometry.")
    k, centers, receivers = (geometry[key] for key in ("k", "centers", "receivers"))
    incident, fit_choices, predict_incident = incident_fit(empty, source_positions, receivers,
                                                          centers, k)
    # Antenna-coincident evaluations are absent for that source by acquisition design.
    with np.errstate(invalid="ignore", divide="ignore"):
        empty_predicted = predict_incident(receivers)
    empty_val = np.zeros_like(mask)
    cyl_gain = np.empty(ns, complex)
    opposite = (source_degrees + 180) % 360
    for source, row in enumerate(fit_choices):
        empty_val[row["validation_channels"], source] = True
        cyl_gain[source] = empty[opposite[source], source] / (
            .25j * hankel1(0, k * (2 * DISTANCE)))
    norm = float(np.linalg.norm(measurement[mask]))
    if not np.isfinite(empty_predicted[empty_val]).all():
        raise ArithmeticError("Nonfinite calibrated empty-field holdout.")
    empty_defect = np.linalg.norm((empty_predicted - empty)[empty_val])
    empty_val_incident = np.linalg.norm(empty[empty_val])
    empty_val_scatter = np.linalg.norm(measurement[empty_val])
    # Convert the empty test RMS to an all-observed-channel sensitivity scale.
    # This is empirical bookkeeping and is not a noise upper bound.
    empty_rms = float(empty_defect / np.sqrt(empty_val.sum()))
    opposite_gain_incident = .25j * hankel1(0, k * np.linalg.norm(
        centers[:, None] - source_positions, axis=2)) * cyl_gain
    cyl_predict = .25j * hankel1(0, k * np.linalg.norm(
        receivers[:, None] - source_positions, axis=2)) * cyl_gain
    record = dict(name=name, frequency_GHz=frequency, id=f"fresnel_2005/{name}/raw/{frequency}GHz",
                  raw_file=str(path.relative_to(ROOT)), measured_norm=norm,
                  field_unit="complex raw receiver units, normalized by all measured scattered Frobenius norm",
                  source_degrees=source_degrees.tolist(), source_count=ns,
                  measured_channels=int(mask.sum()), missing_channels=int(mask.size-mask.sum()),
                  common_full_array_receiver_count=int(mask.all(axis=1).sum()),
                  empty_calibration=dict(sources=fit_choices, gain_cyl_amplitude_min=float(abs(cyl_gain).min()),
                      gain_cyl_amplitude_max=float(abs(cyl_gain).max()),
                      gain_cyl_phase_degrees=np.angle(cyl_gain, deg=True).tolist(),
                      held_channels=int(empty_val.sum()), held_empty_relative=float(empty_defect / empty_val_incident),
                      held_empty_error_over_held_scatter=float(empty_defect / empty_val_scatter),
                      held_empty_error_over_all_scatter=float(empty_defect / norm),
                      rms_transferred_over_all_scatter=float(empty_rms*np.sqrt(mask.sum())/norm),
                      opposite_cylinder_held_empty_relative=float(np.linalg.norm(
                          (cyl_predict-empty)[empty_val]) / empty_val_incident),
                      incident_center_multipole_vs_cylinder=float(np.linalg.norm(
                          incident - opposite_gain_incident) / np.linalg.norm(incident))),
                  preparation=elapsed(start))
    if model_audit:
        audit_start = clock()
        hold_tx = mask & (tx % 4 == 0)[None, :]
        hold_rx = mask & ~hold_tx & ((relative_angle - 60) % 5 == 0)
        train = mask & ~hold_tx & ~hold_rx
        prediction, solve_record = nominal_forward(name, k, receivers, predict_incident, SIDE)
        fine_prediction, fine_record = nominal_forward(name, k, receivers, predict_incident, 120)
        gain = np.vdot(prediction[train], measurement[train]) / np.vdot(prediction[train], prediction[train])
        residuals = {}
        for label, channels in (("all", mask), ("gain_fit", train), ("held_sources", hold_tx),
                                ("held_receivers_on_fit_sources", hold_rx)):
            ref = measurement[channels]
            pred, fine = prediction[channels], fine_prediction[channels]
            denominator = np.linalg.norm(ref)
            residuals[label] = dict(channels=int(channels.sum()),
                                    nominal60=float(np.linalg.norm(pred - ref) / denominator),
                                    nominal120=float(np.linalg.norm(fine - ref) / denominator),
                                    nominal60_with_train_gain=float(np.linalg.norm(gain * pred - ref) / denominator),
                                    discretization_shift=float(np.linalg.norm(fine - pred) / denominator))
        record["model_audit"] = dict(relative_residuals=residuals, direct_solves=[solve_record, fine_record],
             diagnostic_scatter_gain_abs=float(abs(gain)), diagnostic_scatter_gain_phase_deg=float(np.angle(gain, deg=True)),
             diagnostic_scatter_gain_real_imag=[float(gain.real), float(gain.imag)],
             model_error_hard_upper_bound=None, **elapsed(audit_start))
        np.savez_compressed(RESULTS/f"{name}_{frequency}GHz_ports.npz",
                            raw_scattered=measurement, raw_empty=empty, mask=mask,
                            source_degrees=source_degrees, nominal60=prediction, nominal120=fine_prediction,
                            empty_held_mask=empty_val, nominal_held_tx=hold_tx, nominal_held_rx=hold_rx,
                            diagnostic_scatter_gain=gain)
    return dict(B=incident/norm, Y=measurement/norm, mask=mask, record=record, norm=norm)
