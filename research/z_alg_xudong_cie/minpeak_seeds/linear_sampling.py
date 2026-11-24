"""Build two data-only support seeds from the frozen full-aperture port matrix."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve, svd
from scipy.ndimage import gaussian_filter
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1] / "z_alg_optical_inertia"
PUBLIC = STUDY / "synthetic/public"
GEOMETRY = "PW360_10x20"
CASE = "optfresh_01__pw360_10x20__n02"
SIDE = 32
NATIVE_SIDE = 96
START_SIDE = 24
CONTRASTS = (.75, 1.5, 2.4, 3.2)
ROI_FRACTION = .28


def time_since(wall, cpu):
    return dict(wall_seconds=perf_counter() - wall, cpu_seconds=process_time() - cpu)


def unit_columns(matrix):
    norms = np.linalg.norm(matrix, axis=0)
    if np.any(norms == 0):
        raise ArithmeticError("The sampling vectors must be nonzero.")
    return matrix / norms


def construct_indicators(y, sensor, incident, radius):
    """Interrogate the two singular spaces with their matching cell responses."""
    u, singular, vh = svd(y, full_matrices=True, check_finite=False)
    receiver = u.conj().T @ unit_columns(sensor)
    source = vh @ unit_columns(incident.conj().T)
    left_energy = abs(receiver)**2
    right_energy = abs(source)**2
    beta = radius / np.sqrt(y.shape[1])
    observed = np.r_[singular, np.zeros(y.shape[0] - len(singular))]
    left_picard = 1 / np.sqrt((1 / (observed**2 + beta**2)) @ left_energy)
    right_picard = 1 / np.sqrt((1 / (singular**2 + beta**2)) @ right_energy)

    truncation = max(4 * beta, .08 * singular[0])
    retained = max(1, int(np.count_nonzero(singular > truncation)))
    left_range = np.sum(left_energy[:retained], axis=0)
    right_range = np.sum(right_energy[:retained], axis=0)
    return dict(smooth_picard=np.sqrt(left_picard * right_picard),
                truncated_range=np.sqrt(left_range * right_range)), dict(
                    singular_values=singular.tolist(), beta=float(beta),
                    truncation_threshold=float(truncation), retained_rank=retained)


def broad_roi(indicator):
    smoothed = gaussian_filter(indicator.reshape(SIDE, SIDE), sigma=1.0, mode="nearest")
    count = int(np.rint(ROI_FRACTION * smoothed.size))
    selected = np.argsort(smoothed.ravel(), kind="stable")[-count:]
    mask = np.zeros(smoothed.size, dtype=bool)
    mask[selected] = True
    return mask, smoothed


def coarse_trial(green, sensor, incident, y, active, contrast):
    g = green[np.ix_(active, active)]
    state = np.eye(len(active), dtype=complex) - contrast * g
    rhs = contrast * incident[active]
    current = lu_solve(lu_factor(state, check_finite=False), rhs, check_finite=False)
    state_error = np.linalg.norm(state @ current - rhs) / np.linalg.norm(rhs)
    if not state_error < 1e-9:
        raise ArithmeticError("The proposed coarse physical state did not solve its VIE.")
    output = sensor[:, active] @ current
    return dict(contrast=contrast, residual_Frobenius=float(np.linalg.norm(output - y)),
                relative_state_residual=float(state_error))


def run(output):
    wall, cpu = perf_counter(), process_time()
    geometry = json.loads((PUBLIC / "manifest.json").read_text())
    with np.load(PUBLIC / geometry["green_file"], allow_pickle=False) as source:
        green = source["G"]
    with np.load(PUBLIC / geometry["geometries"][GEOMETRY]["file"], allow_pickle=False) as source:
        sensor, incident = source["S"], source["B"]
    with np.load(PUBLIC / f"{CASE}.npz", allow_pickle=False) as source:
        y = source["Y"]
    feasibility = json.loads((STUDY / "results/theory/fine_feasibility.json").read_text())
    radius = next(row["inflated_radius"] for row in feasibility["observations"] if row["id"] == CASE)
    indicators, spectrum = construct_indicators(y, sensor, incident, radius)
    records = {}
    for kind, indicator in indicators.items():
        stage_wall, stage_cpu = perf_counter(), process_time()
        mask, smoothed = broad_roi(indicator)
        active = np.flatnonzero(mask)
        trials = [coarse_trial(green, sensor, incident, y, active, contrast) for contrast in CONTRASTS]
        chosen = min(trials, key=lambda trial: trial["residual_Frobenius"])
        coarse = chosen["contrast"] * mask.reshape(SIDE, SIDE)
        factor = NATIVE_SIDE // SIDE
        native = np.repeat(np.repeat(coarse, factor, axis=0), factor, axis=1)
        width = NATIVE_SIDE // START_SIDE
        grid24 = native.reshape(START_SIDE, width, START_SIDE, width).mean(axis=(1, 3))
        target = output / f"{CASE}__{kind}.npz"
        np.savez_compressed(target, m_seed=grid24.ravel(), m_native=native.ravel(),
                            m_coarse=coarse.ravel(), support_coarse=mask, indicator=smoothed)
        records[kind] = dict(
            path=target.name, n_active_coarse=int(len(active)), chosen=chosen, trials=trials,
            binary64_coarse_LU=len(trials), original_forward_source_rhs=y.shape[1] * len(trials),
            **time_since(stage_wall, stage_cpu))
        print(kind, chosen["contrast"], chosen["residual_Frobenius"], target, flush=True)
    receipt = dict(
        observation_id=CASE, geometry=GEOMETRY, data_boundary="public S, B, G, Y and archived radius",
        coarse_side=SIDE, native_side=NATIVE_SIDE, m_seed_side=START_SIDE, original_radius=radius,
        zero_material_residual_Frobenius=float(np.linalg.norm(y)),
        scope="initial material proposals screened on coarse full VIE; no native feasibility claim",
        roi_smoothing_sigma_coarse_cells=1.0, retained_fraction_before_threshold_rounding=ROI_FRACTION,
        source_and_receiver_sampling=spectrum, variants=records, **time_since(wall, cpu),
    )
    (output / "linear_sampling.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(1):
        run(args.output)
