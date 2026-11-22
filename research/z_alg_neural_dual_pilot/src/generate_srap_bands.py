"""Offline-only rescaling of the original SRAP22 noise directions to 1% and 2%."""

import gc
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits


PILOT = Path(__file__).resolve().parents[1]
PAPER = PILOT.parent / "z_alg_optical_inertia"
SRAP = PILOT.parent / "srap_cie_array_study_20261121"
ARCHIVE = PAPER / "validation" / "srap"
DATA = PILOT / "data"
NOISE_PCT = (1, 2)

# The original generator and verifier use local imports with common names.
# This command runs separately from the proposal and learning processes.
sys.dont_write_bytecode = True
for path in (PILOT.parents[1], SRAP, PAPER / "synthetic", ARCHIVE):
    sys.path.insert(0, str(path))

import data2d
from nonvacuity import ExactFourierGeometry, StableActiveModel
from offline_feasible import ActiveSRAP
from online import ARRAY, FINE, SCOPE, ball, frobenius_squared, high, load_geometry, low


class ExactActive(ActiveSRAP, StableActiveModel):
    """Original SRAP material support and outward correction of native port fields."""


def now():
    return perf_counter(), process_time()


def timing(t):
    return dict(zip(("wall_seconds", "cpu_seconds"), (a - b for a, b in zip(now(), t))))


def norm(matrix):
    return frobenius_squared(matrix).sqrt()


def directional_increment(index, archived_delta, shape, array_index):
    seed = np.random.SeedSequence([data2d.NOISE_SEED, index, array_index, 0])
    rng = np.random.default_rng(seed)
    direction = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    return direction * (archived_delta / np.linalg.norm(direction))


def main():
    if not (PILOT / "results" / "holdout_open.json").is_file():
        raise SystemExit("The frozen holdout has not been opened.")
    ctx.prec, ctx.threads = 112, 1
    t = now()
    generator = data2d.SyntheticOperators(96, data2d.CONFIRM_ARRAYS)
    spectrum, sensor, incident = load_geometry()
    matches = dict(
        spectrum=np.array_equal(spectrum, generator.green.spectrum),
        sensor=np.array_equal(sensor, generator.fine_s[generator.rx_ids[ARRAY]]),
        incident=np.array_equal(incident, generator.fine_e[:, generator.tx_ids[ARRAY]]),
    )
    if not all(matches.values()):
        raise ArithmeticError("The replayed generator does not match the native frozen arrays.")

    geometry = ExactFourierGeometry(FINE)
    geometry.S, geometry.B = sensor, incident
    labels = np.random.default_rng(data2d.FAMILY_SEED).permutation(np.repeat(np.arange(6), 4))
    expected_ids = [f"s2d_{index:03d}__pw360_10x20" for index in SCOPE]
    public = ARCHIVE / "public_observations.npz"
    with np.load(public, allow_pickle=False) as saved:
        if saved["ids"].tolist() != expected_ids:
            raise ValueError("The original SRAP22 order or scope differs from its frozen list.")
        archived_y, archived_delta = saved["Y"], saved["original_delta"]
    if archived_y.shape != (22, 20, 10) or archived_delta.shape != (22,) or np.any(archived_delta <= 0):
        raise ValueError("Invalid original SRAP public observation pack.")
    setup = timing(t)
    obs, radii, receipt = [], [], []
    array_index = data2d.CONFIRM_ARRAYS.index(ARRAY)
    guard = arb(1) / 10**6

    with threadpool_limits(limits=1):
        for n, index in enumerate(SCOPE):
            item_t = now()
            object_id = f"s2d_{index:03d}"
            y10 = archived_y[n]
            increment = directional_increment(index, archived_delta[n], y10.shape, array_index)
            replay_center = y10 - increment

            material = data2d.random_material(generator.fine_centers, index, labels[index - 1])
            model = ExactActive(geometry, material)
            inverse, inv_cost = model.inverse_bound()
            exact_field, field_remainder, field_cost = model.enclose_field(inverse, 128)
            center_error = (norm(ball(replay_center) - exact_field).upper() + field_remainder).upper()

            y_pair, radius_pair, band_receipt = [], [], []
            for pct in NOISE_PCT:
                y = y10 - (1 - pct / 10) * increment
                noisy_center = ball(y)
                nominal = norm(noisy_center - ball(replay_center)).upper()
                radius = high(nominal + guard * norm(noisy_center))
                actual_upper = (norm(noisy_center - exact_field).upper() + field_remainder).upper()
                margin = low(arb(radius) - actual_upper)
                if not margin > 0:
                    raise ArithmeticError(f"{object_id} {pct}% noise ball was not admitted by the native solve.")
                path_defect = norm(
                    ball(y10) - noisy_center - ball(increment) * (1 - arb(pct) / 10)
                ).upper()
                y_pair.append(y)
                radius_pair.append(radius)
                band_receipt.append(dict(
                    noise_pct=pct, sample_noise_norm_upper=high(nominal), padded_radius=radius,
                    padded_ball_slack_lower=margin, affine_rounding_defect_upper=high(path_defect),
                ))
            obs.append(y_pair)
            radii.append(radius_pair)
            receipt.append(dict(
                object_id=object_id, anchor_id=expected_ids[n], center_replay_error_upper=high(center_error),
                forward_solve_timing={k: field_cost[k] for k in ("wall_seconds", "cpu_seconds")},
                inverse_bound_timing={k: inv_cost[k] for k in ("wall_seconds", "cpu_seconds")},
                bands=band_receipt, **timing(item_t),
            ))
            print(object_id, "verified 1%/2%", "seconds", round(receipt[-1]["wall_seconds"], 2), flush=True)
            del material, model, exact_field, replay_center, increment
            gc.collect()

    ids = np.array([f"s2d_{index:03d}" for index in SCOPE])
    path = DATA / "srap_supplement.npz"
    tmp_npz = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp_npz, object_id=ids, Y=np.array(obs), radius=np.array(radii),
                        noise_pct=np.array(NOISE_PCT))
    tmp_npz.replace(path)
    doc = dict(
        scope="Supplemental holdout noise sensitivity after the primary freeze, 22 same objects",
        object_count=len(ids), observation_count=2 * len(ids), noise_percent=list(NOISE_PCT),
        material_generator=str((SRAP / "data2d.py").relative_to(PILOT.parents[1])),
        anchor_observation_pack=str(public.relative_to(PILOT.parents[1])),
        native_operator=str(FINE.relative_to(PILOT.parents[1])),
        geometry_replay_binary64_equal=matches,
        noise_construction=(
            "Recreate the original complex Gaussian draw with the original archive seed; "
            "normalize to the archived 10% increment norm; form Y_p = Y_10 - (1-p/10)*increment."
        ),
        noise_seed_entropy_prefix=data2d.NOISE_SEED, confirmation_array_index=array_index,
        radius_construction="Outward noise norm plus 1e-6 of the noisy field norm.",
        radius_audit="Original isolated ExactActive material solve, Arb inverse and corrected-port enclosures.",
        interval_precision_bits=ctx.prec, interval_threads=ctx.threads,
        setup=setup, cases=receipt, **timing(t),
    )
    tmp_json = path.with_suffix(".tmp.json")
    tmp_json.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n")
    tmp_json.replace(path.with_suffix(".json"))


if __name__ == "__main__":
    main()
