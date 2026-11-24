"""Score threshold area on the original reference arrays after the proofs."""

import argparse
from pathlib import Path

import numpy as np
from flint import arb, ctx

from area_certificate import (HERE, NATIVE_ID, OPTICAL, RESEARCH, RING_ID, as_path, clock,
                              read_json, save_json, timing)

RING_REFERENCE = (
    RESEARCH / "delegated/max_struct_r1/r10_worktree/research/research_max_struct_r10"
    / "sealed_truth/twod/r10_twod_0001.npz"
)
NATIVE_REFERENCE = OPTICAL / "synthetic/evaluation/optfresh_01.npz"
RING_FORWARD_AUDIT = OPTICAL / "submission/diagnostic/interval_and_offline.json"
NATIVE_FORWARD_AUDIT = OPTICAL / "results/theory/fine_feasibility.json"


def score(path, key, threshold, shape):
    with np.load(path, allow_pickle=False) as data:
        m = data[key]
    if m.size != shape or not np.isfinite(m).all() or np.min(m) < 0:
        raise ValueError("The offline reference does not match the native contrast contract.")
    count = int(np.count_nonzero(m > threshold))
    return dict(source=as_path(path), array=key, threshold=threshold, pixels=shape,
                material_peak_binary64=float(np.max(m)), true_high_pixels=count,
                true_high_pixel_fraction=count / shape)


def compare(radius_table, interval):
    floor, ceiling = map(arb, interval)
    result = {}
    for key, row in radius_table.items():
        lower, upper = arb(row["radius_lower"]), arb(row["radius_upper"])
        result[key] = "inside" if ceiling <= lower else "outside" if floor > upper else "unresolved"
    return result


def run(online, output):
    started = clock()
    ctx.prec = 120
    online = read_json(online)
    ring, native = (online["cases"][key] for key in (RING_ID, NATIVE_ID))
    ring_result = score(RING_REFERENCE, "m", ring["threshold"], 576)
    native_result = score(NATIVE_REFERENCE, "m_fine", native["threshold"], 9216)
    audited_ring = read_json(RING_FORWARD_AUDIT)[RING_ID]["offline_reference"]["reference_interval"]
    audited_native = next(r for r in read_json(NATIVE_FORWARD_AUDIT)["observations"] if r["id"] == NATIVE_ID)
    ring_interval = [audited_ring["full_wave_distance_lo"], audited_ring["full_wave_distance_hi"]]
    native_interval = [audited_native["true_forward_to_Y_lower"], audited_native["true_forward_to_Y_upper"]]
    ring_result.update(forward_distance_interval=ring_interval, forward_audit=as_path(RING_FORWARD_AUDIT),
                       radius_compatibility=compare(ring["budgets"], ring_interval))
    native_result.update(forward_distance_interval=native_interval, forward_audit=as_path(NATIVE_FORWARD_AUDIT),
                         radius_compatibility=compare(native["budgets"], native_interval))
    receipt = dict(
        scope="offline count and reuse of archived true-solver audits; no truth used for certificates",
        new_material_forward_solves=0, cases={RING_ID: ring_result, NATIVE_ID: native_result},
        **timing(started),
    )
    save_json(output, receipt)
    print(RING_ID, ring_result["true_high_pixels"], ring_result["radius_compatibility"])
    print(NATIVE_ID, native_result["true_high_pixels"], native_result["radius_compatibility"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--online", type=Path, default=HERE / "area_receipt.json")
    parser.add_argument("--output", type=Path, default=HERE / "truth_offline.json")
    args = parser.parse_args()
    run(args.online, args.output)
