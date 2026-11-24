"""Score sealed native-lifted synthetic predictions on private material and ports."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage
from threadpoolctl import threadpool_limits

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.z_alg_xudong_cie.holdout.data import (
    DEFAULT_ROOT, INVERSE_SIDE, NATIVE_SIDE, OBJECT_IDS, elapsed, get_native_case, planned_cases,
    require_freeze, timers, write_json,
)
from research.z_alg_xudong_cie.holdout.physics import NativeForward, NativeForwardError, nearest_lift


def material_scores(m, truth):
    """SRAP contract: one truth norm and one truth half-peak for both maps."""
    if m.shape != truth.shape or not np.isfinite(truth).all():
        raise ValueError("Candidate and reference must share their declared grid.")
    d = m - truth
    reference_norm = float(np.linalg.norm(truth))
    peak = float(truth.max())
    if not (reference_norm > 0 and peak > 0):
        raise ValueError("The blind material reference is empty.")
    threshold = peak / 2
    a = m >= threshold
    b = truth >= threshold
    shape = (int(np.sqrt(m.size)),) * 2

    def topology(mask):
        bitmap = mask.reshape(shape)
        outer = ndimage.label(bitmap)[1]
        holes = ndimage.label(ndimage.binary_fill_holes(bitmap) & ~bitmap)[1]
        return int(outer), int(holes)

    cn, hn = topology(a)
    ct, ht = topology(b)
    return dict(rel_L2=float(np.linalg.norm(d) / reference_norm),
                inside_L2=float(np.linalg.norm(d[truth > 0]) / reference_norm),
                background_L2=float(np.linalg.norm(d[truth == 0]) / reference_norm),
                amplitude_NRMSE=float(np.sqrt(np.mean(d**2)) / peak),
                half_IoU=float(np.count_nonzero(a & b) / np.count_nonzero(a | b)),
                half_peak_threshold=threshold, material_reference_norm=reference_norm,
                foreground_mean_ratio=float(np.mean(m[b]) / np.mean(truth[b])),
                peak_contrast=float(np.max(m)), truth_peak_contrast=peak,
                holes=hn, true_holes=ht, components=cn, true_components=ct)


def relative_error(estimate, measured):
    dist = float(np.linalg.norm(estimate - measured))
    scale = float(np.linalg.norm(measured))
    if not scale > 0:
        raise ValueError("Held or fitting field must have nonzero norm.")
    return dict(absolute_Frobenius=dist, relative_to_measured=dist / scale, field_norm=scale)


def score(identifier, candidate, *, freeze_record, root=DEFAULT_ROOT):
    freeze = require_freeze(freeze_record)
    start = timers()
    rows = [row for row in planned_cases()
            if row["case"] == identifier or (identifier in OBJECT_IDS and row["object_id"] == identifier)]
    if not rows:
        raise ValueError("Specify a blind material or one of its exact fitting observation IDs.")
    oid = rows[0]["object_id"]
    m = np.asarray(candidate).reshape(-1)
    if (m.size != INVERSE_SIDE**2 or not np.all(np.isfinite(m)) or np.iscomplexobj(m) or np.any(m < 0)):
        raise ValueError("The sealed prediction must contain one real nonnegative 32² material.")
    up = nearest_lift(m)
    offline = Path(root) / "offline"
    with np.load(offline / "objects" / f"{oid}.npz", allow_pickle=False) as stored:
        truth32, truth96, held_response = (stored[k] for k in ("m_coarse", "m_native", "clean_held"))
    with np.load(offline / "hidden_geometry96.npz", allow_pickle=False) as stored:
        hidden_s, hidden_b = stored["S"], stored["B"]
    fits = [(row["case"], row["array"], get_native_case(row["case"], root)[1]) for row in rows]
    named = {ap: (bundle["S"], bundle["E"]) for _, ap, bundle in fits}
    named["held"] = hidden_s, hidden_b
    spectrum = fits[0][2]["spectrum"]
    s, sr = np.unique(np.concatenate([value[0] for value in named.values()]), axis=0, return_inverse=True)
    b, br = np.unique(np.concatenate([value[1].T for value in named.values()]), axis=0, return_inverse=True)
    rx = dict(zip(named, np.split(sr, np.arange(1, len(named))*20), strict=True))
    tx = dict(zip(named, np.split(br, np.arange(1, len(named))*10), strict=True))
    solver = NativeForward(spectrum, NATIVE_SIDE)
    base = dict(
        scope="offline predictions evaluated with the original native physics under 32-to-96 lift",
        material_grid_for_primary_comparison=NATIVE_SIDE,
        identifier=identifier, object_id=oid, freeze_record=str(freeze),
        material32=material_scores(m, truth32), material96=material_scores(up, truth96),
    )
    try:
        ports, work = solver.forward(up, b.T, s)
    except NativeForwardError as exc:
        return dict(**base, native_fitting={}, held_transmit_and_receive=None,
                    native_forward=dict(status="unresolved", reason=str(exc), **exc.cost),
                    **elapsed(start))
    native_fitting = {}
    for name, ap, bundle in fits:
        errors = relative_error(ports[np.ix_(rx[ap], tx[ap])], bundle["Y"])
        delta = float(bundle["delta"])
        native_fitting[name] = dict(**errors, padded_radius=delta,
                                    numerical_margin=delta-errors["absolute_Frobenius"])
    held = relative_error(ports[np.ix_(rx["held"], tx["held"])], held_response)
    return dict(**base, native_fitting=native_fitting, held_transmit_and_receive=held,
                native_forward=dict(status="scored", **work), **elapsed(start))


def score_batch(manifest, output, *, freeze_record, root=DEFAULT_ROOT):
    require_freeze(freeze_record)
    start = timers()
    raw = json.loads(manifest.read_text(encoding="utf-8"))["predictions"]
    records = []
    for row in raw:
        row_start = timers()
        path = Path(row["path"])
        if not path.is_absolute():
            path = manifest.parent / path
        try:
            with np.load(path, allow_pickle=False) as saved:
                candidate = saved[row.get("array", "m")]
            if "index" in row:
                candidate = candidate[row["index"]]
            result = score(row["case"], candidate, freeze_record=freeze_record, root=root)
            status = "scored" if result["native_forward"]["status"] == "scored" else "material_only"
            records.append(dict(status=status, source=row, source_resolved_path=str(path.resolve()), **result))
        except (ArithmeticError, RuntimeError, OSError, ValueError, KeyError) as exc:
            records.append(dict(status="unresolved", source=row, reason=str(exc), **elapsed(row_start)))
        write_json(output, dict(rows=records, **elapsed(start)))
        print(row["case"], row.get("arm", ""), records[-1]["status"], flush=True)
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--freeze-record", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(1):
        score_batch(args.predictions, args.output, freeze_record=args.freeze_record, root=args.root)
