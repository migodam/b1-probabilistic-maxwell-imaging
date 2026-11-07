"""Public-anchor-only inference: no manifest, labels, or truth paths are read."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
import numpy as np
import torch

try:
    from .networks import NC, NS, OperatorCache, ResidualUNet, choose_device, normalized_features, proposal
except ImportError:
    from networks import NC, NS, OperatorCache, ResidualUNet, choose_device, normalized_features, proposal


def _sync(device):
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elif device.type == "mps":
        torch.mps.synchronize()


def _tensor(value, device, complex_value=False):
    value = torch.from_numpy(np.asarray(value, dtype=np.complex64 if complex_value else np.float32).copy())
    if complex_value and device.type == "mps":
        value = torch.view_as_real(value)
    return value.to(device)


def _numpy(value, current=False):
    array = value.detach().cpu().numpy()
    if current and not np.iscomplexobj(array):
        array = array[..., 0]+1j*array[..., 1]
    return array


def _read_optional_json(npz, key):
    if key not in npz:
        return None
    value = npz[key]
    if value.ndim == 0:
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        return json.loads(value)
    return None


def predict_anchors(checkpoint_path, anchor_paths, operators_path, output_dir, device="auto", operators=None):
    """Write one NPZ and cost sidecar per public anchor, returning receipts.

    Setup is invoiced separately once. Incremental preparation includes the
    anchor's explicit feature_receipt (if present), input normalization and the
    network/rowLS proposal, but never the base exact-Maxwell anchor construction.
    Missing feature-cost evidence is marked incomplete rather than guessed.
    """
    target_device = choose_device(device)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint.get("format") != "gn_nn_r3_v1":
        raise ValueError("Unrecognized checkpoint format")
    arm = checkpoint["config"]["arm"]
    passes = checkpoint["config"]["cascade_passes"]
    if checkpoint.get("r4_format") != "gn_nn_r4_preservation_v1":
        raise ValueError("R4 inference refuses old R3 checkpoints")
    alpha = float(checkpoint["config"].get("alpha", .3))
    stats = checkpoint["feature_stats"]
    model = ResidualUNet(arm).to(target_device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    operators_reused = operators is not None
    operators = operators if operators_reused else OperatorCache(operators_path).to(target_device)
    if operators.G_re.device != target_device:
        raise ValueError("Shared operator device mismatch")
    _sync(target_device)
    setup = {"checkpoint_and_operator_setup_seconds": time.perf_counter()-start,
             "V48_source": operators.v48_source, "device": str(target_device),
             "checkpoint": str(checkpoint_path), "operators": str(operators_path),
             "arm": arm, "cascade_passes": passes, "alpha": alpha, "loss": checkpoint["config"]["loss"],
             "projectors": operators.projector_metadata, "operators_reused": operators_reused,
             "base_anchor_cost_included": False}
    setup["dev_only"] = bool(checkpoint["config"].get("dev_only", False))
    (output_dir/"inference_setup.json").write_text(json.dumps(setup, indent=2)+"\n")
    receipts = []
    for path in anchor_paths:
        path = Path(path)
        obs_id = path.stem
        load_start = time.perf_counter()
        with np.load(path, allow_pickle=False) as p:
            anchor = {key: np.array(p[key]) for key in ("m0", "j0", "Q0", "Y", "gb", "noise")}
            feature_receipt = _read_optional_json(p, "feature_receipt")
            base_anchor_cost = _read_optional_json(p, "base_anchor_cost")
        anchor_sidecar = path.with_suffix(".json")
        if anchor_sidecar.exists():
            recorded = json.loads(anchor_sidecar.read_text())
            if feature_receipt is None and "feature_wall_seconds" in recorded and "feature_counts" in recorded:
                feature_receipt = {"wall_seconds": recorded["feature_wall_seconds"],
                                   **recorded["feature_counts"], "complete": True,
                                   "source": str(anchor_sidecar)}
            if base_anchor_cost is None:
                base_anchor_cost = {"solve_wall_seconds": recorded.get("solve_wall_seconds"),
                                    "counts": recorded.get("counts"),
                                    "source": str(anchor_sidecar)}
        anchor_io_seconds = time.perf_counter()-load_start
        if anchor["j0"].shape != (NC, NS) or anchor["Q0"].shape != (NC, NS):
            raise ValueError(f"Expected cell x source current/field: {path}")
        if not all(np.isfinite(value).all() for value in anchor.values()):
            raise ValueError(f"Nonfinite public anchor {path}")
        _sync(target_device)
        feature_start = time.perf_counter()
        x = normalized_features(anchor, stats, arm)[None].to(target_device)
        m0 = _tensor(anchor["m0"].reshape(1, NC), target_device)
        j0 = _tensor(anchor["j0"][None], target_device, complex_value=True)
        _sync(target_device)
        feature_seconds = time.perf_counter()-feature_start
        inference_start = time.perf_counter()
        with torch.no_grad():
            prediction = proposal(model, operators, x, m0, j0, stats, passes, alpha)
        _sync(target_device)
        inference_seconds = time.perf_counter()-inference_start
        transfer_start = time.perf_counter()
        arrays = {"m_proposal": _numpy(prediction["m_proposal"])[0],
                  "raw_delta": _numpy(prediction["raw_delta"])[0],
                  "stage_materials": _numpy(prediction["stage_materials"])[0]}
        for key in ("j_aux", "q_aux", "delta_current", "stage_currents", "stage_fields"):
            if key in prediction:
                arrays[key] = _numpy(prediction[key], current=True)[0]
        transfer_seconds = time.perf_counter()-transfer_start
        if not all(np.isfinite(value).all() for value in arrays.values()):
            raise FloatingPointError(f"Nonfinite prediction for {obs_id}")
        feature_known = bool(feature_receipt and feature_receipt.get("complete", False))
        existing_feature_seconds = float((feature_receipt or {}).get("wall_seconds", 0.))
        cost = {"wall_seconds": existing_feature_seconds+feature_seconds+inference_seconds+transfer_seconds,
                "LU": int((feature_receipt or {}).get("LU", 0)),
                "actual_complex_RHS": int((feature_receipt or {}).get("actual_complex_RHS", 0)),
                "complete": feature_known, "feature_receipt": feature_receipt,
                "feature_normalization_seconds": feature_seconds,
                "inference_seconds": inference_seconds, "output_transfer_seconds": transfer_seconds,
                "base_anchor_included": False, "global_setup_reference": "inference_setup.json"}
        receipt = {"obs_id": obs_id, "arm": arm, "cascade_passes": passes,
                   "alpha": alpha, "loss": checkpoint["config"]["loss"],
                   "anchor_path": str(path), "preparation_cost": cost,
                   "inference_seconds": inference_seconds, "anchor_io_seconds": anchor_io_seconds,
                   "base_anchor_cost": base_anchor_cost,
                   "truth_access": "NONE", "full_vie_verification": "NOT_RUN_EXTERNAL_REQUIRED"}
        # The sidecar stores complete costs, while scalar costs also travel with NPZ.
        np.savez_compressed(output_dir/f"{obs_id}.npz", **arrays,
                            inference_seconds=inference_seconds,
                            preparation_cost_json=json.dumps(cost))
        (output_dir/f"{obs_id}.json").write_text(json.dumps(receipt, indent=2, allow_nan=False)+"\n")
        receipts.append(receipt)
    return receipts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--operators", required=True)
    parser.add_argument("--anchors", nargs="+", required=True,
                        help="Explicit public anchor NPZ files or one directory")
    parser.add_argument("--output", required=True)
    parser.add_argument("--device", default="auto")
    args = parser.parse_args()
    paths = [Path(p) for p in args.anchors]
    if len(paths) == 1 and paths[0].is_dir():
        paths = sorted(paths[0].glob("*.npz"))
    if not paths:
        raise ValueError("No public anchors supplied")
    results = predict_anchors(args.checkpoint, paths, args.operators, args.output, args.device)
    print(json.dumps({"status": "PREDICTED", "count": len(results), "output": args.output,
                      "truth_access": "NONE"}), flush=True)


if __name__ == "__main__":
    main()
