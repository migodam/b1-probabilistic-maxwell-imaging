"""Public-anchor-only recurrent inference; never opens manifests, labels or truth."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np
import torch

try:
    from .networks import NC, NS, OperatorCache, choose_device, normalized_features
    from .recurrent_feedback import RecurrentFeedback, HINGE_CLEAN_DELTA, HINGE_NOISY_DELTA
    from .predict import _sync, _tensor, _numpy, _read_optional_json
except ImportError:
    from networks import NC, NS, OperatorCache, choose_device, normalized_features
    from recurrent_feedback import RecurrentFeedback, HINGE_CLEAN_DELTA, HINGE_NOISY_DELTA
    from predict import _sync, _tensor, _numpy, _read_optional_json


def predict_recurrent_anchors(checkpoint_path, anchor_paths, operators_path, output_dir,
                              device="mps", *, operators=None, stats=None):
    """Write evaluation-compatible candidates and complete per-case cost sidecars.

    Anchor paths are explicit public NPZ files. Optional immutable operators and
    TRAIN normalization are reused from the trainer; neither is refit here.
    Setup cost is recorded once separately; base-anchor solving is excluded.
    """
    target = choose_device(device)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint.get("format") != "gn_nn_r4_recurrent_feedback_v1":
        raise ValueError("Expected an independently trained R4 recurrent checkpoint")
    config = checkpoint["config"]
    stats_reused = stats is not None
    stats = checkpoint["feature_stats"] if stats is None else stats
    model = RecurrentFeedback(config["mode"], config["rho"]).to(target)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    operators_reused = operators is not None
    operators = OperatorCache(operators_path).to(target) if operators is None else operators
    _sync(target)
    setup = {"checkpoint_and_operator_setup_seconds": time.perf_counter()-started,
             "checkpoint": str(checkpoint_path), "operators": str(operators_path),
             "device": str(target), "arm": config["arm"], "mode": config["mode"],
             "rho": config["rho"], "cascade_passes": 2,
             "feedback_step_rule": "hinge_majorizer" if config["mode"] == "feedback_hinge" else "quadratic_exact_line",
             "noise_added_to_NN_channels": False,
             "operators_reused": operators_reused, "TRAIN_feature_stats_reused": stats_reused,
             "new_operator_uploads": int(not operators_reused),
             "projectors": operators.projector_metadata, "base_anchor_cost_included": False,
             "truth_access": "NONE", "class_or_evaluator_access": "NONE"}
    setup_name = f"inference_setup_{time.time_ns()}.json"
    setup["receipt_filename"] = setup_name
    setup_text = json.dumps(setup, indent=2)+"\n"
    (output_dir/setup_name).write_text(setup_text)
    (output_dir/"inference_setup.json").write_text(setup_text)
    receipts = []
    for anchor_path in anchor_paths:
        path = Path(anchor_path)
        obs_id = path.stem
        io_started = time.perf_counter()
        with np.load(path, allow_pickle=False) as p:
            anchor = {key: np.array(p[key]) for key in ("m0", "j0", "Q0", "Y", "gb", "noise")}
            feature_receipt = _read_optional_json(p, "feature_receipt")
            base_anchor_cost = _read_optional_json(p, "base_anchor_cost")
        sidecar = path.with_suffix(".json")
        if sidecar.exists():
            recorded = json.loads(sidecar.read_text())
            if feature_receipt is None and "feature_wall_seconds" in recorded and "feature_counts" in recorded:
                feature_receipt = {"wall_seconds": recorded["feature_wall_seconds"],
                                   **recorded["feature_counts"], "complete": True, "source": str(sidecar)}
            if base_anchor_cost is None:
                base_anchor_cost = {"solve_wall_seconds": recorded.get("solve_wall_seconds"),
                                    "counts": recorded.get("counts"), "source": str(sidecar)}
        anchor_io_seconds = time.perf_counter()-io_started
        if anchor["j0"].shape != (NC, NS) or anchor["Q0"].shape != (NC, NS):
            raise ValueError(f"Expected public cell-by-source current/field: {path}")
        if not all(np.isfinite(value).all() for value in anchor.values()):
            raise ValueError(f"Nonfinite public anchor: {path}")
        _sync(target)
        feature_started = time.perf_counter()
        features = normalized_features(anchor, stats, "j_full")[None].to(target)
        m0 = _tensor(anchor["m0"].reshape(1, NC), target)
        j0 = _tensor(anchor["j0"][None], target, complex_value=True)
        q0 = _tensor(anchor["Q0"][None], target, complex_value=True)
        public_y = _tensor(anchor["Y"][None], target, complex_value=True)
        public_noise = _tensor(np.asarray(anchor["noise"]).reshape(1), target)
        _sync(target)
        feature_seconds = time.perf_counter()-feature_started
        inference_started = time.perf_counter()
        with torch.no_grad():
            output = model(operators, features, m0, j0, public_y, stats, q0, public_noise=public_noise)
        _sync(target)
        inference_seconds = time.perf_counter()-inference_started
        transfer_started = time.perf_counter()
        arrays = {"m_proposal": _numpy(output["m_proposal"])[0],
                  "raw_delta": _numpy(output["raw_delta"])[0],
                  "stage_materials": _numpy(output["stage_materials"])[0]}
        for key in ("j_aux", "q_aux", "delta_current", "stage_currents", "stage_fields"):
            arrays[key] = _numpy(output[key], current=True)[0]
        diagnostics = {}
        if model.feedback_steps:
            for key in ("tau", "gradient_norm2", "curvature", "phi_before", "phi_after",
                        "zero_gradient", "zero_curvature", "zero_y_normalization", "zero_j0_normalization",
                        "data_tolerance", "data_theta", "relative_data_residual_before", "relative_data_residual_after"):
                value = torch.stack([torch.stack([step[key] for step in layer], dim=1)
                                     for layer in output["feedback_diagnostics"]], dim=1)
                diagnostics["feedback_"+key] = _numpy(value)[0]
        transfer_seconds = time.perf_counter()-transfer_started
        if not all(np.isfinite(value).all() for value in {**arrays, **diagnostics}.values()):
            raise FloatingPointError(f"Nonfinite recurrent candidate for {obs_id}")
        fixed_counts = output["preparation_counts"]
        known = bool(feature_receipt and feature_receipt.get("complete", False))
        cost = {"wall_seconds": float((feature_receipt or {}).get("wall_seconds", 0.))
                                +feature_seconds+inference_seconds+transfer_seconds,
                "LU": int((feature_receipt or {}).get("LU", 0))+int(fixed_counts["LU"]),
                "actual_complex_RHS": int((feature_receipt or {}).get("actual_complex_RHS", 0))
                                      +int(fixed_counts["actual_complex_RHS"]),
                "complete": known, "feature_receipt": feature_receipt,
                "feature_normalization_seconds": feature_seconds, "inference_seconds": inference_seconds,
                "output_transfer_seconds": transfer_seconds, "fixed_operator_cost": fixed_counts,
                "base_anchor_included": False, "global_setup_reference": setup_name}
        receipt = {"obs_id": obs_id, "arm": config["arm"], "mode": config["mode"],
                   "rho": config["rho"], "cascade_passes": 2,
             "feedback_step_rule": "hinge_majorizer" if config["mode"] == "feedback_hinge" else "quadratic_exact_line",
             "noise_added_to_NN_channels": False, "anchor_path": str(path),
                   "preparation_cost": cost, "fixed_operator_cost": fixed_counts,
                   "inference_seconds": inference_seconds, "anchor_io_seconds": anchor_io_seconds,
                   "base_anchor_cost": base_anchor_cost, "truth_access": "NONE",
                   "class_or_evaluator_access": "NONE", "full_vie_verification": "NOT_RUN_EXTERNAL_REQUIRED"}
        np.savez_compressed(output_dir/f"{obs_id}.npz", **arrays, **diagnostics,
                            inference_seconds=inference_seconds,
                            fixed_operator_cost_json=json.dumps(fixed_counts),
                            preparation_cost_json=json.dumps(cost))
        (output_dir/f"{obs_id}.json").write_text(json.dumps(receipt, indent=2, allow_nan=False)+"\n")
        receipts.append(receipt)
    return receipts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--operators", type=Path, required=True)
    p.add_argument("--anchors", nargs="+", type=Path, required=True,
                   help="Explicit public anchor NPZ paths; no manifest or evaluator file")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--device", default="mps")
    args = p.parse_args()
    result = predict_recurrent_anchors(args.checkpoint, args.anchors, args.operators, args.output, args.device)
    print(json.dumps({"status": "PREDICTED", "count": len(result), "output": str(args.output),
                      "truth_access": "NONE", "class_or_evaluator_access": "NONE"}), flush=True)


if __name__ == "__main__":
    main()
