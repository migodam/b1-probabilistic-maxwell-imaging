"""TRAIN/VAL-only fitting. Prediction lives in predict.py and never opens GT."""
from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

try:
    from .networks import (ARMS, NC, NS, MATERIAL_SCALE, OperatorCache, ResidualUNet,
                           choose_device, fit_feature_stats, normalized_features,
                           proposal, material_loss, architecture_cost)
except ImportError:
    from networks import (ARMS, NC, NS, MATERIAL_SCALE, OperatorCache, ResidualUNet,
                          choose_device, fit_feature_stats, normalized_features,
                          proposal, material_loss, architecture_cost)


def read_manifest(path):
    path = Path(path)
    if path.suffix == ".jsonl":
        records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    else:
        value = json.loads(path.read_text())
        records = value if isinstance(value, list) else value["records"]
    ids = [r["obs_id"] for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate obs_id in manifest")
    return records


def load_anchor(path):
    with np.load(path, allow_pickle=False) as p:
        result = {k: np.array(p[k]) for k in ("m0", "j0", "Q0", "Y", "gb", "noise")}
    for k in ("m0", "j0", "Q0", "Y", "gb", "noise"):
        if not np.isfinite(result[k]).all():
            raise ValueError(f"Nonfinite anchor {path}: {k}")
    if result["m0"].size != NC or result["j0"].shape != (NC, NS) or result["Q0"].shape != (NC, NS):
        raise ValueError(f"Anchor dimensions invalid: {path}")
    return result


def resolve_path(root, value):
    path = Path(value)
    return path if path.is_absolute() else Path(root)/path


class SupervisedAnchors(Dataset):
    def __init__(self, records, data_root, anchor_root, stats, arm, dev_only=False):
        self.records = list(records)
        self.data_root, self.anchor_root = Path(data_root), Path(anchor_root)
        self.stats, self.arm = stats, arm
        allowed = ("debugtrain", "debugval") if dev_only else ("train", "val", "validation")
        if any(str(r["split"]).lower() not in allowed for r in self.records):
            raise ValueError("Supervised loader refuses every split except TRAIN/VAL")

    def __len__(self):
        return len(self.records)

    def __getitem__(self, index):
        record = self.records[index]
        anchor = load_anchor(self.anchor_root/f"{record['obs_id']}.npz")
        # This is the sole GT read in the training stack, after the split guard.
        path = resolve_path(self.data_root, record.get("truth_path") or
                            f"truth/{str(record['split']).lower()}/{record['parent_id']}.npz")
        with np.load(path, allow_pickle=False) as p:
            m = np.asarray(p["m"], np.float32).reshape(NC).copy()
            j = np.asarray(p["j"], np.complex64).reshape(NC, NS).copy()
            amplitude = float(np.asarray(p["amplitude"]).reshape(-1)[0])
        if not np.isfinite(m).all() or not np.isfinite(j).all() or not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError(f"Nonfinite TRAIN/VAL truth {path}")
        return {"x": normalized_features(anchor, self.stats, self.arm),
                "m0": torch.from_numpy(np.asarray(anchor["m0"], np.float32).reshape(NC)),
                "j0": torch.from_numpy(np.asarray(anchor["j0"], np.complex64)),
                "m_gt": torch.from_numpy(m), "j_gt": torch.from_numpy(j),
                "amplitude": torch.tensor(amplitude, dtype=torch.float32)}


def _json_write(path, value):
    path = Path(path)
    temp = path.with_name(path.name+".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temp.replace(path)


def _sync(device):
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elif device.type == "mps":
        torch.mps.synchronize()


def to_device(value, device):
    if device.type == "mps" and value.is_complex():
        value = torch.view_as_real(value)
    return value.to(device)


def _run_epoch(model, operators, loader, stats, passes, device, optimizer=None):
    training = optimizer is not None
    model.train(training)
    sums = {"loss": 0., "material_l1_normalized": 0., "edge_l1_normalized": 0.,
            "material_mae": 0., "material_relative_l2": 0., "material_amplitude_nrmse": 0.,
            "current_l1_normalized": 0.}
    count = 0
    gradient_receipts = []
    with torch.set_grad_enabled(training):
        for batch in loader:
            batch = {key: to_device(value, device) for key, value in batch.items()}
            if training:
                optimizer.zero_grad(set_to_none=True)
            output = proposal(model, operators, batch["x"], batch["m0"], batch["j0"], stats, passes)
            mloss, terms = material_loss(output["m_proposal"], batch["m_gt"])
            closs = torch.zeros((), device=device)
            if model.arm == "j_null_current":
                target = operators.null_project(batch["j_gt"]-batch["j0"])
                difference = output["delta_current"]-target
                magnitude = difference.abs() if difference.is_complex() else torch.linalg.vector_norm(difference, dim=-1)
                closs = magnitude.mean()/float(stats["current_scale"])
                loss = closs
            else:
                loss = mloss
            if not torch.isfinite(loss):
                raise FloatingPointError("Nonfinite training/validation loss")
            if training:
                loss.backward()
                head_grad = torch.sqrt(sum(p.grad.square().sum() for p in model.head.parameters() if p.grad is not None))
                stem_grad = torch.sqrt(sum(p.grad.square().sum() for p in model.stem.parameters() if p.grad is not None))
                gradient_receipts.append({"head_l2": float(head_grad.detach()), "stem_l2": float(stem_grad.detach())})
                optimizer.step()
            n = batch["m_gt"].shape[0]
            err = output["m_proposal"]-batch["m_gt"]
            nrmse = torch.linalg.vector_norm(err, dim=1)/torch.linalg.vector_norm(batch["m_gt"], dim=1).clamp_min(1e-12)
            amplitude_nrmse = err.square().mean(dim=1).sqrt()/batch["amplitude"].clamp_min(1e-12)
            values = {"loss": loss, **terms, "material_mae": err.abs().mean(),
                      "material_relative_l2": nrmse.mean(), "material_amplitude_nrmse": amplitude_nrmse.mean(),
                      "current_l1_normalized": closs}
            for key, value in values.items():
                sums[key] += float(value.detach())*n
            count += n
    result = {key: value/count for key, value in sums.items()}
    if gradient_receipts:
        result["gradient_first_step"] = gradient_receipts[0]
        result["gradient_last_step"] = gradient_receipts[-1]
        result["optimizer_steps"] = len(gradient_receipts)
    return result


def train_model(data_root, manifest_path, operators_path, output_dir, arm,
                seed=7, epochs=100, batch_size=64, patience=15, device="auto",
                cascade_passes=1, anchor_root=None, num_workers=0, dev_only=False):
    if arm not in ARMS:
        raise ValueError(arm)
    if epochs < 1 or epochs > 100 or batch_size < 1 or patience < 1:
        raise ValueError("Require 1<=epochs<=100, positive batch_size and patience")
    if cascade_passes == 2 and arm not in ("v_physics", "j_null"):
        raise ValueError("Cascade applies only to v_physics/j_null")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    target_device = choose_device(device)
    data_root, output_dir = Path(data_root), Path(output_dir)
    anchor_root = Path(anchor_root) if anchor_root else data_root.parent/"anchors"
    output_dir.mkdir(parents=True, exist_ok=True)
    records = read_manifest(manifest_path)
    train_splits = ("debugtrain",) if dev_only else ("train",)
    val_splits = ("debugval",) if dev_only else ("val", "validation")
    train = [r for r in records if str(r["split"]).lower() in train_splits]
    val = [r for r in records if str(r["split"]).lower() in val_splits]
    if not train or not val:
        raise ValueError("Both nonempty TRAIN and VAL records required")
    if set(r["parent_id"] for r in train) & set(r["parent_id"] for r in val):
        raise ValueError("Parent leakage between TRAIN and VAL")
    start = time.perf_counter()
    stats = fit_feature_stats(load_anchor(anchor_root/f"{r['obs_id']}.npz") for r in train)
    if dev_only:
        stats["fitted_on"] = "DEV_ONLY_debugtrain_public_anchors_only"
    stats_seconds = time.perf_counter()-start
    model = ResidualUNet(arm).to(target_device)
    setup_start = time.perf_counter()
    operators = OperatorCache(operators_path).to(target_device)
    _sync(target_device)
    operator_setup_seconds = time.perf_counter()-setup_start
    generator = torch.Generator().manual_seed(seed)
    train_loader = DataLoader(SupervisedAnchors(train, data_root, anchor_root, stats, arm, dev_only),
                              batch_size=batch_size, shuffle=True, generator=generator,
                              num_workers=num_workers)
    val_loader = DataLoader(SupervisedAnchors(val, data_root, anchor_root, stats, arm, dev_only),
                            batch_size=batch_size, shuffle=False, num_workers=num_workers)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    config = {"arm": arm, "seed": seed, "max_epochs": epochs, "batch_size": batch_size,
              "patience": patience, "lr": 1e-3, "optimizer": "AdamW", "weight_decay": .01,
              "cascade_passes": cascade_passes, "device": str(target_device),
              "requested_device": device, "train_count": len(train), "val_count": len(val),
              "input_channels": 27, "unet_widths": [16, 32, 64],
              "complex_dtype": "float32_real_pair" if target_device.type == "mps" else "complex64", "network_dtype": "float32",
              "selection_metric": "val_current_l1" if arm == "j_null_current" else "val_material_l1_plus_0.1_edge",
              "feature_stats_seconds": stats_seconds, "operator_setup_seconds": operator_setup_seconds,
              "V48_source": operators.v48_source, "truth_access": "DEV_ONLY_EXPLICITLY_AUTHORIZED" if dev_only else "TRAIN_VAL_ONLY",
              "dev_only": dev_only,
              "residual_head_initialization": "all_zero_weights_and_bias",
              "offline_metric_definitions": {"material_relative_l2": "||m_proposal-m_GT||_2/||m_GT||_2",
                                              "material_amplitude_nrmse": "sqrt(mean((m_proposal-m_GT)^2))/GT_amplitude"},
              "arithmetic_estimate": architecture_cost(arm, cascade_passes),
              "full_vie_verification": "NOT_RUN_EXTERNAL_REQUIRED"}
    _json_write(output_dir/"config.json", config)
    best, stale, history = float("inf"), 0, []
    fit_start = time.perf_counter()
    for epoch in range(1, epochs+1):
        epoch_start = time.perf_counter()
        tm = _run_epoch(model, operators, train_loader, stats, cascade_passes, target_device, optimizer)
        vm = _run_epoch(model, operators, val_loader, stats, cascade_passes, target_device)
        _sync(target_device)
        row = {"epoch": epoch, "train": tm, "val": vm, "seconds": time.perf_counter()-epoch_start}
        history.append(row)
        if vm["loss"] < best:
            best, stale = vm["loss"], 0
            checkpoint = {"format": "gn_nn_r3_v1", "state_dict": model.state_dict(),
                          "feature_stats": stats, "config": config, "best_epoch": epoch,
                          "best_val": vm, "parameter_count": sum(p.numel() for p in model.parameters())}
            temp = output_dir/"best.pt.tmp"
            torch.save(checkpoint, temp)
            temp.replace(output_dir/"best.pt")
        else:
            stale += 1
        _json_write(output_dir/"history.json", history)
        print(json.dumps(row), flush=True)
        if stale >= patience:
            break
    summary = {"status": "TRAINED", "checkpoint": str(output_dir/"best.pt"),
               "best_loss": best, "epochs_run": len(history),
               "early_stopped": stale >= patience, "training_seconds": time.perf_counter()-fit_start,
               "feature_stats_seconds": stats_seconds, "operator_setup_seconds": operator_setup_seconds,
               "best_epoch": int(torch.load(output_dir/"best.pt", map_location="cpu", weights_only=False)["best_epoch"]),
               "config": config}
    _json_write(output_dir/"summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--operators", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--anchors")
    parser.add_argument("--arm", choices=ARMS, required=True)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--cascade-passes", type=int, choices=(1, 2), default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--dev-only", action="store_true")
    args = parser.parse_args()
    result = train_model(args.data_root, args.manifest, args.operators, args.output, args.arm,
                         args.seed, args.epochs, args.batch_size, args.patience, args.device,
                         args.cascade_passes, args.anchors, args.num_workers, args.dev_only)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
