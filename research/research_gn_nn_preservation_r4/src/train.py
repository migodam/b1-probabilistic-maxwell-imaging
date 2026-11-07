"""TRAIN/VAL-only fitting. Prediction lives in predict.py and never opens GT."""
from __future__ import annotations

import argparse
import json
import random
import time
from collections import OrderedDict
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

try:
    from .networks import (ARMS, NC, NS, MATERIAL_SCALE, OperatorCache, ResidualUNet,
                           choose_device, fit_feature_stats, normalized_features,
                           proposal, material_loss, architecture_cost, public_regularizer)
except ImportError:
    from networks import (ARMS, NC, NS, MATERIAL_SCALE, OperatorCache, ResidualUNet,
                          choose_device, fit_feature_stats, normalized_features,
                          proposal, material_loss, architecture_cost, public_regularizer)


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
    def __init__(self, records, data_root, anchor_root, stats, arm, dev_only=False, cache_max_bytes=512*1024**2):
        self.records = list(records)
        self.data_root, self.anchor_root = Path(data_root), Path(anchor_root)
        self.stats, self.arm = stats, arm
        allowed = ("debugtrain", "debugval") if dev_only else ("train", "val", "validation")
        if any(str(r["split"]).lower() not in allowed for r in self.records):
            raise ValueError("Supervised loader refuses every split except TRAIN/VAL")
        # CPU tensors only. The initial 640-observation cohort fits comfortably
        # under the shared 512 MiB budget. Larger cohorts use a bounded LRU.
        bytes_per_item = 27*NC*4 + 2*NC*4 + 2*NC*NS*8 + 4
        self.max_cached_items = max(0, int(cache_max_bytes)//bytes_per_item)
        self._cache = OrderedDict()
        start = time.perf_counter()
        eager = len(self.records) <= self.max_cached_items
        if eager:
            for index in range(len(self.records)):
                self._cache[index] = self._read(index)
        used = sum(v.numel()*v.element_size() for item in self._cache.values() for v in item.values())
        self.cache_receipt = {"mode": "eager_cpu_all" if eager else "bounded_cpu_lru",
                              "load_seconds": time.perf_counter()-start,
                              "resident_bytes": used, "budget_bytes": int(cache_max_bytes),
                              "cached_items": len(self._cache), "total_items": len(self.records),
                              "max_cached_items": self.max_cached_items}

    def __len__(self):
        return len(self.records)

    def __getitem__(self, index):
        if index in self._cache:
            self._cache.move_to_end(index)
            return self._cache[index]
        item = self._read(index)
        if self.max_cached_items:
            if len(self._cache) >= self.max_cached_items:
                self._cache.popitem(last=False)
            self._cache[index] = item
        return item

    def _read(self, index):
        record = self.records[index]
        anchor = load_anchor(self.anchor_root/f"{record['obs_id']}.npz")
        # This is the sole GT read in the training stack, after the split guard.
        path = resolve_path(self.data_root, record.get("truth_path") or
                            f"truth/{str(record['split']).lower()}/{record['parent_id']}.npz")
        expected_truth_dir = (self.data_root/"truth"/str(record["split"]).lower()).resolve()
        if not path.resolve().is_relative_to(expected_truth_dir):
            raise ValueError("Supervised truth path must remain inside its TRAIN/VAL split")
        with np.load(path, allow_pickle=False) as p:
            m = np.asarray(p["m"], np.float32).reshape(NC).copy()
            amplitude = float(np.asarray(p["amplitude"]).reshape(-1)[0])
        if not np.isfinite(m).all() or not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError(f"Nonfinite TRAIN/VAL truth {path}")
        return {"x": normalized_features(anchor, self.stats, self.arm),
                "m0": torch.from_numpy(np.asarray(anchor["m0"], np.float32).reshape(NC)),
                "j0": torch.from_numpy(np.asarray(anchor["j0"], np.complex64)),
                "m_gt": torch.from_numpy(m),
                "q0_norm2": torch.tensor(float(np.sum(np.abs(anchor["Q0"])**2)), dtype=torch.float32),
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


def _run_epoch(model, operators, loader, stats, passes, device, optimizer=None,
               alpha=.3, loss_name="material", regularizer_weight=.01):
    training = optimizer is not None
    model.train(training)
    sums = {key: 0. for key in ("loss", "material_l1_normalized", "edge_l1_normalized",
             "material_mae", "material_relative_l2", "material_amplitude_nrmse", "public_regularizer")}
    count, gradient_receipts = 0, []
    with torch.set_grad_enabled(training):
        for batch in loader:
            batch = {key: to_device(value, device) for key, value in batch.items()}
            if training:
                optimizer.zero_grad(set_to_none=True)
            output = proposal(model, operators, batch["x"], batch["m0"], batch["j0"], stats, passes, alpha)
            mloss, terms = material_loss(output["m_proposal"], batch["m_gt"])
            rloss = public_regularizer(output, batch["j0"], batch["q0_norm2"], operators, loss_name)
            loss = mloss+regularizer_weight*rloss
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
            relative_l2 = torch.linalg.vector_norm(err, dim=1)/torch.linalg.vector_norm(batch["m_gt"], dim=1).clamp_min(1e-12)
            amplitude_nrmse = err.square().mean(dim=1).sqrt()/batch["amplitude"].clamp_min(1e-12)
            values = {"loss": loss, **terms, "material_mae": err.abs().mean(),
                      "material_relative_l2": relative_l2.mean(), "material_amplitude_nrmse": amplitude_nrmse.mean(),
                      "public_regularizer": rloss}
            for key, value in values.items():
                sums[key] += float(value.detach())*n
            count += n
    result = {key: value/count for key, value in sums.items()}
    if gradient_receipts:
        result.update(gradient_first_step=gradient_receipts[0], gradient_last_step=gradient_receipts[-1],
                      optimizer_steps=len(gradient_receipts))
    return result


def strict_training_release(study_root):
    path = Path(study_root)/"audit/DATA_READY.json"
    if not path.exists():
        return False
    try:
        r = json.loads(path.read_text())
    except (ValueError, OSError):
        return False
    return (r.get("status") == "TRAIN_VAL_STRICT_ANCHORS_READY" and r.get("train") == 1008
            and r.get("val") == 156 and r.get("all_primary_full24") is True
            and r.get("all_primary_discrepancy_stop") is False)


def released_training_cohort(data_root, manifest_path, anchor_root):
    """One launch-time cohort/provenance check, reused by all arm runs."""
    data_root, anchor_root = Path(data_root), Path(anchor_root)
    records = read_manifest(manifest_path)
    train_splits = ("train",)
    val_splits = ("val", "validation")
    train = [r for r in records if str(r["split"]).lower() in train_splits]
    val = [r for r in records if str(r["split"]).lower() in val_splits]
    if len(train) != 1008 or len(val) != 156:
        raise ValueError("Complete R4 TRAIN1008/VAL156 anchor cohorts required")
    if len(set(r["parent_id"] for r in train)) != 504 or len(set(r["parent_id"] for r in val)) != 78:
        raise ValueError("Complete R4 TRAIN504/VAL78 parent cohorts required")
    if any(r.get("digit") in (4, 7) for r in train):
        raise ValueError("Withheld digit4/7 may never enter TRAIN")
    study_root = data_root.parent
    if not strict_training_release(study_root):
        raise ValueError("Wait for data worker DATA_READY release")
    if any(not (anchor_root/f"{r['obs_id']}.npz").exists() for r in train+val):
        raise ValueError("Full TRAIN/VAL anchors required before training")
    if set(r["parent_id"] for r in train) & set(r["parent_id"] for r in val):
        raise ValueError("Parent leakage between TRAIN and VAL")
    return {"train": train, "val": val}


def train_model(data_root, manifest_path, operators_path, output_dir, arm,
                seed=7, epochs=100, batch_size=64, patience=15, device="auto",
                cascade_passes=2, anchor_root=None, num_workers=0, dev_only=False, cache_ram_mb=1536,
                alpha=.3, loss_name="material", regularizer_weight=.01, operators=None, feature_stats=None, cohort=None):
    if arm not in ARMS:
        raise ValueError(arm)
    if dev_only:
        raise ValueError("R4 forbids partial-cohort training; use synthetic wiring checks")
    if seed not in (7, 17, 27):
        raise ValueError("Frozen R4 seeds are 7/17/27")
    if epochs < 1 or epochs > 100 or batch_size < 1 or patience < 1:
        raise ValueError("Require 1<=epochs<=100, positive batch_size and patience")
    if cascade_passes != 2:
        raise ValueError("R4 requires exactly two shared-weight passes")
    if loss_name not in ("material", "const", "internal"):
        raise ValueError(loss_name)
    if arm == "v_physics" and loss_name != "material":
        raise ValueError("Voxel baseline has no learned auxiliary current")
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
    cohort = cohort or released_training_cohort(data_root, manifest_path, anchor_root)
    train, val = cohort["train"], cohort["val"]
    start = time.perf_counter()
    stats = feature_stats or fit_feature_stats(load_anchor(anchor_root/f"{r['obs_id']}.npz") for r in train)
    if dev_only:
        stats["fitted_on"] = "DEV_ONLY_debugtrain_public_anchors_only"
    stats_seconds = time.perf_counter()-start
    model = ResidualUNet(arm).to(target_device)
    setup_start = time.perf_counter()
    operators_reused = operators is not None
    operators = operators if operators_reused else OperatorCache(operators_path).to(target_device)
    _sync(target_device)
    operator_setup_seconds = time.perf_counter()-setup_start
    generator = torch.Generator().manual_seed(seed)
    total_items = len(train)+len(val)
    cache_bytes = max(0, int(cache_ram_mb*1024**2))
    train_dataset = SupervisedAnchors(train, data_root, anchor_root, stats, arm, dev_only,
                                      cache_bytes*len(train)//total_items)
    val_dataset = SupervisedAnchors(val, data_root, anchor_root, stats, arm, dev_only,
                                    cache_bytes*len(val)//total_items)
    train_loader = DataLoader(train_dataset,
                              batch_size=batch_size, shuffle=True, generator=generator,
                              num_workers=num_workers)
    val_loader = DataLoader(val_dataset,
                            batch_size=batch_size, shuffle=False, num_workers=num_workers)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    config = {"arm": arm, "seed": seed, "max_epochs": epochs, "batch_size": batch_size,
              "patience": patience, "lr": 1e-3, "optimizer": "AdamW", "weight_decay": .01,
              "cascade_passes": cascade_passes, "passes": cascade_passes, "alpha": float(alpha),
              "loss": loss_name, "regularizer_weight": regularizer_weight,
              "train_from_scratch": True, "initial_weights_source": "NONE",
              "old_checkpoint_weights_loaded": False, "operators_reused": operators_reused,
              "projectors": operators.projector_metadata, "device": str(target_device),
              "requested_device": device, "train_count": len(train), "val_count": len(val),
              "input_channels": 27, "unet_widths": [16, 32, 64],
              "complex_dtype": "float32_real_pair" if target_device.type == "mps" else "complex64", "network_dtype": "float32",
              "selection_metric": "val_material_l1_plus_0.1_edge",
              "selection_uses_regularizer": False,
              "feature_stats_seconds": stats_seconds, "operator_setup_seconds": operator_setup_seconds,
              "supervised_cpu_cache": {"train": train_dataset.cache_receipt, "val": val_dataset.cache_receipt},
              "V48_source": operators.v48_source, "truth_access": "TRAIN_VAL_ONLY",
              "split_exposure": {"TRAIN": {"material_targets": True, "optimizer": True, "observations": len(train)},
                                 "VAL": {"material_targets": True, "optimizer": False, "selection": True, "observations": len(val)},
                                 "TEST": {"truth": False, "selection": False}},
              "parameter_count": sum(p.numel() for p in model.parameters()),
              "dev_only": dev_only,
              "residual_head_initialization": "all_zero_weights_and_bias",
              "offline_metric_definitions": {"material_relative_l2": "||m_proposal-m_GT||_2/||m_GT||_2",
                                              "material_amplitude_nrmse": "sqrt(mean((m_proposal-m_GT)^2))/GT_amplitude"},
              "arithmetic_estimate": architecture_cost(arm, cascade_passes, operators.projector_metadata["W_rank"], alpha),
              "full_vie_verification": "NOT_RUN_EXTERNAL_REQUIRED"}
    _json_write(output_dir/"config.json", config)
    best, stale, history = float("inf"), 0, []
    fit_start = time.perf_counter()
    for epoch in range(1, epochs+1):
        epoch_start = time.perf_counter()
        tm = _run_epoch(model, operators, train_loader, stats, cascade_passes, target_device, optimizer, alpha, loss_name, regularizer_weight)
        vm = _run_epoch(model, operators, val_loader, stats, cascade_passes, target_device, alpha=alpha, loss_name=loss_name, regularizer_weight=regularizer_weight)
        _sync(target_device)
        row = {"epoch": epoch, "train": tm, "val": vm, "seconds": time.perf_counter()-epoch_start}
        history.append(row)
        selection = vm["material_l1_normalized"]+.1*vm["edge_l1_normalized"]
        row["selection_metric"] = selection
        if selection < best:
            best, stale = selection, 0
            checkpoint = {"format": "gn_nn_r3_v1", "r4_format": "gn_nn_r4_preservation_v1", "state_dict": model.state_dict(),
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
    parser.add_argument("--cascade-passes", type=int, choices=(2,), default=2)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--dev-only", action="store_true")
    parser.add_argument("--alpha", type=float, default=.3)
    parser.add_argument("--loss", choices=("material", "const", "internal"), default="material")
    args = parser.parse_args()
    result = train_model(args.data_root, args.manifest, args.operators, args.output, args.arm,
                         args.seed, args.epochs, args.batch_size, args.patience, args.device,
                         args.cascade_passes, args.anchors, args.num_workers, args.dev_only,
                         alpha=args.alpha, loss_name=args.loss)
    print(json.dumps(result), flush=True)


class _FullVIEFunction(torch.autograd.Function):
    """Exact retained-LU real-material VJP; CPU complex128 physics, no FD."""
    @staticmethod
    def forward(ctx, material, bridge):
        from scipy.linalg import lu_factor, lu_solve
        started = time.perf_counter()
        transfer = time.perf_counter()
        m = material.detach().cpu().numpy().astype(np.float64, copy=True).reshape(-1, NC)
        input_transfer_seconds = time.perf_counter()-transfer
        ys, retained = [], []
        solve_start = time.perf_counter()
        for mi in m:
            a = np.eye(NC, dtype=np.complex128)-mi[:, None]*bridge.G
            lu, piv = lu_factor(a, check_finite=False)
            j = lu_solve((lu, piv), mi[:, None]*bridge.E, check_finite=False)
            q = bridge.E+bridge.G@j
            y = bridge.S@j
            if not np.isfinite(y).all():
                raise FloatingPointError("Nonfinite exact full-VIE prediction")
            ys.append(y)
            retained.append((lu, piv, q))
        physics_seconds = time.perf_counter()-solve_start
        y = np.stack(ys)
        pair = np.stack((y.real, y.imag), axis=-1).astype(np.float32)
        transfer = time.perf_counter()
        result = torch.from_numpy(pair).to(material.device)
        _sync(material.device)
        output_transfer_seconds = time.perf_counter()-transfer
        ctx.bridge, ctx.retained = bridge, retained
        ctx.material_shape, ctx.material_device, ctx.material_dtype = material.shape, material.device, material.dtype
        bridge._record("forward", len(m), len(m), len(m)*NS, 0,
                       time.perf_counter()-started, input_transfer_seconds+output_transfer_seconds,
                       physics_seconds)
        return result

    @staticmethod
    def backward(ctx, grad_y_pair):
        from scipy.linalg import lu_solve
        started = time.perf_counter()
        transfer = time.perf_counter()
        gy = grad_y_pair.detach().cpu().numpy().astype(np.float64, copy=False)
        gy = gy[..., 0]+1j*gy[..., 1]
        input_transfer_seconds = time.perf_counter()-transfer
        gradients = []
        solve_start = time.perf_counter()
        for gi, (lu, piv, q) in zip(gy, ctx.retained):
            rhs = ctx.bridge.S.conj().T@gi
            adjoint = lu_solve((lu, piv), rhs, trans=2, check_finite=False)
            # grad_y_pair already carries both real-output derivatives: no extra2.
            gradients.append(np.real(np.sum(adjoint.conj()*q, axis=1)))
        physics_seconds = time.perf_counter()-solve_start
        gradient = np.stack(gradients)
        if not np.isfinite(gradient).all():
            raise FloatingPointError("Nonfinite exact full-VIE VJP")
        transfer = time.perf_counter()
        result = torch.from_numpy(gradient).to(device=ctx.material_device, dtype=ctx.material_dtype).reshape(ctx.material_shape)
        _sync(ctx.material_device)
        output_transfer_seconds = time.perf_counter()-transfer
        ctx.bridge._record("backward", len(gradient), 0, 0, len(gradient)*NS,
                          time.perf_counter()-started, input_transfer_seconds+output_transfer_seconds,
                          physics_seconds)
        return result, None


class FullVIEBridge(torch.nn.Module):
    """Optional full-Maxwell differentiable backend, never invoked by base training.

    Input (B,576) real material. Output (B,48,12,2) real/imag receiver pair.
    Each forward stores one complex128 CPU LU per case for its exact adjoint.
    The explicit arm/parent/epoch gate stays with the parent experiment runner.
    """
    def __init__(self, operators_path):
        super().__init__()
        started = time.perf_counter()
        with np.load(operators_path, allow_pickle=False) as p:
            self.G = np.asarray(p["G"], np.complex128).copy()
            self.E = np.asarray(p["E0"] if "E0" in p else p["E"], np.complex128).copy()
            self.S = np.asarray(p["S"], np.complex128).copy()
        if self.G.shape != (NC, NC) or self.E.shape != (NC, NS) or self.S.shape != (48, NC):
            raise ValueError("FullVIEBridge requires frozen R4 scalar-TM operators")
        self.setup_seconds = time.perf_counter()-started
        self.operators_path = str(operators_path)
        self.calls = []
        self.counts = {"LU": 0, "source_complex_RHS": 0, "adjoint_complex_RHS": 0,
                       "actual_complex_RHS": 0, "forward_cases": 0, "backward_cases": 0}

    def _record(self, stage, cases, lu, source_rhs, adjoint_rhs, wall, transfers, physics):
        row = {"stage": stage, "cases": cases, "LU": lu,
               "source_complex_RHS": source_rhs, "adjoint_complex_RHS": adjoint_rhs,
               "actual_complex_RHS": source_rhs+adjoint_rhs, "wall_seconds": wall,
               "transfer_wall_seconds_including_input_producer_wait": transfers,
               "cpu_physics_wall_seconds": physics}
        self.calls.append(row)
        self.counts["LU"] += lu
        self.counts["source_complex_RHS"] += source_rhs
        self.counts["adjoint_complex_RHS"] += adjoint_rhs
        self.counts["actual_complex_RHS"] += source_rhs+adjoint_rhs
        self.counts[stage+"_cases"] += cases

    def forward(self, material):
        if material.is_complex() or material.shape[-1] != NC:
            raise ValueError("Expected real material shape(B,576)")
        return _FullVIEFunction.apply(material, self)

    def receipt(self):
        return {"backend": "scipy_complex128_retained_LU_exact_real_material_VJP",
                "operators": self.operators_path, "setup_seconds": self.setup_seconds,
                "counts": dict(self.counts), "calls": list(self.calls),
                "truth_access": "NONE", "finite_difference_physics": False,
                "runtime_status": "NOT_RUN" if not self.calls else "EXECUTED",
                "total_wall_seconds": sum(r["wall_seconds"] for r in self.calls)}


def full_maxwell_hinge(material, public_y_pair, bridge, tolerance):
    """Squared excess relative public discrepancy; tolerance supplied by protocol.

    No target current, support, material truth, or online GT-dependent choice.
    Parent must authorize the separate fine-tuning stage before calling this.
    """
    predicted = bridge(material)
    if public_y_pair.is_complex():
        raise ValueError("Pass public observations as (...,Re/Im) pairs")
    y = public_y_pair.to(device=material.device, dtype=predicted.dtype)
    if y.shape != predicted.shape:
        raise ValueError("Public Y must have shape(B,48,12,2)")
    numerator = torch.linalg.vector_norm((predicted-y).flatten(1), dim=1)
    denominator = torch.linalg.vector_norm(y.flatten(1), dim=1).clamp_min(1e-20)
    relative = numerator/denominator
    tol = torch.as_tensor(tolerance, dtype=relative.dtype, device=relative.device)
    return torch.relu(relative-tol).square().mean(), {"relative_public_residual": relative}


class PhysicsSupervisedAnchors(SupervisedAnchors):
    """Optional-stage TRAIN/VAL dataset adding only public Y/noise to base data."""
    def _read(self, index):
        item = super()._read(index)
        record = self.records[index]
        with np.load(self.anchor_root/f"{record['obs_id']}.npz", allow_pickle=False) as anchor:
            y = np.asarray(anchor["Y"], np.complex64)
            pair = np.stack((y.real, y.imag), axis=-1).astype(np.float32)
            noise = float(np.asarray(anchor["noise"]).reshape(-1)[0])
        item["public_y_pair"] = torch.from_numpy(pair.copy())
        item["public_noise"] = torch.tensor(noise, dtype=torch.float32)
        return item


def _run_finetune_epoch(model, operators, loader, stats, config, device, optimizer,
                        bridge, contract):
    model.train(True)
    sums = {key: 0. for key in ("loss", "existing_loss", "material_l1_normalized",
            "edge_l1_normalized", "public_regularizer", "full_data_hinge",
            "relative_public_residual")}
    count, updates = 0, 0
    for batch in loader:
        batch = {key: to_device(value, device) for key, value in batch.items()}
        optimizer.zero_grad(set_to_none=True)
        output = proposal(model, operators, batch["x"], batch["m0"], batch["j0"], stats,
                          config["cascade_passes"], config["alpha"])
        mloss, terms = material_loss(output["m_proposal"], batch["m_gt"])
        rloss = public_regularizer(output, batch["j0"], batch["q0_norm2"], operators, config["loss"])
        existing = mloss+config["regularizer_weight"]*rloss
        hinge = torch.zeros((), device=device)
        residual = torch.zeros((), device=device)
        if bridge is not None:
            tolerance = torch.where(batch["public_noise"] > 0,
                                    float(contract["noisy_delta"]), float(contract["clean_delta"]))
            hinge, physical = full_maxwell_hinge(output["m_proposal"], batch["public_y_pair"],
                                                bridge, tolerance)
            residual = physical["relative_public_residual"].mean()
        loss = existing+float(contract["lambda_full_data_hinge"])*hinge
        if not torch.isfinite(loss):
            raise FloatingPointError("Nonfinite optional fine-tuning loss")
        loss.backward()
        optimizer.step()
        n = len(batch["m_gt"])
        values = {"loss": loss, "existing_loss": existing, **terms,
                  "public_regularizer": rloss, "full_data_hinge": hinge,
                  "relative_public_residual": residual}
        for key, value in values.items():
            sums[key] += float(value.detach())*n
        count += n
        updates += 1
    return {**{key: value/count for key, value in sums.items()}, "optimizer_steps": updates,
            "public_residual_computed": bridge is not None}


def finetune_model(checkpoint_path, data_root, manifest_path, operators_path, output_dir,
                   anchor_root=None, contract_path=None, device="mps", *, authorized=False,
                   authorization=None):
    """Run the frozen, matched optional-stage pair only after parent authorization.

    Both arms start from exactly the same R4 best checkpoint and use one saved
    seed7 batch schedule. Both use pure material/edge supervision in this stage;
    the treatment adds the exact full-VIE hinge. All ten
    epochs run; full VAL material loss alone selects each best checkpoint.
    """
    if not authorized:
        raise PermissionError("Optional fine-tuning requires the parent's VAL/time gate authorization")
    checkpoint_path, data_root, output_dir = Path(checkpoint_path), Path(data_root), Path(output_dir)
    anchor_root = Path(anchor_root) if anchor_root else data_root.parent/"anchors"
    contract_path = Path(contract_path) if contract_path else data_root.parent/"audit/OPTIONAL_FT_CONTRACT.json"
    contract = json.loads(contract_path.read_text())
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint.get("r4_format") != "gn_nn_r4_preservation_v1":
        raise ValueError("Optional fine-tuning requires an R4 checkpoint; historical R3 weights are excluded")
    base_config, stats = checkpoint["config"], checkpoint["feature_stats"]
    records = read_manifest(manifest_path)
    parents = set(contract["parent_ids"])
    train = [r for r in records if r["split"].lower() == "train" and r["parent_id"] in parents]
    val = [r for r in records if r["split"].lower() in ("val", "validation")]
    if len(train) != contract["observations"] or len(parents) != 128 or len(val) != 156:
        raise ValueError("Optional stage requires the frozen128 TRAIN mothers/256 observations and full156 VAL")
    target_device = choose_device(device)
    setup_start = time.perf_counter()
    operators = OperatorCache(operators_path).to(target_device)
    train_dataset = PhysicsSupervisedAnchors(train, data_root, anchor_root, stats, base_config["arm"])
    val_dataset = SupervisedAnchors(val, data_root, anchor_root, stats, base_config["arm"])
    generator = torch.Generator().manual_seed(7)
    epochs, batch_size = int(contract["epochs"]), int(contract["batch_size"])
    permutations = [torch.randperm(len(train), generator=generator).tolist() for _ in range(epochs)]
    output_dir.mkdir(parents=True, exist_ok=True)
    schedule_path = output_dir/"batch_schedule.json"
    schedule = {"seed": 7, "parent_ids": contract["parent_ids"],
                "observation_ids": [r["obs_id"] for r in train],
                "epoch_index_permutations": permutations, "batch_size": batch_size,
                "epochs": epochs, "drop_last": False,
                "optimizer_steps_per_arm": sum((len(p)+batch_size-1)//batch_size for p in permutations)}
    _json_write(schedule_path, schedule)
    setup = {"status": "AUTHORIZED_MATCHED_PAIR", "base_checkpoint": str(checkpoint_path.resolve()),
             "contract": str(contract_path.resolve()), "authorization": authorization or {"explicit_parent_authorization": True},
             "shared_setup_seconds": time.perf_counter()-setup_start, "operator_uploads": 1,
             "feature_statistics": "FROZEN_FROM_R4_BASE_CHECKPOINT_TRAIN_ONLY", "device": str(target_device),
             "train_cache": train_dataset.cache_receipt, "val_cache": val_dataset.cache_receipt,
             "truth_access": "FROZEN_TRAIN_SUBSET_AND_FULL_VAL_ONLY", "test_truth_access": "NONE"}
    _json_write(output_dir/"setup.json", setup)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0)
    summaries = []
    pair_start = time.perf_counter()
    for mode in ("material_loss_control", "full_vie_hinge"):
        random.seed(7)
        np.random.seed(7)
        torch.manual_seed(7)
        model = ResidualUNet(base_config["arm"]).to(target_device)
        model.load_state_dict(checkpoint["state_dict"])
        optimizer = torch.optim.AdamW(model.parameters(), lr=float(contract["learning_rate"]),
                                      weight_decay=float(base_config.get("weight_decay", .01)))
        run_dir = output_dir/mode
        run_dir.mkdir(parents=True, exist_ok=True)
        config = json.loads(json.dumps(base_config))
        config.update(seed=7, lr=float(contract["learning_rate"]), max_epochs=epochs, batch_size=batch_size,
                      patience=None, train_from_scratch=False, initial_weights_source=str(checkpoint_path.resolve()),
                      residual_head_initialization="same_R4_best_checkpoint_for_both_optional_arms",
                      device=str(target_device), requested_device=device, train_count=len(train), val_count=len(val))
        config.update(loss='material',regularizer_weight=0.)
        config["fine_tune"] = {"mode": mode, "base_checkpoint": str(checkpoint_path.resolve()),
                               "base_checkpoint_best_epoch": checkpoint["best_epoch"],
                               "base_seed": base_config["seed"], "num_parents": len(parents),
                               "parent_ids": contract["parent_ids"], "observations": len(train),
                               "contract_path": str(contract_path.resolve()), "batch_schedule": str(schedule_path.resolve()),
                               "epochs": epochs, "batch_size": batch_size,
                               "planned_optimizer_steps": schedule["optimizer_steps_per_arm"],
                               "lambda_full_data_hinge": float(contract["lambda_full_data_hinge"]) if mode == "full_vie_hinge" else 0.,
                               "clean_delta": float(contract["clean_delta"]), "noisy_delta": float(contract["noisy_delta"]),
                               "source_checkpoint_loss": base_config["loss"],
                               "existing_loss": "material",
                               "existing_regularizer_weight": 0.,
                               "authorization": setup["authorization"]}
        config["split_exposure"] = {"TRAIN": {"material_targets": True, "optimizer": True,
                                                "parents": len(parents), "observations": len(train)},
                                    "VAL": {"material_targets": True, "optimizer": False,
                                            "selection": True, "observations": len(val)},
                                    "TEST": {"truth": False, "selection": False}}
        _json_write(run_dir/"config.json", config)
        bridge = FullVIEBridge(operators_path) if mode == "full_vie_hinge" else None
        history, best, updates = [], float("inf"), 0
        fit_start = time.perf_counter()
        try:
            for epoch, order in enumerate(permutations, start=1):
                epoch_start = time.perf_counter()
                batches = [order[i:i+batch_size] for i in range(0, len(order), batch_size)]
                loader = DataLoader(train_dataset, batch_sampler=batches, num_workers=0)
                tm = _run_finetune_epoch(model, operators, loader, stats, config, target_device,
                                         optimizer, bridge, contract)
                updates += tm["optimizer_steps"]
                vm = _run_epoch(model, operators, val_loader, stats, config["cascade_passes"], target_device,
                                alpha=config["alpha"], loss_name=config["loss"],
                                regularizer_weight=config["regularizer_weight"])
                _sync(target_device)
                selection = vm["material_l1_normalized"]+.1*vm["edge_l1_normalized"]
                row = {"epoch": epoch, "train": tm, "val": vm, "selection_metric": selection,
                       "optimizer_steps_cumulative": updates, "seconds": time.perf_counter()-epoch_start}
                history.append(row)
                if selection < best:
                    best = selection
                    saved = {"format": "gn_nn_r3_v1", "r4_format": "gn_nn_r4_preservation_v1",
                             "fine_tune_format": "gn_nn_r4_optional_matched_ft_v1",
                             "state_dict": model.state_dict(), "feature_stats": stats, "config": config,
                             "best_epoch": epoch, "best_val": vm,
                             "parameter_count": sum(p.numel() for p in model.parameters())}
                    temporary = run_dir/"best.pt.tmp"
                    torch.save(saved, temporary)
                    temporary.replace(run_dir/"best.pt")
                _json_write(run_dir/"history.json", history)
                if bridge is not None:
                    _json_write(run_dir/"bridge_cost.json", bridge.receipt())
                print(json.dumps({"optional_fine_tune_mode": mode, **row}), flush=True)
        finally:
            costs = bridge.receipt() if bridge is not None else {
                "runtime_status": "NOT_USED_CONTROL_ARM", "counts": {"LU": 0,
                "source_complex_RHS": 0, "adjoint_complex_RHS": 0, "actual_complex_RHS": 0},
                "calls": [], "setup_seconds": 0., "total_wall_seconds": 0.}
            _json_write(run_dir/"bridge_cost.json", costs)
            _json_write(run_dir/"run_progress.json", {"epochs_completed": len(history),
                        "optimizer_steps": updates, "training_wall_seconds": time.perf_counter()-fit_start,
                        "planned_optimizer_steps": schedule["optimizer_steps_per_arm"]})
        summary = {"status": "TRAINED", "mode": mode, "checkpoint": str(run_dir/"best.pt"),
                   "epochs_run": len(history), "optimizer_steps": updates, "best_material_loss": best,
                   "training_wall_seconds": time.perf_counter()-fit_start, "bridge_cost": costs,
                   "config": config}
        _json_write(run_dir/"summary.json", summary)
        summaries.append(summary)
        _json_write(output_dir/"matched_pair_summary.json", {"status": "RUNNING" if len(summaries) < 2 else "TRAINED",
                    "base_checkpoint": str(checkpoint_path.resolve()), "runs": summaries,
                    "pair_wall_seconds": time.perf_counter()-pair_start, "test_truth_access": "NONE"})
    return {"status": "TRAINED", "runs": summaries, "pair_wall_seconds": time.perf_counter()-pair_start,
            "shared_setup": setup}


if __name__ == "__main__":
    main()
