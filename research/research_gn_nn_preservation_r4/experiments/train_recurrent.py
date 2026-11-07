"""Prepared recurrent-feedback ablations; no training starts without --authorized."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import random
import sys
import time
import numpy as np
import torch
from torch.utils.data import DataLoader

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY/"src"))
from networks import (NC, NS, GRID, MATERIAL_SCALE, OperatorCache, choose_device,
                      fit_feature_stats, material_loss, public_regularizer)
from recurrent_feedback import (MODES, RecurrentFeedback, feedback_step,
                                HINGE_CLEAN_DELTA, HINGE_NOISY_DELTA)
from predict_recurrent import predict_recurrent_anchors
from train import (SupervisedAnchors, load_anchor, released_training_cohort,
                   to_device, _sync, _json_write)


class RecurrentAnchors(SupervisedAnchors):
    def _read(self, index):
        item = super()._read(index)
        obs_id = self.records[index]["obs_id"]
        with np.load(self.anchor_root/f"{obs_id}.npz", allow_pickle=False) as p:
            item["public_y"] = torch.from_numpy(np.asarray(p["Y"], np.complex64).copy())
            item["anchor_q"] = torch.from_numpy(np.asarray(p["Q0"], np.complex64).copy())
            item["public_noise"] = torch.tensor(float(np.asarray(p["noise"]).reshape(-1)[0]), dtype=torch.float32)
        return item


def _epoch(model, operators, loader, stats, device, loss_name, regularizer_weight, optimizer=None):
    model.train(optimizer is not None)
    sums = {k: 0. for k in ("loss", "material_l1_normalized", "edge_l1_normalized", "public_regularizer")}
    count, updates, counts = 0, 0, {}
    first_gradient = last_gradient = None
    started = time.perf_counter()
    with torch.set_grad_enabled(optimizer is not None):
        for batch in loader:
            batch = {key: to_device(value, device) for key, value in batch.items()}
            if optimizer is not None:
                optimizer.zero_grad(set_to_none=True)
            out = model(operators, batch["x"], batch["m0"], batch["j0"], batch["public_y"], stats, batch["anchor_q"],
                        public_noise=batch["public_noise"])
            mloss, terms = material_loss(out["m_proposal"], batch["m_gt"])
            rloss = public_regularizer(out, batch["j0"], batch["q0_norm2"], operators, loss_name)
            loss = mloss+regularizer_weight*rloss
            if not torch.isfinite(loss):
                raise FloatingPointError("Nonfinite recurrent loss")
            if optimizer is not None:
                loss.backward()
                grad = {"head_l2": float(torch.sqrt(sum(p.grad.square().sum() for p in model.network.head.parameters())).detach()),
                        "stem_l2": float(torch.sqrt(sum(p.grad.square().sum() for p in model.network.stem.parameters())).detach())}
                first_gradient = first_gradient or grad
                last_gradient = grad
                optimizer.step()
                updates += 1
            n = len(batch["m_gt"])
            for key, value in {"loss": loss, **terms, "public_regularizer": rloss}.items():
                sums[key] += float(value.detach())*n
            for key, value in out["preparation_counts"].items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    # Operator call counts are batched; column/FLOP fields include batch size.
                    if key not in ("shared_parameter_count", "feedback_steps_per_layer", "feedback_steps_total"):
                        counts[key] = counts.get(key, 0)+value
            count += n
    _sync(device)
    return {**{key: value/count for key, value in sums.items()}, "optimizer_steps": updates,
            "observations": count, "wall_seconds": time.perf_counter()-started,
            "forward_operation_totals": counts, "backward_operation_counts": "NOT_ESTIMATED",
            "gradient_first_step": first_gradient, "gradient_last_step": last_gradient}


def train_recurrent(data_root, manifest_path, operators_path, anchor_root, output_dir,
                    *, rho=1., modes=("refreshed_q", "feedback1"), seed=7, device="mps", stats_checkpoint=None,
                    loss_name="material", authorized=False, predict_val=False):
    if not authorized:
        raise PermissionError("Parent must include this optional VAL ablation before training")
    data_root, anchor_root, output_dir = Path(data_root), Path(anchor_root), Path(output_dir)
    cohort = released_training_cohort(data_root, manifest_path, anchor_root)
    train, val = cohort["train"], cohort["val"]
    stats_source = "full_TRAIN_public_anchors_only"
    if stats_checkpoint is None:
        stats = fit_feature_stats(load_anchor(anchor_root/f"{r['obs_id']}.npz") for r in train)
    else:
        c = torch.load(stats_checkpoint, map_location="cpu", weights_only=False)
        if c.get("r4_format") != "gn_nn_r4_preservation_v1":
            raise ValueError("Statistics checkpoint must be from the new R4 study")
        stats, stats_source = c["feature_stats"], str(Path(stats_checkpoint).resolve())
        # Only normalization is read; network state is never loaded.
    target = choose_device(device)
    setup_start = time.perf_counter()
    operators = OperatorCache(operators_path).to(target)
    train_set = RecurrentAnchors(train, data_root, anchor_root, stats, "j_full")
    val_set = RecurrentAnchors(val, data_root, anchor_root, stats, "j_full")
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_receipt = {"operator_uploads": 1, "device": str(target),
                 "wall_seconds": time.perf_counter()-setup_start, "feature_stats_source": stats_source,
                 "normalization_refitted_per_layer": False, "old_network_weights_loaded": False,
                 "training_status": "EXPLICIT_PARENT_AUTHORIZED", "test_truth_access": "NONE",
                 "invocation_unix_ns": time.time_ns()}
    _json_write(output_dir/f"shared_setup_{setup_receipt['invocation_unix_ns']}.json", setup_receipt)
    _json_write(output_dir/"shared_setup.json", setup_receipt)
    summaries = []
    for mode in modes:
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        model = RecurrentFeedback(mode, rho).to(target)
        run_dir = output_dir/f"recurrent_{mode}_seed{seed}_cascade2"
        if (run_dir/"best.pt").exists():
            raise FileExistsError(run_dir/"best.pt")
        run_dir.mkdir(parents=True, exist_ok=True)
        generator = torch.Generator().manual_seed(seed)
        train_loader = DataLoader(train_set, batch_size=64, shuffle=True, generator=generator, num_workers=0)
        val_loader = DataLoader(val_set, batch_size=64, shuffle=False, num_workers=0)
        optimizer = torch.optim.AdamW(model.parameters(), lr=.001)
        config = {"arm": f"recurrent_{mode}", "network_arm": "j_full", "mode": mode,
                  "rho": float(rho), "rho_fixed": True, "cascade_passes": 2,
                  "device": str(target), "requested_device": device,
                  "feedback_data_objective": "noise_ball_squared_hinge" if mode == "feedback_hinge" else "quadratic",
                  "feedback_step_rule": "hinge_majorizer" if mode == "feedback_hinge" else "quadratic_exact_line",
                  "feedback_data_tolerance": {"clean": HINGE_CLEAN_DELTA, "noisy": HINGE_NOISY_DELTA,
                                              "source": "explicit_public_anchor_noise"} if mode == "feedback_hinge" else 0.,
                  "noise_added_to_NN_channels": False,
                  "feedback_steps_per_layer": model.feedback_steps, "seed": seed,
                  "max_epochs": 100, "batch_size": 64, "patience": 15, "lr": .001,
                  "optimizer": "AdamW", "weight_decay": .01, "loss": loss_name,
                  "regularizer_weight": .01, "train_from_scratch": True,
                  "initial_weights_source": "NONE", "old_checkpoint_weights_loaded": False,
                  "input_channels": 27, "unet_widths": [16,32,64],
                  "parameter_count": sum(p.numel() for p in model.parameters()),
                  "feature_stats_source": stats_source, "feature_stats_policy": "same_TRAIN_statistics_all_layers",
                  "Q_policy": "fixed_anchor" if mode == "no_feedback" else "refreshed_from_current",
                  "gb_noise_policy": "fixed_anchor", "physics_inverse_calls_per_layer": 0,
                  "normalization_denominators": "fixed_public_Y_and_anchor_j0_norm_squared",
                  "selection_metric": "VAL_material_L1_plus_0.1_edge",
                  "complex_dtype": "float32_real_pair" if target.type == "mps" else "complex64",
                  "split_exposure": {"TRAIN": {"optimizer": True, "material_targets": True, "observations": len(train)},
                                     "VAL": {"optimizer": False, "material_targets": True, "selection": True, "observations": len(val)},
                                     "TEST": {"truth": False, "selection": False}},
                  "full_VIE_verification": "NOT_RUN_EXTERNAL_REQUIRED"}
        _json_write(run_dir/"config.json", config)
        history, best, stale = [], float("inf"), 0
        started = time.perf_counter()
        for epoch in range(1, 101):
            tm = _epoch(model, operators, train_loader, stats, target, loss_name, .01, optimizer)
            vm = _epoch(model, operators, val_loader, stats, target, loss_name, .01)
            selection = vm["material_l1_normalized"]+.1*vm["edge_l1_normalized"]
            row = {"epoch": epoch, "train": tm, "val": vm, "selection_metric": selection}
            history.append(row)
            if selection < best:
                best, stale = selection, 0
                saved = {"format": "gn_nn_r4_recurrent_feedback_v1", "state_dict": model.state_dict(),
                         "feature_stats": stats, "config": config, "best_epoch": epoch, "best_val": vm}
                temporary = run_dir/"best.pt.tmp"
                torch.save(saved, temporary)
                temporary.replace(run_dir/"best.pt")
            else:
                stale += 1
            _json_write(run_dir/"history.json", history)
            print(json.dumps({"recurrent_mode": mode, **row}), flush=True)
            if stale >= 15:
                break
        summary = {"status": "TRAINED", "mode": mode, "checkpoint": str(run_dir/"best.pt"),
                   "epochs_run": len(history), "best_material_loss": best,
                   "training_wall_seconds": time.perf_counter()-started, "config": config}
        _json_write(run_dir/"summary.json", summary)
        if predict_val:
            # This list is the complete156-VAL cohort already established at launch.
            prediction_dir = STUDY/"predictions"/run_dir.name
            receipts = predict_recurrent_anchors(run_dir/"best.pt",
                [anchor_root/f"{r['obs_id']}.npz" for r in val], operators_path, prediction_dir,
                device, operators=operators, stats=stats)
            summary["public_val_predictions"] = {"count": len(receipts), "output": str(prediction_dir),
                                                   "truth_access": "NONE"}
            _json_write(run_dir/"summary.json", summary)
            print(json.dumps({"status": "RECURRENT_TRAINED_VAL_PREDICTED", "mode": mode,
                              "checkpoint": str(run_dir/"best.pt"), "predictions": str(prediction_dir)}), flush=True)
        summaries.append(summary)
        _json_write(output_dir/"recurrent_runs.json", summaries)
    return summaries


def check_wiring(operators_path, device="mps"):
    """One deterministic synthetic wiring check; no cohort/GT or accuracy-floor test."""
    target = choose_device(device)
    operators = OperatorCache(operators_path).to(target)
    torch.manual_seed(7)
    e = operators._complex("E") if target.type != "mps" else torch.stack((operators.E_re, operators.E_im), dim=-1)
    m0 = torch.full((2,NC), .12, device=target)
    j0 = m0[...,None]*e if e.is_complex() else m0[...,None,None]*e
    q0 = operators.total_field(j0)
    y = operators._matmul("S", j0)*1.01
    stats = {"mean": [0.]*27, "std": [1.]*27, "current_scale": .04}
    features = torch.zeros(2,27,GRID,GRID,device=target)
    features[:,0] = m0.reshape(2,GRID,GRID)/MATERIAL_SCALE
    target_m = m0+.01*torch.sin(torch.arange(NC,device=target)/20.)[None]
    model = RecurrentFeedback("feedback2",rho=1.).to(target)
    optimizer = torch.optim.AdamW(model.parameters(),lr=.001)
    gradients=[]
    for step in range(2):
        optimizer.zero_grad(set_to_none=True)
        out=model(operators,features,m0,j0,y,stats,q0)
        loss,_=material_loss(out["m_proposal"],target_m)
        loss.backward()
        head=float(torch.sqrt(sum(p.grad.square().sum() for p in model.network.head.parameters())).detach())
        stem=float(torch.sqrt(sum(p.grad.square().sum() for p in model.network.stem.parameters())).detach())
        if not np.isfinite(head+stem) or head<=0 or (step==1 and stem<=0):
            raise RuntimeError("Recurrent gradient connection failed")
        gradients.append({"step":step,"head_l2":head,"stem_l2":stem,"loss":float(loss.detach())})
        optimizer.step()
    zero_j=torch.zeros_like(j0[:1])
    zero_y=torch.zeros_like(y[:1])
    zero_m=torch.zeros_like(m0[:1])
    zero=feedback_step(operators,zero_j,zero_m,zero_y,zero_j,1.)
    if not bool(torch.isfinite(zero["current"]).all()) or not bool((zero["tau"]==0).all()):
        raise RuntimeError("Zero-gradient sample produced invalid feedback")
    return {"status":"PASS","device":str(target),"source":"deterministic_synthetic_only",
            "truth_or_cohort_access":"NONE","retained_model_weights":False,
            "check_scope":"two-layer/two-feedback-step autograd plus zero-gradient branch",
            "gradient_steps":gradients,"stage_materials_shape":list(out["stage_materials"].shape),
            "zero_gradient_tau":zero["tau"].detach().cpu().tolist(),
            "preparation_counts":out["preparation_counts"]}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-root",type=Path,default=STUDY/"data")
    p.add_argument("--manifest",type=Path,default=STUDY/"data/manifests/evaluation.json")
    p.add_argument("--operators",type=Path,default=STUDY/"data/operators.npz")
    p.add_argument("--anchors",type=Path,default=STUDY/"anchors")
    p.add_argument("--output",type=Path,default=STUDY/"models/recurrent")
    p.add_argument("--modes",choices=MODES,nargs="+",default=["refreshed_q","feedback1"])
    p.add_argument("--rho",type=float,choices=(1.,),default=1.,help="Frozen rho=1 for the authorized finite VAL ablation")
    p.add_argument("--seed",type=int,choices=(7,17,27),default=7)
    p.add_argument("--device",default="mps")
    p.add_argument("--stats-checkpoint",type=Path,help="Reuse new R4 TRAIN normalization only, never its weights")
    p.add_argument("--loss",choices=("material","const","internal"),default="material")
    p.add_argument("--authorized",action="store_true")
    p.add_argument("--predict-val",action="store_true",help="Write complete156 public VAL candidates after each trained model")
    p.add_argument("--check-wiring",action="store_true")
    args=p.parse_args()
    if args.check_wiring:
        print(json.dumps(check_wiring(args.operators,args.device)),flush=True)
        return
    if not args.authorized:
        p.error("Training awaits parent inclusion/authorization within the existing6h budget")
    result=train_recurrent(args.data_root,args.manifest,args.operators,args.anchors,args.output,
                           rho=args.rho,modes=args.modes,seed=args.seed,device=args.device,
                           stats_checkpoint=args.stats_checkpoint,loss_name=args.loss,authorized=True,predict_val=args.predict_val)
    print(json.dumps({"status":"TRAINED","runs":len(result)}),flush=True)


if __name__=="__main__":
    main()
