"""Compare R4 frozen Voxel/J-full proposals and cold voxel Full-J on the same new data."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from weak_r4 import (
    DEFAULT_ARCHIVE, HERE, R4, archived_material_score, new_sources,
    optimize, ratio, read_npz, rms,
)

CHECKPOINTS = {
    "r4_voxel": "v_physics_seed7_cascade2",
    "r4_j_full": "j_full_seed7_cascade2",
}


def physically_score(Model, ops, m, obs, noise, truth):
    tick = perf_counter()
    model = Model(ops)
    residual = model.forward(m)["y"] - obs
    return dict(**archived_material_score(m, truth),
                relative_residual=ratio(residual, np.linalg.norm(obs)),
                noise_multiple=ratio(residual, noise),
                seconds=perf_counter()-tick, counts=dict(model.counts))


def prepare_proposals(r4, cases, need_replay):
    proposals = {}
    provenance = {}
    if need_replay:
        import torch
        from research.research_gn_nn_preservation_r4.src.networks import OperatorCache
        from research.research_gn_nn_preservation_r4.src.predict import predict_anchors
        torch.set_num_threads(1)
        operator_path = r4/"data/operators.npz"
        operator_cache = OperatorCache(operator_path).to(torch.device("cpu"))
    for label, old_name in CHECKPOINTS.items():
        block = {}
        prediction_dir = r4 / "predictions" / old_name
        anchor_paths = [r4/"anchors"/f'{c["obs_id"]}.npz' for c in cases]
        stored = all((prediction_dir/p.name).exists() for p in anchor_paths)
        if not stored and not need_replay:
            raise FileNotFoundError(f"No archived proposals at {prediction_dir}. Run with --replay-checkpoints.")
        checkpoint = r4/"models"/old_name/"best.pt"
        with TemporaryDirectory(prefix="weak_r4_replay_") as tmp:
            tmp = Path(tmp)
            if need_replay:
                replay_receipts = predict_anchors(checkpoint, anchor_paths, operator_path, tmp,
                                                   "cpu", operators=operator_cache)
            replay_diffs = []
            for p in anchor_paths:
                receipt_file = (prediction_dir if stored else tmp) / p.with_suffix(".json").name
                receipt = json.loads(receipt_file.read_text())
                prediction_file = (prediction_dir if stored else tmp) / p.name
                m = read_npz(prediction_file)["m_proposal"]
                block[p.stem] = dict(
                    material=m, strict_base_cost=receipt["base_anchor_cost"],
                    inference_cost=receipt["preparation_cost"],
                )
                if need_replay and stored:
                    replay = read_npz(tmp/p.name)["m_proposal"]
                    replay_diffs.append(dict(obs_id=p.stem, material_rms=rms(m-replay),
                                             material_max=float(np.abs(m-replay).max())))
            provenance[label] = dict(frozen_arm=old_name, checkpoint=str(checkpoint.relative_to(r4)),
                                     used="archived original predictions" if stored else "CPU checkpoint replay",
                                     replay_comparison=replay_diffs,
                                     replay_inference_seconds=[p["inference_seconds"]
                                                               for p in replay_receipts] if need_replay else [])
        proposals[label] = block
    return proposals, provenance


def summarize(report, output):
    groups = defaultdict(list)
    for c in report["cases"]:
        for k, val in c["original_only"].items():
            groups["old_only", "archived", k].append(val)
        for a in c["arms"]:
            for k, val in a["solutions"].items():
                groups[a["library"], a["selector"], k].append(val)
    with output.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(("library","selector","method","n","median_rmse","median_nrmse","median_iou",
                         "median_noise_multiple","qualifies","median_VIE_LU","median_complex_RHS","median_seconds"))
        for (g,s,k),vals in groups.items():
            med = lambda key: np.median([v[key] for v in vals])
            costs = [v["counts"] for v in vals]
            writer.writerow([g,s,k,len(vals),*[med(key) for key in ("rmse","nrmse","iou","noise_multiple")],
                             sum(v["noise_multiple"]<=1.10 for v in vals),
                             np.median([v["LU"] for v in costs]),
                             np.median([v["actual_complex_RHS"] for v in costs]), med("seconds")])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--budget", type=int, choices=(1,4), default=4)
    ap.add_argument("--replay-checkpoints", action="store_true")
    args = ap.parse_args()
    sys.path.insert(0, str(args.archive.resolve()))
    from research.research_z_alg_r4_diversity.src.physics import Model, Regularizer, material_basis
    r4 = args.archive/R4
    result_dir = HERE/"results"
    source_file = result_dir/("weak_r4" if args.budget==4 else "weak_r4_one")
    report = json.loads(source_file.with_suffix(".json").read_text())
    from_study = read_npz(source_file.with_suffix(".npz"))
    ops = read_npz(r4/"data/operators.npz")
    rows = json.loads((r4/"data/manifests/evaluation.json").read_text())["records"]
    records = {r["parent_id"]: r for r in rows if r["split"]=="test" and r["noise"]==.001}
    reg = Regularizer((24,24),smooth=1e-5,tv=1e-5,epsilon=.01)
    far_xy, far_E0 = new_sources(ops, report["config"]["same_radius"], .5)
    near_xy, near_E0 = new_sources(ops, report["config"]["near_radius"], .5)
    all_E0 = np.column_stack((far_E0,near_E0,ops["E0"]))
    arrays = {}
    saved_inputs = {k:v for k,v in ops.items() if k in
                    ("G","S","E0","k","radius","sources","receivers","centers","grid_shape")}
    saved_inputs.update(supplemental_E0=all_E0,
                        supplemental_sources=np.vstack((far_xy,near_xy,ops["sources"])))
    out = result_dir / ("weak_baselines" if args.budget == 4 else "weak_baselines_one")
    pub = dict(
        base_inputs=source_file.name, new_exposures=args.budget,
        correction="same regularized full-J and 1.10 noise discrepancy as R12 branch warm refit, 14 steps",
        cold="same 8,12,24-grid ladder and discrete bounds as R12 ordinary GN, with full analytic J and zero material regularizer",
        cases=[],
    )
    with threadpool_limits(1):
        proposals, pub["prediction_provenance"] = prepare_proposals(
            r4, report["cases"], args.replay_checkpoints)
        for c in report["cases"]:
            pid, obs_id = c["parent_id"], c["obs_id"]
            row = records[pid]
            truth = read_npz(r4/"data"/row["truth_path"])
            old_Y = from_study[f"{pid}_old_Y"]
            full_new_Y = from_study[f"{pid}_new_Y_all"]
            eps = float(row["noise"] * np.linalg.norm(old_Y))
            per_source = eps/np.sqrt(ops["E0"].shape[1])
            case = dict(parent_id=pid, obs_id=obs_id, original_only={}, arms=[])
            pub["cases"].append(case)
            saved_inputs.update({f"{pid}_truth":truth["m"], f"{pid}_alpha":truth["alpha"],
                                 f"{pid}_amplitude":truth["amplitude"], f"{pid}_old_Y":old_Y,
                                 f"{pid}_new_Y_all":full_new_Y})
            ladder = tuple((material_basis((24,24),k),8) for k in (8,12))+((None,14),)
            cold = np.full(576,.1)
            for tag, inputs in proposals.items():
                m = inputs[obs_id]["material"]
                arrays[f"{pid}_{tag}_raw"] = m
                case["original_only"][tag+"_raw"] = physically_score(Model,ops,m,old_Y,eps,truth)
                case["original_only"][tag+"_raw"].update(
                    archived_strict_GN_cost=inputs[obs_id]["strict_base_cost"],
                    archived_network_cost=inputs[obs_id]["inference_cost"])
                corrected = optimize(Model,reg,ops,old_Y,eps,m,((None,14),),True,tag+"_corrected")
                estimate = corrected.pop("material")
                arrays[f"{pid}_{tag}_old_corrected"] = estimate
                case["original_only"][tag+"_corrected"] = dict(**corrected,**archived_material_score(estimate,truth))
            original_voxel = optimize(Model,reg,ops,old_Y,eps,cold,ladder,False,"cold_voxel_fullJ")
            original_voxel_m = original_voxel.pop("material")
            arrays[f"{pid}_old_cold_voxel"] = original_voxel_m
            case["original_only"]["cold_voxel_fullJ"] = dict(
                **original_voxel,**archived_material_score(original_voxel_m,truth))
            visited = {}
            for arm in c["arms"]:
                ids = tuple(map(int, arm["source_ids"]))
                prefix = f'{pid}_{arm["library"]}_{arm["selector"]}'
                if ids not in visited:
                    arm_ops = dict(ops,E0=np.column_stack((ops["E0"],all_E0[:,ids])))
                    Y = np.column_stack((old_Y,full_new_Y[:,ids]))
                    joint_noise = float(np.hypot(eps,per_source*np.sqrt(len(ids))))
                    solutions = {}
                    fits = {"cold_voxel_fullJ": (cold,ladder,False)}
                    for tag, inputs in proposals.items():
                        m = inputs[obs_id]["material"]
                        solutions[tag+"_raw"] = physically_score(Model,arm_ops,m,Y,joint_noise,truth)
                        fits[tag+"_corrected"] = (m,((None,14),),True)
                    fitted = {}
                    for tag,(m,stages,smooth) in fits.items():
                        fit = optimize(Model,reg,arm_ops,Y,joint_noise,m,stages,smooth,tag)
                        estimate = fit.pop("material")
                        fitted[tag] = estimate
                        solutions[tag] = dict(**fit,**archived_material_score(estimate,truth))
                    visited[ids] = solutions,fitted
                solutions,fitted = visited[ids]
                case["arms"].append(dict(library=arm["library"],selector=arm["selector"],
                                         source_ids=ids,solutions=solutions))
                for tag,m in fitted.items():
                    arrays[f"{prefix}_{tag}"] = m
            out.with_suffix(".json").write_text(json.dumps(pub,indent=2,default=float)+"\n")
            np.savez_compressed(out.with_suffix(".npz"),**arrays)
            print(pid,[(a["library"],a["selector"],
                        *(round(a["solutions"][key]["rmse"],5) for key in
                          ("r4_voxel_corrected","r4_j_full_corrected","cold_voxel_fullJ")))
                       for a in case["arms"] if a["selector"] in ("candidate_gap","equally_spaced")],flush=True)
    summarize(pub,out.with_suffix(".csv"))
    if args.budget == 4:
        np.savez_compressed(result_dir/"weak_inputs.npz",**saved_inputs)


if __name__ == "__main__":
    main()
