"""Use actual A13 Phase B weights on their own existing native observations."""

import argparse
import json
import sys

import numpy as np
from scipy.cluster.hierarchy import linkage, cut_tree
from threadpoolctl import threadpool_limits

from pool_common import RESULTS, FLOW, FullWave, clock, elapsed, write_arrays, write_json

SPLIT = FLOW / "results/phase_b/data/main/test"
RUN = FLOW / "results/phase_b/runs/main/seed_11/material"
OUTPUT = RESULTS / "flow_results"
IDS = (0, 37, 111)
SEED = 11077
COUNT = 64
CLUSTERS = 12
TOP_K = (1, 3, 6)
CAP = 5.


def generation():
    """Generation/refinement phase has no file path to material labels or held-out ports."""
    import torch
    sys.path.insert(0, str(FLOW / "src"))
    import physics_learning as learned

    if (OUTPUT / "online.json").exists():
        raise FileExistsError("Online predictions already exist.")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)
    started = clock()
    spec = json.loads((SPLIT / "manifest.json").read_text())
    solver = learned.make_solver()
    if not (solver.n == spec["n"] and solver.k == spec["k"]
            and solver.n_rx == spec["n_rx"] and solver.n_illum == spec["n_illum"]
            and solver.rx_radius == spec["rx_radius"]):
        raise ValueError("Incompatible learned checkpoint and acquisition geometry.")
    levels = learned.nested_rx_indices(16)
    mask = learned.nested_rx_mask(16)[None]
    flow, decoder = learned.build_model(), learned.build_model()
    for model, stage in ((flow, "flow"), (decoder, "probe")):
        state = torch.load(RUN / stage / "checkpoint.pt", map_location="cpu", weights_only=False)
        model.load_state_dict(state["model"])
        model.eval()
    stats = learned.load_latent_stats(RUN / "latent_stats.npz")
    flow_metrics = json.loads((RUN / "flow" / "done.json").read_text())
    if flow_metrics["noise"] != 0.0 or 16 not in flow_metrics["rx_levels"]:
        raise ArithmeticError("The declared acquisition disagrees with training.")
    setup = elapsed(started)
    # Only training-receiver values leave the observation file for the online function.
    observed = np.load(SPLIT / "y.npy", mmap_mode="r")
    receipts = {}
    for index in IDS:
        y = np.asarray(observed[index, levels, :], dtype=np.complex128)
        object_id = f"a13_main_test_{spec['seed'] + index:06d}__rx16"
        receipts[object_id] = run_case(object_id, index, y, solver, levels, mask, flow,
                                       decoder, stats, learned, torch)

    write_json(OUTPUT / "online.json", dict(
        cases=receipts, setup=setup, weights_seed=11, sample_count=COUNT, cluster_count=CLUSTERS,
        radius="1% of measured training norm; clean archived acquisition, nominal tolerance to stored-field rounding",
        support="Full 32² free, signed contrast; a [-5,+5] policy is declared; clipping counts reported",
        training_mask_levels=flow_metrics["rx_levels"], training_noise=flow_metrics["noise"],
        label_access="online uses only y.npy projected to 16 training receivers",
        source_paths={name: str((RUN / name / "checkpoint.pt").relative_to(FLOW)) for name in ("flow", "probe")},
        physical_certificate="No rigorous full-grid shifted or Neumann proof run for these 1024-pixel signed boxes.",
        elapsed=elapsed(started)))


def run_case(object_id, index, y, solver, used_rx, mask, flow, decoder, stats, learned, torch):
    begin = clock()
    y_input = np.zeros((1, solver.n_rx, solver.n_illum), dtype=np.complex64)
    y_input[0, used_rx] = y
    tokens = learned.measurement_feature_tensor(torch.as_tensor(y_input), torch.as_tensor(mask), 0.,
                                                 solver.receiver_pos, solver.illum_dir)
    with torch.no_grad():
        latent = learned.sample_latent(flow, tokens, samples=COUNT, seed=SEED + index)[0]
        latent = latent * stats["std"] + stats["mean"]
        raw = learned.image_to_real(decoder.decode(latent, learned.MODE_INDEX["MATERIAL"]).cpu().numpy())
    candidates = np.clip(raw.astype(np.float64), -CAP, CAP)
    labels = cut_tree(linkage(candidates, method="ward"), n_clusters=CLUSTERS).ravel()
    lowers = np.stack([np.maximum(-CAP, candidates[labels == i].min(0) - .10) for i in range(CLUSTERS)])
    uppers = np.stack([np.minimum(CAP, candidates[labels == i].max(0) + .10) for i in range(CLUSTERS)])
    sample_time = elapsed(begin)

    model = FullWave(solver.GD, solver.GS[used_rx], solver.incident, y)
    begin_point = clock()
    point = np.stack([model.evaluate(m)[0] for m in candidates])
    point_cost = elapsed(begin_point)
    point_residuals = np.linalg.norm((point - y).reshape(COUNT, -1), axis=1)
    costs, materials, data = [], [], []
    ranked = np.argsort(point_residuals, kind="stable")
    representatives = [next(i for i in ranked if labels[i] == c) for c in range(CLUSTERS)]
    for c, rep in enumerate(representatives):
        m, reading, cost = model.refine(candidates[rep], lowers[c], uppers[c], max_calls=9)
        materials.append(m)
        data.append(reading)
        costs.append(cost)

    radius = .01 * np.linalg.norm(y)
    order = list(map(int, np.argsort(point_residuals[representatives], kind="stable")))
    pipelines = {}
    for name, visited in [("all_cluster", order)] + [(f"top_{k}", order[:k]) for k in TOP_K]:
        winner = min(visited, key=lambda c: costs[c]["residual"])
        pipelines[name] = dict(winner=winner, branches=len(visited), visit=visited,
                               logical_full_wave_lu=COUNT + sum(costs[c]["counts"]["lu"] for c in visited),
                               online_wall_seconds=sample_time["wall_seconds"] + point_cost["wall_seconds"]
                               + sum(costs[c]["wall_seconds"] for c in visited),
                               online_cpu_seconds=sample_time["cpu_seconds"] + point_cost["cpu_seconds"]
                               + sum(costs[c]["cpu_seconds"] for c in visited),
                               feasible_endpoint=bool(costs[winner]["residual"] <= radius))
    write_arrays(OUTPUT / f"{object_id}.npz", train=y, used_rx=used_rx, raw_samples=raw,
                 samples=candidates, sample_predictions=point, cluster=labels, radius=radius,
                 lows=lowers, highs=uppers, refined=np.array(materials), train_prediction=np.array(data))
    print(f"{object_id} 64 samples; best point/r {min(point_residuals)/radius:.3f}; "
          f"best 12-box refinement/r {min(c['residual'] for c in costs)/radius:.3f}", flush=True)
    return dict(index=index, radius=radius, point_residuals=point_residuals.tolist(), cluster_sizes=np.bincount(labels).tolist(),
                clipped_entries=int(np.count_nonzero(raw != candidates)), setup=sample_time, point_cost=point_cost,
                counts=dict(model.counts), representative_indices=list(map(int, representatives)),
                branch_cost=costs, pipelines=pipelines)


def offline():
    sys.path.insert(0, str(FLOW / "src"))
    import physics

    results = json.loads((OUTPUT / "online.json").read_text())
    solver = physics.Solver()
    y_all = np.load(SPLIT / "y.npy", mmap_mode="r")
    truth_all = np.load(SPLIT / "chi.npy", mmap_mode="r")
    scores, port_arrays = {}, {}
    for object_id, row in results["cases"].items():
        i = row["index"]
        clean = y_all[i]
        truth = truth_all[i].astype(float)
        with np.load(OUTPUT / f"{object_id}.npz") as saved:
            test_ix = np.setdiff1d(np.arange(solver.n_rx), saved["used_rx"])
            model = FullWave(solver.GD, solver.GS, solver.incident, clean)
            physical_truth = model.evaluate(truth)[0]
            background = np.abs(truth) < .1 * np.max(np.abs(truth))
            ref, predictions = [], []
            for material in saved["refined"]:
                predicted = model.evaluate(material)[0]
                predictions.append(predicted)
                ref.append(dict(
                    overall=float(np.linalg.norm(material - truth) / np.linalg.norm(truth)),
                    inside=float(np.linalg.norm((material-truth)[~background]) / np.linalg.norm(truth[~background])),
                    background_rms=float(np.sqrt(np.mean((material-truth)[background]**2))),
                    holdout_receiver=float(np.linalg.norm((predicted-clean)[test_ix]) / np.linalg.norm(clean[test_ix]))))
            port_arrays[f"{object_id}__truth"] = truth
            port_arrays[f"{object_id}__heldout_rx"] = test_ix
            port_arrays[f"{object_id}__clean_ports"] = clean
            port_arrays[f"{object_id}__refined_full_ports"] = np.stack(predictions)
            scores[object_id] = dict(
                strict_box_contains_truth=np.flatnonzero(np.all((saved["lows"] <= truth)
                                                                 & (truth <= saved["highs"]), axis=1)).tolist(),
                forward_versus_saved=float(np.linalg.norm(physical_truth - clean) / np.linalg.norm(clean)),
                truth_range=[float(truth.min()), float(truth.max())],
                branch_scores=ref,
                pipelines={key: ref[value["winner"]] for key, value in row["pipelines"].items()})
    write_arrays(OUTPUT / "offline_geometry_scores.npz", **port_arrays)
    write_json(OUTPUT / "offline.json", scores)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("online", "offline"))
    with threadpool_limits(limits=1):
        (generation if parser.parse_args().phase == "online" else offline)()
