"""Measure separate end-to-end pipeline passes and replay a lower-endpoint ablation."""

import argparse
import json
import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from pool_common import FullWave, clock, elapsed, write_arrays, write_json
from box_certificates import Bounds, ball, gershgorin, lower, star
from structured_pool import BRANCHES, MAX_CALLS, N_STARTS, OUT, PUBLIC, RX, TX, filenames, model_inputs, proposals

METHODS = ("all_cluster", "top_best_1", "top_best_3", "top_best_6", "top_early",
           "neumann_standard", "neumann_directional", "shifted")
ORDER_SEED = 33937


def run_one(method, g, s, b, y, radius, seed):
    total = clock()
    model = FullWave(g, s, b, y)
    samples, lo, hi, prep = proposals(g, s, b, y, model, seed)
    starts = samples[::N_STARTS].copy()
    begin_screen = clock()
    bounds = Bounds(g, s, b, y)
    branch_order = list(range(len(BRANCHES)))
    branch_choice = []
    screen_rows = []
    if method.startswith("top"):
        point = np.stack([model.evaluate(m)[0] for m in samples])
        dist = np.linalg.norm((point-y).reshape(len(samples), -1), axis=1).reshape(-1, N_STARTS)
        best = np.argmin(dist, axis=1)
        branch_order = np.argsort(dist[np.arange(len(BRANCHES)), best], kind="stable").tolist()
        selected = N_STARTS*np.arange(len(BRANCHES))+best
        starts = samples[selected]
        if method.startswith("top_best"):
            branch_order = branch_order[:int(method.rsplit("_", maxsplit=1)[1])]
    elif method != "all_cluster":
        kind = "strict_lower" if method == "shifted" else method.split("_")[1]+"_lower"
        verify = bounds.shifted if method == "shifted" else bounds.neumann
        screen_rows = [verify(l, h, radius) for l, h in zip(lo, hi, strict=True)]
        branch_order = [i for i, row in enumerate(screen_rows) if row[kind] <= radius]
    screen_time = elapsed(begin_screen)
    best_residual, winner, selected_arrays = np.inf, None, None
    fits = []
    checked_upper = None
    endpoint = dict(wall_seconds=0., cpu_seconds=0.)
    for i in branch_order:
        branch_choice.append(i)
        material, field, info = model.refine(starts[i], lo[i], hi[i], MAX_CALLS)
        fits.append(info)
        if info["residual"] < best_residual:
            best_residual, winner, selected_arrays = info["residual"], i, (material, field)
        if method == "top_early" and info["residual"] <= radius:
            begin_endpoint = clock()
            checked_upper = bounds.point_upper(material)
            t = elapsed(begin_endpoint)
            endpoint = {name: t[name]+endpoint[name] for name in endpoint}
            if checked_upper <= radius:
                break
    if checked_upper is None and winner is not None:
        begin_endpoint = clock()
        checked_upper = bounds.point_upper(selected_arrays[0])
        endpoint = elapsed(begin_endpoint)
    return dict(retained=len(branch_order), visited=branch_choice, winner=winner,
                selector="first strict feasible" if method == "top_early" else "minimum observed residual",
                residual=best_residual if winner is not None else None,
                verified_residual_upper=checked_upper, preparation=prep,
                screening=screen_time, refinement_fits=fits, endpoint_audit=endpoint,
                dense_optimization_counts=dict(model.counts),
                neumann_numeric_center_factorizations=sum(r.get("center_factorizations", 0) for r in screen_rows),
                neumann_numeric_center_transfer_rhs=sum(r.get("center_transfer_rhs", 0) for r in screen_rows),
                neumann_interval_factorizations=sum(r.get("verified_inverse_factors", 0) for r in screen_rows),
                shifted_proposals=sum(r.get("proposals", 0) for r in screen_rows),
                shifted_interval_successes=sum(r.get("proof") is not None for r in screen_rows),
                **elapsed(total)), selected_arrays


def exact_surjectivity(s, b):
    """Nonsingular row Gram, with its SVD preconditioner rechecked in Arb."""
    result = {}
    combined = np.column_stack((s.T, b)).T
    for axis, matrix in [("receiver", s), ("source", b.T), ("reciprocal", combined)]:
        m = ball(matrix)
        u, sigma, _ = np.linalg.svd(matrix, full_matrices=False)
        preconditioned = star(ball(u / sigma)) * m
        gap = gershgorin(preconditioned * star(preconditioned), upper=False)
        result[axis] = dict(strict_positive_gram=bool(gap > 0),
                            preconditioned_interval_gap_lower=lower(gap))
    return result


def support_only():
    ctx.prec, ctx.threads = 112, 1
    start = clock()
    _, sensor, incident = model_inputs()
    result = exact_surjectivity(sensor[RX], incident[:, TX])
    write_json(OUT / "support_only.json", dict(port_range=result, **elapsed(start)))


def main(timing_only=False):
    if (OUT / "timed_pipelines.json").exists():
        raise FileExistsError(OUT / "timed_pipelines.json")
    clock_all = clock()
    offline_frozen = (OUT / "offline.json").is_file()
    ctx.prec, ctx.threads = 112, 1
    start_setup = clock()
    full_g, full_s, full_b = model_inputs()
    g, s, b = full_g, full_s[RX], full_b[:, TX]
    surjective = exact_surjectivity(s, b)
    setup = elapsed(start_setup)
    receipts, saved = {}, {}
    rng = np.random.default_rng(ORDER_SEED)
    for seed, case in filenames():
        with np.load(PUBLIC / f"{case}.npz", allow_pickle=False) as data:
            y, radius = data["train"], float(data["radius"])
        methods = [str(v) for v in rng.permutation(METHODS)]
        method_receipts = {}
        for method in methods:
            receipt, arrays = run_one(method, g, s, b, y, radius, seed+1007)
            method_receipts[method] = receipt
            if arrays is not None:
                saved[f"{case}__{method}__material"], saved[f"{case}__{method}__train"] = arrays
        ablation = {}
        if not timing_only:
            with np.load(OUT / f"{case}.npz", allow_pickle=False) as source:
                lows, highs = source["lower"], source["upper"]
            bounds = Bounds(g, s, b, y)
            row = {}
            for name, lo, hi in zip(BRANCHES, lows, highs, strict=True):
                cap_only = bounds.shifted(np.zeros_like(lo), hi, radius)
                row[name] = cap_only
            ablation = dict(remove_all_positive_lower_endpoints=row,
                            excluded=sum(v["strict_lower"] > radius for v in row.values()))
        print(case, ", ".join(f"{k}:{method_receipts[k]['wall_seconds']:.3f}s" for k in METHODS),
              ("cap-only="+str(ablation.get("excluded","n/a"))), flush=True)
        receipts[case] = dict(radius=radius, methods=method_receipts, execution_order=methods,
                              peak_only_ablation=ablation)
    write_arrays(OUT / "timed_pipeline_arrays.npz", **saved)
    write_json(OUT / "timed_pipelines.json", dict(geometry_setup=setup, port_range=surjective, cases=receipts,
               retrospective_same_objects=offline_frozen, total=elapsed(clock_all)))


def offline():
    """Score endpoints solely after the independent timed pipeline pass is written."""
    timed = json.loads((OUT / "timed_pipelines.json").read_text())
    g, s, b = model_inputs()
    scores, predictions = {}, {}
    rx_hold = np.setdiff1d(np.arange(s.shape[0]), RX)
    tx_hold = np.setdiff1d(np.arange(b.shape[1]), TX)
    from structured_pool import PRIVATE, SIDE, residual

    with np.load(OUT / "timed_pipeline_arrays.npz", allow_pickle=False) as arrays:
        for _, case in filenames():
            with np.load(PRIVATE / f"{case}.npz", allow_pickle=False) as private:
                truth = private["material"]
                ports = private["noiseless_ports"]
            model = FullWave(g, s, b, ports)
            scores[case] = {}
            for method, row in timed["cases"][case]["methods"].items():
                if row["winner"] is None:
                    scores[case][method] = None
                    continue
                estimated = arrays[f"{case}__{method}__material"]
                y_all = model.evaluate(estimated)[0]
                predictions[f"{case}__{method}"] = y_all
                relative_l2 = float(np.linalg.norm(estimated-truth)/np.linalg.norm(truth))
                scores[case][method] = dict(full_image_relative_l2=relative_l2, inside_relative_l2=relative_l2,
                    background_rms=0., background_cells=SIDE**2-len(truth),
                    heldout_receiver_relative_l2=residual(y_all, ports, rx_hold, TX),
                    heldout_source_relative_l2=residual(y_all, ports, RX, tx_hold))
    write_arrays(OUT / "timed_pipeline_offline_ports.npz", **predictions)
    write_json(OUT / "timed_pipelines_offline.json", scores)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timing-only", action="store_true")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--support-only", action="store_true")
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.offline:
            offline()
        elif args.support_only:
            support_only()
        else:
            main(args.timing_only)
