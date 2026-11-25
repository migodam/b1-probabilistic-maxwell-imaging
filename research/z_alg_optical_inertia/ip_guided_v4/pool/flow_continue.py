"""Additional A13 refinements on previously frozen native samples and boxes."""

import argparse
import json
import sys

import numpy as np
from threadpoolctl import threadpool_limits

from pool_common import FLOW, FullWave, clock, elapsed, write_arrays, write_json
from flow_audit import OUTPUT


def experiment():
    started = clock()
    source = json.loads((OUTPUT / "online.json").read_text())
    sys.path.insert(0, str(FLOW / "src"))
    import physics

    begin_setup = clock()
    solver = physics.Solver()
    setup_cost = elapsed(begin_setup)
    receipts, evidence = {}, {}
    for key, row in source["cases"].items():
        radius = row["radius"]
        with np.load(OUTPUT / f"{key}.npz", allow_pickle=False) as data:
            y, rxs, tests = data["train"], data["used_rx"], data["refined"]
            lower, upper = data["lows"], data["highs"]
        model = FullWave(solver.GD, solver.GS[rxs], solver.incident, y)
        initial = np.argsort([r["residual"] for r in row["branch_cost"]], kind="stable")
        receipt = {"radius": radius, "refined_box_results": {}}
        begin = clock()
        for c in initial:
            fitted, port, rec = model.refine(tests[c], lower[c], upper[c], max_calls=30)
            receipt["refined_box_results"][str(c)] = rec
            evidence[f"{key}__box{c:02d}__material"] = fitted
            evidence[f"{key}__box{c:02d}__port"] = port
        receipt["all_box_continuation"] = elapsed(begin)
        receipt["box_continuation_counts"] = dict(model.counts)
        broad, port, unrestricted = model.refine(tests[initial[0]], np.full(1024, -5.),
                                                 np.full(1024, 5.), max_calls=60)
        evidence[f"{key}__wide_box__material"] = broad
        evidence[f"{key}__wide_box__port"] = port
        receipt["same_support_wide_box"] = unrestricted
        receipt["original_best_branch"] = int(initial[0])
        receipts[key] = receipt
        print(f"{key} best 12 native boxes/r = "
              f"{min(x['residual'] for x in receipt['refined_box_results'].values())/radius:.4f}, "
              f"wide/r = {unrestricted['residual']/radius:.4f}", flush=True)
    write_arrays(OUTPUT / "continuation_arrays.npz", **evidence)
    write_json(OUTPUT / "continuation.json", dict(
        current_pass=receipts, setup=setup_cost, totals=elapsed(started),
        scope="Existing A13 proposals, post-primary audit, historical observations already exposed. "
              "30 more calls for every box, 60 calls for standalone signed [-5,5] full-domain control."))


def offline():
    sys.path.insert(0, str(FLOW / "src"))
    import physics

    results = json.loads((OUTPUT / "continuation.json").read_text())["current_pass"]
    solver = physics.Solver()
    report = {}
    with np.load(OUTPUT / "offline_geometry_scores.npz") as private, \
            np.load(OUTPUT / "continuation_arrays.npz") as data:
        for key, row in results.items():
            clean, truth, heldout = (private[f"{key}__clean_ports"], private[f"{key}__truth"],
                                     private[f"{key}__heldout_rx"])
            branch_id = min(row["refined_box_results"], key=lambda c: row["refined_box_results"][c]["residual"])
            model = FullWave(solver.GD, solver.GS, solver.incident, clean)
            background = np.abs(truth) < .1 * np.max(np.abs(truth))
            by_method = {}
            for method, approach in [("continued_boxes", f"box{int(branch_id):02d}"), ("wide_box", "wide_box")]:
                material = data[f"{key}__{approach}__material"]
                predicted = model.evaluate(material)[0]
                by_method[method] = dict(
                    material_relative_l2=float(np.linalg.norm(material-truth)/np.linalg.norm(truth)),
                    inside_relative_l2=float(np.linalg.norm((material-truth)[~background])/np.linalg.norm(truth[~background])),
                    background_rms=float(np.sqrt(np.mean((material-truth)[background]**2))),
                    heldout_rx=float(np.linalg.norm((predicted-clean)[heldout])/np.linalg.norm(clean[heldout])))
            report[key] = dict(best_box=int(branch_id), scores=by_method)
    write_json(OUTPUT / "continuation_offline.json", report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true")
    with threadpool_limits(limits=1):
        (offline if parser.parse_args().offline else experiment)()
