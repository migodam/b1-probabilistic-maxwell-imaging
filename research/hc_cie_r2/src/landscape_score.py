"""Post-seal material scoring for the free-face curvature intervention."""
import argparse
from pathlib import Path
import time

from imaging_run import LAB, json_write, recorded_input_path
from imaging_score import get_case, material_score, read, write_csv, mark_pareto

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np


ARMS = ["continue", "gauss_newton", "random", "curvature"]
LABELS = dict(gauss_newton="Gauss–Newton", random="Random", curvature="Curvature", continue_="Continue")
LABELS["continue"] = LABELS.pop("continue_")


def read_runs(base):
    seal = read(base / "PREDICTION_SEAL.json")
    if not seal["without_truth"] or seal["expected"] != len(seal["jobs"]):
        raise ValueError("Required prediction seal is incomplete.")
    cases = {}
    for job in seal["jobs"]:
        directory = base / job["directory"]
        trace = read(directory / "TRACE.json")
        assert trace["without_truth"]
        rows = trace["chain"]
        with np.load(directory / "PREDICTIONS.npz", allow_pickle=False) as raw:
            arrays = {name: raw[name] for name in raw.files}
        with np.load(recorded_input_path(trace["input_paths"][-2]), allow_pickle=False) as raw:
            measured = raw["Y"]
        chosen = min([row for row in rows if row["paid"]], key=lambda row: row["Phi"])["i"]
        assert chosen == job["selected"] == int(arrays["selected_index"])
        assert np.isfinite(arrays["material"]).all()
        assert arrays["material"].min() >= .02 - 1e-12
        assert arrays["material"].max() <= 20. + 1e-12
        for index, ports in [(chosen, "selected_port"), (0, "endpoint_port")]:
            assert np.isclose(np.linalg.norm(arrays[ports]-measured)/np.linalg.norm(measured),
                              rows[index]["relative_port_error"], atol=1e-14, rtol=1e-12)
        case = cases.setdefault(job["case"], {})
        case[job["arm"]] = (trace, arrays)
    return seal, cases


def score_row(trace, arrays, row, truth, stage):
    return dict(case=trace["case"], endpoint=trace["endpoint"], arm=trace["arm"], source=stage,
                group=truth["group"], object_id=truth["object_id"], termination=trace["termination"],
                total_cpu=trace["total_cpu"], total_wall=trace["total_wall"], direction_cpu=trace["direction_cpu"],
                total_full_LU=trace["cost"].get("LU", 0),
                total_complex_RHS=trace["cost"].get("actual_complex_RHS", 0),
                **{k: v for k, v in row.items() if k not in {"setup_counts", "original_counts"}},
                **material_score(arrays["material"][row["i"]], truth, (24, 24)))


def decision(rows):
    rng = np.random.default_rng(26112843)
    selected = [r for r in rows if r["source"] == "selected"]
    all_end = {(r["case"], r["arm"]): r for r in rows if r["source"] == "selected"}
    reference = {(r["case"], r["arm"]): r for r in rows if r["source"] == "endpoint"}
    summary = []
    groups = sorted({r["group"] for r in selected})
    for arm in ARMS:
        for group in ["all", *groups]:
            tests = [r for r in selected if r["arm"] == arm and (group == "all" or r["group"] == group)]
            if not tests:
                continue
            delta = np.array([r["NRMSE"] - all_end[(r["case"], "continue")]["NRMSE"] for r in tests])
            deltas_bs = delta[rng.integers(0, len(tests), size=(4000, len(tests)))].mean(axis=1)
            summary.append(dict(arm=arm, group=group, objects=len(tests),
                                mean_NRMSE=float(np.mean([r["NRMSE"] for r in tests])),
                                mean_delta_vs_continue=float(np.mean(delta)),
                                descriptive_bootstrap_95=np.quantile(deltas_bs, [.025, .975]).tolist(),
                                more_physically_accurate_than_continue=sum(
                                    r["Phi"] < all_end[(r["case"], "continue")]["Phi"] for r in tests),
                                new_pass=sum(r["image_pass"] and not reference[(r["case"], arm)]["image_pass"] for r in tests),
                                lost_pass=sum(not r["image_pass"] and reference[(r["case"], arm)]["image_pass"] for r in tests)))
    return summary


def gallery(cases, truths, outcomes, figdir):
    ids = list(cases)
    selected = {(r["case"], r["arm"]): r for r in outcomes if r["source"] == "selected"}
    vmax = max(float(np.abs(cases[c]["curvature"][1]["direction"]).max()) for c in ids)
    for page, start in enumerate(range(0, len(ids), 4), 1):
        subset = ids[start:start+4]
        fig, axes = plt.subplots(len(subset), 7, figsize=(14, len(subset)*2.5),
                                 layout="constrained", squeeze=False)
        for i, case in enumerate(subset):
            records = cases[case]
            trace, arrays = records["curvature"]
            values = [truths[case]["material"], arrays["material"][0]]
            labels = ["Truth", case + "\nSaved start"]
            for arm in ARMS:
                t, a = records[arm]
                values.append(a["material"][t["selected"]])
                labels.append(f"{LABELS[arm]}\nNRMSE {selected[(case, arm)]['NRMSE']:.3f}")
            for j, (material, name) in enumerate(zip(values, labels)):
                contrast = axes[i, j].imshow(material.reshape(24,24), origin="lower",
                                             vmin=0, vmax=20, cmap="viridis")
                axes[i, j].set_title(name, fontsize=9)
            modal = axes[i, 6].imshow(arrays["direction"].reshape(24,24), origin="lower",
                                     vmin=-vmax, vmax=vmax, cmap="RdBu_r")
            eigen = trace["diagnostics"]["free_min_eigenvalue"]
            axes[i, 6].set_title(f"Free direction\nλ={eigen:.2g}", fontsize=9)
            for axis in axes[i]:
                axis.set_axis_off()
        fig.colorbar(contrast, ax=axes[:, :6], shrink=.7, label="Material contrast")
        fig.colorbar(modal, ax=axes[:, 6], shrink=.7, label="Direction, unit L2")
        fig.savefig(figdir/f"landscape_directions_page{page}.png", dpi=145)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sealed", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--figures", type=Path, default=LAB / "figures")
    args = parser.parse_args()
    started = time.process_time()
    seal, cases = read_runs(args.sealed)
    truth = {case: get_case(case) for case in cases}
    outcomes = []
    prefixes = []
    diagnostic = []
    for case, variants in cases.items():
        for arm, (trace, arrays) in variants.items():
            chain = trace["chain"]
            for r in chain:
                if r["stage"] in {"endpoint", "jump"}:
                    outcomes.append(score_row(trace, arrays, r, truth[case], r["stage"]))
            outcomes.append(score_row(trace, arrays, chain[trace["selected"]], truth[case], "selected"))
            for cap in trace["configuration"]["prefixes"]:
                row = min([r for r in chain if r["paid"] and r["cpu"] <= cap], key=lambda r: r["Phi"])
                prefixes.append(dict(budget_cpu=cap, **score_row(trace, arrays, row, truth[case], "prefix")))
            diagnostic.append(dict(case=case, endpoint=trace["endpoint"], arm=arm, direction_cpu=trace["direction_cpu"],
                                   all_cpu=trace["total_cpu"], all_wall=trace["total_wall"],
                                   origin_relative_port_disagreement=trace["endpoint_port_disagreement"],
                                   unique_paid_jumps=sum(r["stage"] == "jump" and r["paid"] for r in chain),
                                   **trace["diagnostics"]))
    args.out.mkdir(parents=True, exist_ok=False)
    mark_pareto(prefixes)
    write_csv(args.out/"OUTCOMES.csv", outcomes)
    write_csv(args.out/"PAID_PREFIXES.csv", prefixes)
    write_csv(args.out/"DIRECTIONS.csv", diagnostic)
    json_write(args.out/"PAIRED_SUMMARY.json", decision(outcomes))
    args.figures.mkdir(parents=True, exist_ok=True)
    gallery(cases, truth, outcomes, args.figures)
    json_write(args.out/"SCORING_RECEIPT.json",
               dict(post_seal_scoring=True, n=seal["expected"], process_cpu=time.process_time()-started))


if __name__ == "__main__":
    main()
