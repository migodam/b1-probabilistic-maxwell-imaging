"""Score sealed predictions. All truth loading and outcome figures live here."""
import argparse
import csv
import json
from pathlib import Path
import time

from imaging_run import LAB, REPO, R13, R6, R6_TREE, json_write, recorded_input_path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
from scipy import ndimage
from scipy.stats import rankdata


FAILURES = {"r13_fresh_0007", "r13_fresh_0016", "r13_fresh_0029",
            "r13_fresh_0037", "r13_fresh_0044"}
ARM_LABEL = dict(continue_="Continue", restart="Restart", single="Single",
                 pair_random="Random pair", pair_born="Born pair",
                 pair_secant="Secant pair", pair_exact="Exact pair")
ARM_LABEL["continue"] = ARM_LABEL.pop("continue_")


def read(path):
    return json.loads(path.read_text())


def get_case(case):
    if case.startswith("r13_fresh_"):
        with np.load(R13 / "sealed_truth/fresh48" / f"{case}.npz", allow_pickle=False) as z:
            material = z["material_coarse"].ravel()
        support = material > 1.
        amplitude = float(np.max(material) - .1)
        object_id = case
        group = "R13_common_failure" if case in FAILURES else "R13_positive_control"
        rule = "contrast > 1"
        data_rule = ("CN99", 1.)
    else:
        cases = read(R6 / "data_public/cases.json")["records"]
        scoring = R6 / "sealed_truth/scoring"
        mapping = scoring / "offline_replay_mapping.json"
        if not mapping.is_file():
            mapping = scoring / "historical_mapping_PRIVATE.json"
        origins = read(mapping)
        public = next(r for r in cases if r["case_id"] == case)
        truth = next(r for r in origins if r["case_id"] == case)
        with np.load(R6_TREE / truth["truth_path"], allow_pickle=False) as z:
            material = z["m"].ravel()
            support = z["mask"].ravel().astype(bool)
        amplitude = float(public["amplitude"])
        object_id = truth.get("object_id", case)
        group = "R6_historical_curved_bar"
        rule = "Otsu256"
        if case == "obs_b20d5a55677d4eca":
            with np.load(R6_TREE / public["observation_file"], allow_pickle=False) as z:
                sigma = float(z["sigma"])
                threshold = float(z["noise_radius"])
            data_rule = ("legacy_whitened", sigma * threshold)
        else:
            scores = R6_TREE / "research/research_fbg_r3_experiments/results/e3_all_scores.csv"
            with scores.open() as f:
                threshold = next(float(r["physical_data_gate"]) for r in csv.DictReader(f)
                                 if r["case_id"] == case and r["physical_data_gate"])
            data_rule = ("legacy_relative", threshold)
    return dict(material=material, support=support, amplitude=amplitude, rule=rule,
                object_id=object_id, group=group, exposed=True, data_rule=data_rule)


def score_joint_shape(segmentation, shape):
    binary = segmentation.reshape(shape)
    components = int(ndimage.label(binary)[1])
    background, n_background = ndimage.label(~binary, structure=np.ones((3, 3), dtype=int))
    boundary = set(np.r_[background[0, :], background[-1, :],
                         background[:, 0], background[:, -1]]) - {0}
    return components, n_background - len(boundary)


def otsu(material):
    if np.ptp(material) <= 1e-12:
        return np.zeros(material.size, dtype=bool)
    hist, edges = np.histogram(material, bins=256)
    center = (edges[:-1] + edges[1:]) / 2
    weight = hist.cumsum()
    total = (hist * center).cumsum()
    den = weight * (weight[-1] - weight)
    score = np.full(256, -np.inf)
    good = den > 0
    score[good] = (weight[good] * total[-1] - total[good] * weight[-1])**2 / den[good]
    return material > center[np.argmax(score)]


def material_score(material, case, shape):
    err = material - case["material"]
    support = case["support"]
    amp = case["amplitude"]
    seg = material > 1. if case["rule"] == "contrast > 1" else otsu(material)
    inter = np.count_nonzero(seg & support)
    union = np.count_nonzero(seg | support)
    nrmse = float(np.sqrt(np.mean(err**2)) / amp)
    iou = float(inter / union) if union else 1.
    c, h = score_joint_shape(seg, shape)
    true_c, true_h = score_joint_shape(support, shape)
    return dict(NRMSE=nrmse, IoU=iou, L2=float(np.linalg.norm(err) / np.linalg.norm(case["material"])),
                RMSE=float(np.sqrt(np.mean(err**2))),
                foreground_NRMSE=float(np.sqrt(np.mean(err[support]**2)) / amp),
                background_NRMSE=float(np.sqrt(np.mean(err[~support]**2)) / amp),
                foreground_bias=float(np.mean(err[support])),
                precision=float(inter / np.count_nonzero(seg)) if seg.any() else 0.,
                recall=float(inter / np.count_nonzero(support)) if support.any() else 0.,
                components=c, holes=h, true_components=true_c, true_holes=true_h,
                matches_structure=c == true_c and h == true_h,
                image_pass=nrmse <= .15 and iou >= .7)


def add_result(recipe, trace_row, material, truth, source):
    res = trace_row["residual_over_delta"]
    score = material_score(material, truth, tuple(recipe["geometry_shape"]))
    data_kind, data_threshold = truth["data_rule"]
    if data_kind == "CN99":
        tested = res
    elif data_kind == "legacy_relative":
        tested = trace_row["relative_residual"]
    else:
        tested = trace_row["relative_residual"] * recipe["measured_port_norm"]
    return dict(variant=recipe["config"]["id"], case=recipe["case"], object=truth["object_id"],
                group=truth["group"], endpoint=recipe["endpoint"], arm=recipe["arm"],
                source=source, index=trace_row["i"], stage=trace_row["stage"], cpu=trace_row["cpu"],
                wall=trace_row["wall"], total_cpu=recipe["total_cpu"],
                total_wall=recipe["total_wall"], stop=recipe["termination"],
                Phi=trace_row["Phi"], D=trace_row["D"], R=trace_row["R"],
                relative_port_error=trace_row["relative_residual"], residual_over_delta=res,
                material_and_CN99=bool(score["image_pass"] and res <= 1.) if res is not None else None,
                historical_data_test=data_kind, historical_data_threshold=data_threshold,
                historical_data_pass=tested <= data_threshold,
                historical_joint_pass=bool(score["image_pass"] and tested <= data_threshold),
                prior_LU=trace_row["cost"].get("LU", 0),
                prior_actual_complex_RHS=trace_row["cost"].get("actual_complex_RHS", 0),
                small_LU=trace_row["proposal_cost"].get("small_LU", 0),
                small_RHS=trace_row["proposal_cost"].get("small_RHS", 0),
                total_LU=recipe["counts"].get("LU", 0),
                total_actual_complex_RHS=recipe["counts"].get("actual_complex_RHS", 0),
                forecast_Phi=trace_row["forecast"], proposal_cpu=recipe.get("proposal_cpu"),
                **score)


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        table = csv.DictWriter(f, fieldnames=keys)
        table.writeheader()
        table.writerows(rows)


def load_sealed(path):
    data = read(path / "PREDICTION_SEAL.json")
    if not data["without_truth"] or data["n"] != len(data["jobs"]):
        raise ValueError(f"Invalid prediction seal in {path}")
    measurements = {}
    for row in data["jobs"]:
        directory = path / row["directory"]
        trace = read(directory / "TRACE.json")
        if trace["status"] != "PREDICTION_COMPLETE" or not trace["without_truth"]:
            raise ValueError(f"Unsealed data {directory}")
        with np.load(directory / "PREDICTIONS.npz", allow_pickle=False) as z:
            materials = z["material"]
            actual_index = int(z["selected_index"])
            selected_port, endpoint_port = z["selected_port"], z["endpoint_port"]
        chain = trace["chain"]
        eligible = [entry for entry in chain if entry["paid"]]
        assert len(chain) == row["candidate_count"] == len(materials)
        assert np.isfinite(materials).all()
        assert materials.min() >= trace["lower"] - 1e-12
        assert materials.max() <= trace["upper"] + 1e-12
        assert min(eligible, key=lambda c: c["Phi"])["i"] == row["selected"] == actual_index
        obs = trace["input_paths"][-2]
        if obs not in measurements:
            with np.load(recorded_input_path(obs), allow_pickle=False) as z:
                measurements[obs] = z["Y"]
        y = measurements[obs]
        trace["measured_port_norm"] = float(np.linalg.norm(y))
        for port, index in [(selected_port, actual_index), (endpoint_port, 0)]:
            residual = np.linalg.norm(port - y) / trace["measured_port_norm"]
            assert np.isclose(residual, chain[index]["relative_residual"], rtol=1e-12, atol=1e-14)
        yield trace, materials


def read_port_reliability(trace):
    chain = trace["chain"]
    jump = [line for line in chain if line["stage"] == "jump"]
    if not jump:
        return None
    actual = np.array([line["Phi"] for line in jump])
    forecast = np.array([line["forecast"] for line in jump])
    if np.ptp(actual) > 0 and np.ptp(forecast) > 0:
        rank = float(np.corrcoef(rankdata(actual), rankdata(forecast))[0, 1])
    else:
        rank = None
    return dict(variant=trace["config"]["id"], case=trace["case"], end=trace["endpoint"],
                arm=trace["arm"], n_audited=len(jump),
                n_with_strict_improvement=int(np.sum(actual < chain[0]["Phi"])),
                endpoint_Phi=chain[0]["Phi"], best_actual_jump_Phi=float(actual.min()),
                best_forecast_Phi=float(forecast.min()), rank_correlation=rank,
                max_Phi_forecast_error=float(np.max(np.abs(actual - forecast))),
                relative_Phi_forecast_L1=float(np.sum(np.abs(actual - forecast)) / np.sum(actual)),
                best_jump_index=int(jump[int(np.argmin(actual))]["i"]))


def cross_validate(rows, folds):
    """Fit only a method selector on public training residuals, hold out whole objects."""
    eligible = sorted({r["arm"] for r in rows})
    cases = {r["case"]: r["object"] for r in rows}
    if len({case for fold in folds for case in fold}) != len(cases):
        raise ValueError("Unexpected CV observation coverage")
    test_counter = set()
    table = []
    detail = []
    lookup = {(r["case"], r["endpoint"], r["arm"]): r for r in rows}
    for i, fold in enumerate(folds, 1):
        test_objects = {cases[t] for t in fold}
        if test_counter & test_objects:
            raise ValueError("Object appears in more than one holdout fold")
        test_counter |= test_objects
        train_objects = set(cases.values()) - test_objects
        baseline = [r for r in rows if r["arm"] == "continue"]
        train_scores = {}
        for arm in eligible:
            by_object = {obj: [] for obj in train_objects}
            for row in baseline:
                obj = row["object"]
                if obj in by_object:
                    trial = lookup[(row["case"], row["endpoint"], arm)]
                    by_object[obj].append(np.log(trial["Phi"] / row["Phi"]))
            train_scores[arm] = float(np.median([np.mean(ds) for ds in by_object.values()]))
        chosen = min(eligible, key=lambda arm: (train_scores[arm], arm != "continue", arm))
        fold_results = []
        for case in fold:
            matches = [r for r in baseline if r["case"] == case]
            improvement = []
            for row in matches:
                trial = lookup[(row["case"], row["endpoint"], chosen)]
                delta_in_image = trial["NRMSE"] - row["NRMSE"]
                improvement.append(delta_in_image)
                detail.append(dict(fold=i, object=row["object"], case=row["case"],
                                   endpoint=row["endpoint"], chosen=chosen, test_Phi=trial["Phi"],
                                   continue_Phi=row["Phi"], test_NRMSE=trial["NRMSE"],
                                   continue_NRMSE=row["NRMSE"], NRMSE_delta=delta_in_image))
            fold_results.append(float(np.mean(improvement)))
        table.append(dict(fold=i, train_materials=sorted(train_objects), held_out_materials=sorted(test_objects),
                          train_selection_statistic="median of object means of log(Phi/Phi_continue)",
                          train_scores=train_scores, chosen=chosen,
                          held_out_object_mean_NRMSE_delta=float(np.mean(fold_results))))
    if test_counter != set(cases.values()):
        raise ValueError("Missing held out object")
    return table, detail


def robustness(rows, comparisons):
    lookup = {(r["case"], r["endpoint"], r["arm"], r["source"]): r for r in rows}
    pairs = []
    for row in rows:
        if row["source"] != "selected" or row["arm"] not in comparisons:
            continue
        for against in ("endpoint", "continue", "restart", "pair_random", "single", "pair_born"):
            source = "endpoint" if against == "endpoint" else "selected"
            baseline_arm = row["arm"] if against == "endpoint" else against
            base = lookup[(row["case"], row["endpoint"], baseline_arm, source)]
            pairs.append(dict(case=row["case"], object=row["object"], group=row["group"],
                              endpoint=row["endpoint"], method=row["arm"], reference=against,
                              NRMSE_diff=row["NRMSE"] - base["NRMSE"], IoU_diff=row["IoU"] - base["IoU"],
                              Phi_ratio=row["Phi"] / base["Phi"], reference_image_pass=base["image_pass"],
                              image_pass=row["image_pass"], delta_cost=row["total_cpu"] - base["total_cpu"]))
    rng = np.random.default_rng(261128)
    summary = []
    groups = sorted({r["group"] for r in pairs})
    for method in comparisons:
        for reference in ("endpoint", "continue", "restart", "pair_random", "single", "pair_born"):
            for group in ["all", *groups]:
                records = [row for row in pairs if row["method"] == method and row["reference"] == reference
                           and (group == "all" or group == row["group"])]
                if not records:
                    continue
                objects = sorted({r["object"] for r in records})
                delta = np.array([np.mean([r["NRMSE_diff"] for r in records if r["object"] == ob]) for ob in objects])
                boots = delta[rng.integers(0, len(delta), size=(4000, len(delta)))].mean(axis=1)
                summary.append(dict(method=method, reference=reference, group=group, objects=len(objects),
                                    endpoint_pairs=len(records), mean_NRMSE_delta=float(delta.mean()),
                                    bootstrap95=[float(x) for x in np.quantile(boots, [.025, .975])],
                                    improved_objects=int(np.count_nonzero(delta < -1e-9)),
                                    rescued_endpoints=sum(r["image_pass"] and not r["reference_image_pass"] for r in records),
                                    lost_endpoints=sum(r["reference_image_pass"] and not r["image_pass"] for r in records)))
    return pairs, summary


def mark_pareto(rows):
    grouped = {}
    for row in rows:
        group = (row["case"], row["endpoint"])
        grouped.setdefault(group, []).append(row)
    for group in grouped.values():
        xyz = np.array([[r["budget_cpu"], r["Phi"], r["NRMSE"]] for r in group])
        for row, costs in zip(group, xyz):
            row["Pareto"] = not bool(np.any(np.all(xyz <= costs, axis=1) & np.any(xyz < costs, axis=1)))


def gallery(cases, selected, pictures, out, variant):
    for case, truth in cases.items():
        arm_order = ["continue", "restart", "single", "pair_random", "pair_born", "pair_secant", "pair_exact"]
        endpoints = list(dict.fromkeys(r["endpoint"] for r in selected if r["case"] == case))
        arms = [a for a in arm_order if any(r["arm"] == a for r in selected)]
        cols = 2 + len(arms)
        fig, axes = plt.subplots(len(endpoints), cols, figsize=(2.15 * cols, 2.25 * len(endpoints)),
                                 squeeze=False, layout="constrained")
        by = {(r["endpoint"], r["arm"], r["source"]): r for r in selected if r["case"] == case}
        for y, endpoint in enumerate(endpoints):
            images = [("Truth", truth["material"], None),
                      (endpoint + " start", pictures[(case, endpoint, arms[0], "endpoint")],
                       by[(endpoint, arms[0], "endpoint")])]
            for arm in arms:
                images.append((ARM_LABEL.get(arm, "Continue"),
                               pictures[(case, endpoint, arm, "selected")], by[(endpoint, arm, "selected")]))
            for x, (label, m, row) in enumerate(images):
                axis = axes[y, x]
                im = axis.imshow(m.reshape(24, 24), origin="lower", cmap="viridis", vmin=0, vmax=20)
                axis.set_title(label if row is None else f"{label}\nNRMSE {row['NRMSE']:.3f} · IoU {row['IoU']:.2f}",
                               fontsize=9)
                axis.set_axis_off()
        fig.suptitle(f"{case} · {variant} · historical exposure", fontsize=12)
        fig.colorbar(im, ax=axes, label="Material contrast", shrink=.68)
        fig.savefig(out / f"imaging_{variant}_{case}.png", dpi=125)
        plt.close(fig)


def trajectories(rows, out, variant):
    cases = sorted({r["case"] for r in rows})
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    arms = sorted({r["arm"] for r in rows})
    for case in cases:
        endpoints = list(dict.fromkeys(r["endpoint"] for r in rows if r["case"] == case))
        fig, axes = plt.subplots(len(endpoints), 2, figsize=(10, 2.3 * len(endpoints)), squeeze=False, layout="constrained")
        for idx, endpoint in enumerate(endpoints):
            for color, arm in zip(colors, arms):
                entries = [r for r in rows if r["arm"] == arm and r["case"] == case and r["endpoint"] == endpoint]
                x = [r["budget_cpu"] for r in entries]
                label = ARM_LABEL.get(arm, "Continue")
                axes[idx, 0].plot(x, [r["Phi"] for r in entries], label=label, color=color, marker=".")
                axes[idx, 1].plot(x, [r["NRMSE"] for r in entries], label=label, color=color, marker=".")
            axes[idx, 0].set_yscale("log")
            axes[idx, 0].set_ylabel(endpoint + " · original Φ")
            axes[idx, 1].set_ylabel("Amplitude NRMSE")
            for ax in axes[idx]:
                ax.grid(alpha=.25)
        for ax in axes[-1]:
            ax.set_xlabel("Paid process CPU (s)")
        axes[0, 0].legend(loc="best", ncol=3, fontsize=7)
        fig.suptitle(case + " · online min-Φ incumbent")
        fig.savefig(out / f"imaging_{variant}_{case}_cost.png", dpi=130)
        plt.close(fig)


def score(frozen, out, figdir, plots):
    receipts = list(load_sealed(frozen))
    out.mkdir(parents=True, exist_ok=False)
    figdir.mkdir(parents=True, exist_ok=True)
    case_data = {name: get_case(name) for name in dict.fromkeys(trace["case"] for trace, _ in receipts)}
    selected = []
    pictures = {}
    prefixes = []
    forecasts = []
    provenance = []
    base_costs = {}
    for trace, materials in receipts:
        case = trace["case"]
        truth = case_data[case]
        chain = trace["chain"]
        reliable = read_port_reliability(trace)
        if reliable:
            forecasts.append(reliable)
        chosen = trace.get("chosen_jump")
        if chosen is None and reliable:
            chosen = reliable["best_jump_index"]
        identifiers = [("endpoint", 0), ("selected", trace["selected"]),
                       ("before_polish", trace["selected_before_polish"])]
        if chosen is not None:
            identifiers.append(("best_raw_jump", chosen))
        for label, position in identifiers:
            picture = materials[position]
            row = add_result(trace, chain[position], picture, truth, label)
            selected.append(row)
            pictures[(case, trace["endpoint"], trace["arm"], label)] = picture
        for cap in trace["config"]["prefix_cpu_seconds"]:
            feasible = [entry for entry in chain if entry["paid"] and entry["cpu"] <= cap]
            if not feasible:
                continue
            entry = min(feasible, key=lambda r: r["Phi"])
            row = add_result(trace, entry, materials[entry["i"]], truth, "prefix")
            prefixes.append(dict(budget_cpu=cap, **row))
        provenance.append(dict(case=case, object=truth["object_id"], observation=trace["input_paths"],
                               endpoint=trace["endpoint"], method=trace["arm"], group=truth["group"],
                               exposure="previously scored public material", IoU_rule=truth["rule"],
                               amplitude_denominator=truth["amplitude"],
                               port_relative_disagreement=trace["endpoint_relative_port_disagreement"]))
        parent_key = (case, trace["endpoint"])
        if parent_key not in base_costs:
            origin = recorded_input_path(trace["input_paths"][-1])
            parent = read(origin.with_name("RECEIPT.json"))
            base_costs[parent_key] = dict(
                case=case, endpoint=trace["endpoint"], file=str(origin.relative_to(REPO)),
                parent_cpu=parent.get("whole_process_CPU_seconds", parent.get("cpu_seconds")),
                parent_wall=parent.get("wall_seconds"),
                timing_type="historical R13 / i9 CPU, separate clock" if case.startswith("r13_")
                else "historical R6 CPU, separate clock")
    methods = ["pair_secant"]
    if "pair_exact" in receipts[0][0]["config"]["arms"]:
        methods.append("pair_exact")
    pairs, summary = robustness(selected, methods)
    folds, predictions = cross_validate([r for r in selected if r["source"] == "selected"],
                                         receipts[0][0]["config"]["folds"])
    mark_pareto(prefixes)
    write_csv(out / "OUTCOMES.csv", selected)
    write_csv(out / "PAID_PREFIXES.csv", prefixes)
    write_csv(out / "FORECAST_AUDIT.csv", forecasts)
    write_csv(out / "HISTORICAL_ENDPOINT_COSTS.csv", list(base_costs.values()))
    write_csv(out / "PAIRED_ENDPOINTS.csv", pairs)
    write_csv(out / "CV_HELD_OUT.csv", predictions)
    json_write(out / "CV_FOLDS.json", folds)
    json_write(out / "PAIRED_OBJECT_SUMMARY.json", summary)
    json_write(out / "PROVENANCE.json", provenance)
    if plots:
        variant = frozen.name
        gallery(case_data, selected, pictures, figdir, variant)
        trajectories(prefixes, figdir, variant)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sealed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--figures", type=Path, default=LAB / "figures")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    started = time.process_time()
    score(args.sealed, args.out, args.figures, not args.no_plots)
    json_write(args.out / "SCORING_RECEIPT.json",
               {"seal": str(args.sealed), "post_seal_scoring": True,
                "process_cpu_seconds": time.process_time() - started, "plots": not args.no_plots})


if __name__ == "__main__":
    main()
