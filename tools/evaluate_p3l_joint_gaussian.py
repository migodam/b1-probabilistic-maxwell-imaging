"""Isolated offline DEV scoring for the completed joint-probit GN revision.

No solver, Target, backend or Book imports. Load and freeze online estimates
before reading DEV chi; never modify proposals, select parents, tune quadrature
or qualify a posterior using truth. IS remains an unqualified diagnostic.
"""
from pathlib import Path
import argparse
import json
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "research/p3l_local_posterior_imaging"
PILOT = ["DEV_000", "DEV_028", "DEV_030", "DEV_026"]
ESTIMATORS = ["historical_native_q_anchor", "decoder_at_x_mode",
              "fixed_property_gaussian_mean_at_x_mode", "joint_GH5_mean",
              "joint_GH9_mean", "joint_GH17_mean", "full_target_IS_unqualified_mean"]
BASELINES = ESTIMATORS[:3]
CAP = 128 * 1024 * 1024


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, allow_nan=False, indent=2) + "\n")
    tmp.replace(path)


def load(path, names):
    with zipfile.ZipFile(path) as archive:
        for name in names:
            if archive.getinfo(name + ".npy").file_size > CAP:
                raise MemoryError("array exceeds frozen 128MiB bound")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in names}
    for value in arrays.values():
        if value.nbytes > CAP or not np.all(np.isfinite(value)):
            raise ValueError("nonfinite or oversized array")
    return arrays


def vector(path, name):
    value = load(path, [name])[name]
    if value.shape != (576,):
        raise ValueError("material estimate must have shape (576,)")
    return value


def configuration(figures):
    return dict(revision="P3L_B_JOINT_GAUSSIAN_OFFLINE_EVALUATION_001",
        parents=[dict(parent=f"DEV_{i:03d}", role="train" if i < 24 else "validation") for i in range(32)],
        estimators=ESTIMATORS, baselines=BASELINES,
        losses=["complex_MSE", "complex_MAE"], masks=["full_grid", "ROI"],
        difference="estimator risk minus baseline risk; negative favors estimator",
        bootstrap=dict(replicates=2000, seed=2026102295, interval=[.025, .975],
                       unit="paired parent resampling within fixed cohort; mean differences",
                       missing="Any missing parent makes cohort statistic null; never remove missing parents"),
        cohorts=["all_development", "train", "validation", "fixed_IS_pilot"],
        IS_pilot=PILOT, figures=bool(figures),
        figure_rules="first four fixed y-only parents; shared epsilon Re and Im ranges; uncalibrated Gaussian modulus SD separately",
        phases="closed online job -> frozen evaluation config -> all estimates loaded -> truth chi reads -> risks",
        qualifications="GN Gaussian-reference images and unqualified IS only; no posterior/teacher/coverage qualification",
        new_physics_RHS=0, online_proposals=False, GT_selection=False, heldout_access=False)


def gather(out, declared, online_rows):
    estimates = {}
    rows = []
    for declaration in declared:
        parent = declaration["parent"]
        folder = out / f"B/joint_gaussian/{parent}"
        online = online_rows.get(parent, {})
        row = dict(**declaration, online_status=online.get("status", "ONLINE_ROW_MISSING"),
                   stationarity=online.get("stationarity"), estimator_failures={}, risks={}, differences={})
        values = {}
        sources = {
            "historical_native_q_anchor": (folder / "NATIVE_ANCHOR_CONTROL.npz", "chi"),
            "decoder_at_x_mode": (folder / "CONTROLS.npz", "decoder_at_probit_mode"),
            "fixed_property_gaussian_mean_at_x_mode": (folder / "CONTROLS.npz", "fixed_property_mode_gaussian_mean"),
            **{f"joint_GH{order}_mean": (folder / f"JOINT_GH_ORDER_{order}.npz", "chi_mean") for order in (5, 9, 17)},
            "full_target_IS_unqualified_mean": (folder / "FULL_TARGET_IS.npz", "weighted_material_mean")}
        for name, (path, key) in sources.items():
            try:
                values[name] = vector(path, key)
            except Exception as error:
                row["estimator_failures"][name] = dict(status="MISSING_OR_FAILED_RETAINED", error=repr(error))
        # The historical optimized anchor is available even if the new paid job
        # failed before copying its control. Reuse only that exact pre-existing
        # A anchor; never choose between candidate anchors using offline risk.
        if "historical_native_q_anchor" not in values:
            try:
                values["historical_native_q_anchor"] = vector(out / f"A/anchors/{parent}.npz", "chi")
                row["historical_anchor_source"] = "A/anchors: exact historical anchor; missing revision copy"
                row["estimator_failures"].pop("historical_native_q_anchor", None)
            except Exception as error:
                row["estimator_failures"]["historical_native_q_anchor"] = dict(status="MISSING_OR_FAILED_RETAINED", error=repr(error))
        for filename, label in (("GH_REFINEMENT.json", "online_refinement"),
                                ("FULL_TARGET_IS.json", "IS_diagnostic_unqualified"),
                                ("MODE_OPTIMIZATION.json", "mode_optimization")):
            path = folder / filename
            if path.exists():
                try:
                    content = json.loads(path.read_text())
                    # Avoid repeating the full optimizer trace in the risk file.
                    if label == "mode_optimization":
                        content.pop("history", None)
                    row[label] = content
                except Exception as error:
                    row[label] = dict(status="METADATA_READ_FAILED", error=repr(error))
        row["saved_array_refinement"] = {}
        for lower, higher in ((5, 9), (9, 17)):
            a, b = values.get(f"joint_GH{lower}_mean"), values.get(f"joint_GH{higher}_mean")
            key = f"GH{higher}_minus_GH{lower}"
            row["saved_array_refinement"][key] = (None if a is None or b is None else
                dict(complex_RMS=float(np.sqrt(np.mean(np.abs(b-a)**2))),
                     complex_MAE=float(np.mean(np.abs(b-a))),
                     complex_max=float(np.max(np.abs(b-a))),
                     interpretation="quadrature refinement only; not Gaussian adequacy or posterior error"))
        try:
            variance = vector(folder / "JOINT_GH_ORDER_17.npz", "chi_variance")
            if np.iscomplexobj(variance) or np.any(variance < 0):
                raise ValueError("Gaussian modulus variance must be real and nonnegative")
            values["uncalibrated_gaussian_modulus_SD"] = np.sqrt(variance)
        except Exception as error:
            row["uncertainty_status"] = dict(status="MISSING_OR_FAILED_RETAINED", error=repr(error))
        estimates[parent] = values
        rows.append(row)
    return estimates, rows


def risks(values, truth, masks):
    return {name: {label: dict(complex_MSE=float(np.mean(np.abs(value[mask]-truth[mask])**2)),
                               complex_MAE=float(np.mean(np.abs(value[mask]-truth[mask]))))
                    for label, mask in masks.items()}
            for name, value in values.items() if name in ESTIMATORS}


def aggregation(rows, config):
    results = {}
    for cohort_name in config["cohorts"]:
        cohort = [r for r in rows if cohort_name == "all_development" or
                  (cohort_name == "fixed_IS_pilot" and r["parent"] in PILOT) or r["role"] == cohort_name]
        bootstrap = config["bootstrap"]
        rng = np.random.default_rng(bootstrap["seed"])
        indices = rng.integers(0, len(cohort), size=(bootstrap["replicates"], len(cohort)))
        summary = dict(parent_count=len(cohort), parents=[r["parent"] for r in cohort],
                       mean_risks={}, paired_differences={})
        for estimator in ESTIMATORS:
            for mask in config["masks"]:
                for metric in config["losses"]:
                    values = [r.get("risks", {}).get(estimator, {}).get(mask, {}).get(metric) for r in cohort]
                    missing = [r["parent"] for r, value in zip(cohort, values) if value is None]
                    summary["mean_risks"][f"{estimator}/{mask}/{metric}"] = dict(
                        mean=None if missing else float(np.mean(values)), missing_parents=missing,
                        status="INCOMPLETE_COHORT_NO_MISSING_PARENT_REMOVAL" if missing else "DESCRIPTIVE_COMPLETE_DEV_COHORT")
                    for baseline in BASELINES:
                        if estimator == baseline:
                            continue
                        values = [r.get("differences", {}).get(estimator, {}).get(baseline, {}).get(mask, {}).get(metric) for r in cohort]
                        missing = [r["parent"] for r, value in zip(cohort, values) if value is None]
                        item = dict(mean_difference=None, paired_bootstrap_95_interval=None, missing_parents=missing,
                                    status="INCOMPLETE_COHORT_NO_MISSING_PARENT_REMOVAL" if missing else "DESCRIPTIVE_COMPLETE_DEV_COHORT")
                        if not missing:
                            values = np.asarray(values)
                            item.update(mean_difference=float(values.mean()),
                                        paired_bootstrap_95_interval=np.quantile(values[indices].mean(axis=1), bootstrap["interval"]).tolist())
                        summary["paired_differences"][f"{estimator}_minus_{baseline}/{mask}/{metric}"] = item
        results[cohort_name] = summary
    return results


def figures(out, estimates, truth_by_parent):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    columns = [("truth", "DEV truth"), (ESTIMATORS[0], "Historical native-q anchor"),
               (ESTIMATORS[1], "Probit coordinate mode"), (ESTIMATORS[2], "Fixed-property Gaussian mean"),
               ("joint_GH17_mean", "Joint GN Gaussian mean (GH17)"),
               (ESTIMATORS[-1], "IS diagnostic (unqualified)")]
    values = {p: dict(estimates[p], truth=truth_by_parent[p]) for p in PILOT if p in truth_by_parent}
    destination = out / "B/joint_gaussian_figures"
    destination.mkdir(parents=True, exist_ok=True)
    report = {}
    for component in ("real", "imag"):
        displayed = [getattr(v[name], component) + (1 if component == "real" else 0)
                     for v in values.values() for name, _ in columns if name in v]
        if not displayed:
            report[component] = dict(status="NO_ARRAYS_AVAILABLE")
            continue
        lo, hi = min(float(a.min()) for a in displayed), max(float(a.max()) for a in displayed)
        fig, axes = plt.subplots(4, len(columns), figsize=(19, 12), squeeze=False, constrained_layout=True)
        image = None
        for i, parent in enumerate(PILOT):
            for j, (name, title) in enumerate(columns):
                ax = axes[i, j]
                ax.set_title(title if i == 0 else "", fontsize=9)
                if name in values.get(parent, {}):
                    field = getattr(values[parent][name], component) + (1 if component == "real" else 0)
                    image = ax.imshow(field.reshape(24, 24), origin="lower", vmin=lo, vmax=hi, cmap="viridis")
                else:
                    ax.text(.5, .5, "MISSING / FAILED", ha="center", va="center", transform=ax.transAxes, fontsize=8)
                ax.set_xticks([]); ax.set_yticks([])
                if j == 0:
                    ax.set_ylabel(parent)
        if image is not None:
            fig.colorbar(image, ax=axes.ravel().tolist(), shrink=.7, label=f"epsilon {component} (shared scale)")
        fig.suptitle("Fixed DEV panel: Gaussian-reference images; no full-posterior qualification")
        path = destination / f"fixed_panel_epsilon_{component}.png"
        fig.savefig(path, dpi=150); plt.close(fig)
        report[component] = dict(path=str(path.relative_to(out)), shared_min=lo, shared_max=hi)
    fig, axes = plt.subplots(1, 4, figsize=(14, 4), squeeze=False, constrained_layout=True)
    spreads = [estimates[p]["uncalibrated_gaussian_modulus_SD"] for p in PILOT
               if "uncalibrated_gaussian_modulus_SD" in estimates[p]]
    upper = max((float(v.max()) for v in spreads), default=1.)
    image = None
    for ax, parent in zip(axes[0], PILOT):
        value = estimates[parent].get("uncalibrated_gaussian_modulus_SD")
        if value is not None:
            image = ax.imshow(value.reshape(24, 24), origin="lower", vmin=0, vmax=upper, cmap="magma")
        else:
            ax.text(.5, .5, "MISSING / FAILED", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(parent); ax.set_xticks([]); ax.set_yticks([])
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=.8, label="sqrt(E|chi-Echi|²), uncalibrated")
    fig.suptitle("GN Gaussian-reference modulus SD (GH17): uncalibrated; no credible-region claim")
    path = destination / "fixed_panel_gaussian_uncalibrated_SD.png"
    fig.savefig(path, dpi=150); plt.close(fig)
    report["uncertainty"] = dict(path=str(path.relative_to(out)), shared_min=0., shared_max=upper,
                                 interpretation="Gaussian-reference complex modulus SD, not component SD or calibrated posterior uncertainty")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--figures", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    # A final closure must exist: a still-running online job cannot be scored.
    closure = out / "B/jobs/joint_probit_gaussian/CLOSED.json"
    if not closure.exists() or json.loads(closure.read_text()).get("status") not in ("COMPLETE", "ERROR"):
        raise RuntimeError("online joint Gaussian job has not finished; no offline reads permitted")
    online_config = json.loads((out / "configs/B_JOINT_GAUSSIAN_REVISION.json").read_text())
    config = configuration(args.figures)
    if online_config["parents"] != [r["parent"] for r in config["parents"]]:
        raise RuntimeError("online all32 cohort differs")
    declarations = json.loads((out / "configs/PARENTS.json").read_text())
    development = [dict(parent=r["parent"], role=r["role"]) for r in declarations if r["group"] == "development"]
    if development != config["parents"]:
        raise RuntimeError("development identities/order/roles differ from frozen evaluation")
    frozen_path = out / "configs/B_JOINT_GAUSSIAN_EVALUATION.json"
    if frozen_path.exists():
        if json.loads(frozen_path.read_text()) != config:
            raise RuntimeError("evaluation config differs from frozen revision")
    else:
        with frozen_path.open("x") as stream:
            stream.write(json.dumps(config, allow_nan=False, indent=2) + "\n")
    online = json.loads((out / "B/JOINT_GAUSSIAN_REFERENCE.json").read_text())
    online_rows = {r["parent"]: r for r in online["rows"]}
    estimates, rows = gather(out, config["parents"], online_rows)
    write(out / "B/JOINT_GAUSSIAN_EVALUATION_INPUTS.json", dict(
        phase="ALL_ONLINE_ESTIMATES_LOADED_BEFORE_TRUTH", truth_read=False, rows=rows,
        online_closure=json.loads(closure.read_text()), all32_retained=True))
    cfg = json.loads((out / "configs/TARGET.json").read_text())
    points = np.asarray(json.loads((ROOT / "research/p1_material_target_20261018/inputs/geometry.json").read_text())["points"])
    masks = dict(full_grid=np.ones(576, bool), ROI=np.linalg.norm(points, axis=1) <= cfg["roi_radius"])
    truth_by_parent = {}
    for row in rows:
        parent = row["parent"]
        try:
            truth = vector(out / f"inputs/offline_only/{parent}.npz", "chi")
            row["risks"] = risks(estimates[parent], truth, masks)
            row["differences"] = {name: {baseline: {mask: {
                metric: row["risks"][name][mask][metric] - row["risks"][baseline][mask][metric]
                for metric in config["losses"]} for mask in config["masks"]}
                for baseline in BASELINES if baseline != name and baseline in row["risks"]}
                for name in row["risks"]}
            row["evaluation_status"] = "OFFLINE_DEV_RISKS_RECORDED_NO_QUALIFICATION"
            if args.figures and parent in PILOT:
                truth_by_parent[parent] = truth.copy()
        except Exception as error:
            row.update(evaluation_status="TRUTH_EVALUATION_FAILED_PARENT_RETAINED", evaluation_error=repr(error))
    report = dict(config="configs/B_JOINT_GAUSSIAN_EVALUATION.json", all32_retained=True,
                  online_closure_status=json.loads(closure.read_text())["status"],
                  rows=rows, cohorts=aggregation(rows, config),
                  new_physics_RHS=0, GT_selection=False, online_proposals=False,
                  reference="JOINT_PROBIT_GN_GAUSSIAN_REFERENCE_NOT_FULL_POSTERIOR",
                  IS="DISTINCT_UNQUALIFIED_FULL_TARGET_IMPORTANCE_APPROXIMATION",
                  interpretation="Descriptive DEV risks; coordinate optimization and Gaussian averaging controls separated; no posterior/teacher/coverage qualification from GT")
    if args.figures:
        try:
            report["figures"] = figures(out, estimates, truth_by_parent)
        except Exception as error:
            report["figures"] = dict(status="FIGURE_FAILURE_RISKS_RETAINED", error=repr(error))
    destination = out / "B/JOINT_GAUSSIAN_IMAGE_RISK.json"
    write(destination, report)
    print(json.dumps(dict(status="OFFLINE_DEV_EVALUATION_RECORDED", parents=len(rows),
                         report=str(destination), new_physics_RHS=0)))


if __name__ == "__main__":
    main()
