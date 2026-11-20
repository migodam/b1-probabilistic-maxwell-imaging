"""Summarize paired short runs and the separately archived R10 long schedule."""
from native_current import ADDITION_ARMS, ARMS, CASE_COHORT, ROOT

import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cohorts import historical_study, load_truth
from score_materials import material_metrics

LABEL = {
    "native_fft": "Native FFT",
    "small_fft": "Narrow FFT",
    "small_swap": "Narrow FFT + six swaps",
    "small_add_fourier": "Add six Fourier",
    "small_add_residual": "Add six residual",
    "archived_auxiliary": "R10 long FFT-CIE",
    "archived_corrector": "R10 long + corrector",
}
COLOR = {"native_fft": "#1a5973", "small_fft": "#bd8029", "small_swap": "#963b78"}


def quartiles(number):
    return {key: float(v) for key, v in zip(("q25", "median", "q75"), np.percentile(number, [25, 50, 75]))}


def historical_record(case_id, truth, stage):
    study = historical_study(10) / "RESULTS" / "twod" / f"{case_id}__global_dual__noise_ball_99"
    with np.load(study / f"{stage}.npz") as arrays:
        material = material_metrics(arrays["m"], truth)
    candidate = json.loads((study / f"{stage}.json").read_text())
    return {
        "true_data_residual": candidate["physical_relative_residual"],
        "physical_objective": candidate["objective"],
        "materials": material,
    }


def aggregate(cases, methods=(*ARMS, "archived_auxiliary", "archived_corrector")):
    scorecard = {}
    for method in methods:
        records = [record["arms"][method] for record in cases]
        residual = [r["true_data_residual"] for r in records]
        success = [r["materials"]["material_success"] for r in records]
        noise = [r["relative_rms_noise"] for r in cases]
        scorecard[method] = {
            "N": len(records),
            "true_data_residual": quartiles(residual),
            "count_within_2x_rms_noise": int(np.count_nonzero(np.divide(residual, noise) <= 2)),
            "count_material_success": int(sum(success)),
            **{
                metric: quartiles([r["materials"][metric] for r in records])
                for metric in ("RMSE", "amplitude_NRMSE", "fixed_support_IoU", "SSIM_fixed_range")
            },
        }
        if method in (*ARMS, *ADDITION_ARMS):
            for clock in ("wall_s", "cpu_s"):
                value = [r["timing"][clock] for r in records]
                scorecard[method][f"{clock}_per_case"] = quartiles(value)
                scorecard[method][f"{clock}_sum"] = float(sum(value))
            scorecard[method]["actual_stage_iterations_values"] = sorted({
                tuple(s["iterations"] for s in r["stages"]) for r in records
            })
            scorecard[method]["material_inner_evaluations"] = quartiles([
                r["operations"]["material_inner_evaluations"] for r in records
            ])
            scorecard[method]["full_wave_candidate_LU_sum"] = sum(r["full_wave_candidate_lu"] for r in records)
            scorecard[method]["full_wave_candidate_rhs_sum"] = sum(r["full_wave_candidate_rhs"] for r in records)
            scorecard[method]["best_stage_counts"] = {
                str(i): sum(r["selected_stage"] == i for r in records) for i in range(-1, 4)
            }
    return scorecard


def comparison(cases, arm, baseline):
    paired = []
    for record in cases:
        left, right = [record["arms"][key] for key in (arm, baseline)]
        paired.append([
            left["true_data_residual"] - right["true_data_residual"],
            left["materials"]["RMSE"] - right["materials"]["RMSE"],
            left["materials"]["fixed_support_IoU"] - right["materials"]["fixed_support_IoU"],
        ])
    paired = np.asarray(paired)
    return {
        key: dict(
            delta=quartiles(paired[:, idx]),
            negative=int(np.count_nonzero(paired[:, idx] < -1e-12)),
            near_tie=int(np.count_nonzero(abs(paired[:, idx]) <= 1e-12)),
            positive=int(np.count_nonzero(paired[:, idx] > 1e-12)),
        )
        for idx, key in enumerate(("true_data_residual", "RMSE", "IoU"))
    }


def make_figure(cases, stats):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(12.1, 3.9), layout="constrained")
    ax = axes[0]
    native = np.array([case["arms"]["native_fft"]["true_data_residual"] for case in cases])
    for method in ARMS[1:]:
        residual = [case["arms"][method]["true_data_residual"] for case in cases]
        ax.scatter(native, residual, s=30, color=COLOR[method], alpha=.8,
                   marker="o" if method == "small_fft" else "^", label=LABEL[method])
    ax.plot([.008, 1.2], [.008, 1.2], color="#8e9eaa", lw=1, ls="--")
    ax.set(xscale="log", yscale="log", xlim=(.008, 1.2), ylim=(.008, 1.2),
           xlabel="Native FFT residual", ylabel="Small-model residual",
           title="a   Full-wave fit after inversion")
    ax.legend(frameon=False, loc="upper left", fontsize=8)

    ax = axes[1]
    for y, name in enumerate(ARMS):
        arm = stats[name]
        q = arm["cpu_s_per_case"]
        ax.errorbar(q["median"], y, xerr=[[q["median"] - q["q25"]], [q["q75"] - q["median"]]],
                    ms=7, fmt="o", color=COLOR[name], capsize=3)
    ax.set(yticks=range(len(ARMS)), yticklabels=[LABEL[n] for n in ARMS],
           xlabel="Total CPU seconds per case", title="b   Complete inversion time")
    ax.invert_yaxis()

    ax = axes[2]
    for marker, name in zip(("o", "^", "s"), ARMS):
        ax.plot(range(1, 5), [
            np.median([case["arms"][name]["stages"][stage]["true_data_residual"] for case in cases])
            for stage in range(4)
        ], marker=marker, color=COLOR[name], label=LABEL[name], lw=1.8)
    ax.set(xticks=[1, 2, 3, 4], xlabel="Stage", ylabel="Median full-wave residual",
           title="c   Fit by continuation stage")

    for ax in axes:
        ax.grid(alpha=.2, color="#8c99a8")
        ax.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"native_current_comparison.{suffix}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def additions_analysis(original, additions):
    together = []
    stage_agreement = []
    for case_id, row in additions.items():
        old = original[case_id]["arms"]
        added = row["arms"]
        anchors = [r["stages"][0]["true_data_residual"] for r in (*old.values(), *added.values())
                   if "stages" in r]
        stage_agreement.append(max(anchors) - min(anchors))
        together.append({**row, "arms": {**old, **added}})
    methods = (*ARMS, *ADDITION_ARMS)
    groups = {
        "all24": together,
        "weak8": [c for c in together if c["truth_peak"] <= 2.5 + 1e-6],
        "strong16": [c for c in together if c["truth_peak"] > 2.5 + 1e-6],
    }
    pairs = [
        (new_arm, reference) for new_arm in ADDITION_ARMS for reference in ("small_fft", "small_swap", "native_fft")
    ] + [("small_add_residual", "small_add_fourier")]
    rank_report = {key: {} for key in ADDITION_ARMS}
    for key in ADDITION_ARMS:
        rank_records = [c["arms"][key]["final_basis_rank"] for c in together]
        rank_report[key] = {
            "observed_final_ranks": sorted({r["effective_dark_rank"] for r in rank_records}),
            "smallest_basis_singular_value": quartiles([r["smallest_singular_value"] for r in rank_records]),
        }
    results = {
        "scorecard": {group: aggregate(cases, methods) for group, cases in groups.items()},
        "paired": {f"{a}_against_{b}": comparison(together, a, b) for a, b in pairs},
        "observed_final_basis_ranks": rank_report,
        "shared_first_stage_max_physical_residual_delta": max(stage_agreement),
    }
    path = ROOT / "results" / "native_current_additions_summary.json"
    path.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9), layout="constrained")
    for i, (label, bounds) in enumerate((("Full-wave residual", (.009, 1.1)), ("Material RMSE", (.2, 10)))):
        ax = axes[i]
        for method, color, marker in zip(ADDITION_ARMS, ("#bd8029", "#963b78"), ("o", "^")):
            narrow = [c["arms"]["small_fft"] for c in together]
            improved = [c["arms"][method] for c in together]
            x, y = (
                ([r["true_data_residual"] for r in narrow], [r["true_data_residual"] for r in improved])
                if i == 0 else ([r["materials"]["RMSE"] for r in narrow], [r["materials"]["RMSE"] for r in improved])
            )
            ax.scatter(x, y, alpha=.75, s=34, marker=marker, color=color, label=LABEL[method])
        ax.plot(bounds, bounds, ls="--", lw=1, color="#8e9eaa")
        ax.set(xscale="log", yscale="log", xlim=bounds, ylim=bounds,
               xlabel=f"{label} · original narrow FFT", ylabel=f"{label} · with added modes")
        ax.grid(alpha=.2, color="#8c99a8")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(loc="upper left", frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"native_current_additions.{suffix}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    result = json.loads((ROOT / "results" / "native_current.json").read_text())
    original = result["cases"]
    for case_id, case in original.items():
        truth = load_truth(CASE_COHORT, case_id)
        case["arms"]["archived_auxiliary"] = historical_record(case_id, truth, "selected_before_corrector")
        case["arms"]["archived_corrector"] = historical_record(case_id, truth, "prediction")
    cases = list(original.values())
    stats = {
        "all24": aggregate(cases),
        "weak8": aggregate([c for c in cases if c["truth_peak"] <= 2.5 + 1e-6]),
        "strong16": aggregate([c for c in cases if c["truth_peak"] > 2.5 + 1e-6]),
    }
    paired = {f"{arm}_against_{baseline}": comparison(cases, arm, baseline)
              for arm, baseline in (("small_fft", "native_fft"), ("small_swap", "native_fft"), ("small_swap", "small_fft"))}
    diagnostics = {
        "shared_first_stage_max_physical_residual_delta": max(
            c["identical_initial_stage_data_residual_max_difference"] for c in cases
        ),
        "relative_projected_gradient_energy_retained": quartiles([
            c["arms"]["small_swap"]["basis_replacement"]["retained_singular_energy"] for c in cases
        ]),
        "retained_fourier_overlap_max": max(
            c["arms"]["small_swap"]["basis_replacement"]["retained_dark_overlap"] for c in cases
        ),
        "receiver_overlap_max": max(c["arms"]["small_swap"]["basis_replacement"]["receiver_overlap"] for c in cases),
    }
    archived = {case_id: {arm: case["arms"][arm] for arm in ("archived_auxiliary", "archived_corrector")}
                for case_id, case in original.items()}
    report = {"scorecard": stats, "paired": paired, "diagnostics": diagnostics, "historical_references": archived}
    (ROOT / "results" / "native_current_summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    make_figure(cases, stats["all24"])
    additions = ROOT / "results" / "native_current_additions.json"
    if additions.exists():
        additions_analysis(original, json.loads(additions.read_text())["cases"])


if __name__ == "__main__":
    main()
