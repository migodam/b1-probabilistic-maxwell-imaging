"""Score the current-state approximation at a common material and draw figures."""
from state_rom import ROOT, COHORTS, build_models

import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from cohorts import load_cohort, load_truth

ARCHIVES = (*COHORTS, "r11_crossgrid8_matched")
ARMS = ("full", "tsom36", "tsom80", "tsom144", "fft80", "tsom80_enriched")
NAMES = {
    "full": "Full",
    "tsom36": "TSOM 36",
    "tsom80": "TSOM 80",
    "tsom144": "TSOM 144",
    "fft80": "FFT 80",
    "tsom80_enriched": "TSOM 80→92",
}
COLORS = {
    "full": "#172b40",
    "tsom36": "#60aeb1",
    "tsom80": "#207a8c",
    "tsom144": "#11485d",
    "fft80": "#e69a47",
    "tsom80_enriched": "#87367b",
}
SHORT = {
    "r10_twod": "R10 · 24",
    "r11_independent48": "R11 independent · 48",
    "r11_crossgrid8": "R11 cross-grid · 8",
    "r11_bothblocks16": "R11 both-blocks · 16",
    "r11_crossgrid8_matched": "R11 matched pair · 8",
    "r10_weak": "R10 weak contrast · 8",
    "r10_strong": "R10 strong contrast · 16",
}


def statistics(array):
    a = np.asarray(array)
    return dict(zip(("min", "q25", "median", "q75", "max"), map(float, np.percentile(a, [0, 25, 50, 75, 100]))))


def fixed_state_comparison(models, material, Y):
    full = models["full"]
    J_true, Y_true, _ = full.forward(material)
    enriched, _ = models["tsom80"].enrich(material)
    operators = {**models, "tsom80_enriched": enriched}
    report = {}
    for key in ARMS[1:]:
        arm = operators[key]
        z, Y_rom, _ = arm.forward(material)
        J_rom = arm.lift(z)
        rhs = material[:, None] * arm.E0
        defect = J_rom - material[:, None] * arm.field(z)
        report[key] = {
            "relative_current_error": float(np.linalg.norm(J_rom - J_true) / np.linalg.norm(J_true)),
            "relative_state_defect": float(np.linalg.norm(defect) / np.linalg.norm(rhs)),
            "relative_receiver_error": float(np.linalg.norm(Y_rom - Y_true) / np.linalg.norm(Y)),
        }
    return report


def aggregate(cases):
    summary = {}
    for key in ARMS:
        values = [c["roles"][key] for c in cases]
        speed = [c["roles"]["full"]["warm_reconstruction"]["wall_s"] / c["roles"][key]["warm_reconstruction"]["wall_s"] for c in cases]
        speed_check = [c["roles"]["full"]["warm_including_terminal_validation"]["wall_s"] / c["roles"][key]["warm_including_terminal_validation"]["wall_s"] for c in cases]
        noise_ratio = [value["true_data_residual"] / c["relative_rms_noise"] for c, value in zip(cases, values)]
        summary[key] = {
            key: statistics([value[key] for value in values])
            for key in (
                "true_data_residual", "model_data_residual", "state_defect_rel_rhs",
                "material_rmse", "material_amplitude_nrmse", "material_support_iou",
                "objective_evaluations",
            )
        }
        summary[key].update(
            n_cases=len(cases),
            physical_data_residual_le_2x_rms_noise=int(np.count_nonzero(np.asarray(noise_ratio) <= 2)),
            ratio_to_rms_noise=statistics(noise_ratio),
            paired_warm_wall_speedup=statistics(speed),
            paired_warm_wall_speedup_with_validation=statistics(speed_check),
        )
        for timing in (
            "warm_reconstruction", "cold_reconstruction",
            "warm_including_terminal_validation", "cold_including_terminal_validation",
            "shared_operator_and_method_setup",
        ):
            for clock in ("wall_s", "cpu_s"):
                summary[key][f"{timing}_{clock}"] = statistics([value[timing][clock] for value in values])
        if key != "full":
            fixed = [c["fixed_material_state_comparison"][key] for c in cases]
            summary[key]["fixed_material_state_comparison"] = {
                k: statistics([value[k] for value in fixed]) for k in fixed[0]
            }
    return summary


def read_all():
    by_archive = {}
    for name in ARCHIVES:
        out = ROOT / "results" / f"state_{name}"
        saved = json.loads(out.with_suffix(".json").read_text())
        model, cases = load_cohort(name)
        models, _, _, _ = build_models(model)
        full = models["full"]
        with np.load(out.with_suffix(".npz")) as predictions:
            for case in cases:
                case_id, Y = case["id"], case["Y"]
                record = saved["cases"][case_id]
                truth = np.asarray(load_truth(name, case_id)).ravel()
                _, Y_truth, _ = full.forward(truth)
                expected_noise = case["sigma"] * np.sqrt(Y.size)
                record["relative_rms_noise"] = float(expected_noise / np.linalg.norm(Y))
                record["truth_max"] = float(truth.max())
                record["archive_truth_full_residual_over_rms_noise"] = float(np.linalg.norm(Y_truth - Y) / expected_noise)
                record["fixed_material_state_comparison"] = fixed_state_comparison(
                    models, predictions[f"{case_id}::tsom80"], Y,
                )
        by_archive[name] = saved
    return by_archive


def suite_cost(by_archive):
    main = [by_archive[name] for name in COHORTS]
    cost = {}
    for arm in ARMS:
        rows = [case["roles"][arm] for archived in main for case in archived["cases"].values()]
        cost[arm] = {"warm_factorizations": sum(row["warm_factorizations"] for row in rows)}
        for clock in ("wall_s", "cpu_s"):
            setup = sum(
                archived["shared_cohort_archive_load"][clock]
                + next(iter(archived["cases"].values()))["roles"][arm]["shared_operator_and_method_setup"][clock]
                for archived in main
            )
            for scope in ("warm_reconstruction", "warm_including_terminal_validation"):
                cost[arm][f"suite_{scope}_{clock}"] = setup + sum(row[scope][clock] for row in rows)
    return cost


def save_figure(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def compare_plot(by_archive, summaries):
    cases = [c for name in COHORTS for c in by_archive[name]["cases"].values()]
    fig, axes = plt.subplots(2, 2, figsize=(11.7, 7.7), layout="constrained")
    ax = axes[0, 0]
    for key in ARMS[1:]:
        fixed = [c["fixed_material_state_comparison"][key]["relative_current_error"] for c in cases]
        ax.plot(np.sort(fixed), np.arange(1, len(fixed) + 1) / len(fixed), color=COLORS[key], lw=2)
    ax.set(xscale="log", xlabel="Relative error of the lifted current", ylabel="Fraction of cases",
           title="a    Current error at a common material")

    ax = axes[0, 1]
    for i, key in enumerate(ARMS[1:]):
        s = summaries["main_96"][key]
        t = s["paired_warm_wall_speedup"]
        tv = s["paired_warm_wall_speedup_with_validation"]
        ax.errorbar(t["median"], i, xerr=[[t["median"] - t["q25"]], [t["q75"] - t["median"]]],
                    fmt="o", color=COLORS[key], ms=7, capsize=3)
        ax.plot(tv["median"], i, marker="D", mfc="white", mec=COLORS[key], ms=6)
    ax.axvline(1, color="#aab3bf", lw=1)
    ax.set(yticks=range(len(ARMS) - 1), yticklabels=[NAMES[key] for key in ARMS[1:]],
           xlabel="Speedup over full inversion (×)", title="b    End-to-end inversion time")
    ax.invert_yaxis()
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="#4e5a66", label="Warm reconstruction (IQR)"),
            Line2D([0], [0], marker="D", mfc="white", color="#4e5a66", linestyle="None", label="Including full-wave check"),
        ],
        fontsize=8, frameon=False, loc="upper left",
    )

    ax = axes[1, 0]
    x = np.arange(1, len(ARMS) + 1)
    groups = ("r10_weak", "r10_strong", "r11_independent48", "r11_bothblocks16", "r11_crossgrid8", "r11_crossgrid8_matched")
    styles = ("-", "-", "--", ":", "-.", (0, (3, 1, 1, 1, 1, 1)))
    for group, style in zip(groups, styles):
        medians = [summaries[group][key]["material_rmse"]["median"] for key in ARMS]
        color = "#258d76" if group == "r10_weak" else "#563f84" if group == "r10_strong" else "#65758b"
        ax.plot(x, medians, linestyle=style, lw=2, color=color, alpha=.85, label=SHORT[group])
    ax.set(title="c    Recovered material error", ylabel="Material RMSE",
           xticks=x, xticklabels=[NAMES[k] for k in ARMS])
    ax.tick_params(axis="x", labelrotation=20)
    ax.legend(fontsize=7, frameon=False, ncols=2, loc="center left")

    ax = axes[1, 1]
    for i, key in enumerate(ARMS):
        q = [c["roles"][key]["true_data_residual"] for c in cases]
        ax.plot(np.sort(q), np.arange(1, len(q) + 1) / len(q), color=COLORS[key],
                lw=2, ls=["-", "-", "--", ":", "-.", "-"][i], label=NAMES[key])
    ax.set(xscale="log", xlabel="Relative full-wave data residual",
           ylabel="Fraction of cases", title="d    Data residual after full-wave recalculation")
    ax.legend(fontsize=8, frameon=False, loc="upper left", ncols=2)
    for ax in axes.ravel():
        ax.grid(alpha=.2, color="#8c99a8")
        ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, "state_rom_accuracy_time")


def recon_plot(by_archive):
    examples = (
        ("r10_twod", "r10_twod_0000", "R10 0000 · weak annulus"),
        ("r10_twod", "r10_twod_0006", "R10 0006 · weak islands"),
        ("r10_twod", "r10_twod_0001", "R10 0001 · strong annulus"),
        ("r11_crossgrid8", "r11_crossgrid8_0000", "R11 0000 · cross-grid"),
    )
    methods = ("full", "tsom36", "tsom80", "tsom144", "tsom80_enriched")
    fig, axes = plt.subplots(len(examples), len(methods) + 1, figsize=(11.9, 7.7), layout="constrained")
    for row, (name, case_id, label) in enumerate(examples):
        shape = by_archive[name]["grid_shape"]
        truth = load_truth(name, case_id).reshape(shape)
        roles = by_archive[name]["cases"][case_id]["roles"]
        axes[row, 0].set_ylabel(label, fontsize=9)
        with np.load(ROOT / "results" / f"state_{name}.npz") as predictions:
            mats = [truth] + [predictions[f"{case_id}::{key}"].reshape(shape) for key in methods]
        vmax = max(truth.max(), *[np.quantile(m, .995) for m in mats])
        for col, mat in enumerate(mats):
            ax = axes[row, col]
            im = ax.imshow(mat, origin="lower", vmin=0, vmax=vmax, cmap="viridis", interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            if row == 0:
                ax.set_title("Archived truth" if col == 0 else NAMES[methods[col - 1]], fontsize=10)
            if col:
                scored = roles[methods[col - 1]]
                ax.set_xlabel(f"RMSE {scored['material_rmse']:.2f}   IoU {scored['material_support_iou']:.2f}", fontsize=7.5)
        fig.colorbar(im, ax=axes[row], shrink=.87, fraction=.03, pad=.01, label="Contrast")
    save_figure(fig, "state_rom_material_examples")


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    by_archive = read_all()
    summaries = {name: aggregate(list(saved["cases"].values())) for name, saved in by_archive.items()}
    r10 = by_archive["r10_twod"]["cases"].values()
    summaries["r10_weak"] = aggregate([c for c in r10 if c["truth_max"] <= 2.5 + 1e-6])
    summaries["r10_strong"] = aggregate([c for c in r10 if c["truth_max"] > 2.5 + 1e-6])
    summaries["main_96"] = aggregate([c for name in COHORTS for c in by_archive[name]["cases"].values()])
    compact = {}
    for name, saved in by_archive.items():
        compact[name] = {
            case_id: {key: value for key, value in case.items() if key != "roles"}
            for case_id, case in saved["cases"].items()
        }
    output = {
        "summary": summaries,
        "main_96_suite_cost_including_archive_and_setup": suite_cost(by_archive),
        "fixed_material_and_noise_diagnostics": compact,
    }
    (ROOT / "results" / "state_summary.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    compare_plot(by_archive, summaries)
    recon_plot(by_archive)


if __name__ == "__main__":
    main()
