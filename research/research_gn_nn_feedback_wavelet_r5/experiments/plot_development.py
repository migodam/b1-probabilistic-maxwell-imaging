"""Development-only R5 figures; never opens TRAIN truth or any TEST artifact.

Reference VAL truth is opened only after all five prediction files exist.
Run once on completion, rather than polling: python plot_development.py
Use --architecture-only while the finite training request is running.
Phase II: --phase noiseaware routes primary inputs to results/noiseaware and
models_noiseaware, and writes only figures/noiseaware. A separately labelled
comparison may read both phases' existing raw summaries/case identities.
Both phases remain DEVELOPMENT VAL, never fresh TEST.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import torch

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY.parents[1]))
from research.research_gn_nn_feedback_wavelet_r5.src.closed_loop import FeedbackCNN
from research.research_gn_nn_feedback_wavelet_r5.src.wavelets import WaveletPyramid2D

ARMS = ("pixel_closed_coupled", "haar_closed_coupled", "db2_closed_coupled",
        "haar_closed_no_fullfield", "haar_free_coupled")
LABEL = {"pixel_closed_coupled": "Pixel closed", "haar_closed_coupled": "Haar closed",
         "db2_closed_coupled": "db2 closed", "haar_closed_no_fullfield": "Haar no field loss",
         "haar_free_coupled": "Haar free control"}
COLORS = dict(zip(ARMS, ("#244b78", "#1b896c", "#b55a2f", "#9c6db2", "#a28725")))
SCOPE = "DEVELOPMENT VAL used for checkpoint/model selection"
EXPOSURE = "Image/mother absent TRAIN; classes/fonts may be seen. No fresh TEST."
OUTPUTS: list[str] = []
ACTIVE_PHASE = "strict"


def results_root(study):
    return study / "results" / "noiseaware" if ACTIVE_PHASE == "noiseaware" else study / "results"


def models_root(study):
    return study / "models_noiseaware" if ACTIVE_PHASE == "noiseaware" else study / "models"


def phase_label():
    return "Phase II: common noise-aware GN inputs" if ACTIVE_PHASE == "noiseaware" else "Phase I"


def write_table(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    OUTPUTS.append(path.name)


def read_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def truthy(value):
    return str(value).lower() in ("true", "1")


def raw_summary_table(rows, out):
    """Full raw cohort only; corrected subset is exported separately."""
    summary = []
    for arm in ARMS:
        for noisy in (False, True):
            selected = [row for row in rows if row["arm"] == arm and (float(row["noise"]) > 0) == noisy]
            if not selected:
                continue
            summary.append({"phase": phase_label(), "stage": "RAW_NN_FULL_COHORT",
                            "arm": arm, "noise_group": "noisy" if noisy else "clean",
                            "observations": len(selected), "mothers": len({r["parent_id"] for r in selected}),
                            **{key: float(np.mean([float(row[key]) for row in selected]))
                               for key in ("material_NRMSE", "material_RMSE", "IoU", "data_residual")},
                            "qualified": sum(truthy(row["qualified"]) for row in selected),
                            "joint": sum(truthy(row["joint"]) for row in selected)})
    write_table(out / "RAW_NN_FULL_COHORT_SUMMARY.csv", summary)
    return {"rows": len(rows), "observations_per_arm": {arm: sum(row["arm"] == arm for row in rows) for arm in ARMS},
            "table": "RAW_NN_FULL_COHORT_SUMMARY.csv", "corrected_rows_included": 0}


def title(fig, text):
    fig.suptitle(f"{phase_label()} | {text}\n{SCOPE}\n{EXPOSURE}", fontsize=13, y=.995)


def save(fig, out, name, pdf=False):
    path = out / (name + ".png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    OUTPUTS.append(path.name)
    if pdf:
        path = out / (name + ".pdf")
        fig.savefig(path, bbox_inches="tight", facecolor="white")
        OUTPUTS.append(path.name)
    plt.close(fig)


def architecture(study, out):
    models = {rep: FeedbackCNN(rep) for rep in ("pixel", "haar", "db2")}
    counts = {rep: model.parameter_counts() for rep, model in models.items()}
    pixel_width = models["pixel"].prior_network.cnn.layers[0].out_channels
    wavelet_width = models["haar"].prior_network.cnn.low.layers[0].out_channels
    metric_width = models["haar"].metric_network.cnn.low.layers[0].out_channels
    fig, ax = plt.subplots(figsize=(18, 11))
    ax.set_xlim(0, 18); ax.set_ylim(0, 9); ax.axis("off")
    title(fig, "Implemented architecture: conditional prior + bounded SPD current metric")

    def box(x, y, width, height, text, color="#edf2f7", size=11):
        ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=.10",
                                   linewidth=1.2, edgecolor="#486173", facecolor=color))
        ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=size)

    def arrow(start, end, color="#3d5363", style="-", connection="arc3"):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                                    linewidth=1.5, color=color, linestyle=style, connectionstyle=connection))

    anchor_label = "Common noise-aware GN anchor" if ACTIVE_PHASE == "noiseaware" else "Fixed strict-GN public anchor"
    box(.25, 7.6, 3.5, 1.15, anchor_label + "\n$m_{GN}, j_0, Q_0, Y$\n27 normalized feature maps", size=10.5)
    box(4.5, 7.6, 4.9, 1.15, r"""PriorCNN: material residual
$c=clip(m_{GN}+0.4\,CNN(x), .02,20)$
Closed: computed once; fixed across steps""", "#e5f3ea")
    arrow((3.85, 8.2), (4.4, 8.2))
    box(10.2, 7.6, 7.1, 1.15,
        r"""Fixed public weights / normalizations
$a=\|Y\|^2,\ b=\|j_0\|^2,\ \lambda_i=\rho\|Q_{0,i}\|^2/b$
$\rho=1$; no GT-dependent online setting""")
    box(.25, 4.65, 4.4, 1.8,
        "Exact shared-real material row solve\n$q=E+Gj$; minimize energy over $m_i$\n$e=j-mq$\nbox constraint $[.02,20]$", "#fff1d6")
    box(5.45, 4.65, 5.2, 1.8,
        r"""MetricCNN on 7 physical feature maps
$p=1+0.5\tanh(raw)$, $0.5\leq p\leq1.5$
$P=W^T diag(p)W$
Same spatial gains for all sources / real-imag""", "#e4eef9")
    box(11.45, 4.65, 5.85, 1.8,
        r"""Fixed-operator current update
$g=S^*(Sj-Y)/a+\rho[e-G^*(me)]/b$
$d=Pg$; per-case Armijo $j^+=j-\tau d$
Initial step + up to 6 halvings; rejected: keep $j$""", "#e4eef9")
    arrow((6.9, 7.5), (2.5, 6.55)); arrow((13.7, 7.5), (14.3, 6.55))
    arrow((4.75, 5.55), (5.35, 5.55)); arrow((10.75, 5.55), (11.35, 5.55))
    arrow((14.4, 4.55), (2.5, 4.55), connection="arc3,rad=-.08")
    ax.text(3.8, 3.75, "Repeat the SAME CNN weights\nRow solve again after accepted current update",
            ha="center", fontsize=10)
    box(.25, 1.35, 5.3, 1.9,
        f"Pixel / complete wavelet coordinates\nPixel prior width {pixel_width}; wavelet prior width {wavelet_width}\nMetric width {metric_width}; circular 3x3 / SiLU\nHaar or db2: 24→12→6; ALL 7 bands\n36 + 3×36 + 3×144 = 576 coefficients", "#edf2f7", 10.5)
    box(6.15, 1.35, 5.25, 1.9,
        "Coupled TRAIN loss on selected stages\nMaterial + 0.5 prior supervision\n+ actual nonlinear Maxwell field hinge\n+ constitutive and tail penalties\nRetained full LU + analytic adjoint/VJP", "#f7e7df", 10.5)
    box(12.0, 1.35, 5.3, 1.9,
        "Declared controls\nNo-field-loss: same inverse architecture\nFree cascade: recompute c from current m/Q\nFree control: NO shared fixed-energy guarantee\nNo LSTM/GRU; no inverse-layer LU", "#f6e9ee", 10.5)
    arrow((14.5, 4.55), (8.8, 3.4), color="#a45236", style="--")
    arrow((8.0, 3.35), (8.0, 4.55), color="#a45236", style="--")
    ax.text(13.4, 3.45, "TRAIN gradients reach prior and metric CNNs\nthrough actual physics and the inverse path",
            ha="center", fontsize=9, color="#a45236")
    ax.text(.3, .65, "Actual class parameter counts: " + "; ".join(
        f"{rep}: {counts[rep]['total']:,}" for rep in models) +
        "\nClosed fixed-energy descent is conditional; no global convergence or GT-recovery guarantee.", fontsize=11)
    fig.subplots_adjust(top=.86, bottom=.03)
    save(fig, out, "01_actual_architecture", pdf=True)
    return counts


def curves(study, out):
    histories = {arm: read_json(models_root(study) / arm / "history.json")
                 for arm in ARMS if (models_root(study) / arm / "history.json").exists()}
    if not histories:
        return False
    metrics = (("material_NRMSE", "Material NRMSE"),
               ("relative_public_residual", "Actual Maxwell residual, mean of evaluated stages"),
               ("selection_score", "Common composite selection score"))
    fig, axs = plt.subplots(3, 2, figsize=(17, 12), squeeze=False)
    title(fig, "TRAIN and DEVELOPMENT VAL histories (separate panels)")
    for column, split in enumerate(("TRAIN", "VAL")):
        for row, (metric, label) in enumerate(metrics):
            ax = axs[row, column]
            for arm, history in histories.items():
                epochs = [item["epoch"] for item in history]
                values = [item.get(split, {}).get(metric, np.nan) for item in history]
                if np.isfinite(values).any():
                    ax.plot(epochs, values, marker=".", color=COLORS[arm], label=LABEL[arm])
            ax.set_title(f"{split}: {label}"); ax.set_xlabel("Completed epoch")
            ax.grid(alpha=.2)
            if row:
                ax.set_yscale("log")
            if column == 0 and row == 0:
                ax.legend(fontsize=9)
    fig.text(.5, .01, "TRAIN field residual is unavailable for the no-field-loss arm; missing values are omitted, never shown as zero.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, .035, 1, .91))
    save(fig, out, "02_train_val_curves")
    return True


def metadata_label(record, obs, parent, noise):
    def seen(key):
        value = record.get(key)
        return "unknown" if value is None else ("seen" if value else "absent")
    identity = (f"digit {record.get('digit')}" if record.get("family") == "digits"
                else str(record.get("family", "unknown family")))
    font = str(record.get("font") or "not applicable")
    font = textwrap.fill(font, 22)
    return (f"{identity}; noise={noise:g}\n{obs}\n{parent}\n"
            f"image:{seen('seen_image')}; mother:{seen('seen_parent')}\n"
            f"class:{seen('seen_class')}; font:{seen('seen_font')}\n{font}")


def reconstruction_panel(reference, predictions, manifest, indices, out, name, headline):
    methods = ("GT", "GN", *ARMS)
    n = len(indices)
    height = 2.6 * n + 2.0
    fig, axs = plt.subplots(n, len(methods), figsize=(20, height), squeeze=False)
    title(fig, headline)
    truth = reference["m_gt"].reshape(-1, 576)
    amplitude = reference["amplitude"].reshape(-1)
    obs_ids = reference["obs_id"].astype(str)
    parents = reference["parent_id"].astype(str)
    for row, index in enumerate(indices):
        for column, method in enumerate(methods):
            ax = axs[row, column]
            value = (truth[index] if method == "GT" else reference["m_gn"][index]
                     if method == "GN" else predictions[method]["m_final"][index]).reshape(24, 24)
            image = ax.imshow(value, vmin=.02, vmax=.55, cmap="viridis", interpolation="nearest")
            lower, upper = int((value < .02).sum()), int((value > .55).sum())
            label = method if method in ("GT", "GN") else LABEL[method]
            extra = f"sat <.02:{lower}, >.55:{upper}"
            if method != "GT":
                nrmse = np.sqrt(np.mean((value.ravel() - truth[index]) ** 2)) / amplitude[index]
                extra = f"NRMSE={nrmse:.3f}; {extra}"
            if method in ARMS:
                residual = float(predictions[method]["data_residual"][index])
                qualified = bool(predictions[method]["qualified"][index])
                extra += f"\nr={residual:.2e}; {'QUALIFIED' if qualified else 'FAIL'}"
            ax.set_title(label + "\n" + extra, fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
        record = manifest.get(obs_ids[index], {})
        axs[row, 0].text(-.12, .5, metadata_label(record, obs_ids[index], parents[index],
                                               float(reference["noise"][index])),
                        transform=axs[row, 0].transAxes, ha="right", va="center", fontsize=8)
    fig.subplots_adjust(left=.23, right=.96, top=1 - 1.65 / height, bottom=.04,
                        hspace=.55, wspace=.35)
    color_ax = fig.add_axes((.975, .15, .007, .55))
    fig.colorbar(image, cax=color_ax, label="Material; fixed display [.02,.55]")
    save(fig, out, name)


def choose_prespecified(reference, manifest):
    chosen, missing = [], []
    ids = reference["obs_id"].astype(str)
    for group in ("digit2", "digit8", "islands", "primitive_composites"):
        for noisy in (False, True):
            matches = []
            for index, obs in enumerate(ids):
                record = manifest.get(obs, {})
                identity = (str(record.get("digit")) == group[-1] and record.get("family") == "digits"
                            if group.startswith("digit") else record.get("family") == group)
                if (identity and np.isclose(float(reference["amplitude"][index]), .2)
                        and (float(reference["noise"][index]) > 0) == noisy):
                    matches.append(index)
            if matches:
                chosen.append(matches[0])
            else:
                missing.append(f"{group}/{'noisy' if noisy else 'clean'}/A=.2")
    return chosen, missing


def band_energies(reference, predictions, out):
    bank = WaveletPyramid2D("haar")
    methods = {"GN": reference["m_gn"], **{arm: predictions[arm]["m_final"] for arm in ARMS}}
    energies, max_difference = {}, 0.
    rows = []
    for method, values in methods.items():
        error = torch.as_tensor(np.asarray(values - reference["m_gt"], np.float32).reshape(-1, 24, 24))
        bands = {key: value.numpy() / 576 for key, value in bank.energy_by_band(error).items()}
        total = sum(bands.values())
        mse = error.square().mean(dim=(-2, -1)).numpy()
        max_difference = max(max_difference, float(np.max(np.abs(total - mse))))
        energies[method] = bands
        for index, obs in enumerate(reference["obs_id"].astype(str)):
            rows.append({"method": method, "obs_id": obs, "noise": float(reference["noise"][index]),
                         **{key: float(value[index]) for key, value in bands.items()},
                         "sum_band_MSE": float(total[index]), "voxel_MSE": float(mse[index])})
    with (out / "HAAR_ERROR_BAND_ENERGIES.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    OUTPUTS.append("HAAR_ERROR_BAND_ENERGIES.csv")
    fig, axs = plt.subplots(1, 2, figsize=(17, 7))
    title(fig, "Absolute Haar error-band energies: summed band contributions equal voxel MSE")
    palette = plt.get_cmap("tab20").colors
    names = list(methods)
    for ax, noisy in zip(axs, (False, True)):
        mask = (reference["noise"].reshape(-1) > 0) == noisy
        bottom = np.zeros(len(names))
        for color_index, band in enumerate(bank.band_order):
            values = np.array([np.mean(energies[method][band][mask]) for method in names])
            ax.bar(np.arange(len(names)), values, bottom=bottom, label=band, color=palette[color_index])
            bottom += values
        ax.set_xticks(np.arange(len(names)), ["GN" if name == "GN" else LABEL[name] for name in names],
                      rotation=25, ha="right")
        ax.set_ylabel("Mean squared material error (absolute, not a percentage)")
        ax.set_title(f"{'Noisy' if noisy else 'Clean'} DEVELOPMENT VAL; n={mask.sum()}")
        ax.grid(axis="y", alpha=.2)
    axs[1].legend(ncol=4, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, .87)); save(fig, out, "04_absolute_haar_error_energies")
    return {"max_abs_band_sum_minus_voxel_MSE": max_difference, "normalization": "sum squares / 576"}


def cost_summary(study, rows, predictions, reference, out):
    fig, axs = plt.subplots(2, 2, figsize=(18, 11))
    title(fig, "Material quality, actual full-forward residual and charged computation")
    x = np.arange(len(ARMS)); width = .36
    for offset, noisy, color in ((-.5, False, "#365f86"), (.5, True, "#b77a43")):
        group = "Noisy" if noisy else "Clean"
        mean_quality, mean_residual, qualifications = [], [], []
        for arm in ARMS:
            selected = [r for r in rows if r["arm"] == arm and (float(r["noise"]) > 0) == noisy]
            mean_quality.append(np.mean([float(r["material_NRMSE"]) for r in selected]) if selected else np.nan)
            mask = (reference["noise"].reshape(-1) > 0) == noisy
            mean_residual.append(np.mean(predictions[arm]["data_residual"][mask]))
            qualifications.append(np.mean(predictions[arm]["qualified"][mask].astype(bool)))
        axs[0, 0].bar(x + offset * width, mean_quality, width, label=group, color=color)
        bars = axs[0, 1].bar(x + offset * width, mean_residual, width, label=group, color=color)
        for bar, fraction in zip(bars, qualifications):
            axs[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                           f"{fraction:.0%} Q", ha="center", va="bottom", fontsize=8)
    axs[0, 0].set_title("Material NRMSE; separate clean/noisy means")
    axs[0, 1].set_title("Actual nonlinear Maxwell residual; labels = qualified fraction")
    axs[0, 1].set_yscale("log"); axs[0, 1].axhline(.001, color="#365f86", linestyle=":", label="Clean tolerance .001")
    axs[0, 1].axhline(.003, color="#b77a43", linestyle=":", label="Noisy tolerance .003")
    charges = []
    for arm in ARMS:
        path = models_root(study) / arm / "physics_cost.json"
        if path.exists():
            receipt = read_json(path); count = receipt.get("counts", {})
            charges.append({"arm": arm, "LU": count.get("LU"),
                            "source_complex_RHS": count.get("source_complex_RHS"),
                            "adjoint_complex_RHS": count.get("adjoint_complex_RHS"),
                            "physics_wall_seconds": receipt.get("total_wall_seconds")})
        else:
            charges.append({"arm": arm, "LU": None, "status": "MISSING_RECEIPT"})
    lus = [item["LU"] if item.get("LU") is not None else np.nan for item in charges]
    axs[1, 0].bar(x, lus, color=[COLORS[arm] for arm in ARMS])
    axs[1, 0].set_title("Charged full LU: completed TRAIN + VAL loss/evaluation")
    axs[1, 0].set_ylabel("Actual retained full-wave LU count")
    for index, charge in enumerate(charges):
        if charge.get("LU") is not None:
            axs[1, 0].text(index, charge["LU"], f"{charge['LU']:,}\n"
                           f"src RHS {charge.get('source_complex_RHS', '?')}\n"
                           f"adj RHS {charge.get('adjoint_complex_RHS', '?')}",
                           ha="center", va="bottom", fontsize=8)
    if any(np.isfinite(lus)):
        axs[1, 0].set_ylim(0, np.nanmax(lus) * 1.30)
    means = [np.mean([float(r["forward_wall_seconds"]) for r in rows if r["arm"] == arm]) for arm in ARMS]
    axs[1, 1].bar(x, means, color=[COLORS[arm] for arm in ARMS])
    axs[1, 1].set_title("Reported prediction forward wall time per observation")
    axs[1, 1].set_ylabel("Seconds; CSV reporting scope")
    for ax in axs.flat:
        ax.set_xticks(x, [LABEL[arm] for arm in ARMS], rotation=20, ha="right")
        ax.grid(axis="y", alpha=.2)
    axs[0, 0].legend(); axs[0, 1].legend(fontsize=8)
    fig.text(.5, .01, "Inverse iterations use LU=0. Full nonlinear TRAIN/VAL loss/evaluation uses charged retained LU and adjoints. "
             "Prediction timings and training receipts have different scopes; no speed or recovery claim follows from cost alone.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .045, 1, .9)); save(fig, out, "05_quality_residual_actual_cost")
    return charges


def corrected_comparison(study, raw_rows, out):
    """Compare chosen online corrected outputs to exactly matched raw cases.

    Never uses best-K/oracle candidate scores. This 360-row subset is not the
    full 780-row raw cohort, and the two stages retain separate qualification.
    """
    path = results_root(study) / "CORRECTED_BY_CASE.csv"
    if not path.exists():
        return {"status": "NOT_AVAILABLE"}
    corrected = read_csv(path)
    raw = {(row["arm"], row["obs_id"]): row for row in raw_rows}
    if len(raw) != len(raw_rows):
        raise ValueError("duplicate raw case keys")
    pairs = []
    seen = set()
    for row in corrected:
        key = (row["arm"], row["obs_id"])
        if key in seen or key not in raw:
            raise ValueError(f"corrected/raw alignment failure: {key}")
        seen.add(key); pairs.append((raw[key], row))
    matched_summary = []
    for arm in ARMS:
        for noisy in (False, True):
            selected = [(r, c) for r, c in pairs if c["arm"] == arm and (float(c["noise"]) > 0) == noisy]
            if not selected:
                continue
            matched_summary.append({"phase": phase_label(), "stage": "CORRECTED_MATCHED_SUBSET",
                                    "arm": arm, "noise_group": "noisy" if noisy else "clean",
                                    "observations": len(selected),
                                    "mothers": len({c["parent_id"] for _, c in selected}),
                                    "raw_material_NRMSE": float(np.mean([float(r["material_NRMSE"]) for r, _ in selected])),
                                    "corrected_material_NRMSE": float(np.mean([float(c["NRMSE"]) for _, c in selected])),
                                    "raw_data_residual": float(np.mean([float(r["data_residual"]) for r, _ in selected])),
                                    "corrected_data_residual": float(np.mean([float(c["data_residual"]) for _, c in selected])),
                                    "raw_qualified": sum(truthy(r["qualified"]) for r, _ in selected),
                                    "corrected_qualified": sum(truthy(c["qualified"]) for _, c in selected),
                                    "raw_joint": sum(truthy(r["joint"]) for r, _ in selected),
                                    "corrected_joint": sum(truthy(c["joint"]) for _, c in selected)})
    write_table(out / "CORRECTED_MATCHED_SUBSET_SUMMARY.csv", matched_summary)
    fig, axs = plt.subplots(2, 2, figsize=(19, 12))
    title(fig, f"Raw NN versus separate physical corrector: exactly matched {len(pairs)} outputs")
    x = np.arange(len(ARMS)); width = .36
    counts = []
    for ax, noisy in zip(axs[0], (False, True)):
        raw_quality, corrected_quality = [], []
        for arm in ARMS:
            selected = [(r, c) for r, c in pairs if c["arm"] == arm and (float(c["noise"]) > 0) == noisy]
            raw_quality.append(np.mean([float(r["material_NRMSE"]) for r, _ in selected]))
            corrected_quality.append(np.mean([float(c["NRMSE"]) for _, c in selected]))
        ax.bar(x-width/2, raw_quality, width, color="#8c9aa5", label="Raw NN, matched subset")
        ax.bar(x+width/2, corrected_quality, width, color="#34778b", label="Separate physical corrector")
        ax.set_title(f"{'Noisy' if noisy else 'Clean'} matched DEVELOPMENT VAL: material NRMSE")
        ax.set_ylabel("Mean material NRMSE"); ax.legend(fontsize=9)
    raw_residual, corrected_residual = [], []
    incremental_wall, deployment_wall, incremental_lu, deployment_lu = [], [], [], []
    for arm in ARMS:
        selected = [(r, c) for r, c in pairs if c["arm"] == arm]
        if not selected:
            raise ValueError(f"corrected cohort missing arm {arm}")
        tolerances = [.003 if float(c["noise"]) > 0 else .001 for _, c in selected]
        raw_residual.append(np.mean([float(r["data_residual"])/tol for (r, _), tol in zip(selected, tolerances)]))
        corrected_residual.append(np.mean([float(c["data_residual"])/tol for (_, c), tol in zip(selected, tolerances)]))
        count = {"arm": arm, "cases": len(selected),
                 "raw_qualified": sum(truthy(r["qualified"]) for r, _ in selected),
                 "corrected_qualified": sum(truthy(c["qualified"]) for _, c in selected),
                 "raw_joint": sum(truthy(r["joint"]) for r, _ in selected),
                 "corrected_joint": sum(truthy(c["joint"]) for _, c in selected)}
        counts.append(count)
        incremental_wall.append(np.mean([float(c["wall"]) for _, c in selected]))
        deployment_wall.append(np.mean([float(c["deployment_wall"]) for _, c in selected]))
        incremental_lu.append(np.mean([float(c["LU"]) for _, c in selected]))
        deployment_lu.append(np.mean([float(c["deployment_LU"]) for _, c in selected]))
    ax = axs[1, 0]
    for sign, values, color, prefix in ((-1, raw_residual, "#8c9aa5", "raw"),
                                        (1, corrected_residual, "#34778b", "corrected")):
        bars = ax.bar(x+sign*width/2, values, width, color=color,
                      label="Raw NN" if prefix == "raw" else "Separate physical corrector")
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height(),
                    f"Q {count[prefix+'_qualified']}/{count['cases']}\nJ {count[prefix+'_joint']}/{count['cases']}",
                    ha="center", va="bottom", fontsize=8)
    ax.axhline(1, linestyle=":", color="#a14d35", label="Qualification tolerance")
    ax.set_yscale("log"); ax.set_ylim(bottom=min(corrected_residual)*.3, top=max(raw_residual)*2.5)
    ax.set_title("Mean actual residual / per-case tolerance; Q=qualified, J=recorded joint")
    ax.set_ylabel("Residual / tolerance (per-case normalization)"); ax.legend(fontsize=8)
    ax = axs[1, 1]
    for sign, values, lus, color, label in ((-1, incremental_wall, incremental_lu, "#c59a59", "Additional corrector wall"),
                                          (1, deployment_wall, deployment_lu, "#697b91", "Reported full deployment wall")):
        bars = ax.bar(x+sign*width/2, values, width, color=color, label=label)
        for bar, lu in zip(bars, lus):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height(), f"LU {lu:.1f}",
                    ha="center", va="bottom", fontsize=8)
    ax.set_title("Charged corrected-stage timing; labels = mean full-LU counts")
    ax.set_ylabel("Mean reported seconds per case"); ax.legend(fontsize=9)
    ax.set_ylim(0, max(deployment_wall)*1.2)
    for ax in axs.flat:
        ax.set_xticks(x, [LABEL[arm] for arm in ARMS], rotation=20, ha="right")
        ax.grid(axis="y", alpha=.2)
    fig.text(.5, .01, "Same selected online outputs; no best-K/oracle replacement. Physical qualification and material quality are distinct. "
             "Corrector qualification is not raw NN qualification. Full deployment costs include the reused GN anchor as charged by the result CSV.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .90)); save(fig, out, "07_raw_vs_corrected_matched_subset")
    return {"status": "PLOTTED_MATCHED_SUBSET", "counts": counts,
            "rows": len(pairs), "table": "CORRECTED_MATCHED_SUBSET_SUMMARY.csv",
            "matched_obs_ids_by_arm": {arm: [c["obs_id"] for _, c in pairs if c["arm"] == arm] for arm in ARMS},
            "selection": "recorded online chosen outputs only; no best-K/oracle fields used"}


def input_mode_comparison(study, out):
    """Observed whole-pipeline comparison, not a causal wavelet/loss ablation."""
    if ACTIVE_PHASE != "noiseaware":
        return {"status": "NOT_REQUESTED_FOR_PHASE_I"}
    roots = {"Phase I / strict GN": study / "results",
             "Phase II / common noise-aware GN": study / "results" / "noiseaware"}
    paths = [root / name for root in roots.values()
             for name in ("DEVELOPMENT_SUMMARY.json", "RESULTS_BY_CASE.csv")]
    if not all(path.exists() for path in paths):
        return {"status": "WAITING_FOR_BOTH_RAW_SUMMARIES", "missing": [str(path.relative_to(study)) for path in paths if not path.exists()]}
    summaries, cases = {}, {}
    for label, root in roots.items():
        document = read_json(root / "DEVELOPMENT_SUMMARY.json")
        summaries[label] = document["rows"] if isinstance(document, dict) else document
        cases[label] = read_csv(root / "RESULTS_BY_CASE.csv")
    def case_keys(rows):
        return {(r["arm"], r["obs_id"], r["parent_id"], float(r["noise"]), float(r["amplitude"])) for r in rows}
    phase_names = list(roots)
    if any(len(case_keys(rows)) != len(rows) for rows in cases.values()):
        return {"status": "DUPLICATE_RAW_CASE_IDENTITIES_NO_COMPARISON_PLOTTED"}
    if case_keys(cases[phase_names[0]]) != case_keys(cases[phase_names[1]]):
        return {"status": "COHORT_IDENTITIES_DIFFER_NO_COMPARISON_PLOTTED"}
    if any(len(rows) != 780 for rows in cases.values()):
        return {"status": "INCOMPLETE_RAW_COHORT_NO_COMPARISON_PLOTTED"}
    tables, indexed = [], {}
    for label, rows in summaries.items():
        indexed[label] = {(row["arm"], float(row["noise"]) > 0): row for row in rows}
        if any((arm, noisy) not in indexed[label] for arm in ARMS for noisy in (False, True)):
            return {"status": "INCOMPLETE_RAW_SUMMARY_NO_COMPARISON_PLOTTED"}
        for row in rows:
            tables.append({"phase_and_anchor_mode": label, "stage": "RAW_NN_FULL_COHORT",
                           **{key: row.get(key) for key in ("arm", "noise", "observations", "mothers",
                                                         "NRMSE", "RMSE", "IoU", "data_residual", "qualified", "joint")}})
    write_table(out / "INPUT_MODE_RAW_PIPELINE_COMPARISON.csv", tables)
    fig, axs = plt.subplots(2, 2, figsize=(18, 11))
    fig.suptitle("Phase I versus Phase II | observed GN-input pipeline comparison\n" + SCOPE + "\n" + EXPOSURE,
                 fontsize=13, y=.995)
    x, width = np.arange(len(ARMS)), .36
    for column, noisy in enumerate((False, True)):
        for row_index, metric in enumerate(("NRMSE", "data_residual")):
            ax = axs[row_index, column]
            for phase_index, (label, color) in enumerate(zip(phase_names, ("#8c9aa5", "#34778b"))):
                values = [float(indexed[label][(arm, noisy)][metric]) for arm in ARMS]
                bars = ax.bar(x+(phase_index-.5)*width, values, width, color=color, label=label)
                if metric == "data_residual":
                    for bar, arm in zip(bars, ARMS):
                        record = indexed[label][(arm, noisy)]
                        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height(),
                                f"Q {int(record['qualified'])}/{int(record['observations'])}",
                                ha="center", va="bottom", fontsize=8)
            ax.set_title(f"{'Noisy' if noisy else 'Clean'} raw DEVELOPMENT VAL: " +
                         ("material NRMSE" if metric == "NRMSE" else "actual nonlinear Maxwell residual"))
            ax.set_xticks(x, [LABEL[arm] for arm in ARMS], rotation=20, ha="right")
            ax.grid(axis="y", alpha=.2)
            if metric == "data_residual":
                ax.set_yscale("log")
                ax.axhline(.003 if noisy else .001, linestyle=":", color="#a14d35")
                ax.margins(y=.2)
            if row_index == 0:
                ax.legend(fontsize=9)
    fig.text(.5, .01, "Same declared cohort, architecture/loss and caps; GN inputs and separate from-scratch fits differ. "
             "This observed pipeline difference does not isolate causal wavelet benefit or a loss-only effect. Corrected subsets are excluded.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .90)); save(fig, out, "08_gn_input_mode_raw_pipeline_comparison")
    return {"status": "PLOTTED_OBSERVED_PIPELINE_DIFFERENCE", "raw_rows_per_phase": 780,
            "identity_match": True, "corrected_rows_included": 0,
            "attribution": "GN inputs and independently trained weights differ; no isolated causal wavelet/loss claim",
            "sources": [str(path.relative_to(study)) for path in paths]}


def no_field_material_physics_gap(raw_rows, out):
    """Existing VAL material scores versus measured nonlinear forward residual.

    These two raw arms share the selected phase's GN input protocol. Keep the
    no-full-field ablation distinct from the extra noise-ball mechanism arm.
    """
    if ACTIVE_PHASE != "noiseaware":
        return {"status": "NOT_REQUESTED_FOR_PHASE_I"}
    coupled, ablated = "haar_closed_coupled", "haar_closed_no_fullfield"
    indexed = {arm: {row["obs_id"]: row for row in raw_rows if row["arm"] == arm}
               for arm in (coupled, ablated)}
    if set(indexed[coupled]) != set(indexed[ablated]):
        raise ValueError("no-field ablation case identities do not match coupled Haar")
    matched = []
    for obs, base in indexed[coupled].items():
        control = indexed[ablated][obs]
        if (base["parent_id"] != control["parent_id"] or float(base["noise"]) != float(control["noise"])
                or float(base["amplitude"]) != float(control["amplitude"])):
            raise ValueError(f"no-field ablation metadata mismatch: {obs}")
        matched.append({"phase": phase_label(), "stage": "MATCHED_RAW_NN_ONLY", "obs_id": obs,
                        "parent_id": base["parent_id"], "noise": float(base["noise"]),
                        "amplitude": float(base["amplitude"]),
                        "coupled_material_NRMSE": float(base["material_NRMSE"]),
                        "no_field_material_NRMSE": float(control["material_NRMSE"]),
                        "coupled_data_residual": float(base["data_residual"]),
                        "no_field_data_residual": float(control["data_residual"]),
                        "coupled_qualified": truthy(base["qualified"]),
                        "no_field_qualified": truthy(control["qualified"]),
                        "delta_material_NRMSE_no_field_minus_coupled": float(control["material_NRMSE"])-float(base["material_NRMSE"]),
                        "delta_data_residual_no_field_minus_coupled": float(control["data_residual"])-float(base["data_residual"])})
    write_table(out / "NO_FIELD_MATERIAL_PHYSICS_GAP.csv", matched)
    fig, axs = plt.subplots(1, 2, figsize=(17, 8))
    title(fig, "No-field-loss ablation: VAL material error versus actual physical discrepancy")
    summaries = []
    for ax, noisy in zip(axs, (False, True)):
        group = [row for row in matched if (row["noise"] > 0) == noisy]
        tolerance = .003 if noisy else .001
        summary = {"noise_group": "noisy" if noisy else "clean", "observations": len(group)}
        for prefix, color, marker, label in (("coupled", "#246787", "o", "Haar coupled"),
                                             ("no_field", "#b55438", "D", "Haar no field loss")):
            material = np.array([row[prefix+"_material_NRMSE"] for row in group])
            residual = np.array([row[prefix+"_data_residual"] for row in group])
            ax.scatter(material, residual, s=23, alpha=.5, color=color, marker=marker,
                       label=f"{label}, n={len(group)}")
            ax.scatter([material.mean()], [residual.mean()], s=155, color=color, marker=marker,
                       edgecolor="black", linewidth=1.2)
            summary[prefix+"_mean_material_NRMSE"] = float(material.mean())
            summary[prefix+"_mean_data_residual"] = float(residual.mean())
            summary[prefix+"_qualified"] = sum(row[prefix+"_qualified"] for row in group)
        ax.axhline(tolerance, color="#707070", linestyle=":", label=f"Physical tolerance {tolerance:g}")
        ax.set_title(f"{'Noisy' if noisy else 'Clean'} raw DEVELOPMENT VAL; large symbols = means")
        ax.set_xlabel("Material NRMSE against existing VAL truth")
        ax.set_ylabel("Actual nonlinear Maxwell relative data residual")
        ax.set_yscale("log"); ax.grid(alpha=.2); ax.legend(fontsize=9)
        summaries.append(summary)
    fig.text(.5, .02, "Two separate outcomes: material closeness and original-data qualification. "
             "Both are raw inverse outputs using common noise-aware GN inputs; corrected results and the extra noise-ball CNN are excluded.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .05, 1, .86)); save(fig, out, "09_no_field_material_vs_physical_gap")
    return {"status": "PLOTTED_MATCHED_RAW_HAAR_ABLATION", "rows": len(matched), "groups": summaries,
            "extra_noise_ball_arm_included": False, "corrected_outputs_included": False,
            "table": "NO_FIELD_MATERIAL_PHYSICS_GAP.csv"}


def depth(study, out):
    path = results_root(study) / "DEPTH_DIAGNOSTICS.csv"
    if not path.exists():
        return False
    rows = read_csv(path)
    metrics = ("material_NRMSE", "data_residual", "constitutive_defect", "current_update", "energy")
    groups = [("", "Clean + noisy", rows)]
    if rows and "noise" in rows[0]:
        groups += [("_clean", "Clean", [row for row in rows if float(row["noise"]) == 0]),
                   ("_noisy", "Noisy", [row for row in rows if float(row["noise"]) > 0])]
    scope = {"source": str(path.relative_to(study)), "phase": phase_label(), "rows": len(rows), "groups": {}}
    for suffix, group, subset in groups:
        if not subset:
            continue
        counts = {arm: len({item["obs_id"] for item in subset if item["arm"] == arm}) for arm in ARMS}
        sizes = sorted(set(counts.values()))
        size_label = str(sizes[0]) if len(sizes) == 1 else f"{min(sizes)}–{max(sizes)}"
        scope["groups"][group] = {"case_counts": counts,
                                  "obs_ids": sorted({item["obs_id"] for item in subset}),
                                  "depths": sorted({int(item["depth"]) for item in subset})}
        fig, axs = plt.subplots(2, 3, figsize=(18, 10))
        title(fig, f"{group} depth diagnostics on declared subset: {size_label} observations per arm")
        for ax, metric in zip(axs.flat, metrics):
            for arm in ARMS:
                selected = [item for item in subset if item["arm"] == arm]
                depths = sorted({int(item["depth"]) for item in selected})
                values = [np.mean([float(item[metric]) for item in selected if int(item["depth"]) == level])
                          for level in depths]
                if depths:
                    ax.plot(depths, values, marker="o", color=COLORS[arm], label=LABEL[arm],
                            linestyle="--" if "free" in arm else "-")
            ax.set_title("Mean " + metric.replace("_", " ")); ax.set_xlabel("Iteration depth"); ax.grid(alpha=.2)
            if metric in ("data_residual", "constitutive_defect", "current_update"):
                ax.set_yscale("symlog", linthresh=1e-8)
        axs[1, 2].axis("off")
        axs[1, 2].text(0, .9, "Declared diagnostic subset, not all 156 VAL cases.\n\nEnergy descent is conditional on the declared\nfixed prior and accepted updates.\n\nFree cascade refreshes c; its energies are not\na shared fixed-energy convergence certificate.\n\nSmaller state changes do not establish\ncloser material or global convergence.", fontsize=11, va="top")
        axs[0, 0].legend(fontsize=8)
        fig.tight_layout(rect=(0, 0, 1, .9)); save(fig, out, "06_depth_diagnostics" + suffix)
    (out / "DEPTH_PLOT_SCOPE.json").write_text(json.dumps(scope, indent=2) + "\n")
    OUTPUTS.append("DEPTH_PLOT_SCOPE.json")
    return True


def main():
    global ACTIVE_PHASE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=STUDY)
    parser.add_argument("--phase", choices=("strict", "noiseaware"), default="strict",
                        help="strict=immutable Phase I roots; noiseaware=separate Phase II result/model/figure roots")
    parser.add_argument("--architecture-only", action="store_true")
    parser.add_argument("--require-results", action="store_true")
    parser.add_argument("--no-gallery", action="store_true")
    parser.add_argument("--depth-only", action="store_true", help="Render only the selected phase's optional depth CSV; preserve other figures")
    args = parser.parse_args(); ACTIVE_PHASE = args.phase
    study = args.study.resolve()
    result_dir, model_dir = results_root(study), models_root(study)
    out = study / "figures" / "noiseaware" if ACTIVE_PHASE == "noiseaware" else study / "figures"
    phase_id = "PHASE_II_COMMON_NOISEAWARE_GN_INPUTS" if ACTIVE_PHASE == "noiseaware" else "PHASE_I_STRICT_GN_INPUTS"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white"})
    if args.depth_only:
        plotted = depth(study, out)
        path = out / "FIGURES.json"
        receipt = read_json(path) if path.exists() else {"phase": phase_id, "scope": SCOPE}
        receipt.update(depth_plotted=plotted, depth_status="PLOTTED" if plotted else "NOT_AVAILABLE")
        receipt["files"] = list(dict.fromkeys(receipt.get("files", []) + OUTPUTS))
        path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        print(json.dumps({"status": receipt["depth_status"], "files": OUTPUTS})); return
    receipt = {"scope": SCOPE, "exposure": EXPOSURE, "parameter_counts": architecture(study, out),
               "phase": phase_id, "Phase_II_inputs_used": ACTIVE_PHASE == "noiseaware",
               "results_root": str(result_dir.relative_to(study)),
               "models_root": str(model_dir.relative_to(study)),
               "figures_root": str(out.relative_to(study)),
               "input_mode_attribution": "anchor conditioning and independent from-scratch fits; not isolated wavelet/loss causality",
               "TRAIN_truth_opened": False, "TEST_opened": False, "VAL_truth_opened": False}
    if args.architecture_only:
        receipt["status"] = "ARCHITECTURE_READY_RESULTS_NOT_READ"
    else:
        receipt["histories_plotted"] = curves(study, out)
        prediction_paths = {arm: result_dir / "VAL_PREDICTIONS" / (arm + ".npz") for arm in ARMS}
        missing = [str(path.relative_to(study)) for path in prediction_paths.values() if not path.exists()]
        required = ("VAL_REFERENCE.npz", "VAL_MANIFEST.json", "RESULTS_BY_CASE.csv")
        missing += [str((result_dir / name).relative_to(study)) for name in required if not (result_dir / name).exists()]
        if missing:
            receipt.update(status="WAITING_FOR_COMPLETE_PREDICTIONS", missing=missing)
        else:
            # Predictions are checked before opening the reference containing VAL GT.
            predictions = {arm: dict(np.load(path, allow_pickle=False)) for arm, path in prediction_paths.items()}
            reference = dict(np.load(result_dir / "VAL_REFERENCE.npz", allow_pickle=False))
            receipt["VAL_truth_opened"] = True
            n = len(reference["obs_id"])
            if n != 156:
                raise ValueError(f"Expected ordered156 DEVELOPMENT VAL rows, got {n}")
            for arm, values in predictions.items():
                if values["m_final"].shape != (n, 576):
                    raise ValueError(f"Unexpected prediction shape for {arm}")
                if "obs_id" in values and not np.array_equal(values["obs_id"].astype(str), reference["obs_id"].astype(str)):
                    raise ValueError(f"Prediction order mismatch for {arm}")
            records = read_json(result_dir / "VAL_MANIFEST.json")
            records = records.get("records", records) if isinstance(records, dict) else records
            manifest = {item["obs_id"]: item for item in records}
            selected, absent = choose_prespecified(reference, manifest)
            receipt.update(prespecified_indices=selected, prespecified_missing=absent,
                           prespecified_rule="first manifest/reference ordered match, A=.2, digit2/digit8/islands/primitive_composites, clean then noisy; no outcome selection")
            if selected:
                reconstruction_panel(reference, predictions, manifest, selected, out,
                                     "03_prespecified_val_reconstructions", "Prespecified A=.2 reconstruction examples; all five arms")
            if not args.no_gallery:
                for start in range(0, n, 6):
                    reconstruction_panel(reference, predictions, manifest, list(range(start, min(start + 6, n))), out,
                                         f"gallery_{start//6+1:02d}", f"Full ordered DEVELOPMENT VAL gallery: rows {start+1}–{min(start+6,n)} of {n}; failures retained")
            receipt["haar_energy"] = band_energies(reference, predictions, out)
            rows = read_csv(result_dir / "RESULTS_BY_CASE.csv")
            receipt["raw_full_cohort_table"] = raw_summary_table(rows, out)
            receipt["training_physics_receipts"] = cost_summary(study, rows, predictions, reference, out)
            receipt["raw_vs_corrected"] = corrected_comparison(study, rows, out)
            receipt["gn_input_mode_comparison"] = input_mode_comparison(study, out)
            receipt["no_field_material_physical_gap"] = no_field_material_physics_gap(rows, out)
            receipt["depth_plotted"] = depth(study, out)
            receipt["history_completed_epochs"] = {arm: len(read_json(model_dir/arm/"history.json")) for arm in ARMS}
            input_paths = list(prediction_paths.values()) + [result_dir/name for name in required]
            input_paths += [model_dir/arm/"history.json" for arm in ARMS]
            input_paths += [model_dir/arm/"physics_cost.json" for arm in ARMS]
            for name in ("CORRECTED_BY_CASE.csv", "DEPTH_DIAGNOSTICS.csv"):
                if (result_dir/name).exists():
                    input_paths.append(result_dir/name)
            input_paths += [study/path for path in receipt["gn_input_mode_comparison"].get("sources", [])]
            input_paths = list(dict.fromkeys(input_paths))
            receipt["input_provenance"] = [{"path": str(path.relative_to(study)), "bytes": path.stat().st_size,
                                             "mtime_ns": path.stat().st_mtime_ns} for path in input_paths]
            receipt["status"] = "COMPLETE_DEVELOPMENT_FIGURES"
    receipt["files"] = OUTPUTS
    (out / "FIGURES.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": receipt["status"], "files": len(OUTPUTS), "missing": receipt.get("missing", [])}))
    if args.require_results and receipt["status"] != "COMPLETE_DEVELOPMENT_FIGURES":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
