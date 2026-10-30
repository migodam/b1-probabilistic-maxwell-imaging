"""Post-seal Stage 4 figures, seedless truth release, and paired descriptions.

Requires a closed prediction seal and its already-written offline scorer output.
It never imports a solver or performs a Maxwell solve. Failed images remain in
every panel. The release mapping is built from an explicit field whitelist.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np

HERE = Path(__file__).resolve().parents[1]
METHODS = ("native_trf", "reciprocal_trf", "data_lbfgs", "gd_all",
           "feedback_beta:0.2", "feedback_qbeta")
LABELS = ("Native", "Reciprocal", "Data L-BFGS", "GD all", "β = 0.2", "qβ")
COLORS = {"native_trf": "#4c78a8", "reciprocal_trf": "#72b7b2",
          "data_lbfgs": "#9c755f", "gd_all": "#777777",
          "feedback_beta:0.2": "#f28e2b", "feedback_qbeta": "#b279a2"}
PUBLIC_KEYS = ("observation_id", "group_id", "shape_slot", "shape_family",
               "amplitude", "background", "noise", "noise_realization",
               "source_layout", "truth_file")
ARRAY_KEYS = ("m", "alpha", "mask", "Y_clean")


def read(path: Path):
    return json.loads(path.read_text())


def write(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def verify_closed(seal: dict, scores: dict) -> None:
    if (seal.get("status") != "ALL_PLANNED_PREDICTIONS_CLOSED" or
            seal.get("schema") != "fbg-r2-prediction-seal-v2" or
            scores.get("schema") != "fbg-r2-blind-scores-v2" or
            not scores.get("truth_disclosed_only_after_seal") or
            seal.get("planned_observations") != 22 or
            seal.get("planned_predictions") != 132 or
            len(scores.get("rows", [])) != 132):
        raise ValueError("Stage 4 closure and 132 scored method observations required")
    if ({r["observation_id"] for r in seal["rows"]} !=
            {r["observation_id"] for r in scores["rows"]}):
        raise ValueError("Scored IDs differ from sealed IDs")


def export_sanitized(private_dir: Path, release_dir: Path,
                     seal: dict) -> dict:
    if release_dir.exists() and any(release_dir.iterdir()):
        raise FileExistsError("Preserve existing truth-after-freeze release")
    private = read(private_dir / "mapping_PRIVATE.json")
    records = private["records"]
    if (len(records) != 22 or {r["observation_id"] for r in records} !=
            {r["observation_id"] for r in seal["rows"]}):
        raise ValueError("Private mapping and closed public IDs disagree")
    release_dir.mkdir(parents=True, exist_ok=True)
    sanitized = []
    for r in records:
        clean = {key: r[key] for key in PUBLIC_KEYS}
        if any("seed" in key.lower() or "entropy" in key.lower() for key in clean):
            raise RuntimeError("Metadata whitelist violated")
        with np.load(private_dir / r["truth_file"], allow_pickle=False) as f:
            if set(f.files) != set(ARRAY_KEYS):
                raise ValueError("Private truth NPZ contains unexpected arrays")
            arrays = {key: f[key] for key in ARRAY_KEYS}
        np.savez_compressed(release_dir / r["truth_file"], **arrays)
        sanitized.append(clean)
    output = {"schema": "fbg-r2-sanitized-truth-after-freeze-v1",
              "truth_policy": "Released only after all 132 predictions were closed and scored",
              "records": sanitized}
    write(release_dir / "mapping_PRIVATE.json", output)
    return output


def sealed_images(seal: dict) -> dict:
    out = {}
    for row in seal["rows"]:
        ident = row["observation_id"]
        for method in row["methods"]:
            selected = method.get("online_selected_start_id")
            start = next((s for s in method["starts"] if s["start_id"] == selected), None)
            m = (start or {}).get("prediction") or {}
            out[(ident, method["method"])] = m.get("m")
    return out


def map_scores(scores: dict) -> dict:
    return {(r["observation_id"], r["method"]): r for r in scores["rows"]}


def render_grid(name: str, selected: list[dict], release_dir: Path,
                images: dict, scored: dict, figures: Path) -> None:
    rows = len(selected)
    fig, axes = plt.subplots(rows, 7, figsize=(23, max(7, rows*2.55)),
                             constrained_layout=True, squeeze=False)
    norm = PowerNorm(gamma=.46, vmin=.02, vmax=20.)
    image = None
    for i, rec in enumerate(selected):
        ident = rec["observation_id"]
        with np.load(release_dir / rec["truth_file"], allow_pickle=False) as f:
            truth = f["m"].reshape(24, 24)
        image = axes[i, 0].imshow(truth, cmap="magma", norm=norm,
                                  origin="upper", interpolation="nearest")
        title = (f"Truth · {rec['shape_family']} · A={rec['amplitude']:g}\n"
                 f"{rec['noise']} {rec['noise_realization'] if rec['noise'] == 'noisy' else ''} · {rec['source_layout']}")
        axes[i, 0].set_title(title, fontsize=7.8)
        axes[i, 0].set_ylabel(f"Object {rec['shape_slot']+1}", fontsize=9)
        for j, (method, label) in enumerate(zip(METHODS, LABELS), start=1):
            ax = axes[i, j]
            row = scored.get((ident, method))
            online = (row or {}).get("online") or {}
            material = images.get((ident, method))
            if material is not None and len(material) == 576 and np.isfinite(material).all():
                ax.imshow(np.asarray(material).reshape(24, 24), cmap="magma", norm=norm,
                          origin="upper", interpolation="nearest")
                gate = online.get("joint_gate", "FAIL")
                nrmse = online.get("NRMSE")
                resid = online.get("relative_data_residual")
                ax.set_title(f"{label} · {gate}\nN={nrmse:.3f}  R={resid:.1e}" if nrmse is not None and resid is not None
                             else f"{label} · {gate}\nN/R missing",
                             fontsize=7.6, color="#157a45" if gate == "PASS" else "#8b2727")
            else:
                ax.set_facecolor("#eeeeee")
                ax.text(.5, .5, "NO FINITE\nENDPOINT", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{label} · FAILED/NOT_RUN", fontsize=7.6)
        for ax in axes[i]:
            ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(image, ax=axes, shrink=.68, pad=.008)
    cb.set_label("Material m · shared [0.02, 20], power γ=0.46")
    fig.suptitle(f"Stage 4 fresh blind imaging after seal · {name.replace('_', ' ')}\n"
                 "Full nonlinear 24×24 scalar-TM material; every failure retained",
                 fontsize=14, fontweight="bold")
    figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / (name + ".png"), dpi=135)
    plt.close(fig)


def paired_analysis(scores: dict) -> dict:
    rows = scores["rows"]
    by_id = {(r["observation_id"], r["method"]): r for r in rows}
    clean = [r for r in rows if r["noise"] == "clean" and r["method"] == METHODS[0]]
    clean.sort(key=lambda r: r["shape_slot"])
    case_rows = []
    for base in clean:
        ident = base["observation_id"]
        methods = {m: by_id[(ident, m)]["online"] for m in METHODS}
        case_rows.append({"shape_slot": base["shape_slot"],
                          "shape_family_OFFLINE": base["shape_family_OFFLINE"],
                          "amplitude": base["amplitude"],
                          "group_id": base["group_id"],
                          "clean_observation_id": ident,
                          "online_NRMSE": {m: methods[m]["NRMSE"] for m in METHODS},
                          "online_joint_gate": {m: methods[m]["joint_gate"] for m in METHODS},
                          "delta_NRMSE_beta_minus_GD": methods["feedback_beta:0.2"]["NRMSE"] - methods["gd_all"]["NRMSE"],
                          "delta_NRMSE_qbeta_minus_GD": methods["feedback_qbeta"]["NRMSE"] - methods["gd_all"]["NRMSE"]})
    # Equal-object cluster means include all available clean/noisy variants.
    clusters = []
    for group_id in sorted({r["group_id"] for r in rows}):
        group = [r for r in rows if r["group_id"] == group_id]
        method_mean = {m: float(np.mean([r["online"]["NRMSE"] for r in group if r["method"] == m]))
                       for m in METHODS}
        clusters.append({"group_id": group_id,
                         "shape_slot": group[0]["shape_slot"],
                         "variants_per_method": len([r for r in group if r["method"] == METHODS[0]]),
                         "mean_NRMSE": method_mean,
                         "delta_beta_minus_GD": method_mean["feedback_beta:0.2"] - method_mean["gd_all"]})
    delta = np.array([c["delta_beta_minus_GD"] for c in clusters])
    rng = np.random.default_rng(20261030)
    boot = [float(np.mean(rng.choice(delta, size=len(delta), replace=True))) for _ in range(1000)]
    high = [r for r in case_rows if r["amplitude"] == 12.]
    return {"schema": "fbg-r2-stage4-paired-descriptive-v1",
            "status": "OFFLINE_AFTER_SEAL", "independent_objects": 14,
            "observation_variants": 22,
            "best_at_K_scope": "K=1 for every method; offline best@K is the same endpoint as online top-1, not multibranch evidence",
            "clean_per_shape": case_rows, "equal_object_variant_clusters": clusters,
            "high_A12_clean_delta_beta_minus_GD_mean": float(np.mean([r["delta_NRMSE_beta_minus_GD"] for r in high])),
            "equal_object_mean_delta_beta_minus_GD": float(np.mean(delta)),
            "equal_object_bootstrap_1000_percentile_95_descriptive": [float(x) for x in np.quantile(boot, [.025, .975])],
            "bootstrap_scope": "resamples 14 equal-weight observed object clusters; descriptive stability only, not population inference or a p-value"}


def cost_quality(scores: dict, figures: Path) -> None:
    rows = scores["rows"]
    categories = [
        ("High A=12 clean · 12 objects", [r for r in rows if r["noise"] == "clean" and r["amplitude"] == 12.]),
        ("Contrast controls clean · 2 objects", [r for r in rows if r["noise"] == "clean" and r["amplitude"] != 12.]),
        ("Noisy · 8 observations / 4 objects", [r for r in rows if r["noise"] == "noisy"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for ax, (title, subset) in zip(axes, categories):
        for method in METHODS:
            group = [r for r in subset if r["method"] == method]
            xs = [r["measured_costs_across_starts"]["wall_seconds_sum_known"] for r in group]
            ys = [r["online"]["NRMSE"] if r.get("online") else np.nan for r in group]
            passes = [r["online"]["joint_gate"] == "PASS" if r.get("online") else False for r in group]
            for x, y, passed in zip(xs, ys, passes):
                if x is not None and np.isfinite(y):
                    ax.scatter(x, y, marker="o" if passed else "x", s=55,
                               color=COLORS[method], alpha=.8)
        ax.axhline(.15, color="#2f8c61", ls="--", lw=1)
        ax.set_title(title)
        ax.set_yscale("log")
        ax.set_xlabel("Charged task wall (s), per method/observation")
        ax.set_ylabel("Material NRMSE · ○ joint PASS, × FAIL")
        ax.grid(alpha=.2)
    handles = [plt.Line2D([], [], marker="o", lw=0, color=COLORS[m], label=l)
               for m, l in zip(METHODS, LABELS)]
    axes[0].legend(handles=handles, fontsize=8, ncol=2)
    fig.suptitle("Fresh Stage 4 · physical image quality versus observed solver work\n"
                 "Methods are paired within observations; noisy variants share four objects",
                 fontsize=14, fontweight="bold")
    fig.savefig(figures / "fresh_cost_quality.png", dpi=160)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seal", type=Path, default=HERE / "audit/STAGE4_PREDICTIONS_CLOSED.json")
    ap.add_argument("--scores", type=Path, default=HERE / "results/STAGE4_SCORES.json")
    ap.add_argument("--private-dir", type=Path, required=True)
    ap.add_argument("--release-out", type=Path, default=HERE / "truth_after_freeze")
    ap.add_argument("--figures", type=Path, default=HERE / "figures")
    a = ap.parse_args()
    seal, scores = read(a.seal), read(a.scores)
    verify_closed(seal, scores)
    release = export_sanitized(a.private_dir, a.release_out, seal)
    images, scored = sealed_images(seal), map_scores(scores)
    records = release["records"]
    clean = sorted((r for r in records if r["noise"] == "clean"), key=lambda r: r["shape_slot"])
    digits = [r for r in clean if r["shape_slot"] < 9]
    noisy = sorted((r for r in records if r["noise"] == "noisy"),
                   key=lambda r: (r["shape_slot"], r["noise_realization"]))
    if (len(clean), len(digits), len(noisy)) != (14, 9, 8):
        raise ValueError("Post-seal groups must be 14 clean, nine digits, eight noisy")
    render_grid("fresh_all_shapes", clean, a.release_out, images, scored, a.figures)
    render_grid("fresh_digits_1to9", digits, a.release_out, images, scored, a.figures)
    render_grid("fresh_noise", noisy, a.release_out, images, scored, a.figures)
    cost_quality(scores, a.figures)
    paired = paired_analysis(scores)
    write(HERE / "results/STAGE4_PAIRED_ANALYSIS.json", paired)
    print(json.dumps({"status": "POST_SEAL_FIGURES_AND_SEEDLESS_RELEASE_COMPLETE",
                      "clean_objects": len(clean), "noisy_observations": len(noisy),
                      "released_npz": len(records),
                      "figures": [str(a.figures / (name + ".png")) for name in
                                  ("fresh_all_shapes", "fresh_digits_1to9", "fresh_noise", "fresh_cost_quality")]}))


if __name__ == "__main__":
    main()
