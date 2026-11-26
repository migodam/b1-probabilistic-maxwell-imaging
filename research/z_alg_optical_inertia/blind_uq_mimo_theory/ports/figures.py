"""Figure and tables for the discrete MIMO sensitivity experiment."""
from collections import defaultdict
import csv
import json
from pathlib import Path
import statistics

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = Path(__file__).resolve().parent
FIG = HERE/"figures"
RES = HERE/"results"
STRENGTHS = (16, 48, 144)
COUNTS = (6,8,12,16)


def values(row, column):
    return float(row[column])


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    rows = read_csv(RES/"mechanism_witness.csv")
    lookup = {(int(r["frequency_MHz"]), int(r["q"]), float(r["radiation_scale_ohms"]),
               float(r["calibration_eta"])): r for r in rows}
    dual = {}
    for (f, q, strength, eta), row in lookup.items():
        if f != 400:
            continue
        extra = lookup[460,q,strength,eta]
        r = np.hypot(float(row["noise_radius"]), float(extra["noise_radius"]))
        # Readout sqrt(W): split 1 W available equally over two actual frequencies.
        # Fixed absolute instrument floor is incurred on both bands.
        def stacked(field):
            return np.hypot(values(row,field), values(extra,field))/np.sqrt(2)
        nominal_dist = stacked("nominal_readout_frobenius_separation")
        noisy_dist = stacked("readout_frobenius_separation")
        beta = np.hypot(
            float(row["noise_radius"]) + values(row,"calibration_radius_full")/np.sqrt(2),
            float(extra["noise_radius"]) + values(extra,"calibration_radius_full")/np.sqrt(2))
        dual[q,strength,eta] = dict(q=q, radiation_scale_ohms=strength, calibration_eta=eta,
                                    pair_distance=float(noisy_dist), joint_noise_radius=float(r),
                                    nominal_pair_distance=float(nominal_dist),
                                    optional_full_transfer_joint_radius=float(beta),
                                    actual_pair_separation_over_2r=float(noisy_dist/(2*r)),
                                    nominal_pair_separation_over_inflated_2r=float(nominal_dist/(2*beta)))
    with (RES/"dual_band_grid.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(next(iter(dual.values()))))
        w.writeheader()
        w.writerows(dual.values())

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.3), layout="constrained")
    color = TwoSlopeNorm(vmin=0, vcenter=1, vmax=5)
    panels = []
    for eta in (0.0, 0.005, 0.025):
        panels.append((f"400 MHz, δ/R = {eta:.3f}", [values(lookup[400,q,float(s),eta],
                        "nominal_pair_sep_over_two_inflated_r") for s in STRENGTHS for q in COUNTS]))
    for eta in (0.0, 0.005, 0.025):
        panels.append((f"400 + 460 MHz, δ/R = {eta:.3f}", [dual[q,s,eta][
                        "nominal_pair_separation_over_inflated_2r"] for s in STRENGTHS for q in COUNTS]))
    for i, (ax, (title, data)) in enumerate(zip(axes.ravel(), panels, strict=True)):
        matrix = np.asarray(data).reshape(len(STRENGTHS), len(COUNTS))
        rendered = ax.imshow(np.clip(matrix,0,5), norm=color, cmap="RdYlBu", origin="lower", aspect="auto")
        for r, c in np.ndindex(matrix.shape):
            val = matrix[r,c]
            ax.text(c, r, f"{val:.2f}", ha="center", va="center",
                    color="white" if val < .27 or val > 4.0 else "#112031", fontweight="semibold")
        ax.set_title(title, loc="left", fontsize=10, fontweight="semibold")
        ax.set_xticks(range(len(COUNTS)), [f"{q}×{q}" for q in COUNTS])
        ax.set_yticks(range(len(STRENGTHS)), [str(a) for a in STRENGTHS])
        if i % 3 == 0:
            ax.set_ylabel("Radiation scale (Ω)")
    fig.suptitle("Separation of the same two full-wave materials", fontsize=15, fontweight="semibold")
    cbar = fig.colorbar(rendered, ax=axes.ravel().tolist(), shrink=.85, extend="max")
    cbar.set_label("Pair distance / combined radii")
    fig.savefig(FIG/"mimo_phase.png", dpi=195)
    plt.close(fig)

    with np.load(RES/"mechanism_witness.npz") as w:
        with_roi, without_roi, roi = w["with_region"], w["without_region"], w["roi"]
    fig, axes = plt.subplots(1,3, figsize=(11,3.2), layout="constrained")
    im = axes[0].imshow(with_roi, origin="lower", vmin=0, vmax=4.5, cmap="viridis",
                       extent=(-1,1,-1,1))
    axes[0].set_title("Mean 0.75 in outlined ROI")
    axes[1].imshow(without_roi, origin="lower", vmin=0, vmax=4.5, cmap="viridis",
                   extent=(-1,1,-1,1))
    axes[1].set_title("Exactly zero in outlined ROI")
    for ax in axes[:2]:
        axis = (np.arange(len(roi))+.5)*2/len(roi)-1
        ax.contour(axis, axis, roi, levels=[0.5], colors="#ff9a4f", linewidths=1.3)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
    axes[2].plot(COUNTS,[values(lookup[400,q,48.0,0.0],"actual_pair_sep_over_two_r") for q in COUNTS],
                 marker="o", label="400 MHz", color="#18629B")
    axes[2].plot(COUNTS,[dual[q,48.0,0.0]["actual_pair_separation_over_2r"] for q in COUNTS],
                 marker="s", label="400 + 460 MHz", color="#7C3E87")
    axes[2].axhline(1, color="#686868", linestyle="--", linewidth=1)
    axes[2].set_ylabel("Pair distance / combined radii")
    axes[2].set_xticks(COUNTS)
    axes[2].set_xlabel("Same transmit and receive port count")
    axes[2].legend(frameon=False, loc="best")
    fig.colorbar(im, ax=axes[:2].tolist(), shrink=.77, label="Real contrast")
    fig.savefig(FIG/"mimo_mechanism_pair.png", dpi=195)
    plt.close(fig)
    summary = {}
    for n in (32,96):
        path = RES/f"array_grid_n{n}.csv"
        if not path.exists():
            continue
        rows = read_csv(path)
        groups = defaultdict(list)
        for r in rows:
            groups[int(r["frequency_MHz"]),int(r["q"]),int(float(r["radiation_scale_ohms"])),
                   round(values(r,"normalized_calibration_error"),5)].append(r)
        summary[n] = dict(rows=len(rows), materials=len({r["case_id"] for r in rows}),
                          max_transfer_identity_error=max(values(r,"transfer_identity_max_abs") for r in rows),
                          max_energy_balance_error_W=max(values(r,"energy_balance_max_W") for r in rows),
                          min_loaded_antenna_passivity_ohm=min(values(r,"loaded_resistance_min_ohm") for r in rows),
                          rows_by_setting={"_".join(map(str,k)):{
                              field: statistics.median(values(r,field) for r in v) for field in (
                                  "off_diagonal_Z_operator_ohms", "loaded_total_condition", "accepted_total_W",
                                  "receiver_10pct_modes", "illumination_10pct_modes",
                                  "Y_train_calibration_difference", "Y_train_calibration_bound_full",
                              )} for k,v in groups.items()})
    (RES/"mimo_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
