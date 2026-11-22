"""Summarize the masked Fresnel scans and their model/calibration error ledgers."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from physical import FREQUENCIES, HERE, NAMES, RESULTS
from scan import CAPS, RADII

FAMILIES = ("real", "optical", "two_zone")
LABELS = {"real": "Real sector", "optical": "Global nonzero optical",
          "two_zone": "Two-zone weights"}
COLORS = {"real": "#747982", "optical": "#05668d", "two_zone": "#e36414"}


def csv_write(path, data):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)


def effective_radius(delta, eta):
    return (delta + eta)/(1.-eta) if eta < 1. else np.inf


def empty_scan_repeats(output):
    rows = []
    for frequency in FREQUENCIES:
        records = {}
        for name in NAMES:
            path = output/f"{name}_{frequency}GHz_ports.npz"
            if path.exists():
                with np.load(path) as item:
                    records[name] = tuple(item[k] for k in ("source_degrees", "raw_empty", "raw_scattered"))
        for i, (name, (angles, empty, y)) in enumerate(records.items()):
            for other, (other_angles, other_empty, other_y) in tuple(records.items())[i+1:]:
                degrees, ia, ib = np.intersect1d(angles, other_angles, return_indices=True)
                first = empty[:, ia]
                second = other_empty[:, ib]
                available = np.isfinite(first) & np.isfinite(second)
                absolute = np.linalg.norm((first-second)[available])
                rows.append(dict(first_case=name, second_case=other, GHz=frequency,
                    shared_source_degrees=",".join(map(str,degrees)), same_physical_ports=int(available.sum()),
                    difference_over_first_empty=float(absolute/np.linalg.norm(first[available])),
                    difference_over_first_scatter=float(absolute/np.linalg.norm(y[:,ia][available])),
                    difference_over_second_scatter=float(absolute/np.linalg.norm(other_y[:,ib][available]))))
    if rows:
        csv_write(output/"REPEATED_EMPTY_PORTS.csv", rows)


def sensitivity_table(scan, verified):
    rows, overview = [], []
    for record in scan["cases"]:
        delta = record["empty_calibration"]["rms_transferred_over_all_scatter"]
        name = record["name"]
        frequency = record["frequency_GHz"]
        model = record["model_audit"]
        eta = abs(complex(*model["diagnostic_scatter_gain_real_imag"]) - 1.)
        with_gain = effective_radius(delta, eta)
        v = verified.get(name, {})
        rows_this_case = []
        for cap in CAPS:
            for family in FAMILIES:
                candidates = [entry for entry in record["grid"]
                              if entry["u"] == cap and entry["family"] == family]
                best = max(candidates, key=lambda row: row["best"]["numerical_bound"])
                winner = best["best"]
                bound = winner["numerical_bound"]
                indef = max((entry["best_indefinite"] for entry in candidates
                             if entry["best_indefinite"] is not None),
                             key=lambda entry: entry["numerical_bound"], default=None)
                strict = max((row["certificate_distance_lower"] for row in v.get("results", [])
                              if row["u"] == cap and row["family"] == family
                              and "certificate_distance_lower" in row), default=None)
                rows_this_case.append(dict(
                    case=name, GHz=frequency, cap_chi=cap, family=family, profile=best["profile"],
                    lower_numeric=bound, lower_outward=strict,
                    partition=winner["partition"], mu=winner["mu"],
                    retained=winner["retained_channels"], measured=record["measured_channels"],
                    informative_channels=winner["informative_channels"],
                    inertia=",".join(map(str, winner["inertia"])),
                    numerical_joint_indefinite_best=None if indef is None else indef["numerical_bound"],
                    numerical_joint_indefinite_groups=None if indef is None else indef["informative_indefinite_groups"],
                    numerical_joint_indefinite_partition=None if indef is None else indef["partition"],
                    numerical_joint_indefinite_inertia=None if indef is None else ",".join(map(str,indef["inertia"])),
                    bcrit_empty_rms_no_gain=bound-delta,
                    bcrit_empty_rms_and_diagnostic_gain=bound-with_gain if np.isfinite(with_gain) else None,
                    bcrit_outward_empty_rms_no_gain=None if strict is None else strict-delta,
                    bcrit_outward_empty_rms_and_diagnostic_gain=None if strict is None
                        or not np.isfinite(with_gain) else strict-with_gain,
                    **{f"margin_if_radius_{int(radius*100)}pct":bound-radius for radius in RADII}
                ))
        def largest(family, threshold, field):
            return max((row["cap_chi"] for row in rows_this_case if row["family"] == family
                        and row[field] is not None and row[field] > threshold), default=None)

        nominal = model["relative_residuals"]
        overview.append(dict(
            case=name, GHz=frequency, sources=record["source_count"], measured=record["measured_channels"],
            calibration_hold_relative_incident=record["empty_calibration"]["held_empty_relative"],
            calibration_hold_relative_scatter=record["empty_calibration"]["held_empty_error_over_held_scatter"],
            calibration_empty_rms_transfer_all_scatter=delta,
            diagnostic_gain_abs=model["diagnostic_scatter_gain_abs"],
            diagnostic_gain_phase_deg=model["diagnostic_scatter_gain_phase_deg"],
            diagnostic_gain_deviation_from_unit=eta,
            nominal_model_relative_residual_all=nominal["all"]["nominal60"],
            nominal_model_relative_residual_held_tx=nominal["held_sources"]["nominal60"],
            nominal_model_relative_residual_held_rx=nominal["held_receivers_on_fit_sources"]["nominal60"],
            gained_nominal_model_relative_residual_held_tx=nominal["held_sources"]["nominal60_with_train_gain"],
            gained_nominal_model_relative_residual_held_rx=nominal["held_receivers_on_fit_sources"]["nominal60_with_train_gain"],
            nominal_fine_coarse_shift_relative_all=nominal["all"]["discretization_shift"],
            air_only_outward=v.get("air_only_discrete_distance_outward_lower"),
            **{f"{f}_largest_numeric_over_empty_proxy": largest(f, delta, "lower_numeric") for f in FAMILIES},
            **{f"{f}_largest_outward_over_empty_proxy": largest(f, delta, "lower_outward") for f in FAMILIES},
            **{f"{f}_largest_numeric_over_empty_plus_gain": largest(f, with_gain, "lower_numeric")
               for f in FAMILIES}
        ))
        rows.extend(rows_this_case)
    return rows, overview


def plot(rows, overview, indefinite):
    fig, axes = plt.subplots(3, 3, figsize=(12.6, 9.2), sharex=True, sharey=True)
    for record in overview:
        name, ghz = record["case"], record["GHz"]
        ax = axes[NAMES.index(name), FREQUENCIES.index(ghz)]
        for family in FAMILIES:
            selected = [r for r in rows if r["case"] == name and r["GHz"] == ghz and r["family"] == family]
            ax.plot([r["cap_chi"] for r in selected], [100*r["lower_numeric"] for r in selected],
                    marker="o", color=COLORS[family], label=LABELS[family], lw=1.45, ms=3.)
            strict = [r for r in selected if r["lower_outward"] is not None]
            ax.scatter([r["cap_chi"] for r in strict], [100*r["lower_outward"] for r in strict],
                       facecolors="none", edgecolors=COLORS[family], s=63, zorder=5)
        delta = record["calibration_empty_rms_transfer_all_scatter"]
        eta = record["diagnostic_gain_deviation_from_unit"]
        ax.axhline(100*delta, color="#317b41", ls="--", lw=1.2, label="Held-empty RMS proxy")
        if eta < 1.:
            ax.axhline(100*effective_radius(delta, eta), color="#913d51", ls=":", lw=1.35,
                       label="Held-empty plus diagnostic gain")
        proof = indefinite.get(ghz)
        if proof and proof.get("lower_outward", 0.) > 0. and proof["proposal"]["case"] == name:
            ax.scatter([proof["proposal"]["u"]], [100*proof["lower_outward"]],
                       marker="^", facecolors="white", edgecolors="#1f1f1f", s=35, zorder=6)
        ax.set_title(f"{name.removesuffix('TM')} • {ghz} GHz", fontsize=10)
        ax.set_xscale("log")
        ax.set_xlim(.0038, 3.)
        ax.set_ylim(0., 100.)
        ax.grid(alpha=.20)
    axes[0, 0].legend(loc="upper right", fontsize=6.8)
    fig.supxlabel("Tested cap on real contrast χ = εᵣ − 1 (log scale)", y=.038)
    fig.supylabel("Lower bound on measured-channel residual (%)", x=.021)
    fig.suptitle("Fresnel 2005 TM • calibrated data masked by physical ports", y=.985, fontsize=15)
    fig.text(.5, .012, "Rings and triangles are outward discrete checks (positive and indefinite inertia). Horizontal lines are sensitivity scenarios.",
             ha="center", fontsize=8.5)
    fig.tight_layout(rect=(.035, .045, 1., .966), h_pad=1.15)
    fig.savefig(HERE.parent / "figures/fresnel_margins.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def write_tables(summary, rows, indefinite):
    def pct(x):
        return "—" if x is None else f"{100*x:.2f}%"

    def capstr(x):
        return "—" if x is None else f"{x:g}"

    text = [
        "# Fresnel 2005 误差和余量",
        "",
        "| 观测 | 空场留出相对空场 | 空场留出相对同口散射 | 空场 RMS 转至全观测 | 名义模型留出发射 / 接收 | 60→120 网格差 | 数值通过的最大 χ | 外向复核通过的最大 χ |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for record in summary:
        label = f"{record['case'].removesuffix('TM')} {record['GHz']} GHz"
        text.append("| " + " | ".join((label,
            pct(record["calibration_hold_relative_incident"]),
            pct(record["calibration_hold_relative_scatter"]),
            pct(record["calibration_empty_rms_transfer_all_scatter"]),
            pct(record["nominal_model_relative_residual_held_tx"]) + " / "
                + pct(record["nominal_model_relative_residual_held_rx"]),
            pct(record["nominal_fine_coarse_shift_relative_all"]),
            " / ".join(capstr(record[f"{f}_largest_numeric_over_empty_proxy"]) for f in FAMILIES),
            " / ".join(capstr(record[f"{f}_largest_outward_over_empty_proxy"]) for f in FAMILIES))) + " |")
    text.extend(["", "后两列依次为全局实扇区、全局非零光学、双区乘子。通过指下界超过空场 RMS 代理，模型与标定上界仍需单独给定。",
                 "", "## 九例统一检验 χ≤0.5", "",
                 "| 观测 | 数值下界（实 / 光 / 双区） | 空场后余量（实 / 光 / 双区） | 空场加诊断增益后余量（实 / 光 / 双区） |",
                 "|---|---:|---:|---:|"])
    for record in summary:
        same = [row for row in rows if row["case"] == record["case"] and row["GHz"] == record["GHz"]
                and row["cap_chi"] == .5]
        label = f"{record['case'].removesuffix('TM')} {record['GHz']} GHz"
        text.append("| " + label + " | " + " | ".join(
            " / ".join(pct(r[key]) for r in same)
            for key in ("lower_numeric", "bcrit_empty_rms_no_gain", "bcrit_empty_rms_and_diagnostic_gain")) + " |")
    text.extend(["", "以每例全部真实散射读数的 Frobenius 范数为分母。负值表示即使取模型误差 b=0，该敏感性场景也超出了算得的下界。",
                 "", "## 分组不定证书", ""])
    for freq, result in indefinite.items():
        if result["status"] == "discrete_indefinite_outward":
            p = result["proposal"]
            counts = [str(row.get("K_negative_verified", "未闭合")) for row in result["groups"]]
            text.append(f"- {p['case']} {freq} GHz，χ≤{p['u']:g}，逐组已证负方向 "
                        f"{' / '.join(counts)}，外向残差下界 {pct(result['lower_outward'])}，"
                        f"空场代理后余量 {pct(result['empirical_empty_RMS_margin'])}。")
        elif "proposal" in result:
            text.append(f"- {freq} GHz 的不定候选未闭合。{result.get('reason', '')}")
        else:
            text.append(f"- {freq} GHz 在外向复核采用的 μ≥10⁻⁷ 候选中没有正的多源不定下界。"
                        "较小 μ 的数值结果见敏感性表。")
    text.extend(["", "[全部临界余量](SENSITIVITY.csv) · [逐例误差账本](ERROR_LEDGER.csv) · [九例图](../../figures/fresnel_margins.png)", ""])
    (RESULTS / "TABLES_ZH.md").write_text("\n".join(text))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS)
    args = parser.parse_args()
    rows, overview, indefinite = [], [], {}
    for frequency in FREQUENCIES:
        scan = json.loads((args.output / f"scan_{frequency}GHz.json").read_text())
        path = args.output / f"verified_{frequency}GHz.json"
        verified = json.loads(path.read_text())["cases"] if path.exists() else {}
        per_scan, per_overview = sensitivity_table(scan, verified)
        rows.extend(per_scan)
        overview.extend(per_overview)
        indef_path = args.output / f"indefinite_{frequency}GHz.json"
        if indef_path.exists():
            indefinite[frequency] = json.loads(indef_path.read_text())
    csv_write(args.output / "SENSITIVITY.csv", rows)
    csv_write(args.output / "ERROR_LEDGER.csv", overview)
    empty_scan_repeats(args.output)
    write_tables(overview, rows, indefinite)
    plot(rows, overview, indefinite)


if __name__ == "__main__":
    main()
