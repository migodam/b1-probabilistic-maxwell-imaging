"""Offline N2 graphics; run only after the main reader's blind receipt and scoring.

The sealed online selections are displayed unchanged. GT is used only for
offline panels and numerical descriptions; it never selects a displayed action.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESEARCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESEARCH / "agent_c"))
import offline_score as o

METHODS = ["H0_GN", "H0_DATA_RML", "R9_STYLE_FIXED8", "H1_FULL16", "H2_FULL16", "HA_FULL16"]
LABELS = ["H0 ordinary GN", "H0 data RML", "Fixed8", "H1 Full16", "H2 Full16", "HA Full16"]


def full_occ(value, roi):
    if value is None:
        return None
    a = np.asarray(value, float)
    if a.shape == (256,):
        b = np.zeros(576)
        b[roi] = a
        return b
    return a.reshape(576)


def topo_fields(material, occ, gt, cfg):
    out = {"material_" + k: v for k, v in o.structural(material.real, gt["material"].real, cfg["material_topology_threshold"]).items()}
    if occ is not None and gt["occupancy"] is not None:
        out.update({"occupancy_" + k: v for k, v in o.structural(occ, gt["occupancy"], cfg["occupancy_topology_threshold"]).items()})
    return out


def values(material, occ, view, cfg):
    if view.startswith("OCCUPANCY"):
        a = occ
    else:
        a = material.real
    if a is None:
        return None
    if view.endswith("MASK"):
        threshold = cfg["occupancy_topology_threshold"] if view.startswith("OCCUPANCY") else cfg["material_topology_threshold"]
        a = a >= threshold
    return np.asarray(a).reshape(24, 24)


def image(ax, data, view, cfg):
    ax.set_xticks([])
    ax.set_yticks([])
    if data is None:
        ax.text(.5, .5, "UNAVAILABLE", ha="center", va="center", transform=ax.transAxes, fontsize=7)
        return None
    if view.endswith("MASK") or view.startswith("OCCUPANCY"):
        return ax.imshow(data, cmap="gray_r", vmin=0, vmax=1, interpolation="none")
    return ax.imshow(data, cmap="viridis", vmin=cfg["host_real"], vmax=1.9, interpolation="none")


def topology_label(material, occ, gt, view, cfg):
    use_occ = view.startswith("OCCUPANCY")
    a, target = (occ, gt["occupancy"]) if use_occ else (material.real, gt["material"].real)
    if a is None or target is None:
        return "occupancy unavailable"
    threshold = cfg["occupancy_topology_threshold"] if use_occ else cfg["material_topology_threshold"]
    t = o.structural(a, target, threshold)
    return f"IoU {t['IoU']:.2f}  C{t['predicted_components']} H{t['predicted_holes']}"


def reconstruction(groups, truths, rows, out, cfg):
    parents = sorted({g["parent"] for g in groups}, key=lambda p: truths[("N", p)]["digit"])
    by = {(g["parent"], g["method"]): g for g in groups}
    files = []
    for view in ["MATERIAL_REAL", "OCCUPANCY", "MATERIAL_MASK", "OCCUPANCY_MASK"]:
        fig, axes = plt.subplots(len(parents), 8, figsize=(16.4, 2.12 * len(parents)), squeeze=False)
        last = None
        for ri, p in enumerate(parents):
            gt = truths[("N", p)]
            g0 = by[(p, METHODS[0])]
            a = g0["anchor"]
            anchor_occ = full_occ(a.get("occupancy"), g0["roi"])
            panels = [(gt["material"], gt["occupancy"], "GT", True, None),
                      (a["material"], anchor_occ, "Common anchor", a["legal"] and a["phi"] <= g0["noise_max"], a["phi"])]
            for method, label in zip(METHODS, LABELS):
                g = by[(p, method)]
                mat, occ, sid, qualified = o.selected_material(g)
                phi = g["phi"][g["selected"]] if g["selected"] is not None else g["anchor"]["phi"]
                panels.append((mat, occ, label, qualified, phi))
            for ci, (mat, occ, label, qualified, phi) in enumerate(panels):
                ax = axes[ri, ci]
                last = image(ax, values(mat, occ, view, cfg), view, cfg)
                if ri == 0:
                    ax.set_title(label, fontsize=10)
                if ci == 0:
                    ax.set_ylabel(f"Digit {gt['digit']}", fontsize=11)
                if ci == 0:
                    t = o.topology((gt["occupancy"] >= .5) if view.startswith("OCCUPANCY") else (gt["material"].real >= .75))
                    label2 = f"GT C{t['components']} H{t['holes']}"
                else:
                    label2 = topology_label(mat, occ, gt, view, cfg)
                    label2 += f"\nQ{int(qualified)} Phi {phi:.1f}"
                ax.set_xlabel(label2, fontsize=7)
        description = {"MATERIAL_REAL": "Real material (common scale 0.1 to 1.9)", "OCCUPANCY": "Occupancy (common scale 0 to 1)",
                       "MATERIAL_MASK": "Material mask Re chi >= 0.75", "OCCUPANCY_MASK": "Occupancy mask >= 0.5"}[view]
        fig.suptitle(f"N2: digits 1–9; sealed online actions | {description}\nRaw 24 x 24; no interpolation; topology FG4/BG8; Q = engine qualification", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, .97), h_pad=.9, w_pad=.3)
        if view == "MATERIAL_REAL":
            name = "DIGITS_123456789_RECONSTRUCTION.png"
        else:
            name = f"DIGITS_123456789_{view}.png"
        path = out / name
        fig.savefig(path, dpi=150)
        plt.close(fig)
        files.append(str(path))
    return files


def candidate_galleries(groups, truths, out, cfg):
    files = []
    parent = next(p for (ds, p), t in truths.items() if ds == "N" and t["digit"] == 8)
    gt = truths[("N", parent)]
    for method, label in zip(METHODS, LABELS):
        g = next(g for g in groups if g["parent"] == parent and g["method"] == method)
        a = g["anchor"]
        noise = np.asarray(g["z"].get("noise_qualified", g["finite"] & (g["phi"] <= g["noise_max"])), bool)
        legal = np.asarray(g["z"].get("legal", np.ones(len(g["materials"]), bool)), bool)
        panels = [(gt["material"], gt["occupancy"], "GT", ""),
                  (a["material"], full_occ(a.get("occupancy"), g["roi"]), "Original anchor", f"Phi {a['phi']:.1f}")]
        for i in range(len(g["materials"])):
            occ = None if g["occupancy"] is None else g["occupancy"][i]
            q = f"Q{int(g['qualified'][i])} N{int(noise[i])} L{int(legal[i])} Phi {g['phi'][i]:.1f}"
            panels.append((g["materials"][i], occ, f"#{i} {str(g['ids'][i])[:24]}", q))
        for view in ["MATERIAL_REAL", "OCCUPANCY", "MATERIAL_MASK", "OCCUPANCY_MASK"]:
            nc = 6
            nr = math.ceil(len(panels) / nc)
            height = 3.6 * nr + .9
            fig, axes = plt.subplots(nr, nc, figsize=(15, height), squeeze=False)
            for j, ax in enumerate(axes.ravel()):
                if j >= len(panels):
                    ax.axis("off")
                    continue
                mat, occ, title, flags = panels[j]
                image(ax, values(mat, occ, view, cfg), view, cfg)
                ax.set_title(title, fontsize=8)
                desc = topology_label(mat, occ, gt, view, cfg)
                ax.set_xlabel(desc + ("\n" + flags if flags else ""), fontsize=7)
                if j >= 2 and not g["qualified"][j-2]:
                    for spine in ax.spines.values():
                        spine.set_edgecolor("#c94435")
                        spine.set_linewidth(2)
            threshold = .5 if view.startswith("OCCUPANCY") else .75
            fig.suptitle(f"Digit 8 | {label} | {view.replace('_', ' ').lower()}\nAll actual candidates in immutable order; red frame Q0.\nQ = legal + noise + greedy distinct (.002); N = actual noise ball; L = legal | threshold {threshold}, FG4/BG8", fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, 1 - .9 / height), h_pad=3.0, w_pad=.8)
            path = out / f"DIGIT8_CLOUD_{method}_{view}.png"
            fig.savefig(path, dpi=140)
            plt.close(fig)
            files.append(str(path))
    return files


def tables(groups, truths, score_rows, out, cfg):
    selected_rows, candidate_rows, cloud_rows = [], [], []
    score_by = {(r["parent"], r["method"]): r for r in score_rows}
    for g in groups:
        gt = truths[(g["dataset"], g["parent"])]
        identity = {"dataset": g["dataset"], "parent": g["parent"], "digit": gt["digit"], "method": g["method"]}
        mat, occ, sid, qualified = o.selected_material(g)
        sr = score_by[(g["parent"], g["method"])]
        i = g["selected"]
        noise_pass = bool(g["finite"][i] and g["phi"][i] <= g["noise_max"]) if i is not None else bool(g["anchor"]["phi"] <= g["noise_max"])
        selected_rows.append(dict(**identity, selected_id=sid, qualified=qualified,
                                  selected_noise_pass=noise_pass,
                                  qualified_semantics="legal AND noise AND greedy material-RMS .002 distinct flag; Q0 may be duplicate",
                                  material_MSE=sr["online_selected_MSE"], ROI_material_MSE=sr["online_selected_ROI_MSE"],
                                  **topo_fields(mat, occ, gt, cfg)))
        noise = np.asarray(g["z"].get("noise_qualified", g["finite"] & (g["phi"] <= g["noise_max"])), bool)
        legal = np.asarray(g["z"].get("legal", np.ones(len(g["materials"]), bool)), bool)
        qualified_topology = []
        for i, mat in enumerate(g["materials"]):
            occ = None if g["occupancy"] is None else g["occupancy"][i]
            row = dict(**identity, index=i, candidate_id=str(g["ids"][i]), layer=int(g["layers"][i]),
                       qualified=bool(g["qualified"][i]), noise_qualified=bool(noise[i]), legal=bool(legal[i]),
                       Phi=float(g["phi"][i]), finite=bool(g["finite"][i]))
            if g["finite"][i]:
                row.update(topo_fields(mat, occ, gt, cfg))
                row["material_MSE"] = float(np.mean(abs(mat - gt["material"])**2))
            candidate_rows.append(row)
            if g["qualified"][i]:
                qualified_topology.append(row)
        cr = dict(**identity, raw_candidates=len(g["materials"]), qualified_raw_count=int(g["qualified"].sum()),
                  noise_failed_raw_count=int((~noise).sum()), primary_coverage=sr["coverage"][str(cfg["epsilon_primary"])],
                  qualified_unique_total=sr["qualified_unique_total"], nearest_material_MSE=sr["nearest_material_MSE"],
                  online_selected_MSE=sr["online_selected_MSE"], anchor_MSE=sr.get("anchor_MSE"))
        for name in ["material", "occupancy"]:
            rr = [r for r in qualified_topology if name + "_topology_exact" in r]
            cr[name + "_topology_exact_qualified_raw_count"] = sum(r[name + "_topology_exact"] for r in rr)
            cr[name + "_IoU_max_qualified_raw"] = max((r[name + "_IoU"] for r in rr), default=None)
        cloud_rows.append(cr)
    o.csv_write(out / "N2_SELECTED_STRUCTURE.csv", selected_rows)
    o.csv_write(out / "N2_CANDIDATE_STRUCTURE.csv", candidate_rows)
    o.csv_write(out / "N2_CLOUD_STRUCTURE.csv", cloud_rows)
    o.write(out / "N2_SELECTED_STRUCTURE.json", selected_rows)
    o.write(out / "N2_CLOUD_STRUCTURE.json", cloud_rows)
    return selected_rows, cloud_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--truth-manifest", required=True)
    ap.add_argument("--score-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    manifest, groups, cfg = o.load_predictions(args.manifest)
    score_rows = o.read(Path(args.score_dir) / "GROUP_SCORES.json")
    truths = {}
    for identity, path in o.truth_paths(args.truth_manifest).items():
        with np.load(path, allow_pickle=False) as z:
            truths[identity] = dict(material=np.asarray(z["material"], complex),
                                    occupancy=np.asarray(z["occupancy"], float) if "occupancy" in z.files else None,
                                    digit=int(np.asarray(z["target_digit"]).item()))
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    selected, clouds = tables(groups, truths, score_rows, out, cfg)
    files = reconstruction(groups, truths, score_rows, out, cfg)
    files += candidate_galleries(groups, truths, out, cfg)
    o.write(out / "FIGURE_CONTENTS.json", dict(files=files, new_physics_calls=0, online_actions="sealed selections unchanged",
            image_grid="raw24; no interpolation", topology="FG4/BG8; material .75, occupancy .5; all holes",
            candidate_display="all raw digit8 candidates in immutable order; all qualified and all failed flags retained",
            offline_only="GT panels and numerical descriptions; no GT-best replacement of online selection"))
    print(json.dumps({"figures": len(files), "selected_objects": len(selected), "cloud_objects": len(clouds), "out": str(out)}))


if __name__ == "__main__":
    main()
