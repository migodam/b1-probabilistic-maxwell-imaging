"""Compact case-level public receipts from the blind experiment."""
import csv
import json

from shared import OUT, write_json


def read(case_path, name):
    p = case_path/name
    return json.loads(p.read_text()) if p.exists() else {}


def main():
    output = {}
    for root in sorted(OUT.glob("n*_c??")):
        base = read(root, "baseline.json")
        if not base:
            continue
        outputs = {"dual":root}
        if (root/"stationary"/"baseline.json").exists():
            outputs["stationary"] = root/"stationary"
        for family, path in outputs.items():
            b = read(path, "baseline.json")
            c = read(path, "continuation.json")
            b1 = read(path, "b1.json")
            b2 = read(path, "b2.json")
            h = read(path,"held.json").get("predictions", {})
            verified = read(path,"one_sided_witness.json")
            selected_roi = b1.get("selected_roi")
            if selected_roi and not verified:
                verified = {"selected":read(path,f"native_one_sided_{selected_roi}.json")}
            row = dict(case=root.name, family=family, side=b["side"],
                       baseline_qualified=b.get("converged", b.get("stable")),
                       baseline_updates=sum(s["updates"] for s in b.get("stages", [])) or b.get("updates"),
                       baseline_wall=b["wall_s"], baseline_cpu=b["cpu_s"],
                       baseline_r=b["original_physics"]["residual_over_radius"],
                       baseline_held=h.get("baseline",{}).get("held_relative_residual"),
                       rois=";".join(str(p["cells"]) for p in b["segmentation"]["components"]),
                       b1_zero_roi=selected_roi,
                       b1_zero_r=b1.get("selected", {}).get("residual_over_radius"),
                       b1_held=h.get("b1",{}).get("held_relative_residual"),
                       b1_wall_control_r=b1.get("matched_wall_full_wave",{}).get("residual_over_radius"),
                       b1_wall_control_held=h.get("b1_wall_control",{}).get("held_relative_residual"),
                       b1_cie_match_r=b1.get("matched_cie_selected",{}).get("residual_over_radius"),
                       b1_one_sided_r=verified.get("selected",{}).get("residual_upper_over_radius"),
                       b1_target_fit_wall=b1.get("compared_paid_budget",{}).get("wall_s"),
                       b1_target_fit_cpu=b1.get("compared_paid_budget",{}).get("cpu_s"),
                       tighter_r=c.get("tightened",{}).get("residual_over_radius"),
                       tighter_delta_x=c.get("relative_material_change"),
                       tighter_updates=c.get("updates"),
                       tightened_qualified=c.get("stop")=="COEFFICIENT_RMS" if c else None,
                       k_neg_modes=b2.get("modes",{}).get("negative_kept"),
                       b2_mode_wall=b2.get("modes",{}).get("wall_s"),
                       b2_matched_wall=b2.get("common_extra_wall_s"))
            for arm in ("plain","ridge","random","dark"):
                raw = b2.get("arms",{}).get(arm,{})
                row.update({
                    f"b2_{arm}_r": raw.get("selected",{}).get("residual_over_radius"),
                    f"b2_{arm}_held":h.get("b2_"+arm,{}).get("held_relative_residual"),
                    f"b2_{arm}_energy_fraction":raw.get("actual_physical_k_negative_energy_fraction"),
                    f"b2_{arm}_time_r":raw.get("matched_wall_selected",{}).get("residual_over_radius"),
                    f"b2_{arm}_time_held":h.get("b2_"+arm+"_time",{}).get("held_relative_residual"),
                    f"b2_{arm}_wall":raw.get("end_to_end",{}).get("wall_s"),
                    f"b2_{arm}_cpu":raw.get("end_to_end",{}).get("cpu_s"),
                })
            output[root.name+"/"+family] = row
    write_json(OUT/"summary.json", output)
    with (OUT/"summary.csv").open("w") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(next(iter(output.values()))))
        writer.writeheader()
        writer.writerows(output.values())


if __name__ == "__main__":
    main()
