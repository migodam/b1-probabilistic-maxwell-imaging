"""Summarize strictly checked Manitoba floors and separate calibration proxies."""
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from physical import HERE, TARGETS

LABELS = dict(zip(TARGETS, ("Large nylon cylinder", "Small nylon cylinder",
                            "Two nylon cylinders", "Polyethylene E phantom"), strict=True))
PROXY_NAMES = (
    ("fixed_gain_empty", "drift_over_scatter"),
    ("recalibrated_before_after", "refit_before_after_empty_difference_relative"),
    ("two_metal_references", "two_reference_columns_difference_relative"),
    ("empty_plus_held_metal", "drift_plus_holdout_over_scatter"),
)
CAPS = (.25, .5, .75, 1.)


def load(name):
    return json.loads((HERE/name).read_text())


def comparison_rows(physical, verified):
    diagnostics = {row["id"]:row for row in physical["cases"]}
    compared = []
    for x in verified["rows"]:
        d = diagnostics[x["case_id"]]
        relative = x["verified_lower_relative"]
        row = {key:x[key] for key in ("case_id","target","frequency_GHz","cap","arm","field_norm",
                                     "verified_lower","verified_lower_relative","numeric_bound_relative")}
        row["positive_groups"] = sum(any(q.get("verified_lower",0)>0 for q in g["attempts"])
                                     for g in x["groups"])
        for name,key in PROXY_NAMES:
            row[name+"_relative"] = d[key]
            row["margin_after_"+name] = relative - d[key]
        row["margin_after_both_cross_checks"] = (
            relative - max(row["two_metal_references_relative"],row["empty_plus_held_metal_relative"]))
        row["nominal80_prediction_relative"] = d["nominal80_rel"]
        row["nominal_grid_difference_relative"] = d["nominal_grid_difference_over_scatter"]
        compared.append(row)
    return compared


def assemble(physical, numeric, verified, compared):
    by_arm = defaultdict(list)
    for cap in CAPS:
        for arm in ("real","positive_optical","free_optical"):
            population = [r for r in compared if r["cap"] == cap and r["arm"] == arm]
            per_target = {}
            for target in TARGETS:
                target_rows = [i for i in population if i["target"] == target]
                per_target[target] = {
                    "verified_relative_max": max(r["verified_lower_relative"] for r in target_rows),
                    "exceeds_both_empirical_proxies_frequencies_GHz": [r["frequency_GHz"]
                        for r in target_rows if r["margin_after_both_cross_checks"] > 0],
                }
            by_arm[arm].append(dict(
                cap=cap, observation_count=len(population), by_target=per_target,
                strict_positive=sum(r["verified_lower"]>0 for r in population),
                exceeds_both_cross_checks=sum(r["margin_after_both_cross_checks"]>0 for r in population),
                exceeds_each_check={name:sum(r["margin_after_"+name]>0 for r in population)
                                    for name,_ in PROXY_NAMES},
                numeric_but_unverified=sum(r["numeric_bound_relative"]>0 and r["verified_lower"]==0
                                           for r in population)))
    profiles = [p for freq in verified["costs"] for p in freq["proofs"]]
    return dict(
        bounds_are_outward_for_frozen_discrete_operator=True,
        physical_model_error_hard_upper=None,
        profiles_verified=len(profiles),
        largest_compact_proof_dimension=max(p["proof_dimension"] for p in profiles),
        physical_preparation_seconds=physical["wall_seconds"],
        numeric_scan_seconds=numeric["wall_seconds"],
        verification_seconds=verified["wall_seconds"],
        by_arm=by_arm,
        comparative_empirical_examples=[r for r in compared if r["target"]=="polyethylene_e" and
                                        r["cap"]==.5])


def plot(compared):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"pdf.fonttype":42})
    fig, axes = plt.subplots(2,2,sharex=True,sharey=True,figsize=(8.3,5.35))
    fig.subplots_adjust(left=.098,right=.985,bottom=.20,top=.9,hspace=.24,wspace=.12)
    palette=("#1f5068","#197e6b","#b65f15","#af4269")
    for ax,target in zip(axes.flat,TARGETS,strict=True):
        for cap,color,symbol in zip(CAPS,palette,("o","s","^","D"),strict=True):
            rows=[r for r in compared if r["target"]==target and r["cap"]==cap
                  and r["arm"]=="positive_optical"]
            xx=[r["frequency_GHz"] for r in rows]
            ax.plot(xx,[100*r["verified_lower_relative"] for r in rows],marker=symbol,
                    ms=4,lw=1.3,label=f"Cap {cap:g}",color=color)
        ax.plot(xx,[100*r["two_metal_references_relative"] for r in rows],ls="--",lw=1.25,
                label="Two metal references",color="#4d5157")
        ax.plot(xx,[100*r["empty_plus_held_metal_relative"] for r in rows],ls=":",lw=1.65,
                label="Empty drift + held metal",color="#7d6260")
        ax.set_title(LABELS[target],fontweight="bold",fontsize=9.5,loc="left")
        ax.set_ylim(-2,100)
        ax.set_yticks((0,20,40,60,80,100))
        ax.set_xticks((3,3.5,4,4.5,5,5.5,6))
        ax.grid(axis="y",color="#d7dde0",lw=.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Positive optical certificates on Manitoba measurements",
                 fontsize=13,fontweight="bold",y=.99)
    fig.supxlabel("Frequency (GHz)",y=.128)
    fig.supylabel("Percent of measured scattering norm",x=.014)
    lines, labels=axes[0,0].get_legend_handles_labels()
    fig.legend(lines,labels,loc="lower center",ncol=3,frameon=False,fontsize=8,
               bbox_to_anchor=(.52,.012),columnspacing=1.65,handlelength=2.1)
    fig.savefig(HERE/"calibration_vs_certificate.png",dpi=230)
    fig.savefig(HERE/"calibration_vs_certificate.pdf")


if __name__=="__main__":
    physical=load("physical_audit.json")
    numeric=load("numeric_scan.json")
    verified=load("verified_scan.json")
    compared=comparison_rows(physical,verified)
    with (HERE/"calibration_margins.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(compared[0]))
        writer.writeheader()
        writer.writerows(compared)
    (HERE/"summary.json").write_text(json.dumps(assemble(physical,numeric,verified,compared),
                                                indent=2,allow_nan=False)+"\n")
    plot(compared)
