"""Summarize independently timed CIE, CIE+corrector and original direct replays."""

import csv
import json
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
METHODS = ("all_cluster","point_to_fit","neumann_directional","shifted_threshold","direct_point_to_fit")
NAMES = ("All boxes", "Point to first fit", "Directional Neumann", "Shifted threshold", "Direct point + 42")
COLORS = ("#b2825d","#666970","#3c6b98","#2e877f","#8b58a2")


def table(path):
    return json.loads(path.read_text())


def merge():
    original = table(HERE/"results"/"online.json")["cases"]
    rows, checks = [], {}
    for schedule, directory in [("cie_only","results"),("cie_then_corrector","corrected_results")]:
        path = HERE/directory
        online, scored = table(path/"online.json")["cases"], table(path/"offline.json")
        checked = [row for scored_item in scored.values() for row in scored_item["arb_checked_shifted_exclusions"]]
        checks[schedule] = dict(arbitrary_precision_arb_exclusions=len(checked),
            contained_truth_total=sum(len(s["contains_truth"]) for s in scored.values()),
            arb_independent_wall=sum(x["independent_check"]["wall_seconds"] for x in checked),
            arb_independent_cpu=sum(x["independent_check"]["cpu_seconds"] for x in checked))
        with np.load(path/"predictions.npz",allow_pickle=False) as materials:
            for case, item in online.items():
                assert np.array_equal(materials[case+"__all_cluster__material"],materials[case+"__shifted_threshold__material"])
                raw = item["methods"]
                for method,a in raw.items():
                    rows.append(serialize(schedule,method,case,item["cohort"],a,scored[case]["methods"][method]))
    direct = table(HERE/"direct_results"/"online.json")["cases"]
    known = table(HERE/"direct_results"/"offline.json")
    for case,a in direct.items():
        rows.append(serialize("direct_corrector","direct_point_to_fit",case,original[case]["cohort"],a,known[case]))

    output = HERE/"per_object.csv"
    with output.open("w",newline="") as fd:
        writer=csv.DictWriter(fd,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summarized = {}
    for cohort in ("all","development","confirmation"):
        values = [row for row in rows if cohort in ("all",row["cohort"])]
        group = {}
        for phase in ("cie_only","cie_then_corrector","direct_corrector"):
            subset = [row for row in values if row["schedule"]==phase]
            by_method = {}
            for method in METHODS:
                picked = [row for row in subset if row["method"]==method]
                if not picked:
                    continue
                by_method[method] = aggregate(picked)
            if "shifted_threshold" in by_method:
                fast = {row["case"]:row["wall_seconds"] for row in subset if row["method"]=="shifted_threshold"}
                by_method["shifted_threshold"]["median_all_over_shifted"] = median(row["wall_seconds"]/fast[row["case"]]
                    for row in subset if row["method"]=="all_cluster")
                by_method["shifted_threshold"]["median_point_over_shifted"] = median(row["wall_seconds"]/fast[row["case"]]
                    for row in subset if row["method"]=="point_to_fit")
            group[phase] = by_method
        summarized[cohort] = group
    receipt = dict(scores=summarized, independent_checks=checks,
                   exact_equal_all_shifted_predicted_materials_each_phase=12)
    (HERE/"summary.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    return receipt


def serialize(phase,method,case,cohort,info,eval_):
    is_direct = method=="direct_point_to_fit"
    radius=info["radius"]
    return dict(schedule=phase,method=method,case=case,cohort=cohort,winner=info["winner"],
                visits=len(info["visited"]),wall_seconds=info["wall_seconds"],cpu_seconds=info["cpu_seconds"],
                physical_lu=(info["dense_optimization_counts"] if is_direct else info["physical_counts"])["lu"],
                pr_steps=0 if is_direct else info["cie_counts"]["PR_CG_steps"],
                residual_upper_over_radius=(info["verified_residual_upper"] if is_direct else info["verified_winner_upper"])/radius,
                feasible=(eval_["feasible"] if is_direct else eval_["strict_data_feasible"]),
                full_relative_l2=eval_["full_relative_l2"],inside_relative_l2=eval_["inside_relative_l2"],
                background_rms=eval_["background_rms"],heldout_receivers=eval_["heldout_receivers"],
                heldout_sources=eval_["heldout_sources"])


def aggregate(picked):
    metrics=("full_relative_l2","inside_relative_l2","background_rms","heldout_receivers","heldout_sources")
    totals=("wall_seconds","cpu_seconds","physical_lu","pr_steps","visits","feasible")
    return dict(materials=len(picked),
                totals={name:sum(row[name] for row in picked) for name in totals},
                residual_upper_over_radius_range=[min(row["residual_upper_over_radius"] for row in picked),
                                                  max(row["residual_upper_over_radius"] for row in picked)],
                median={name:median(row[name] for row in picked) for name in metrics})


def figure(data):
    tally=data["scores"]["all"]
    fig,axes=plt.subplots(1,2,figsize=(10.7,4.),layout="constrained",sharex=True)
    for ax,stage,title in zip(axes,("cie_only","cie_then_corrector"),("CIE alone","CIE followed by full-wave correction"),strict=True):
        scored=tally[stage]
        if stage=="cie_then_corrector":
            scored={**scored,**tally["direct_corrector"]}
        for i,(ident,color) in enumerate(zip(METHODS,COLORS,strict=True)):
            if ident not in scored:
                continue
            record=scored[ident]["totals"]
            ax.barh(i,record["wall_seconds"],color=color,height=.61)
            ax.text(record["wall_seconds"]+.11,i,f"{record['feasible']}/12 feasible",va="center",fontsize=8.5)
        ax.set_yticks(range(len(NAMES)),NAMES)
        ax.set_ylim(4.65,-.65)
        ax.set_title(title,loc="left",fontsize=11.5)
        ax.set_xlim(0,8.4)
        ax.set_xlabel("Paid online wall time over all 12 materials (s)")
        ax.grid(axis="x",alpha=.18)
        ax.set_axisbelow(True)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="y",length=0)
    axes[0].set_yticks(range(4),NAMES[:4])
    fig.suptitle("Candidate screening reduces CIE work, while point ranking finds feasible materials first",
                 fontsize=12.3)
    fig.savefig(HERE/"paired_cost_fit.png",dpi=185)
    fig.savefig(HERE/"paired_cost_fit.pdf")


if __name__=="__main__":
    figure(merge())
