"""Create the phase, material, and cost figures from frozen decision records."""
from collections import Counter
import json

from common import BGS, PRIVATE, PUBLIC, RATIOS, RESULTS, load, write_json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


def figures():
    s = json.loads((RESULTS/"summary.json").read_text())
    zero = json.loads((RESULTS/"known_zero_control.json").read_text())
    replay = [json.loads(p.read_text()) for p in sorted((RESULTS/"replay_inertia").glob("*.json"))]
    if len(replay)!=24 or not all(r["pass_all"] for r in replay):
        raise ArithmeticError("The figure needs the complete timed independent replay.")
    timings = Counter()
    counts = {"C":Counter(),"D":Counter()}
    limits = Counter()
    legal_divergence_examples = set()
    raw = load(PUBLIC)
    for rep in replay:
        timings.update(rep["seconds_by_method"])
        for row in rep["signed_state_inertia"]:
            key = ("unverified" if row["positive_negative_zero"] is None else
                   "indefinite" if row["positive_negative_zero"][1]>0 else "positive_definite")
            counts[row["method"]][key] += 1
        for row in rep["confirmed_positive_endpoint_divergence_examples"]:
            _,b,w,slot = row["box_positions"]
            legal_divergence_examples.add((rep["pair"],b,w,slot))
        pair = rep["pair"]
        rec = json.loads((RESULTS/"raw"/f"pair_{pair:02d}.json").read_text())
        for record in rec["rows"]:
            _,b,w,slot=record["box_positions"]
            col = raw["boxes"][pair-1,b,w,slot]
            for task in record["methods"]["E"].get("attempts",[]):
                limits["total_attempts"] += 1
                width = col[1]-col[0]
                unknown = bool(np.any(width[task["local_rows"]]>0.))
                limits["attempts_with_positive_width_rows"] += unknown
                if task["excluded"]:
                    limits["successful_attempts_with_positive_width_rows"] += unknown
                    limits["successful_attempts_with_only_zero_width_rows"] += not unknown
    result = dict(timed_independent_replay_wall_seconds=sum(r["seconds"] for r in replay),
                  independent_replay_wall_seconds_by_method=dict(timings),
                  verified_signed_C_D_state_inertia=counts,
                  E_local_interval_attempts=dict(limits),
                  distinct_demonstrated_legal_endpoint_divergences=len(legal_divergence_examples),
                  matched_known_zero_inertia=zero)
    write_json(RESULTS/"diagnostics_and_costs.json",result)

    plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,
                         "axes.spines.right":False,"font.size":10})
    palette = LinearSegmentedColormap.from_list("confidence",["#f8f3e7","#a6cccc","#368f92","#193a58"])
    fig,axs = plt.subplots(2,5,figsize=(15.2,6.2),sharex=True,sharey=True,layout="constrained")
    labels = ["A  Neumann","B  PSD/QCQP","C  Full inertia","D  Exact trace","E  Shared rows"]
    for i,family in enumerate(("amplitude","structure")):
        eligible = [r for r in s["phase"] if r["family"]==family]
        for j,m in enumerate("ABCDE"):
            ax=axs[i,j]
            values = np.array([r["states"][m]["EXCLUDED"]/r["opposite_noncover_boxes"]
                               if r["opposite_noncover_boxes"] else np.nan for r in eligible]).reshape(4,7)
            heat = ax.pcolormesh(np.arange(8)-.5,np.arange(5)-.5,values,vmin=0,vmax=1,
                                 cmap=palette,shading="flat")
            ax.set_xlim(-.5,6.5)
            ax.set_ylim(3.5,-.5)
            for k,col in enumerate(eligible):
                q,p = divmod(k,7)
                value=values[q,p]
                if np.isfinite(value):
                    ax.text(p,q,str(col["states"][m]["EXCLUDED"]),ha="center",va="center",fontsize=8,
                            color="white" if value>.73 else "#1b3445")
            ax.set_xticks(np.arange(7),[str(r) for r in RATIOS])
            ax.set_yticks(np.arange(4),["0","1/32","1/16","1/8"])
            ax.tick_params(length=0,pad=5,labelsize=8)
            if i==0:
                ax.set_title(labels[j],fontsize=11,loc="left")
            if i==1:
                ax.set_xlabel("Relative width")
            if j==0:
                ax.set_ylabel(("Amplitude pairs" if i==0 else "Structure pairs")+"\nBackground cap")
    fig.suptitle("Strict rejection of opposite boxes that do not contain truth",fontsize=14,
                 fontweight="bold",x=.015,ha="left")
    cbar=fig.colorbar(heat,ax=axs,location="right",shrink=.85,pad=.015)
    cbar.solids.set_rasterized(False)
    cbar.set_label("Fraction of eligible boxes excluded")
    fig.savefig(RESULTS/"phase_diagram.png",dpi=185,facecolor="white")
    fig.savefig(RESULTS/"phase_diagram.pdf",facecolor="white")
    plt.close(fig)

    plt.rcParams.update({"axes.spines.left":False})
    fig,axs = plt.subplots(1,2,figsize=(11.8,4.1),layout="constrained")
    rows = s["raw_run_timing_seconds"]["method_components"]
    shared = s["raw_run_timing_seconds"]["shared"]
    bank = (shared["trace"]/5 + shared["reduced_ball"]+shared["reduced_numeric"])
    colors=["#2e7390","#4ba89c","#d89057","#364763"]
    slots = np.arange(5)
    labels = ["A","B","C","D","E"]
    solve = [rows[m]["proposal"]+rows[m]["reused"]+(rows["D"]["proposal"] if m=="E" else 0.) for m in labels]
    strict = [rows[m]["strict"]+(rows["D"]["strict"] if m=="E" else 0.) for m in labels]
    prep = [0,0,0,bank,bank]
    played = [timings[m]+(timings["D"] if m=="E" else 0.) for m in labels]
    x0 = np.zeros(5)
    for series,color,legend in zip((solve,strict,prep,played),colors,
                                   ("Proposal","Online proof","Trace setup","Independent replay"),strict=True):
        axs[0].barh(slots,series,left=x0,color=color,height=.65,label=legend)
        x0+=series
    axs[0].set_yticks(slots,labels)
    axs[0].invert_yaxis()
    axs[0].set_xlabel("Sum of measured wall seconds")
    axs[0].set_title("Work on 2,688 boxes",loc="left")
    axs[0].legend(ncols=2,frameon=False,fontsize=7,loc="upper right")
    projections = s["projected_wall_amortization_seconds"]
    for role,color in (("C","#2e7390"),("D","#d89057")):
        x=np.array([int(i) for i in projections]); y=[projections[str(i)][role] for i in x]
        axs[1].plot(x,y,marker="o",label=role,color=color)
    axs[1].set_xscale("log")
    axs[1].set_xticks([100,1000,2688],["100","1,000","2,688"])
    axs[1].set_ylabel("Projected wall seconds")
    axs[1].set_xlabel("Boxes with repeated observation reuse")
    axs[1].set_title("Technical amortization of exact trace reuse",loc="left")
    axs[1].legend(frameon=False)
    fig.savefig(RESULTS/"amortization.png",dpi=185,facecolor="white")
    fig.savefig(RESULTS/"amortization.pdf",facecolor="white")
    plt.close(fig)

    truth = load(PRIVATE)["material"]
    pair_ids = [7,8,14]
    fig,axs = plt.subplots(3,2,figsize=(6.4,7.2),layout="constrained")
    vmax=max(float(np.max(truth[p,:,3])) for p in pair_ids)
    for i,pair in enumerate(pair_ids):
        for j,m in enumerate(truth[pair,:,3]):
            ax=axs[i,j]
            heat = ax.pcolormesh(np.arange(6),np.arange(6),m.reshape(5,5),
                                 vmin=0.,vmax=vmax,cmap="cividis",shading="flat")
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"Pair {pair+1} · acquisition {j}",fontsize=9,loc="left")
    fig.suptitle("Same-support generating pairs at background cap 1/8",fontsize=12,
                 fontweight="bold",x=.03,ha="left")
    cbar=fig.colorbar(heat,ax=axs,location="right",shrink=.85,pad=.02)
    cbar.solids.set_rasterized(False)
    cbar.set_label("True contrast")
    fig.savefig(RESULTS/"material_pairs.png",dpi=185,facecolor="white")
    fig.savefig(RESULTS/"material_pairs.pdf",facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    figures()
