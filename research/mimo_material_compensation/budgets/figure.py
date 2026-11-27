"""Render the budget comparison from the raw receipts."""
from pathlib import Path
import json
import os

os.environ["MPLBACKEND"]="Agg"

import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE
COLORS={"A":"#176b87","B":"#9b362d"}
PHASE=["hideable","unresolved","unhideable"]
PHASE_COLORS=["#2c7f92","#f4e9d9","#a54636"]


def bundle(key):
    return json.loads((HERE/f"budget_{key}.json").read_text())


def scale_sort(rows,lay,contract,t,kind="max"):
    return sorted([r for r in rows if (r["layout"],r["contract"],r["t"],r["kind"])==(lay,contract,t,kind)],
                  key=lambda q:q["beta"])


def plot_phase(ax, rows, lay):
    from matplotlib.colors import ListedColormap
    t_values=[.22,.23,.24,.25,.252,.255,.26,.27,.28]
    b_values=[0.,.0001,.0003,.0005,.001,.0015,.002,.0025,.003]
    s={(r["t"],r["beta"]):r["status"] for r in rows if (r["layout"],r["contract"],r["kind"])==(lay,"prospective","max")}
    phase=np.array([[PHASE.index(s[t,b]) if (t,b) in s else np.nan for b in b_values] for t in t_values])
    cmap=ListedColormap(PHASE_COLORS)
    cmap.set_bad(color="#ececec")
    ax.imshow(phase,origin="lower",aspect="auto",cmap=cmap,vmin=0,vmax=2)
    ax.set_yticks(range(len(t_values)),[str(t) for t in t_values])
    ax.set_xticks(range(len(b_values)),[f"{b*1000:g}" for b in b_values],rotation=35)
    ax.set_xlabel(r"Exterior maximum budget $\beta_\infty$ ($10^{-3}$)")
    ax.set_title(f"{lay}   |   Fixed-reference phase")
    for a,spine in ax.spines.items():
        spine.set_visible(False)


def main():
    plt.rcParams.update({"font.size":10,"legend.fontsize":8,"axes.titlesize":11,"axes.labelsize":9,
                          "pdf.fonttype":42,"ps.fonttype":42,"figure.dpi":130})
    joint=bundle("joint32")
    focus=bundle("focus_joint32")
    actual=joint["records"]
    arr=joint["records"]+focus["records"]
    fig,axs=plt.subplots(2,2,figsize=(10.5,7.1),layout="constrained")
    ax=axs[0,0]
    x=np.array([r["beta"] for r in scale_sort(arr,"B","prospective",.27)])*1000
    b=scale_sort(arr,"B","prospective",.27)
    a=scale_sort(arr,"A","prospective",.27)
    ax.plot(x,[r["lower_over_threshold"] for r in b],color=COLORS["B"],marker="o",
            markersize=4,label="B directional lower bound")
    ax.plot(x,[r["empirical_direction_lower_over_threshold"] for r in b],color=COLORS["B"],ls="--",
            label="B empirical-direction lower bound")
    ax.plot(x,[r["norm_control_lower_over_threshold"] for r in b],color="#a96e61",ls=":",
            label="B norm-subtraction lower bound")
    ax.plot(x,[r["witness_upper_over_threshold"] for r in a],color=COLORS["A"],marker="s",
            markersize=4,label="A witness upper bound")
    ax.axhline(1.,color="black",lw=.8)
    ax.axvspan(1.,3.,color="#e4f1f3",zorder=0)
    ax.set_ylim(0,1.38)
    ax.set_title(r"$32^2$ joint-fit reference, $t=0.27$")
    ax.set_ylabel(r"Separation relative to $2r$")
    ax.set_xlabel(r"Exterior maximum budget $\beta_\infty$ ($10^{-3}$)")
    ax.legend(loc="lower left",frameon=False)
    ax=axs[0,1]
    b=scale_sort(actual,"B","actual",.045)
    a=scale_sort(actual,"A","actual",.045)
    pairs=[(float(w["beta"])*1000,z,w) for z,w in zip(b,a) if w["beta"]<=.0003]
    ax.plot([p[0] for p in pairs],[p[2]["lower_over_threshold"] for p in pairs],
            color=COLORS["A"],marker="o",label="A directional lower bound")
    ax.plot([p[0] for p in pairs],[p[1]["witness_upper_over_threshold"] for p in pairs],
            color=COLORS["B"],marker="s",label="B witness upper bound")
    ax.axhline(1.,color="black",lw=.8)
    ax.axvspan(0.,.1,color="#f7e8e4",zorder=0)
    ax.set_title(r"Recorded noisy readings, $t=0.045$")
    ax.set_xlabel(r"Exterior maximum budget $\beta_\infty$ ($10^{-3}$)")
    ax.set_ylabel(r"Observed residual relative to $r$")
    ax.legend(loc="upper right",frameon=False)
    for i,a in enumerate(("A","B")):
        plot_phase(axs[1,i],arr,a)
    axs[1,0].set_ylabel(r"ROI contrast increase $t$")
    from matplotlib.patches import Patch
    axs[1,1].legend(handles=[Patch(facecolor=c,label=l.capitalize()) for l,c in zip(PHASE,PHASE_COLORS)]+
                    [Patch(facecolor="#ececec",label="Grid cell absent")],
                    loc="upper left",bbox_to_anchor=(0,-.38),ncols=2,frameon=False,borderaxespad=0)
    fig.savefig(OUT/"budget_phase.pdf",bbox_inches="tight")
    fig.savefig(OUT/"budget_phase.png",dpi=200,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    main()
