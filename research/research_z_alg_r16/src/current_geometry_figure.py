"""Plot archived current-capture geometry and independently scored materials."""
from collections import Counter
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"results"

BLUE = "#2d7490"
ORANGE = "#a9651f"
GREEN = "#388866"
GRAY = "#777779"
STAGES = {
    "background": ("Fixed background", BLUE, "o"),
    "voxel_fit": ("After voxel fit", ORANGE, "s"),
}
COHORTS = {
    "r10_twod": ("R10 (24)", BLUE, "o"),
    "r11_independent48": ("R11 independent (48)", GREEN, "^"),
    "r12_confirm24": ("R12 confirmation (24)", ORANGE, "s"),
}
COVER_NAMES = {
    "fixed48": "Fixed 48 dark",
    "optimal48_dark": "Constitutive-optimal 48 dark\n= all-dark ceiling",
    "updated48_soft": "Step-adapted 48\nwith bright allowed",
}


def read_geometries():
    records = json.loads((RESULTS/"current_geometry.json").read_text())["rows"]
    by_stage = {s:{(r["cohort"],r["case"]):r for r in records if r["stage"]==s}
                for s in STAGES}
    keys = set(by_stage["background"])
    assert len(keys)==56 and set(by_stage["voxel_fit"])==keys
    omitted = {s:[k for k,r in values.items() if r["capture"]["data_gradient"] is None]
               for s,values in by_stage.items()}
    assert not omitted["background"] and len(omitted["voxel_fit"])==1
    coincidence = max(
        abs(c["all_dark_ceiling"]-c["optimal48_dark"]) for r in records
        if (c:=r["capture"]["data_gradient"]) is not None)
    assert coincidence < 2e-12
    return by_stage,omitted,coincidence


def capture_panel(ax, geometry):
    percentiles = {}
    for stage_idx,(stage,(label,color,marker)) in enumerate(STAGES.items()):
        offset = (.13,-.13)[stage_idx]
        rows = [r["capture"]["data_gradient"] for r in geometry[stage].values()
                if r["capture"]["data_gradient"] is not None]
        stage_summary = {}
        for i,key in enumerate(COVER_NAMES):
            values = np.array([r[key] for r in rows])
            assert np.min(values)>-1e-12 and np.max(values)<1+1e-12
            values = np.clip(values,0,1)
            lo,mid,hi = np.percentile(values,(1,50,99))
            stage_summary[key] = [float(p) for p in (lo,mid,hi)]
            y = len(COVER_NAMES)-1-i+offset
            ax.hlines(y,lo,hi,color=color,lw=2.4,zorder=2)
            ax.vlines([lo,hi],y-.035,y+.035,color=color,lw=1.2,zorder=2)
            ax.scatter(mid,y,color=color,marker=marker,s=55,
                       label=f"{label} ({len(rows)})" if i==0 else None,
                       edgecolor="white",linewidth=.75,zorder=3)
            if key!="updated48_soft":
                y_label = y+(.13,-.13)[stage_idx]
                ax.text(mid,y_label,f"{100*mid:.1f}%",ha="center",va="center",fontsize=7.5,color=color)
        percentiles[stage] = stage_summary
    ax.axvline(1.,color="#9c9c9c",linestyle=":",linewidth=.85,zorder=1)
    ax.text(.977,-.42,"100% on the selected direction",ha="right",va="center",fontsize=7.5,color=GRAY)
    ax.set(title="(a) Fraction of the material tangent captured",
           xlabel=r"Energy retained, $\|\Pi_{\mathcal{A}}t\|_F^2/\|t\|_F^2$",
           xlim=(-.025,1.04),ylim=(-.54,2.57),yticks=np.arange(len(COVER_NAMES))[::-1],
           yticklabels=list(COVER_NAMES.values()))
    ax.set_xticks(np.arange(0,1.01,.2))
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.tick_params(axis="y",length=0,pad=5)
    ax.legend(loc="center left",bbox_to_anchor=(.005,-.27),fontsize=7.5,
              ncol=2,frameon=False,columnspacing=1,handletextpad=.3)
    return percentiles


def spectrum_panel(ax, geometry):
    med_counts = {}
    rank = np.arange(1,577)
    for stage,(label,color,_) in STAGES.items():
        rows = list(geometry[stage].values())
        eig = np.stack([r["eigenvalues"][::-1] for r in rows])
        assert eig.shape==(56,576) and eig.min()>-1e-11 and eig.max()<1+1e-11
        lo,med,hi = np.percentile(np.clip(eig,0,1),(1,50,99),axis=0)
        if stage=="voxel_fit":
            ax.fill_between(rank,lo,hi,color=color,alpha=.19,lw=0.)
        ax.plot(rank,med,color=color,lw=1.8,ls="--" if stage=="background" else "-",
                label=label+(" · 1st–99th pct. band" if stage=="voxel_fit" else ""))
        med_counts[stage] = {name:float(np.median([r[name] for r in rows]))
                             for name in ("above_half","above_nine_tenths","top","trace")}
    ax.axhline(.5,color="#bcbcbc",linestyle=":",linewidth=.8,zorder=0)
    ax.axhline(.9,color="#bcbcbc",linestyle=":",linewidth=.8,zorder=0)
    ax.set(title="(b) Fixed-space material geometry",
           xlabel="Eigenvalue rank (descending)",ylabel=r"Eigenvalue of $K$",
           xscale="log",xlim=(1,576),ylim=(0,1.02))
    ax.set_xticks([1,4,16,48,144,576],labels=["1","4","16","48","144","576"])
    ax.legend(loc="upper right",fontsize=7.3,frameon=False,handlelength=1.8)
    ax.text(.04,.08,"After voxel fit\nMedian 44 modes ≥ 0.50, 6 ≥ 0.90",fontsize=7.5,color=GRAY,
            transform=ax.transAxes,ha="left",va="bottom")
    return med_counts


def quality_panel(ax):
    fixed,adapted = "shared48","shared36_tangent_soft"
    with (RESULTS/"current_cases.csv").open() as stream:
        selected = [r for r in csv.DictReader(stream)
                    if r["cohort"] in COHORTS and r["arm"] in (fixed,adapted)]
    by_arm = {a:{(r["cohort"],r["case"]):r for r in selected if r["arm"]==a} for a in (fixed,adapted)}
    keys = sorted(by_arm[fixed])
    assert len(keys)==96 and set(by_arm[adapted])==set(keys) and len(selected)==192
    nx,ny = ([float(by_arm[arm][k]["amplitude_NRMSE"]) for k in keys] for arm in (fixed,adapted))
    nx,ny = np.asarray(nx),np.asarray(ny)
    xmax = 1.055*max(nx.max(),ny.max())
    ax.plot([0,xmax],[0,xmax],color="#9c9c9c",lw=.9,linestyle="--",zorder=0)
    counts = Counter(key[0] for key in keys)
    for group,(label,color,marker) in COHORTS.items():
        mask = np.array([k[0]==group for k in keys])
        ax.scatter(nx[mask],ny[mask],marker=marker,color=color,alpha=.72,s=18,linewidth=.25,
                   edgecolor="white",label=label,zorder=2)
    ax.set(title="(c) Final material fidelity on 96 paired cases",
           xlabel="Fixed 48 · amplitude NRMSE",ylabel="Step-adapted 48 · amplitude NRMSE",
           xlim=(0,xmax),ylim=(0,xmax))
    ax.set_aspect("equal",adjustable="box")
    ax.legend(loc="upper left",fontsize=7.2,frameon=False,handletextpad=.25,borderpad=.15)
    successful = {arm:sum(by_arm[arm][k]["material_success"]=="True" for k in keys) for arm in (fixed,adapted)}
    ax.text(.98,.055,f"Material successes\nFixed {successful[fixed]}/96 · Adapted {successful[adapted]}/96",
            transform=ax.transAxes,ha="right",va="bottom",fontsize=7.6,color=GRAY)
    return dict(cases=dict(counts),median_fixed=float(np.median(nx)),median_adapted=float(np.median(ny)),
                median_paired_change=float(np.median(ny-nx)),successful=successful)


def run():
    plt.rcParams.update({"font.size":9,"axes.titlesize":9.5,"axes.labelsize":9})
    data,omitted,coincidence=read_geometries()
    fig,axes = plt.subplots(1,3,figsize=(14.4,3.8),layout="constrained",
                            gridspec_kw={"width_ratios":[1.35,1.10,1.14]})
    percentiles = capture_panel(axes[0],data)
    med_counts = spectrum_panel(axes[1],data)
    quality = quality_panel(axes[2])
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
    fig.supxlabel("Coverage intervals show 1st–99th percentiles. Of 56 paired observations, one fitted feasible "
                  "gradient is zero and omitted in (a).",fontsize=8,color="#4a4a4c")
    for extension in ("png","pdf"):
        fig.savefig(ROOT/"figures"/f"current_geometry.{extension}",dpi=220,bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(dict(intervals_and_median=percentiles,spectral_medians=med_counts,
                          omitted=omitted,equal_dark_within=coincidence,quality=quality),indent=2))


if __name__=="__main__":
    run()
