"""Plot saved full-wave calibration and interval summaries without forward solves."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def render():
    folder=ROOT/"results"
    before=json.loads((folder/"geometry_mismatch.json").read_text())
    after=json.loads((folder/"geometry_calibrated.json").read_text())
    outer=json.loads((folder/"host_error_bound.json").read_text())
    teal, gold, grey="#227b72", "#edc675", "#a9b4bf"
    blue, rust="#305d83", "#a65432"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,
                         "axes.spines.right":False,"pdf.fonttype":42,"ps.fonttype":42})
    fig=plt.figure(figsize=(8.65,5.65))
    gs=fig.add_gridspec(2,2,left=.085,right=.975,top=.885,bottom=.21,hspace=.37,wspace=.16,
                        height_ratios=[.81,1.14])
    counts=fig.add_subplot(gs[0,:])
    primary=fig.add_subplot(gs[1,0])
    confirmation=fig.add_subplot(gs[1,1],sharey=primary)
    fig.text(.025,.969,"A",fontsize=13,weight="bold",va="center")
    fig.text(.065,.969,"Coverage with a 1% host-contrast error",fontsize=12,weight="bold",va="center")
    entries=[
        ("Nominal host subtraction",before["combined"]["wrong_host_baseline"]),
        ("Actual empty-host subtraction",before["combined"]["actual_host_baseline"]),
        ("Calibrated shots, noise-only interval",after["combined"]["one_calibration"]),
        ("Calibrated shots, host-aware interval",outer["outer"]),
    ]
    # Bars occupy x = 0..32; the negative coordinate space carries the method labels.
    counts.set_xlim(-20.7,32.4)
    counts.set_ylim(-.45,len(entries)-.48)
    for y,(name,rec) in zip(range(len(entries)-1,-1,-1),entries):
        n,covered,empty=rec["n"],rec["covered"],rec["empty"]
        sizes=[covered,n-covered-empty,empty]
        left=0
        counts.text(-.66,y,name,ha="right",va="center",fontsize=9)
        for j,(size,color) in enumerate(zip(sizes,[teal,gold,grey])):
            if size:
                counts.barh(y,size,left=left,height=.67,color=color,linewidth=.7,edgecolor="white")
                counts.text(left+size/2,y,str(size),ha="center",va="center",size=9,
                            color="white" if j==0 else "#242729",weight="bold")
            left+=size
    counts.set_axis_off()
    handles=[Line2D([],[],ls="",marker="s",markersize=7,markeredgewidth=0,color=c,label=l)
             for c,l in zip([teal,gold,grey],["Contains truth","Excludes truth, nonempty","Empty"])]
    fig.legend(handles=handles,loc="upper left",bbox_to_anchor=(.361,.931),ncol=3,frameon=False,
               handletextpad=.35,columnspacing=1.15,fontsize=8.2)
    for ax,key,letter,title in [(primary,"primary","B","Primary fixture"),(confirmation,"confirmation","C","Confirmation fixture")]:
        sample=sorted([row for row in outer["rows"] if row["layout"]==key],key=lambda row:row["case"])
        x=np.array([row["case"]+1 for row in sample])
        nominal=np.array([row["nominal"]["width"] for row in sample])
        robust=np.array([row["outer"]["width"] for row in sample])
        ax.vlines(x,nominal,robust,color="#b6c1c6",lw=1.35,zorder=1)
        ax.scatter(x,nominal,color=blue,s=22,edgecolors="white",linewidths=.4,zorder=3)
        ax.scatter(x,robust,color=rust,marker="D",s=24,edgecolors="white",linewidths=.4,zorder=3)
        ax.set_title(f"{letter}   {title}",loc="left",fontsize=11.2,weight="bold",pad=7)
        ax.set_yscale("log")
        ax.set(xlim=(.45,16.55),ylim=(.025,17),xlabel="Target case")
        ax.set_xticks([1,4,8,12,16])
        ax.set_yticks([.03,.1,.3,1,3,10],["0.03","0.1","0.3","1","3","10"])
        ax.minorticks_off()
        ax.tick_params(labelsize=8.3)
        ax.grid(axis="y",color="#dcdfe1",lw=.65,zorder=0)
        ax.axhline(14,ls=(0,(3,2)),color="#727b80",lw=.9)
    primary.set_ylabel("Hull width in material units",fontsize=9)
    confirmation.tick_params(axis="y",labelleft=False)
    legend=[Line2D([],[],color=blue,ls="",marker="o",ms=5,label="Noise-only interval"),
            Line2D([],[],color=rust,ls="",marker="D",ms=5,label="Host-aware outer interval"),
            Line2D([],[],color="#727b80",ls=(0,(3,2)),label="Full material box (width 14)")]
    fig.legend(handles=legend,loc="lower center",bbox_to_anchor=(.52,.063),ncol=3,frameon=False,
               columnspacing=1.35,handletextpad=.45,fontsize=8.4)
    fig.text(.5,.029,"Three target shots per material (total energy 2), plus one calibration shot per fixture (energy 1) shared by 16 materials.",
             ha="center",va="center",fontsize=8.25,color="#293139")
    fig.savefig(ROOT/"figures/host_calibration.pdf")
    fig.savefig(ROOT/"figures/host_calibration.png",dpi=260)
    plt.close(fig)


if __name__=="__main__":
    render()
