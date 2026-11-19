"""Paper plots for matched-data stationarity and the mirror full-box control."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"results"
FIGURES = ROOT/"figures"

BLUE = "#176b9e"
ORANGE = "#c65a28"
GREEN = "#34856a"
GRAY = "#575759"


def read(name):
    return json.loads((RESULTS/f"geometry-{name}.json").read_text())


def save(fig,name):
    for extension in ("png","pdf"):
        fig.savefig(FIGURES/f"geometry_{name}.{extension}",dpi=190,bbox_inches="tight")
    plt.close(fig)


def stationarity():
    data = read("fold")
    known, shared = data[1],data[3]
    fig, axes = plt.subplots(1,3,figsize=(11.5,3.4),layout="constrained")

    ax = axes[0]
    good = known["start_candidates"][1]
    bad = known["start_candidates"][0]
    ax.scatter([good["material_rms"],bad["material_rms"]],[good["score"],bad["score"]],
               facecolors="white",edgecolors=[BLUE,ORANGE],marker="o",linewidths=1.7,s=58,zorder=3)
    bad_end = known["fits"][known["certified_wrong_minima"][0]["fit_index"]]
    good_end = known["fits"][known["multistart"]["selected_index"]]
    for start,end,color in ((good,good_end,BLUE),(bad,bad_end,ORANGE)):
        a = (start["material_rms"],start["score"])
        b = (end["material_rms"],end["score"])
        ax.scatter(*b,color=color,s=31,zorder=4)
        ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="->",color=color,shrinkA=5,shrinkB=5))
    ax.axhline(1,color=GRAY,ls="--",lw=.9)
    ax.annotate("noise contract",(3.5e-4,1.22),fontsize=8,color=GRAY)
    ax.annotate("strict local\nminimum",(.6,45),(.024,19),
                color=ORANGE,fontsize=8,arrowprops=dict(arrowstyle="-",color=ORANGE))
    ax.set(xscale="log",yscale="log",xlabel="Material RMS error",ylabel="Residual / noise radius",
           title="(a) Same new read, two starts",ylim=(.34,115))

    ax = axes[1]
    for row,key,label,color,style in (
        (known,"r13_eta","R13 absolute",ORANGE,"--"),
        (known,"relative_eta","Scaled, calibrated",BLUE,"-"),
        (shared,"relative_eta","Scaled, shared gain",GREEN,"-."),
    ):
        rows = [r for r in row["certificates"]["domain_checks"] if r[key] is not None]
        ax.plot([r["material_half_width"] for r in rows],[r[key] for r in rows],
                color=color,ls=style,marker="o",markersize=3.2,label=label)
    ax.axhline(1/np.sqrt(2),color=GRAY,ls=":",lw=1,label=r"Positive bound if $\eta<1/\sqrt{2}$")
    ax.set(xscale="log",yscale="log",xlim=(.0008,.4),ylim=(.012,80),
           xlabel="Material box half-width",ylabel=r"Upper bound on $\eta$",title="(b) Certified direction range")
    ax.legend(frameon=False,fontsize=6.8,loc="lower right")

    ax = axes[2]
    natural = read("natural")
    confirm = read("confirm")
    for group,label,color,marker in ((natural,"Natural, 18",BLUE,"o"),(confirm,"Finite aperture, 8",ORANGE,"s")):
        half = [r["certificates"]["relative_width"] for r in group]
        error = [np.linalg.norm(
            np.array(r["fits"][r["multistart"]["selected_index"]]["material"])-np.array(r["truth"]),
            np.inf) for r in group]
        ax.scatter(half,error,marker=marker,color=color,s=24,alpha=.75,label=label)
    ax.plot([.0025,22],[.0025,22],ls="--",lw=.9,color=GRAY)
    ax.set(xscale="log",yscale="log",xlim=(.0027,.35),ylim=(.0027,23),
           xlabel="Certified material half-width",ylabel=r"Actual center error ($\ell_\infty$)",
           title="(c) Does the box contain truth?")
    ax.legend(frameon=False,fontsize=7,loc="lower right")
    save(fig,"stationarity")


def mirror():
    result = read("mirror-box-gain")
    scored = read("mirror-box-gain-scores")
    names = ("symmetric_broadband","asymmetric_broadband","asymmetric_same_frequency")
    labels = ("Symmetric, broadband","Asymmetric, broadband","Asymmetric, one frequency")
    colors = ("#8795a1",BLUE,ORANGE)
    fig, axes = plt.subplots(1,3,figsize=(11.1,3.4),layout="constrained")
    with np.load(RESULTS/result["cells_file"],allow_pickle=False) as inputs:
        for i,ax in enumerate(axes[:2]):
            for key,color in zip(names,colors,strict=True):
                if i==0 and key==names[0]:
                    boxes = inputs[key]
                    x,y = boxes.mean(axis=1).T
                    ax.scatter(x,y,s=2,color=color,alpha=.23,rasterized=True)
                elif key!=names[0]:
                    lo,hi = result["arms"][key]["hull"]
                    ax.add_patch(Rectangle(lo,*(np.asarray(hi)-lo),fill=False,ec=color,
                                           lw=1.6,linestyle="--" if key==names[2] else "-"))
            truth,other = np.asarray(scored["truth"]),np.asarray(scored["other"])
            ax.scatter(*truth,marker="+",color="black",s=80,zorder=5)
            if i==0:
                ax.scatter(*other,marker="x",color="black",s=45,zorder=5)
                low,hi = (9.89,11.17)
                ax.plot((low,hi),(low,hi),ls=":",color=GRAY,lw=.8)
                ax.set(xlim=(low,hi),ylim=(low,hi),title="(a) Full-box survivors")
            else:
                ax.set(xlim=(10.904,10.945),ylim=(10.052,10.097),title="(b) Asymmetric hulls")
            ax.set_aspect("equal")
            ax.set(xlabel="Material 1",ylabel="Material 2")
            ax.tick_params(axis="both",labelsize=8)
    ax=axes[2]
    diam = [result["arms"][k]["rms_diameter_upper"] for k in names]
    bars=ax.barh(np.arange(3),diam,color=colors,height=.55)
    ax.axvline(scored["pair_rms"],ls=":",color="black",lw=1)
    for bar,value in zip(bars,diam,strict=True):
        ax.annotate(f"{value:.3f}",(value,bar.get_y()+bar.get_height()/2),xytext=(4,0),
                    textcoords="offset points",va="center",fontsize=8)
    ax.set(xscale="log",xlim=(.009,2.5),xlabel="Certified RMS diameter upper bound",
           yticks=np.arange(3),yticklabels=labels,title="(c) Same energy and shared gain")
    ax.tick_params(axis="y",labelsize=7.5)
    ax.invert_yaxis()
    fig.legend(handles=[Patch(color=colors[i],label=labels[i]) for i in range(3)],
               loc="outside lower center",ncol=3,frameon=False,fontsize=8)
    save(fig,"mirror")


if __name__=="__main__":
    plt.rcParams.update({"axes.spines.top":False,"axes.spines.right":False,"font.size":9,
                         "axes.titlesize":9,"axes.labelsize":9})
    stationarity()
    mirror()
