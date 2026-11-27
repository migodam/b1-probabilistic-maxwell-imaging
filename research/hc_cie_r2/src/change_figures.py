"""Plot sealed grouped outcomes, physical exclusion bounds and position/amplitude examples."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LAB=Path(__file__).resolve().parents[1]; OUT=LAB/"figures"
COL=["#385988","#cd6e48","#2b796a","#8c679b"]


def read(path):
    with (LAB/"results"/path).open() as f: return list(csv.DictReader(f))


def panel(rows,scores,ax,labels,noise,methods,norm):
    for j,m in enumerate(methods):
        nums=[]
        for cond in noise:
            z=[r for r in rows if (r.get("condition") or r.get("noise"))==cond and
               (r.get("mode") or r.get("method"))==m]
            nums.append(sum(int(r["exact_support"]) for r in z))
        t=np.arange(len(noise))-.34+(j+.5)*.68/len(methods)
        bars=ax.bar(t,nums,.67/len(methods),color=COL[j],label=scores[j])
        ax.bar_label(bars,fontsize=7,padding=1)
    ax.set_xticks(np.arange(len(noise)),labels);ax.set_ylim(0,norm*1.14)
    ax.set_ylabel("Exact two-site supports");ax.spines[["right","top"]].set_visible(False)
    ax.legend(fontsize=7,loc="upper right",bbox_to_anchor=(1.01,-.19),ncol=2)


def main():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titlelocation":"left"})
    r=read("interaction/offline.csv");p=read("support/monitor_offline.csv")
    x=read("cross_changes/offline.csv");vis=read("interaction/port_visibility.csv")
    fig,axes=plt.subplots(2,2,figsize=(11,8.8),layout="constrained")
    panel(r+[t for t in p if t["mode"]=="conditional_readout"],
          ["MUSIC + fit","Diagonal + fit","Interaction + fit","Interaction + readout"],axes[0,0],
          ["Clean","10% of change","10% of total"],["clean","ten_percent_change","ten_percent_total"],
          ["music","diagonal","interaction","conditional_readout"],48)
    axes[0,0].set_title("A   Known background, 24 materials × 2 changes")
    panel(x,["MUSIC + readout","Interaction + readout","Interaction + fit"],axes[0,1],
          ["32 clean","32 noisy","96 → 32\nclean","96 → 32\nnoisy"],
          ["matched_clean","matched_noisy","cross_clean","cross_noisy"],
          ["music_closed","interaction_closed","interaction_fit"],24)
    axes[0,1].set_title("B   Six archived families, paired acquisition layouts")
    for j,dist in enumerate(["adjacent","separated"]):
        vs=sorted(float(t["relative_cross"])*100 for t in vis if t["noise"]=="clean" and t["distance"]==dist)
        axes[1,0].plot(range(1,len(vs)+1),vs,"o-",ms=3,color=COL[j],label=dist.capitalize())
    ax=axes[1,0];ax.axhline(10,ls="--",c=".35",label="10% differential-noise radius")
    ax.set_yscale("log");ax.set_xlabel("Ordered parent within each change type")
    ax.set_ylabel("Distinguishable interaction / change field (%)")
    ax.set_title("C   Joint scattering after removing single-site signals");ax.legend(fontsize=8)
    ax=axes[1,1];q=read("support/all_pairs.csv")
    ids=[f"{int(t['parent']):04d}" for t in q];t=np.arange(len(ids))
    ax.scatter(t,[float(v["rank2_by_delta"]) for v in q],s=32,color=COL[0],label="Any rank-two port response")
    ax.scatter(t,[float(v["span_by_delta"]) for v in q],s=30,marker="D",color=COL[1],label="Best of all 165,600 symmetric pairs")
    ax.axhline(1,ls="--",c=".35");ax.set_yscale("log");ax.set_xticks(t,ids,rotation=35);ax.set_xlabel("Historical Dual-CIE endpoint")
    ax.set_ylabel("Minimum distance / noise radius");ax.legend(fontsize=7);ax.set_title("D   Which two-pixel corrections can reach the data?")
    fig.savefig(OUT/"material_changes.png",dpi=180);plt.close(fig)
    fig,axes=plt.subplots(4,4,figsize=(11,10),layout="constrained")
    # Same parent and same change under both acquisition layouts; example chosen by ID, not score.
    for i,(parent,distance) in enumerate([("s2d_001","adjacent"),("s2d_001","separated"),
            ("s2d_004","adjacent"),("s2d_004","separated")]):
        files=LAB/"results/cross_changes"
        base=np.load(files/"public"/(parent+"__pw360_10x20__"+distance+".npz"))["known"].reshape(32,32)
        q=np.load(files/"scoring"/(parent+"__pw360_10x20__"+distance+".npz"));m=np.zeros(1024);m[q["pair"]]=q["steps"]
        panels=[base,m.reshape(32,32)]
        for layout in ["pw360_10x20","ls_line_same_10x20"]:
            row=next(t for t in x if t["parent"]==parent and t["distance"]==distance and t["layout"]==layout and t["condition"]=="cross_noisy" and t["method"]=="interaction_closed")
            pred=np.zeros(1024);pred[[int(row["i"]),int(row["j"])]]=[float(row["di"]),float(row["dj"])]
            panels.append(pred.reshape(32,32))
        for j,img in enumerate(panels):
            a=axes[i,j];im=a.imshow(img,vmin=0,vmax=8,cmap="viridis",interpolation="nearest")
            a.set_xticks([]);a.set_yticks([])
            if i==0:a.set_title(["Known reference","True change","96 → 32, full ring","96 → 32, same-side line"][j])
            if j==0:a.set_ylabel(parent+"\n"+distance)
    fig.colorbar(im,ax=axes,shrink=.6,label="Real contrast (shared scale; changes in columns 2–4)")
    fig.savefig(OUT/"cross_changes_maps.png",dpi=160);plt.close(fig)


if __name__=="__main__":main()
