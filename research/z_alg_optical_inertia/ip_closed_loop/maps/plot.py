"""Reader figures from archived score arrays."""
import json
from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from common import HERE

RES=HERE/"results"
FIG=HERE/"figures"
CMAP=ListedColormap(["#17675f","#ff613d","#ffdf93"])
STAT=[Patch(facecolor=CMAP.colors[0],label="Excluded"),
      Patch(facecolor=CMAP.colors[1],label="Physical witness"),
      Patch(facecolor=CMAP.colors[2],label="Unresolved"),
      Line2D([0],[0],color="#262d39",marker="o",markerfacecolor="none",markersize=7,linewidth=0,label="True centre")]


def read(path):
    return np.load(path,allow_pickle=False)


def tidy(ax,n):
    ax.set_xlim(-.5,n-.5)
    ax.set_ylim(-.5,n-.5)
    if n==32:
        ax.set_xticks([0,15,31])
        ax.set_yticks([0,15,31])
    else:
        ax.set_xticks([0,47,95])
        ax.set_yticks([0,47,95])
    ax.tick_params(length=1.5,labelsize=8,width=.45)
    for spine in ax.spines.values():
        spine.set_linewidth(.5)


def display_state(ax,arr,truth=None):
    ax.imshow(arr,origin="lower",cmap=CMAP,vmin=-.5,vmax=2.5,interpolation="nearest")
    if truth is not None:
        ax.plot(*truth,"o",color="#292831",markerfacecolor="none",markersize=6,mew=1.05)
    tidy(ax,len(arr))


def out(fig,names):
    for ext in ("png","pdf"):
        fig.savefig(FIG/f"{names}.{ext}",dpi=185,facecolor="white",bbox_inches="tight")
    plt.close(fig)


def atlas(rows):
    fig,ax=plt.subplots(4,4,figsize=(9.3,8.6))
    for a,row in zip(ax.ravel(),rows,strict=True):
        oid=row["obs_id"]
        shapes=read(RES/f"{oid}__b0.03125__library_all_union.npz")
        cover=shapes["status"].copy()
        if "witness_position" in row:
            x,y=row["witness_position"]
            assert cover[y,x]!=0
            cover[y,x]=1
        display_state(a,cover,row["true_centre"])
        a.set_title(oid.split("__")[1]+f"  {np.count_nonzero(cover==0)}/1024",fontsize=9,pad=3)
    fig.suptitle("Full field of view with a four-shape prior",fontsize=14,fontweight="bold",y=.98)
    fig.text(.5,.935,"32 × 32 grid  |  same 7 illuminations and 9 receivers  |  3% noise  |  actual weak background",
             ha="center",fontsize=9)
    fig.legend(handles=STAT,loc="lower center",bbox_to_anchor=(.5,.006),ncol=4,frameon=False)
    fig.subplots_adjust(top=.9,bottom=.07,left=.07,right=.97,hspace=.4,wspace=.26)
    out(fig,"weak-background-library")


def native(rows):
    fig,axes=plt.subplots(2,4,figsize=(11.1,5.55),gridspec_kw={"width_ratios":[1,1,1,1]})
    im=None
    for axs,row in zip(axes,rows,strict=True):
        oid=row["obs_id"]
        n=row["grid"]
        truth=read(HERE/"offline"/f"{oid}.npz")["truth"].reshape(n,n)
        ctl=read(RES/f"{oid}__b0.03125__primary.npz")
        labels=read(RES/f"{oid}__b0.03125__score_map.npz")["status"]
        wit=read(RES/f"{oid}__b0.03125__witness.npz")["material"].reshape(n,n)
        xy=row["true_centre"]
        im=axs[0].imshow(truth,origin="lower",cmap="magma",vmin=0,vmax=4.25,interpolation="nearest")
        axs[0].set_title("Generating material",fontsize=10)
        axs[1].imshow(ctl["osm"]/ctl["osm"].max(),origin="lower",cmap="Blues",
                       vmin=0,vmax=1,interpolation="nearest")
        axs[1].plot(*xy,"+",color="#a82a21",mew=1.15,markersize=7)
        axs[1].set_title("Backpropagation",fontsize=10)
        display_state(axs[2],labels,xy)
        axs[2].set_title(f"Certified states  {row['excluded']}/9216",fontsize=10)
        axs[3].imshow(wit,origin="lower",cmap="magma",vmin=0,vmax=4.25,interpolation="nearest")
        axs[3].set_title(f"Physical witness  {row['witness_residual_ratio']:.3f} r",fontsize=10)
        axs[0].set_ylabel(oid.split("__")[1],fontsize=10)
        for a in axs:
            tidy(a,n)
    fig.suptitle("Native 96 × 96 full-field maps under actual weak background",
                 fontsize=14,fontweight="bold",y=.982)
    fig.legend(handles=STAT,loc="lower center",bbox_to_anchor=(.49,.005),ncol=4,frameon=False)
    bar_ax=fig.add_axes([.92,.22,.012,.56])
    c=fig.colorbar(im,cax=bar_ax)
    c.set_label("Permittivity contrast",fontsize=8)
    c.ax.tick_params(labelsize=8)
    fig.subplots_adjust(top=.88,bottom=.12,left=.08,right=.89,hspace=.29,wspace=.28)
    out(fig,"native-weak-background")


def broad(rows):
    fig,axes=plt.subplots(1,len(rows),figsize=(2.3*len(rows),3))
    for a,row in zip(axes,rows,strict=True):
        oid,n=row["obs_id"],row["grid"]
        r=read(RES/f"{oid}__b0.03125__broad_sector.npz")
        m=r["margin_sector"].reshape(n,n)
        state=np.where(m>0,0,2)
        if "witness_position" in row:
            x,y=row["witness_position"]
            assert state[y,x]!=0
            state[y,x]=1
        n_nu=int(np.count_nonzero((m>0)&(r["negative_index"].reshape(m.shape)>0)))
        assert n_nu==np.count_nonzero(m>0)
        display_state(a,state,row["true_centre"])
        a.set_title(oid.split("__")[1]+f"  {n_nu}/{n*n}",fontsize=9,pad=3)
    fig.suptitle("Broad annulus certificate  |  core [3, 5], free shell [0, 5]",
                 fontsize=13,fontweight="bold",y=.975)
    fig.text(.5,.855,"Unknown exterior [0, 1/32]. Fractions show exclusions certified with negative state inertia.",
             ha="center",fontsize=8.5)
    fig.legend(handles=STAT,loc="lower center",bbox_to_anchor=(.5,.005),ncol=4,frameon=False)
    fig.subplots_adjust(top=.76,bottom=.17,left=.04,right=.98,wspace=.35)
    out(fig,"broad-annulus")


def phase():
    files=sorted((RES/"stress").glob("stress_*.json"))
    if not files:
        return
    grouped={}
    for path in files:
        entries=list(json.loads(path.read_text())["rows"].values())
        for row in entries:
            if row["n_centres"]!=1024:
                continue
            b=row["inference_background"]
            noise=row["noise"]
            grouped.setdefault((b,noise),[]).append(row["excluded"]/1024)
    if len(grouped)!=6:
        return
    fig,ax=plt.subplots(figsize=(6.4,2.7))
    ns=(.01,.03,.10)
    labels=(0.,.03125)
    values=np.array([[np.mean(grouped[(b,p)]) for p in ns] for b in labels])
    im=ax.imshow(values,vmin=.96,vmax=1,cmap="GnBu",aspect="auto")
    for y,x in np.ndindex(values.shape):
        v=grouped[labels[y],ns[x]]
        ax.text(x,y,f"{100*values[y,x]:.2f}%\n{min(v)*100:.2f}–{max(v)*100:.2f}%",
                ha="center",va="center",color="white" if values[y,x]>.985 else "#18333e",fontsize=9)
    ax.set_xticks(range(3),["1%","3%","10%"])
    ax.set_yticks(range(2),["Actual vacuum, b = 0","Actual weak background, b = 1/32"])
    ax.set_xlabel("Relative joint complex noise")
    ax.set_title("Strict rejection over every centre  |  four paired 32 × 32 objects",
                 fontsize=11,fontweight="bold",pad=9)
    fig.colorbar(im,ax=ax,shrink=.76).ax.tick_params(labelsize=8)
    out(fig,"noise-background-phase")


def main():
    FIG.mkdir(exist_ok=True)
    summary=json.loads((RES/"summary.json").read_text())
    rows=[r for r in summary["rows"] if r["actual_scene"]=="weak"]
    atlas([r for r in rows if r["grid"]==32])
    native([r for r in rows if r["grid"]==96])
    broad([r for r in rows if r["shape"]=="ring"])
    phase()


if __name__=="__main__":
    main()
