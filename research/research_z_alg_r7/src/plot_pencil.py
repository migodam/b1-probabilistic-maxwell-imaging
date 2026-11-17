"""Summarize and plot the physical patch pencil comparisons."""
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import numpy as np
from numpy.polynomial import polynomial as poly

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/"results"
FIGURES=ROOT/"figures"
COLOR=dict(pencil="#245e85", anchor="#b75539", random="#c58d2b", completion="#2a8371",
           interval="#6e567f")
LABEL=dict(pencil="Coordinated", anchor="Local anchors", random="Random span",
           completion="All ten feeds", interval="Coordinated interval")
LAYOUTS=("primary","confirmation")
REGIONS=("L","square")
ARMS=("pencil","anchor","random","completion")
NOISE=(0.,2e-5,2e-4)
LAST=NOISE[-1]


def descr(values):
    v = np.asarray([a for a in values if a is not None])
    if v.size==0:
        return None
    return dict(n=v.size,min=float(v.min()),q25=float(np.quantile(v,.25)),
                median=float(np.median(v)),q75=float(np.quantile(v,.75)),max=float(v.max()))


def summarize(raw):
    groups=defaultdict(list)
    for row in raw["observations"]:
        groups[row["layout"],row["region"],row["epsilon"],row["arm"]].append(row)
    combined=[]
    for (layout,region,noise,arm),rows in groups.items():
        iv=[r["interval"] for r in rows if r["interval"] is not None]
        err=[abs(r["trf"]["best"]["x"][0]-r["truth"]) for r in rows]
        wit=[r["witness"]["spread"] for r in rows if r["witness"] is not None]
        combined.append(dict(
            layout=layout,region=region,noise=noise,arm=arm,n=len(rows),
            trf_error=descr(err),trf_cpu_s=descr([r["trf"]["cpu_s"] for r in rows]),
            trf_reduced_lu=descr([r["trf"]["reduced_lu"] for r in rows]),
            trf_fs_count=sum(r["trf"]["best"]["residual_l2"]<=r["noise_radius_with_rounding"] for r in rows),
            trf_start_evaluation_limit=sum(u["status"]==0 for r in rows for u in r["trf"]["starts"]),
            trf_residual=descr([r["trf"]["best"]["residual_l2"] for r in rows]),
            interval_width=descr([i["width"] for i in iv]),interval_coverage=sum(i["covers"] for i in iv),
            interval_error=descr([abs(r["interval"]["estimate"]-r["truth"])
                                  for r in rows if r["interval"] and r["interval"]["estimate"] is not None]),
            interval_disjoint=sum(len(i["components"])>1 for i in iv),
            interval_empty=sum(len(i["components"])==0 for i in iv),
            interval_full_hull=sum(i["full_box_hull"] for i in iv),
            interval_cpu_s=descr([i["scalar_cpu_s"] for i in iv]),
            clean_nontruth_real_root_count=sum(len(i["clean_nontruth_real_roots"]) for i in iv),
            exam_min_true_visibility=descr([i["true_clean_visibility"] for i in iv]),
            off_anchor_rejects=descr([r["profile_off_anchor_rejections"] for r in rows]),
            off_anchor_null_rank_nonzero=descr([r["profile_off_anchor_null_dim_nonzero"] for r in rows]),
            witness_spread=descr(wit),maximum_feed_norm=max(r["maximum_feed_norm"] for r in rows),
            bounded_SNR_dB=descr([r["bounded_SNR_dB"] for r in rows]),
        ))
    out = dict(
        group_summary=combined,
        enrollment_cpu_s=raw["cpu_seconds"],common_total_energy=raw["common_total_energy"],
        max_feed_norm=max(r["maximum_feed_norm"] for r in raw["observations"]),
        max_reduced_full_field_norm_abs=max(r["full_reduced_diagnostic"] for r in raw["observations"]),
        median_born_response_disagreement=descr([r["relative_Born_response_disagreement"] for r in raw["observations"]]),
        numerical_checks=raw["numerical_checks"],simple_interval_checks=raw["simple_interval_checks"],
        chain_rank_diagnostics={
            g:{r:{k:v[k] for k in ("degree", "regional_observability_rank",
                                   "regional_observability_singular_values", "smallest_zero_degree_chain_singular",
                                   "coefficient_span_diagnostics")} for r,v in l["regions"].items()}
            for g,l in raw["design"].items()
        },
    )
    (RESULTS/"pencil_summary.json").write_text(json.dumps(out,indent=2)+"\n")
    return {(v["layout"],v["region"],v["noise"],v["arm"]): v for v in combined}


def geometry_panel(ax,d):
    a=np.asarray(d["active"])
    pos=np.asarray(d["cells"])
    sources=np.asarray(d["sources"])
    host=np.asarray(d["host"])
    T_L=a[d["regions"]["L"]["target"]]
    T_S=a[d["regions"]["square"]["target"]]
    ax.add_patch(Circle((0,0),1.76,ls=":",fc="none",ec="#b0b6c2",lw=.8))
    for j,(xy,v) in enumerate(zip(pos,host)):
        if j in T_L:
            color="#bf6b41"
        elif j in T_S:
            color="#d9ab55"
        elif j in a:
            color="white"
        else:
            color=plt.cm.Blues(.23+.1*v)
        ax.add_patch(Circle(xy,d["disk_radius"]*1.32,fc=color,ec="#537389",lw=.75))
    ax.scatter(*sources.T,s=21,facecolor="white",edgecolor=COLOR["pencil"],lw=1.15,label="10 feeds")
    ax.scatter(*np.asarray(d["receiver"]).ravel(),s=28,marker="D",c="#21252e",label="Receiver")
    ax.set(xlim=(-1.92,1.92),ylim=(-1.92,1.92),xlabel=r"$x/\lambda$",ylabel=r"$y/\lambda$")
    ax.set_aspect("equal")
    legend=[
        Line2D([0],[0],marker="o",color="white",markeredgecolor=COLOR["pencil"],ms=5,label="10 feeds"),
        Line2D([0],[0],marker="D",color="white",markerfacecolor="#21252e",ms=6,label="Receiver"),
        Line2D([0],[0],marker="o",color="white",markerfacecolor="#bf6b41",ms=7,label="L patch"),
        Line2D([0],[0],marker="o",color="white",markerfacecolor="#d9ab55",ms=7,label="Square adds cell"),
    ]
    ax.legend(handles=legend,loc="upper center",bbox_to_anchor=(.5,-.19),ncol=2,
              frameon=False,fontsize=7.8,columnspacing=.85)


def profiles_panel(ax,raw,region):
    d={r["arm"]:r for r in raw["observations"] if (
        r["region"]==region and r["layout"]=="confirmation" and r["epsilon"]==LAST and r["case"]==0
    )}
    q=np.asarray(raw["q_test_grid"])
    for arm in ("completion","pencil","anchor","random"):
        y=np.asarray(d[arm]["profile_grid"])[:,0]
        if arm=="anchor" and region=="square":
            keep=y>1e-10
            ax.plot(q[~keep],y[~keep],"--",color=COLOR[arm],linewidth=1.8)
            ax.scatter(q[keep],y[keep],marker="D",s=20,color=COLOR[arm],label=LABEL[arm],zorder=4)
        elif arm=="random":
            ax.scatter(q[::5],y[::5],marker="x",s=13,color=COLOR[arm],label=LABEL[arm],zorder=4)
        else:
            style="--" if arm=="anchor" else "-"
            ax.plot(q,y,style,color=COLOR[arm],lw=1.7 if arm=="anchor" else 2,label=LABEL[arm])
    ax.axhline(LAST,lw=1.0,ls=":",color="#6d6d73")
    ax.axvline(d["pencil"]["truth"],lw=1.0,ls=":",color="#444448")
    ax.text(1.1,LAST*1.24,"Noise radius",color="#66666d",fontsize=8)
    ax.set(xlim=(.8,15.25),ylim=(-.35e-5,.16),yscale="symlog",xlabel="Candidate patch contrast",
           ylabel="Complex nuisance profile distance")
    ax.set_yscale("symlog",linthresh=LAST/4)
    ax.legend(frameon=False,loc="upper center",fontsize=7.8,ncol=2,columnspacing=.8)


def error_panel(ax,summary):
    group=[(g,r) for g in LAYOUTS for r in REGIONS]
    order=("interval","pencil","anchor","random","completion")
    offset=np.linspace(-.32,.32,len(order))
    error_labels=dict(interval="Pencil interval",pencil="Pencil + TRF",anchor="Anchors + TRF",
                      random="Random + TRF",completion="Ten feeds + TRF")
    for i,arm in enumerate(order):
        vals=[summary[g,r,LAST,"pencil"]["interval_error"]["median"] if arm=="interval"
              else summary[g,r,LAST,arm]["trf_error"]["median"] for g,r in group]
        ax.bar(np.arange(len(group))+offset[i],vals,width=.145,color=COLOR[arm],label=error_labels[arm],zorder=3)
    ax.set_xticks(np.arange(len(group)),["Primary\nL","Primary\nsquare","Confirm.\nL","Confirm.\nsquare"])
    ax.set(yscale="log",ylim=(.003,10),ylabel="Median absolute contrast error")
    ax.legend(frameon=False,fontsize=7,loc="upper right",ncol=2,columnspacing=.8)


def rejection_panel(ax,summary):
    group=[(g,r) for g in LAYOUTS for r in REGIONS]
    for i,arm in enumerate(ARMS):
        vals=[summary[g,r,LAST,arm]["off_anchor_rejects"]["median"] for g,r in group]
        x=np.arange(len(group))+(i-1.5)*.18
        ax.bar(x,vals,width=.15,color=COLOR[arm],zorder=3)
        zeros=np.array(vals)==0
        ax.scatter(x[zeros],np.zeros(np.count_nonzero(zeros)),s=13,c=COLOR[arm],zorder=5,clip_on=False)
    ax.set_xticks(np.arange(len(group)),["Primary\nL","Primary\nsquare","Confirm.\nL","Confirm.\nsquare"])
    ax.set(ylim=(-.045,1.04),ylabel="Fraction of off-anchor candidates rejected")


def witness_panel(ax,raw):
    tickpos, tickname=[],[]
    short=dict(pencil="pencil",anchor="anchors",random="random",completion="ten feeds")
    for i,g in enumerate(LAYOUTS):
        v={o["arm"]:o for o in raw["observations"] if o["layout"]==g and o["region"]=="square"
            and o["case"]==0 and o["epsilon"]==LAST}
        off=(1-i)*5
        ax.axvline(v["pencil"]["truth"],ymin=(off-.2)/8.5,ymax=(off+3.25)/8.5,
                   color="#6d6d73",ls=":",lw=1)
        for k,arm in enumerate(ARMS):
            row=v[arm]
            y=off+3-k
            tickpos.append(y)
            tickname.append(("Primary" if i==0 else "Confirm.")+" "+short[arm])
            if arm=="pencil":
                for lo,hi in row["interval"]["components"]:
                    ax.hlines(y,lo,hi,color=COLOR["interval"],lw=5,alpha=.42)
            wit=row["witness"]["extrema"]
            if len(wit)==2:
                edges=[e["q"] for e in wit]
                ax.hlines(y,*edges,color=COLOR[arm],lw=1.4)
                ax.scatter(edges,[y]*2,c=COLOR[arm],s=17,zorder=5)
    ax.set_yticks(tickpos,tickname,fontsize=7.0)
    ax.set(xlim=(.7,15.3),ylim=(-.5,8.45),xlabel="Feasible contrast endpoints")
    ax.text(7.6,8.13,"Dots are full-VIE witnesses",ha="center",va="center",fontsize=7.6,color="#53535b")


def draw(raw,summary):
    plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
                         "axes.titlesize":11,"axes.labelsize":9,"xtick.labelsize":8.4,"ytick.labelsize":8.4,
                         "pdf.fonttype":42,"axes.labelpad":4})
    f,ax=plt.subplots(2,3,figsize=(13.4,7.9),gridspec_kw=dict(wspace=.35,hspace=.55))
    geometry_panel(ax[0,0],raw["design"]["primary"])
    profiles_panel(ax[0,1],raw,"L")
    profiles_panel(ax[0,2],raw,"square")
    rejection_panel(ax[1,0],summary)
    error_panel(ax[1,1],summary)
    witness_panel(ax[1,2],raw)
    titles=("Known host and variable patch","L patch, four shots","Square patch, three shots",
            "Off-anchor rejection","Patch contrast errors","Square patch, feasible witnesses")
    for i,(p,title) in enumerate(zip(ax.flat,titles)):
        p.set_title(f"({chr(i+97)}) {title}",loc="left",pad=9,fontweight="bold")
        if i:
            p.grid(axis="y" if i!=5 else "x",color="#d8dbe1",alpha=.63,lw=.7)
            p.set_axisbelow(True)
    f.subplots_adjust(left=.08,right=.985,top=.929,bottom=.10)
    f.savefig(FIGURES/"pencil_comparison.pdf",facecolor="white")
    f.savefig(FIGURES/"pencil_comparison.png",dpi=240,facecolor="white")
    plt.close(f)




def disconnected_figure(raw, followup):
    probe=followup["disjoint_probes"][0]
    row=next(r for r in raw["observations"] if r["arm"]=="pencil"
             and (r["layout"],r["region"],r["case"])==(probe["layout"],probe["region"],probe["case"])
             and len(r["interval"]["components"])>1)
    data=row["interval"]["g"]
    coords=raw["design"][probe["layout"]]["regions"][probe["region"]]["template"]["coord"]
    g=np.asarray(data["real"])+1j*np.asarray(data["imag"])
    c=np.asarray(coords["real"])+1j*np.asarray(coords["imag"])
    c/=np.sqrt(raw["common_total_energy"]/row["shots"])
    score_q=np.linspace(1.,15.,2001)
    t=(score_q-8.)/7.
    min_res=abs(poly.polyval(t,g))/np.linalg.norm(poly.polyval(t,c.T),axis=0)
    eps=row["noise_radius_with_rounding"]
    f,ax=plt.subplots(figsize=(7.1,3.5))
    ax.semilogy(score_q,min_res/eps,color=COLOR["pencil"],lw=1.8,label="Focused residual")
    for j,(lo,hi) in enumerate(row["interval"]["components"]):
        ax.axvspan(lo,hi,color=COLOR["interval"],alpha=.2,label="Outer interval components" if j==0 else None)
    ax.axhline(1,color="#62636a",ls=":",lw=1)
    ax.text(14.95,1.06,"Noise radius",ha="right",va="bottom",fontsize=8,color="#55555d")
    ax.axvline(row["truth"],color="#252832",ls=":",lw=1,label="True contrast")
    q=[x["q"] for x in probe["attempts"] if x["independently_feasible"]]
    v=[x["full_residual"]/eps for x in probe["attempts"] if x["independently_feasible"]]
    ax.scatter(q,v,marker="o",facecolors="white",edgecolors=COLOR["anchor"],s=22,
               zorder=5,label="Full-VIE feasible material")
    ax.set(xlabel="Candidate L-patch contrast",ylabel="Residual / noise radius",
           xlim=(.8,15.2),ylim=(.04,105),title="One noisy observation admits two separate contrast intervals")
    ax.grid(axis="y",alpha=.22)
    ax.set_axisbelow(True)
    ax.legend(frameon=False,ncol=1,fontsize=8,loc="lower right")
    f.tight_layout()
    f.savefig(FIGURES/"pencil_disconnected.pdf",facecolor="white")
    f.savefig(FIGURES/"pencil_disconnected.png",dpi=240,facecolor="white")
    plt.close(f)


if __name__=="__main__":
    raw=json.loads((RESULTS/"pencil_raw.json").read_text())
    summary=summarize(raw)
    draw(raw,summary)
    followup=RESULTS/"pencil_feasible_followup.json"
    if followup.is_file():
        disconnected_figure(raw,json.loads(followup.read_text()))
