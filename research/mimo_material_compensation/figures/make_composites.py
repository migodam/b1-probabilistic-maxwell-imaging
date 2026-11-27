"""Render two print figures directly from the frozen field and certificate arrays."""
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as pe
from matplotlib.ticker import MaxNLocator
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ENV, BUDGET, PATHS = (BASE / p for p in ("environment", "budgets", "paths"))
INK = "#152733"
MUTED = "#4b5962"
GRID = "#dfe5e6"
TEAL = "#146b72"
TEAL_PALE = "#dcefed"
RUST = "#ae4635"
RUST_PALE = "#f6e3dd"
GRAY = "#efeee9"
PURPLE = "#535184"
CORE = "#f28145"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 9,
    "axes.titlesize": 10.0, "axes.titleweight": "bold", "axes.titlepad": 7,
    "axes.edgecolor": "#7b8990", "axes.linewidth": .7, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.major.size": 3, "ytick.major.size": 3, "legend.fontsize": 8,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "mathtext.fontset": "dejavusans", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def read(folder, name):
    return json.loads((folder / name).read_text())


def outward_label(value, is_lower, places=3):
    quantum = Decimal(1).scaleb(-places)
    rounding = ROUND_FLOOR if is_lower else ROUND_CEILING
    return format(Decimal.from_float(float(value)).quantize(quantum,rounding=rounding),"f")


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.pdf", metadata={"Creator": "Matplotlib certificate plotting"})
    fig.savefig(HERE / f"{stem}.png", dpi=300)
    plt.close(fig)


def heading(fig, title, subtitle):
    fig.text(.065, .958, title, fontsize=17, fontweight="bold", va="top")
    fig.text(.065, .912, subtitle, fontsize=9, color=MUTED, va="top")


def chart_axis(ax, x=False):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="both" if x else "y", color=GRID, linewidth=.55)


def material(ax, array, mask, side, title, dosage):
    cmap = ax.imshow(array.reshape(side, side), origin="lower", extent=(-1, 1, -1, 1),
                     interpolation="nearest", cmap="cividis", norm=Normalize(0, 4.5))
    edge = np.linspace(-1., 1., side+1)
    roi = np.pad(mask.reshape(side, side), 1)
    outline = []
    for j,i in np.argwhere(mask.reshape(side, side)):
        x0,x1,y0,y1 = edge[i],edge[i+1],edge[j],edge[j+1]
        for dj,di,segment in [(-1,0,[(x0,y0),(x1,y0)]),(1,0,[(x0,y1),(x1,y1)]),
                              (0,-1,[(x0,y0),(x0,y1)]),(0,1,[(x1,y0),(x1,y1)])]:
            if not roi[j+1+dj,i+1+di]:
                outline.append(segment)
    ax.add_collection(LineCollection(outline,colors=CORE,linewidths=1.12,
        path_effects=[pe.Stroke(linewidth=2.55,foreground="#25212b"),pe.Normal()]))
    ax.text(-.92, .86, dosage, color="white", fontsize=8.4, fontweight="bold",
            bbox={"facecolor":"#1c2831", "edgecolor":"none", "pad":3.1, "alpha":.91})
    ax.set_title(title, loc="left", fontsize=9.2)
    ax.set_xticks((-1, 0, 1))
    ax.set_yticks((-1, 0, 1))
    ax.set_aspect("equal")
    ax.tick_params(length=2.5, pad=2, labelsize=7.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return cmap


def response_bounds(ax, values, upp, midpoint, title, area):
    beta = np.array([row["beta"] for row in values]) * 1000
    lower = np.array([row["open_lower_over_two_r"] for row in values])
    max_x = 27.5
    ax.axvspan(beta[-1], max_x, color=GRAY, alpha=.79, zorder=0)
    # A class bound at a tested budget also holds for every nested, smaller box.
    x = np.r_[0., beta[beta > 0]]
    effective = lower[beta > 0]
    for left, right, value in zip(x[:-1], x[1:], effective):
        ax.fill_between([left, right], [1.,1.], [value,value],
                        color=RUST_PALE, alpha=.82, linewidth=0, zorder=1)
        ax.hlines(value, left, right, color=RUST, linewidth=1.65, zorder=3)
    ax.plot(beta, lower, "o", color=RUST, ms=3.65, zorder=4)
    ax.hlines(upp, 0, max_x, color=TEAL, linewidth=1.55)
    ax.hlines(midpoint, 0, max_x, color=TEAL, linewidth=1.05, linestyles="dotted")
    ax.axhline(1., color=INK, linewidth=.85, linestyle=(0, (4, 3)))
    ax.annotate(outward_label(lower[-1],True,4), (beta[-1], lower[-1]), xytext=(0, 9),
                textcoords="offset points", ha="right", color=RUST, fontsize=8.3)
    ax.text(max_x-.45, upp-.07, outward_label(upp,False,4), color=TEAL, size=7.7, va="top", ha="right")
    ax.text(max_x-.45, midpoint-.07, outward_label(midpoint,False,4), color=TEAL, size=7.7, va="top", ha="right")
    ax.set(xlim=(0,max_x), ylim=(0,4.02), xticks=[0,5,10,15,20,25],
           ylabel=r"Separation bound $d/(2r)$", title=title)
    ax.text(.98,.94, area, transform=ax.transAxes, va="top", ha="right",
            color=MUTED, fontsize=8)
    chart_axis(ax)


def matched():
    c32 = read(ENV,"frequency_classes_open.json")["records"][0]["bands"]
    c96 = read(ENV,"matched_amplitude.json")
    energy = read(ENV/"low_energy_32","proposal.json")
    configs = [
        (32,"frequency_classes_32.npz",c32["460"]["budgets"],
         c32["400"]["open_400_common_midpoint_residual_over_r_upper"],
         read(ENV,"kyp_verified_responses_32_2.0_phase_32_two_shells_400.json")["upper_over_2r"]),
        (96,"matched_amplitude.npz",c96["outputs"]["460"]["budgets"],
         c96["outputs"]["400"]["open_common_midpoint_residual_over_r_upper"],
         read(ENV,"kyp_verified_responses_96_2.0_wave_96_2.0_1_400.json")["upper_over_2r"]),
    ]
    fig = plt.figure(figsize=(11.5,8.25))
    heading(fig,"The same core increase on two native grids",
            "For each grid, the public reference specifications share possible 400 MHz data. A 460 MHz action separates the two material classes.")
    gs = fig.add_gridspec(2,4,left=.069,right=.972,bottom=.205,top=.855,
                          height_ratios=(1,1),width_ratios=(1.06,1.06,.09,2.02),
                          hspace=.33,wspace=.31)
    for row,(side,file,values,mid,upp) in enumerate(configs):
        with np.load(ENV/file,allow_pickle=False) as raw:
            low,high,core = (raw[key] for key in ("minus","plus","core"))
        if not np.all(low[core]==0) or not np.all(high[core]==2):
            raise ValueError("The matched-amplitude core is not exact")
        lab = "abc" if row==0 else "def"
        for col,(image,title,dose) in enumerate([(low,"Low reference",r"$\chi_R=0$"),
                                                 (high,"High reference",r"$\chi_R=2$")]):
            a=fig.add_subplot(gs[row,col])
            mat=material(a,image,core,side,f"{lab[col]}   {title}",dose)
            if col:
                a.set_yticklabels([])
            else:
                a.set_ylabel("Position [m]")
            if row:
                a.set_xlabel("Position [m]")
        bound=fig.add_subplot(gs[row,3])
        cells=int(core.sum())
        response_bounds(bound,values,upp,mid,f"{lab[2]}   Native {side} × {side}",
                        f"{cells} core cells · {cells*(2/side)**2:.5f} m²")
        bound.set_xlabel(r"All-pixel box radius $\beta$ ($\times 10^{-3}$)")
    cb = fig.colorbar(mat,cax=fig.add_subplot(gs[:,2]),ticks=[0,1,2,3,4.5])
    cb.ax.set_title(r"$\chi$", pad=6, fontsize=10)
    cb.ax.yaxis.set_ticks_position("left")
    cb.ax.tick_params(labelsize=7.4,length=2,pad=3)
    legend = [
        Line2D([],[],color=RUST,marker="o",markersize=4,label="460 MHz lower bound, every pair (added port open)"),
        Line2D([],[],color=TEAL,label="400 MHz reference-pair upper bound, all passive loads and coherent drives"),
        Line2D([],[],color=TEAL,linestyle="dotted",label="400 MHz open-port midpoint, reference pair"),
        Patch(facecolor=GRAY,edgecolor="none",label="Beyond the verified 460 MHz range"),
    ]
    fig.legend(handles=legend,loc="upper left",bbox_to_anchor=(.066,.143),ncol=2,
               frameon=False,borderaxespad=0,columnspacing=1.6,handlelength=2.25)
    fig.text(.069,.086,"Full action at either frequency: 8 shots, 64 complex reads, 1 W summed available power. "
                      "The two grids use separate native reference specimens.",fontsize=8.15,color=MUTED)
    fig.text(.069,.059,f"At β = 0.01 on 32 × 32 and fixed absolute noise, a 460 MHz action resolved the calibrated classes "
                      f"at {energy['operating_energy']:.3f} E₀. A 400 MHz separator needs >{outward_label(energy['necessary_400'],True)} E₀.",
             fontsize=8.15,color=MUTED)
    fig.text(.069,.033,"E₀ denotes the full-action available energy at equal dwell. The hollow core outline marks the same physical task in each material pair.",
             fontsize=7.8,color=MUTED)
    save(fig,"matched_core_amplitude")


def key_row(items, archive, name):
    return [v for v in items if v["source"]==f"{archive} / {name}"]


def paths_panel(ax, data, side, labels, xlim, ylim, title):
    curves=[]
    for archive,name,text,color in labels:
        raw=key_row(data,archive,name)
        if not raw:
            raise ValueError(f"No certified stripes in {archive} / {name}")
        steps = np.array([[part["task_t_from"],part["task_t_to"],part["sufficient_exterior_budget"]["L2_m"]]
                          for part in raw])
        curves.append((steps,text,color))
    # Shade each certified step. A faint second colour can overlap where both tasks are feasible.
    for steps,text,color in curves:
        for x0,x1,y in steps:
            ax.fill_between([x0,x1],[y,y],[ylim,ylim],color=color,alpha=.09,linewidth=0)
        line_x=np.r_[steps[0,0],steps[:,1]]
        ax.stairs(steps[:,2],line_x,color=color,lw=1.7,baseline=None,label=text)
        ax.plot(steps[-1,1],steps[-1,2],"o",ms=4,color=color,zorder=4)
    if side==32:
        ax.axvline(2.,color=INK,lw=.7,ls=(0,(3,3)))
        ax.text(1.993,.10,"2",ha="right",va="bottom",color=MUTED,fontsize=7.5)
    ax.set(xlim=(0,xlim),ylim=(0,ylim),title=title,xlabel="Mean increase in the task region",
           ylabel=r"Exterior budget $\|e\|_{2,\mathrm{area}}$ [m]")
    ax.legend(loc="lower right",frameon=True,facecolor="white",edgecolor="none",
              framealpha=.88,handlelength=1.65,fontsize=7.6)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    chart_axis(ax,True)


def bridge_panel(ax, data):
    rows=data["records"]
    xs=np.array([u["beta"] for u in rows])*1000
    excluded=np.array([u["certified_alpha_through"] for u in rows])
    hidden=np.array([u["witness"]["alpha_upper"] for u in rows])
    top=.263
    ax.set_facecolor(GRAY)
    for i in range(len(xs)-1):
        # Inclusion of the hypothesis sets makes both step implications rigorous.
        w=xs[i+1]-xs[i]
        ax.add_patch(Rectangle((xs[i],0.),w,excluded[i+1],fc=RUST_PALE,ec="none"))
        ax.add_patch(Rectangle((xs[i],hidden[i]),w,top-hidden[i],fc=TEAL_PALE,ec="none"))
        ax.hlines(excluded[i+1],xs[i],xs[i+1],lw=1.35,color=RUST)
        ax.hlines(hidden[i],xs[i],xs[i+1],lw=1.35,color=TEAL)
    ax.scatter(xs,excluded,s=14,color=RUST,zorder=5,clip_on=False)
    ax.scatter(xs,hidden,s=14,color=TEAL,zorder=5,clip_on=False)
    ax.axvline(2.5,color="#566570",lw=.7,ls=(0,(3,3)),zorder=1)
    ax.text(1.34,.104,"Unresolved",color=MUTED,size=8.2,ha="center",va="center")
    ax.text(.10,.016,"Cannot hide",color=RUST,size=7.8,va="bottom")
    ax.text(2.94,.235,"Can hide",color=TEAL,size=8.1,ha="right")
    ax.set(xlim=(0.,3.),ylim=(0,top),xticks=xs,title="c   Interleaved B with variable ROI shape",
           xlabel=r"Exterior maximum budget $\beta_\infty$ ($\times 10^{-3}$)",
           ylabel=r"Allowed per-cell deviation $\alpha$ from mean")
    ax.spines[["top","right"]].set_visible(False)


def reversal_panel(ax, rows, bridge):
    flat_a=next(u for u in rows if u["t"]==.27 and u["beta"]==.0025 and u["layout"]=="A" and u["contract"]=="prospective" and u["kind"]=="max")
    flat_b=next(u for u in rows if u["t"]==.27 and u["beta"]==.0025 and u["layout"]=="B" and u["contract"]=="prospective" and u["kind"]=="max")
    free=next(u for u in bridge["records"] if u["beta"]==0.)["witness"]
    entries=[(2,flat_a["witness_upper_over_threshold"],TEAL,False),
             (1,flat_b["lower_over_threshold"],RUST,True),
             (0,free["strict_W_over_two_r"],TEAL,False)]
    lo,hi=.79,1.13
    ax.axvspan(lo,1.,color=TEAL_PALE,alpha=.65)
    ax.axvspan(1.,hi,color=RUST_PALE,alpha=.55)
    ax.axvline(1.,color=INK,lw=.75,ls=(0,(3,3)))
    for y,value,color,lower in entries:
        ax.scatter([value],[y],c=color,s=31,zorder=3)
        edge=value+.048 if lower else value-.063
        ax.annotate("",xy=(edge,y),xytext=(value,y),
            arrowprops={"arrowstyle":"-|>","lw":1.35,"color":color})
        ax.text(value,y+.19,f"{'≥' if lower else '≤'} {outward_label(value,lower)}",
                color=color,ha="center",va="bottom",fontsize=8.8,fontweight="bold")
    ax.set(xlim=(lo,hi),ylim=(-.46,2.56),yticks=[0,1,2],
           yticklabels=["B  Mean only\nβ∞ = 0", "B  Uniform\nβ∞ = 0.0025", "A  Uniform\nβ∞ = 0.0025"],
           title="d   Same mean increase, different conclusion",
           xlabel=r"Certified separation relative to $2r$")
    ax.tick_params(axis="y",length=0,pad=6,labelsize=8)
    ax.spines[["top","right","left"]].set_visible(False)


def compensation():
    phase=read(PATHS,"phase_bands.json")["proved_hideable_segments"]
    bridge=read(BUDGET/"bridge","receipt.json")
    focus=read(BUDGET,"budget_focus_joint32.json")["records"]
    fig=plt.figure(figsize=(11.5,8.5))
    heading(fig,"What can exterior compensation hide?",
            "Top: certified material paths under the original noisy A measurements. Bottom: a fixed-reference prospective test.")
    outer=fig.add_gridspec(2,2,left=.082,right=.974,top=.848,bottom=.172,hspace=.44,wspace=.34)
    paths_panel(fig.add_subplot(outer[0,0]),phase,32,[
        ("n32_c01_paths_broad.json","ball_min_margin","Uniform 43-cell patch",TEAL),
        ("n32_c01_paths_whole_region.json","ball_shape_free","Free shape in 172-cell ROI",PURPLE)],
        2.03,3.61,"a   Original A · native 32 × 32")
    paths_panel(fig.add_subplot(outer[0,1]),phase,96,[
        ("n96_c01_paths_broad.json","ball_min_margin","Uniform 156-cell patch",TEAL),
        ("n96_c01_paths_ball.json","ball_shape_free","Free shape in 156-cell patch",PURPLE)],
        .49,1.04,"b   Original A · native 96 × 96")
    bridge_panel(fig.add_subplot(outer[1,0]),bridge)
    right=outer[1,1].subgridspec(1,2,width_ratios=[.12,.88],wspace=0)
    reversal_panel(fig.add_subplot(right[0,1]),focus,bridge)
    legend=[Patch(facecolor=TEAL_PALE,edgecolor="none",label="Proven hideable"),
            Patch(facecolor=RUST_PALE,edgecolor="none",label="Proven unhideable"),
            Patch(facecolor=GRAY,edgecolor="none",label="Unresolved")]
    fig.legend(handles=legend,loc="upper left",bbox_to_anchor=(.082,.112),ncol=3,
               frameon=False,borderaxespad=0,handlelength=1.25,columnspacing=1.45)
    fig.text(.082,.063,"Panels a–b shade sufficient budgets along each verified material path. Panels c–d use the same six-cell 32² c02 region, "
                      "mean increase 0.27,",color=MUTED,fontsize=8.15)
    fig.text(.082,.037,"and a fixed public 400 MHz reference. In d, each arrow points to the side allowed by the certified bound.",
             color=MUTED,fontsize=8.15)
    save(fig,"compensation_phase")


if __name__ == "__main__":
    matched()
    compensation()
