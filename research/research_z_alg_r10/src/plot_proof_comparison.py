"""Compare archived whole-region proofs and the joint-error counterexample."""

import math

from plot_certified import ALTERNATIVE, GRID, MUTED, NEW, OLD, PAPER, ROOT, TRUTH, clean_axes, objective, plt, read, setup
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np


def proof_runs(comparison):
    runs = []
    for row in comparison["regions"]:
        pivot = row["r10_pivot_bernstein"]
        direct = row["direct_full_bernstein"]
        fallback = row.get("r10_chart_fallback", {})
        if not direct["excluded_all"] or not (pivot["excluded_all"] or fallback.get("excluded_all")):
            continue
        shared = row["shared_witness_cpu_seconds"]
        pivot_time = shared + fallback.get("total_cpu_seconds_including_first_chart", pivot["cpu_seconds"])
        direct_time = shared + direct["cpu_seconds"]
        runs.append((bool(fallback), row["sample"], row["old_witness_id"], pivot_time*1000, direct_time*1000))
    return sorted(runs)


def proof_panel(ax, runs):
    for i, (_, _, _, reduced, direct) in enumerate(runs):
        ax.plot([reduced,direct], [i,i], color=OLD, linewidth=2.2, zorder=2)
        ax.scatter(reduced,i,s=76,color=NEW,edgecolor=PAPER,linewidth=.9,zorder=4)
        ax.scatter(direct,i,s=64,marker="s",color=TRUTH,edgecolor=PAPER,linewidth=.9,zorder=4)
        ax.annotate(f"{reduced:.2f}", (reduced,i), textcoords="offset points", xytext=(0,-13),
                    ha="center", va="top", fontsize=8.7, color=NEW)
        ax.annotate(f"{direct:.2f}", (direct,i), textcoords="offset points", xytext=(0,-13),
                    ha="center", va="top", fontsize=8.7, color=TRUTH)
        ax.annotate(f"{direct/reduced:.1f} 倍", (math.sqrt(reduced*direct),i),
                    textcoords="offset points", xytext=(0,7), ha="center", va="bottom", fontsize=10.5)
    labels = [f"{name} · 见证 {index}" for _,name,index,_,_ in runs]
    ax.set(xscale="log", xlim=(1.45, 310), ylim=(len(runs)-1+.64,-.46),
           yticks=list(range(len(runs))), yticklabels=labels,
           xlabel="证明同一区域无解的单核 CPU 时间（毫秒，对数轴）")
    ax.xaxis.set_major_locator(FixedLocator([2,5,10,20,50,100,200]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
    ax.minorticks_off()
    ax.grid(axis="x",color=GRID,linewidth=.9)
    clean_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y",left=False,pad=8)


def counterexample_panel(ax, example):
    cert = example["certificate"]
    assert cert["strict_counterexample"]
    singleton_count, pair_count = cert["strict_singletons"], cert["strict_pairs"]
    singleton_max = objective({"objective_enclosure":cert["bounds"]["max_single_score_upper"]})
    pair_max = objective({"objective_enclosure":cert["bounds"]["max_pair_score_upper"]})
    forward = objective(cert["full_forward_witness"])
    with np.load(ROOT/"results"/"natural_inputs.npz") as data:
        read_count = len(data[f"{example['case']}/active/y"])
    labels = [
        f"{singleton_count} 条关系逐条检验\n最高分",
        f"{pair_count} 对关系两两检验\n最高分",
        f"全 {read_count} 束共用同一组误差\n完整正演分数",
    ]
    colors = [NEW, NEW, ALTERNATIVE]
    levels = [singleton_max,pair_max,forward]
    ax.axvspan(0,1,color="#ebf5f1",zorder=0)
    ax.axvline(1,color="#657988",linewidth=1.25,linestyle="--",zorder=3)
    ax.barh(np.arange(3),levels,color=colors,height=.29,zorder=2)
    for i, value in enumerate(levels):
        ax.text(value-.025,i,f"{value:.4f}",color=PAPER,ha="right",va="center",fontsize=10.5)
    ax.annotate("误差上限 1", (1,-.35), textcoords="offset points", xytext=(0,3),
                ha="center",va="bottom",color=MUTED,fontsize=9.5)
    ax.set(xlim=(0,1.36), ylim=(2.55,-.64), yticks=[0,1,2], yticklabels=labels,
           xticks=[0,.25,.5,.75,1,1.25], xlabel="解释读数所需的最小归一化联合误差平方")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
    ax.grid(axis="x",color=GRID,linewidth=.9)
    clean_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y",left=False,pad=8)


def main():
    setup()
    comparison = read("direct_box_unselected.json")
    example = read("all_circuit_counterexample.json")["example"]
    runs = proof_runs(comparison)
    fig = plt.figure(figsize=(13.8,5.25))
    grid = fig.add_gridspec(1,2,left=.148,right=.97,top=.752,bottom=.249,wspace=.53,width_ratios=(1.08,1))
    proof_ax, counter_ax = (fig.add_subplot(grid[0,i]) for i in range(2))
    proof_panel(proof_ax,runs)
    counterexample_panel(counter_ax,example)
    fig.text(.048,.948,"全盒排除的证明成本与局部误差检验的局限",fontsize=22,va="top")
    fig.text(.049,.853,f"a   同一误差界下严格排除 {len(runs)} 个区域",fontsize=13)
    fig.text(.593,.853,f"b   {example['case']} · 同一份读数与同一试探材料",fontsize=13)
    handles = [
        Line2D([],[],color=NEW,marker="o",ls="None",label="降维外包证明"),
        Line2D([],[],color=TRUTH,marker="s",ls="None",label="直接六维证明"),
    ]
    fig.legend(handles=handles,ncol=2,frameon=False,loc="upper left",bbox_to_anchor=(.146,.154),
               fontsize=9,handletextpad=.35,columnspacing=1.4,borderaxespad=0)
    fig.text(.049,.078,"两种时间均计入同一次见证核验，降维法还计入更换坐标对的时间。",fontsize=9,color=MUTED)
    fig.text(.593,.140,"每对关系内部共用一次误差解释。所有单项与双项均在限内。",
             fontsize=9,color=MUTED)
    fig.text(.593,.078,"完整正演要求所有读数共用同一组误差，严格排除该试探材料。",
             fontsize=9,color=MUTED)
    for ext in ("png","pdf"):
        fig.savefig(ROOT/"figures"/f"proof_comparison.{ext}",dpi=230)
    plt.close(fig)


if __name__ == "__main__":
    main()
