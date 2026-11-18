"""Plot the conditional noise-risk prediction and the recorded Maxwell checks."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.font_manager import get_font_names
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np

from risk_design import conditional_risk


ROOT = Path(__file__).resolve().parents[1]
INK = "#162938"
MUTED = "#475761"
GRID = "#e2e9ec"
BLUE = "#12769a"
ORANGE = "#aa5e13"

METHODS = ("noise_free_plus_one_percent", "ten_percent_risk", "double_band")
GEOMETRIES = {"mirror": ("镜像几何", BLUE, "o"),
              "next_balanced_layout": ("新几何（第七组）", ORANGE, "s")}


def prediction(r):
    one = np.ones(1)
    return conditional_risk(one, one, np.ones((1, 1)), mu=0., var=r*r/2,
                            sigma=.45, P=2., side=1)


def wilson(count, total, z=1.95996398454):
    value = count/total
    correction = 1 + z*z/total
    center = (value + z*z/(2*total))/correction
    width = z*np.sqrt((value*(1-value) + z*z/(4*total))/total)/correction
    return max(0.,center-width), min(1.,center+width)


def appearance(ax, title, xlabel, ylabel=None, grid="both"):
    ax.set_title(title, loc="left", pad=11, fontsize=14, fontweight="bold", color=INK)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.grid(axis=grid, color=GRID, lw=.8)
    ax.set_axisbelow(True)
    ax.spines[["top","right"]].set_visible(False)
    ax.spines[["left","bottom"]].set_color("#a4b5bf")
    ax.tick_params(colors=MUTED, length=3)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def plotted_counts(row, method):
    samples = [s for s in row["samples"] if s["method"] == method]
    return (len(samples), sum(s["weak_minima"]>=2 for s in samples),
            sum(s["nearer_generating_branch"] for s in samples))


def make_figure(records):
    fig = plt.figure(figsize=(11.8,5.18),facecolor="white")
    gs = fig.add_gridspec(1,2,left=.08,right=.965,bottom=.23,top=.71,
                          width_ratios=(1.26,1.),wspace=.55)
    fig.text(.053,.954,"从预测的存活概率选择带宽",
             ha="left",va="top",fontsize=22,fontweight="bold",color=INK)
    fig.text(.055,.844,"基频两侧各加测一次。增益已知，归一高斯噪声标准差为 0.45，存活概率在采样前计算。",
             ha="left",va="top",fontsize=10.7,color=MUTED)

    ax = fig.add_subplot(gs[0,0])
    r = np.linspace(.72,2.11,200)
    ax.plot(r,[prediction(x) for x in r],color=INK,lw=2.2,zorder=3)
    xs = np.array([1.01,records[0]["chosen_scale"],2.])
    for j,x in enumerate(xs):
        ax.axvline(x,color="#b3c1c8",lw=.85,ls=":",zorder=1)
        if j<2:
            px = max(records[0]["risk"][METHODS[j]])
            ax.scatter(x,px,marker="D",s=32,color=INK,zorder=5)
            offset=(-7,31) if j==0 else (-11,22)
            ax.annotate(f"预测 {100*px:.1f}%",(x,px),xytext=offset,textcoords="offset points",
                        fontsize=9.0,color=INK,ha="right",
                        arrowprops={"arrowstyle":"-","color":MUTED,"lw":.8},
                        bbox={"facecolor":"white","edgecolor":"none","pad":.2})

    for i,row in enumerate(records):
        _,color,marker = GEOMETRIES[row["label"]]
        for j,(method,x) in enumerate(zip(METHODS,xs)):
            total, count, _ = plotted_counts(row,method)
            risk = count/total
            low,high = wilson(count,total)
            # Horizontal jitter separates the two estimates at the same bandwidth.
            plot_x = x + (-.021 if i==0 else .021)
            ax.errorbar(plot_x,risk,yerr=[[max(0,risk-low)],[max(0,high-risk)]],
                        marker=marker,ls="none",color=color,markersize=5.,
                        markeredgecolor="white",markeredgewidth=.55,
                        elinewidth=1.25,capsize=3.,zorder=4)
            if j==1:
                ax.annotate(f"{count}/{total}",(plot_x,risk),
                            xytext=(8,8 if i==0 else -9),textcoords="offset points",
                            color=color,fontsize=9.4,ha="left",va="center",zorder=6)

    ax.set(xlim=(.73,2.095),ylim=(-.033,.90),
           xticks=xs,xticklabels=("1.01",f"{xs[1]:.3f}","2.00"),
           yticks=np.arange(0,.81,.2))
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1,decimals=0))
    appearance(ax,"a  第二个局部极小仍存在","归一半带宽 r（无噪声局部阈值为 1）","存活概率")
    handles=[Line2D([],[],color=INK,lw=2.2)]
    names=["采样前理论预测"]
    for name,color,marker in GEOMETRIES.values():
        handles.append(Line2D([],[],color=color,marker=marker,ls=""))
        names.append(name)
    ax.legend(handles,names,ncol=3,loc="upper left",bbox_to_anchor=(-.04,-.185),
              frameon=False,fontsize=8.9,handlelength=1.15,handletextpad=.35,columnspacing=.85)

    ax = fig.add_subplot(gs[0,1])
    y = np.arange(3)
    for i,row in enumerate(records):
        _,color,_ = GEOMETRIES[row["label"]]
        bar_y = y+(-.16 if i==0 else .16)
        for j,method in enumerate(METHODS):
            total,_,count=plotted_counts(row,method)
            ax.barh(bar_y[j],total,height=.282,color="#e7eef0",zorder=2)
            ax.barh(bar_y[j],count,height=.282,color=color,zorder=3)
            ax.text(count-.72,bar_y[j],f"{count}/{total}",
                    ha="right",va="center",fontsize=9.8,fontweight="bold",color="white")
    ax.set(ylim=(-.48,2.48),xlim=(0,40.12),yticks=y,
           yticklabels=("r = 1.01",f"r = {xs[1]:.3f}","r = 2.00"),
           xticks=np.arange(0,41,10))
    ax.invert_yaxis()
    appearance(ax,"b  从旧错误支起步的 TRF","恢复到更接近生成支的次数",grid="x")

    fig.text(.055,.042,"每个几何、每档带宽各 40 次，两支各生成 20 次。a 中竖线为 95% Wilson 二项近似区间，三档带宽使用配对的噪声抽样。",
             color=MUTED,fontsize=9.15)
    target=ROOT/"figures/noise_risk"
    fig.savefig(target.with_suffix(".png"),dpi=230,facecolor="white")
    fig.savefig(target.with_suffix(".pdf"),facecolor="white")
    plt.close(fig)


def main():
    records = json.loads((ROOT/"results/risk_design.json").read_text())
    fonts = set(get_font_names())
    family = next((x for x in ("Noto Sans CJK SC","Noto Sans SC","PingFang SC","STHeiti")
                   if x in fonts),None)
    if not family:
        raise RuntimeError("Install a Simplified Chinese font such as Noto Sans CJK SC.")
    with plt.rc_context({"font.family":["DejaVu Sans",family],"font.size":10.5,
                         "axes.unicode_minus":False,"pdf.fonttype":3}):
        make_figure(records)


if __name__ == "__main__":
    main()
