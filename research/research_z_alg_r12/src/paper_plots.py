"""Make the Chinese spectral-moment and box-face figures from saved results."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.font_manager import get_font_names
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
INK = "#162938"
MUTED = "#475761"
GRID = "#e2e9ec"
BLUE = "#12769a"
GOLD = "#aa5e13"
PURPLE = "#83558d"
GRAY = "#6a7880"

SCHEDULES = ("opposite_edges", "two_positive", "two_negative")
SCHEDULE_LABELS = ("正负各一次", "两次正偏", "两次负偏")
SCHEDULE_COLORS = (BLUE, GOLD, PURPLE)
FACIAL_CONTROLS = ("profile", "fixed_face_ignored", "largest_residual", "random")
FACIAL_LABELS = ("剖面选源", "忽略联动选源", "最大残差", "固定随机源")


def read(filename):
    return json.loads((RESULTS / filename).read_text())


def style(ax, title, xlabel=None, ylabel=None, grid="x"):
    ax.set_title(title, loc="left", pad=12, fontweight="bold", color=INK, fontsize=13)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=GRID, linewidth=.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#a8bac2")
    ax.tick_params(colors=MUTED, length=3, width=.7)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def header(fig, title, subtitle):
    fig.text(.06, .965, title, ha="left", va="top", fontsize=21, fontweight="bold", color=INK)
    fig.text(.061, .897, subtitle, ha="left", va="top", fontsize=10.2, color=MUTED)


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=230, facecolor="white")
    fig.savefig(FIGURES / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def worst_band(row, pattern, eta=1.):
    runs = [r for r in row["runs"] if r["pattern"] == pattern and r["eta"] == eta]
    if len(runs) != 2 or not all(r["collision_valid"] for r in runs):
        raise ValueError(f"Incomplete paired collision for {row['label']}, {pattern}, {eta}")
    return max(r["numerical_band"] for r in runs)


def spectral_figure():
    layouts = read("spectral_new_layouts.json")
    layouts = [x for x in layouts if x["included"]]
    fig = plt.figure(figsize=(10.7, 7.45), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=.12, right=.96, bottom=.17, top=.81,
                          width_ratios=(1.22, 1), height_ratios=(1.08, 1), wspace=.39, hspace=.68)
    header(fig, "双侧频偏较早消除较难一支的错误极小",
           "六组独立物理布局均含局部折叠。a、b 各方案新增两次半能量曝光，几何位置固定。")

    ax = fig.add_subplot(gs[0,0])
    bands = np.array([[worst_band(row,p) for p in SCHEDULES] for row in layouts])
    normalized = bands / bands[:,:1]
    for i, row in enumerate(normalized):
        ax.hlines(i, 1., max(row), color="#cbd7dc", lw=2)
        ax.scatter(row[0], i, marker="D", s=42, color=BLUE, zorder=5)
        ax.scatter(row[1], i-.13, color=GOLD, marker="o", s=35, zorder=5)
        ax.scatter(row[2], i+.13, color=PURPLE, marker="o", s=35, zorder=5)
    ax.set(ylim=(-.57, len(layouts)-.45), xlim=(.93, 1.93),
           yticks=np.arange(len(layouts)), yticklabels=[f"布局 {i+1}" for i in range(len(layouts))],
           xticks=[1., 1.2, 1.4, 1.6, 1.8])
    ax.invert_yaxis()
    style(ax, "a  两支中较大的碰撞半带宽", "相对各自的双侧半带宽")
    ax.legend([Line2D([],[],ls="",marker=m,color=c) for m,c in zip(("D","o","o"),SCHEDULE_COLORS)],
              SCHEDULE_LABELS, loc="upper left", bbox_to_anchor=(-.03,-.19),
              ncol=3, fontsize=8.3, frameon=False, columnspacing=1.00, handletextpad=.18)

    ax = fig.add_subplot(gs[0,1])
    counts = [sum(r["decision"]["conditional_TRF_correct_branch"] for row in layouts
                  for r in row["runs"] if r["pattern"] == p and r["eta"] == 1.)
              for p in SCHEDULES]
    total = 2*len(layouts)
    y = np.arange(len(counts))
    ax.barh(y, [total]*len(counts), color="#e4ecef", height=.52)
    ax.barh(y, counts, color=SCHEDULE_COLORS, height=.52)
    for j,n in enumerate(counts):
        ax.text(n-.24,j,f"{n}/{total}",va="center",ha="right",color="white",fontweight="bold",fontsize=11)
    ax.set(yticks=y, yticklabels=SCHEDULE_LABELS, xlim=(0,total+.05), xticks=np.arange(0,total+1,3))
    ax.invert_yaxis()
    style(ax, "b  同半带宽下恢复到生成支", "成功数（两支分别生成读数）")

    ax = fig.add_subplot(gs[1,0])
    locations = [np.array([-1/np.sqrt(2),1/np.sqrt(2)]), np.array([-1.,0.,1.])]
    powers = [np.array([.5,.5]), np.array([.25,.50,.25])]
    colors = (BLUE,PURPLE)
    names = ("带内双侧", "带边双侧加重复旧频")
    for x,w,c in zip(locations,powers,colors):
        ax.vlines(x, 0, w, color=c, linewidth=2.6)
        ax.scatter(x,w,color=c,s=44,zorder=4)
    ax.set(xlim=(-1.17,1.17), ylim=(0,.58), xticks=(-1.,-.707,0.,.707,1.),
           xticklabels=("-1","-0.707","0","0.707","1"), yticks=(0.,.25,.5))
    style(ax, "c  两种谱具有相同的一阶、二阶矩", "频偏 / 最大允许频偏", "新增能量占比", grid="y")
    ax.legend([Line2D([],[],color=c,marker="o",lw=1.5) for c in colors], names,
              loc="upper left", bbox_to_anchor=(-.03,-.26), ncol=2, fontsize=8.6,
              frameon=False, handlelength=1.15, columnspacing=1.03)

    ax = fig.add_subplot(gs[1,1])
    proof = next(row for row in read("spectral_moments.json") if row["label"]=="fold_0.4")
    patterns = ("inner_symmetric","edges_and_repeats")
    for iy,eta in enumerate((1.,.25)):
        points = [worst_band(proof, p, eta)*1e4 for p in patterns]
        for j,(x,c) in enumerate(zip(points,colors)):
            yval = iy+(-.085 if j==0 else .085)
            ax.scatter(x,yval,s=48,color=c,zorder=4)
        ax.vlines(np.mean(points), iy-.085, iy+.085, color=GRAY, lw=1.2)
        ax.text(np.mean(points)+.16,iy,f"{np.mean(points):.2f}",ha="left",va="center",color=INK,fontsize=10)
    ax.set(ylim=(-.31,1.33),xlim=(0,8.5),yticks=(0,1),yticklabels=("增益已知","剖面后 η=0.25"),
           xticks=(0,2,4,6,8))
    ax.invert_yaxis()
    style(ax, "d  等矩频谱的数值碰撞半带宽", "两支最坏阈值（万分之一）")
    fig.text(.065,.028,"图 c 合并显示四次新增曝光的同频能量。图 d 只改变增益校准条件。",
             fontsize=9.2,color=MUTED)
    save(fig,"spectral_moments")


def basin_figure():
    ring,shell = [read(f"face_{case}.json") for case in ("ring","shell")]
    fig = plt.figure(figsize=(10.7,7.45),facecolor="white")
    gs = fig.add_gridspec(2,2,left=.17,right=.965,top=.81,bottom=.17,width_ratios=(1.03,1),
                          hspace=.69,wspace=.32)
    header(fig, "首个边界释放可以预测，材料恢复仍依赖后续路径",
           "576 个材料参数。源权重延拓沿同一驻点支追踪，材料误差只对已完成的终点评分。")

    ax = fig.add_subplot(gs[:,0])
    for block,case in enumerate((ring,shell)):
        for i,key in enumerate(FACIAL_CONTROLS):
            r = case["rows"][key]
            p = r["pred_t"]*1e5
            obs = np.mean(r["observed_t"])*1e5
            yval=block*5+i
            color = BLUE if key=="profile" else GRAY
            ax.hlines(yval-.04,min(obs,p),max(obs,p),color=color,lw=1.5)
            ax.scatter(obs,yval-.13,s=38,color=color,zorder=5)
            ax.scatter(p,yval+.13,s=36,facecolors="white",edgecolors=color,lw=1.4,zorder=5)
            if key=="profile":
                ax.annotate(f"索引 {r['observed_index']}",
                            (obs,yval),xytext=(19,0),textcoords="offset points",va="center",
                            fontsize=9.4,color=BLUE)
    y = [i+offset for offset in (0,5) for i in range(len(FACIAL_CONTROLS))]
    ax.set(xlim=(.49,16.),ylim=(-.7,8.7),xscale="log",yticks=y,
           yticklabels=[f"{case} · {text}" for case in ("ring","shell") for text in FACIAL_LABELS],
           xticks=(.5,1,2,5,10))
    ax.xaxis.set_major_locator(FixedLocator([.5,1,2,5,10]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}"))
    ax.invert_yaxis()
    ax.axhline(4,color=GRID,lw=.8)
    style(ax, "a  首次下界释放", "源权重混合参数 t（十万分之一）")
    ax.legend([Line2D([],[],ls="",marker="o",mfc="white",color=GRAY),
               Line2D([],[],ls="",marker="o",color=GRAY)],
              ["一阶预测","数值延拓"],loc="upper left",bbox_to_anchor=(.03,-.125),
              ncol=2,fontsize=9,frameon=False,handlelength=1.2,columnspacing=1.0)

    ax = fig.add_subplot(gs[0,1])
    controls = read("basin_scores.json")
    scores = {(r["case"],r["method"]):r for r in controls if r["file"].startswith("basin_controls_")}
    entries = [
        ("ring","initial","o",GOLD, None, (0,0)),
        ("ring","profile_long","D",GOLD,None,(0,0)),
        ("ring","cold_midbox_long","^",GOLD,"ring 中点冷启长跑",(8,-3)),
        ("shell","initial","o",BLUE,None,(0,0)),
        ("shell","profile","s",BLUE,"剖面短校正",(0,-13)),
        ("shell","profile_long","D",BLUE,"剖面长校正",(0,12)),
    ]
    for case,method,marker,color,label,offset in entries:
        result=scores[case,method]
        x = result["phi"]/scores[case,"initial"]["phi"]
        y = result["NRMSE"]
        hollow = method == "initial"
        ax.scatter(x,y,s=54 if hollow else 33,marker=marker,
                   facecolors="white" if hollow else color,edgecolors=color,
                   linewidth=1.2,zorder=4 if hollow else 5)
        if label:
            ax.annotate(label,(x,y),xytext=offset,textcoords="offset points",
                        color=color,ha="left" if offset[0] else "center",va="center",fontsize=8.5)
    ax.set(xscale="log",xlim=(.024,1.28),ylim=(.33,1.065),yticks=(.4,.6,.8,1.))
    ax.xaxis.set_major_locator(FixedLocator([.03,.1,.3,1.]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}"))
    style(ax, "b  旧数据拟合与材料误差", "拟合目标 / 旧起点目标", "材料 NRMSE", grid="both")
    ax.legend([Line2D([],[],ls="",marker="o",mfc="white",color=c) for c in (GOLD,BLUE)],
              ["ring 起点","shell 起点"],loc="upper right",fontsize=8.4,frameon=False,
              handletextpad=.25, borderaxespad=.65,labelspacing=.35)

    ax = fig.add_subplot(gs[1,1])
    cases = ("ring","shell","bridge")
    methods = ("wide_direct_bad","wide_cumulative_bad","narrow_cumulative_bad")
    labels = ("宽频同时","宽频累计","窄频累计")
    for i,case in enumerate(cases):
        result = read(f"freq_high/{case}_wide_7200_scores.json")
        initial = read(f"freq_high/{case}_scores.json")["bad"]["relative_L2"]
        ax.hlines(i,0,initial,color="#becbd1",lw=1.3,zorder=1)
        ax.scatter(initial,i,s=42,marker="|",c=INK,zorder=3)
        for j,(key,c,m) in enumerate(zip(methods,(GOLD,BLUE,GRAY),("o","D","x"))):
            ax.scatter(result[key],i+(j-1)*.13,s=42,color=c,marker=m,zorder=4)
    ax.set(yticks=np.arange(len(cases)),yticklabels=cases,ylim=(-.5,2.57),xlim=(0,1.13),
           xticks=(0,.25,.5,.75,1.))
    ax.invert_yaxis()
    style(ax, "c  旧坏起点之后的材料恢复", "相对材料误差（7200 个复 RHS 上限）")
    handles=[Line2D([],[],ls="",marker=m,color=c) for m,c in zip(("o","D","x","|"),(GOLD,BLUE,GRAY,INK))]
    ax.legend(handles,[*labels,"旧起点"],loc="upper left",bbox_to_anchor=(-.025,-.275),
              ncol=4,fontsize=8.1,frameon=False,columnspacing=.72,handletextpad=.20,
              handlelength=.90)
    fig.text(.065,.028,"图 b 使用同一批旧频数据，NRMSE 按材料幅值归一。图 c 加测两次半功率曝光，宽频累计按 -50%、-25%、基频推进。",
             color=MUTED,fontsize=9.1)
    save(fig,"basin_release")


if __name__ == "__main__":
    preferred = ("Noto Sans CJK SC", "Noto Sans SC", "PingFang SC", "STHeiti", "Hiragino Sans GB")
    fonts=set(get_font_names())
    chinese = next((x for x in preferred if x in fonts),None)
    if chinese is None:
        raise RuntimeError("Install a Simplified Chinese font, for example Noto Sans CJK SC.")
    with plt.rc_context({"font.family":["DejaVu Sans",chinese],"font.size":10,
                         "mathtext.fontset":"dejavusans","axes.unicode_minus":False,
                         "pdf.fonttype":3}):
        spectral_figure()
        basin_figure()
