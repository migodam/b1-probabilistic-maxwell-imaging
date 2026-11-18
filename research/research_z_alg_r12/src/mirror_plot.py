"""Render the mirror-cell phase control from the recorded Maxwell results."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.font_manager import get_font_names
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INK = "#162938"
MUTED = "#475761"
GRID = "#e2e9ec"
BLUE = "#12769a"
GOLD = "#aa5e13"
PURPLE = "#83558d"
GRAY = "#6a7880"


def style(ax, title, xlabel, ylabel=None, grid="both"):
    ax.set_title(title, loc="left", pad=11, fontweight="bold", color=INK, fontsize=12)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=GRID, lw=.7)
    ax.spines[["top","right"]].set_visible(False)
    ax.spines[["bottom","left"]].set_color("#a8bac2")
    ax.tick_params(length=3,width=.65,colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def geometry_panel(ax, geo):
    cells = np.array(geo["centers"])
    sources = np.array(geo["sources"])
    receivers = np.array(geo["receivers"])
    ax.axvline(0, color="#97abb6", ls="--", lw=1.1)
    for i in range(0,len(sources),2):
        beam=sources[i:i+2]
        ax.plot(*beam.T, color=BLUE, alpha=.55, lw=1.0)
    for p in cells[2:]:
        ax.add_patch(Circle(p,geo["radius"],facecolor="#c6d1d6",edgecolor=GRAY,lw=.8,zorder=5))
    for p in cells[:2]:
        ax.add_patch(Circle(p,geo["radius"],facecolor="#e7a451",edgecolor=GOLD,lw=1.0,zorder=5))
    ax.scatter(*sources.T,marker="^",s=34,c=BLUE,zorder=6)
    ax.scatter(*receivers.T,marker="s",s=27,c=PURPLE,zorder=6)
    ax.set(xlim=(-1.09,1.09),ylim=(-1.21,1.31),xticks=(-1,0,1),yticks=(-1,0,1),aspect="equal")
    style(ax,"a  两个活动材料单元","横坐标（模型长度）","纵坐标")
    handles = [Patch(fc="#e7a451",ec=GOLD,label="未知材料"),
               Patch(fc="#c6d1d6",ec=GRAY,label="已知宿主"),
               Line2D([],[],marker="^",color=BLUE,ls="-",label="成对源端口"),
               Line2D([],[],marker="s",color=PURPLE,ls="",label="轴上接收")]
    ax.legend(handles=handles,ncol=2,loc="upper center",bbox_to_anchor=(.49,-.20),
              fontsize=8.2,frameon=False,columnspacing=1.05,handlelength=1.0,handletextpad=.35)


def signal_panel(ax, data):
    for key,color,label,marker in [
        ("equal_feeds",GRAY,"等幅同相馈电","o"),
        ("physical_delays",BLUE,"馈电相位随频偏变化","D"),
    ]:
        points = sorted((100*float(b),gap) for b,gap in data[key]["absolute_gaps"].items())
        d,v=np.array(points).T
        ax.scatter(d,v,marker=marker,color=color,s=27,zorder=3,label=label)
    ax.set(xlim=(-55,55),ylim=(9e-18,1e-1),yscale="log",
           xticks=(-50,-25,0,25,50),yticks=(1e-16,1e-12,1e-8,1e-4,1e-2))
    ax.axvline(0,color="#a8bac2",ls=":",lw=1.)
    style(ax,"b  交换材料后的读数差","相对频偏（%）","复读数差的二范数")
    ceiling = max(data["equal_feeds"]["absolute_gaps"].values())
    ax.annotate(f"等馈最大 {ceiling/1e-17:.1f}×10⁻¹⁷",(-50,1.3e-16),color=MUTED,fontsize=8.6)
    ax.legend(loc="upper center",bbox_to_anchor=(.51,-.20),ncol=1,frameon=False,
              fontsize=8.4,handletextpad=.28,labelspacing=.25)


def collision_panel(ax, data):
    codes = ("opposite_edges","two_negative","two_positive")
    labels = ("正负各一次","两次负偏","两次正偏")
    entries = {(r["experiment"],r["pattern"],r["side"]):r
               for r in data["runs"] if r["gain"]=="known" and r["a"]==.2}
    plans = (("完整公式",BLUE,"o"),("略去附加法向项",GOLD,"s"),
             ("固定传播，仅改相位",PURPLE,"D"))
    full=[]; omitted=[]; phase=[]
    for i,code in enumerate(codes):
        for side in (-1,1):
            phys=entries["physical_delays",code,side]
            control=entries["fixed_frequency_phase_control",code,side]
            b=phys["numerical_band"]
            p=100*np.abs(np.r_[phys["predicted_band"],
                               phys["predicted_omitting_extra_normal_components"],
                               control["numerical_band"]]/b-1.)
            full.append(p[0]); omitted.append(p[1]); phase.append(p[2])
            for j, (v,(_,c,marker)) in enumerate(zip(p,plans)):
                ax.scatter(v,i+(j-1)*.18+side*.043,color=c,marker=marker,s=26,zorder=4)
    major=(.2,.5,1.,5.,10.,40.)
    ax.set(xscale="log",xlim=(.18,61.),ylim=(-.48,2.62),
           yticks=np.arange(len(codes)),yticklabels=labels)
    ax.xaxis.set_major_locator(FixedLocator(major))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}"))
    ax.axvline(1,color="#97abb6",ls=":",lw=1.15)
    ax.invert_yaxis()
    style(ax,"c  碰撞阈值的数值对照","相对数值阈值的绝对偏差（%）",grid="x")
    ax.annotate(f"最大 {max(omitted):.1f}%",(42,.67),ha="right",fontsize=8.7,color=GOLD)
    ax.annotate(f"完整 ≤{max(full):.2f}%，相位对照 ≤{max(phase):.2f}%",
                (.24,2.51),fontsize=8.4,color=INK)
    handles=[Line2D([],[],marker=m,color=c,ls="",label=n) for n,c,m in plans]
    ax.legend(handles=handles,loc="upper center",bbox_to_anchor=(.56,-.20),
              frameon=False,ncol=2,fontsize=8.1,columnspacing=.8,handletextpad=.24)


def main():
    data = json.loads((ROOT/"results/mirror_spectrum.json").read_text())
    allowed = ("Noto Sans CJK SC","Noto Sans SC","PingFang SC","STHeiti","Hiragino Sans GB")
    available = set(get_font_names())
    font = next((x for x in allowed if x in available),None)
    if font is None:
        raise RuntimeError("Install a Simplified Chinese font such as Noto Sans CJK SC.")
    with plt.rc_context({"font.family":["DejaVu Sans",font],"font.size":9.7,
                         "axes.unicode_minus":False,"pdf.fonttype":3}):
        fig = plt.figure(figsize=(12.6,4.15),facecolor="white")
        grid = fig.add_gridspec(1,3,left=.07,right=.985,top=.70,bottom=.306,
                               width_ratios=(1.,1.12,1.53),wspace=.43)
        fig.text(.049,.957,"镜像材料的扫频简并与馈电相位",va="top",fontsize=19,weight="bold",color=INK)
        fig.text(.05,.854,"六个 Maxwell 离散单元中有两个未知材料参数。两组镜像源端口相干馈电，三个接收器位于对称轴。",
                 va="top",fontsize=10.3,color=MUTED)
        geometry_panel(fig.add_subplot(grid[0,0]),data["geometry"])
        signal_panel(fig.add_subplot(grid[0,1]),data)
        collision_panel(fig.add_subplot(grid[0,2]),data)
        fig.text(.05,.035,"图 b 使用 a=0.60 的镜像材料。图 c 取 a=0.20、已知增益，两支分别测试，以完整变频阈值为参照。馈电保持等功率。",
                 fontsize=9.1,color=MUTED)
        target = ROOT/"figures"/"mirror_spectrum"
        fig.savefig(target.with_suffix(".png"),dpi=240,facecolor="white")
        fig.savefig(target.with_suffix(".pdf"),facecolor="white")
        plt.close(fig)


if __name__ == "__main__":
    main()
