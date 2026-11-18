"""Plot the stored scalar-TM uncertainty certificates and point witnesses."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
INK = "#162e40"
MUTED = "#496172"
OLD = "#bacbd6"
NEW = "#168480"
ALTERNATIVE = "#a35032"
TRUTH = "#16466a"
PAPER = "#ffffff"
GRID = "#ecf0f3"


def setup():
    font_paths = [
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    ]
    for path in font_paths:
        if Path(path).exists():
            font_manager.fontManager.addfont(path)
            family = font_manager.FontProperties(fname=path).get_name()
            break
    else:
        family = "Noto Sans CJK SC"
    plt.rcParams.update({
        "font.family": family, "font.size": 10, "text.color": INK,
        "axes.labelcolor": INK, "axes.edgecolor": "#bcc9d1", "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "pdf.fonttype": 42, "axes.unicode_minus": False,
        "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    })


def read(name):
    return json.loads((RESULTS / name).read_text())


def clean_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=3)
    ax.set_axisbelow(True)


def shade_mask(ax, cells, box, levels, color):
    mask = np.zeros((levels, levels), dtype=bool)
    mask[cells[:, 1], cells[:, 0]] = True
    plt_color = ListedColormap([color])
    ax.imshow(np.ma.masked_where(~mask, mask), extent=[*box, *box], origin="lower",
              interpolation="none", cmap=plt_color, vmin=0, vmax=1, zorder=1)


def two_coordinate_case(fig, axes, data, cells, truth):
    case = "fresh6-00"
    before, after = (data[f"{case}/{key}"] for key in ("old", "active"))
    assert before["pivots"][:2] == after["pivots"] and before["depth"] == after["depth"]
    box, levels = after["box"], 2**after["depth"]
    pivots = after["pivots"]
    witness = np.asarray(after["witness_coordinates"])
    truth_point = truth[pivots]
    assert np.allclose(witness[0], truth_point)
    for ax, key, color, title in zip(axes, ("old", "active"), (OLD, NEW),
                                      ("3 束旧照明", "4 束 · 追加主动照明")):
        projected = np.unique(cells[f"{case}/{key}"][:, :2], axis=0)
        shade_mask(ax, projected, box, levels, color)
        for j, point in enumerate(witness[1:], 1):
            removed = key == "active" and after["witnesses"][j]["outside"]
            ax.scatter(*point, s=63, marker="x" if removed else "D",
                       color=ALTERNATIVE, edgecolor=PAPER if not removed else None,
                       linewidth=1.65 if removed else .85, zorder=4)
        ax.scatter(*truth_point, marker="*", s=210, color=TRUTH,
                   edgecolor=PAPER, linewidth=1, zorder=5)
        ticks = (5, 10, 15, 20, 25, 30)
        ax.set(xlim=box, ylim=box, xticks=ticks, yticks=ticks, title=title,
               xlabel=f"单元 {pivots[0]+1} 的对比度")
        ax.set_aspect("equal")
        ax.grid(color=GRID, linewidth=.8, zorder=0)
        clean_axes(ax)
    axes[0].set_ylabel(f"单元 {pivots[1]+1} 的对比度")
    axes[1].tick_params(labelleft=False, left=False)
    handles = [
        Line2D([], [], marker="*", linestyle="None", ms=12, color=TRUTH, label="生成真值"),
        Line2D([], [], marker="D", linestyle="None", ms=6, color=ALTERNATIVE, label="旧可行见证"),
        Line2D([], [], marker="x", linestyle="None", ms=7, color=ALTERNATIVE, label="被第 4 束排除的见证"),
    ]
    fig.legend(handles=handles, ncol=3, frameon=False, fontsize=9, handletextpad=.25,
               columnspacing=1.2, loc="upper left", bbox_to_anchor=(.059, .467), borderaxespad=0)


def twenty_four_case(ax, cert):
    choices = [("branch", "按分支选源"), ("geometry_projection", "按几何投影选源"), ("random", "随机选源")]
    control = cert["24/branch"]
    truth, alternative = np.array(control["witness_coordinates"]).ravel()
    assert all(w["inside"] for w in control["old_witnesses"])
    ax.scatter(truth, 3, marker="*", s=180, color=TRUTH, linewidth=.9, edgecolor=PAPER, zorder=5)
    ax.scatter(alternative, 3, marker="D", s=52, color=ALTERNATIVE, linewidth=.8, edgecolor=PAPER, zorder=5)
    for x in (alternative, truth):
        ax.text(x, 3.23, f"{x:.3f}", ha="center", fontsize=9)
    for row, (key, _) in enumerate(choices):
        y = 2 - row
        for lo, hi in cert[f"24/{key}"]["outer"][0]:
            ax.plot((lo, hi), (y, y), color=NEW, lw=6, solid_capstyle="butt")
            ax.vlines((lo, hi), y-.075, y+.075, color=NEW, lw=1.0)
            lo_display = np.floor(lo*1000)/1000
            hi_display = np.ceil(hi*1000)/1000
            ax.text(hi+.26, y, f"{lo_display:.3f}–{hi_display:.3f}", va="center", fontsize=8.1, color=MUTED)
    ax.set(xlim=(.5, 12), ylim=(-.42, 3.59), xticks=(1, 3, 5, 7, 9, 11), yticks=(3, 2, 1, 0),
           yticklabels=["22 束 · 旧见证", *(label for _, label in choices)],
           xlabel=f"单元 {control['pivots'][0]+1} 的对比度")
    ax.grid(axis="x", color=GRID, linewidth=.8)
    clean_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False)


def objective(witness):
    return float(witness["objective_enclosure"].removeprefix("[").split()[0])


def confirmation(ax, records, case):
    row = records[f"{case}/active"]
    ax.axvspan(.56, 1, color="#e1f0eb", zorder=0)
    ax.axvline(1, color="#617786", linewidth=1.2, linestyle="--", zorder=1)
    for i, (previous, subsequent) in enumerate(zip(row["old_witnesses"], row["new_witnesses"])):
        old, new = objective(previous), objective(subsequent)
        y = 1-i
        color = TRUTH if i == 0 else ALTERNATIVE
        marker, size = ("*", 160) if i == 0 else ("D", 62)
        ax.plot((old, new), (y, y), color=color, lw=1.4, alpha=.66, zorder=2)
        ax.scatter(old, y, marker=marker, s=size, facecolor=PAPER, edgecolor=color, linewidth=1.0, zorder=3)
        ax.scatter(new, y, marker=marker, s=size, facecolor=color, edgecolor=PAPER, linewidth=.85, zorder=4)
        ax.annotate(f"{new:.3f}", xy=(new, y), xytext=(8, 6), textcoords="offset points",
                    fontsize=9, color=color, va="bottom")
    alternative = row["new_witnesses"][1]
    gap = row["all_charts_linf_radii"][1]
    title = "仍有两个可行见证" if alternative["inside"] else f"另一见证被排除 · 邻域半径 {gap:.3f}"
    ax.set(xscale="log", xlim=(.56, 11.0), ylim=(-.46, 1.48), yticks=(1, 0),
           yticklabels=("生成真值", "另一旧见证"), title=f"{case} · {title}")
    ax.xaxis.set_major_locator(FixedLocator([.6, 1, 2, 4, 8]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
    ax.minorticks_off()
    ax.grid(axis="x", color=GRID, linewidth=.8)
    clean_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False)


def main():
    setup()
    fig = plt.figure(figsize=(13.25, 8.35))
    grid = fig.add_gridspec(2, 2, left=.070, right=.973, top=.838, bottom=.113,
                            width_ratios=(1.40, 1), height_ratios=(1.56, 1),
                            hspace=.78, wspace=.26)
    top = grid[0, 0].subgridspec(1, 2, wspace=.16)
    ax_old, ax_new = (fig.add_subplot(top[0, i]) for i in range(2))
    ax_short = fig.add_subplot(grid[0, 1])
    confirm = grid[1, :].subgridspec(1, 2, wspace=.26)
    confirm_axes = [fig.add_subplot(confirm[0, i]) for i in range(2)]
    fig.text(.059, .955, "追加一次照明后，哪些材料仍与读数相容？", fontsize=22,
             ha="left", va="top", fontweight="bold")
    fig.text(.060, .905,
             "标量 TM 全波模拟 · 联合误差界  |  复加性误差半径 2 × 10⁻⁴，公共接收增益误差半径 1%",
             fontsize=10, color=MUTED)
    fig.text(.059, .871, "a   自然六单元 · fresh6-00", fontsize=13, fontweight="bold")
    fig.text(.610, .871, "b   二十四单元 · 22 → 23 束", fontsize=13, fontweight="bold")
    fig.text(.610, .447, "上行仅标已知可行点。下三行色条为追加后的严格外包。", fontsize=9, color=MUTED)
    fig.text(.059, .429,
             "色块为严格外包，两图均为 64 × 64 格。色块内可能有不可行材料。",
             fontsize=9, color=MUTED)
    fig.text(.059, .382, "c   独立六单元确认 · 同一误差界", fontsize=13, fontweight="bold")
    fig.text(.07, .066,
             "横轴：解释所示读数需要的最小归一化联合误差平方；虚线为上限 1。空心为 3 束旧数据，实心为追加一束后。",
             fontsize=9, color=MUTED)
    fig.text(.07, .038,
             "确认组的第 4 束使用随机后备照明。邻域半径按材料最大分量距离计，由全部 15 对坐标的外包共同验证。",
             fontsize=9, color=MUTED)

    with np.load(RESULTS/"natural_old_certificate_cells.npz") as cells, \
            np.load(RESULTS/"natural_inputs.npz") as natural:
        two_coordinate_case(fig, (ax_old, ax_new), read("natural_old_certificate.json"),
                            cells, natural["fresh6-00/truth"])
    twenty_four_case(ax_short, read("short_certificate.json"))
    records = read("confirm_all_pairs.json")["records"]
    for ax, case in zip(confirm_axes, ("confirm6-00", "confirm6-05")):
        confirmation(ax, records, case)
    for ext in ("pdf", "png"):
        fig.savefig(ROOT/"figures"/f"certified_regions.{ext}", dpi=230)
    plt.close(fig)


if __name__ == "__main__":
    main()
