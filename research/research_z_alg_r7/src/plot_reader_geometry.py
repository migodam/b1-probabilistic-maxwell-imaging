"""Draw the R7 primary square-target geometry from the archived design."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle


ROOT = Path(__file__).resolve().parents[1]
INK = "#192530"
MUTED = "#475968"
SOURCE = "#535B9B"
RECEIVER = "#192530"
TARGET = "#AD4C25"
NEIGHBOR = "#176B68"
BACKGROUND = "#788895"


def chinese_font(weight="regular"):
    mac_font = Path("/System/Library/Fonts") / {
        "regular": "STHeiti Light.ttc",
        "semibold": "STHeiti Medium.ttc",
    }[weight]
    if mac_font.is_file():
        return FontProperties(fname=mac_font)
    for family in ("Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Zen Hei"):
        font = FontProperties(family=family, weight=weight)
        try:
            findfont(font, fallback_to_default=False)
        except ValueError:
            continue
        return font
    raise RuntimeError("Install a Simplified Chinese font such as Noto Sans CJK SC.")


def draw_cells(ax, design, target, neighbors, *, annotate):
    radius = design["disk_radius"]
    for index, (x, y) in enumerate(design["cells"]):
        if index in target:
            color, text = TARGET, "q"
        elif index in neighbors:
            color, text = NEIGHBOR, str(neighbors[index])
        else:
            color, text = BACKGROUND, ""
        ax.add_patch(Circle((x, y), radius,
                            facecolor=color if text else "white",
                            edgecolor=color, linewidth=.85, zorder=3))
        if annotate and text:
            ax.text(x, y, text, ha="center", va="center", color="white", fontsize=8.1,
                    fontweight="bold", fontstyle="italic" if text == "q" else "normal",
                    zorder=4)


def key(ax, entries, font, x, rows):
    for (marker, fill, edge, label), y in zip(entries, rows):
        handle = Line2D([x], [y], marker=marker, markersize=7.5,
                        markerfacecolor=fill, markeredgecolor=edge,
                        markeredgewidth=1.0, linestyle="none",
                        transform=ax.transAxes)
        ax.add_line(handle)
        ax.text(x + .028, y, label, transform=ax.transAxes, va="center",
                fontproperties=font, fontsize=12.3, color=INK)


def main():
    design = json.loads((ROOT / "results/pencil_raw.json").read_text())["design"]["primary"]
    active = design["active"]
    target = {active[index] for index in design["regions"]["square"]["target"]}
    neighbors = {index: label for label, index in enumerate(
        (index for index in active if index not in target), start=1)}

    regular = chinese_font()
    semibold = chinese_font("semibold")
    plt.rcParams.update({"pdf.fonttype": 42, "figure.facecolor": "white"})
    fig = plt.figure(figsize=(10.5, 5.15))
    fig.text(.057, .925, "外部发射、单点接收与待测材料",
             fontproperties=semibold, fontsize=24, color=INK)

    overview = fig.add_axes([.088, .227, .345, .614])
    detail = fig.add_axes([.558, .227, .345, .614])
    for ax, limit in ((overview, 1.99), (detail, .96)):
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal")
        ax.axis("off")

    fig.text(.256, .833, "主实验的发射与接收位置", ha="center",
             fontproperties=semibold, fontsize=14, color=INK)
    fig.text(.731, .833, f'{len(design["cells"])} 个计算单元（放大）', ha="center",
             fontproperties=semibold, fontsize=14, color=INK)

    draw_cells(overview, design, target, neighbors, annotate=False)
    draw_cells(detail, design, target, neighbors, annotate=True)
    overview.scatter(*zip(*design["sources"]), marker="^", s=66, color=SOURCE, zorder=5)
    overview.scatter(*zip(*design["receiver"]), marker="D", s=66, color=RECEIVER, zorder=5)

    x, y = zip(*design["cells"])
    margin = design["disk_radius"] + .06
    overview.add_patch(Rectangle(
        (min(x) - margin, min(y) - margin),
        max(x) - min(x) + 2 * margin, max(y) - min(y) + 2 * margin,
        fill=False, edgecolor=BACKGROUND, linewidth=.95, linestyle=(0, (4, 4)), zorder=2))

    legend = fig.add_axes([0, 0, 1, 1])
    legend.axis("off")
    key(legend, [
        ("^", SOURCE, SOURCE, f'{len(design["sources"])} 个外部发射位置'),
        ("D", RECEIVER, RECEIVER, f'{len(design["receiver"])} 个接收器'),
    ], regular, x=.127, rows=(.165, .111))
    key(legend, [
        ("o", TARGET, TARGET, f"{len(target)} 个待测单元，共用一个材料值 q"),
        ("o", NEIGHBOR, NEIGHBOR, f"{len(neighbors)} 个周围单元，各有独立的未知材料值"),
        ("o", "white", BACKGROUND, f'{len(design["cells"]) - len(active)} 个背景单元，材料值已知'),
    ], regular, x=.57, rows=(.165, .111, .057))

    for extension in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"reader_geometry.{extension}", dpi=250)
    plt.close(fig)


if __name__ == "__main__":
    main()
