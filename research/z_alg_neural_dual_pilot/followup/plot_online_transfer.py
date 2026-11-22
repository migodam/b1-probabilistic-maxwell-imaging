"""Plot the paired native transfer comparison from its saved scoring receipt."""

import json
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib.patches import Patch, Rectangle


HERE = Path(__file__).resolve().parent
FIXED = "fixed_archived_five"
GRADIENT = "native"
INK = "#21313d"
COLORS = {
    (True, True): ("#2b609e", "＝", "共同过线"),
    (False, True): ("#258273", "＋", "仅梯度过线"),
    (False, False): ("#e7edf1", "", "均未过线"),
}


def setup_font():
    available = {font.name for font in font_manager.fontManager.ttflist}
    preferred = ("Heiti SC", "Microsoft YaHei", "Noto Sans CJK SC", "Hiragino Sans GB", "PingFang SC")
    font = next((name for name in preferred if name in available), None)
    if font is None:
        raise RuntimeError("A Simplified Chinese sans-serif font is required for this figure.")
    heading = 400 if font == "Heiti SC" else "bold"
    plt.rcParams.update({
        "font.family": [font, "DejaVu Sans"], "font.weight": 400,
        "pdf.fonttype": 42 if font in ("Heiti SC", "Microsoft YaHei") else 3, "text.color": INK,
    })
    return heading


def view_matrix(ax, rows, objects, noise_levels):
    for i, oid in enumerate(objects):
        for j, noise in enumerate(noise_levels):
            row = rows[oid, noise]
            if any(row[method]["verdict"] not in ("PASS", "BELOW_OR_EQUAL") for method in (FIXED, GRADIENT)):
                raise ValueError("The paired view requires a separate unresolved category.")
            status = tuple(row[method]["verdict"] == "PASS" for method in (FIXED, GRADIENT))
            if status == (True, False):
                raise ValueError("A fixed-only pass requires a separate paired category.")
            color, mark, _ = COLORS[status]
            ax.add_patch(Rectangle((j - 0.46, i - 0.42), 0.92, 0.84, facecolor=color, linewidth=0))
            if mark:
                ax.text(j, i, mark, fontsize=9, color="white", ha="center", va="center")
    ax.set_xlim(-0.5, len(noise_levels) - 0.5)
    ax.set_ylim(len(objects) - 0.5, -0.5)
    ax.set_xticks(range(len(noise_levels)), [f"{noise:g}%" for noise in noise_levels], fontsize=10)
    ax.set_yticks(range(len(objects)), [oid.rsplit("_", 1)[-1] for oid in objects], fontsize=8.5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=6, labelcolor=INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(
        handles=[Patch(facecolor=color, label=label) for color, _, label in COLORS.values()],
        loc="lower left", bbox_to_anchor=(-0.024, -0.142), ncols=3,
        frameon=False, fontsize=8.5, handlelength=1.4, columnspacing=1.1,
    )


def comparison(ax, receipt, rows, noise_levels, heading):
    ax.set_axis_off()
    all_views = receipt["summary"]["all_views"]
    n_views = all_views["views"]
    pass_rows = [row for row in rows.values() if row[GRADIENT]["verdict"] == "PASS"]
    n_objects = all_views["independent_objects"]
    fixed_passes = len(all_views["pass_views"][FIXED])
    fixed_objects = len(all_views["pass_objects"][FIXED])
    gradient_passes = len(pass_rows)
    gradient_objects = len(all_views["pass_objects"][GRADIENT])
    fixed_med = all_views["median_online"][FIXED]["wall_seconds"]
    gradient_med = all_views["median_online"][GRADIENT]["wall_seconds"]
    success_med = median(row[GRADIENT]["wall_seconds"] for row in pass_rows)
    fixed_prep = receipt["fixed_prep"]["matched_response_preparation"]["wall_seconds"]

    ax.text(0, 0.997, "严格过线数", fontsize=12, weight=heading)
    ax.text(0.53, 0.929, "旧五固定权重", ha="right", fontsize=9.2)
    ax.text(1, 0.929, "64 区普通梯度", ha="right", fontsize=9.2, color=COLORS[(False, True)][0])
    entries = [(f"全部 {n_views} 个视图", f"{fixed_passes}/{n_views}", f"{gradient_passes}/{n_views}"),
               (f"至少一档过线的材料", f"{fixed_objects}/{n_objects}", f"{gradient_objects}/{n_objects}")]
    for noise in noise_levels:
        section = receipt["summary"][f"{noise:g}_pct"]
        n = section["views"]
        entries.append((f"{noise:g}% 噪声", f"{len(section['pass_views'][FIXED])}/{n}",
                        f"{len(section['pass_views'][GRADIENT])}/{n}"))
    for i, (label, fixed, gradient) in enumerate(entries):
        y = 0.846 - 0.078 * i
        ax.text(0, y, label, fontsize=9.2, va="center")
        ax.text(0.53, y, fixed, ha="right", fontsize=10, va="center")
        ax.text(1, y, gradient, ha="right", fontsize=10, va="center", color=COLORS[(False, True)][0])
        ax.plot([0, 1], [y - 0.035, y - 0.035], color="#dce5e9", linewidth=0.65)

    ax.text(0, 0.385, "单视图在线墙钟中位数", fontsize=11, weight=heading)
    costs = [
        (f"固定表评分 · 全 {n_views} 视图", f"{fixed_med * 1000:.2f} 毫秒"),
        (f"梯度搜索＋核验 · 全 {n_views} 视图", f"{gradient_med:.2f} 秒"),
        (f"梯度搜索＋核验 · 仅 {gradient_passes} 个过线视图", f"{success_med:.2f} 秒"),
    ]
    for i, (label, cost) in enumerate(costs):
        y = 0.306 - 0.071 * i
        ax.text(0, y, label, fontsize=9.2, va="center")
        ax.text(1, y, cost, ha="right", fontsize=9.2, va="center")
    ax.text(0, 0.083, f"固定表另有 {fixed_prep:.2f} 秒一次性准备。", fontsize=8.7)
    ax.text(0, 0.02, "1 秒是粗网搜索预算。在线费用计入原生核验和补证。", fontsize=8.7)
    ax.set(xlim=(0, 1), ylim=(0, 1))


def main():
    receipt = json.loads((HERE / "results" / "online_transfer_paired.json").read_text())
    rows = {(row["object_id"], row["noise_pct"]): row for row in receipt["rows"]}
    objects = sorted({oid for oid, _ in rows})
    noise_levels = sorted({noise for _, noise in rows})
    every_pair = {(oid, noise) for oid in objects for noise in noise_levels}
    if set(rows) != every_pair or len(rows) != len(receipt["rows"]):
        raise ValueError("The receipt must contain exactly one observation for each material and noise level.")
    heading = setup_font()
    fig = plt.figure(figsize=(9.5, 5.65), facecolor="white")
    gains = len(receipt["summary"]["all_views"]["gain_over_fixed_bound"])
    fig.text(0.048, 0.956, f"普通梯度多证实 {gains} 个视图可以排除材料上限",
             fontsize=17, weight=heading, va="top")
    fig.text(0.048, 0.881,
             f"{len(objects)} 件新材料 · 每件 {len(noise_levels)} 个噪声档"
             f" · 96² 原生核验 · 峰值上限 {receipt['cap']:g}", fontsize=10)
    matrix = fig.add_axes((0.068, 0.175, 0.372, 0.618))
    fig.text(0.048, 0.825, "材料", fontsize=9)
    fig.text(0.254, 0.844, "噪声档", fontsize=9, ha="center")
    view_matrix(matrix, rows, objects, noise_levels)
    comparison(fig.add_axes((0.5, 0.15, 0.453, 0.692)), receipt, rows, noise_levels, heading)
    fig.text(0.048, 0.039, "过线表示证书距离的严格下界超过该视图的原噪声半径。",
             fontsize=9)
    for ext in ("png", "pdf"):
        fig.savefig(HERE / "figures" / f"online_transfer.{ext}", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
