"""Render figures and concise machine-readable tables from the frequency results."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np

from freq_fold import DATA, GAIN_SD, ROBUST_SIGMA, visible_gap
from freq_model import FrequencyTM

ROOT = Path(__file__).resolve().parents[1]
C_NEG = "#2166ac"
C_POS = "#b2182b"
COLORS = {"-1": C_NEG, "1": C_POS}
METHODS = ("direct", "pair_multistart", "generic_multistart", "same_read_continuation", "cumulative")


def fork_data(case, sign):
    return case["exact"][sign]["trajectory"], case["exact"][sign]["fold"]


def compare_magnitude(case):
    return next(c["delta"] for c in case["comparisons"][2:] if c["delta"] > 0)


def paper_figure(cases, natural, zh=False):
    def title(en, cn):
        return cn if zh else en
    ref = cases[2]
    model = FrequencyTM(json.loads(DATA.read_text())["geometry"])
    pair = np.array(ref["candidates"])
    gap_curve = np.linspace(.0001,.0022,100)
    positive_gaps = np.array([visible_gap(model,pair,d)[0] for d in gap_curve])
    negative_gaps = np.array([visible_gap(model,pair,-d)[0] for d in gap_curve])

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), layout="constrained")
    ax = axes[0,0]
    for gaps, sign in [(negative_gaps, "-1"), (positive_gaps, "1")]:
        ax.plot(1e3*gap_curve, 1e5*gaps, color=COLORS[sign], lw=2,
                ls="-" if sign == "-1" else "--")
        fold = ref["exact"][sign]["fold"]
        x = fold["delta"]
        ax.scatter(abs(x)*1e3, visible_gap(model,pair,x)[0]*1e5,
                   c=COLORS[sign], zorder=4, marker="X", s=54)
    ax.plot(1e3*gap_curve,1e5*gap_curve*ref["pair_frequency_jet_norm"],
            color="0.38", ls=":", lw=1.5)
    ax.axvline(1e3*compare_magnitude(ref), color="0.5", ls=":", lw=1., zorder=0)
    ax.set(xlabel=title(r"Absolute frequency offset $|\delta|$ ($10^{-3}$)", "频偏绝对值 |δ|（千分之一）"),
           ylabel=title(r"Candidate response gap ($10^{-5}$)", "两候选读数差（十万分之一）"),
           title=title("a   Equal-power, signed measurements", "a   等功率正负频偏"))
    ax.legend(handles=[Line2D([0],[0],color=C_NEG,label=title(r"$\delta<0$","负频偏")),
                       Line2D([0],[0],color=C_POS,ls="--",label=title(r"$\delta>0$","正频偏")),
                       Line2D([0],[0],color="0.38",ls=":",label=title("First-order response", "读数差一阶预测")),
                       Line2D([0],[0],color="0.38",lw=0,marker="X",label=title("Numerical fold", "数值碰撞"))],
              fontsize=9, frameon=False, loc="upper left")

    ax = axes[0,1]
    for sign in ("-1","1"):
        path, fold = fork_data(ref, sign)
        for name in ("minimum","saddle"):
            x = [1e3*abs(p["delta"]) for p in path] + [1e3*abs(fold["delta"])]
            y = [p[name]["t"] for p in path] + [fold["t"]]
            ax.plot(x,y, color=COLORS[sign], ls="-" if name=="minimum" else "--",
                    lw=1.8, marker=".", ms=4)
        ax.scatter(1e3*abs(fold["delta"]),fold["t"], c=COLORS[sign],marker="X",zorder=4,s=54)
        ax.axvline(1e3*abs(ref["predicted"][sign]), color=COLORS[sign], ls=":", lw=1)
    ax.axvline(1e3*compare_magnitude(ref), color="0.5", ls=":", lw=1., zorder=0)
    ax.set(xlabel=title(r"Absolute frequency offset $|\delta|$ ($10^{-3}$)", "频偏绝对值 |δ|（千分之一）"),
           ylabel=title(r"Weak material coordinate $t/a$", "弱方向材料坐标 t/a"), title=title("b   Minimum and saddle collide", "b   错误极小与鞍点相撞"))
    ax.legend(handles=[Line2D([0],[0],color="0.38",label=title("Wrong minimum", "错误极小")),
                       Line2D([0],[0],color="0.38",ls="--",label=title("Saddle", "鞍点")),
                       Line2D([0],[0],color="0.38",ls=":",label=title("Jet prediction", "局部展开预测"))],
              fontsize=9, frameon=False, loc="center left")

    ax = axes[1,0]
    labels = ["R11 .23 +", "R11 .47 -", "R11 .83 +", title("Known gain + noise", "增益已知，有噪声"),
              "R10 fresh6-02", "R10 fresh6-14"]
    display = [*cases[:4],*(case for case in natural if case.get("exact"))]
    for i, case in enumerate(display):
        for s in ("-1","1"):
            _, fold = fork_data(case,s)
            x = abs(fold["delta"]); pred = abs(case["predicted"][s])
            y = i+.10*int(s)
            ax.plot([pred,x],[y,y],color=COLORS[s],lw=1.4)
            ax.plot(pred,y,marker="o",mec=COLORS[s],mfc="white",ms=5)
            ax.plot(x,y,marker="X",mec=COLORS[s],mfc=COLORS[s],ms=6)
    ax.set(yticks=np.arange(len(labels)), yticklabels=labels,
           xscale="log", xlabel=title(r"Absolute fold offset $|\delta_*|$", "碰撞频偏绝对值 |δ*|"),
           title=title("c   Candidate jet versus numerical fold", "c   局部展开预测与数值碰撞"))
    ax.invert_yaxis()
    ax.legend(handles=[Line2D([0],[0],color="0.38",mfc="white",lw=0,marker="o",label=title("Prediction", "预测")),
                       Line2D([0],[0],color="0.38",lw=0,marker="X",label=title("Numerical fold", "数值碰撞"))],
              fontsize=9, frameon=False, loc="lower left")

    ax = axes[1,1]
    short = {"direct":"Direct TRF", "pair_multistart":"Old pair",
             "generic_multistart":"Generic multi-start", "same_read_continuation":"Weight continuation",
             "cumulative":"Cumulative frequencies"}
    if zh:
        short = dict(zip(METHODS, ["直接 TRF", "旧候选双起点", "普通双起点", "同读数增权延拓", "累计频率"]))
    for i, method in enumerate(METHODS):
        for arm in ref["comparisons"][2:]:
            s = str(int(np.sign(arm["delta"])))
            ax.scatter(max(arm[method]["material_RMS"], 1e-12), i+.12*int(s), c=COLORS[s], s=28)
    ax.set(yticks=np.arange(len(METHODS)), yticklabels=[short[m] for m in METHODS],
           xscale="log", xlabel=title("Contrast RMS error", "材料对比度均方根误差"),
           title=title("d   Same offset magnitude, different basin", "d   等幅频偏下的材料恢复"))
    ax.invert_yaxis()
    ax.set_xlim(4e-13,2)
    ax.axvline(1e-12,color=".65",lw=.8,ls=":")
    ax.set_xticks([1e-12,1e-9,1e-6,1e-3,1.])

    for ax in axes.flat:
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(alpha=.18,axis="x")
        panel_title = ax.get_title()
        ax.set_title("")
        ax.set_title(panel_title, loc="left", fontweight="bold")
    output = ROOT/"figures"/("freq_signed_folds_zh" if zh else "freq_signed_folds")
    fig.savefig(output.with_suffix(".pdf"),bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"),bbox_inches="tight",dpi=235)


def summarized(cases, natural):
    model = FrequencyTM(json.loads(DATA.read_text())["geometry"])
    rows = []
    for case in [*cases,*natural]:
        if "skipped" in case:
            rows.append({"tag":case["tag"],"skip":case["skipped"]})
            continue
        row = {"tag":case["tag"],"candidate_RMS":case["candidate_material_RMS"],
               "selected_signed_offset":case["selected_offset"], "fold_center_found":case["estimated_fold"],
               "predicted":case["predicted"],"conditional_swapped_hypothesis":case["predicted_reversed"],
               "baseline_at_test_offset":{}}
        if case["exact"]:
            row["actual_fold"] = {}
            for s in ("-1", "1"):
                f = case["exact"][s]["fold"]
                if f is None:
                    continue
                estimate = case["predicted"][s]
                entry = dict(offset=f["delta"],
                             relative_prediction_error=abs((estimate-f["delta"])/f["delta"]) if estimate else None,
                             material_RMS=f["material_RMS"], min_H=f["eig"][0], next_H=f["eig"][1],
                             gradient=f["gradient"], profile_third=f["profile_third"],
                             profile_gradient_delta=f["profile_gradient_delta"])
                if not case["natural"]:
                    c = (ROBUST_SIGMA**2/GAIN_SD**2, 1.) if case["noisy"] == "profiled" else None
                    entry["raw_gap"], entry["joint_calibrated_gap"] = visible_gap(
                        model,np.array(case["candidates"]),f["delta"],c)
                row["actual_fold"][s] = entry
        for arm in case["comparisons"]:
            row["baseline_at_test_offset"][str(arm["delta"])] = {
                method:{"material_RMS":arm[method]["material_RMS"],"factorizations":arm[method]["factorizations"]}
                for method in METHODS}
        rows.append(row)
    (ROOT/"results"/"freq_summary.json").write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')


if __name__ == "__main__":
    cases = json.loads((ROOT/"results"/"freq_signed.json").read_text())
    natural = json.loads((ROOT/"results"/"freq_natural.json").read_text())
    paper_figure(cases, natural)
    preferred = ["Noto Sans CJK SC", "Noto Sans SC", "PingFang SC", "STHeiti", "Hiragino Sans GB"]
    available = set(font_manager.get_font_names())
    chinese = next((f for f in preferred if f in available), None)
    if chinese:
        with plt.rc_context({"font.family":["DejaVu Sans", chinese], "mathtext.fontset":"dejavusans", "axes.unicode_minus":False}):
            paper_figure(cases,natural,zh=True)
    else:
        print("Install a Chinese font (e.g. Noto Sans CJK SC) to render the Chinese variant.")
    summarized(cases, natural)
