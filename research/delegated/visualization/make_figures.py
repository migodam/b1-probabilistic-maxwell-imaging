#!/usr/bin/env python
"""B1 visualization: high-quality static scientific figures + a self-contained
local HTML Top-K dashboard.

READ-ONLY over results/. This script:
  * reads only the declared saved derived arrays/JSON (no physics, no budget),
  * rebuilds per-topology representative material images from the saved
    conditional-mean material parameters + the frozen protocol geometry,
    and verifies that the event-weighted mixture reproduces the saved
    ``posterior_mean`` / ``void_probability`` arrays exactly (audit check),
  * writes figures to ``figures/``, the dashboard to ``dashboard.html``, and
    its own script / audit / receipt under ``research/delegated/visualization``.

It never runs a Maxwell forward, never trains, never consumes budget, never
touches locked/historical arrays, and never edits src/config/results.

Interpreter (explicit full path):
  /Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python
"""
from __future__ import annotations

import base64
import io
import json
import logging
import os
import platform
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

WALL_T0 = time.perf_counter()
CPU_T0 = time.process_time()

ROOT = Path(__file__).resolve().parents[3]
FIG_DIR = ROOT / "figures"
DELEG = ROOT / "research" / "delegated" / "visualization"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------------
# Style / fonts
# ----------------------------------------------------------------------------
_avail = {f.name for f in font_manager.fontManager.ttflist}
_cjk = [n for n in ["PingFang SC", "Hiragino Sans GB", "Heiti SC", "Songti SC",
                    "Arial Unicode MS"] if n in _avail]
plt.rcParams.update({
    "font.family": _cjk + ["DejaVu Sans"],
    "axes.unicode_minus": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.titlesize": 11.5,
    "axes.labelsize": 10.5,
    "axes.titleweight": "bold",
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
})

TOP_NAME = {0: "实心 (无空腔)", 1: "中心空腔", 2: "中心+环状空腔"}
TOP_SHORT = {0: "实心", 1: "中心空腔", 2: "双空腔"}
TOP_COLORS = {0: "#496b9b", 1: "#147d92", 2: "#db9245"}
CASE_ORDER = ["competitive", "distinguishable", "information_limited",
              "gain_mismatch", "gain_marginalized", "representation_mismatch"]
CASE_ZH = {
    "competitive": "竞争情形",
    "distinguishable": "可分辨",
    "information_limited": "信息受限",
    "gain_mismatch": "增益失配",
    "gain_marginalized": "增益边缘化",
    "representation_mismatch": "表示失配 (负对照)",
}


def base_name(name: str) -> str:
    return name.split("refined_", 1)[1] if name.startswith("refined_") else name


def case_label(name: str) -> str:
    if name in CASE_ZH:
        return CASE_ZH[name]
    return CASE_ZH.get(base_name(name), base_name(name)) + "（精化）"

# fixed, comparable color scales (real background = 1, imag background = 0)
REAL_VMIN, REAL_VMAX = 1.0, 2.8
IMAG_VMIN, IMAG_VMAX = 0.0, 0.16


def zh_ok() -> bool:
    return bool(_cjk)


# ----------------------------------------------------------------------------
# Loaders (read-only)
# ----------------------------------------------------------------------------
def load_json(path: Path):
    return json.loads(Path(path).read_text())


def material_image(H, theta, config, grid=151):
    """Exact replica of src/b1/inference.material_image (deterministic, no physics)."""
    z = np.linspace(-1.1, 1.1, grid)
    X, Y = np.meshgrid(z, z)
    r = np.hypot(X, Y)
    inside = r <= config["outer_radius"]
    void = np.zeros_like(inside)
    for a, b in config["topologies"][int(H)]["void_intervals"]:
        void |= (r >= a) & (r <= b) & inside
    eps = np.ones(r.shape, complex)
    eps[inside & ~void] = complex(*theta)
    return z, eps, void


def gather():
    config = load_json(ROOT / "configs" / "protocol.json")
    summary = load_json(ROOT / "results" / "PRIMARY_SUMMARY.json")
    mcm = load_json(ROOT / "results" / "MODEL_CORRECTION_MECHANISM.json")
    prior_sens = load_json(ROOT / "results" / "PRIOR_SENSITIVITY.json")
    importance = load_json(ROOT / "results" / "importance" / "RESULT.json")
    calib = load_json(ROOT / "results" / "CALIBRATION_RESULTS.json")
    cases = {}
    for pj_path in sorted((ROOT / "results" / "cases").glob("*/probabilities.json")):
        name = pj_path.parent.name
        pj = load_json(pj_path)
        mf = np.load(pj_path.parent / "material_fields.npz")
        cases[name] = {
            "probabilities": pj,
            "axis": mf["axis"],
            "posterior_mean": mf["posterior_mean"],
            "material_complex_sd": mf["material_complex_sd"],
            "void_probability": mf["void_probability"],
        }
    return config, summary, mcm, prior_sens, importance, calib, cases


def per_topology_images(case_pj, config):
    """Rebuild the three per-topology representative material fields (region means)."""
    ims, voids = [], []
    for h in range(3):
        mean = case_pj["result"]["conditional_material"][h]["mean"]
        _, im, void = material_image(h, mean, config)
        ims.append(im)
        voids.append(void)
    return np.asarray(ims), np.asarray(voids)


def validate(config, cases):
    """Audit: event-weighted mixture of rebuilt representatives must equal saved arrays."""
    checks = []
    for name, c in cases.items():
        p = np.asarray(c["probabilities"]["result"]["topology_probabilities"])
        ims, voids = per_topology_images(c["probabilities"], config)
        mix = np.sum(p[:, None, None] * ims, axis=0)
        vp = np.sum(p[:, None, None] * voids, axis=0)
        e1 = float(np.max(np.abs(mix - c["posterior_mean"])))
        e2 = float(np.max(np.abs(vp - c["void_probability"])))
        checks.append(dict(case=name, mixture_max_abs_err=e1, void_max_abs_err=e2,
                           ok=bool(max(e1, e2) < 1e-9)))
    return checks


# ----------------------------------------------------------------------------
# Figure helpers
# ----------------------------------------------------------------------------
def fixed_geometry_circles(ax, config, H=None, lw=1.0):
    ax.add_patch(plt.Circle((0, 0), config["outer_radius"], fill=False,
                            ec="white", lw=lw, ls="-", alpha=0.9))
    ax.add_patch(plt.Circle((0, 0), config["outer_radius"], fill=False,
                            ec="black", lw=lw, ls=":", alpha=0.9))
    if H is not None:
        for a, b in config["topologies"][int(H)]["void_intervals"]:
            if a > 0:
                ax.add_patch(plt.Circle((0, 0), a, fill=False, ec="cyan",
                                        lw=lw, ls="--", alpha=0.9))
            if b < config["outer_radius"]:
                ax.add_patch(plt.Circle((0, 0), b, fill=False, ec="cyan",
                                        lw=lw, ls="--", alpha=0.9))


def fig_to_b64(fig, dpi):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def save_high(fig, name):
    fig.savefig(FIG_DIR / name, dpi=200, bbox_inches="tight", facecolor="white")


def p_label(p, digits=2):
    """Percent formatter without rounding tiny nonzero masses to zero."""
    return f"{100*p:.{digits}f}%" if p >= 0.00005 else f"{100*p:.2e}%"


def coverage_marginal(calib):
    # The parent has corrected the key: these intervals use H-marginal weights.
    # Legacy fallback lets the renderer run while the JSON rename is in progress.
    return calib.get("marginal_material_90_coverage", calib.get("conditional_material_90_coverage"))


# ----------------------------------------------------------------------------
# Figure 1: competitive representative material fields (3 topologies + mixture)
# ----------------------------------------------------------------------------
def fig_representatives(config, case_pj, case_name, dpi=200, high=True,
                        selected=None, include_mix=True):
    ims, voids = per_topology_images(case_pj, config)
    p = np.asarray(case_pj["result"]["topology_probabilities"])
    mix = np.sum(p[:, None, None] * ims, axis=0)
    z = np.linspace(-1.1, 1.1, ims.shape[1])
    ext = [z[0], z[-1], z[0], z[-1]]
    cols = (list(range(3)) if selected is None else list(selected)) + (["mix"] if include_mix else [])
    nc = len(cols)
    fig = plt.figure(figsize=(3.7*nc+1.0, 8.7))
    gs = fig.add_gridspec(2, nc+1, width_ratios=[1]*nc + [0.065], wspace=0.12,
                          hspace=0.19, left=0.055, right=0.935, top=0.80, bottom=0.095)
    im_real = im_imag = None
    for row, (q, cmap, vmin, vmax) in enumerate(
            [("real", "viridis", REAL_VMIN, REAL_VMAX),
             ("imag", "magma", IMAG_VMIN, IMAG_VMAX)]):
        for col, key in enumerate(cols):
            ax = fig.add_subplot(gs[row, col])
            if key == "mix":
                field = mix
                ttl = ("事件边际均值\n可不对应任一模板" if row == 0 else "")
            else:
                field = ims[key]
                ttl = (f"H{key} · {TOP_NAME[key]}\n拓扑事件质量 {p_label(p[key], 4)}" if row == 0 else "")
            data = field.real if q == "real" else field.imag
            m = ax.imshow(data, extent=ext, origin="lower", cmap=cmap,
                          vmin=vmin, vmax=vmax)
            if q == "real":
                im_real = m
            else:
                im_imag = m
            fixed_geometry_circles(ax, config, None if key == "mix" else key, lw=0.9)
            ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
            ax.tick_params(labelsize=8)
            ax.grid(False)
            if row == 0:
                ax.set_title(ttl, fontsize=10.6, pad=10)
            if col == 0:
                ax.set_ylabel(("实部 Re ε" if q == "real" else "虚部 Im ε"),
                              fontsize=11)
            if row == 1:
                ax.set_xlabel("x / R", fontsize=9)
    cax1 = fig.add_subplot(gs[0, nc]); fig.colorbar(im_real, cax=cax1)
    cax1.set_title("Re ε", fontsize=9)
    cax2 = fig.add_subplot(gs[1, nc]); fig.colorbar(im_imag, cax=cax2)
    cax2.set_title("Im ε", fontsize=9)
    display_name = "三个拓扑" if nc > 2 else "Top-2 拓扑"
    fig.suptitle(f"B1 · {case_label(case_name)}：{display_name}的代表材料场",
                 fontsize=14, fontweight="bold", y=0.975)
    if nc > 2:
        subtitle = ("三维同心球模型的二维截面；各分支使用条件均值材料参数。标注为拓扑事件质量，连续先验下单幅精确图像概率为 0。\n"
                    "上排 Re ε 描述材料响应，下排 Im ε 描述损耗；各排共用色标。背景 Re ε=1、Im ε=0；青色虚线为固定空腔边界。")
    else:
        subtitle = ("按拓扑事件质量排序；两分支保持原始概率，不重新归一化。\n"
                    "三维同心球的二维截面；条件均值代表场，背景 Re ε=1、Im ε=0。")
    fig.text(0.055, 0.922, subtitle, fontsize=9.1, va="top", linespacing=1.7)
    fig.text(0.055, 0.028,
             "圆心与边界固定；概率条件于三个声明模板。更广泛物理几何后验未量化；代表图不等于唯一恢复图像。",
             fontsize=8.8, color="#8a4b00")
    if high:
        save_high(fig, f"fig01_{case_name}_representatives.png")
    return fig


# ----------------------------------------------------------------------------
# Figure 2: Top-2 + remaining modeled mass
# ----------------------------------------------------------------------------
def fig_top2(summary, prior_strip=True, dpi=200):
    p = np.asarray(summary["cases"]["competitive"]["result"]["topology_probabilities"])
    order = np.argsort(-p)
    fig = plt.figure(figsize=(15.0, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.65], wspace=0.17,
                          left=0.065, right=0.98, top=0.78, bottom=0.17)
    ax = fig.add_subplot(gs[0, 0])
    cols = [TOP_COLORS[int(h)] for h in order]
    ax.bar(range(3), p[order], color=cols, edgecolor="#344254", linewidth=0.7, width=0.64)
    for i, h in enumerate(order):
        ax.text(i, p[h]+0.014, p_label(p[h],4), ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"Top-1 · H{order[0]}\n{TOP_SHORT[int(order[0])]}",
                        f"Top-2 · H{order[1]}\n{TOP_SHORT[int(order[1])]}",
                        f"剩余 · H{order[2]}\n{TOP_SHORT[int(order[2])]}"])
    ax.set_ylim(0,0.62); ax.set_ylabel("原始拓扑后验质量")
    ax.set_title("竞争情形：Top-2 与保留质量",fontsize=11)
    ax.grid(axis="x",visible=False)
    ax.text(0.5,-0.27,"Top-2 总质量 "+p_label(p[order[:2]].sum(),4)+"；未重新归一化",
            transform=ax.transAxes,ha="center",fontsize=9,color="#344254")

    # An exact table avoids labels overflowing almost-zero stacked segments.
    ar = fig.add_subplot(gs[0,1]); ar.axis("off")
    ar.set_xlim(0,1); ar.set_ylim(0,1)
    ar.set_title("六个主情形的排序与原始 5/7 参考状态",fontsize=11,pad=12)
    widths=[0.25,0.19,0.19,0.19,0.18]
    headers=["主情形","Top-1","Top-2","保留第三","原始 5/7"]
    x=0
    for w,hdr in zip(widths,headers):
        _grid_cell(ar,x,0.84,w,0.12,"#e8eef6",hdr,fs=9,bold=True);x+=w
    names=[n for n in CASE_ORDER if n in summary["cases"]]
    for i,n in enumerate(names):
        pp=np.asarray(summary["cases"][n]["result"]["topology_probabilities"])
        o=np.argsort(-pp); st=summary["cases"][n]["refinement_5_vs_7"]["status"]
        vals=[CASE_ZH[n]]+[f"H{int(h)}\n{p_label(pp[h])}" for h in o]+[
            "数值 PASS" if st=="NUMERICAL_REFERENCE_PASS" else "UNRESOLVED"]
        x=0;y=0.84-(i+1)*0.12
        for j,(w,t) in enumerate(zip(widths,vals)):
            bg="#f4f7fb" if i%2==0 else "white"
            if j==4: bg=_status_color(st)
            _grid_cell(ar,x,y,w,0.12,bg,t,fs=8.8,bold=j==4,
                       color="#9c3020" if j==4 and "UNRESOLVED" in t else "#16202b");x+=w
    fig.suptitle("B1 · Top-2 概率与剩余建模质量",fontsize=14,fontweight="bold",y=0.975)
    fig.text(0.065,0.91,
             "每个百分比是整个拓扑分支的后验质量；Top-2 展示保留第三分支，不把两项重新缩放到 100%。\n"
             "声明先验在模板外为 0；更广泛物理几何后验未定义、未量化。数值 PASS 仅指细化诊断，均无严格认证。",
             fontsize=9.2,va="top",linespacing=1.6)
    save_high(fig,"fig02_top2_remaining.png")
    return fig


# ----------------------------------------------------------------------------
# Figure 3: cavity occupancy + complex material uncertainty (event-marginal)
# ----------------------------------------------------------------------------
def fig_cavity_uncertainty(config, case, case_name, dpi=200, high=True):
    z=case["axis"]; ext=[z[0],z[-1],z[0],z[-1]]
    vp=case["void_probability"]; sd=case["material_complex_sd"]
    fig=plt.figure(figsize=(13.6,6.9))
    gs=fig.add_gridspec(1,2,wspace=0.26,left=0.065,right=0.945,top=0.78,bottom=0.17)
    ax0=fig.add_subplot(gs[0,0])
    m0=ax0.imshow(vp,extent=ext,origin="lower",cmap="cividis",vmin=0,vmax=1)
    fixed_geometry_circles(ax0,config,None)
    for h,color in [(1,"#00dfe8"),(2,"#ef6f58")]:
        for a,b in config["topologies"][h]["void_intervals"]:
            for rr in (a,b):
                if rr>0 and rr<config["outer_radius"]:
                    ax0.add_patch(plt.Circle((0,0),rr,fill=False,ec=color,lw=0.9,ls="--"))
    ax0.set_title("空腔占据概率\n"+r"$P(\mathrm{void}\mid y)=\sum_H p_H\,v_H$",fontsize=11,pad=10)
    cb0=fig.colorbar(m0,ax=ax0,fraction=0.048,pad=0.035);cb0.set_label("概率")
    ax1=fig.add_subplot(gs[0,1])
    m1=ax1.imshow(sd,extent=ext,origin="lower",cmap="inferno",vmin=0,vmax=float(sd.max()))
    fixed_geometry_circles(ax1,config,None)
    ax1.set_title("复介电常数的边际标准差\n"+r"$\sigma_\varepsilon=\sqrt{\mathbb{E}|\varepsilon-\mathbb{E}\varepsilon|^2}$",fontsize=11,pad=10)
    cb1=fig.colorbar(m1,ax=ax1,fraction=0.048,pad=0.035);cb1.set_label(r"$\sigma_\varepsilon$ (complex)")
    for ax in [ax0,ax1]:
        ax.grid(False);ax.set_xlabel("x / R");ax.set_ylabel("y / R")
        ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1])
    fig.suptitle(f"B1 · {case_label(case_name)}：空腔占据与材料不确定度",fontsize=14,fontweight="bold",y=0.975)
    fig.text(0.065,0.915,
             "左图把三个固定空腔模板按拓扑概率相加；右图同时包含拓扑差异和拓扑内部材料参数方差。\n"
             "两图均为三维同心球的二维截面。概率与方差都依赖所选模板和先验。",
             fontsize=9.2,va="top",linespacing=1.6)
    fig.text(0.065,0.062,
             "青虚线：H1 半径 0.5；红虚线：H2 半径 0.24、0.46、0.6；外边界 R=1。背景真空 Re ε=1、Im ε=0。\n"
             "圆心固定在原点，未进行定位推断；更广泛物理几何后验未定义。环状空腔不具有唯一的点中心。",
             fontsize=9.0,va="top",linespacing=1.6,color="#8a4b00")
    if high: save_high(fig,f"fig03_{case_name}_cavity_uncertainty.png")
    return fig


# ----------------------------------------------------------------------------
# Figure 4: cases probabilities / prior / negative controls
# ----------------------------------------------------------------------------
def fig_cases_prior_negative(summary, prior_sens):
    fig=plt.figure(figsize=(15.2,10.0))
    gs=fig.add_gridspec(2,2,height_ratios=[1.1,1],hspace=0.48,wspace=0.25,
                         left=0.065,right=0.97,top=0.78,bottom=0.15)
    colors=[TOP_COLORS[h] for h in range(3)]
    ax=fig.add_subplot(gs[0,:]);names=[n for n in CASE_ORDER if n in summary["cases"]]
    x=np.arange(len(names));w=0.24
    for h in range(3):
        vals=[summary["cases"][n]["result"]["topology_probabilities"][h] for n in names]
        ax.bar(x+(h-1)*w,vals,w,color=colors[h],label=f"H{h} · {TOP_NAME[h]}",edgecolor="#344254",linewidth=0.5)
    for i,n in enumerate(names):
        t=summary["cases"][n].get("truth",{})
        if n=="gain_marginalized":t=summary["cases"]["gain_mismatch"].get("truth",{})
        truth=f"真值 H{t.get('H','—')}"+(" · 支持外" if t.get("within_inference_support") is False else "")
        st=summary["cases"][n]["refinement_5_vs_7"]["status"]
        ax.text(i,1.10,truth,ha="center",fontsize=9,transform=ax.get_xaxis_transform(),clip_on=False)
        ax.text(i,1.025,"5/7 数值 PASS" if st=="NUMERICAL_REFERENCE_PASS" else "5/7 UNRESOLVED",
                ha="center",fontsize=8.5,color="#216a37" if st=="NUMERICAL_REFERENCE_PASS" else "#a13c26",
                transform=ax.get_xaxis_transform(),clip_on=False)
    ax.set_xticks(x);ax.set_xticklabels([CASE_ZH[n] for n in names]);ax.set_ylim(0,1.0)
    ax.set_ylabel("拓扑后验质量");ax.grid(axis="x",visible=False)
    fig.legend(*ax.get_legend_handles_labels(),ncol=3,loc="upper center",bbox_to_anchor=(0.52,0.875),frameon=False)

    ab=fig.add_subplot(gs[1,0])
    sets=[("竞争\n均匀先验",summary["cases"]["competitive"]["result"]["topology_probabilities"]),
          ("竞争\nBeta 材料先验",prior_sens["material_beta_prior"]["topology"]),
          ("信息受限\n均匀先验",summary["cases"]["information_limited"]["result"]["topology_probabilities"]),
          ("信息受限\n实心偏好",prior_sens["solid_favored"]["topology"]),
          ("信息受限\n双空腔偏好",prior_sens["double_void_favored"]["topology"])]
    xx=np.arange(len(sets))
    for h in range(3):
        ab.bar(xx+(h-1)*w,[v[1][h] for v in sets],w,color=colors[h],edgecolor="#344254",linewidth=0.5)
    ab.set_xticks(xx);ab.set_xticklabels([v[0] for v in sets],fontsize=9)
    ab.set_ylim(0,1);ab.set_ylabel("拓扑后验质量");ab.grid(axis="x",visible=False)
    ab.set_title("先验敏感性：材料先验与拓扑先验",fontsize=11,pad=12)

    ac=fig.add_subplot(gs[1,1]);rm=summary["cases"]["representation_mismatch"]
    pp=rm["result"]["topology_probabilities"]
    ac.bar(range(3),pp,color=colors,edgecolor="#344254",linewidth=0.5,width=0.58)
    for i,v in enumerate(pp):
        ac.text(i,v+0.025,p_label(v),ha="center",fontsize=9.5,fontweight="bold")
    ac.set_xticks(range(3));ac.set_xticklabels([f"H{h} · {TOP_SHORT[h]}" for h in range(3)])
    ac.set_ylim(0,1);ac.set_ylabel("拓扑后验质量");ac.grid(axis="x",visible=False)
    ac.set_title("负对照：真实几何不在三个模板内",fontsize=11,pad=12)
    fig.suptitle("B1 · 情形概率、先验敏感性与负对照",fontsize=14,fontweight="bold",y=0.975)
    fig.text(0.065,0.915,"主情形采用原始 order7 全网格结果；原始 5/7 参考失败保留。所有概率均条件于三个声明模板。",
             fontsize=9.2,va="top")
    fig.text(0.065,0.077,
             "材料 Beta 先验仅在 order7 评估，未独立细化；低信息条件下，更改拓扑先验可显著改变排序。\n"
             f"负对照的预测尾部概率为 {rm['result']['posterior_predictive_tail_probability']:.4g}；模型仍被迫把全部质量分配给模板内分支。\n"
             "声明先验在三个模板外为 0；更广泛物理几何的后验质量未定义、未量化，不能据此断言模板外几何不可能。",
             fontsize=9.0,va="top",linespacing=1.6,color="#8a4b00")
    save_high(fig,"fig04_cases_prior_negative.png");return fig


# ----------------------------------------------------------------------------
# Figure 5: calibration diagnostics (24 parent clusters / 720 noises)
# ----------------------------------------------------------------------------
def fig_calibration(calib, review=None):
    fig=plt.figure(figsize=(15.4,8.3))
    gs=fig.add_gridspec(1,2,width_ratios=[1,1.12],wspace=0.38,left=0.07,right=0.97,top=0.78,bottom=0.28)
    rv=review['calibration'] if review else None
    bidx={(b['H'],round(b['low'],6)):b for b in rv['calibration_bins_with_cluster_uncertainty']} if rv else {}
    ax=fig.add_subplot(gs[0,0]);ax.plot([0,1],[0,1],color="#888",ls="--",lw=1,label="理想校准")
    for h in range(3):
        bb=[b for b in calib["calibration_bins"] if b["H"]==h and b["count"]>0]
        mp=np.array([b["mean_probability"] for b in bb]);of=np.array([b["observed_frequency"] for b in bb]);cnt=np.array([b["count"] for b in bb])
        edges=[];lws=[]
        for b,xx,yy in zip(bb,mp,of):
            r=bidx.get((h,round(b['low'],6)))
            flag=False
            if r:
                lo,hi=r['observed_frequency_cluster95']
                ax.errorbar(xx,yy,yerr=[[max(0,yy-lo)],[max(0,hi-yy)]],fmt='none',ecolor=TOP_COLORS[h],capsize=3,lw=.8,alpha=.6)
                d0,d1=r['observed_minus_reported_cluster95'];flag=(d0>0 or d1<0)
            edges.append('#a13c26' if flag else '#263646');lws.append(1.8 if flag else .6)
        ax.scatter(mp,of,s=20+220*cnt/cnt.max(),color=TOP_COLORS[h],edgecolor=edges,linewidth=lws,
                   zorder=3,label=f"H{h} · {TOP_SHORT[h]}")
        ax.plot(mp,of,color=TOP_COLORS[h],lw=1,alpha=0.65)
        for xx,yy,n in zip(mp,of,cnt):
            if n<10:ax.annotate(f"n={int(n)}",(xx,yy),xytext=(8,-8),textcoords="offset points",fontsize=8.5,color="#8a4b00")
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_xlabel("报告概率（分箱均值）");ax.set_ylabel("事件观测频率")
    ax.set_title("可靠性图：母体聚类 95% 区间",fontsize=11,pad=12)
    ax.legend(loc="upper left",frameon=True,fontsize=9)
    ar=fig.add_subplot(gs[0,1])
    cv=coverage_marginal(calib)
    items=[("90% 拓扑可信集合覆盖",calib["topology_credible_90_coverage"],.90),
           ("Top-2 经验覆盖",calib["top2_empirical_coverage"],None),
           ("Top-2 平均报告质量",calib["top2_mean_reported_mass"],None),
           ("材料边际 90% 覆盖 Re ε",cv[0],.90),
           ("材料边际 90% 覆盖 Im ε",cv[1],.90)]
    metric_keys=['topology_credible90_coverage','top2_coverage','top2_reported_mass','marginal_real90_coverage','marginal_imag90_coverage']
    if rv:
        items += [('已知真 H 条件 90% 覆盖 Re ε',rv['metrics']['truth_H_conditional_real90_coverage']['estimate'],.90),
                  ('已知真 H 条件 90% 覆盖 Im ε',rv['metrics']['truth_H_conditional_imag90_coverage']['estimate'],.90)]
        metric_keys += ['truth_H_conditional_real90_coverage','truth_H_conditional_imag90_coverage']
    yy=np.arange(len(items));ar.barh(yy,[v[1] for v in items],height=.56,color="#d3e3ef",edgecolor="#506377",linewidth=.5)
    for i,(label,v,nom) in enumerate(items):
        if rv:
            lo,hi=rv['metrics'][metric_keys[i]]['parent_stratified_bootstrap95']
            ar.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='none',ecolor='#637784',capsize=3,lw=1,zorder=3)
        if nom is not None:ar.plot([nom,nom],[i-.34,i+.34],color="#b74432",lw=2)
        ar.text(1.025,i,p_label(v),va="center",fontsize=9.4,fontweight="bold")
    ar.set_yticks(yy);ar.set_yticklabels([v[0] for v in items],fontsize=9.3)
    ar.set_xlim(0,1.18);ar.invert_yaxis();ar.grid(axis="y",visible=False)
    ar.set_title("覆盖率：母体聚类 95% 区间；红短线=90%",fontsize=11,pad=12)
    ar.set_xlabel("覆盖比例 / 报告质量")
    fig.suptitle("B1 · 探索性校准诊断：24 个独立母体 / 720 次噪声",fontsize=14,fontweight="bold",y=.975)
    fig.text(.07,.915,
             "720 次噪声来自 24 个连续材料母体，不能视为 720 个独立物体。材料区间使用 H 边际权重。\n"
             "聚类自助区间不包含求解器偏差或表示偏差；这是声明模型内的开发诊断。",
             fontsize=9.2,va="top",linespacing=1.6,color="#8a4b00")
    bci=calib['Brier_cluster_bootstrap_95'];lci=calib['log_score_cluster_bootstrap_95']
    fig.text(.07,.18,f"多类 Brier = {calib['Brier_multiclass']:.3f}\n母体聚类自助 95% 区间 [{bci[0]:.3f}, {bci[1]:.3f}]",
             fontsize=10,va="top",linespacing=1.6)
    fig.text(.50,.18,f"对数损失 = {calib['log_score']:.3f}\n母体聚类自助 95% 区间 [{lci[0]:.3f}, {lci[1]:.3f}]",
             fontsize=10,va="top",linespacing=1.6)
    fig.text(.07,.082,"点面积=分箱噪声数；红外圈=观测−报告差值 95% 区间不含 0（描述性，未校正多重比较）。\n"
             "已知真 H 的条件覆盖仅是 oracle 描述性对照。此母体集细化通过不构成校准通过或严格认证。",
             fontsize=9.0,color="#8a4b00",va='top',linespacing=1.6)
    save_high(fig,"fig05_calibration.png");return fig


# ----------------------------------------------------------------------------
# Figure 6: reference status
# ----------------------------------------------------------------------------
def _status_color(s):
    return "#d9f0dd" if s == "NUMERICAL_REFERENCE_PASS" else "#fde2c8"


def _grid_cell(ax, x0, y0, w, h, fc, txt, fs=8.8, bold=False, color="#16202b"):
    ax.add_patch(Rectangle((x0, y0), w, h, facecolor=fc, edgecolor="#8a97a6", lw=0.6))
    ax.text(x0 + w / 2, y0 + h / 2, txt, ha="center", va="center",
            fontsize=fs, fontweight=("bold" if bold else "normal"), color=color)


def fig_reference(summary, importance):
    fig=plt.figure(figsize=(15.4,9.0))
    gs=fig.add_gridspec(2,2,height_ratios=[1.35,1],hspace=.35,wspace=.29,
                         left=.055,right=.975,top=.83,bottom=.15)
    ax=fig.add_subplot(gs[0,:]);ax.axis("off");ax.set_xlim(0,1);ax.set_ylim(0,1)
    names=[n for n in CASE_ORDER if n in summary["cases"]]
    headers=["主情形","3/5 Δp","3/5 Δlog BF","5/7 Δp","5/7 Δlog BF","5/7 状态","3/5 状态"]
    widths=[.17,.13,.145,.13,.145,.14,.14]
    x=0
    for w,h in zip(widths,headers):_grid_cell(ax,x,.87,w,.11,"#e8eef6",h,fs=9.4,bold=True);x+=w
    for i,n in enumerate(names):
        ss=summary["cases"][n];r5=ss["refinement_5_vs_7"];r3=ss.get("refinement_3_vs_5")
        row=[CASE_ZH[n],f"{r3['max_probability_difference']:.2e}" if r3 else "未计算",
             f"{r3['max_log_BF_difference']:.2e}" if r3 else "未计算",
             f"{r5['max_probability_difference']:.2e}",f"{r5['max_log_BF_difference']:.2e}",
             "数值 PASS" if r5['status']=="NUMERICAL_REFERENCE_PASS" else "REFERENCE\nUNRESOLVED",
             "未计算" if r3 is None else ("数值 PASS" if r3['status']=="NUMERICAL_REFERENCE_PASS" else "REFERENCE\nUNRESOLVED")]
        x=0;y=.87-(i+1)*.125
        for j,(w,t) in enumerate(zip(widths,row)):
            ref=r3 if j in [1,2,6] else r5
            fc="#f0f4f8" if j==0 else (_status_color(ref['status']) if ref else "#efefef")
            _grid_cell(ax,x,y,w,.125,fc,t,fs=9.1,bold=j>=5,
                       color="#96301d" if "UNRESOLVED" in t else "#16202b");x+=w
    ap=fig.add_subplot(gs[1,0]);p0=np.asarray(summary['cases']['competitive']['result']['topology_probabilities'])
    pi=np.asarray(importance['topology_probabilities']);ci=np.asarray(importance['bootstrap_95_intervals']);xx=np.arange(3)
    ap.errorbar(xx-.055,pi,yerr=np.stack([pi-ci[:,0],ci[:,1]-pi]),fmt='o',color='#a76e17',capsize=5,
                label="重要性采样 + 自助 95% 区间")
    ap.scatter(xx+.055,p0,marker='D',s=40,color='#263e63',label="原始 order7 全网格",zorder=3)
    ap.set_xticks(xx);ap.set_xticklabels([f"H{h} · {TOP_SHORT[h]}" for h in range(3)])
    ap.set_ylim(0,.63);ap.set_ylabel("拓扑后验质量");ap.set_title("竞争情形：两种数值近似",fontsize=11)
    ap.legend(loc='upper left',fontsize=9);ap.grid(axis='x',visible=False)
    ae=fig.add_subplot(gs[1,1]);err=importance['probability_error_vs_order7'];tol=.02
    ae.barh([0],[err],height=.44,color='#c89241',edgecolor='#48525e')
    ae.axvline(tol,color='#b74432',ls='--',lw=1.5)
    ae.text(tol,.36,"容差 0.02",ha='center',color='#a13c26',fontsize=9.5)
    ae.text(err+.0008,0,f"{err:.5f}",va='center',fontsize=11,fontweight='bold')
    ae.set_xlim(0,.038);ae.set_ylim(-.5,.5);ae.set_yticks([]);ae.grid(axis='y',visible=False)
    ae.set_xlabel("最大拓扑概率差");ae.set_title("采样误差超过容差：不认证",fontsize=11)
    fig.suptitle("B1 · 原始全网格参考状态与采样对照",fontsize=14,fontweight='bold',y=.975)
    fig.text(.055,.925,
             "参考 PASS 要求最大 Δp ≤ 0.02 且最大 Δlog BF ≤ 0.1；四个 REFERENCE_UNRESOLVED 原样保留。\n"
             "这里检查的是数值细化，所有 certified_bound 均为 false；新增经验细化结果在另一张图中独立展示。",
             fontsize=9.2,va='top',linespacing=1.6,color='#8a4b00')
    fig.text(.055,.07,
             "重要性采样仅 16 样本/分支（共 48）；有限样本归一化比值有偏，自助区间不能认证未发现的分支内模态。\n"
             "三个质量均条件于声明模板；更广泛物理几何后验未定义、未量化。",
             fontsize=9.1,va='top',linespacing=1.6,color='#8a4b00')
    save_high(fig,'fig06_reference_status.png');return fig


def fig_refined_comparison(summary, refinement):
    fig=plt.figure(figsize=(15.4,10.4))
    gs=fig.add_gridspec(2,1,height_ratios=[1.35,1.0],hspace=.31,left=.14,right=.96,top=.82,bottom=.12)
    ax=fig.add_subplot(gs[0,0]);names=[n for n in CASE_ORDER if n in refinement['cases']]
    for i,n in enumerate(names):
        yy=len(names)-1-i
        for shift,source,tag in [(+.15,summary,'原'),(-.15,refinement,'细')]:
            pp=np.asarray(source['cases'][n]['result']['topology_probabilities']);left=0
            for h,pv in enumerate(pp):
                ax.barh(yy+shift,pv,left=left,height=.23,color=TOP_COLORS[h],edgecolor='white',linewidth=.3)
                if pv>.14:ax.text(left+pv/2,yy+shift,f"H{h} · {p_label(pv)}",ha='center',va='center',color='white',fontsize=8.7,fontweight='bold')
                left+=pv
            ax.text(-.012,yy+shift,tag,ha='right',va='center',fontsize=8.7)
        change=refinement['cases'][n]['topology_change_from_original_order7']
        ax.text(1.025,yy,f"最大变化\n{100*change:.4g} 个百分点",va='center',fontsize=8.7,color='#5b4a36')
    ax.set_yticks(range(len(names)));ax.set_yticklabels([CASE_ZH[n] for n in reversed(names)],fontsize=9.5)
    ax.tick_params(axis='y',pad=26)
    ax.set_xlim(0,1.22);ax.set_xticks(np.arange(0,1.01,.2));ax.set_xlabel('模板内拓扑后验质量（原：全网格 order7；细：poly7 + order25 积分）')
    ax.grid(axis='y',visible=False);ax.axvline(1,ls=':',color='#77818d',lw=.8)
    handles=[Rectangle((0,0),1,1,color=TOP_COLORS[h]) for h in range(3)]
    ax.legend(handles,[f"H{h} · {TOP_SHORT[h]}" for h in range(3)],ncol=3,loc='lower left',bbox_to_anchor=(0,1.01),frameon=False)
    ar=fig.add_subplot(gs[1,0]);ar.axis('off');ar.set_xlim(0,1);ar.set_ylim(0,1)
    widths=[.23,.18,.20,.18,.21]
    headers=['情形','积分 17/25 Δp','积分 17/25 Δlog BF','插值 poly5/7 Δp','插值 poly5/7 Δlog BF']
    x=0
    for w,hdr in zip(widths,headers):_grid_cell(ar,x,.85,w,.11,'#e8eef6',hdr,fs=8.8,bold=True);x+=w
    for i,n in enumerate(names):
        alt=refinement['cases'][n]['alternative_reference'];integ=alt['integration_17_vs_25'];poly=alt['forward_interpolation_order5_vs7']
        vals=[CASE_ZH[n],f"{integ['max_probability_difference']:.2e}",f"{integ['max_log_BF_difference']:.2e}",
              f"{poly['max_probability_difference']:.2e}",f"{poly['max_log_BF_difference']:.2e}"]
        x=0;y=.85-(i+1)*.115
        for w,txt in zip(widths,vals):_grid_cell(ar,x,y,w,.115,'#f0f6f5' if i%2==0 else 'white',txt,fs=9.0);x+=w
    fig.suptitle('B1 · 独立展示的经验细化结果',fontsize=14,fontweight='bold',y=.975)
    fig.text(.065,.92,
             '六个替代参考状态均为 EMPIRICALLY_REFINED_REFERENCE；积分和插值细化诊断在声明容差内。\n'
             '它近似原始连续全 Mie 后验；没有统一插值余项界，也没有未见尾部误差认证。原始四个参考失败不被改写。',
             fontsize=9.2,va='top',linespacing=1.6,color='#8a4b00')
    fig.text(.065,.062,
             '圆心和边界固定；所有质量条件于三个模板，不含更广泛物理几何。增益边际化情形在细化后 Top-1 可翻转，分支仍接近各半。\n'
             '保存的 24 次全正演插值审计与 4 次独立 treams 对照另图展示；经验点检不等于严格后验认证。',
             fontsize=9.0,va='top',linespacing=1.6,color='#8a4b00')
    save_high(fig,'fig07_refined_reference_comparison.png');return fig


def fig_independent_forward(refinement, forward):
    fig=plt.figure(figsize=(14.8,7.0))
    gs=fig.add_gridspec(1,2,width_ratios=[1.5,1],wspace=.25,left=.07,right=.97,top=.77,bottom=.23)
    ax=fig.add_subplot(gs[0,0]);checks=forward['checks']
    for i,c in enumerate(checks):
        marker='s' if c['type']=='prior_box_corner' else 'o'
        ax.plot([i+1,i+1],[c['order5_relative_field_error'],c['relative_field_error']],color='#bcc6d0',lw=.8)
        ax.scatter(i+1,c['order5_relative_field_error'],marker=marker,s=40,facecolors='none',edgecolors='#687887')
        ax.scatter(i+1,c['relative_field_error'],marker=marker,s=45,color=TOP_COLORS[c['H']],edgecolors='#344254',linewidth=.5)
    ax.set_yscale('log');ax.set_ylim(1e-8,2e-4);ax.set_xlim(.3,len(checks)+.7)
    ax.set_xticks([1,4,8,9,12,16,17,20,24]);ax.set_xlabel('已保存的审计点编号：H0 [1–8] / H1 [9–16] / H2 [17–24]')
    ax.set_ylabel('相对散射场误差');ax.set_title('24 点：多项式插值 vs 完整向量 Mie',fontsize=11,pad=10)
    for x in [8.5,16.5]:ax.axvline(x,color='#9facb9',ls=':',lw=.9)
    ax.scatter([],[],facecolors='none',edgecolors='#687887',label='poly5（空心）')
    ax.scatter([],[],color='#147d92',label='poly7（实心，颜色=拓扑）')
    ax.scatter([],[],marker='s',color='#7f8c99',label='先验边角')
    ax.scatter([],[],marker='o',color='#7f8c99',label='新连续材料点')
    ax.legend(loc='upper right',fontsize=8.5,ncol=2)
    ar=fig.add_subplot(gs[0,1]);ind=forward['independent_treams']
    xx=np.arange(len(ind));errs=[c['relative_field_error'] for c in ind]
    ar.scatter(xx,errs,s=65,c=[TOP_COLORS[c['H']] for c in ind],edgecolors='#344254',zorder=3)
    for i,e in enumerate(errs):ar.annotate(f'{e:.2e}',(i,e),xytext=(0,12),textcoords='offset points',ha='center',fontsize=8.8)
    ar.set_xticks(xx);ar.set_xticklabels([f"H{c['H']}\nε={c['theta'][0]:.1f}+{c['theta'][1]:.2f}i" for c in ind],fontsize=8.5)
    ar.set_yscale('log');ar.set_ylim(1e-16,1e-14);ar.set_xlim(-.6,len(ind)-.4)
    ar.set_title('4 个边角点：独立 treams 实现对照',fontsize=11,pad=10)
    ar.set_ylabel('相对散射场误差');ar.grid(axis='x',visible=False)
    fig.suptitle('B1 · 已保存的独立正演检查',fontsize=14,fontweight='bold',y=.975)
    fig.text(.07,.915,
             '左图仅检查 24 个指定材料点上的插值误差；右图仅检查 4 个指定边角点上的实现一致性。\n'
             '这些有限点测试不认证全部先验支持上的误差，也不认证模型外几何或后验尾部。',
             fontsize=9.2,va='top',linespacing=1.6,color='#8a4b00')
    fig.text(.07,.14,
             f"保存的最大 poly7 场误差 = {refinement['max_empirical_relative_field_error']:.3e}；"
             f"独立实现最大场误差 = {refinement['max_independent_treams_relative_field_error']:.3e}。\n"
             '24 + 4 次检查已计入 B1 总预算 384/384。本次制图仅读取结果，新物理调用为 0。',
             fontsize=10,va='top',linespacing=1.6)
    save_high(fig,'fig08_independent_forward_checks.png');return fig


# ----------------------------------------------------------------------------
# Dashboard
# ----------------------------------------------------------------------------
def build_dashboard(config, summary, prior_sens, importance, calib, cases,
                    refinement, diagnostics, mcm):
    primary=[n for n in CASE_ORDER if n in cases]
    refined=["refined_"+n for n in CASE_ORDER if "refined_"+n in cases]
    order_all=primary+refined
    data={}
    for name in order_all:
        c=cases[name];pj=c["probabilities"];res=pj["result"];p=np.asarray(res["topology_probabilities"])
        oh=list(map(int,np.argsort(-p)));bn=base_name(name);ss=summary["cases"].get(bn,{})
        observation=summary["cases"].get("gain_mismatch",{}) if bn=="gain_marginalized" else ss
        is_ref=name.startswith('refined_')
        f_all=fig_representatives(config,pj,name,high=False)
        f_top=fig_representatives(config,pj,name,high=False,selected=oh[:2],include_mix=False)
        f_unc=fig_cavity_uncertainty(config,c,name,high=False)
        data[name]=dict(zh=case_label(name),p=list(map(float,p)),order=oh,
                       is_refined=is_ref,truth=observation.get("truth",{}),noise_sd=observation.get("noise_real_sd"),
                       ref5=(pj.get('reference',{}).get('status') if is_ref else ss['refinement_5_vs_7']['status']),
                       original_ref5=ss['refinement_5_vs_7']['status'],
                       tail=res.get('posterior_predictive_tail_probability'),
                       remainder_modeled=float(pj['top_k']['remaining_modeled_mass']),
                       rep_all=fig_to_b64(f_all,105),rep_top2=fig_to_b64(f_top,110),unc=fig_to_b64(f_unc,100))
    model_rows=[(label,mcm['competitive'][key]) for label,key in [
        ('完整向量 Mie（主结果）','full_topology'),
        ('单锚点全导数线性化','single_anchor_full_derivative_linear'),
        ('单锚点偶极导数校准线性化','single_anchor_dipole_derivative_calibrated_linear'),
        ('偶极截断模型','truncated_dipole_model')]]
    payload=dict(cases=data,images=diagnostics,model_rows=model_rows,
                 calibration_scope=dict(parents=calib['independent_continuous_material_parents'],
                                        noises=calib['noise_replicates'],marginal_material_90_coverage=coverage_marginal(calib)),
                 refined_status=refinement['status'])
    blob=json.dumps(payload,ensure_ascii=False).replace('</',r'<\/')
    doc=_DASH_TEMPLATE.replace('__BLOB__',blob).replace('__ORDER__',json.dumps(order_all,ensure_ascii=False))
    (ROOT/'dashboard.html').write_text(doc,encoding='utf-8')
    return ROOT/'dashboard.html'


_DASH_TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>B1 · 概率 Maxwell 成像</title>
<style>
:root{--bg:#f3f6f8;--paper:#fff;--ink:#1b2c3b;--mut:#5c7182;--line:#d6e0e7;--teal:#147d92;--amber:#ad7227;--red:#a13c26}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.65 -apple-system,"PingFang SC","Hiragino Sans GB",sans-serif}
header{background:#113542;color:white;padding:34px max(24px,calc((100vw - 1200px)/2));}h1{font-size:28px;line-height:1.3;margin:8px 0 14px}h2{font-size:20px;margin:0 0 12px;line-height:1.4}h3{font-size:14px;margin:12px 0 8px;color:var(--mut)}p{margin:8px 0}.eyebrow{letter-spacing:.08em;font-size:12px;color:#b3d8df}.intro{max-width:850px;color:#d3e4e9;font-size:15px}.wrap{max-width:1248px;margin:auto;padding:20px 24px 50px}.card{background:var(--paper);border:1px solid var(--line);border-radius:12px;padding:22px;margin:20px 0;box-shadow:0 3px 15px #19314206}
.scope{border-left:4px solid var(--amber);background:#fff9ef;padding:15px 18px;border-radius:8px;font-size:13px;margin:16px 0}.scope strong{color:#805113}.help{font-size:13px;color:var(--mut);margin-top:12px}.cases{display:flex;flex-wrap:wrap;gap:8px}.case,.view-toggle,.figure-open{appearance:none;cursor:pointer;font:inherit;border:1px solid var(--line);background:white;color:var(--ink);padding:6px 14px;border-radius:20px;font-size:13px}.case.active,.view-toggle.active{background:var(--teal);color:white;border-color:var(--teal)}.case.refined{border-style:dashed}.case.refined.active{background:#ad7227;border-color:#ad7227}.cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin:16px 0}.kv{background:#f4f7f9;border:1px solid var(--line);padding:12px 14px;border-radius:8px}.k{font-size:12px;color:var(--mut)}.v{font-size:15px;font-weight:650;margin-top:4px;overflow-wrap:anywhere}.pill{display:inline-block;padding:2px 8px;border-radius:20px;font-size:11px;line-height:1.5;border:1px solid var(--line)}.pill.ok{color:#1c6656;border-color:#9bbfb4;background:#eff8f5}.pill.no{color:var(--red);background:#fff1e9;border-color:#dfb4a3}.pill.empirical{color:#855516;background:#fff6e6;border-color:#d9bb91}
.prob-row{display:grid;grid-template-columns:190px 1fr 100px;align-items:center;gap:16px;margin:15px 0}.track{height:24px;background:#e9eef3;border-radius:5px;overflow:hidden}.fill{height:100%;border-radius:4px}.value{font-size:19px;font-weight:700;text-align:right;font-variant-numeric:tabular-nums}.prob-label{font-size:14px}.remainder{margin:18px 0 0;padding:12px 16px;background:#f3f5f8;border:1px solid var(--line);border-radius:8px;font-size:14px}.remainder b{font-size:18px}.mass-scope{font-size:12px;color:var(--mut);margin:10px 0}.figure-head{display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;margin:8px 0}.figure-open{font-size:12px;border-radius:6px;color:var(--teal)}img.fig{width:100%;display:block;background:white;margin-top:12px;border:1px solid var(--line);border-radius:6px}.two{display:grid;grid-template-columns:1fr 1fr;gap:18px}.two .card{margin:0}code{font:11px/1.6 ui-monospace,SFMono-Regular,monospace;overflow-wrap:anywhere;background:#eff3f7;padding:2px 4px;border-radius:3px}details{margin:14px 0;font-size:13px;color:var(--mut)}summary{cursor:pointer;color:var(--ink)}table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13px}td,th{text-align:left;padding:9px;border-bottom:1px solid var(--line);font-variant-numeric:tabular-nums}th{background:#f3f6f8}footer{font-size:12px;color:var(--mut);margin-top:30px}dialog{border:0;border-radius:10px;max-width:96vw;width:1500px;max-height:95vh;padding:16px;overflow:auto}dialog::backdrop{background:#112b3ac9}.zoom-head{display:flex;justify-content:space-between;align-items:center;margin:0 0 12px;gap:20px}#zoomImg{width:100%;height:auto;display:block}.table-scroll{overflow-x:auto}button:focus-visible{outline:3px solid #67afc3;outline-offset:3px}
@media(max-width:800px){header{padding:24px}h1{font-size:24px}.wrap{padding:12px}.card{padding:16px}.cards{grid-template-columns:1fr}.prob-row{grid-template-columns:135px 1fr 82px;gap:8px}.value{font-size:16px}.two{grid-template-columns:1fr}.prob-label{font-size:12px}.case{padding:6px 11px;font-size:12px}}
</style></head><body>
<header><div class="eyebrow">B1 · 3D VECTOR MIE · 概率成像</div>
<h1>同一散射观测，保留多种材料解释</h1>
<p class="intro">本面板比较实心、中心空腔、中心与环状空腔三个固定同心模板。散射数据给出各模板的后验质量；每张代表图展示该模板内材料参数的条件均值。</p>
</header>
<main class="wrap">
<div class="scope"><strong>概率的范围。</strong>所有质量都条件于声明的三个模板，三项合计为 100%。该先验在模板外为 0；更广泛物理几何的后验质量未定义、未量化，不能当作第四分支，也不能据此断言模板外几何不可能。圆心与边界固定，未进行定位推断。</div>
<section class="card"><h2>选择观测情形</h2><h3>原始全网格结果（默认）</h3><div class="cases" id="primaryButtons"></div>
<h3>经验细化结果（单独展示）</h3><div class="cases" id="refinedButtons"></div>
<div class="cards" id="caseFacts"></div><p class="help" id="caseStatus"></p></section>
<section class="card"><h2 id="probTitle">Top-2 拓扑事件质量</h2><div id="probRows"></div>
<div class="remainder" id="remainder"></div><p class="mass-scope" id="massScope"></p>
<p class="help">“拓扑事件”包含该模板的整个连续材料参数区域。连续先验下，单幅精确图像的概率为 0；这些百分比属于分支，不能理解为某幅图像的单点概率。</p></section>
<section class="card"><div class="figure-head"><h2 id="repTitle">Top-2 代表材料场</h2><div><button class="view-toggle active" data-view="top2">Top-2（默认）</button> <button class="view-toggle" data-view="all">全部 3 模板与均值</button></div></div>
<div class="figure-head"><span class="help">上排：Re ε；下排：Im ε。同一数量共用色标。三维球模型以二维截面展示。</span><button class="figure-open" data-image="repImg">放大检查</button></div>
<img class="fig" id="repImg" alt="当前情形的拓扑代表材料场"><p class="help">外部与空腔均为真空（Re ε=1、Im ε=0）；青色虚线表示固定空腔边界。事件边际均值可能不对应任何单一模板。</p></section>
<section class="card"><div class="figure-head"><h2 id="uncTitle">空腔占据与材料不确定度</h2><button class="figure-open" data-image="uncImg">放大检查</button></div><img class="fig" id="uncImg" alt="当前情形空腔占据概率与复介电常数边际标准差"><p class="help">左图按分支概率混合固定空腔模板；右图的复材料标准差同时包含模板之间的差异和模板内部材料参数方差。两项都依赖所选模板与先验。</p></section>
<section class="card"><div class="figure-head"><h2>主情形、先验与负对照</h2><button class="figure-open" data-image="casesImg">放大检查</button></div><img class="fig" id="casesImg" alt="主情形拓扑概率、先验敏感性、负对照"><p class="help">改变先验可改变低信息条件下的排序。负对照的真实空腔半径不属于任何声明模板；即便如此，模板内概率仍被归一化到 100%，因此还需检查观测残差和预测诊断。</p></section>
<section class="card"><div class="figure-head"><h2>原始参考失败与采样误差</h2><button class="figure-open" data-image="referenceImg">放大检查</button></div><img class="fig" id="referenceImg" alt="保留原始5/7参考失败和独立重要性采样误差"><p class="help">数值参考要求 Δp ≤ 0.02 且 Δlog BF ≤ 0.1。四个原始 <code>REFERENCE_UNRESOLVED</code> 保留；竞争情形的重要性采样误差为 0.02565，超出容差。有限样本自助区间不能认证未发现的模态。</p></section>
<section class="card"><div class="figure-head"><h2>经验细化：与原始概率对照</h2><button class="figure-open" data-image="refinedImg">放大检查</button></div><img class="fig" id="refinedImg" alt="独立展示经验细化概率和积分及插值收敛诊断"><p class="help"><code>EMPIRICALLY_REFINED_REFERENCE</code> 表示 order17/25 积分与 poly5/7 插值诊断在声明容差内。它近似连续全 Mie 后验；没有统一余项界或未见尾部认证，不改写原始网格失败。</p></section>
<section class="card"><div class="figure-head"><h2>独立正演点检</h2><button class="figure-open" data-image="forwardImg">放大检查</button></div><img class="fig" id="forwardImg" alt="保存的24个插值审计点和4个独立treams实现对照"><p class="help">这里只展示已保存的 24 + 4 次检查；有限点的一致性不保证全部先验支持上的误差。B1 物理预算已用满 384/384，本面板生成没有增加正演。</p></section>
<section class="card"><div class="figure-head"><h2>探索性校准诊断</h2><button class="figure-open" data-image="calibrationImg">放大检查</button></div><img class="fig" id="calibrationImg" alt="24独立连续材料母体、720噪声实现的探索性校准诊断"><p class="help"><b>24 个独立连续材料母体 / 720 次噪声实现</b>。材料边际区间使用 H 边际权重；“已知真 H 条件”行仅作 oracle 描述性对照。红色外圈表示分箱偏差区间不含 0（未校正多重比较）。母体聚类自助区间不包含求解器偏差或表示偏差；结果属于声明模型内的开发诊断。</p></section>
<section class="card"><h2>竞争情形：正演近似对概率的影响</h2><div class="table-scroll"><table id="modelTable"></table></div><p class="help">各列仍是同一三个模板内的事件质量。对照用于查看模型近似如何改变推断；余项来自付费离线全正演网格，未证明降维或加速。</p></section>
<details><summary>数据来源与生成范围</summary><p><code>configs/protocol.json</code>；<code>results/cases/*/probabilities.json</code> 与 <code>material_fields.npz</code>；<code>PRIMARY_SUMMARY.json</code>、<code>PRIOR_SENSITIVITY.json</code>、<code>MODEL_CORRECTION_MECHANISM.json</code>、<code>CALIBRATION_RESULTS.json</code>、<code>importance/RESULT.json</code>、<code>refinement/SUMMARY.json</code>、<code>refinement/FULL_FORWARD_CHECKS.json</code>；本地审阅生成的 <code>CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json</code>。</p><p>这是自包含本地文件：数据和图片均内联，无远程脚本。生成时仅读取已保存的结果；新物理调用 0、训练 0、锁定或历史数组读取 0。A32 仅是来源名称。</p></details>
<footer>本面板提供探索性数值结果与证据范围；数值细化或点检通过均不等于严格后验认证。</footer>
</main>
<dialog id="figureDialog"><div class="zoom-head"><b id="zoomTitle">图像检查</b><button class="figure-open" id="closeZoom">关闭</button></div><img id="zoomImg" alt="放大科学图"></dialog>
<script>
const D=__BLOB__;
const ORDER=__ORDER__;
const TN=['实心','中心空腔','中心+环状空腔'];
const TC=['#496b9b','#147d92','#db9245'];
let current='competitive',view='top2';
function pct(v,n=4){return v<0.00005?(v*100).toExponential(2)+'%':(v*100).toFixed(n)+'%';}
function status(s){if(s==='NUMERICAL_REFERENCE_PASS')return '<span class="pill ok">数值细化 PASS · 无严格证书</span>';if(s==='EMPIRICALLY_REFINED_REFERENCE')return '<span class="pill empirical">经验细化参考 · 无严格证书</span>';return '<span class="pill no">REFERENCE_UNRESOLVED</span>';}
function render(name){
current=name;const c=D.cases[name];const t=c.truth||{};document.querySelectorAll('.case').forEach(b=>b.classList.toggle('active',b.dataset.case===name));
const facts=[['真值模板',(t.H===undefined?'未记录':'H'+t.H+' · '+TN[t.H])+(t.within_inference_support===false?'（真实几何在支持外）':'')],['观测噪声标准差',c.noise_sd==null?'未记录':c.noise_sd.toFixed(4)],[c.is_refined?'经验细化状态':'原始全网格 5/7 状态',status(c.ref5)]];
document.getElementById('caseFacts').innerHTML=facts.map(([k,v])=>`<div class="kv"><div class="k">${k}</div><div class="v">${v}</div></div>`).join('');
document.getElementById('caseStatus').innerHTML=c.is_refined?`此处单独展示经验细化。原始 5/7 状态仍为 ${status(c.original_ref5)}；没有统一插值余项证书。`:'主结果来自原始 order7 全网格积分；参考状态独立保留。';
document.getElementById('probTitle').textContent=c.zh+' · Top-2 拓扑事件质量';
document.getElementById('probRows').innerHTML=c.order.slice(0,2).map((h,i)=>`<div class="prob-row"><div class="prob-label">Top-${i+1} · H${h}<br>${TN[h]}</div><div class="track"><div class="fill" style="width:${c.p[h]*100}%;background:${TC[h]}"></div></div><div class="value">${pct(c.p[h])}</div></div>`).join('');
const h3=c.order[2];document.getElementById('remainder').innerHTML=`保留第三模板 · H${h3} ${TN[h3]}：<b>${pct(c.p[h3])}</b>　<span class="help">实际建模余量，未重归一化</span>`;
document.getElementById('massScope').textContent='Top-2 合计 '+pct(c.p[c.order[0]]+c.p[c.order[1]])+' + 第三模板 '+pct(c.p[h3])+' = 三个声明模板内 100%。更广泛几何质量未量化。';
document.getElementById('repImg').src='data:image/png;base64,'+(view==='top2'?c.rep_top2:c.rep_all);
document.getElementById('repTitle').textContent=c.zh+(view==='top2'?' · Top-2 代表材料场':' · 三个代表场与事件边际均值');
document.getElementById('uncImg').src='data:image/png;base64,'+c.unc;document.getElementById('uncTitle').textContent=c.zh+' · 空腔占据与材料不确定度';
document.body.dataset.currentCase=name;document.body.dataset.representativeView=view;
}
for(const name of ORDER){const b=document.createElement('button');b.className='case'+(D.cases[name].is_refined?' refined':'');b.dataset.case=name;b.textContent=D.cases[name].zh;b.onclick=()=>render(name);document.getElementById(D.cases[name].is_refined?'refinedButtons':'primaryButtons').appendChild(b);}
document.querySelectorAll('[data-view]').forEach(b=>b.onclick=()=>{view=b.dataset.view;document.querySelectorAll('[data-view]').forEach(x=>x.classList.toggle('active',x===b));render(current);});
const maps={casesImg:'cases',referenceImg:'reference',refinedImg:'refined',forwardImg:'forward',calibrationImg:'calibration'};
for(const [id,key] of Object.entries(maps))document.getElementById(id).src='data:image/png;base64,'+D.images[key];
document.getElementById('modelTable').innerHTML='<tr><th>正演模型</th>'+TN.map((t,h)=>'<th>H'+h+' · '+t+'</th>').join('')+'</tr>'+D.model_rows.map(([label,p])=>'<tr><td>'+label+'</td>'+p.map(v=>'<td>'+pct(v)+'</td>').join('')+'</tr>').join('');
const dlg=document.getElementById('figureDialog');document.querySelectorAll('[data-image]').forEach(b=>b.onclick=()=>{const img=document.getElementById(b.dataset.image);document.getElementById('zoomImg').src=img.src;document.getElementById('zoomTitle').textContent=img.alt;dlg.showModal();});document.getElementById('closeZoom').onclick=()=>dlg.close();dlg.addEventListener('click',e=>{if(e.target===dlg)dlg.close();});
render('competitive');
</script></body></html>
"""


def main():
    config, summary, mcm, prior_sens, importance, calib, cases = gather()
    checks = validate(config, cases)
    if not all(c["ok"] for c in checks):
        raise RuntimeError(f"representative reconstruction mismatch: {checks}")

    refinement=load_json(ROOT/'results/refinement/SUMMARY.json')
    forward=load_json(ROOT/'results/refinement/FULL_FORWARD_CHECKS.json')
    review_path=ROOT/'research/delegated/probability_review/CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json'
    review=load_json(review_path) if review_path.exists() else None
    # Each static figure is also rendered once at dashboard resolution.
    diagnostics={}
    for key,fig in [
        ('representatives',fig_representatives(config,cases['competitive']['probabilities'],'competitive')),
        ('top2',fig_top2(summary)),
        ('uncertainty',fig_cavity_uncertainty(config,cases['competitive'],'competitive')),
        ('cases',fig_cases_prior_negative(summary,prior_sens)),
        ('calibration',fig_calibration(calib,review)),
        ('reference',fig_reference(summary,importance)),
        ('refined',fig_refined_comparison(summary,refinement)),
        ('forward',fig_independent_forward(refinement,forward))]:
        diagnostics[key]=fig_to_b64(fig,110)
    dash=build_dashboard(config,summary,prior_sens,importance,calib,cases,refinement,diagnostics,mcm)

    wall = time.perf_counter() - WALL_T0
    cpu = time.process_time() - CPU_T0
    figs = sorted(p.name for p in FIG_DIR.glob("*.png"))
    primary = [n for n in CASE_ORDER if n in cases]
    refined = sorted(n for n in cases if n not in primary)
    receipt = dict(
        task="B1 static scientific figures + self-contained Top-K dashboard",
        title="B1",
        a32_role="source name only (not used as a quantity)",
        executed_at=datetime.now(timezone.utc).isoformat(),
        host=socket.gethostname(), platform=platform.platform(),
        interpreter=str(Path(__import__("sys").executable)),
        python=platform.python_version(),
        wall_seconds=round(wall, 4), cpu_seconds=round(cpu, 4),
        physics_calls=0, budget_consumed=0, nn_training=0,
        locked_or_historical_arrays_read=0,
        forward_or_mie_executed=False,
        numpy_used=str(np.__version__),
        matplotlib_used=matplotlib.__version__,
        cjk_fonts=_cjk,
        inputs_read=[
            "configs/protocol.json",
            "results/cases/*/probabilities.json",
            "results/cases/*/material_fields.npz",
            "results/PRIMARY_SUMMARY.json",
            "results/MODEL_CORRECTION_MECHANISM.json",
            "results/PRIOR_SENSITIVITY.json",
            "results/importance/RESULT.json",
            "results/CALIBRATION_RESULTS.json",
            "results/refinement/SUMMARY.json",
            "results/refinement/FULL_FORWARD_CHECKS.json",
            "src/b1/pilot.py (read to confirm gain_marginalized shares gain_mismatch observation)",
            "research/delegated/probability_review/CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json",
        ],
        outputs_figures=["figures/" + f for f in figs],
        output_dashboard="dashboard.html",
        default_display=dict(case='competitive',top_k=2,renormalized=False),
        calibration_material_coverage_semantics='H-marginal weights; true-H conditional coverage displayed separately as oracle descriptive diagnostic',
        calibration_parent_cluster_intervals_used=bool(review),
        thread_environment={k:os.getenv(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},
        representation_reconstruction_check=checks,
        primary_cases=primary,
        secondary_refined_cases=refined,
        secondary_reference_status={n: cases[n]["probabilities"].get("reference", {}).get("status")
                                    for n in refined},
        preserved_status=dict(
            reference_unresolved_preserved=True,
            unmodeled_support_mass="undefined/unquantified in broader physical universe",
            declared_prior_outside_support_mass=0,
            note="All reported probabilities sum to one conditional on the declared three-template universe; no fourth-branch mass estimate"),
    )
    (DELEG / "receipt.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
                                        encoding="utf-8")

    findings = [
        ("ARTIFACT", f"{len(figs)} static figures under figures/ ({', '.join(figs)})"),
        ("ARTIFACT", f"self-contained dashboard: dashboard.html ({dash.stat().st_size/1e6:.2f} MB, inline base64, no remote JS)"),
        ("CHECK", "rebuilt representatives reproduce saved posterior_mean/void_probability (max err 0.0); "
                  f"ingested {len(primary)} primary + {len(refined)} secondary refined cases"),
        ("FINDING", "primary competitive Top-2 = H1 50.6418% / H2 40.4265%, remaining modeled H0 8.9317%; not renormalized"),
        ("FINDING", "reference: competitive + information_limited PASS; distinguishable, gain_mismatch, gain_marginalized, representation_mismatch REFERENCE_UNRESOLVED (preserved)"),
        ("FINDING", "importance error vs order7 = 0.0257 > 0.02 tolerance (not certified)"),
        ("FINDING", "calibration EXPLORATORY: 24 parent clusters / 720 noises; Brier 0.518, top2 coverage 0.865"),
        ("UNCERTAINTY", "three-template prior assigns zero outside support by definition; broader physical geometry posterior undefined/unquantified"),
        ("UNCERTAINTY", "fixed concentric geometry: center fixed at origin, not inferred; annular void has no unique center"),
        ("NO-PHYSICS", f"physics_calls=0 budget_consumed=0 nn_training=0; wall={wall:.2f}s cpu={cpu:.2f}s"),
    ]
    for tag, msg in findings:
        print(f"[{tag}] {msg}")


if __name__ == "__main__":
    main()
