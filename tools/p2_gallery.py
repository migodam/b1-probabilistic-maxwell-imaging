"""Render existing Joint P2 probability records; performs no physics or inference.

The parent meters this command through tools/p2_job.py. Every displayed mass
comes from a saved package/reference; Top-K is never renormalized here. The
only project import is the pure material decoder/image implementation.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.b1_joint.material import image, decode  # pure, no solver or accounting

BLUE, RED, GOLD, GREY = "#236c93", "#b5454c", "#bd8736", "#63727e"
CASE_LABELS = {
    "solid": "实心观测", "small_cavity": "小中心空腔观测",
    "large_cavity": "大中心空腔观测", "annular_cavity": "环形空腔观测",
    "combined": "中心与环形空腔观测", "prior_boundary": "先验边界压力案例",
    "gain_negative_control": "增益失配负控制",
}
SHORT_LABELS = {
    "solid": "Solid", "small_cavity": "Small central", "large_cavity": "Large central",
    "annular_cavity": "Annular", "combined": "Combined", "prior_boundary": "Prior boundary",
    "gain_negative_control": "Gain mismatch",
}


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#7c8993", "axes.labelcolor": "#24323c",
        "xtick.color": "#344751", "ytick.color": "#344751",
        "text.color": "#24323c", "figure.facecolor": "white",
        "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def save(fig, directory: Path, name: str):
    fig.savefig(directory / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(directory / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    return f"figures/{name}.png"


def event_label(e):
    e = int(e)
    return f"C{e // 3}A{e % 3}"


def case_key(path: Path, package):
    return str(package.get("case_id") or path.stem)


def mass_error_number(value):
    if isinstance(value, (float, int)):
        return float(value)
    if isinstance(value, dict):
        for key in ("absolute", "absolute_error", "event_mass", "numerical_error", "estimated_error"):
            if isinstance(value.get(key), (float, int)):
                return float(value[key])
    return None


def compact(value, length=180):
    text = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return text if len(text) <= length else text[:length] + "…"


def percent(value):
    return "未量化" if value is None else f"{100 * float(value):.2f}%"


def topk_figure(package, key, directory):
    hypotheses = package["hypotheses"]
    if not hypotheses:
        raise ValueError(f"{key}: empty hypotheses cannot be shown as a Top-K result")
    ncols = min(6, len(hypotheses))
    blocks = math.ceil(len(hypotheses) / ncols)
    fig, axes = plt.subplots(2 * blocks, ncols, figsize=(2.5 * ncols + 1.0, 4.65 * blocks + 1.1),
                             squeeze=False, layout="constrained")
    re_norm, im_norm = Normalize(1, 2.8), Normalize(0, .16)
    re_artist = im_artist = None
    for index, hyp in enumerate(hypotheses):
        latent = np.asarray(hyp["representative_latent"], dtype=float)
        axis, field = image(latent, n=161)
        geometry = decode(latent)
        block, column = divmod(index, ncols)
        real_ax, loss_ax = axes[2 * block, column], axes[2 * block + 1, column]
        re_artist = real_ax.imshow(field.real, origin="lower", extent=(axis[0], axis[-1], axis[0], axis[-1]),
                                   cmap="viridis", norm=re_norm, interpolation="nearest")
        im_artist = loss_ax.imshow(field.imag, origin="lower", extent=(axis[0], axis[-1], axis[0], axis[-1]),
                                   cmap="magma", norm=im_norm, interpolation="nearest")
        mass_error=hyp.get('mass_error',{});error=mass_error_number(mass_error)
        if isinstance(mass_error,dict) and isinstance(mass_error.get('MC_residual_SE'),(int,float)):
            error_text=f" · MC SE {100*mass_error['MC_residual_SE']:.3f}pp"
        elif isinstance(mass_error,dict) and isinstance(mass_error.get('quadrature_refinement_max'),(int,float)):
            error_text=f" · quad Δ {100*mass_error['quadrature_refinement_max']:.3f}pp"
        else:error_text="" if error is None else f" ±{100 * error:.2f}pp"
        real_ax.set_title(f"{event_label(hyp['event_id'])} · {100 * hyp['posterior_mass']:.2f}%{error_text}\n"
                          f"r={geometry['central_radius']:.3f}R, h={geometry['annular_halfwidth']:.3f}R", fontsize=10)
        loss_ax.set_title(f"Re ε={latent[2]:.3f}, Im ε={latent[3]:.3f}", fontsize=9)
        for ax in (real_ax, loss_ax):
            ax.set_xticks([-1, 0, 1])
            ax.set_yticks([-1, 0, 1])
            ax.tick_params(labelsize=8)
            ax.set_aspect("equal")
            ax.set_xlabel("x/R", fontsize=8)
            if column == 0:
                ax.set_ylabel("y/R", fontsize=8)
        if column == 0:
            real_ax.text(-.28, .5, "Re ε", transform=real_ax.transAxes, va="center", rotation=90, color=BLUE)
            loss_ax.text(-.28, .5, "Im ε", transform=loss_ax.transAxes, va="center", rotation=90, color=RED)
    for index in range(len(hypotheses), blocks * ncols):
        block, column = divmod(index, ncols)
        axes[2 * block, column].axis("off")
        axes[2 * block + 1, column].axis("off")
    fig.colorbar(re_artist, ax=axes[::2, :].ravel().tolist(), shrink=.8, label="Re ε")
    fig.colorbar(im_artist, ax=axes[1::2, :].ravel().tolist(), shrink=.8, label="Im ε")
    fig.suptitle(f"{SHORT_LABELS.get(key, key)} · material-event representatives, not point probabilities\n"
                 f"Shown mass {percent(package.get('represented_mass'))}; unshown modeled mass "
                 f"{percent(package.get('unshown_modeled_mass'))}; outside support unquantified", fontsize=12)
    return save(fig, directory, f"topk_{key}")


def event_mass_figure(reference, directory):
    keys = [k for k in CASE_LABELS if k in reference.get("cases", {})]
    keys.extend(k for k in reference.get("cases", {}) if k not in keys)
    masses, errors, labels = [], [], []
    for key in keys:
        target = reference["cases"][key].get("C")
        if not target or not target.get("orders"):
            continue
        masses.append(target["orders"][-1]["event_mass"])
        errors.append(target.get("refinement", {}).get("event_mass"))
        labels.append(SHORT_LABELS.get(key, key))
    if not masses:
        return None
    data = np.asarray(masses)
    fig, ax = plt.subplots(figsize=(12.8, .65 * len(data) + 2.4), layout="constrained")
    artist = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(12), [event_label(e) for e in range(12)], fontsize=9)
    ax.set_yticks(range(len(labels)), labels)
    for row in range(len(data)):
        for col in range(12):
            value = data[row, col]
            label = f"{100 * value:.1f}" if value >= .0005 else "<0.05"
            ax.text(col, row, label, ha="center", va="center", fontsize=8,
                    color="white" if value > .45 else "#193544")
    ax.set_xlabel("C: central-radius bin 0 / (0,.10] / (.10,.20] / (.20,rmax]; "
                  "A: annular halfwidth 0 / (0,.06] / (.06,hmax]")
    ax.set_title("Target C · 12 predeclared, disjoint material events\n"
                 "Percent of the full modeled posterior; rows sum to 100%, no Top-K normalization")
    fig.colorbar(artist, ax=ax, shrink=.7, label="Posterior event mass")
    return save(fig, directory, "event_mass_C")


def compensation_figure(compensation, directory):
    records = compensation["records"]
    geometry_keys = list(dict.fromkeys(r["geometry"] for r in records))
    fig, axes = plt.subplots(1, 2, figsize=(12.3, 5.3), layout="constrained")
    width = .24
    for panel, freq in enumerate(([1.6], [1.6, 2.4])):
        ax = axes[panel]
        x = np.arange(len(geometry_keys))
        selected = [r for r in records if r["frequencies"] == freq]
        raw, eps, gain = [], [], []
        for key in geometry_keys:
            by_kind = {r["compensation"]: r for r in selected if r["geometry"] == key}
            eps_record, gain_record = by_kind["epsilon"], by_kind["epsilon_gain"]
            raw.append(eps_record["raw_distance"])
            eps.append(eps_record["full_compensated_distance"])
            gain.append(gain_record["full_compensated_distance"])
        for offset, values, label, color in ((-width, raw, "Raw secant", GREY), (0, eps, "ε compensation", BLUE),
                                              (width, gain, "ε + gain compensation", RED)):
            ax.bar(x + offset, values, width, label=label, color=color)
        ax.axhline(1, color=GOLD, linestyle="--", lw=1, label="unit whitened distance")
        ax.set_yscale("log")
        ax.set_xticks(x, [SHORT_LABELS.get(k, k).replace(" ", "\n", 1) for k in geometry_keys], fontsize=8)
        ax.set_ylabel("Noise-whitened full-Mie distance")
        ax.set_title("kR = " + ", ".join(map(str, freq)))
        ax.grid(axis="y", alpha=.15)
        if panel == 0:
            ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("Structure changes can be absorbed; extra frequency can reveal the residual\n"
                 "Bounded shared gain [.7,1.3]; best-found feasible endpoints, not posterior evidence", fontsize=12)
    return save(fig, directory, "compensation")


def forward_audit_figure(audit, directory):
    checks = audit.get("checks", [])
    if not checks:
        return None
    values = np.asarray([c["whitened_error"] for c in checks])
    fig, ax = plt.subplots(figsize=(7.1, 3.8), layout="constrained")
    ax.hist(values, bins=18, color=BLUE, alpha=.85)
    ax.axvline(audit["max_whitened"], color=RED, lw=1.5, label=f"max = {audit['max_whitened']:.5f}")
    ax.set_xlabel("Noise-whitened surrogate vs independent full-Mie field error")
    ax.set_ylabel("Independent prior probes")
    ax.set_title(f"{audit.get('points', len(checks))} independent probes · empirical forward audit\n"
                 "No uniform bound; probes do not certify every posterior mode")
    ax.legend(fontsize=9)
    return save(fig, directory, "forward_audit")


def downstream_figure(data, directory):
    """Declared figure schema: cases mapping with event_mass_before/after arrays."""
    rows = data.get("cases", {})
    if isinstance(rows, list):
        rows = {r["case_id"]: r for r in rows}
    valid = {k: r for k, r in rows.items() if "event_mass_before" in r and "event_mass_after" in r}
    key = "large_cavity" if "large_cavity" in valid else next(iter(valid), None)
    if key is None:
        return None
    record = valid[key]
    before, after = np.asarray(record["event_mass_before"]), np.asarray(record["event_mass_after"])
    fig, ax = plt.subplots(figsize=(11.2, 4.3), layout="constrained")
    x = np.arange(12)
    ax.bar(x - .18, before, .36, color=BLUE, label="Maxwell posterior")
    ax.bar(x + .18, after, .36, color=RED, label="After independent noisy void-volume sensor")
    ax.set_xticks(x, [event_label(i) for i in x])
    ax.set_ylabel("Full modeled posterior event mass")
    ax.set_title(f"{SHORT_LABELS.get(key, key)} · downstream reweighting, not repeated scattering evidence\n"
                 "Context likelihood: c = normalized void volume + N(0, .01²)")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=.15)
    return save(fig, directory, "downstream")


def method_figure(data, directory):
    """Keep fresh proxy densities and one full-Mie development case separate."""
    methods = data.get("methods", [])
    if isinstance(methods, dict):
        methods = [dict(value, method=key) for key, value in methods.items()]
    raw = [r for r in methods if r.get("target") == "fresh common proxy target"
           and "event_mass_MAE" in r]
    full = [r for r in methods if r.get("target") == "full-Mie numerical C target"
            and "event_mass_MAE" in r]
    if not raw or not full:
        return None

    raw_order = ["full_distill", "diag_distill", "gaussian_distill", "full_joint", "diag_joint"]
    raw_labels = ["Full-cov.\ndistillation", "Diagonal\ndistillation", "Single Gaussian\ndistillation",
                  "Full-cov.\njoint training", "Diagonal\njoint training"]
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.7))
    seed_colors = [BLUE, GOLD, GREY]
    left_means = []
    for i, method in enumerate(raw_order):
        rows = sorted((r for r in raw if r["method"].split(":")[0] == method), key=lambda r: r["method"])
        if len(rows) != 3:
            raise ValueError(f"Fresh density figure needs the three reported seeds for {method}")
        values = [r["event_mass_MAE"] for r in rows]
        mean = float(np.mean(values))
        left_means.append(mean)
        axes[0].bar(i, mean, width=.7, color=BLUE if "distill" in method else GREY, alpha=.15)
        for j, value in enumerate(values):
            axes[0].scatter(i + (j - 1) * .13, value, color=seed_colors[j], s=31, zorder=4,
                            label=f"Seed {rows[j]['method'].split(':')[-1]}" if i == 0 else None)
        axes[0].plot([i - .25, i + .25], [mean, mean], color="#263d4a", lw=1.7, zorder=5)
        axes[0].text(i, max(values) + .002, f"mean {mean:.4f}", ha="center", fontsize=8)
    axes[0].set_xticks(np.arange(5), raw_labels, fontsize=8)
    axes[0].set_ylim(0, max(r["event_mass_MAE"] for r in raw) * 1.25)
    axes[0].set_ylabel("Event-mass MAE vs matched proxy reference")
    axes[0].set_title("Fresh proxy posterior · 64 parents × 2 noise realizations\n"
                      "Raw density; 3 seed dots and their mean", fontsize=11)
    axes[0].legend(fontsize=7, loc="upper left")

    full_order = ["full_distill full-IS64", "diag_distill full-IS64", "gaussian_distill full-IS64",
                  "generic SMC64:261013211", "generic SMC64:261013212"]
    right_labels = ["Full-cov.\nIS64", "Diagonal\nIS64", "Gaussian\nIS64", "SMC64\nseed 1", "SMC64\nseed 2"]
    ordered = [{r["method"]: r for r in full}[method] for method in full_order]
    right_values = [r["event_mass_MAE"] for r in ordered]
    axes[1].bar(np.arange(5), right_values, width=.66, color=[BLUE, GOLD, GREY, RED, "#c86c73"])
    for i, row in enumerate(ordered):
        axes[1].text(i, row["event_mass_MAE"] + .0012,
                     f"{row['event_mass_MAE']:.4f}\n{row['full_calls']} full calls", ha="center", fontsize=8)
    axes[1].set_xticks(np.arange(5), right_labels, fontsize=8)
    axes[1].set_ylim(0, max(right_values) * 1.38)
    axes[1].set_ylabel("Event-mass MAE vs matched full-Mie reference")
    axes[1].set_title("Full-Mie C posterior · large_cavity development case\n"
                      "Finite correction / sampling; actual physics budgets differ", fontsize=11)
    for ax in axes:
        ax.grid(axis="y", alpha=.15)

    pool = data.get("matched_corrected_proposals", {})
    cost_notes = []
    if isinstance(pool.get("wall_seconds"), (int, float)):
        cost_notes.append(f"Three IS64 runs pooled: {pool['wall_seconds']:.2f}s / "
                          f"{pool.get('actual_total_physics', 'unreported')} full calls")
    ledger = directory.parent / "results/JOB_LEDGER.jsonl"
    if ledger.exists():
        smc_jobs = [json.loads(line) for line in ledger.read_text(encoding="utf-8").splitlines() if line.strip()]
        smc_jobs = [r for r in smc_jobs if r.get("id") == "p2-generic_smc" and r.get("status") == "PASS"
                    and isinstance(r.get("wall_seconds"), (int, float))]
        if smc_jobs:
            cost_notes.append(f"Two SMC64 runs pooled: {sum(r['wall_seconds'] for r in smc_jobs):.2f}s / "
                              f"{sum(r.get('new_physics', 0) for r in smc_jobs)} full calls; no per-seed wall receipt")
    fig.suptitle("Comparisons are within each panel's target, never between panels", fontsize=13, y=.99)
    fig.text(.06, .09, " · ".join(cost_notes) if cost_notes else "Pooled online timing not reported; consult original cost ledger.",
             fontsize=8, color=GREY)
    fig.text(.06, .055, "Teacher, paid response bank and density training are additional amortized costs. "
             "One development case cannot establish global neural superiority.", fontsize=8, color=GREY)
    fig.tight_layout(rect=(0, .15, 1, .94), w_pad=3)
    return save(fig, directory, "method_comparison")


def proposal_timing_figure(data, directory):
    """Saved 64-draw acquisition timing, excluding physical likelihood correction."""
    methods = data.get("methods", [])
    medians = data.get("medians", {})
    if not methods or not medians:
        return None
    labels = {"full_distill": "Full-cov. distillation", "diag_distill": "Diagonal distillation",
              "gaussian_distill": "Single Gaussian distillation"}
    rows = data.get("rows", [])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    for ax, key, label, color in ((axes[0], "cold_load_infer_sample_logq_seconds", "Cold model reload + condition + 64 draws + log q", BLUE),
                                  (axes[1], "warm_sample_logq_seconds", "Warm condition reuse + fresh 64 draws + log q", GOLD)):
        values = [1000 * medians[method][key] for method in methods]
        ax.barh(np.arange(len(methods)), values, color=color, alpha=.24)
        for i, method in enumerate(methods):
            repeated = [1000 * r[key] for r in rows if r["method"] == method]
            offsets = np.linspace(-.12, .12, len(repeated))
            ax.scatter(repeated, i + offsets, color=color, s=20, zorder=3)
            ax.plot([values[i], values[i]], [i - .22, i + .22], color="#263d4a", lw=1.6)
        ax.set_yticks(np.arange(len(methods)), [labels.get(m, m) for m in methods], fontsize=8)
        ax.set_xlabel("Measured milliseconds; dots = 5 rotated-order repetitions")
        ax.set_title(label, fontsize=10)
        ax.grid(axis="x", alpha=.15)
        ax.invert_yaxis()
    fig.suptitle("Proposal acquisition only · large_cavity · frozen model seed 26101321", fontsize=12)
    fig.text(.06, .055, "No Maxwell or likelihood correction in these timings; no OS-cache flush. "
             "Teacher/bank/training and full correction remain separate costs.", fontsize=8, color=GREY)
    fig.tight_layout(rect=(0, .10, 1, .94), w_pad=2)
    return save(fig, directory, "proposal_timing")

def material_fields_figure(data, directory):
    names=list(data['case_ids']);axis=np.linspace(-1.05,1.05,161);X,Y=np.meshgrid(axis,axis);radius=np.hypot(X,Y)
    fig,axes=plt.subplots(2,3,figsize=(11.2,7),layout='constrained')
    for row,name in enumerate(['solid','large_cavity']):
        i=names.index(name);r=data['radius'];values=[data['vacuum_probability'][i],data['epsilon_mean'][i].real,data['epsilon_mean'][i].imag]
        for col,(v,title,cmap,lo,hi) in enumerate(zip(values,['P(exact vacuum)','Posterior mean Re ε','Posterior mean Im ε'],['Greys','viridis','magma'],[0,1,0],[1,2.8,.16])):
            field=np.interp(radius,r,v);field[radius>1]=np.nan;ax=axes[row,col];art=ax.imshow(field,origin='lower',extent=(-1.05,1.05,-1.05,1.05),cmap=cmap,vmin=lo,vmax=hi);ax.set_title(f'{SHORT_LABELS[name]} · {title}');ax.set_xlabel('x/R');ax.set_ylabel('y/R');fig.colorbar(art,ax=ax,shrink=.7)
    fig.suptitle('Position probabilities and material means · weighted radial posterior\nEvent masses corrected; finite proxy conditional particles, no open-world certificate',fontsize=12)
    return save(fig,directory,'material_probability_fields')


def uncertainty_table(uncertainty):
    if not isinstance(uncertainty, dict):
        return f"<pre>{html.escape(json.dumps(uncertainty, ensure_ascii=False, indent=2))}</pre>"
    names = ["中心空腔半径 / R", "环形空腔半宽 / R", "Re ε", "Im ε"]
    mean, sd = uncertainty.get("mean"), uncertainty.get("SD")
    quantiles = uncertainty.get("quantiles_05_50_95")
    if not isinstance(mean, list) or not isinstance(sd, list) or len(mean) != 4 or len(sd) != 4:
        return f"<pre>{html.escape(json.dumps(uncertainty, ensure_ascii=False, indent=2))}</pre>"
    rows = []
    for i, name in enumerate(names):
        q = "—"
        if isinstance(quantiles, list) and len(quantiles) == 4:
            q = " / ".join(f"{float(v):.4g}" for v in quantiles[i])
        rows.append(f"<tr><td>{name}</td><td>{float(mean[i]):.4g}</td><td>{float(sd[i]):.4g}</td><td>{q}</td></tr>")
    return "<table><thead><tr><th>事件内材料量</th><th>均值</th><th>SD</th><th>q05 / q50 / q95</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"


def hypothesis_table(package):
    rows = []
    for hyp in package["hypotheses"]:
        mass, prior = hyp["posterior_mass"], hyp.get("prior_mass")
        odds = "—" if prior is None or prior <= 0 or prior >= 1 or mass >= 1 else f"{(mass / (1 - mass)) / (prior / (1 - prior)):.3g}"
        latent = np.asarray(hyp["representative_latent"])
        geo = decode(latent)
        details = uncertainty_table(hyp.get("material_uncertainty", {}))
        details += "<details><summary>Forward prediction、补偿与数值误差</summary><pre>" + html.escape(json.dumps({
            k: hyp.get(k) for k in ("forward_prediction", "compensation_diagnostics", "indistinguishable_from", "mass_error")
        }, ensure_ascii=False, indent=2)) + "</pre></details>"
        rows.append(f"<tr><td><b>{event_label(hyp['event_id'])}</b><br><span class='quiet'>{html.escape(compact(hyp.get('event_definition'), 190))}</span></td>"
                    f"<td><b>{percent(mass)}</b><br><span class='quiet'>误差：{html.escape(compact(hyp.get('mass_error'), 100))}</span></td>"
                    f"<td>{percent(prior)}<br><span class='quiet'>事件 odds 倍率 {odds}</span></td>"
                    f"<td>r={geo['central_radius']:.4f}R<br>h={geo['annular_halfwidth']:.4f}R<br>ε={latent[2]:.4f}+{latent[3]:.4f}i</td>"
                    f"<td><details><summary>条件材料分布与预测</summary>{details}</details></td></tr>")
    return "<div class='table-wrap'><table><thead><tr><th>不重叠材料事件</th><th>后验质量</th><th>先验与 evidence 信息</th><th>代表材料</th><th>事件内部的不确定性</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"


def write_html(directory, packages, assets, reference, compensation, qualification, audit, optional):
    panels, options = [], []
    keys = sorted(packages, key=lambda k: (k != "large_cavity", k != "solid", list(CASE_LABELS).index(k) if k in CASE_LABELS else 99))
    for key in keys:
        package, package_path = packages[key]
        options.append(f"<option value='{html.escape(key)}'>{html.escape(CASE_LABELS.get(key, key))}</option>")
        figure = assets[f"topk_{key}"]
        stem = Path(figure).stem
        numerical = html.escape(compact(package.get("numerical_status"), 800))
        panels.append(f"""<section class="case-panel" id="case-{html.escape(key)}">
          <div class="case-top"><div><span class="eyebrow">{html.escape(str(package['observation_id']))}</span><h3>{html.escape(CASE_LABELS.get(key, key))}</h3></div>
          <a class="button" href="results/packages/{html.escape(package_path.name)}">机器可读 JSON ↗</a></div>
          <p class="target"><b>Target：</b>{html.escape(compact(package['target_identity'], 1000))}</p>
          <div class="mass-strip"><span>已展示 <b>{percent(package.get('represented_mass'))}</b></span><span>已建模、未展示 <b>{percent(package.get('unshown_modeled_mass'))}</b></span><span>模型外 <b>无法量化</b></span></div>
          <a href="figures/{stem}.pdf"><img class="science" src="{figure}" alt="{html.escape(key)}：各材料事件的代表性实部与损耗截面，概率不重新归一化"></a>
          <p class="caption">概率属于整个材料事件；图中材料仅是该事件的代表。C 表示中心空腔半径分箱，A 表示环形空腔半宽分箱；不同分箱不自动等于分离的 posterior mode。所有图共用实部与损耗色标。</p>
          {hypothesis_table(package)}
          <p class="numerical"><b>数值状态：</b>{numerical}<br><b>未探索质量：</b>{html.escape(compact(package.get('unexplored_mass_status'), 500))}<br><b>模型适用性：</b>{html.escape(compact(package.get('model_adequacy_status'), 500))}</p>
          <p><a href="results/packages/{html.escape(package_path.stem)}.npz">加权 posterior particles / NPZ</a> · <a href="results/packages/{html.escape(package_path.name)}">完整 forward、事件定义与 context 接口</a></p>
        </section>""")
    ref_converged = sum(bool(v.get("C", {}).get("independent_adaptive", {}).get("converged")) for v in reference.get("cases", {}).values())
    checks = qualification.get("checks", [])
    max_multipole = max((c["n14_n18_relative"] for c in checks), default=float("nan"))
    max_recursion = max((c["recursive_relative"] for c in checks), default=float("nan"))
    compensation_rows = []
    for r in compensation.get("records", []):
        if r["compensation"] != "epsilon_gain":
            continue
        compensation_rows.append(f"<tr><td>{html.escape(CASE_LABELS.get(r['geometry'], r['geometry']))}</td><td>{', '.join(map(str,r['frequencies']))}</td>"
                                 f"<td>{r['raw_distance']:.4f}</td><td>{r['full_compensated_distance']:.4f}</td>"
                                 f"<td>{100*r['absorbed_fraction']:.1f}%</td><td>{r['gain']:.4f}</td><td>{r['optimized_epsilon'][0]:.4f}+{r['optimized_epsilon'][1]:.4f}i</td></tr>")
    def fig(name, caption):
        source = assets.get(name)
        if not source:
            return ""
        return f"<figure><a href='figures/{name}.pdf'><img class='science' src='{source}' alt='{html.escape(caption)}'></a><figcaption>{caption}</figcaption></figure>"
    downstream_text = "<p>下游结果尚未写入本次输入；不得由图库推断效果。请查 <a href='DOWNSTREAM_DECISION_TEST.md'>下游实验记录</a>。</p>"
    if assets.get("downstream"):
        downstream_text = ("<p>额外上下文是声明似然的独立、带噪声空腔体积传感器，而非直接揭示 truth。"
                           "更新使用事件内样本；只给代表图会遗漏条件变化，未展示质量也必须保留。</p>" +
                           fig("downstream", "相同 Maxwell 数据；独立几何传感器通过 posterior reweighting 改变材料事件质量。"))
    method_text = ("<p class='small'>方法图按 target 分开：左侧是 64 个 fresh parents × 两次噪声的 proxy posterior，"
                   "比较五种原始密度、三个训练 seeds；右侧是在单个 large_cavity development 案例上的 full-Mie C posterior，"
                   "比较三种 IS64 修正和两次 SMC64。只允许同一面板内部比较，不能跨面板或据此宣称整体 NN 优越性。</p>" +
                   fig("method_comparison", "左：三个 seed 点与均值，监督信息量差异单独解释。右：有限 full-likelihood 修正/SMC 的事件误差及真实 64/295/301 调用数。时间仅显示有原始回执的合计；teacher、bank 和训练 setup 另计。") +
                   fig("proposal_timing", "仅测量 frozen proposal 的 64 次 draw 与 log q：五次轮换顺序，冷加载与 warm condition 复用分列；没有 Maxwell solve 或 likelihood correction，不能直接当作完整 inference 的速度。"))
    content = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>B1 × C1 Joint P2 · 保留、解释并传递材料歧义</title>
    <style>
    :root{{--ink:#1c303c;--muted:#5e727e;--blue:#236c93;--red:#b5454c;--line:#d7e0e5;--paper:#f7f9fa}}*{{box-sizing:border-box}}
    body{{margin:0;color:var(--ink);background:var(--paper);font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Noto Sans CJK SC",sans-serif;line-height:1.75}}
    main{{max-width:1240px;margin:auto;padding:44px 28px 70px}}a{{color:var(--blue);text-underline-offset:3px}}header{{max-width:940px;margin-bottom:40px}}
    h1{{font-size:clamp(28px,4vw,43px);line-height:1.3;letter-spacing:-.03em;margin:14px 0 22px}}h2{{font-size:25px;margin:35px 0 12px}}h3{{font-size:23px;margin:5px 0}}
    .eyebrow{{font-size:12px;font-weight:650;letter-spacing:.06em;color:var(--blue)}}.lede{{font-size:19px;color:#405762}}.quiet,.caption,figcaption{{color:var(--muted);font-size:13px}}
    .nav{{display:flex;gap:18px;flex-wrap:wrap;font-size:14px}}.notice{{padding:18px 23px;background:#eaf1f5;border-left:4px solid var(--blue);border-radius:4px;margin:20px 0}}
    .case-panel{{background:white;border:1px solid var(--line);padding:25px;border-radius:10px;margin:18px 0}}.case-panel[hidden]{{display:none}}
    .case-top{{display:flex;gap:18px;justify-content:space-between;align-items:center}}.button{{font-size:13px;padding:8px 14px;border:1px solid var(--line);border-radius:6px;text-decoration:none;white-space:nowrap}}
    select{{font:inherit;padding:9px 14px;min-width:300px;background:white;border:1px solid #a7bbc7;border-radius:6px;color:var(--ink)}}.target{{font-size:13px;overflow-wrap:anywhere}}
    .mass-strip{{display:flex;gap:26px;flex-wrap:wrap;background:#f3f6f8;padding:12px 15px;border-radius:6px;font-size:14px}}.mass-strip b{{font-variant-numeric:tabular-nums}}
    .science{{display:block;width:100%;height:auto;max-width:1180px;margin:18px auto;border:1px solid #eef1f3;background:white}}figure{{margin:20px 0 28px}}
    table{{border-collapse:collapse;width:100%;font-size:13px;min-width:650px}}th,td{{border-bottom:1px solid var(--line);padding:11px 10px;text-align:left;vertical-align:top}}th{{background:#eef3f6;color:#355463;font-weight:650}}
    .table-wrap{{overflow-x:auto}}td b{{font-variant-numeric:tabular-nums}}details{{font-size:13px}}summary{{color:var(--blue);cursor:pointer}}details table{{min-width:390px}}
    pre{{white-space:pre-wrap;word-break:break-word;font-size:12px;background:#f2f5f7;padding:12px;max-height:450px;overflow:auto}}.numerical{{font-size:13px;background:#faf3eb;padding:14px;border-left:3px solid #bd8736;overflow-wrap:anywhere}}
    .columns{{display:grid;grid-template-columns:1fr 1fr;gap:24px}}.small{{font-size:14px}}footer{{border-top:1px solid var(--line);padding-top:22px;margin-top:45px;font-size:13px;color:var(--muted)}}
    @media(max-width:700px){{main{{padding:25px 14px}}.case-panel{{padding:16px}}.case-top{{align-items:flex-start;flex-direction:column}}.columns{{grid-template-columns:1fr}}select{{width:100%;min-width:0}}}}
    </style></head><body><main>
    <header><span class="eyebrow">B1 × C1 · JOINT P2 · 2026-10-13</span><h1>同一份散射测量，<br>可能对应哪些材料世界？</h1>
    <p class="lede">保留具有显著后验质量的物理解释，揭示材料与校准自由度如何吸收结构差异，并把不确定性传递给下游判断。</p>
    <nav class="nav"><a href="RESEARCH_REPORT_ZH.md">中文研究报告</a><a href="START_HERE.md">复现入口</a><a href="METHOD_COMPARISON.md">方法比较</a><a href="RESOURCE_SUMMARY.md">完整费用</a><a href="../../">B1 项目首页</a></nav></header>
    <div class="notice"><b>概率的含义：</b>每个数字对应预先定义的、互不重叠的材料事件；它不是一张连续图像的单点概率。图中候选不会重新归一化到 100%。当前物理支持是三维矢量 Maxwell 的同心分层球，展示球对称材料的中心截面；一般偏心、多物体材料世界的概率尚无法量化。</div>
    <h2>一份观测，多个材料解释</h2><p>选择观测案例。各候选必须同时保留事件内部的几何与复介电常数分布；不同半径分箱也可能属于同一条连续 ridge。</p>
    <label for="case-picker" class="quiet">观测案例</label><br><select id="case-picker">{''.join(options)}</select>{''.join(panels)}
    <h2>重要材料事件是否被保留？</h2>
    <p>下图展示 target C：ε 与有界、跨频共享 gain 均未知。每一行是十二个完整事件的数值后验质量，而非 Top-K 条件概率。</p>
    {fig('event_mass_C','C0：无中心空腔；C1/C2/C3：中心半径逐级增加。A0：无环形空腔；A1/A2：环形半宽逐级增加。事件划分固定，不能用一个巨大事件冒充覆盖。')}
    <h2>结构变化被什么吸收？</h2><p>原始与补偿后距离采用相同噪声 whitening。ε 补偿改变材料本构，gain 补偿改变仪器响应，两者分别记录；低 residual 只是辨识诊断，不能代替 proper prior 下的边缘 evidence。</p>
    {fig('compensation','灰：raw secant；蓝：ε 补偿；红：ε + gain 补偿。右图追加 kR=2.4 频率。距离 1 为单位 whitening 参考尺度，不是通用检验阈值；最小值为已验证可行点的 best-found upper bound。')}
    <details><summary>查看已验证的 full-Mie ε + gain 补偿数据</summary><div class="table-wrap"><table><thead><tr><th>材料变化</th><th>kR</th><th>raw 距离</th><th>comp 距离</th><th>吸收比例</th><th>gain</th><th>ε</th></tr></thead><tbody>{''.join(compensation_rows)}</tbody></table></div></details>
    <h2>额外情景信息如何帮助下游？</h2>{downstream_text}
    <p class="small"><a href="DOWNSTREAM_DECISION_TEST.md">MAP、Top-K、完整 weighted samples 与条件样本的决策对照</a> · <a href="results/DOWNSTREAM.json">原始更新数据</a></p>
    <h2>概率参考与方法的实际限制</h2><div class="columns"><div><p>显示概率使用保存的连续 prior 与数值 response target。材料结构的原子质量由 inactive latent 的先验测度产生；CDF、均值、SD 与 quantiles 保持同一 target。</p>
    <p class="small">独立自适应 reference 的内部收敛标记：{ref_converged}/{len(reference.get('cases',{}))} 个 target-C 案例。该标记与逐阶差、full-forward fidelity 分开报告，不能把小差异自动称为全域认证。</p></div>
    <div><p class="small">物理维度：{html.escape(qualification.get('physical_dimension','未声明'))}<br>最高 n14/n18 相对场差：{max_multipole:.3g}<br>独立递推/边界解最大相对差：{max_recursion:.3g}<br>独立 prior probes 最大 whitening 误差：{audit.get('max_whitened',float('nan')):.5g}<br>Uniform certificate：<b>没有</b></p></div></div>
    {fig('forward_audit','独立 prior probes 的 surrogate / full-Mie 响应误差。该经验核验不保证所有 posterior modes 均被发现，也不构成现实材料模型认证。')}{method_text}
    <p><a href="POSTERIOR_REFERENCE.md">参考误差与连续积分</a> · <a href="MODE_COVERAGE.md">false exclusion 与 false confidence</a> · <a href="FAILURES.md">失败与 NOT_RUN</a> · <a href="THEORY_STATUS.md">理论可声明范围</a> · <a href="GATE_DECISION.md">最终科研裁决</a></p>
    <footer>所有图仅由已保存的实验数据与纯材料 decoder 渲染。Galer​​y 不调用 forward solver、不做 inference、不产生新的候选概率。原始 B1 384 次物理账本与 C1-P1 保留只读。<br><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging">B1 仓库</a> · <a href="https://github.com/migodam/c1-physics-state-generative-imaging">C1 仓库</a></footer>
    </main><script>const picker=document.getElementById('case-picker');function showCase(){{document.querySelectorAll('.case-panel').forEach(el=>{{el.hidden=el.id!=='case-'+picker.value;}});}}picker.addEventListener('change',showCase);showCase();</script></body></html>"""
    (directory / "index.html").write_text(content, encoding="utf-8")


def main():
    started = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT / "research/joint_p2")
    args = parser.parse_args()
    directory, results = args.root, args.root / "results"
    figure_dir = directory / "figures"
    paths = sorted((results / "packages").glob("*.json"))
    if not paths:
        raise SystemExit("No packages available. Parent must create actual posterior packages before rendering.")
    packages = {}
    for path in paths:
        package = read(path)
        if not all(k in package for k in ("observation_id", "target_identity", "hypotheses")):
            continue
        key = case_key(path, package)
        if not re.fullmatch(r"[A-Za-z0-9_-]+", key):
            raise ValueError("case_id must be a filename-safe explicit identifier")
        packages[key] = (package, path)
    if not packages:
        raise SystemExit("No hypothesis package records matched the declared interface.")
    packages={key:packages[key] for key in CASE_LABELS if key in packages}
    figure_dir.mkdir(parents=True, exist_ok=True)
    style()
    reference = read(results / "PILOT_REFERENCE.json")
    compensation = read(results / "COMPENSATION.json")
    qualification = read(results / "PHYSICS_QUALIFICATION.json")
    audit = read(results / "INDEPENDENT_FORWARD_AUDIT.json")
    assets, optional = {}, {}
    for key, (package, _) in packages.items():
        assets[f"topk_{key}"] = topk_figure(package, key, figure_dir)
    assets["event_mass_C"] = event_mass_figure(reference, figure_dir)
    assets["compensation"] = compensation_figure(compensation, figure_dir)
    assets["forward_audit"] = forward_audit_figure(audit, figure_dir)
    if (results/'MATERIAL_FIELDS.npz').exists():
        with np.load(results/'MATERIAL_FIELDS.npz') as data:assets['material_probability_fields']=material_fields_figure(data,figure_dir)
    for filename, name, renderer in (("DOWNSTREAM.json", "downstream", downstream_figure),
                                       ("METHOD_COMPARISON.json", "method_comparison", method_figure),
                                       ("PROPOSAL_TIMING.json", "proposal_timing", proposal_timing_figure)):
        path = results / filename
        if path.exists():
            data = read(path)
            assets[name] = renderer(data, figure_dir)
            optional[name] = "rendered" if assets[name] else "input lacks declared figure schema; read original JSON"
        else:
            optional[name] = "NOT_RUN: result not available to renderer"
    write_html(directory, packages, assets, reference, compensation, qualification, audit, optional)
    if assets.get('material_probability_fields'):
        target=directory/'index.html';content=target.read_text();content=content.replace('<footer>',"<h2>位置概率与条件材料场</h2><p>严格真空的径向位置概率，以及材料实部/损耗后验均值。均值图是分布摘要，不是另一个具有点概率的候选；事件内仍保留不确定性，数值与支持范围见原始文件。</p><figure><img class='science' src='figures/material_probability_fields.png' alt='位置概率与复材料后验均值'><figcaption>全部12事件的质量与有限event-conditional proxy particles；模型外质量未量化。</figcaption></figure><footer>");target.write_text(content)
    manifest = dict(renderer="tools/p2_gallery.py", physics_evaluations=0, inference_evaluations=0,
                    masses_renormalized=False, cases=list(packages), images={k:v for k,v in assets.items() if v},
                    optional_inputs=optional, wall_seconds=time.perf_counter() - started)
    (figure_dir / "GALLERY_RENDER.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
