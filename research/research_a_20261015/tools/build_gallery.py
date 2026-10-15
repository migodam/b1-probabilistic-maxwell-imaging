#!/usr/bin/env python3
"""Render Research A from saved arrays only; no solver imports or evaluations.

Truth is read solely for offline pictures/scores.  The second displayed case is
chosen by online structural count, then full-data residual, then original order.
Missing experiments remain absent; this program cannot manufacture candidates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ARMS = {
    "finite_static": "Finite edit + compensation",
    "multistart": "Independent multistart",
    "population": "Repulsive population",
    "mdn": "Frozen C1 proposal + compensation",
    "adaptive": "Adaptive next solution",
}
COLORS = dict(zip(ARMS, ["#2075b0", "#d87c24", "#3a8d63", "#9660a4", "#cb4549"]))


def load(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def finite(x):
    return isinstance(x, (int, float)) and np.isfinite(x)


def mean(values):
    values = [float(v) for v in values if finite(v)]
    return float(np.mean(values)) if values else None


def decoded_image(p, n=193):
    r, h, er, ei = map(float, p)
    axis = np.linspace(-1.04, 1.04, n)
    xx, yy = np.meshgrid(axis, axis)
    radius = np.hypot(xx, yy)
    material = np.full(radius.shape, complex(er, ei))
    material[radius < r] = 1
    if h > 0:
        material[np.abs(radius - .65) < h] = 1
    material[radius > 1] = complex(np.nan, np.nan)
    return axis, material


def structural_count(archive):
    representatives = []
    for c in archive:
        p = c["p"]
        if all(abs(p[0] - a[0]) >= .05 or abs(p[1] - a[1]) >= .025
               for a in representatives):
            representatives.append(p)
    return len(representatives)


def distance(a, b):
    # Exact volume-weighted field distance in the declared radial support.
    edges = sorted({0., 1., float(a[0]), float(b[0]),
                    *([.65-a[1], .65+a[1]] if a[1] > 0 else []),
                    *([.65-b[1], .65+b[1]] if b[1] > 0 else [])})
    total = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (lo+hi)/2
        def eps(p):
            return 1 if s < p[0] or (p[1] > 0 and abs(s-.65) < p[1]) else complex(p[2], p[3])
        total += (hi**3-lo**3) * abs(eps(a)-eps(b))**2
    return float(np.sqrt(total)/abs(1.5+.1j))


def save(fig, folder, name):
    fig.savefig(folder/(name+".png"), dpi=165, bbox_inches="tight", facecolor="white")
    fig.savefig(folder/(name+".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def topk_figure(case, method, truth, figures, name):
    archive = method["archive"]
    selected = method.get("selected_indices", list(range(len(archive))))[:8]
    items = []
    if truth is not None:
        items.append(("Truth: offline only", truth["p"], None))
    if case.get("base"):
        items.append(("Common coarse base", case["base"]["p"], case["base"].get("residual")))
    for j, i in enumerate(selected):
        c = archive[i]
        items.append((f"Candidate {j+1} / archive {i}", c["p"], c.get("residual")))
    if not items:
        return
    fig, axes = plt.subplots(2, len(items), figsize=(2.15*len(items), 5.4), squeeze=False)
    for col, (label, p, residual) in enumerate(items):
        _, image = decoded_image(p)
        for row, part, lo, hi, cmap in ((0, image.real, 1., 2.8, "viridis"),
                                       (1, image.imag, 0., .16, "magma")):
            im = axes[row, col].imshow(part, origin="lower", extent=(-1.04,1.04,-1.04,1.04),
                                      vmin=lo, vmax=hi, cmap=cmap)
            axes[row, col].set_xticks([]); axes[row, col].set_yticks([])
            if col == 0:
                axes[row, col].set_ylabel("Re permittivity" if row == 0 else "Im permittivity")
        label += f"\nr={p[0]:.3f}, h={p[1]:.3f}\nε={p[2]:.3f}+{p[3]:.3f}i"
        if finite(residual): label += f"\nwhitened residual={residual:.3f}"
        axes[0, col].set_title(label, fontsize=8)
    fig.suptitle(f"{case['id']} | {ARMS.get(method['arm'], method['arm'])}\n"
                 "3D concentric vector-Mie model, central cross-section; candidate-only, probabilities null",
                 fontsize=11, y=1.035)
    fig.subplots_adjust(left=.025,right=.94,top=.72,bottom=.025,wspace=.07,hspace=.12)
    fig.colorbar(axes[0,0].images[0], cax=fig.add_axes([.952,.407,.009,.295]))
    fig.colorbar(axes[1,0].images[0], cax=fig.add_axes([.952,.04,.009,.295]))
    save(fig, figures, name)


def row_groups(rows):
    return {arm: [r for r in rows if r.get("arm") == arm] for arm in ARMS
            if any(r.get("arm") == arm for r in rows)}


def comparison_figures(campaign, rows, figures):
    groups = row_groups(rows)
    if not groups:
        return {}
    summary = {}
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.5))
    for arm, records in groups.items():
        ks, maes, coverages = [], [], []
        for k in (1,2,4,8):
            entries = [r.get("by_k",{}).get(str(k)) for r in records]
            entries = [e for e in entries if isinstance(e,dict)]
            maes.append(mean([e.get("best",{}).get("normalized_mae") for e in entries]))
            coverages.append(mean([int(e["coverage"]) for e in entries if "coverage" in e]))
            ks.append(k)
        c = COLORS.get(arm)
        axes[0,0].plot(ks, [np.nan if x is None else x for x in maes], "o-", color=c, label=ARMS[arm])
        axes[0,1].plot(ks, [np.nan if x is None else x for x in coverages], "o-", color=c)
        valid = mean([r.get("distinct_archive") for r in records])
        structural = mean([structural_count(m["archive"]) for item in campaign
                           for m in item["methods"] if m["arm"] == arm])
        calls = mean([r.get("calls") for r in records])
        wall = mean([r.get("wall") for r in records])
        summary[arm] = dict(parents=len(records), best_of_k_mae=dict(zip(map(str,ks),maes)),
                            coverage_at_k=dict(zip(map(str,ks),coverages)),
                            mean_valid_materials=valid, mean_structural_worlds=structural,
                            mean_search_endpoint_calls=calls, mean_search_wall_seconds=wall,
                            selected_mae=mean([r.get("selected",{}).get("normalized_mae") for r in records if r.get("selected")]),
                            common_base_mae=mean([r.get("base",{}).get("normalized_mae") for r in records]))
    arms = list(groups); positions=np.arange(len(arms)); width=.36
    axes[1,0].bar(positions-width/2,[summary[a]["mean_valid_materials"] for a in arms],width,label="Distinct material fields")
    axes[1,0].bar(positions+width/2,[summary[a]["mean_structural_worlds"] for a in arms],width,label="Structural worlds")
    axes[1,0].set_xticks(positions,[a.replace("_","\n") for a in arms],fontsize=8)
    axes[1,0].legend(fontsize=8)
    axes[1,1].bar(positions,[summary[a]["mean_search_wall_seconds"] for a in arms],color=[COLORS[a] for a in arms])
    axes[1,1].set_xticks(positions,[a.replace("_","\n") for a in arms],fontsize=8)
    for i,a in enumerate(arms):
        w=summary[a]["mean_search_wall_seconds"]; calls=summary[a]["mean_search_endpoint_calls"]
        axes[1,1].text(i,w,f"{calls:g} full endpoints",ha="center",va="bottom",fontsize=8)
    for ax in axes.flat: ax.grid(axis="y",alpha=.18)
    axes[0,0].set(xlabel="K, nested online material selector",ylabel="Offline best-of-K normalized material MAE",xticks=[1,2,4,8])
    axes[0,0].legend(fontsize=7)
    axes[0,1].set(xlabel="K",ylabel="Offline coverage fraction",xticks=[1,2,4,8],ylim=(-.02,1.02))
    axes[1,0].set_ylabel("Mean full-valid archive count per parent")
    axes[1,1].set_ylabel("Actual search wall seconds per parent")
    fig.suptitle("Confirmation panel: quality, coverage, valid diversity and measured search cost\n"
                 "Shared setup/coarse reconstruction and historical teacher cost reported separately",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.93)); save(fig,figures,"comparison_bestk_worlds_cost")

    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.bar(positions,[summary[a]["selected_mae"] for a in arms],color=[COLORS[a] for a in arms])
    baseline=summary[arms[0]]["common_base_mae"]
    ax.axhline(baseline,color="black",linestyle="--",label=f"Common base MAE={baseline:.6f}")
    ax.set_xticks(positions,[a.replace("_","\n") for a in arms],fontsize=9)
    ax.set_ylabel("Truth-scored material MAE of online lowest-residual selection")
    ax.set_title("Online selection versus the common base: discovery is not oracle selection")
    ax.legend(fontsize=9);ax.grid(axis="y",alpha=.2)
    fig.tight_layout();save(fig,figures,"selected_mae_and_common_base")

    fig,ax=plt.subplots(figsize=(8.5,5.4))
    for arm,record in summary.items():
        for k in (1,2,4,8):
            x=record["best_of_k_mae"][str(k)]; y=record["coverage_at_k"][str(k)]
            if x is not None and y is not None:
                ax.scatter(x,y,color=COLORS[arm],s=45)
                ax.annotate(f"{arm}, K={k}",(x,y),xytext=(4,4),textcoords="offset points",fontsize=7)
    ax.set(xlabel="Offline best-of-K normalized material MAE (lower better)",ylabel="Coverage@K (higher better)",ylim=(-.03,1.05))
    ax.grid(alpha=.2); ax.set_title("Coverage–MAE frontier, frozen acceptance and material scale")
    fig.tight_layout();save(fig,figures,"coverage_mae_frontier")

    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    for arm in groups:
        histories=[m.get("history",[]) for c in campaign for m in c["methods"] if m["arm"]==arm]
        maximum=max((len(x) for x in histories),default=0)
        xx=list(range(1,maximum+1))
        for ax,key in zip(axes,("distinct_worlds","structural_worlds")):
            yy=[mean([h[t-1].get(key) for h in histories if len(h)>=t]) for t in xx]
            ax.plot(xx,yy,"o-",color=COLORS[arm],label=ARMS[arm],markersize=3)
            ax.set_xlabel("Charged new full-Mie search endpoints")
            ax.set_ylabel("Mean distinct material fields" if key=="distinct_worlds" else "Mean structural worlds")
            ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle("Sequential discovery curves: history is online, truth is not used")
    fig.tight_layout();save(fig,figures,"discovery_history")
    return summary


def fit_and_structure(case, method, observation, figures):
    selected=method.get("selected_indices",[])[:8]
    if not selected: return
    candidates=[method["archive"][i] for i in selected]
    matrix=np.array([[distance(a["p"],b["p"]) for b in candidates] for a in candidates])
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.3))
    image=axes[0].imshow(matrix,cmap="cividis",vmin=0)
    axes[0].set(title="Volume-weighted material-field distance",xlabel="Candidate",ylabel="Candidate")
    axes[0].set_xticks(range(len(candidates)),range(1,len(candidates)+1));axes[0].set_yticks(range(len(candidates)),range(1,len(candidates)+1))
    fig.colorbar(image,ax=axes[0],fraction=.046)
    for i,c in enumerate(candidates):
        p=c["p"]; axes[1].scatter(p[0],p[1],s=70)
        axes[1].annotate(str(i+1),(p[0],p[1]),xytext=(4,4),textcoords="offset points")
    axes[1].set(xlabel="Central cavity radius / R",ylabel="Annular cavity half-width / R",xlim=(-.01,.305),ylim=(-.005,.125),title="Physical geometry, not latent duplication")
    axes[1].grid(alpha=.2);fig.suptitle(f"{case['id']} | {method['arm']} | online-only selected demo")
    fig.tight_layout();save(fig,figures,"structural_differences")
    if observation is None:return
    y=np.asarray(observation["y"],float); sigma=np.asarray(observation["sigma"],float)
    if sigma.ndim==0:sigma=np.full(y.shape,float(sigma))
    fig,axes=plt.subplots(2,1,figsize=(11,6.2),sharex=True)
    axes[0].errorbar(np.arange(len(y)),y,yerr=sigma,color="black",fmt=".",alpha=.65,label="Same observed y ± 1σ")
    for i,c in enumerate(candidates):
        f=np.asarray(c.get("fitted_prediction",c.get("prediction",[])),float)
        if f.shape!=y.shape:continue
        axes[0].plot(f,alpha=.65,label=f"Candidate {i+1}, g={c.get('gain',1):.3f}")
        axes[1].plot((f-y)/sigma,alpha=.7)
    axes[0].set_ylabel("Packed scattering (Re then Im)");axes[0].legend(fontsize=7,ncol=3)
    axes[1].axhline(0,color="black",linewidth=.7);axes[1].set(xlabel="Measurement channel",ylabel="Whitened residual per channel")
    axes[0].set_title("Original full-Maxwell endpoints; one bounded shared real gain per candidate")
    for a in axes:a.grid(alpha=.2)
    fig.tight_layout();save(fig,figures,"full_maxwell_scattering_fit")


def e2_saved_panel(results, figures):
    """Render only saved legal E2 candidates, with all failed cases retained."""
    folder=results/"E2"
    required=[folder/name for name in ("METRICS.json","SEARCH_HISTORY.json","MATERIAL_AND_PREDICTION_ARRAYS.npz","PROTOCOL.json")]
    if not all(p.exists() for p in required):return dict(files=[],available=False)
    rows=load(required[0]);history=load(required[1]);protocol=load(required[3])
    arrays=np.load(required[2],allow_pickle=False)
    cases=protocol["cases"];threshold=protocol["compatibility_radius"]
    finite_rows={r["case_id"]:r for r in rows if r["arm"]=="finite_history"}
    base_rows={r["case_id"]:r for r in history if r["arm"]=="base"}
    report=[]
    for case in cases:
        record=finite_rows[case];names=record.get("Top_8_selected_names",[])
        additions=[n for n in names if n!="base"]
        base=arrays[case+"__base"];background=arrays[case+"__background"]
        details=[]
        for name in additions:
            candidate=arrays[name+"__chi"]
            changed=int(np.sum((abs(base-background)>.18)!=(abs(candidate-background)>.18)))
            details.append(dict(name=name,mask_changed_cells=changed,new_structural_world=changed>=6,
                                material_rms=float(np.sqrt(np.mean(abs(candidate-base)**2)))))
        report.append(dict(case_id=case,base_compatible=bool(base_rows[case]["compatible"]),
                           valid_materials=record["valid_distinct_material_count"],
                           structural_worlds=record["structurally_distinct_masks_at_least6changedcells"],
                           additions=details,Best4=record["Best_of_4_complex_MAE"]))
    for group,start in enumerate(range(0,len(cases),4),1):
        part=cases[start:start+4]
        fig,axes=plt.subplots(len(part),6,figsize=(15,3.15*len(part)),squeeze=False)
        for row,case in enumerate(part):
            entry=report[cases.index(case)];addition=entry["additions"][:1]
            materials=[arrays[case+"__OFFLINE_TRUTH"],arrays[case+"__base"],
                       arrays[addition[0]["name"]+"__chi"] if addition else None]
            titles=["Truth: offline only", "Base: full-compatible" if entry["base_compatible"] else "Base: DATA-INCOMPATIBLE",
                    (f"Additional full-valid material\nmask change: {addition[0]['mask_changed_cells']} / required 6" if addition
                     else "NO ACCEPTED NEW CANDIDATE")]
            for col,material in enumerate(materials):
                for component in range(2):
                    ax=axes[row,col*2+component]
                    ax.set_xticks([]);ax.set_yticks([])
                    if material is None:
                        ax.set_facecolor("#f2f2f2");ax.text(.5,.5,"NO ACCEPTED\nNEW CANDIDATE",ha="center",va="center",fontsize=9,transform=ax.transAxes,color="#a23e37")
                    else:
                        values=(1+material.real if component==0 else material.imag).reshape(24,24)
                        ax.imshow(values,origin="lower",vmin=.5 if component==0 else 0.,vmax=5. if component==0 else .5,
                                  cmap="viridis" if component==0 else "magma",interpolation="nearest")
                    ax.set_title(titles[col]+("\nRe ε=1+Re χ" if component==0 else "\nIm ε=Im χ"),fontsize=8)
                    if col==0 and component==0:ax.set_ylabel(case.replace("_a1_snr20",""),fontsize=8)
        fig.suptitle("E2: 2D scalar TM/VIE, 24×24; fixed original case order, not 3D evidence\n"
                     "Shared declared material colors: Re ε=[0.5,5], Im ε=[0,0.5]; failed bases are labelled",fontsize=11)
        fig.tight_layout(rect=(0,0,1,.945));save(fig,figures,f"e2_material_cases_{group}")
    arms=list(dict.fromkeys(r["arm"] for r in rows))
    fig,axes=plt.subplots(1,2,figsize=(12,4.4));xx=np.arange(len(cases));width=.19
    for i,arm in enumerate(arms):
        records=[next(r for r in rows if r['case_id']==c and r['arm']==arm) for c in cases]
        axes[0].bar(xx+(i-(len(arms)-1)/2)*width,[r["valid_distinct_material_count"] for r in records],width,label=arm)
    axes[0].set_xticks(xx,[c.replace("_a1_snr20","").replace("p2500","") for c in cases],rotation=45,ha="right",fontsize=8)
    axes[0].set_ylabel("Actual full-compatible distinct materials (includes base)")
    axes[0].set_title("All eight conditions retained; zero means no valid solution")
    axes[0].legend(fontsize=7)
    values=[];eligible=[]
    for arm in arms:
        scores=[r["Best_of_4_complex_MAE"] for r in rows if r["arm"]==arm and r["Best_of_4_complex_MAE"] is not None]
        values.append(mean(scores));eligible.append(len(scores))
    axes[1].bar(np.arange(len(arms)),values,color=["#2075b0","#d87c24","#3a8d63","#9660a4"])
    axes[1].set_xticks(np.arange(len(arms)),[a.replace("_","\n") for a in arms],fontsize=8)
    for i,(v,n) in enumerate(zip(values,eligible)):axes[1].text(i,v,f"{v:.6f}\n{n}/8 eligible",ha="center",va="bottom",fontsize=8)
    axes[1].set_ylabel("Best-of-4 complex material MAE, eligible cases only")
    axes[1].set_title("Equal eligible-case means; five failed cases are not averaged away")
    axes[1].set_ylim(0,max(values)*1.22)
    for ax in axes:ax.grid(axis="y",alpha=.18)
    fig.suptitle("Hard counterexample: compensation does not improve E2 material coverage or Best-of-4",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.92));save(fig,figures,"e2_negative_comparison")
    files=[str(p.relative_to(results.parent)) for p in required]
    result=dict(available=True,files=files,cases=report,
                no_valid_solution_cases=sum(r["valid_materials"]==0 for r in report),
                total_cases=len(cases),base_valid_cases=sum(r["base_compatible"] for r in report),
                color_ranges=dict(real_epsilon=[.5,5.],imag_epsilon=[0.,.5]),
                structural_mask="abs(chi-background)>.18; at least 6 changed cells",
                matched_eligible_best4=dict(zip(arms,values)),posterior_probability=None,
                scope="2D scalar TM/VIE, historically exposed development; separate from radial3D evidence")
    (results/"E2_GALLERY_SEMANTICS.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    return result


def make_html(root, campaign, truth_map, summary, figures, selection, e2):
    gallery=[]
    for c in campaign:
        methods=[]
        for m in c["methods"]:
            indices=m.get("selected_indices",[])[:8]
            methods.append(dict(arm=m["arm"],label=ARMS.get(m["arm"],m["arm"]),
                                full_valid_count=len(m["archive"]),structural_count=structural_count(m["archive"]),
                                calls=m.get("endpoint_calls"),wall=m.get("wall_seconds"),
                                candidates=[dict(archive_index=i,**m["archive"][i]) for i in indices]))
        gallery.append(dict(id=c["id"],base=c.get("base"),truth=truth_map.get(c["id"]),methods=methods))
    data=dict(cases=gallery,summary=summary,selection=selection,numerical_status="candidate_only",
              posterior_probability=None,dimension_scope="3D vector-Mie concentric radial support; center sections only",
              unexplored_structures="unquantified",outside_model_support_mass=None,E2_saved_results=e2.get("files",[]),E2=e2)
    (root/"GALLERY_DATA.json").write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    script=json.dumps(data,ensure_ascii=False,allow_nan=False).replace("</","<\\/")
    source='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research A｜下一种材料解释</title>
<style>body{margin:0;background:#f3f5f7;color:#182736;font:16px/1.6 system-ui,sans-serif}main{max-width:1440px;margin:auto;padding:28px}header{padding:24px;background:#132e40;color:white;border-radius:14px}h1{font-size:30px;line-height:1.3}h2{font-size:23px}.note{background:#fff3d8;padding:14px;border-radius:9px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(235px,1fr));gap:14px}.card{padding:14px;border-radius:10px;background:white;border:1px solid #dbe3e9}.pair{display:grid;grid-template-columns:1fr 1fr;gap:5px}canvas{width:100%;aspect-ratio:1}select{font:inherit;padding:8px;margin-right:12px;max-width:98%}.muted{color:#5c6b77;font-size:14px}img.plot{width:100%;background:white;border-radius:10px;margin:8px 0}a{color:#116ba0}.bar{height:13px;border-radius:3px}.real{background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)}.imag{background:linear-gradient(90deg,#000004,#51127c,#b73779,#fc8961,#fcfdbf)}table{border-collapse:collapse;background:white;width:100%}td,th{padding:9px;border:1px solid #dde5ea;text-align:left}.overflow{overflow:auto}small{font-size:12px}</style>
<main><header><div>Research A · frozen confirmation experiment</div><h1>同一份散射测量，继续寻找下一种材料世界</h1><p>有限结构编辑 → 有界非线性补偿 → 原始 full-Maxwell endpoint → 材料空间去重</p><p>这里展示真实保存的候选；概率不是研究主目标，全部 posterior probability = null。</p></header>
<p class="note">范围：三维矢量 Maxwell 的<strong>同心径向分层球</strong>，图像是中心截面，不代表一般三维未知形状。候选均使用原始 full-Mie 检验；合法材料数量、显著结构数量与后验模式数量是不同概念。Best-of-K 与 truth 仅用于离线评分。</p>
<p><a href="START_HERE.md">阅读入口</a> · <a href="RESEARCH_REPORT_ZH.md">中文研究报告</a> · <a href="results/METRICS.json">原始评分</a> · <a href="GALLERY_DATA.json">机器可读展示数据</a> · <a href="COST_LEDGER.jsonl">成本账本</a></p>
<h2>完整确认面板的实测结果</h2><p>Best-of-4 是离线真值评分。Selected 是在线只按 full 残差选出的单张材料；两者不能混同。候选数含共同基础解。</p><div id="summary" class="overflow"></div>
<h2>材料候选</h2><p>默认展示预注册的第一 confirmation parent。其他案例和方法保持同一材料色标。</p>
<p><label>观测 <select id="case"></select></label><label>方法 <select id="arm"></select></label></p><div id="facts" class="note"></div>
<div class="grid" style="margin:14px 0"><div>实部 Re ε：1.0（真空）→ 2.8<div class="bar real"></div></div><div>损耗 Im ε：0 → 0.16<div class="bar imag"></div></div></div><div id="cards" class="grid"></div>
<h2>Coverage–MAE–Cost 与发现历史</h2><p>每个方法使用相同接受条件和材料距离。图中成本来自本轮实际保存的搜索费用；共同 setup、基础重建、历史 bank/teacher 等费用另外记录，不能据此宣称端到端加速。</p>
<img class="plot" src="figures/comparison_bestk_worlds_cost.png" alt="Best-of-K、coverage、材料与结构数量、实际搜索成本"><img class="plot" src="figures/selected_mae_and_common_base.png" alt="在线选中材料并未改善共同基础解的MAE"><img class="plot" src="figures/coverage_mae_frontier.png" alt="MAE与coverage前沿"><img class="plot" src="figures/discovery_history.png" alt="随收费fullendpoint增长的发现曲线">
<h2>在线规则选出的结构展示案例</h2><p id="demo"></p><img class="plot" src="figures/online_selected_demo.png" alt="按在线结构数量与残差选择的材料案例"><img class="plot" src="figures/structural_differences.png" alt="材料距离和几何差异"><img class="plot" src="figures/full_maxwell_scattering_fit.png" alt="fullMaxwell散射拟合和白化残差">
<h2>二维困难条件</h2><div id="e2"></div><h2>解释边界</h2><p>同一观测下多个合法且材料上不同的解，不能自动证明 posterior calibration、全局唯一性结论或重要解族全部发现。未探索结构状态未量化；模型外概率质量为 null。候选 online 排序不访问 truth。没有用 residual softmax 生成虚假概率。</p><p>原始失败、重复候选、保护结构丢失和费用应与成功一起阅读。科学结论以父 Codex 最终研究报告为准。</p></main>
<script id="data" type="application/json">__DATA__</script><script>
const data=JSON.parse(document.getElementById('data').textContent),cs=document.getElementById('case'),as=document.getElementById('arm');
const esc=x=>String(x).replace(/[&<>\"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;'}[c]));
for(const c of data.cases){let o=document.createElement('option');o.value=c.id;o.textContent=c.id;cs.appendChild(o)}
const palettes={real:[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]],imag:[[0,0,4],[81,18,124],[183,55,121],[252,137,97],[252,253,191]]};
function color(v,part){let t=Math.max(0,Math.min(1,part==='real'?(v-1)/1.8:v/.16)),p=palettes[part],q=t*(p.length-1),i=Math.min(p.length-2,Math.floor(q)),f=q-i;return p[i].map((x,j)=>Math.round(x+(p[i+1][j]-x)*f))}
function draw(canvas,p,part){let n=193;canvas.width=canvas.height=n;let ctx=canvas.getContext('2d'),im=ctx.createImageData(n,n);for(let y=0;y<n;y++)for(let x=0;x<n;x++){let r=Math.hypot((x/(n-1)*2-1)*1.04,(y/(n-1)*2-1)*1.04),vac=r<p[0]||(p[1]>0&&Math.abs(r-.65)<p[1]),v=part==='real'?(vac?1:p[2]):(vac?0:p[3]),c=r>1?[255,255,255]:color(v,part),j=4*(y*n+x);im.data.set([...c,255],j)}ctx.putImageData(im,0,0)}
function card(label,c,truth=false){if(!c||!c.p)return;let d=document.createElement('div');d.className='card';let p=c.p;d.innerHTML='<strong>'+esc(label)+'</strong><div class="pair"><div><small>Re ε</small><canvas></canvas></div><div><small>Im ε</small><canvas></canvas></div></div><div class="muted">r/R='+p[0].toFixed(4)+'；h/R='+p[1].toFixed(4)+'<br>ε='+p[2].toFixed(4)+'+'+p[3].toFixed(4)+'i'+(c.residual!==undefined?'<br>full 白化残差 '+c.residual.toFixed(4):'')+(c.gain!==undefined?'<br>shared gain '+c.gain.toFixed(4):'')+'<br>'+(truth?'仅离线展示，不参与选择':'candidate_only；probability = null')+'</div>';document.getElementById('cards').appendChild(d);draw(d.querySelectorAll('canvas')[0],p,'real');draw(d.querySelectorAll('canvas')[1],p,'imag')}
function fillArms(){let c=data.cases.find(c=>c.id===cs.value),old=as.value;as.innerHTML='';for(let m of c.methods){let o=document.createElement('option');o.value=m.arm;o.textContent=m.label;as.appendChild(o)}if(c.methods.some(m=>m.arm===old))as.value=old;render()}
function render(){let c=data.cases.find(c=>c.id===cs.value),m=c.methods.find(m=>m.arm===as.value);document.getElementById('cards').innerHTML='';document.getElementById('facts').textContent='观测 '+c.id+' · full-valid 材料数 '+m.full_valid_count+' · 显著结构数 '+m.structural_count+' · 新搜索 full endpoints '+m.calls+' · 实测搜索 wall '+Number(m.wall).toFixed(3)+' s。Top-K 使用在线材料距离嵌套选择；不是 posterior mass。';card('真实材料（离线）',c.truth,true);card('共同基础重建',c.base);m.candidates.forEach((a,i)=>card('候选 '+(i+1)+' / archive '+a.archive_index,a))}
cs.onchange=fillArms;as.onchange=render;fillArms();document.getElementById('demo').textContent=data.selection.description;
const num=x=>x===null?'NOT_RUN':Number(x).toFixed(5);document.getElementById('summary').innerHTML='<table><tr><th>方法</th><th>Best-of-4 MAE ↓</th><th>Coverage@4 ↑</th><th>材料数 / 结构数</th><th>Selected MAE ↓</th><th>共同 base MAE</th><th>搜索 wall 秒/parent</th></tr>'+Object.entries(data.summary).map(([a,s])=>'<tr><td>'+esc(a)+'</td><td>'+num(s.best_of_k_mae['4'])+'</td><td>'+num(s.coverage_at_k['4'])+'</td><td>'+num(s.mean_valid_materials)+' / '+num(s.mean_structural_worlds)+'</td><td>'+num(s.selected_mae)+'</td><td>'+num(s.common_base_mae)+'</td><td>'+num(s.mean_search_wall_seconds)+'</td></tr>').join('')+'</table>';
document.getElementById('e2').innerHTML=data.E2.available?'<p class="note"><strong>E2 是有效负结果：'+data.E2.no_valid_solution_cases+'/'+data.E2.total_cases+' 条件没有合法解，只有 '+data.E2.base_valid_cases+' 条件的基础解通过。</strong>四方法 Best-of-4 在3个有效条件上完全相同。唯一额外材料仅改变1个结构 mask cell，未达到≥6 cell标准，并且同样出现在不补偿对照中；不构成新结构族或补偿收益。不能把5个失败条件从结论中删去。</p><p>独立二维标量 TM/VIE，24×24，原始噪声和gain=1；不是三维矢量成像。所有图保持固定顺序和统一色标，未通过的数据相容性基础解明确标注。</p><img class="plot" src="figures/e2_negative_comparison.png" alt="E2八条件失败和四方法相同Best4"><img class="plot" src="figures/e2_material_cases_1.png" alt="E2前四固定条件的真实材料基础解和合法新增材料"><img class="plot" src="figures/e2_material_cases_2.png" alt="E2后四固定条件，无额外合法材料时明确留空"><ul>'+data.E2_saved_results.map(x=>'<li><a href="'+esc(x)+'">'+esc(x)+'</a></li>').join('')+'</ul>':'<p>当前 gallery 输入中尚无 E2 实测结果。此处不建立虚构二维成绩；是否运行及原因见最终报告。</p>';
</script></html>'''
    (root/"index.html").write_text(source.replace("__DATA__",script))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--e2-only",action="store_true",help="Append saved E2 panel without redrawing completed E1 figures")
    args=parser.parse_args();root=args.root.resolve();results=root/"results"
    campaign=load(results/"CONFIRM_CAMPAIGN.json")
    if not campaign:
        parser.exit(2,"No completed CONFIRM_CAMPAIGN.json; no gallery or claims generated.\n")
    if isinstance(campaign,dict):campaign=campaign.get("cases",campaign.get("campaign",[]))
    if not campaign:parser.exit(2,"Confirmation campaign contains no cases.\n")
    if not (results/"METRICS.json").exists():
        parser.exit(2,"No final METRICS.json; campaign may still be running. No gallery generated.\n")
    truth=load(results/"offline_truth/CONFIRM_TRUTH.json",load(root/"offline_truth/CONFIRM_TRUTH.json",[]))
    truths={x["id"]:x for x in truth}
    rows=load(results/"METRICS.json",[])
    if isinstance(rows,dict):rows=rows.get("rows",rows.get("offline_rows",rows.get("metrics",[])))
    if not isinstance(rows,list):raise ValueError("METRICS.json must contain offline row list")
    observations=load(results/"CONFIRM_OBSERVATIONS.json",[])
    obs={x["id"]:x for x in observations}
    figures=root/"figures";figures.mkdir(exist_ok=True)
    first=campaign[0]
    if not args.e2_only:
        for method in first["methods"]:
            topk_figure(first,method,truths.get(first["id"]),figures,"first_confirm_"+method["arm"])
    ranked=[]
    for ci,case in enumerate(campaign):
        for mi,method in enumerate(case["methods"]):
            archive=method["archive"]
            if archive:
                ranked.append((-structural_count(archive),min(a["residual"] for a in archive),ci,mi))
    if not ranked:raise ValueError("No physically valid archived candidates; no positive demo fabricated")
    _,_,ci,mi=min(ranked);demo=campaign[ci];method=demo["methods"][mi]
    selection=dict(case=demo["id"],arm=method["arm"],uses_truth=False,
                   rule="max online structural-world count, then lowest archived full residual, then campaign/method order",
                   description=f"{demo['id']} / {method['arm']}。预声明在线规则：最大显著结构数→最低 full 残差→原始顺序；未使用 truth、MAE 或 posterior mass 选择。")
    if args.e2_only:
        summary=load(results/"GALLERY_PROVENANCE.json")["summary"]
    else:
        topk_figure(demo,method,truths.get(demo["id"]),figures,"online_selected_demo")
        summary=comparison_figures(campaign,rows,figures)
        fit_and_structure(demo,method,obs.get(demo["id"]),figures)
    e2=e2_saved_panel(results,figures)
    make_html(root,campaign,truths,summary,figures,selection,e2)
    (results/"GALLERY_PROVENANCE.json").write_text(json.dumps(dict(
        source="saved arrays only; no physics imports, no generated experiments",first_confirm=first["id"],
        demonstration=selection,summary=summary,probability_semantics="candidate_only; posterior probabilities null",
        new_physics_calls=0,color_ranges=dict(real=[1.,2.8],imag=[0.,.16]),
        E2_files=e2),ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(status="RENDERED_SAVED_EVIDENCE",root=str(root),parents=len(campaign),
                          first=first["id"],demo=selection,new_physics_calls=0),ensure_ascii=False))


if __name__=="__main__":main()
