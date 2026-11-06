"""Seal all66 confirmation tasks, then score/plot offline. No work on import."""
from pathlib import Path
import argparse
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r9.experiments import scoring_common as common

STATUS="ALL66_CONFIRMATION_PUBLIC_OUTPUTS_SEALED"
SCORED="SCORED_ALL66_CONFIRMATION_TASKS"
METHODS=(("cie_single","fixed_001"),("cie_dual","fixed_001"),
         ("global_dual","fixed_001"),("a_plus","fixed_001"),
         ("cie_single","noise_ball_99"),("cie_dual","noise_ball_99"),("global_dual","noise_ball_99"),
         ("a_plus","noise_ball_99"),("physical_lbfgs","fixed_001"),
         ("profile_homotopy",None),("weak_profile_handoff",None))
LABELS=("Fixed single CIE","Fixed Dual","Fixed GlobalDual","Fixed global A+",
        "Noise99 single CIE","Noise99 Dual","Noise99 GlobalDual","Noise99 global A+","Direct physics","Profile continuation","Weak-profile handoff")


def inputs(args):
    public=common.resolve(args.public_manifest)
    cases=[r["case_id"] for r in common.load(public)["records"]]
    manifests=[ROOT/"configs/confirmation_baselines.json",ROOT/"configs/confirmation_landscape.json",
               ROOT/"configs/confirmation_single_noise.json",
               ROOT/"configs/confirmation_weak_handoff.json"]
    tasks=common.tasks_from(manifests,cases)
    if len(cases)!=6 or len(set(cases))!=6 or len(tasks)!=66:
        raise ValueError("Frozen six-case66-task confirmation required")
    for ident in cases:
        keys=[(t.get("algorithm",t.get("arm")),t.get("dpcs_rule")) for t in tasks if t["case_id"]==ident]
        if len(keys)!=11 or set(keys)!=set(METHODS):
            raise ValueError("Frozen eleven-arm coverage required for every case")
    return public,cases,manifests


def seal(args):
    public,cases,manifests=inputs(args)
    extra=[public,Path(__file__)]
    for task in common.tasks_from(manifests,cases):
        if task.get("algorithm")=="weak_profile_handoff":
            extra.extend(common.resolve(task[key]) for key in
                         ("warmup_receipt","warmup_budget","reused_initialization_provenance"))
    sealed=common.make_seal(common.resolve(args.seal),manifests,cases,common.resolve(args.r8_repo),STATUS,
        extra_evidence=extra,portfolio_rule=ROOT/"configs/portfolio_rule.json")
    print(dict(status=sealed["status"],tasks=sealed["requested_tasks"],truth_access="NONE"))


def truths_after_seal(cases):
    base=ROOT/"sealed_truth/confirmation"
    meta={r["case_id"]:r for r in common.load(base/"score_metadata.json")["records"]}
    result={}
    for ident in cases:
        path=base/(ident+".npz")
        with np.load(path,allow_pickle=False) as z:
            result[ident]=dict(m=z["m"].ravel().copy(),mask=z["mask"].ravel().astype(bool),
                amplitude=float(meta[ident]["amplitude"]),same_grid=True,truth_path=str(path))
    return result


def weak_prefix(task):
    """Public paid first-stage prefix; excludes later unneeded detail evaluation."""
    receipt=common.load(task["warmup_receipt"])
    stages=[r for r in common.load(task["warmup_budget"])
            if r["label"]=="profile_0_0.01_endpoint"]
    if len(stages)!=1:raise ValueError("Unique paid weak first-stage endpoint required")
    stage=stages[0];startup=receipt.get("pre_solver_CPU_seconds")
    semantics="EXACT_RECORDED_PRE_SOLVER_PREFIX"
    if startup is None:
        startup=max(0.,receipt["process_CPU_seconds"]-receipt["cost"]["CPU_seconds"])
        semantics="UPPER_SHIFT_including_output_tail"
    with np.load(task["initialization_file"],allow_pickle=False) as z:initial=z["m_init"].ravel()
    with np.load(task["reused_initialization_provenance"],allow_pickle=False) as z:source=z["m"].ravel()
    if not np.array_equal(initial,source):raise ValueError("Weak initializer differs from declared paid stage")
    return dict(CPU_seconds=stage["CPU_seconds"]+startup,CPU_semantics=semantics,
        startup_CPU_seconds=startup,raw_stage_prefix_CPU_seconds=stage["CPU_seconds"],
        source_receipt=task["warmup_receipt"],source_budget=task["warmup_budget"],
        source_label=stage["label"],physical=stage["physical_cost"],auxiliary=stage["auxiliary_cost"],
        count_scope="At paid endpoint, before poststage details/normal-residual extra current factorization")


def add_weak_costs(rows,budgets,prefixes):
    for row in rows:
        prefix=prefixes.get(row["task_id"])
        if prefix is None:continue
        row.update(reused_initialization_cost=prefix,
            standalone_CPU_seconds=row["cpu_seconds"]+prefix["CPU_seconds"] if row["cpu_seconds"] is not None else None,
            standalone_LU=row["LU"]+prefix["physical"].get("LU",0) if row["LU"] is not None else None,
            standalone_complex_RHS=row["complex_RHS"]+prefix["physical"].get("actual_complex_RHS",0) if row["complex_RHS"] is not None else None,
            standalone_current_Cholesky=row["current_Cholesky"]+prefix["auxiliary"].get("current_Cholesky",0),
            standalone_current_RHS=row["current_RHS"]+prefix["auxiliary"].get("current_RHS",0),
            cost_scope="cpu_seconds/counts are actual new continuation only; standalone_* adds reused first-stage cost. Do not doublecharge prefix in actual66 ledger.")
    for row in budgets:
        prefix=prefixes.get(row["task_id"])
        if prefix is None:continue
        row.update(standalone_CPU_seconds=row["CPU_seconds"]+prefix["CPU_seconds"] if row["CPU_seconds"] is not None else None,
            standalone_LU=row["LU"]+prefix["physical"].get("LU",0) if row["LU"] is not None else None,
            standalone_complex_RHS=row["complex_RHS"]+prefix["physical"].get("actual_complex_RHS",0) if row["complex_RHS"] is not None else None,
            standalone_prefix_semantics=prefix["CPU_semantics"])


def score(args):
    begin,cpu=time.perf_counter(),time.process_time()
    _,cases,manifests=inputs(args)
    sealed=common.verify(common.resolve(args.seal),manifests,cases,66,STATUS)
    out=common.resolve(args.output)
    if out.exists():raise FileExistsError(out)
    prefixes={e["task"]["task_id"]:weak_prefix(e["task"]) for e in sealed["entries"]
              if e["task"].get("algorithm")=="weak_profile_handoff" and e["status"] in common.dev.SUCCESS}
    # First confirmation private reads occur strictly below complete seal checks.
    truths=truths_after_seal(cases)
    rows,stages,budgets,materials,metric_path=common.score_entries(sealed,truths,"FROZEN_SIX_CASE_CONFIRMATION")
    add_weak_costs(rows,budgets,prefixes)
    rows+=common.add_portfolios(sealed,rows,materials)
    common.write_bundle(out,rows,stages,budgets,materials,truths,sealed,metric_path,SCORED,(cpu,begin))
    print(dict(status=SCORED,physical_tasks=66,portfolio_rows=len(sealed["portfolios"]),new_physical_calls=0))


def figures(args):
    _,cases,manifests=inputs(args)
    common.verify(common.resolve(args.seal),manifests,cases,66,STATUS)
    out=common.resolve(args.output)
    if common.load(out/"RECEIPT.json")["status"]!=SCORED:
        raise ValueError("Complete authorized offline scoring required before figures")
    rows=common.load(out/"SCORES.json")["rows"]
    index=common.load(out/"MATERIAL_INDEX.json")
    with np.load(out/"MATERIALS.npz",allow_pickle=False) as z:
        material={tid:m for tid,m in zip(index["task_ids"],z["prediction"])}
        truth={cid:m for cid,m in zip(index["case_ids"],z["truth"])}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder=ROOT/"FIGURES/confirmation";folder.mkdir(parents=True,exist_ok=True)
    table={(r["case_id"],r["method"],r.get("dpcs_rule")):r for r in rows}
    fig,axes=plt.subplots(6,12,figsize=(38,18),constrained_layout=True)
    for i,cid in enumerate(cases):
        selected=[table[(cid,)+key] for key in METHODS]
        mats=[truth[cid]]+[material[r["task_id"]] for r in selected if r["task_id"] in material]
        vmax=max(float(np.max(m)) for m in mats);vmin=min(.02,min(float(np.min(m)) for m in mats))
        im=axes[i,0].imshow(truth[cid].reshape(24,24),origin="lower",vmin=vmin,vmax=vmax,cmap="viridis")
        axes[i,0].set_title(cid+" GT (offline)",fontsize=8)
        for j,(label,row) in enumerate(zip(LABELS,selected),1):
            ax=axes[i,j]
            if row["task_id"] in material:
                ax.imshow(material[row["task_id"]].reshape(24,24),origin="lower",vmin=vmin,vmax=vmax,cmap="viridis")
                display_cpu=row.get("standalone_CPU_seconds",row["cpu_seconds"])
                suffix="*" if "standalone_CPU_seconds" in row else ""
                ax.set_title(f"{label}\nN={row['material_nrmse']:.3f} IoU={row['otsu_iou']:.2f}\nCPU{suffix}={display_cpu:.1f}s",fontsize=8)
            else:
                ax.set_facecolor("lightgray");ax.text(.5,.5,"FAILED",ha="center",va="center",transform=ax.transAxes)
                ax.set_title(label,fontsize=8)
        for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
        fig.colorbar(im,ax=list(axes[i]),shrink=.65,label="Full common row material range")
    fig.suptitle("All six frozen confirmation objects | eleven arms | full amplitudes | weak CPU* includes reused first-stage cost")
    fig.savefig(folder/"all_six_materials.png",dpi=160);plt.close(fig)
    # Two frozen portfolios are shown alongside all actual arms, with both costs.
    extra=[(r["name"],"two_frozen_candidates") for r in common.load(ROOT/"configs/portfolio_rule.json")["portfolios"]]
    keys=list(METHODS)+extra;labels=list(LABELS)+[k[0] for k in extra]
    for row in rows:row["method_CPU_seconds"]=row.get("standalone_CPU_seconds",row.get("cpu_seconds"))
    fields=("material_nrmse","otsu_iou","method_CPU_seconds")
    fig,axes=plt.subplots(1,3,figsize=(22,9),constrained_layout=True)
    for ax,field in zip(axes,fields):
        a=np.array([[table[(cid,)+key].get(field) if table[(cid,)+key].get(field) is not None else np.nan for cid in cases] for key in keys],float)
        im=ax.imshow(a,aspect="auto",cmap="viridis_r" if field!="otsu_iou" else "viridis")
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                val=a[i,j];txt="FAILED" if not np.isfinite(val) else f"{val:.1f}" if field=="method_CPU_seconds" else f"{val:.3f}"
                ax.text(j,i,txt,ha="center",va="center",fontsize=7,color="black",bbox=dict(facecolor="white",alpha=.65,edgecolor="none",pad=.5))
        ax.set_xticks(range(6),[c[-4:] for c in cases]);ax.set_yticks(range(len(labels)),labels,fontsize=8)
        ax.set_title(field);fig.colorbar(im,ax=ax,shrink=.7)
    fig.suptitle("All six cases | weak handoff CPU includes init | unchanged two portfolios pay both candidates")
    fig.savefig(folder/"paired_cost_matrix.png",dpi=160);plt.close(fig)
    common.dump(folder/"FIGURE_RECEIPT.json",dict(status="GENERATED_FROM_SEALED_OFFLINE_SCORES",cases=cases,
        material_panels=72,actual_arms=11,portfolio_matrix_rows=2,full_common_row_range=True,
        weak_CPU_includes_standalone_initialization=True,actual_ledger_never_doublecharges_reuse=True,
        new_physical_calls=0,truth_source="Postseal MATERIALS.npz only; no new private reads"))
    print(dict(status="FIGURES_WRITTEN",folder=str(folder)))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("mode",choices=("seal","score","figures"))
    ap.add_argument("--r8-repo",default=str(REPO.parent/"r8_worktree"))
    ap.add_argument("--public-manifest",default=str(ROOT/"configs/confirmation_public_manifest.json"))
    ap.add_argument("--seal",default=str(ROOT/"audit/CONFIRMATION_PUBLIC_SEAL.json"))
    ap.add_argument("--output",default=str(ROOT/"RESULTS/confirmation_scoring"))
    args=ap.parse_args();{"seal":seal,"score":score,"figures":figures}[args.mode](args)


if __name__=="__main__":main()
