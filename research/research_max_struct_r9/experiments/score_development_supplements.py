"""Separate postdevelopment three-task follow-up scoring; original43 untouched."""
from pathlib import Path
import argparse
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r9.experiments import scoring_common as common

CASES=("obs_b20d5a55677d4eca","obs_40ac5b79e8a74ac3")
STATUS="ALL3_DEVELOPMENT_SUPPLEMENT_PUBLIC_OUTPUTS_SEALED"
MANIFESTS=(ROOT/"configs/profile_restart_control.json",ROOT/"configs/curvature_followup.json")
DIAGNOSTIC=ROOT/"RESULTS/curvature/obs_b20d5a55677d4eca__global_dual/RECEIPT.json"


def inputs():
    tasks=common.tasks_from(MANIFESTS,CASES)
    expected={CASES[0]+"__curvature_anchor",CASES[0]+"__curvature_proposal",CASES[1]+"__profile_restarted"}
    if {t["task_id"] for t in tasks}!=expected or len(tasks)!=3:
        raise ValueError("Only frozen restart/anchor/proposal three-task follow-up")
    return tasks


def seal(args):
    inputs()
    result=common.make_seal(common.resolve(args.seal),MANIFESTS,CASES,common.resolve(args.r8_repo),STATUS,
        extra_evidence=(Path(__file__),DIAGNOSTIC,DIAGNOSTIC.parent/"curvature.npz"))
    print(dict(status=result["status"],tasks=3,truth_access="NONE"))


def truth_after_seal(r8):
    historical=common.load(r8/"research/research_max_struct_r6/sealed_truth/scoring/historical_mapping_PRIVATE.json")
    historical={r["case_id"]:r for r in historical}
    confirmation=common.load(r8/"research/research_max_struct_r6/sealed_truth/confirmation/mapping_PRIVATE.json")
    confirmation={r["case_id"]:r for r in confirmation["records"]}
    public=common.load(r8/"research/research_max_struct_r6/data_public/cases.json")
    public={r["case_id"]:r for r in public.get("records",public.get("cases",[]))}
    result={}
    for ident in CASES:
        meta=historical[ident] if ident==CASES[0] else confirmation[ident]
        p=common.resolve(meta["truth_path"],r8)
        with np.load(p,allow_pickle=False) as z:
            result[ident]=dict(m=z["m"].ravel().copy(),mask=z["mask"].ravel().astype(bool),
                amplitude=float(public[ident]["amplitude"]),same_grid=meta.get("inverse_crime"),truth_path=str(p))
    return result


def score(args):
    begin,cpu=time.perf_counter(),time.process_time();inputs()
    sealed=common.verify(common.resolve(args.seal),MANIFESTS,CASES,3,STATUS)
    out=common.resolve(args.output)
    if out.exists():raise FileExistsError(out)
    truths=truth_after_seal(Path(sealed["r8_repo"]))
    rows,stages,budgets,materials,metric_path=common.score_entries(sealed,truths,"POSTHOC_DEVELOPMENT_SUPPLEMENT_ONLY")
    diagnostic=common.load(DIAGNOSTIC);dcost=diagnostic["cost"];physical=dcost["physical"]
    proposal=CASES[0]+"__curvature_proposal"
    for row in rows:
        row["upstream_endpoint_scope"]="Incremental continuation from already paid exposed development endpoint; upstream original43 run not charged again here"
        row["solver_only_CPU_seconds"]=row["cpu_seconds"]
        if row["task_id"]==proposal:
            row["charged_diagnostic"]=dict(receipt_path=str(DIAGNOSTIC),cost=dcost,
                scope="Whole Hessian/eigendecomposition/probe work paid once on proposal; anchor has no such prefix")
            row["cpu_seconds"]+=dcost["CPU_seconds"]
            row["wall_seconds"]+=dcost["wall_seconds"]
            row["LU"]+=physical["LU"];row["complex_RHS"]+=physical["actual_complex_RHS"]
            row["physical_model_cost"]={k:row["physical_model_cost"].get(k,0)+physical.get(k,0)
                for k in row["physical_model_cost"].keys()|physical.keys()}
    for row in budgets:
        if row["task_id"]==proposal:
            row["solver_only_CPU_seconds"]=row["CPU_seconds"]
            if row["CPU_seconds"] is not None:row["CPU_seconds"]+=dcost["CPU_seconds"]
            row["CPU_semantics"]+=" + CHARGED_HESSIAN_AND_PROBES_PREFIX"
            if row["LU"] is not None:row["LU"]+=physical["LU"]
            if row["complex_RHS"] is not None:row["complex_RHS"]+=physical["actual_complex_RHS"]
    common.write_bundle(out,rows,stages,budgets,materials,truths,sealed,metric_path,
        "SCORED_ALL3_DEVELOPMENT_SUPPLEMENTS",(cpu,begin))
    print(dict(status="SCORED_ALL3_DEVELOPMENT_SUPPLEMENTS",new_physical_calls=0,original43_changed=False))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("mode",choices=("seal","score"))
    ap.add_argument("--r8-repo",default=str(REPO.parent/"r8_worktree"))
    ap.add_argument("--seal",default=str(ROOT/"audit/DEVELOPMENT_SUPPLEMENT_PUBLIC_SEAL.json"))
    ap.add_argument("--output",default=str(ROOT/"RESULTS/development_supplement_scoring"))
    args=ap.parse_args();(seal if args.mode=="seal" else score)(args)


if __name__=="__main__":main()
