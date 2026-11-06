"""Separate exposed-DEV weak-stage handoff grading; no old result rewrites."""
from pathlib import Path
import argparse
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r9.experiments import scoring_common as common

CASE="obs_40ac5b79e8a74ac3"
MANIFEST=ROOT/"configs/weak_handoff_control.json"
SOURCE=ROOT/"RESULTS/development"/(CASE+"__profile_homotopy")
STATUS="WEAK_HANDOFF_PUBLIC_PREDICTION_SEALED"


def inputs():
    tasks=common.tasks_from([MANIFEST],[CASE])
    if len(tasks)!=1 or tasks[0]["task_id"]!=CASE+"__weak_profile_handoff":
        raise ValueError("Only the frozen one-case weak handoff control")
    return tasks[0]


def seal(args):
    task=inputs()
    result=common.make_seal(common.resolve(args.seal),[MANIFEST],[CASE],common.resolve(args.r8_repo),STATUS,
        extra_evidence=(Path(__file__),SOURCE/"RECEIPT.json",SOURCE/"BUDGET.json",
                        common.resolve(task["reused_initialization_provenance"])))
    print(dict(status=result["status"],truth_access="NONE"))


def public_prefix(task):
    rec=common.load(SOURCE/"RECEIPT.json")
    endpoints=[r for r in common.load(SOURCE/"BUDGET.json") if r["label"]=="profile_0_0.01_endpoint"]
    if len(endpoints)!=1:raise ValueError("Unique paid first-stage endpoint required")
    endpoint=endpoints[0];startup=rec.get("pre_solver_CPU_seconds")
    kind="EXACT_RECORDED_PRE_SOLVER_PREFIX"
    if startup is None:
        startup=max(0.,rec["process_CPU_seconds"]-rec["cost"]["CPU_seconds"])
        kind="UPPER_SHIFT_including_output_tail"
    with np.load(task["initialization_file"],allow_pickle=False) as z:initial=z["m_init"].ravel()
    with np.load(task["reused_initialization_provenance"],allow_pickle=False) as z:source=z["m"].ravel()
    if not np.array_equal(initial,source):raise ValueError("Handoff initializer differs from declared paid first stage")
    return dict(source_receipt=str(SOURCE/"RECEIPT.json"),source_budget=str(SOURCE/"BUDGET.json"),
        CPU_seconds=endpoint["CPU_seconds"]+startup,CPU_semantics=kind,
        startup_CPU_seconds=startup,raw_stage_prefix_CPU_seconds=endpoint["CPU_seconds"],
        physical=endpoint["physical_cost"],auxiliary=endpoint["auxiliary_cost"],
        source_label=endpoint["label"],original_Phi=endpoint["Phi"])


def score(args):
    begin,cpu=time.perf_counter(),time.process_time();task=inputs()
    sealed=common.verify(common.resolve(args.seal),[MANIFEST],[CASE],1,STATUS)
    out=common.resolve(args.output)
    if out.exists():raise FileExistsError(out)
    prefix=public_prefix(task)
    # Only this existing disclosed label is opened after complete one-task seal.
    r8=Path(sealed["r8_repo"])
    mapping=common.load(r8/"research/research_max_struct_r6/sealed_truth/confirmation/mapping_PRIVATE.json")
    meta=next(r for r in mapping["records"] if r["case_id"]==CASE)
    path=common.resolve(meta["truth_path"],r8)
    with np.load(path,allow_pickle=False) as z:
        truth=dict(m=z["m"].ravel().copy(),mask=z["mask"].ravel().astype(bool),amplitude=12.,
                   same_grid=meta.get("inverse_crime"),truth_path=str(path))
    truths={CASE:truth}
    rows,stages,budgets,materials,metric_path=common.score_entries(sealed,truths,"POSTHOC_DEV_WEAK_HANDOFF_ONLY")
    row=rows[0]
    row.update(reused_initialization_cost=prefix,
        standalone_CPU_seconds=row["cpu_seconds"]+prefix["CPU_seconds"],
        standalone_LU=row["LU"]+prefix["physical"].get("LU",0),
        standalone_complex_RHS=row["complex_RHS"]+prefix["physical"].get("actual_complex_RHS",0),
        standalone_current_Cholesky=row["current_Cholesky"]+prefix["auxiliary"].get("current_Cholesky",0),
        standalone_current_RHS=row["current_RHS"]+prefix["auxiliary"].get("current_RHS",0),
        cost_scope="cpu_seconds/LU/RHS are actual new continuation only; standalone_* additionally charges reused first-stage cost. Never add reused prefix again to actual research ledger.")
    for b in budgets:
        b.update(standalone_CPU_seconds=b["CPU_seconds"]+prefix["CPU_seconds"] if b["CPU_seconds"] is not None else None,
            standalone_LU=b["LU"]+prefix["physical"].get("LU",0) if b["LU"] is not None else None,
            standalone_complex_RHS=b["complex_RHS"]+prefix["physical"].get("actual_complex_RHS",0) if b["complex_RHS"] is not None else None,
            standalone_prefix_semantics=prefix["CPU_semantics"])
    common.write_bundle(out,rows,stages,budgets,materials,truths,sealed,metric_path,"SCORED_WEAK_HANDOFF_CONTROL",(cpu,begin))
    print(dict(status="SCORED_WEAK_HANDOFF_CONTROL",new_physical_calls=0,old43_and3_unchanged=True))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("mode",choices=("seal","score"))
    ap.add_argument("--r8-repo",default=str(REPO.parent/"r8_worktree"))
    ap.add_argument("--seal",default=str(ROOT/"audit/WEAK_HANDOFF_PUBLIC_SEAL.json"))
    ap.add_argument("--output",default=str(ROOT/"RESULTS/weak_handoff_scoring"))
    args=ap.parse_args();(seal if args.mode=="seal" else score)(args)


if __name__=="__main__":main()
