"""Compile actual R9 receipts only: no arrays, scores, GT, solvers or plotting."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import json
import time

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
MANIFESTS={
    "development":["baseline_development.json","baseline_historical_fixed.json","landscape_development.json",
                   "baseline_cost_controls.json","profile_restart_control.json","curvature_followup.json","weak_handoff_control.json"],
    "confirmation":["confirmation_baselines.json","confirmation_landscape.json","confirmation_single_noise.json","confirmation_weak_handoff.json"]}


def load(p):return json.loads(Path(p).read_text())
def resolve(p):
    p=Path(p);return p if p.is_absolute() else REPO/p
def relative(p):
    p=Path(p)
    return str(p.relative_to(REPO)) if p.is_relative_to(REPO) else str(p)
def number(x):return isinstance(x,(int,float)) and not isinstance(x,bool)
def sum_known(values):return sum(v for v in values if number(v))
def utc(s):return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()


def summarize(records):
    physical,auxiliary,peaks=Counter(),Counter(),{}
    for r in records:
        for k,v in (r.get("physical") or {}).items():
            if number(v):physical[k]+=v
        for k,v in (r.get("auxiliary") or {}).items():
            if not number(v):continue
            if k.startswith("peak_"):
                peaks[k]=max(peaks.get(k,0),v)
            elif k not in ("wall_seconds","cpu_seconds"):
                auxiliary[k]+=v
    return dict(records=len(records),known_CPU_seconds=sum_known(r.get("CPU_seconds") for r in records),
        CPU_unmeasured_records=[r["source"] for r in records if not number(r.get("CPU_seconds"))],
        task_wall_seconds_sum=sum_known(r.get("wall_seconds") for r in records),physical=dict(physical),
        wall_unmeasured_records=[r["source"] for r in records if not number(r.get("wall_seconds"))],
        auxiliary=dict(auxiliary),peak_max_across_receipts=peaks,
        adapter_rank_selection_receiver_SVD=sum_known(r.get("adapter_receiver_SVD") for r in records),
        FFT_batch_calls=sum_known(r.get("FFT_batch_calls") for r in records))


def inverse_record(task,group):
    path=resolve(task["output_dir"])/"RECEIPT.json";d=load(path)
    if d.get("status") not in ("SEALED_PUBLIC_PREDICTION","COMPLETED_PUBLIC_PREDICTION","FAILED"):
        raise RuntimeError("All requested inversions must be terminal before cost compilation: "+str(path))
    landscape="algorithm" in task
    physical=(d.get("cost") or {}).get("physical",{}) if landscape else d.get("physical_model_cost",{})
    execution=d.get("auxiliary_execution_receipt") or {}
    auxiliary=(d.get("cost") or {}).get("auxiliary",{}) if landscape else execution.get("operation_counts",{})
    CPU=d.get("process_CPU_seconds") if landscape else d.get("r9_total_cpu_seconds",d.get("cpu_seconds"))
    wall=d.get("process_wall_seconds") if landscape else d.get("r9_total_wall_seconds",d.get("wall_seconds"))
    decision=d.get("dpcs_rule") or {}
    if not isinstance(decision,dict):decision={}
    timing=dict(inclusive_actual_CPU_seconds=CPU,
        inherited_wrapper_CPU_seconds=d.get("cpu_seconds") if not landscape else None,
        rank_adapter_preparation_CPU_seconds=d.get("adapter_preparation_cpu_seconds"),
        setup_CPU_seconds=d.get("setup_cpu_seconds"),
        auxiliary_observation_inclusive_CPU_seconds=d.get("auxiliary_cpu_seconds"),
        auxiliary_execution_CPU_seconds=execution.get("cpu_seconds"),
        corrector_observation_inclusive_CPU_seconds=(d.get("corrector") or {}).get("cpu_seconds"),
        physical_operation_CPU=d.get("physical_operation_cpu"),
        Aplus_readout_measured_subset_CPU_seconds=auxiliary.get("aplus_readout_cpu_microseconds",0)/1e6 if "aplus_readout_cpu_microseconds" in auxiliary else None,
        landscape_inner_CPU_seconds=(d.get("cost") or {}).get("CPU_seconds") if landscape else None,
        landscape_pre_solver_CPU_seconds=d.get("pre_solver_CPU_seconds"),
        isolated_logging_CPU_seconds=None,
        interpretation="Inclusive baselines retain iteration callbacks, snapshots and serialization; timed compute subsets overlap and do not sum to total. Logging alone is not measured.")
    return dict(task_id=task["task_id"],case_id=task["case_id"],group=group,status=d["status"],
        method=task.get("algorithm",task.get("arm")),dpcs_rule=task.get("dpcs_rule"),source=relative(path),
        CPU_seconds=CPU,wall_seconds=wall,physical=physical,auxiliary=auxiliary,
        FFT_batch_calls=execution.get("FFT"),
        adapter_receiver_SVD=decision.get("new_receiver_SVD",0),
        LU_matrix_dimension=576,current_Cholesky_matrix_dimension=576,
        timing_components=timing,
        cost_scope="Raw executed task only. No reused first-stage, historical R8 endpoint, portfolio sum or already-attributed diagnostic prefix added.")


def recipes(tasks,records):
    """Stand-alone costs are comparisons, not additional actual research work."""
    by_id={r["task_id"]:r for r in records};result=[]
    for task in tasks:
        if task.get("algorithm")!="weak_profile_handoff" and not task["task_id"].endswith("__weak_profile_handoff"):
            continue
        path=resolve(task.get("warmup_receipt",ROOT/"RESULTS/development"/(task["case_id"]+"__profile_homotopy")/"RECEIPT.json"))
        budget=resolve(task.get("warmup_budget",path.parent/"BUDGET.json"))
        receipt=load(path);point=next(p for p in load(budget) if p["label"]=="profile_0_0.01_endpoint")
        prefix=receipt.get("pre_solver_CPU_seconds");semantics="EXACT_RECORDED_PRE_SOLVER_PREFIX"
        if prefix is None:
            prefix=max(0.,receipt["process_CPU_seconds"]-receipt["cost"]["CPU_seconds"])
            semantics="UPPER_SHIFT_including_output_tail"
        cost=point["CPU_seconds"]+prefix;r=by_id[task["task_id"]]
        result.append(dict(task_id=task["task_id"],group=r["group"],
            actual_new_CPU_seconds=r["CPU_seconds"],reused_initialization_CPU_seconds=cost,
            standalone_CPU_seconds=r["CPU_seconds"]+cost,
            standalone_LU=r["physical"].get("LU",0)+point["physical_cost"].get("LU",0),
            standalone_complex_RHS=r["physical"].get("actual_complex_RHS",0)+point["physical_cost"].get("actual_complex_RHS",0),
            standalone_current_Cholesky=r["auxiliary"].get("current_Cholesky",0)+point["auxiliary_cost"].get("current_Cholesky",0),
            standalone_current_RHS=r["auxiliary"].get("current_RHS",0)+point["auxiliary_cost"].get("current_RHS",0),
            initializer_source=relative(budget),prefix_semantics=semantics,
            count_scope="Paid endpoint before unnecessary poststage details factorization. Reused work is already in full profile actual receipt; do not add to grand total."))
    return result


def queue_records():
    result=[]
    for path in sorted(ROOT.glob("delegated/*queue*/RECEIPT.json")):
        d=load(path)
        if d.get("status") not in ("COMPLETED","COMPLETED_WITH_ERRORS"):
            raise RuntimeError("Queue still active: "+str(path))
        start=d.get("started_utc");end=d.get("finished_utc",d.get("ended_utc"))
        result.append(dict(source=relative(path),status=d["status"],jobs=len(d.get("jobs",[])),
            started_utc=start,finished_utc=end,measured_makespan_seconds=d.get("wall_seconds"),
            UTC_interval_seconds=utc(end)-utc(start) if start and end else None,
            workers=d.get("workers",d.get("maximum_numeric_processes")),
            CPU_seconds=d.get("CPU_seconds"),
            note="Controller elapsed interval; child task CPU/wall already counted elsewhere. Do not add as another compute charge."))
    intervals=[(utc(r["started_utc"]),utc(r["finished_utc"]),r["source"]) for r in result if r["started_utc"] and r["finished_utc"]]
    events=sorted({v for a,b,_ in intervals for v in (a,b)});segments=[]
    for a,b in zip(events,events[1:]):
        active=[name for x,y,name in intervals if x<b and y>a]
        if active:segments.append(dict(start_utc=datetime.fromtimestamp(a,timezone.utc).isoformat(),end_utc=datetime.fromtimestamp(b,timezone.utc).isoformat(),seconds=b-a,active_queues=active))
    union=sum(s["seconds"] for s in segments);overlap=sum(s["seconds"] for s in segments if len(s["active_queues"])>1)
    for r in result:
        r["exclusive_nonoverlap_seconds"]=sum(s["seconds"] for s in segments if s["active_queues"]==[r["source"]])
        r["overlap_with_other_queue_seconds"]=sum(s["seconds"] for s in segments if len(s["active_queues"])>1 and r["source"] in s["active_queues"])
    return dict(records=result,interval_segments=segments,queue_interval_union_seconds=union,
        time_with_multiple_queues_active_seconds=overlap,
        exactly_one_queue_active_seconds=union-overlap,
        first_to_last_queue_span_seconds=max(b for _,b,_ in intervals)-min(a for a,_,_ in intervals),
        sum_queue_makespans_seconds=sum_known(r["measured_makespan_seconds"] for r in result),
        scope="UTC interval union distinguishes overlaps/gaps. It is queue-active elapsed time, not total research elapsed, CPU utilization or summed task wall.")


def timed_analysis():
    paths=[ROOT/"delegated/baseline_audit/ANALYSIS_RECEIPT.json"]
    paths+=sorted(ROOT.glob("RESULTS/*scoring/RECEIPT.json"))
    paths+=sorted(ROOT.glob("audit/*plot_receipt.json"))
    result=[]
    for path in paths:
        d=load(path);attempts=d if isinstance(d,list) else [d]
        # Development plot array retains both the failed layout and repaired run.
        for i,r in enumerate(attempts):
            CPU=r.get("CPU_seconds",r.get("scoring_CPU_seconds"))
            wall=r.get("wall_seconds",r.get("scoring_wall_seconds"))
            result.append(dict(source=relative(path),attempt=i,status=r.get("status"),CPU_seconds=CPU,wall_seconds=wall,
                physical={},auxiliary={"public_operator_SVD":r.get("operator_SVDs",0)},
                timing_scope=r.get("scope","Recorded analysis/scoring function timing; unmeasured imports/export tails remain excluded")))
    return result


def main():
    begin,cpu=time.perf_counter(),time.process_time();records=[];all_tasks=[];sources=[]
    for group,names in MANIFESTS.items():
        tasks=[]
        for name in names:
            path=ROOT/"configs"/name;sources.append(relative(path));tasks+=load(path)["tasks"]
        expected=47 if group=="development" else 66
        if len(tasks)!=expected or len({t["task_id"] for t in tasks})!=expected:
            raise ValueError("Expected DEV47 and CONF66 actual inverse tasks")
        all_tasks+=tasks;records += [inverse_record(t,group) for t in tasks]
    if len({r["source"] for r in records})!=113:raise ValueError("Duplicate executed receipt")
    diagnostics=[]
    for path in sorted(ROOT.glob("RESULTS/curvature/*/RECEIPT.json")):
        d=load(path);cost=d["cost"]
        diagnostics.append(dict(source=relative(path),status=d["status"],CPU_seconds=cost.get("CPU_seconds"),wall_seconds=cost.get("wall_seconds"),
            physical=cost["physical"],auxiliary={"dense_eigendecompositions":cost.get("dense_eigendecompositions",0)},
            note="Count once here; curvature proposal inversion raw receipt excludes this diagnostic prefix"))
    if len(diagnostics)!=2:raise ValueError("Expected two curvature diagnostics")
    path=ROOT/"audit/gradient_wiring/RECEIPT.json";d=load(path);cost=d["cost"]
    gradient=[dict(source=relative(path),CPU_seconds=cost.get("CPU_seconds"),wall_seconds=cost.get("wall_seconds"),
        physical=cost["physical"],auxiliary=cost["auxiliary"],timing_scope="Inner wiring experiment CPU; process imports/wall unmeasured")]
    path=ROOT/"audit/generation_confirmation.json";d=load(path)
    generation=[dict(source=relative(path),CPU_seconds=d.get("CPU_seconds"),wall_seconds=d.get("wall_seconds"),objects=d["objects_completed"],
        physical=dict(Counter({k:sum(c["physical_cost"].get(k,0) for c in d["costs"]) for k in {k for c in d["costs"] for k in c["physical_cost"]}})),
        auxiliary={"operator_builds":d.get("operator_builds",0)},
        timing_scope="Whole generation receipt, including nested setup and per-object solves; do not add those subtimers again")]
    archive=load(ROOT/"delegated/development_queue/ARCHIVED_FAILED_OUTPUTS.json");failed=[]
    for item in archive["maps"]:
        path=resolve(item["preserved_path"])/"RECEIPT.json"
        if path.exists():
            d=load(path);failed.append(dict(source=relative(path),CPU_seconds=d.get("process_CPU_seconds"),wall_seconds=d.get("process_wall_seconds"),
                physical={},auxiliary={},status=d.get("status"),physical_calls=0))
    analysis=timed_analysis();queues=queue_records()
    groups={"development_inversions":summarize([r for r in records if r["group"]=="development"]),
        "confirmation_inversions":summarize([r for r in records if r["group"]=="confirmation"]),
        "curvature_diagnostics":summarize(diagnostics),"gradient_wiring":summarize(gradient),
        "confirmation_generation":summarize(generation),"analysis_scoring_plots":summarize(analysis),
        "measured_archived_import_failure_subset":summarize(failed)}
    measured=records+diagnostics+gradient+generation+analysis+failed
    actual=summarize(measured)
    standalone=recipes(all_tasks,records)
    fail_schema=ROOT/"delegated/scoring/FAILED_SCORING_ATTEMPT_01.json"
    unmeasured=[dict(kind="import_failures_unmeasured_baseline_children",attempts=16,CPU_seconds=None,wall_seconds=None,physical_calls=0),
        dict(kind="queue_launch_path_typo",attempts=1,CPU_seconds=None,wall_seconds=None,physical_calls=0,
             source=relative(ROOT/"delegated/development_queue_launch_path_typo_FAILED.log")),
        dict(kind="first_scoring_schema_error",attempts=1,CPU_seconds=None,wall_seconds=None,physical_calls=0,
             source=relative(fail_schema)),
        dict(kind="unmetered_orchestration_literature_seal_metadata_exports_plot_imports",CPU_seconds=None,wall_seconds=None,
             note="No complete CPU/LLM/token/network/power/money metering; never infer zero")]
    invocation=dict(created_utc=datetime.now(timezone.utc).isoformat(),CPU_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-begin,source=relative(Path(__file__)),scope="Cost aggregation through pre-export point; exports and final receipt write not metered")
    log=ROOT/"audit/cost_compiler_invocations.jsonl";log.parent.mkdir(parents=True,exist_ok=True)
    with log.open("a") as f:f.write(json.dumps(invocation)+"\n")
    invocations=[json.loads(s) for s in log.read_text().splitlines() if s.strip()]
    actual["known_CPU_seconds"]+=sum(r["CPU_seconds"] for r in invocations)
    result=dict(status="COMPLETE_RECEIPT_AGGREGATION_KNOWN_COST_LOWER_BOUND",created_utc=invocation["created_utc"],
        task_counts=dict(development=47,confirmation=66,total_inversions=113,curvature_diagnostics=2,gradient_wiring=1,generation_objects=6,
                         import_failure_attempts=28,launch_path_failure_attempts=1,scoring_schema_failure_attempts=1),
        actual_measured_lower_bound=actual,groups=groups,inverse_receipts=records,curvature_receipts=diagnostics,
        gradient_receipts=gradient,generation_receipts=generation,analysis_scoring_plot_receipts=analysis,
        archived_import_failure_receipts=failed,failed_import_archive_index=relative(ROOT/"delegated/development_queue/ARCHIVED_FAILED_OUTPUTS.json"),
        unmeasured_costs=unmeasured,queue_makespans_and_nonoverlap=queues,
        standalone_weak_recipes_NOT_ADDITIONAL_ACTUAL_COST=standalone,
        portfolio_scope="Two frozen portfolios pay both branches as method recipes. Their executed branches are already in113 inversions; do not add portfolio totals again.",
        historical_R8_scope="Reused historical endpoints and R8 training excluded from actual new R9 ledger. No R9 network training occurred.",
        counter_semantics="All physical LU here are576-cell systems; full-current Cholesky576 is separate. Batch FFT calls and per-source2D transforms are separate. Peak storage max is not simultaneous task RSS.",
        timing_semantics="Baseline source wrappers include logging/serialization; available A+ readout and physical-operation timings are overlapping subsets. No isolated logging time or pure-compute speedup invented.",
        monetary_cost=None,energy_cost=None,compiler_invocations=invocations,source_manifests=sources,
        source_helper=relative(Path(__file__)),truth_array_access="NONE",new_scientific_solves=0,hash_checks="NONE")
    target=ROOT/"audit/ALL_COSTS.json";target.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    write_markdown(result)
    print(json.dumps(dict(status=result["status"],inversions=113,known_CPU_lower_bound=actual["known_CPU_seconds"],
        LU=actual["physical"].get("LU"),complex_RHS=actual["physical"].get("actual_complex_RHS"),output=relative(target))))


def write_markdown(result):
    total=result["actual_measured_lower_bound"];q=result["queue_makespans_and_nonoverlap"]
    lines=["# R9 actual cost ledger", "", "All113 inverse receipts (DEV47 +CONF66) are terminal. This ledger reads public cost receipts only: no GT arrays, scores, new solves or plot reruns. It reports a known CPU lower bound, not exact total research cost.", "",
        f"Known measured CPU including retained failed subsets and cost-compiler invocations: **{total['known_CPU_seconds']:.6f}s**. Physical LU **{total['physical'].get('LU',0)}**, complex RHS **{total['physical'].get('actual_complex_RHS',0)}**. Currency and energy remain null.", "",
        "| Category | Receipt count | Known CPU(s) | Known task-wall sum(s) | Physical LU | Complex RHS | Current Cholesky | FFT batch calls |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name,g in result["groups"].items():
        lines.append(f"| {name} | {g['records']} | {g['known_CPU_seconds']:.6f} | {g['task_wall_seconds_sum']:.6f} | {g['physical'].get('LU',0)} | {g['physical'].get('actual_complex_RHS',0)} | {g['auxiliary'].get('current_Cholesky',0)} | {g['FFT_batch_calls']} |")
    compilers=result["compiler_invocations"]
    lines.append(f"| cost_compiler_measured_subsets | {len(compilers)} | {sum(r['CPU_seconds'] for r in compilers):.6f} | {sum(r['wall_seconds'] for r in compilers):.6f} | 0 | 0 | 0 | 0 |")
    lines += ["", "Nested setup, auxiliary/corrector timers, per-object generation timers and scoring cost ledgers are not added again. Inverse totals come directly from raw execution receipts, so curvature diagnostic cost is counted once in its own category rather than twice via attributed supplement summaries. Profile first-stage cost is already paid within full profile runs; weak standalone reuse and portfolio recipe sums are excluded from actual totals. Historical R8 endpoints are not new R9 jobs.", "",
        "## Timing and logging boundaries", "", "Inherited baselines use inclusive r9_total_cpu_seconds (imports, preparation, solver callbacks, paid snapshots and serialization); landscape uses process_CPU_seconds. A+ readout microsecond counters, physical-operation timers, inner auxiliary/corrector CPU and landscape startup are retained per receipt in ALL_COSTS.json. They are overlapping subsets. Isolated logging CPU is unknown; the inclusive baseline is not mislabeled pure arithmetic. Peak counters are maxima across receipts, not summed memory occupancy.", "",
        "## Failures and missing metering", "", "The28 initial wrong-interpreter imports used zero physics. Only12 landscape children recorded CPU; their measured subtotal is included, while16 baseline failure CPU entries remain null. The launcher-path typo and first scorer-schema error each retain a source artifact and unknown CPU/wall. The layout-repair plot attempt is included when present in the plot-attempt array. Literature/orchestration, seal checks, warm-initializer exports, plotting imports explicitly outside timers, and other unmetered tooling are not free: total cost remains a lower bound. New training count is0; monetary/energy cost is null, not0.", "",
        "## Queue elapsed time and overlaps", "",
        f"Queue interval union={q['queue_interval_union_seconds']:.6f}s; exactly one queue active={q['exactly_one_queue_active_seconds']:.6f}s; multiple queues active={q['time_with_multiple_queues_active_seconds']:.6f}s. First-to-last queue span={q['first_to_last_queue_span_seconds']:.6f}s, including orchestration gaps. Sum of measured queue makespans={q['sum_queue_makespans_seconds']:.6f}s is separate from UTC intervals; the two clock boundaries are not asserted identical. The recorded queues are nonoverlapping in this run. Task-wall sums above are concurrent workload, not makespan or total research elapsed.", "",
        "| Queue | Jobs | Measured makespan(s) | Exclusive nonoverlap(s) | Overlap with another queue(s) |",
        "|---|---:|---:|---:|---:|"]
    for r in q["records"]:
        lines.append(f"| {Path(r['source']).parent.name} | {r['jobs']} | {r['measured_makespan_seconds']:.6f} | {r['exclusive_nonoverlap_seconds']:.6f} | {r['overlap_with_other_queue_seconds']:.6f} |")
    lines += ["", "UTC interval segments and actual timestamps are preserved in ALL_COSTS.json. Queue-active overlap is not a claim about SM/core utilization. Queue-controller wall or child CPU is not charged twice.", "",
        "## Weak handoff standalone recipes", "", "These are method costs only; actual new continuation cost is already in the inverse table. Initialization includes exact paid first-stage prefix plus measured startup, before unneeded poststage details factorization. Necessary current Cholesky/RHS are retained.", "",
        "| Task | Actual new CPU(s) | Reused init CPU(s) | Standalone CPU(s) | Standalone LU | Standalone complex RHS | Standalone current Cholesky |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for r in result["standalone_weak_recipes_NOT_ADDITIONAL_ACTUAL_COST"]:
        lines.append(f"| {r['task_id']} | {r['actual_new_CPU_seconds']:.6f} | {r['reused_initialization_CPU_seconds']:.6f} | {r['standalone_CPU_seconds']:.6f} | {r['standalone_LU']} | {r['standalone_complex_RHS']} | {r['standalone_current_Cholesky']} |")
    lines += ["", "Full additive physical/auxiliary counters, FFT transforms, extra rank-selection SVD counts, timing subsets, unknown costs, source receipt paths and source helper are retained in audit/ALL_COSTS.json. Rebuild with experiments/compile_cost_ledger.py; rebuilding only re-reads metadata, appends its measured compiler invocation, and refreshes these derived reports. No raw receipt is changed.", ""]
    (ROOT/"COST_LEDGER.md").write_text("\n".join(lines))


if __name__=="__main__":main()
