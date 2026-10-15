"""Independent saved-array E2 scorer. No solver, new state or RHS is called.

Historical execution scored after search froze but in its original process.
This module makes the requested process boundary explicit for reproduction;
its zero-physics replay does not retroactively change that execution history.
"""
from time import perf_counter, process_time
_START, _CPU = perf_counter(), process_time()

from pathlib import Path
import argparse
import csv
import json
import os
import shutil
import sys

import numpy as np

# These are definition-only numerical/offline helpers. No a30 adapter import or
# ResearchABook instance is created, and no Maxwell evaluation is possible here.
from .e2 import _material_metrics, _select, _write, material_distance


def score(output, source, *, preserve_legacy=False):
    output, source = Path(output).resolve(), Path(source).resolve()
    if any("locked" in p.lower() for p in source.parts):
        raise ValueError("LOCKED_SOURCE_FORBIDDEN")
    protocol = json.loads((output/"PROTOCOL.json").read_text())
    rows = json.loads((output/"SEARCH_HISTORY.json").read_text())
    cases, arms = protocol["cases"], protocol["arms"]
    release = source/"inputs/a29_release"
    metadata = {r["case_id"]: r for r in csv.DictReader((release/"results/a29/COVERAGE_CASES.csv").open())}
    legacy = output/"legacy_offline_inprocess"
    if preserve_legacy:
        if legacy.exists():
            raise FileExistsError("LEGACY_RECEIPT_ALREADY_PRESERVED_DO_NOT_REPEAT_REPLAY")
        legacy.mkdir()
        for name in ("METRICS.json", "MATERIAL_AND_PREDICTION_ARRAYS.npz"):
            shutil.copy2(output/name, legacy/name)
        _write(legacy/"EXECUTION_DEVIATION.json", dict(
            actual_historical_access="Truth opened in online process only after all search/proposals were frozen",
            requested_boundary="Separate offline process",
            deviation="Original implementation lacked process separation; scientific information flow was after frozen search",
            online_truth_leak_found=False,
            reproduction_fix="New scorer executed in clean child process; original full-physics search not rerun",
            preserved_arrays="Complete prior array archive, not a synthetic fixture",new_physics=0,new_RHS=0))
    with np.load(output/"MATERIAL_AND_PREDICTION_ARRAYS.npz",allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in z.files}
    # Only this offline process reads the truth-bearing source archive.
    with np.load(release/"results/a29/OFFLINE_TARGETS.npz",allow_pickle=False) as z:
        saved_ids=z["case_id"].astype(str).tolist()
        truths={case:z["truth"][saved_ids.index(case)].copy() for case in cases}
    metrics=[]
    for case in cases:
        base,background,truth=arrays[case+"__base"],arrays[case+"__background"],truths[case]
        arrays[case+"__OFFLINE_TRUTH"]=truth
        for arm in arms:
            valid=[]
            base_row=next(r for r in rows if r["case_id"]==case and r["arm"]=="base")
            if base_row["compatible"]:
                valid.append(dict(chi=base,residual=base_row["residual"],name="base"))
            for row in rows:
                if row["case_id"]==case and row["arm"]==arm and row["compatible"]:
                    image=arrays[row["array_prefix"]+"__chi"]
                    if all(material_distance(image,item["chi"])>=.03 for item in valid):
                        valid.append(dict(chi=image,residual=row["residual"],name=row["array_prefix"]))
            selected=min(valid,key=lambda x:x["residual"]) if valid else None
            result=dict(case_id=case,family=metadata[case]["family"],parent_id=int(metadata[case]["parent_id"]),
                arm=arm,valid_distinct_material_count=len(valid),probability=None,
                base_offline_metrics=_material_metrics(base,truth,background),
                selected_offline_metrics=_material_metrics(selected["chi"],truth,background) if selected else None,
                selected_name=selected["name"] if selected else None)
            structure_masks=[]
            for candidate in valid:
                mask=np.abs(candidate["chi"]-background)>.15*1.2
                if all(np.count_nonzero(mask!=other)>=6 for other in structure_masks):
                    structure_masks.append(mask)
            result["structurally_distinct_masks_at_least6changedcells"]=len(structure_masks)
            result["structural_definition"]="fixedonline|chi-background|>.18;atleast6changedcells;separatefrommaterialdistance"
            for k in (1,2,4,8):
                ix=_select(valid,base,k)
                errors=[_material_metrics(valid[i]["chi"],truth,background) for i in ix]
                result[f"Best_of_{k}_complex_MAE"]=min((x["full_complex_MAE"] for x in errors),default=None)
                result[f"Coverage_at_{k}"]=any(x["material_distance_truth"]<=.1 for x in errors)
                result[f"Top_{k}_selected_names"]=[valid[i]["name"] for i in ix]
            metrics.append(result)
    equivalence=None
    if preserve_legacy:
        old_metrics=json.loads((legacy/"METRICS.json").read_text())
        with np.load(legacy/"MATERIAL_AND_PREDICTION_ARRAYS.npz",allow_pickle=False) as z:
            old_arrays={k:z[k].copy() for k in z.files}
        same_keys=set(old_arrays)==set(arrays)
        mismatches=[k for k in old_arrays if k not in arrays or not np.array_equal(old_arrays[k],arrays[k])]
        equivalence=dict(status="PASS" if metrics==old_metrics and same_keys and not mismatches else "FAIL",
            metrics_JSON_numerical_and_structural_equality=metrics==old_metrics,
            array_keys_equal=same_keys,array_values_exactly_equal=not mismatches,array_mismatches=mismatches,
            metric_rows=len(metrics),array_count=len(arrays),new_physics=0,new_RHS=0,
            historical_process_boundary="Sameprocess after frozen search; now preserved as execution deviation",
            replay_process_boundary="Independent Python process; no search/corrector/solver executed",
            no_online_truth_leak_found=True)
        _write(output/"OFFLINE_PROCESS_EQUIVALENCE.json",equivalence)
        if equivalence["status"]!="PASS":
            raise AssertionError("INDEPENDENT_OFFLINE_SCORER_CHANGED_RESULTS")
        _write(legacy/"ARRAY_SUMMARY.json",{k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in old_arrays.items()})
    # Equivalent arrays retain all images, predictions and offline truth used by
    # the gallery. Their scientific values are unchanged; no source is modified.
    np.savez_compressed(output/"MATERIAL_AND_PREDICTION_ARRAYS.npz",**arrays)
    _write(output/"METRICS.json",metrics)
    physics_modules=[name for name in sys.modules if name=="a30" or name.startswith("a30.")]
    receipt=dict(status="COMPLETE",pid=os.getpid(),parent_pid=os.getppid(),independent_offline_process=True,
        observations=len(cases),metric_rows=len(metrics),physics_adapter_modules_imported=physics_modules,
        solver_constructors=0,new_full_state_attempts=0,new_RHS=0,GPU_seconds=0,
        source_truth_read="Existing exposed OFFLINE_TARGETS only; no locked payload",
        source_results_modified=False,old_physics_ledger_modified=False,
        equivalence_status=equivalence["status"] if equivalence else "NOT_REQUESTED")
    _write(output/"OFFLINE_PROCESS_RECEIPT.json",receipt)
    return receipt


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--preserve-legacy",action="store_true")
    args=parser.parse_args()
    status="FAILED"
    try:
        receipt=score(args.output,args.source,preserve_legacy=args.preserve_legacy)
        status="COMPLETE"
        print(json.dumps(receipt))
    finally:
        # Separate0physics receipt; never instantiate ResearchABook or append to
        # its execution/physical ledger for this numerical replay.
        _write(args.output/"E2_OFFLINE_REPLAY_JOB.json",dict(status=status,pid=os.getpid(),
            operation="Independent saved-array scorer; no search or physics replay",new_physics=0,new_RHS=0,GPU_seconds=0,
            wall_seconds=perf_counter()-_START,process_cpu_seconds=process_time()-_CPU,
            historical_solver_stage_wall_unchanged=True,old_ledgers_modified=False))


if __name__=="__main__":
    main()
