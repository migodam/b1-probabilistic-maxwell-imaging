"""Score all sealed three-arm outputs and compute fixed paired contrasts."""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
if not __package__:
    sys.path.insert(0, str(HERE.parents[2]))
from research.z_alg_xudong_cie.mechanism.design import DATA, FREEZE, configure_registry
from research.z_alg_xudong_cie.holdout import data as public_registry
from research.z_alg_xudong_cie.holdout.data import elapsed, timers, write_json
from research.z_alg_xudong_cie.holdout.scoring import score
from research.z_alg_xudong_cie.holdout.verify_truth import audit

SCORED = HERE/"scores.json"
MEASURES = ("inside_L2", "background_L2", "half_IoU", "rel_L2", "held", "native_training", "updates", "wall_seconds")
COMPARISONS = (("profile","pixel"),("pixel","reset"),("profile","reset"))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sealed_inputs(root, freeze):
    pairs = [(row,method,root/f"{row['case']}__{method}") for row in public_registry.planned_cases()
             for method in freeze["methods"]]
    missing = [str(stamp) for _,_,stamp in pairs
               if not stamp.with_suffix(".json").is_file() or not stamp.with_suffix(".npz").is_file()]
    if missing:
        raise FileNotFoundError(f"All predictions must seal before scoring starts; missing {len(missing)} files.")
    result = []
    for row,method,stamp in pairs:
        receipt = read_json(stamp.with_suffix(".json"))
        if (receipt["observation"] != row["case"] or receipt["method"] != method or
                not receipt.get("public_only") or not receipt.get("sealed_at") or
                receipt["sealed_at"] <= freeze["frozen_utc"]):
            raise ValueError(f"The training data seal failed for {stamp}.")
        result.append((row,method,stamp,receipt))
    return result


def online_cost(receipt):
    return {key: receipt[key] for key in
            ("wall_seconds","cpu_seconds","budget_seconds","coarse_counts","native_counts",
             "search_runs","native_rank_cost","selected","sealed_at","radius_status",
             "relative_native_training_residual")}


def flat_metrics(row):
    return {**{k: row["material96"][k] for k in ("inside_L2","background_L2","half_IoU","rel_L2")},
            "held": row["held_transmit_and_receive"]["relative_to_measured"],
            "native_training": row["native_fitting"][row["case"]]["relative_to_measured"],
            "updates": sum(x.get("steps",0) for x in row["online"]["search_runs"]),
            "wall_seconds": row["online"]["wall_seconds"]}


def paired(rows):
    scores = {(row["case"],row["method"]): row for row in rows
              if row.get("native_forward",{}).get("status")=="scored"}
    output, summaries = [], []
    for left,right in COMPARISONS:
        for case in public_registry.planned_cases():
            name = case["case"]
            if (name,left) not in scores or (name,right) not in scores:
                continue
            a,b = (flat_metrics(scores[name,method]) for method in (left,right))
            output.append(dict(case=name,object_id=case["object_id"], array=case["array"],
                               left=left,right=right,
                               difference={key:a[key]-b[key] for key in MEASURES},
                               inside_background_win=a["inside_L2"]<b["inside_L2"] and
                               a["background_L2"]<b["background_L2"]))
        for aperture in public_registry.APERTURES:
            subset=[r for r in output if r["left"]==left and r["right"]==right and r["array"]==aperture]
            if not subset:
                continue
            physical_metrics=tuple(k for k in MEASURES if k not in ("updates","wall_seconds"))
            wins={key:sum((r["difference"][key] > 0 if key=="half_IoU" else r["difference"][key] < 0)
                          for r in subset) for key in physical_metrics}
            summaries.append(dict(left=left,right=right,array=aperture,n_independent_materials=len(subset),
                                  median={k:float(np.median([r["difference"][k] for r in subset])) for k in MEASURES},
                                  wins=wins, inside_background_joint_wins=sum(r["inside_background_win"] for r in subset)))
    return output,summaries


def run(root, verify):
    freeze=configure_registry()
    inputs=sealed_inputs(root,freeze)
    if verify:
        truth=audit(DATA, FREEZE, HERE/"truth_feasibility.json", freeze["objects"])
        expected=len(list(public_registry.planned_cases()))
        truth["all_requested_fitting_balls_verified"]=(len(truth["observations"])==expected and
             all(row["true_material_in_padded_ball"] for row in truth["observations"]))
        truth["registry_expected_observations"]=expected
        write_json(HERE/"truth_feasibility.json",truth)
    began=timers()
    result=[]
    for row,method,stamp,receipt in inputs:
        with np.load(stamp.with_suffix(".npz"),allow_pickle=False) as src:
            m=src["m"]
        start=timers()
        try:
            out=score(row["case"],m,freeze_record=FREEZE,root=DATA)
            result.append(dict(case=row["case"],method=method,array=row["array"],
                               material_file=str(stamp.with_suffix(".npz").relative_to(HERE)),
                               online=online_cost(receipt),**out))
        except (RuntimeError,ValueError,ArithmeticError) as err:
            result.append(dict(case=row["case"],method=method,array=row["array"],
                               online=online_cost(receipt),failure=str(err),**elapsed(start)))
        print(row["case"],method,result[-1].get("native_forward",{}).get("status","failed"),flush=True)
    dp, summaries=paired(result)
    output=dict(scope="independent after-seal scoring on the two aperture views per independent material",
                rows=result,paired=dp,summaries=summaries,**elapsed(began))
    write_json(SCORED,output)
    with (HERE/"metrics.csv").open("w",newline="",encoding="utf-8") as stream:
        csv_writer=csv.DictWriter(stream,fieldnames=["case","method","native_score_status",*MEASURES])
        csv_writer.writeheader()
        for row in result:
            state=row.get("native_forward",{}).get("status","failed")
            csv_writer.writerow(dict(case=row["case"],method=row["method"],native_score_status=state,
                                     **(flat_metrics(row) if state=="scored" else {})))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--predictions",type=Path,default=HERE/"predictions")
    p.add_argument("--verify-truth",action="store_true")
    args=p.parse_args()
    with threadpool_limits(1):
        run(args.predictions.resolve(),args.verify_truth)
