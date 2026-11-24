"""Select and score the exposed two-branch portfolios after forecasts are saved."""

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from threadpoolctl import threadpool_limits

from cross_engine import mature

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"diagnosis"))
from diagnose import CASES, LA, Twin, get_case, read_npz, write
from score import error_values

PORTFOLIOS = {"with_pixel":"pixel025_accept_only", "with_constant_065":"dual_start065"}


def read(path):
    return json.loads(path.read_text())


def best_fit(name, arm, budget):
    arm_path = HERE/"development"/name / f"{arm}.json"
    log = read(arm_path)
    allowed = [(i,r) for i,r in enumerate(log["material_audits"]) if r["available_wall"]<=budget]
    index,selected = min(allowed,key=lambda pair:pair[1]["merit"])
    return dict(**selected,arm=arm,record_index=index,planned_branch_wall=budget,
                charged_wall=log["wall_seconds"],
                wall_overrun=max(0.,log["wall_seconds"]-budget),
                measured_cpu=log["cpu_seconds"],received_step=log["reached_step"],
                accepted_corrections=sum(q["accepted"] for q in log["corrections"]))


def evaluate_holdout(model, y, valid, material):
    t0,c0 = perf_counter(),process_time()
    before = model.clock.counts.copy()
    value = model.ports.receive(model.state(material)["j"])
    return dict(held_source_relative=float(np.linalg.norm((value-y)[valid])/np.linalg.norm(y[valid])),
                held_source_count=y.shape[1],held_ports=int(valid.sum()),
                paid_held_wall=perf_counter()-t0,paid_held_cpu=process_time()-c0,
                held_operations=dict(model.clock.counts-before))


def main():
    predictions={}
    for case in CASES:
        budget=24. if case==Twin else 15.
        air=best_fit(case,"dual",budget)
        predictions[case]={}
        for name,arm in PORTFOLIOS.items():
            second=best_fit(case,arm,budget)
            winner=min((air,second),key=lambda row:row["merit"])
            predictions[case][name]=dict(winner=winner,air=air,secondary=second,
                                         charged_wall=air["charged_wall"]+second["charged_wall"],
                                         total_budget=2*budget)
    seal=HERE/"portfolio_selections.json"
    write(seal,dict(selection="minimum original physical training objective",
                    choices=predictions,utc=datetime.now(timezone.utc).isoformat()))
    summary={}
    for case in CASES:
        summary[case]={}
        if case==Twin:
            options,pub=get_case(case)
            raw=read_npz(LA/"heldout"/f"{Twin}.npz")
            y=raw["Y"];valid=np.isfinite(y)
            model=mature.legacy.PaperModel(
                dict(E=raw["E"],S=raw["S"],centers=pub["centers"],radius=pub["radius"],k=pub["k"]),
                options,mature.legacy.Clock(float("inf")))
            prior={}
        else:
            truth=read_npz(HERE.parent/"diagnosis"/"trajectories"/case/"offline_reference.npz")["m"]
        for name,row in predictions[case].items():
            selected=row["winner"]
            arm,idx=selected["arm"],selected["record_index"]
            material=read_npz(HERE/"development"/case/f"{arm}.npz")["audited_material"][idx]
            if case==Twin:
                key=arm,idx
                if key not in prior:
                    prior[key]=evaluate_holdout(model,y,valid,material)
                quality=prior[key]
            else:
                quality={key:float(v[0]) for key,v in error_values(material[None,:],truth).items()}
            summary[case][name]=dict(**row,score=quality)
            print(case,name,arm,selected["step"],selected["physical_train_relative"],
                  quality.get("held_source_relative",quality.get("relative_material_L2")),
                  row["charged_wall"],flush=True)
    write(HERE/"portfolio_scores.json",dict(scope="Exposed development cases only",
                                           cases=summary))


if __name__=="__main__":
    with threadpool_limits(1):
        main()
