#!/usr/bin/env python3
"""Compare candidate law corrections only on the exposed synthetic diagnoses."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from threadpoolctl import threadpool_limits

from cross_engine import iterate, mature

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/"diagnosis"))
from diagnose import CASES, LA, Twin, get_case, read_npz, write
from score import error_values

ARMS = {
    "dual": ("dual",{}),
    "diag16": ("diag",{"strength":16.}),
    "full_ot16": ("full_ot",{"strength":16.}),
    "full16": ("full",{"strength":16.}),
    "regional4_4": ("regional4",{"strength":4.}),
    "regional16_1": ("regional16",{"strength":1.}),
    "pixel025": ("pixel",{"strength":.25}),
    "pixel1": ("pixel",{"strength":1.}),
    "reset": ("reset",{}),
    "diag16_guard": ("diag",{"strength":16.,"stationary_stages":True,"unconditional_pr_reset":True}),
    "full16_guard": ("full",{"strength":16.,"stationary_stages":True,"unconditional_pr_reset":True}),
    "regional16_1_guard": ("regional16",{"strength":1.,"stationary_stages":True,"unconditional_pr_reset":True}),
    "pixel025_guard": ("pixel",{"strength":.25,"stationary_stages":True,"unconditional_pr_reset":True}),
    "pixel025_accept_only": ("pixel",{"strength":.25,"stationary_stages":False}),
    "dual_start065": ("dual",{"initial_level":.65}),
}


def run(proposal):
    name, arm, limit = proposal
    method, settings = ARMS[arm]
    options = dict(settings)
    initial = options.pop("initial_level", 0.)
    if name == Twin:
        options = {**options, "passive": True}
    case,pub = get_case(name)
    folder = HERE/"development"/name
    folder.mkdir(parents=True,exist_ok=True)
    clock = mature.legacy.Clock(limit)
    model = mature.legacy.PaperModel(pub,case,clock)
    Y = pub["Y"]
    m0 = np.full(model.n,initial)
    states = []
    samples = []
    def audit(m, step, stage):
        phi,error,_ = model.objective(m,Y)
        t = clock.wall()
        if t <= limit:
            row=dict(available_wall=t, cpu_seconds=clock.cpu(), step=step, stage=stage,
                     merit=phi, physical_train_relative=error, peak=float(m.max()))
            samples.append(row)
            states.append(m.copy())
        return t+.6

    with threadpool_limits(1):
        status = "SCHEDULE_COMPLETE"
        end_step = 0
        previous = -1
        try:
            due = audit(m0,0,-1)
            for m,step,stage,event in iterate(model,Y,case,method,m0,20.,**options):
                end_step=step
                if clock.wall() >= due or stage != previous:
                    due = audit(m,step,stage)
                previous=stage
            if end_step and samples[-1]["step"]!=end_step:
                audit(m,end_step,stage)
        except mature.legacy.TimeLimit:
            status = "WALL_LIMIT"
        if not samples:
            raise ValueError("No timed physical state was scored.")
        best=min(range(len(samples)),key=lambda i:samples[i]["merit"])
        desc = dict(name=name,arm=arm,method=method,options=options,initial_level=initial,status=status,
                    wall_limit=limit, reached_step=end_step, wall_seconds=clock.wall(),
                    cpu_seconds=clock.cpu(), operations=dict(clock.counts),
                    material_audits=samples, best_index=best,
                    solver_stages=model.trace_route, corrections=model.trace_cross_ot)
        np.savez_compressed(folder/f"{arm}.npz",audited_material=np.asarray(states))
        write(folder/f"{arm}.json",desc)
        print(name,arm,status,end_step,"res",samples[best]["physical_train_relative"],
              "updates",clock.counts["cross_corrections_accepted"],"wall",desc["wall_seconds"],flush=True)


def score_finished(names, arms):
    rows=[]
    for name in names:
        root=HERE/"development"/name
        references=HERE.parent/"diagnosis"/"trajectories"/name/"offline_reference.npz"
        truth=read_npz(references)["m"] if name != Twin else None
        if name == Twin:
            # The held illuminations enter only after all solver arms have returned.
            case, public = get_case(name)
            held = read_npz(LA/"heldout"/f"{Twin}.npz")
            held_y = held["Y"]
            measured = np.isfinite(held_y)
            held_model = mature.legacy.PaperModel(
                dict(E=held["E"], S=held["S"], centers=public["centers"],
                     radius=public["radius"], k=public["k"]),
                case, mature.legacy.Clock(float("inf")))
        for arm in arms:
            file=root/f"{arm}.json"
            log=json.loads(file.read_text())
            materials=read_npz(file.with_suffix(".npz"))["audited_material"]
            i=log["best_index"]
            control = {}
            if truth is not None:
                metrics=error_values(materials,truth)
                score={k:float(v[i]) for k,v in metrics.items()}
                endpoint={k:float(v[-1]) for k,v in metrics.items()}
            else:
                started = perf_counter(),process_time()
                before = held_model.clock.counts.copy()
                scores = {}
                for index in sorted({i,len(materials)-1}):
                    current = held_model.state(materials[index])["j"]
                    error = held_model.ports.receive(current)-held_y
                    scores[index] = float(np.linalg.norm(error[measured])/np.linalg.norm(held_y[measured]))
                score = {"held_source_relative":scores[i]}
                endpoint = {"held_source_relative":scores[len(materials)-1]}
                control = dict(held_source_count=held_y.shape[1],held_complex_ports=int(measured.sum()),
                               held_wall=perf_counter()-started[0],held_cpu=process_time()-started[1],
                               held_operations=dict(held_model.clock.counts-before))
            rows.append(dict(name=name,arm=arm,fit=log["material_audits"][i],shape_scores=score,
                             final=log["material_audits"][-1],final_shape=endpoint,
                             wall=log["wall_seconds"],cpu=log["cpu_seconds"],
                             steps=log["reached_step"],
                             attempted_repairs=len(log["corrections"]),
                             accepted_repairs=sum(r["accepted"] for r in log["corrections"]),**control))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wall",type=float,default=15.)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--case",action="append",choices=CASES)
    parser.add_argument("--arm",action="append",choices=ARMS)
    parser.add_argument("--summary",type=Path,default=HERE/"development_scores.json")
    opts=parser.parse_args()
    names=opts.case or CASES[:3]
    arms=opts.arm or list(ARMS)
    work=[(name,arm,opts.wall) for name in names for arm in arms]
    with ProcessPoolExecutor(max_workers=opts.workers) as ex:
        list(ex.map(run,work))
    with threadpool_limits(1):
        write(opts.summary,dict(wall_limit=opts.wall, rows=score_finished(names,arms)))


if __name__=="__main__":
    main()
