"""Diagnostic rerun of C with the exact vacuum rows eliminated before screening."""
from collections import Counter
import json
from time import perf_counter

from common import PUBLIC, RESULTS, digest, load, save, write_json
from bounds import Acquisition, sector, tail_proposal, threshold_options, verify_inertia, weight_vector
import numpy as np
from run import ArrayStore, PARAMS
from replay import original, replay_terms


def check(acq, box, frozen, store):
    start = perf_counter()
    ell, upp = box
    keep = np.flatnonzero(upp != 0.)
    reduced = Acquisition(acq["G"][np.ix_(keep,keep)],acq["S"][:,keep],acq["B"][keep],
                          acq["Y_block"],acq["radius"])
    trials = []
    for key in PARAMS:
        w,o = [r[keep] for r in weight_vector(len(ell),*key)]
        nominal = reduced.receiver + sector(reduced.J,reduced.E,ell[keep],upp[keep],w,o)
        candidate = tail_proposal(nominal,len(keep),reduced.radius)
        trials.append((nominal,candidate,w,o,key))
    trials.sort(key=lambda item:item[1]["strength"],reverse=True)
    proposal_end = perf_counter()
    verify, independently_replayed = None, None
    for nominal,candidate,w,o,key in trials[:4]:
        opts = list(threshold_options(candidate,reduced.radius))
        if not opts:
            continue
        exact = original(reduced,ell[keep],upp[keep],w,o)
        for terms in opts:
            verify = verify_inertia(exact,nominal,len(keep),terms,reduced.radius,store.add)
            if verify:
                verify["weights"] = w.tolist()
                verify["optical"] = o.tolist()
                break
        if verify:
            independently_replayed = replay_terms(exact,verify,store.arrays,reduced.radius)
            if not independently_replayed:
                raise ArithmeticError("Exact-known-zero control failed replay")
            break
    if verify is not None and frozen["witness"]["audit"]["proved"]:
        raise ArithmeticError("Exact-known-zero control contradicts verified original witness")
    return dict(box_positions=frozen["box_positions"], C=frozen["status"]["C"],
                E=frozen["status"]["E"],
                E_new=frozen["methods"]["E"]["stage"]=="added_physical_rows",
                C_known_zero=("EXCLUDED" if verify is not None else "WITNESS_FOUND" if
                              frozen["witness"]["audit"]["proved"] else "UNRESOLVED"),
                kept_native_rows=keep.tolist(), verified=verify, independently_replayed=independently_replayed,
                proposal_seconds=proposal_end-start, strict_including_replay_seconds=perf_counter()-proposal_end)


def run_control():
    root = RESULTS/"known_zero_control"
    root.mkdir(exist_ok=True)
    arrays = load(PUBLIC)
    for pair in range(24):
        stem = f"pair_{pair+1:02d}"
        output = root/f"{stem}.json"
        if output.exists():
            continue
        primary = RESULTS/"raw"/f"{stem}.json"
        frozen = json.loads(primary.read_text())
        store = ArrayStore()
        nodes = []
        for row in frozen["rows"]:
            t,bg,k,slot = row["box_positions"]
            input_box = arrays["boxes"][pair,bg,k,slot]
            data = arrays | dict(Y_block=arrays["Y"][pair,t,bg],radius=float(arrays["radii"][pair,t,bg]))
            nodes.append(check(data,input_box,row,store))
        save(root/f"{stem}_frames.npz",**store.arrays)
        write_json(output,dict(pair=pair+1,nodes=nodes,original_digest=digest(primary)))
        print(stem,Counter(row["C_known_zero"] for row in nodes),flush=True)
    summarize(root)


def summarize(root):
    true_box = {}
    with (RESULTS/"all_nodes.csv").open() as handle:
        import csv
        for row in csv.DictReader(handle):
            pair = int(row["pair_id"][-2:])
            bg = {"0.0":0,"0.03125":1,"0.0625":2,"0.125":3}[row["background_cap"]]
            ratio = {"0.5":0,"0.9":1,"1.0":2,"1.5":3,"2.0":4,"3.0":5,"5.0":6}[row["width_ratio"]]
            true_box[pair,int(row["target_idx"]),bg,ratio,int(row["box_idx"])] = (
                row["contains_generating_material"]=="True")
    counts, e_counts, c_counts, timings = Counter(), Counter(), Counter(), Counter()
    passed, missed = [], []
    for p in sorted(root.glob("pair_*.json")):
        data = json.loads(p.read_text())
        pair = data["pair"]
        for row in data["nodes"]:
            status = row["C_known_zero"]
            counts[status] += 1
            e_counts[row["E"],status] += 1
            c_counts[row["C"],status] += 1
            timings["proposal"] += row["proposal_seconds"]
            timings["strict_and_replay"] += row["strict_including_replay_seconds"]
            key = (pair,*row["box_positions"])
            if true_box[key] and status=="EXCLUDED":
                raise ArithmeticError(f"Truth was excluded in exact-known-row control: {key}")
            if row["E_new"]:
                (passed if row["independently_replayed"] else missed).append(list(key))
    output = dict(research_question="Can ordinary full inertia give the same added E strict decisions after globally applying known vacuum rows?",
                  registered_primary_arms_unchanged=True, diagnosis_after_raw_primary_decisions=True,
                  uses_only_public_material_endpoints_and_Y_in_test=True, cohort_boxes=sum(counts.values()),
                  states=dict(counts),
                  C_vs_control={f"{original}/{control}":n for (original,control),n in c_counts.items()},
                  E_vs_control={f"{original}/{control}":n for (original,control),n in e_counts.items()},
                  E_new_reproduced_with_known_zero_elimination=passed,
                  E_new_missed_with_known_zero_elimination=missed,
                  timing_seconds=dict(timings), truth_false_exclusions=0)
    write_json(RESULTS/"known_zero_control.json",output)


if __name__=="__main__":
    run_control()
