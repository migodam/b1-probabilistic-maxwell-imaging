"""Sealed-pair scoring at the original independent receiving arc."""
import json
from time import perf_counter, process_time

from continuation import HERE, FROZEN, Model, VerifiedWave, compensation, np, threadpool_limits, write_arrays, write_json
from joint_recheck import physical_receipt


def held(score_shapes):
    """Open the held A-source receive arc only after input shapes are sealed."""
    with np.load(FROZEN/"ports"/"public"/"cases"/"n32_c01.npz",allow_pickle=False) as a:
        y, radius = a["Y_hold"], float(a["radius_hold"])
    with np.load(FROZEN/"ports"/"public"/"operators"/"n32_f400_a08.npz",allow_pickle=False) as a:
        op = dict(a)
    model_op = dict(op,S=op["S_hold"])
    model = Model(model_op,y,radius,"n32_f400_a08 / independent held receive arc",32)
    verify = VerifiedWave(model_op)
    score = {}
    raw = dict(Y_hold=y, radius_hold=radius, port_indices=op["hold_index"])
    for name,x in score_shapes.items():
        score[name],raw[name+"_current"],raw[name+"_prediction"] = physical_receipt(model,verify,x)
    return score,raw


def b_scoring():
    tick, cpu = perf_counter(), process_time()
    source = HERE/"B_independent_trajectory.npz"
    receipt = json.loads(source.with_suffix(".json").read_text())
    if max(receipt["independent_endpoint_audits"][name]["B"]["residual_upper_over_r"]
           for name in ("vacuum","high_patch")) >=1:
        raise ValueError("B endpoints were not sealed as compatible")
    with np.load(source, allow_pickle=False) as a:
        pair = {name:a[name] for name in ("vacuum","high_patch")}
    scores,raw = held(pair)
    result = dict(source=source.name, scoring_only=True, acquisition="original A source columns and 8 independent held receivers",
                  frequency_MHz=400, reads_complex=64, held_radius=float(raw["radius_hold"]), scores=scores,
                  independent_original_pair_distance_over_two_hold_r=float(np.linalg.norm(
                      raw["high_patch_prediction"]-raw["vacuum_prediction"])/(2*raw["radius_hold"])),
                  wall_s=perf_counter()-tick,cpu_s=process_time()-cpu)
    write_json(HERE/"B_independent_held.json",result)
    write_arrays(HERE/"B_independent_held.npz",**raw)


def combine():
    files = {k:HERE/f"joint_original_AB_{k}.json" for k in ("vacuum","high")}
    if not all(f.exists() for f in files.values()):
        return
    reports = {k:json.loads(f.read_text()) for k,f in files.items()}
    decision = {key:reports[name]["selected"][key] for name,key in (("vacuum","vacuum"),("high","high_patch"))}
    pair, raw = {}, {}
    for name, key in (("vacuum","vacuum"),("high","high_patch")):
        with np.load(files[name].with_suffix(".npz"), allow_pickle=False) as a:
            pair[key] = a[key+"_"+decision[key]["seed"]+"_material"]
            raw.update(roi=a["roi"],patch=a["patch"],Y_A=a["Y_A"],Y_B=a["Y_B"])
    cost = compensation((pair["high_patch"]-pair["vacuum"])[~raw["roi"]],(2/32)**2)
    strict = all(d["jointly_verified"] for d in decision.values())
    joint = dict(selected=decision, both_families_strictly_feasible=strict, case="n32_c01",
                 original_A_and_B=True, pulses=16, reads_complex=128,
                 sources={k:v.name for k,v in files.items()},
                 exterior_cost_of_selected_pair=cost, per_run_costs={k:v["cost"] for k,v in reports.items()},
                 heldout=None)
    raw.update(pair)
    # Write the selected material arrays before the sealed held-field scoring.
    write_json(HERE/"joint_original_AB.json",joint)
    write_arrays(HERE/"joint_original_AB.npz",**raw)
    feasible = {k:v for k,v in pair.items() if decision[k]["jointly_verified"]}
    if feasible:
        joint["heldout"],held_raw=held(feasible)
        raw.update({"held_"+k:v for k,v in held_raw.items()})
        write_json(HERE/"joint_original_AB.json",joint)
        write_arrays(HERE/"joint_original_AB.npz",**raw)


if __name__=="__main__":
    with threadpool_limits(1):
        b_scoring()
        combine()
