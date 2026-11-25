"""Matched first-feasible stop with and without prior whole-box pruning."""

import argparse
import json
from statistics import median

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from confirm import OUT, WEIGHTS, geometry
from run import HERE, Border, clock, elapsed
from pool_common import FullWave, write_json
from box_certificates import Bounds
from structured_pool import MAX_CALLS, N_STARTS, RX, TX, proposals

ARMS = ("same_order_unfiltered", "same_order_shifted")
SEED = 453917


def add_times(items):
    return {key: sum(t[key] for t in items) for key in ("wall_seconds", "cpu_seconds")}


def arm(g, s, b, y, radius, seed, prune):
    beginning = clock()
    model = FullWave(g,s,b,y)
    samples, lo, hi, generation = proposals(g,s,b,y,model,seed)

    beginning_screen = clock()
    decisions = []
    if prune:
        certificate = Border(g,s,b)
        certificate.measurement(y)
        decisions = [certificate.check(l,h,y,radius,weights=WEIGHTS,jvalues=(1,))
                     for l,h in zip(lo,hi,strict=True)]
    excluded = {i for i,row in enumerate(decisions) if row["rejected"]}
    screening = elapsed(beginning_screen)

    beginning_rank = clock()
    readings = np.stack([model.evaluate(x)[0] for x in samples])
    residuals = np.linalg.norm((readings-y).reshape(len(samples),-1),axis=1).reshape(-1,N_STARTS)
    best_local = np.argmin(residuals,axis=1)
    per_box = residuals[np.arange(len(lo)),best_local]
    order = np.argsort(per_box,kind="stable").tolist()
    sample_indices = (N_STARTS*np.arange(len(lo))+best_local).tolist()
    rank_time = elapsed(beginning_rank)

    audit, visited, verify_costs = Bounds(g,s,b,y), [], []
    best_material, best_prediction, best_branch, best_residual, best_upper = None,None,None,np.inf,None
    for i in order:
        if i in excluded:
            continue
        branch_start = sample_indices[i]
        material, predicted, result = model.refine(samples[branch_start],lo[i],hi[i],MAX_CALLS)
        residual = result["residual"]
        strict = None
        if residual<=radius:
            start_verify = clock()
            strict = audit.point_upper(material)
            verify_costs.append(elapsed(start_verify))
        visited.append(dict(branch=i,sample_index=branch_start,strict_upper=strict,**result))
        if residual<best_residual:
            best_material,best_prediction,best_branch,best_residual,best_upper = material,predicted,i,residual,strict
        if strict is not None and strict<=radius:
            best_material,best_prediction,best_branch,best_residual,best_upper = material,predicted,i,residual,strict
            break
    if best_material is not None and best_upper is None:
        start_verify = clock()
        best_upper = audit.point_upper(best_material)
        verify_costs.append(elapsed(start_verify))
    return dict(first_feasible=(best_upper is not None and best_upper<=radius),
                reported_branch=best_branch, original_residual=best_residual if best_branch is not None else None,
                verified_residual_upper=best_upper,
                ranking=order, ranked_sample_per_box=sample_indices, excluded=sorted(excluded),
                visited=visited, generation=generation, online_screening=screening, ranking_time=rank_time,
                endpoint_interval_checks=len(verify_costs), endpoint_interval_cost=add_times(verify_costs),
                ordinary_full_wave_counts=dict(model.counts), **elapsed(beginning)), (
                samples,lo,hi,best_material,best_prediction),decisions


def same_optional(a,b):
    return (a is None and b is None) or (a is not None and b is not None and np.array_equal(a,b))


def score(est, truth, clean, g,s,b):
    if est is None:
        return None
    model=FullWave(g,s,b,clean)
    port=model.evaluate(est)[0]
    hold_rx=np.ix_(np.setdiff1d(np.arange(s.shape[0]),RX),TX)
    hold_tx=np.ix_(RX,np.setdiff1d(np.arange(b.shape[1]),TX))
    rel=float(np.linalg.norm(est-truth)/np.linalg.norm(truth))
    return dict(full_relative_l2=rel,active_region_relative_l2=rel,background_rms=0.,
                heldout_receiver=float(np.linalg.norm((port-clean)[hold_rx])/np.linalg.norm(clean[hold_rx])),
                heldout_source=float(np.linalg.norm((port-clean)[hold_tx])/np.linalg.norm(clean[hold_tx])))


def measured_medians(pair):
    values={}
    for method in ARMS:
        a=[v["arms"][method] for v in pair]
        values[method]=dict(wall_seconds=median(r["wall_seconds"] for r in a),
            cpu_seconds=median(r["cpu_seconds"] for r in a),
            proof_wall_seconds=median(r["online_screening"]["wall_seconds"] for r in a),
            proof_cpu_seconds=median(r["online_screening"]["cpu_seconds"] for r in a),
            numeric_full_wave_lu=median(r["ordinary_full_wave_counts"]["lu"] for r in a),
            physical_forward_rhs=median(r["ordinary_full_wave_counts"]["forward_rhs"] for r in a),
            physical_adjoint_rhs=median(r["ordinary_full_wave_counts"].get("adjoint_rhs",0) for r in a),
            refined_boxes=median(len(r["visited"]) for r in a),
            endpoint_interval_checks=median(r["endpoint_interval_checks"] for r in a))
    values["paired_shifted_over_unfiltered_wall"]=median(
        r["arms"][ARMS[1]]["wall_seconds"]/r["arms"][ARMS[0]]["wall_seconds"] for r in pair)
    values["paired_shifted_over_unfiltered_cpu"]=median(
        r["arms"][ARMS[1]]["cpu_seconds"]/r["arms"][ARMS[0]]["cpu_seconds"] for r in pair)
    return values


def main(output,repetitions):
    ctx.prec,ctx.threads=112,1
    overall=clock()
    g,fs,fb=geometry()
    s,b=fs[RX],fb[:,TX]
    cases=json.loads((OUT/"public/cases.json").read_text())
    with np.load(OUT/"predictions.npz",allow_pickle=False) as frozen:
        snapshots={r["id"]:(frozen[r["id"]+"__lower"],frozen[r["id"]+"__upper"]) for r in cases}
    prior=json.loads((OUT/"online.json").read_text())["cases"]
    setup=elapsed(overall)
    rng=np.random.default_rng(SEED)
    arrangements={}
    for case in cases:
        a=list(np.arange(repetitions)%2)
        rng.shuffle(a)
        arrangements[case["id"]]=a
    records={}

    # Do all timed replays before loading any material labels or held-out readings.
    valid_points={}
    for case in cases:
        name=case["id"]
        with np.load(OUT/"public"/f"{name}.npz",allow_pickle=False) as p:
            y,r=p["Y"],float(p["radius"])
        measured=[]
        snapshots_agree=[]
        # The next loop keeps arms separate so both pay the same candidate generation and ranking.
        for first in arrangements[name]:
            execution_order=list(ARMS if first==0 else ARMS[::-1])
            group, material, proofs={},{},{}
            for name_arm in execution_order:
                group[name_arm],material[name_arm],proofs[name_arm]=arm(
                    g,s,b,y,r,case["seed"]+1007,prune=(name_arm==ARMS[1]))
            without,withcut=[group[k] for k in ARMS]
            m0,m1=[material[k] for k in ARMS]
            assert np.array_equal(m0[0],m1[0]) and np.array_equal(m0[1],m1[1]) and np.array_equal(m0[2],m1[2])
            assert np.array_equal(m1[1],snapshots[name][0]) and np.array_equal(m1[2],snapshots[name][1])
            assert without["ranking"]==withcut["ranking"]
            assert without["ranked_sample_per_box"]==withcut["ranked_sample_per_box"]
            assert without["first_feasible"]==withcut["first_feasible"]
            assert without["first_feasible"] is False or without["reported_branch"]==withcut["reported_branch"]
            previously_excluded=[i for i,a in enumerate(prior[name]["methods"]["fast_boxes"]["decisions"]) if a["rejected"]]
            assert previously_excluded==withcut["excluded"]
            assert [a["branch"] for a in withcut["visited"]]==[
                   a["branch"] for a in without["visited"] if a["branch"] not in withcut["excluded"]]
            equality=dict(samples_both_arms=True,boxes_equal_frozen=True,
                          material_arrays_bitwise_equal=same_optional(m0[3],m1[3]),
                          train_prediction_bitwise_equal=same_optional(m0[4],m1[4]))
            snapshots_agree.append(equality)
            measured.append(dict(execution=execution_order,arms=group,**equality))
            if name not in valid_points:
                valid_points[name]=(y,r,m0,m1,proofs[ARMS[1]])
            else:
                assert same_optional(m0[3],valid_points[name][2][3])
                assert same_optional(m1[3],valid_points[name][3][3])

        values=measured_medians(measured)
        records[name]=dict(radius=r,repeats=measured,medians=values,
                           identical_reported_material_every_repeat=all(x["material_arrays_bitwise_equal"] for x in snapshots_agree),
                           identical_reported_training_prediction_every_repeat=all(x["train_prediction_bitwise_equal"] for x in snapshots_agree))
        print(f"{name}: visits {values[ARMS[0]]['refined_boxes']:g} → {values[ARMS[1]]['refined_boxes']:g}; "
              f"LU {values[ARMS[0]]['numeric_full_wave_lu']:g} → {values[ARMS[1]]['numeric_full_wave_lu']:g}; "
              f"wall {values[ARMS[0]]['wall_seconds']:.6f} → {values[ARMS[1]]['wall_seconds']:.6f}; "
              f"fit {measured[-1]['arms'][ARMS[1]]['first_feasible']}",flush=True)

    # Re-evaluate each frozen exclusion on the same exact stored arrays using Arb.
    for name,(y,r,m0,m1,decisions) in valid_points.items():
        start=clock()
        reference=Border(g,s,b)
        reference.measurement(y)
        verified=[dict(branch=i,**reference.check(m1[1][i],m1[2][i],y,r,verify=True,weights=WEIGHTS,jvalues=(1,)))
                  for i,result in enumerate(decisions) if result["rejected"]]
        complete=elapsed(start)
        assert all(v["rejected"] and v["arb_congruence_lower"]>0. for v in verified)
        with np.load(OUT/"scoring"/f"{name}.npz",allow_pickle=False) as labels:
            truth, clean=labels["truth"],labels["ports"]
        contains=np.all((m1[1]<=truth)&(truth<=m1[2]),axis=1)
        assert not any(contains[a["branch"]] for a in verified)
        records[name].update(strict_independent_recheck=dict(
            records=verified,whole_replay_cost=complete,
            arb_congruence_only=add_times([v["independent_check"] for v in verified])),
            truth_containing_boxes=np.flatnonzero(contains).tolist(),
            offline_material_and_ports={k:score(w[3],truth,clean,g,fs,fb) for k,w in zip(ARMS,(m0,m1),strict=True)},
            reported_active_material={k:w[3].tolist() if w[3] is not None else None
                                       for k,w in zip(ARMS,(m0,m1),strict=True)})

    med=[r["medians"] for r in records.values()]
    columns=("wall_seconds","cpu_seconds","numeric_full_wave_lu","physical_forward_rhs","physical_adjoint_rhs","refined_boxes")
    total={method:{field:sum(v[method][field] for v in med) for field in columns} for method in ARMS}
    summary=dict(sum_of_object_medians=total,
        ratio_of_sum_median_walls=total[ARMS[1]]["wall_seconds"]/total[ARMS[0]]["wall_seconds"],
        median_paired_object_ratio_wall=median(x["paired_shifted_over_unfiltered_wall"] for x in med),
        median_paired_object_ratio_cpu=median(x["paired_shifted_over_unfiltered_cpu"] for x in med),
        first_feasible={k:sum(v["repeats"][0]["arms"][k]["first_feasible"] for v in records.values()) for k in ARMS},
        choices_bitwise_equal=sum(x["identical_reported_material_every_repeat"] for x in records.values()),
        strict_exclusions_total=sum(len(x["strict_independent_recheck"]["records"]) for x in records.values()),
        truth_containing_boxes_total=sum(len(x["truth_containing_boxes"]) for x in records.values()),
        arb_congruence_only=add_times([x["strict_independent_recheck"]["arb_congruence_only"] for x in records.values()]),
        whole_independent_replay=add_times([x["strict_independent_recheck"]["whole_replay_cost"] for x in records.values()]))
    write_json(output,dict(cases=records,summary=summary,repetitions=repetitions,shared_geometry_read=setup,
                           total=elapsed(overall),scope="Retrospective, exposed confirmation inputs. All 64 points are ranked in both arms."
                           " Common within-box best point initialization and 42-call full-wave refinement; "
                           "first outward-verified feasible wins, or the numerically best audited endpoint is reported as a failure.",
                           order_seed=SEED))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=type(HERE),default=HERE/"same_stop.json")
    parser.add_argument("--repetitions",type=int,default=10)
    ns=parser.parse_args()
    if ns.repetitions<2 or ns.repetitions%2:
        parser.error("Use an even number of repetitions for equal paired execution orders.")
    with threadpool_limits(limits=1):
        main(ns.output,ns.repetitions)
