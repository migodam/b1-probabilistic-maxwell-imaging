"""Freeze methods, then audit sealed truth, evidence, and paired material outcomes."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import json
from time import perf_counter

from common import BGS, PAPER, PRIVATE, PUBLIC, RATIOS, REGION, RESULTS, SEAL
from common import check_preregistration, digest, load, write_json
import numpy as np
from bounds import Acquisition

DECISIONS = RESULTS / "DECISION_SEAL.json"
METHODS = "ABCDE"
NAMES = ("A directional-Neumann", "B spatial/local/matrix PSD duals",
         "C native inertia", "D identical compressed inertia", "E shared material rows")


def files_to_seal():
    for k in range(1,25):
        stem = f"pair_{k:02d}"
        for path in (RESULTS / "raw" / f"{stem}.json", RESULTS / "raw" / f"{stem}_arrays.npz",
                     RESULTS / "replay" / f"{stem}.json"):
            yield path


def seal_decisions():
    if DECISIONS.exists():
        raise FileExistsError("Decision manifest already sealed.")
    content = dict()
    for path in files_to_seal():
        if path.parent.name == "replay":
            raw = json.loads(path.read_text())
            if not raw["pass_all"]:
                raise ArithmeticError(f"Independent replay failed: {path}")
        elif path.suffix==".json":
            raw = json.loads(path.read_text())
            if len(raw["rows"])!=112:
                raise ArithmeticError(f"Incomplete pair: {path}")
        content[str(path.relative_to(RESULTS))] = digest(path)
    write_json(DECISIONS, dict(freeze_utc=datetime.now(timezone.utc).isoformat(),
                               files=content, pairs=24, proof_records_before_unseal=True))


def historical_controls():
    repo = PAPER.parent
    v5 = repo / "ip_final_v5" / "strong_control"
    older = repo / "ip_unified" / "numerics" / "feasibility" / "strict_endpoints"
    name = "control_disk075__pw360_32x32__pair2pct"
    low = load(older / f"{name}.npz")
    ring = load(v5 / "strong_ring_arrays.npz")
    exact = json.loads((PAPER / "results" / "native_replay.json").read_text())
    illum = ring["ids"]
    cases = list(exact) if isinstance(exact,list) else exact.get("cases",[])
    strict_replays = {}
    for r in cases:
        cap = r.get("cap") or {}
        strict_replays[r["id"]] = dict(nu=cap.get("nu"),
                                      distance_over_radius_lower=cap.get("distance_over_radius_lower"),
                                      legal_endpoint_rho_lower=cap.get("homogeneous_endpoint_moment",{}).get("rho_lower"),
                                      fits_original_ball=r.get("witness",{}).get("fits_noise_ball"))
    return dict(grid="96x96", common_candidate_allowed_cells=len(illum),
                candidate_center_actual_nonzero=int(np.count_nonzero(ring["c_solid"])),
                low_disk_actual_nonzero=int(np.count_nonzero(low["material"])),
                strong_ring_actual_nonzero=int(np.count_nonzero(ring["chi_active"])),
                low_disk_actual_nonzero_outside_allowed_region=int(np.count_nonzero(
                    np.setdiff1d(np.flatnonzero(low["material"]),illum))),
                actual_nonzero_support_overlap_low_disk_with_ring=int(np.count_nonzero(
                    (low["material"][illum]!=0) & (ring["chi_active"]!=0))),
                same_actual_nonzero_support_as_candidate=False,
                native_replay_controls=strict_replays,
                status="exposed historical same allowed-region controls; not prospective same nonzero support repetitions")


def material_digest(x):
    import hashlib
    return hashlib.blake2b(np.ascontiguousarray(x).tobytes(),digest_size=16).hexdigest()


def status_record(rows):
    record = {}
    for method in METHODS:
        c = Counter(row["status"][method] for row in rows)
        record[method] = {name:c[name] for name in ("EXCLUDED","WITNESS_FOUND","UNRESOLVED")}
    return record


def strict_ratios(row, method):
    r = row["methods"][method]
    if not r["excluded"]:
        return None
    if method=="A":
        return r["strict_distance_over_radius_lower"]
    if method=="E":
        if r["stage"]=="inherited_D":
            return strict_ratios(row,"D")
        return r["strict_distance_over_radius_lower"]
    return r["proof"]["strict_distance_over_radius_lower"]


def summary():
    began = perf_counter()
    check_preregistration()
    seal = json.loads(DECISIONS.read_text())
    for name, digest_string in seal["files"].items():
        assert digest(RESULTS/name)==digest_string
    data_seal = json.loads(SEAL.read_text())
    assert digest(PRIVATE)==data_seal["files"][PRIVATE.name]
    assert digest(PUBLIC)==data_seal["files"][PUBLIC.name]
    a, gt = load(PUBLIC),load(PRIVATE)
    records, false_exclusions, truth_failures = [], [], []
    joint_peak_hulls, results, replays = [], [], []
    charge, evidence_times = defaultdict(Counter), Counter()
    setup, between, checks, initial_failures = Counter(), Counter(), Counter(), Counter()
    witness_sources, pair_scored = Counter(), []
    limits = Counter()
    max_rss = 0
    true_forward_ratio = []
    all_bvals = defaultdict(lambda: defaultdict(list))
    for pair in range(24):
        suffix = f"pair_{pair+1:02d}.json"
        pair_rec = json.loads((RESULTS/"raw"/suffix).read_text())
        rep = json.loads((RESULTS/"replay"/suffix).read_text())
        assert rep["pass_all"]
        replays.append(rep)
        for name,value in rep["checked"].items():
            checks[name] += value
        for name,value in rep.get("seconds_by_method",{}).items():
            evidence_times[name] += value
        setup["wall"] += pair_rec["elapsed_wall_seconds"]
        setup["trace"] += sum(row["seconds"] for row in pair_rec["shared_trace_preprocess"])
        max_rss = max(max_rss,pair_rec["peak_resident_memory_bytes"])
        truth, branch = gt["material"][pair],gt["generating_branch"][pair]
        family = str(a["families"][pair])
        for bg, b in enumerate(BGS):
            active = [m != 0 for m in truth[:,bg]]
            same = np.array_equal(*active)
            if same!=(family=="amplitude"):
                truth_failures.append(dict(pair=pair+1,bg=b,reason="incorrect actual support class"))
            if family=="amplitude":
                field1, field2 = truth[:,bg][:,REGION]
                if np.linalg.norm(field1*max(field2)-field2*max(field1))>1e-12:
                    truth_failures.append(dict(pair=pair+1,bg=b,reason="non-amplitude variation inside region"))
        acquisitions = tracker = {}
        by_condition = defaultdict(list)
        for row in pair_rec["rows"]:
            t,bg,w,slt = row["box_positions"]
            key = (t,bg)
            if key not in acquisitions:
                acq = Acquisition(a["G"],a["S"],a["B"],a["Y"][pair,t,bg],float(a["radii"][pair,t,bg]))
                acquisitions[key] = acq
                original_truth_verdict = acq.strict_witness(truth[t,bg],np.zeros(acq.n),np.full(acq.n,32.))
                if not original_truth_verdict["proved"]:
                    truth_failures.append(dict(pair=pair+1,target=t,bg=bg,
                                               reason="generating material does not fit declared data ball"))
                if "residual_upper" in original_truth_verdict:
                    true_forward_ratio.append(original_truth_verdict["residual_upper"]/acq.radius)
                setup["original_acquisition"] += row["acquisition_setup_seconds"]
                setup["reduced_numeric"] += row["reduced_base_numeric_seconds"]
                setup["reduced_ball"] += row["reduced_base_ball_seconds"]
            ell, upp = a["boxes"][pair,bg,w,slt]
            covered = bool(np.all(ell<=truth[t,bg]) and np.all(truth[t,bg]<=upp))
            generating = bool(slt==branch[t])
            if generating and not covered:
                truth_failures.append(dict(pair=pair+1,target=t,bg=bg,w=w,slot=slt,
                                           reason="generating candidate no longer covers material"))
            for method in METHODS:
                raw = row["methods"][method]
                charge[method]["proposal"] += raw.get("proposal_seconds",0.)
                charge[method]["strict"] += raw.get("strict_seconds",0.)
                charge[method]["reused"] += raw.get("reused_full_proposal_seconds",0.)
                if row["status"][method]=="EXCLUDED" and covered:
                    false_exclusions.append(dict(pair=pair+1,row=row["box_positions"],method=method))
                if method=="E":
                    for attempt in raw.get("attempts", []):
                        for leaf in attempt["leaves"]:
                            limits[leaf.get("reason","leaf_certified")] += 1
            witness_sources[row["witness"].get("search",{}).get("stop_reason","failed")] += 1
            setup["witness_proposal"] += row["witness"]["proposal_seconds"]
            setup["witness_strict"] += row["witness"]["strict_seconds"]
            between[row["status"]["C"],row["status"]["D"]] += 1
            for trial in row["methods"]["B"].get("attempted",[]):
                if trial.get("excluded"):
                    all_bvals[family][trial["variant"]].append((pair,t,bg,w,slt))
            by_condition[t,bg,w].append(row)
            records.append(row | dict(independent_pair=pair+1,family=family,
                                      generating=generating, contains_generating_material=covered,
                                      material_hash=material_digest(truth[t,bg]),
                                      strict_ratios={method:strict_ratios(row,method) for method in METHODS}))
        for (t,bg,w), two in by_condition.items():
            hulls = {}
            for method in METHODS:
                allowed = [row for row in two if row["status"][method]!="EXCLUDED"]
                if not allowed:
                    hulls[method] = None
                    continue
                low = min(row["candidate_peak_lower"] for row in allowed)
                high = max(row["candidate_peak_upper"] for row in allowed)
                hulls[method] = [low,high,high-low]
            joint_peak_hulls.append(dict(pair=pair+1,target=t,bg_idx=bg,w_idx=w,
                                        family=family, truth_peak=float(max(truth[t,bg,REGION])),
                                        hull=hulls))
    if false_exclusions or truth_failures:
        write_json(RESULTS/"SCORE_FAILURE.json",{"false_exclusions":false_exclusions,
                                                 "truth_failures":truth_failures})
        raise ArithmeticError("Sealed material coverage audit failed")
    wrong = [row for row in records if not row["generating"]]
    against = [row for row in wrong if not row["contains_generating_material"]]
    matched_pairs = []
    for pair in range(1,25):
        group = [row for row in against if row["independent_pair"]==pair]
        pair_stat = status_record(group)
        matched_pairs.append(dict(pair=pair,family=group[0]["family"],true_opposite_noncover_boxes=len(group),
                                  e_minus_c=pair_stat["E"]["EXCLUDED"]-pair_stat["C"]["EXCLUDED"],
                                  e_minus_b=pair_stat["E"]["EXCLUDED"]-pair_stat["B"]["EXCLUDED"],
                                  states=pair_stat))
    phase = []
    for family in ("amplitude","structure"):
        for bg,b in enumerate(BGS):
            for w, width in enumerate(RATIOS):
                group = [row for row in against if row["family"]==family and row["bg_idx"]==bg
                         and row["width_idx"]==w]
                hulls = [row for row in joint_peak_hulls if row["family"]==family and row["bg_idx"]==bg
                         and row["w_idx"]==w]
                phase.append(dict(family=family,background=b,ratio=width,opposite_noncover_boxes=len(group),
                                  states=status_record(group), peak_hull_median={
                                    m:float(np.median([r["hull"][m][2] for r in hulls])) for m in METHODS}))
    unit_setup = (setup["reduced_numeric"]+setup["reduced_ball"])/192
    constant_trace = (setup["trace"]/
                      max(1,sum(bool(p["shared_trace_preprocess"]) for p in (json.loads(
                          (RESULTS/"raw"/f"pair_{k:02d}.json").read_text()) for k in range(1,25)))))
    amortization = {}
    for n in (100,1000,2688):
        # A stabilized projection on mean measured per-node work, with one
        # common bank and mean observation reuse of 14 boxes.
        amortization[str(n)] = {}
        for m in ("C","D"):
            total = n*(charge[m]["proposal"]+charge[m]["strict"])/2688
            if m=="D":
                total += constant_trace+np.ceil(n/14)*unit_setup
            amortization[str(n)][m] = float(total)
    diff = [
        dict(pair=row["independent_pair"],condition=row["box_positions"],
             family=row["family"], background=row["background_cap"], width=row["width_ratio"],
             C=row["status"]["C"],D=row["status"]["D"],
             numeric_tail_abs_error=row["eligible_C_D_max_numeric_difference"])
        for row in records if row["status"]["C"]!=row["status"]["D"]]
    qualifiers = {
        "inertia_positive_K_only": sum(row["status"]["C"]=="EXCLUDED" and
                                            row["methods"]["C"]["proof"]["numeric_K_negative_count"]==0
                                            for row in records),
        "inertia_with_negative_state_directions": sum(row["status"]["C"]=="EXCLUDED" and
                                            row["methods"]["C"]["proof"]["numeric_K_negative_count"]>0
                                            for row in records),
    }
    extremes = {m: dict(median=float(np.median([strict_ratios(row,m) for row in records
                                               if row["status"][m]=="EXCLUDED"])))
                if any(row["status"][m]=="EXCLUDED" for row in records) else None for m in METHODS}
    added = [
        dict(pair=row["independent_pair"],condition=row["box_positions"],family=row["family"],
             background=row["background_cap"], width=row["width_ratio"],
             q=row["methods"]["E"]["rows_at_first_exclusion"],
             ratio=row["methods"]["E"]["strict_distance_over_radius_lower"],
             method_status=row["status"], material_hash=row["material_hash"])
        for row in records if row["methods"]["E"]["stage"]=="added_physical_rows"]
    deliver = dict(contract_id="fixed-tm-A-25-prospective-v1", pairs=24, materials=48, acquisitions=192,
                   own_material_audits=192, generating_box_paths=1344,
                   truth_contained_box_rows=sum(row["contains_generating_material"] for row in records),
                   maximum_generating_truth_residual_over_radius_upper=max(true_forward_ratio),
                   false_exclusions=false_exclusions, truth_failures=truth_failures, actual_boxes=len(records),
                   full_states=status_record(records), opposite_slot_total=len(wrong),
                   opposite_noncover_total=len(against), opposite_noncover_states=status_record(against),
                   C_D_state_agreement=1-len(diff)/len(records), C_D_disagreements=diff,
                   C_D_max_numeric_trial_ratio_abs_error=max(
                       r["eligible_C_D_max_numeric_difference"] for r in records
                       if r["eligible_C_D_max_numeric_difference"] is not None),
                   E_strict_added=added, E_interval_leaf_diagnostics=dict(limits),
                   certified_neumann_divergence_examples=sum(len(r["confirmed_positive_endpoint_divergence_examples"])
                                                             for r in replays),
                   certified_inertia_negative_count_proposals=qualifiers,
                   strict_exclusion_distance_over_radius=extremes,
                   paired_independent_cluster_totals=matched_pairs, phase=phase,
                   witness_search_termination=dict(witness_sources),
                   B_subvariant_unique_exclusions={family:{variant:len(set(values)) for variant,values in metrics.items()}
                                                  for family,metrics in all_bvals.items()},
                   raw_run_timing_seconds={"method_components":dict(charge),"shared":dict(setup)},
                   projected_wall_amortization_seconds=amortization,
                   measured_process_peak_resident_bytes=max_rss,
                   independent_replay_counts=dict(checks), independent_replay_seconds_by_method=dict(evidence_times),
                   independent_replay_total_seconds=sum(r["seconds"] for r in replays),
                   historical_controls=historical_controls(),
                   score_seconds=perf_counter()-began)
    write_json(RESULTS/"summary.json",deliver)
    write_json(RESULTS/"peak_hulls.json",joint_peak_hulls)
    with (RESULTS/"all_nodes.csv").open("w",newline="") as handle:
        columns = ["pair_id","target_id","material_hash","family","condition_id","target_idx",
                   "box_idx","background_cap","width_ratio","generating","contains_generating_material",
                   "candidate_peak_lower","candidate_peak_upper",*METHODS,*[f"margin_{m}" for m in METHODS]]
        out = csv.DictWriter(handle,fieldnames=columns,extrasaction="ignore")
        out.writeheader()
        for row in records:
            margins = {f"margin_{m}":row["strict_ratios"][m] for m in METHODS}
            out.writerow(row | row["status"] | margins)
    print(json.dumps({k:deliver[k] for k in ("actual_boxes","full_states","opposite_noncover_states",
                                            "C_D_state_agreement","E_strict_added","independent_replay_counts")},
                     indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=["seal","score"])
    args=parser.parse_args()
    if args.stage=="seal":
        seal_decisions()
    else:
        summary()


if __name__=="__main__":
    main()
