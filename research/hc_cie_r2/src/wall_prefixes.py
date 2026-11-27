"""Rescore v2 at recorded job-wall prefixes, using only sealed candidates."""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np

from imaging_run import LAB, json_write
from imaging_score import ARM_LABEL, add_result, get_case, load_sealed, read, write_csv
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator


PREFIXES = (.25, .5, 1., 1.5, 2.)
ARMS = ("continue", "restart", "single", "pair_random", "pair_born", "pair_secant", "pair_exact")
START_LABEL = dict(cie="CIE", dual_cie="Dual-CIE", a_plus_arc_r7="A+-Arc/R7", g19="G19")
WALL_BOUNDARY = (
    "Per-job perf_counter from entry to run_one in an already dispatched worker. "
    "The clock includes in-job imports, input reads, setup, physical solves, and scheduling delays "
    "until each stamp. The stamp is taken after the physical state returns, before audit residual/"
    "prior bookkeeping and candidate insertion, for which a separate timestamp is unavailable. "
    "Worker creation, shared interpreter/module startup, scheduler queue time before run_one, "
    "and final serialization are outside candidate-wall stamps. Existing CPU-paid flags still apply."
)
CONSERVATIVE_BOUNDARY = (
    "Use the next completed-state stamp as an upper bound on when the preceding candidate "
    "finished audit bookkeeping and became available. For the last candidate use recorded "
    "total_wall, which already includes writing PREDICTIONS.npz and excludes the final TRACE.json write."
)


def wall_rows(frozen):
    traces = list(load_sealed(frozen))
    truths = {case: get_case(case) for case in dict.fromkeys(t["case"] for t, _ in traces)}
    for trace, arrays in traces:
        truth = truths[trace["case"]]
        chain = trace["chain"]
        bounds = [entry["wall"] for entry in chain] + [trace["total_wall"]]
        assert all(lo <= hi for lo, hi in zip(bounds, bounds[1:]))
        start = add_result(trace, chain[0], arrays[0], truth, "endpoint")
        for cutoff in PREFIXES:
            base = dict(case=trace["case"], object=truth["object_id"], group=truth["group"],
                        endpoint=trace["endpoint"], arm=trace["arm"], wall_budget=cutoff,
                        expected_job=True, total_completed_full_audits_in_trace=len(chain),
                        first_full_audit_wall=chain[0]["wall"], original_endpoint_Phi=start["Phi"],
                        first_candidate_wall_upper=bounds[1],
                        original_endpoint_NRMSE=start["NRMSE"],
                        original_endpoint_image_pass=start["image_pass"],
                        original_endpoint_joint_pass=start["historical_joint_pass"])
            results = []
            for timestamps in (bounds[:-1], bounds[1:]):
                at_time = [entry for entry, stamp in zip(chain, timestamps) if stamp <= cutoff]
                eligible = [entry for entry in at_time if entry["paid"]]
                row = dict(**base, full_audits_by_wall=len(at_time), paid_full_audits_by_wall=len(eligible),
                           availability="AVAILABLE" if eligible else "NO_PAID_FULL_AUDIT_AT_WALL")
                if eligible:
                    selected = min(eligible, key=lambda entry: (entry["Phi"], entry["i"]))
                    row.update(add_result(trace, selected, arrays[selected["i"]], truth, "wall_prefix"))
                    row["candidate_wall_upper"] = bounds[selected["i"]+1]
                results.append(row)
            assert results[0]["paid_full_audits_by_wall"] - results[1]["paid_full_audits_by_wall"] in (0, 1)
            yield tuple(results)


def available(row):
    return row["availability"] == "AVAILABLE"


def summarize(rows, seal):
    coverage = []
    comparisons = []
    strata = ("all", "R13_common_failure", "R13_positive_control", "R6_historical_curved_bar")
    references = [(mode, "continue") for mode in ARMS if mode != "continue"]
    references += [(mode, reference) for mode in ("pair_secant", "pair_exact")
                   for reference in ("pair_born", "pair_random", "single")]
    for cap in PREFIXES:
        for stratum in strata:
            candidates = [r for r in rows if r["wall_budget"] == cap
                          and (stratum == "all" or stratum == r["group"])]
            for arm in ("all", *ARMS):
                asked = [r for r in candidates if arm == "all" or arm == r["arm"]]
                seen = [r for r in asked if available(r)]
                coverage.append(dict(wall_budget=cap, group=stratum, arm=arm, expected_jobs=len(asked),
                                     available_jobs=len(seen), missing_jobs=len(asked)-len(seen),
                                     zero_audits=sum(r["full_audits_by_wall"] == 0 for r in asked),
                                     endpoint_incumbents=sum(r["index"] == 0 for r in seen),
                                     new_incumbents=sum(r["index"] > 0 for r in seen),
                                     objects_with_available_jobs=len({r["object"] for r in seen}),
                                     new_material_passes=sum(r["image_pass"] and not r["original_endpoint_image_pass"] for r in seen),
                                     new_joint_passes=sum(r["historical_joint_pass"] and
                                                          not r["original_endpoint_joint_pass"] for r in seen)))
            matched = defaultdict(dict)
            for row in candidates:
                matched[row["case"], row["endpoint"]][row["arm"]] = row
            for arm, reference in references:
                objects = defaultdict(list)
                if_matched = 0
                absent_method = absent_reference = 0
                for pair in matched.values():
                    a, b = pair[arm], pair[reference]
                    ok_a, ok_b = available(a), available(b)
                    absent_method += not ok_a
                    absent_reference += not ok_b
                    if not (ok_a and ok_b):
                        continue
                    if_matched += 1
                    objects[a["object"]].append(dict(
                        endpoint=a["endpoint"], NRMSE_delta=a["NRMSE"]-b["NRMSE"],
                        IoU_delta=a["IoU"]-b["IoU"], log_Phi_ratio=np.log(a["Phi"]/b["Phi"]),
                        relative_port_error_delta=a["relative_port_error"]-b["relative_port_error"]))
                pairs = []
                fields = ("NRMSE_delta", "IoU_delta", "log_Phi_ratio", "relative_port_error_delta")
                for case, observations in sorted(objects.items()):
                    pairs.append(dict(material=case, endpoints=[p["endpoint"] for p in observations],
                                      **{key: float(np.mean([p[key] for p in observations])) for key in fields}))
                comparison = dict(wall_budget=cap, group=stratum, arm=arm, reference=reference,
                                  expected_endpoint_pairs=len(matched), available_endpoint_pairs=if_matched,
                                  method_missing_endpoints=absent_method, reference_missing_endpoints=absent_reference,
                                  expected_materials=len({r["object"] for r in candidates}),
                                  available_materials=len(pairs), objects=pairs)
                for key in fields:
                    comparison["mean_paired_"+key] = float(np.mean([x[key] for x in pairs])) if pairs else None
                comparison["objects_NRMSE_better"] = sum(x["NRMSE_delta"] < -1e-9 for x in pairs)
                comparison["objects_NRMSE_worse"] = sum(x["NRMSE_delta"] > 1e-9 for x in pairs)
                comparisons.append(comparison)
    return dict(
        source_seal=seal["id"], jobs=seal["n"], scores=len(rows),
        prefixes=list(PREFIXES), independent_materials=len({r["object"] for r in rows}),
        saved_candidates=sum(job["candidate_count"] for job in seal["jobs"]),
        per_job_wall_scope=WALL_BOUNDARY, historical_batch_wall=seal["batch_wall"],
        historic_workers=seal["workers"], hardware=seal["cpu_info"],
        aggregation="Average endpoint differences within each material, then average over materials with a paired endpoint. Missing candidates remain missing.",
        coverage=coverage, paired_by_material=comparisons,
    )


def two_second_bounds(brackets):
    """Bounds allow a different admissible completion time in each paired job."""
    terminal = {(early["case"], early["endpoint"], early["arm"]): (early, late)
                for early, late in brackets
                if early["wall_budget"] == 2. and early["group"] == "R13_common_failure"}
    pairs = [("pair_exact", reference) for reference in ("continue", "pair_born", "pair_random")]
    pairs.append(("restart", "continue"))
    rows = []
    for arm, reference in pairs:
        by_material = defaultdict(list)
        for (case, endpoint, mode), a in terminal.items():
            if mode != arm:
                continue
            b = terminal[case, endpoint, reference]
            if not all(available(row) for row in (*a, *b)):
                continue
            aa, bb = [r["NRMSE"] for r in a], [r["NRMSE"] for r in b]
            pp, qq = [r["Phi"] for r in a], [r["Phi"] for r in b]
            by_material[a[0]["object"]].append((
                min(aa) - max(bb), max(aa) - min(bb),
                np.log(min(pp) / max(qq)), np.log(max(pp) / min(qq)),
            ))
        ranges = np.mean([np.mean(values, axis=0) for values in by_material.values()], axis=0)
        rows.append(dict(arm=arm, reference=reference, materials=len(by_material),
                         endpoints=sum(map(len, by_material.values())),
                         mean_NRMSE_difference_interval=ranges[:2].tolist(),
                         geometric_mean_Phi_ratio_interval=np.exp(ranges[2:]).tolist()))
    return rows


def plot_failure_endpoints(rows, target):
    failures = [r for r in rows if r["group"] == "R13_common_failure"]
    ids = list(dict.fromkeys(r["case"] for r in failures))
    ends = list(dict.fromkeys(r["endpoint"] for r in failures))
    lines = defaultdict(list)
    for row in failures:
        lines[row["case"], row["endpoint"], row["arm"]].append(row)
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axes = plt.subplots(len(ids), len(ends), figsize=(15, 10.8), sharex=True,
                             squeeze=False, layout="constrained")
    for i, case in enumerate(ids):
        for j, endpoint in enumerate(ends):
            ax = axes[i, j]
            start = lines[case, endpoint, "continue"][0]["original_endpoint_NRMSE"]
            ax.axhline(start, color=".62", ls=":", lw=.9, zorder=1)
            for hue, arm in zip(color, ARMS):
                r = lines[case, endpoint, arm]
                ax.plot([a["wall_budget"] for a in r],
                        [a["NRMSE"] if available(a) else np.nan for a in r],
                        label=ARM_LABEL[arm], marker=".", lw=1.2, color=hue, ms=4)
            ax.grid(alpha=.22)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
            ax.set_title(f"{case[-4:]} · {START_LABEL[endpoint]}", fontsize=9)
            if i == len(ids)-1:
                ax.set_xlabel("Recorded job wall (s)")
            if j == 0:
                ax.set_ylabel("Amplitude NRMSE")
    axes[-1, 0].set_xticks(PREFIXES)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=7, fontsize=9)
    fig.suptitle("Five historical failures · 20 paired endpoints · wall-prefix incumbents\n"
                 "Dotted line marks the saved endpoint · Vertical scales vary by panel", fontsize=13)
    fig.savefig(target, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sealed", type=Path, default=LAB/"results/imaging/v2")
    p.add_argument("--out", type=Path, default=LAB/"results/imaging/score_v2")
    p.add_argument("--figure", type=Path, default=LAB/"figures/wall_prefixes.png")
    a = p.parse_args()
    seal = read(a.sealed/"PREDICTION_SEAL.json")
    brackets = list(wall_rows(a.sealed))
    rows, conservative = map(list, zip(*brackets))
    assert len(rows) == len(PREFIXES)*seal["n"]
    fields = ("availability", "full_audits_by_wall", "paid_full_audits_by_wall", "index",
              "wall", "candidate_wall_upper", "Phi", "relative_port_error", "NRMSE",
              "IoU", "image_pass", "historical_joint_pass")
    table = []
    for early, late in brackets:
        table.append(dict(**early, selected_index_changes=early.get("index") != late.get("index"),
                          **{f"conservative_{key}": late.get(key) for key in fields}))
    write_csv(a.out/"WALL_PREFIXES.csv", table)
    early_summary = summarize(rows, seal)
    late_summary = summarize(conservative, seal)
    early_summary.update(
        conservative_wall_scope=CONSERVATIVE_BOUNDARY,
        max_paid_candidates_between_bounds=max(
            early["paid_full_audits_by_wall"]-late["paid_full_audits_by_wall"]
            for early, late in brackets),
        conservative_coverage=late_summary["coverage"],
        conservative_paired_by_material=late_summary["paired_by_material"],
        mixed_eligibility_bounds_2s_common_failures=two_second_bounds(brackets),
    )
    early_summary["selected_index_changes"] = [
        dict(wall_budget=cap, jobs=sum(r["selected_index_changes"] for r in table if r["wall_budget"]==cap))
        for cap in PREFIXES
    ]
    json_write(a.out/"WALL_SUMMARY.json", early_summary)
    plot_failure_endpoints(rows, a.figure)


if __name__ == "__main__":
    main()
