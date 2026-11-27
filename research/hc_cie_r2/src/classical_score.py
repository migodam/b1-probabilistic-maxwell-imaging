"""Post-seal Austria scoring and an audit of held-out Fresnel sources."""
import argparse
import csv
import json
import os
from pathlib import Path
import time

from imaging_run import LAB, REPO, json_write
from imaging_score import read, write_csv
from classical_run import RUNTIME, TOP, inherited

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np


def read_seal(base):
    seal = read(base/"PREDICTION_SEAL.json")
    if not seal["no_truth_or_holdout_read"] or seal["expected"] != len(seal["jobs"]):
        raise ValueError("Require every paid prediction to be sealed")
    results = {}
    for row in seal["jobs"]:
        directory = base/row["directory"]
        trace = read(directory/"TRACE.json")
        with np.load(directory/"AUDITS.npz", allow_pickle=False) as z:
            arrays = {k:z[k] for k in z.files}
        actual = min((r for r in trace["chain"] if r["paid"]), key=lambda r:r["phi"])
        assert actual["index"] == int(arrays["selected_index"]) == row["selected"]
        assert trace["without_truth_or_holdout"]
        results[row["case"], row["parent"], row["arm"]] = (trace, arrays)
    return seal, results


def raw_audit(case, public, held):
    from research.research_max_struct_r8.src.fresnel_adapter import load_file

    raw_dir = Path(os.environ.get("FRESNEL_2005_DIR",
                                  REPO/"research/delegated/a2_literature/data/2005_iop_21_6_S09"))
    raw_file = raw_dir/(case["object"]+".exp")
    header_rows, raw_y, _ = load_file(raw_file, case["frequency"]/1e9)
    train = np.asarray(case["train_sources"])-1
    test = np.asarray(case["held_sources"])-1
    all_ids = np.r_[train, test]
    assert len(np.unique(all_ids)) == len(all_ids)
    rx = public["rx_ids"]
    training = np.take_along_axis(raw_y[:,train], rx.T, axis=0)
    expected = raw_y[:,test]
    valid = np.isfinite(expected)
    assert np.array_equal(np.isfinite(held["Y"]), valid)
    assert np.array_equal(training, public["Y"])
    assert np.array_equal(expected[valid], held["Y"][valid])
    assert len(header_rows) == training.size + np.count_nonzero(valid)
    return dict(dataset="Institut Fresnel FoamTwinDielTM", frequency_hz=case["frequency"],
                original_rows=len(header_rows), training_source_ids_1_based=(train+1).tolist(),
                held_source_ids_1_based=(test+1).tolist(), actual_training_ports=training.size,
                actual_held_source_ports=int(valid.sum()), all_actual_ports_disjoint=True,
                distinct_train_receiver_masks=len(set(map(tuple,rx.tolist()))),
                raw_train_value_max_abs_difference=float(np.max(np.abs(training-public["Y"]))),
                raw_held_value_max_abs_difference=float(np.max(np.abs(expected[valid]-held["Y"][valid]))),
                empty_field_calibration_median=case["calibration_median_validation"],
                raw_measurement_redistributed=False)


def score(base, out, figures):
    wall0,cpu0 = time.perf_counter(),time.process_time()
    seal, jobs = read_seal(base)
    top = inherited()
    out.mkdir(parents=True, exist_ok=False)
    summary = []
    saved_pred = {}
    audit_receipts = {}
    cost_rows = []
    with (TOP/"metrics.csv").open() as f:
        old = {(r["case"],r["method"],r["key"]):r for r in csv.DictReader(f)}
    for case in dict.fromkeys(k[0] for k in jobs):
        selected = {k[1:]:v for k,v in jobs.items() if k[0] == case}
        public = top.array_file(RUNTIME/"public"/f"{case}.npz")
        meta = json.loads(str(public["case_json"]))
        measured = "rx_ids" in public
        cache = []
        if measured:
            held = top.array_file(RUNTIME/"heldout"/f"{case}.npz")
            audit_receipts[case] = raw_audit(meta,public,held)
            hmodel = top.PaperModel(dict(public,E=held["E"],S=held["S"],rx_ids=None),meta,top.Clock(np.inf))
            valid = np.isfinite(held["Y"])
            holdnorm = np.linalg.norm(held["Y"][valid])
            saved_pred[f"{case}__observed_held_mask"] = valid
        else:
            ref = top.coarse_truth(meta, public["centers"])

        def get_metrics(m, name):
            if not measured:
                return top.errors(m,ref)
            found = next((p for q,p in cache if np.array_equal(q,m)),None)
            if found is None:
                t,c = time.perf_counter(),time.process_time()
                before = dict(hmodel.clock.counts)
                found = hmodel.state(m)["Y"].copy()
                cache.append((m.copy(),found))
                cost_rows.append(dict(case=case, output=name, postseal_scoring_CPU=time.process_time()-c,
                                      postseal_scoring_wall=time.perf_counter()-t,
                                      full_LU=hmodel.clock.counts["LU"]-before.get("LU",0),
                                      RHS=hmodel.clock.counts["complex_RHS"]-before.get("complex_RHS",0)))
            saved_pred[name] = found
            return {"heldout_relative_complex_error":float(np.linalg.norm((found-held["Y"])[valid])/holdnorm)}

        for parent, arm in selected:
            trace, arrays = selected[parent,arm]
            training = public["Y"]
            best = trace["best"]
            original = old[case,parent,trace["historical_receipt"]["key"]]
            assert np.isclose(trace["chain"][0]["relative_port_error"],float(original["relative_data_residual"]))
            for port,index in [(arrays["endpoint_port"],0),(arrays["selected_port"],best)]:
                assert np.isclose(np.linalg.norm(port-training)/np.linalg.norm(training),
                                  trace["chain"][index]["relative_port_error"],rtol=1e-12)
            jump_trace = selected[parent,"finite_two_pixels"][0]
            finite_end = min(jump_trace["budget_cpu"],jump_trace["online_calculation_cpu"])
            paired = [r for r in trace["chain"] if r["paid"] and r["cpu"]<=finite_end]
            stages = [("endpoint",trace["chain"][0]),("budget_selected",trace["chain"][best])]
            if paired:
                stages.append(("finite_actual_cpu_prefix",min(paired,key=lambda r:r["phi"])))
            for name,row in stages:
                key = f"{case}__{parent}__{arm}__{name}"
                metrics = get_metrics(arrays["material"][row["index"]],key)
                if measured and name=="endpoint":
                    assert np.isclose(metrics["heldout_relative_complex_error"],
                                      float(original["heldout_relative_complex_error"]),rtol=1e-7)
                summary.append(dict(case=case, parent=parent, arm=arm, source=name, selected_index=row["index"],
                                    comparison_cap_cpu=trace["budget_cpu"], finite_actual_prefix_cpu=finite_end,
                                    candidate_cpu=row["cpu"], observed_relative_error=row["relative_port_error"],
                                    Phi=row["phi"], online_total_cpu=trace["total_cpu"],
                                    online_excluding_write_cpu=trace["online_calculation_cpu"],
                                    online_wall=trace["total_wall"], full_LU=trace["full_counts"]["LU"],
                                    full_rhs=trace["full_counts"]["complex_RHS"], original_cpu_s=original["cpu_s"],
                                    original_wall_s=original["wall_s"], original_LU=original["LU"],
                                    original_Green_FFT=original["FFT"], original_maximum_wall_budget=original["budget"],
                                    original_ref_data_residual=original["reference_data_residual"],
                                    original_heldout_scoring_wall_s=original["heldout_scoring_wall_s"],
                                    stop=trace["stop"], **metrics))
    write_csv(out/"OUTCOMES.csv",summary)
    write_csv(out/"POSTSEAL_PHYSICS_COSTS.csv",cost_rows)
    json_write(out/"RAW_PORT_ALIGNMENT.json",audit_receipts)
    np.savez_compressed(out/"HELDOUT_PREDICTIONS.npz", **saved_pred)
    figures.mkdir(parents=True,exist_ok=True)
    plot(figures, summary)
    json_write(out/"SCORING_RECEIPT.json",dict(source=str(base), jobs=seal["expected"],
               post_seal_only=True, whole_process_CPU=time.process_time()-cpu0,
               whole_process_wall=time.perf_counter()-wall0))


def plot(figures, rows):
    shows = [r for r in rows if "heldout_relative_complex_error" in r and r["source"]=="budget_selected"]
    fig,axes = plt.subplots(1,2,figsize=(8.0,3.25),layout="constrained")
    for ax,(field,caption) in zip(axes,[("observed_relative_error","13 training sources"),
                                       ("heldout_relative_complex_error","5 unseen sources")]):
        names = ["Paper-NIE","Dual-CIE"]
        parents = ["paper_nie","cie_dual"]
        for i,(arm,label) in enumerate([("finite_two_pixels","Two pixels"),("continue","Continuation")]):
            expected = [next(r for r in shows if r["arm"]==arm and r["parent"]==parent)[field] for parent in parents]
            ax.bar(np.arange(2)+(i-.5)*.32,expected,.30,label=label)
        ax.set_xticks([0,1],names)
        ax.set_ylabel("Relative complex-field error")
        ax.set_title(caption)
        ax.grid(axis="y",alpha=.3)
    axes[0].legend(fontsize=8)
    fig.suptitle("FoamTwinDielTM 8 GHz · 4 CPU-s comparison")
    fig.savefig(figures/"classical_foam_ports.png",dpi=135)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sealed",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--figures",type=Path,default=LAB/"figures")
    a=parser.parse_args()
    score(a.sealed,a.out,a.figures)


if __name__ == "__main__":
    main()
