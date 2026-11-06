"""Use already sealed endpoint norms to diagnose constitutive-to-port transport.

No arrays, truth, new forward solves or candidate selection are read/performed.
The numerical interval is not a certified operator-norm or discretization bound.
"""
import time
START_CPU, START_WALL = time.process_time(), time.perf_counter()
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT / "RESULTS/original_vector3d/NINE_OFFLINE_SCORES.json"
    stored = json.loads(source.read_text())
    rows = []
    for row in stored["primary_rows"]:
        physical = math.sqrt(row["Phi"])
        auxiliary_data = math.sqrt(row["final_auxiliary_data"])
        auxiliary_state = math.sqrt(row["final_auxiliary_state"])
        rows.append(dict(
            cohort=row["cohort"], arm=row["arm"],
            physical_data_RMS=physical, auxiliary_data_RMS=auxiliary_data,
            ordinary_current_state_RMS=auxiliary_state,
            transported_port_norm_lower=max(0., physical-auxiliary_data),
            transported_port_norm_upper=physical+auxiliary_data,
            realized_transport_gain_lower=(max(0., physical-auxiliary_data)/auxiliary_state
                                          if auxiliary_state else None),
            realized_transport_gain_upper=((physical+auxiliary_data)/auxiliary_state
                                          if auxiliary_state else None)))
    output = ROOT / "RESULTS/original_residual_transport"
    output.mkdir(exist_ok=True)
    result = dict(
        status="POSTHOC_DIAGNOSTIC_FROM_ALREADY_PAID_NORMS",
        source=str(source.relative_to(ROOT)), rows=rows,
        identity="r_physical = d_auxiliary - S L(m)^(-1) e; e=L(m)J_aux-D_m E",
        data_norm="sqrt(mean_s ||Z_s||^2/||Y_s||^2)",
        state_norm="sqrt(mean_s ||p(m)e_s||^2/||E_s||^2)",
        interpretation="Interval for the gain on this realized residual direction, under the finite-model identity. Its upper endpoint is NOT an operator-norm upper bound.",
        numerical_limit="Uses computed physical and auxiliary norms; no rigorous linear-solve/discretization error enclosure, so numerical values are diagnostic rather than certified bounds.",
        novelty="Existing Maxwell defect-transport identity plus triangle inequalities, used to interpret new endpoint evidence; not a new theorem.",
        truth_access=False, large_array_access=False, new_physical_solves=0,
        changes_to_predictions_or_selection=False)
    (output / "DIAGNOSTIC.json").write_text(json.dumps(result, indent=2)+"\n")
    with (output / "DIAGNOSTIC.csv").open("w", newline="") as f:
        writer=csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    receipt=dict(platform="Mac", CPU_seconds=time.process_time()-START_CPU,
                 wall_seconds=time.perf_counter()-START_WALL,
                 scope="One metadata-only residual-transport analysis, including imports and exports except receipt tail",
                 new_physical_solves=0, truth_access=False, additive=True, included_in=None)
    (output / "COST.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(rows))


if __name__ == "__main__":
    main()
