#!/usr/bin/env python3
"""Join the four sealed Dual-CIE traces with reference and held-source scores."""

from collections import defaultdict
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from diagnose import (Austria, CASES, DATA, HERE, LA, SRAP, Shell, Slotted, Twin, engine,
                      get_case, read_npz, write)

TRAJECTORIES = HERE / "trajectories"
PROBES = ("optical_diag", "optical_offdiag", "optical_diag_norm", "optical_offdiag_norm",
          "optical_diag_rms", "optical_offdiag_rms", "optical_diag_pair_normalized_rms",
          "optical_offdiag_pair_normalized_rms", "sector_positive_Fro",
          "sector_pixel_positive_trace", "constitutive_relative",
          "aux_data_relative", "J_frob2")


def read(path):
    return json.loads(path.read_text())


def trace(fold):
    with (fold / "updates.jsonl").open() as fd:
        rows = [json.loads(line) for line in fd]
    dpcs, updates = rows[0], rows[1:]
    return dpcs, updates


def references(name, case, public):
    if name == Austria:
        return engine.legacy.coarse_truth(case, public["centers"]), "archived Austria model, 8×8 area quadrature per inverse cell"
    if name == Slotted:
        return read_npz(DATA/"evaluation/optfresh_01.npz")["m"], "archived 96→32 cell-average map"
    if name == Shell:
        return read_npz(SRAP/f"data/truth/{Shell}.npz")["m"], "archived SRAP 96→32 cell-average map"
    return None, "measured target; scored on nine held-out transmitters"


def error_values(materials, truth):
    error = materials - truth[None, :]
    foreground = truth > 0
    half = truth >= np.max(truth)/2
    chosen = materials >= np.max(truth)/2
    return {
        "relative_material_L2": np.linalg.norm(error, axis=1)/np.linalg.norm(truth),
        "inside_relative_L2": np.linalg.norm(error[:, foreground], axis=1)/np.linalg.norm(truth[foreground]),
        "background_RMSE": np.sqrt(np.mean(error[:, ~foreground]**2, axis=1)),
        "background_relative_L2_to_truth_norm": np.linalg.norm(error[:, ~foreground], axis=1)/np.linalg.norm(truth),
        "half_peak_IoU": np.sum(chosen & half[None, :], axis=1)/np.maximum(1, np.sum(chosen | half[None, :], axis=1)),
    }


def heldout(case, public, audit_steps, materials):
    data = read_npz(LA/f"heldout/{Twin}.npz")
    Y = data["Y"]
    saved_valid = np.isfinite(Y)
    model = engine.legacy.PaperModel(dict(E=data["E"], S=data["S"], centers=public["centers"],
                                          k=public["k"], radius=public["radius"]),
                                     case, engine.legacy.Clock(float("inf")))
    start = perf_counter(), process_time()
    out = {}
    for step in audit_steps:
        m = materials[step-1] if step else np.zeros(model.n)
        f = model.ports.receive(model.state(m)["j"])
        residual = float(np.linalg.norm((f-Y)[saved_valid]))
        out[step] = dict(residual=residual, relative=residual/np.linalg.norm(Y[saved_valid]))
    return out, dict(counts=dict(model.clock.counts),
                     source_count=Y.shape[1], complex_channels=int(np.sum(saved_valid)),
                     wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


def correlation(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    valid = np.isfinite(a) & np.isfinite(b)
    n = int(np.sum(valid))
    if n < 3 or np.ptp(a[valid]) == 0. or np.ptp(b[valid]) == 0.:
        return dict(n=n, rho=None)
    return dict(n=n, rho=float(spearmanr(a[valid], b[valid]).statistic))


def short_step(row):
    fields = ("step", "stage", "point", "optical_diag", "optical_offdiag", "optical_diag_rms",
              "optical_offdiag_rms", "optical_diag_norm", "optical_offdiag_norm",
              "optical_diag_pair_normalized_rms", "optical_offdiag_pair_normalized_rms",
              "optical_diag_energy_fraction", "optical_rank_relative_1e_9", "optical_operator_norm",
              "sector_positive_Fro", "sector_positive_rank", "sector_pixel_positive_trace",
              "constitutive_relative", "aux_data_relative", "m_peak")
    return {k: row[k] for k in fields}


def score_case(name):
    fold = TRAJECTORIES / name
    receipt = read(fold/"receipt.json")
    dpcs, replay = trace(fold)
    data = read_npz(fold/"materials.npz")
    if len(replay) != receipt["final_step"] or not np.array_equal(data["step"], list(range(1, len(replay)+1))):
        raise ValueError(f"Incomplete trace for {name}.")
    case, public = get_case(name)
    truth, reference = references(name, case, public)
    quality = error_values(data["materials"], truth) if truth is not None else {}
    if len(quality):
        initial = {k: float(v[0]) for k, v in error_values(np.zeros((1, len(truth))), truth).items()}
        np.savez_compressed(fold/"offline_reference.npz", m=truth, kind=reference)
    else:
        initial = {}
    train = {row["step"]: row for row in receipt["accepted_state_audits"]}
    held, held_costs = heldout(case, public, train, data["materials"]) if name == Twin else ({}, None)
    series = []
    for i, row in enumerate(replay):
        step = row["step"]
        clean = {k: row[k] for k in ("step", "stage", "beta_user_material", "beta_user_current",
                                     "engine_wall_seconds", "m_peak", "optical_rank_relative_1e_9", *PROBES)}
        clean.update({k: float(v[i]) for k,v in quality.items()})
        clean["physical_train_rel"] = train[step]["fit_relative"] if step in train else None
        clean["physical_held_source_rel"] = held[step]["relative"] if step in held else None
        series.append(clean)
    with (fold/"scored_trajectory.csv").open("w", newline="") as fd:
        w = csv.DictWriter(fd, fieldnames=list(series[0]))
        w.writeheader()
        w.writerows(series)
    stages = defaultdict(list)
    for j,row in enumerate(replay):
        stages[row["stage"]].append(j)
    outcomes = [*quality.keys(), "physical_train_rel"]
    if name == Twin:
        outcomes.append("physical_held_source_rel")
    assoc = {}
    for stage, indexes in [("all", list(range(len(replay)))), *stages.items()]:
        assoc[str(stage)] = {probe: {k: correlation([series[i][probe] for i in indexes],
                                                   [series[i][k] for i in indexes])
                                    for k in outcomes} for probe in PROBES}
    q_over_d = np.array([r["optical_offdiag_pair_normalized_rms"]/r["optical_diag_pair_normalized_rms"] for r in replay])
    raw_over_d = np.array([r["optical_offdiag"]/r["optical_diag"] for r in replay])
    finals = {k: series[-1][k] for k in outcomes}
    physical_norm = [r["physical_optical_offdiag_normalized"] for r in train.values()
                     if r["physical_optical_offdiag_normalized"] is not None]
    summary = dict(object_id=name, stage_steps={str(k): len(v) for k,v in stages.items()},
                   dpcs=short_step(dpcs), final=short_step(replay[-1]), quality_initial=initial,
                   quality_final=finals, ref_description=reference, correlations=assoc,
                   train_audits=list(train.values()), held_source_audits=held, held_costs=held_costs,
                   median_pairwise_normalized_offdiag_to_diag=float(np.median(q_over_d)),
                   median_raw_offdiag_to_diag=float(np.median(raw_over_d)),
                   fraction_pairwise_offdiag_larger=float(np.mean(q_over_d > 1.)),
                   fraction_extensive_offdiag_larger=float(np.mean(raw_over_d > 1.)),
                   maximum_physical_optical_offdiag_normalized=max(physical_norm, default=None),
                   optimizer_wall=series[-1]["engine_wall_seconds"],
                   measured_full_wall=receipt["measured_full_wall_seconds"])
    return summary


def main():
    for name in CASES:
        path = TRAJECTORIES / name / "receipt.json"
        if not path.exists() or read(path).get("truth_access") != "NONE":
            raise FileNotFoundError(f"Complete and seal every online trajectory first: {name}")
    seal = HERE/"online_seal.json"
    if not seal.is_file():
        write(seal, dict(cases=CASES, truth_access="NONE",
                         timestamp=datetime.now(timezone.utc).isoformat(),
                         mode="complete trajectories, material snapshots, actual field snapshots and receipts"))
    started = perf_counter(), process_time()
    out = {}
    for name in CASES:
        result = score_case(name)
        out[name] = result
        print(name, "initial/final", result["quality_initial"], result["quality_final"],
              "median raw/entry off-diagonal ratio", result["median_raw_offdiag_to_diag"],
              result["median_pairwise_normalized_offdiag_to_diag"], flush=True)
        write(HERE/"scores.json", dict(cases=out, wall_seconds=perf_counter()-started[0],
                                       cpu_seconds=process_time()-started[1]))


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
