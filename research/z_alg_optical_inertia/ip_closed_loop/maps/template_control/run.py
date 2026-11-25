"""Fit two physical template coefficients at five pre-ranked candidate centres."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parent))

from common import FREEZE, HERE, Geometry, clock, elapsed, geom_path, pattern, place, save_json
from sector_map import SectorGeometry
from witness import certify_witness

B = FREEZE["model"]["background_upper"]
LOW, HIGH = FREEZE["material_classes"]["primary_core_range"]
CANDIDATE_COUNT = 5
ANCHORS = (0., B / 2, B)


def native_point(g, cells, profile, amplitude, background, measured):
    chi = np.full(g.n, background)
    chi[cells] = amplitude * profile
    predicted, current, code = g.forward(chi)
    if code:
        raise ArithmeticError(f"Physical forward failed ({code})")
    field = g.B + g.apply(current)
    state = current - chi[:, None] * field
    error = predicted - measured
    return dict(material=chi, current=current, predicted=predicted, error=error,
                residual=float(np.linalg.norm(error)),
                state_eq_residual=float(np.linalg.norm(state)),
                amplitude=float(amplitude), background=float(background))


def fit_at_centre(g, cells, profile, measured, radius):
    started = clock()
    cache = {}

    def query(z):
        physical = (LOW + (HIGH - LOW) * float(z[0]), B * float(z[1]))
        if physical not in cache:
            cache[physical] = native_point(g, cells, profile, *physical, measured)
        return cache[physical]

    anchors = [query((0.5, q / B)) for q in ANCHORS]
    initial = min(anchors, key=lambda q: q["residual"])

    def residual(z):
        difference = query(z)["error"]
        return np.r_[difference.real.ravel(), difference.imag.ravel()] / radius

    optimum = least_squares(
        residual, [(initial["amplitude"]-LOW)/(HIGH-LOW), initial["background"]/B],
        bounds=(np.zeros(2), np.ones(2)), jac="3-point", diff_step=1e-4,
        xtol=1e-9, ftol=1e-9, gtol=1e-9, max_nfev=24,
    )
    selected = min([initial, query(optimum.x)], key=lambda q: q["residual"])
    row = dict(
        anchor_residual_over_radius=[point["residual"] / radius for point in anchors],
        selected_residual_over_radius=selected["residual"] / radius,
        amplitude=selected["amplitude"], background=selected["background"],
        full_forward_evaluations=len(cache), optimizer_exit=int(optimum.status),
        optimizer_optimality=float(optimum.optimality), optimizer_residual_evaluations=optimum.nfev,
        max_raw_state_equation_residual=max(point["state_eq_residual"] for point in cache.values()),
        **elapsed(started),
    )
    return anchors, selected, row


def fit():
    began = clock()
    by_side = {}
    setup = {}
    for side in (32, 96):
        started = clock()
        by_side[side] = Geometry(geom_path(side))
        setup[str(side)] = elapsed(started)
    roster = json.loads((HERE / "public" / "acquisition.json").read_text())["observations"]
    selected_dir = OUT / "selected"
    selected_dir.mkdir(exist_ok=True)
    records = {}
    for oid, meta in roster.items():
        if meta["scene"] != "weak":
            continue
        started = clock()
        g = by_side[meta["side"]]
        with np.load(HERE / "public" / f"{oid}.npz", allow_pickle=False) as obs:
            measured, radius = obs["Y"], float(obs["r"])
        rank_source = HERE / "results" / f"{oid}__b{B:g}__primary.npz"
        with np.load(rank_source, allow_pickle=False) as old:
            errors = old["reference_error"]
        if len(errors) != g.n or np.any(~np.isfinite(errors)) or np.any(errors < 0):
            raise ArithmeticError(f"Missing finite census ranking for {oid}")
        shortlist = np.argsort(errors, kind="stable")[:CANDIDATE_COUNT]
        offsets, template = pattern(g.side, meta["shape"])
        evaluated = []
        for node in shortlist:
            centre = g.xy[node]
            cells, profile = place(g.side, centre, offsets, template)
            anchors, best, row = fit_at_centre(g, cells, profile, measured, radius)
            if not np.isclose(errors[node], anchors[0]["residual"], rtol=1e-8, atol=1e-10):
                raise ArithmeticError(f"Zero-background nominal residual changed for {oid}")
            row.update(centre=centre.tolist(), nominal_rank_in_shortlist=len(evaluated)+1)
            evaluated.append((row, best, cells, profile))
        chosen, result, cells, profile = min(evaluated, key=lambda value: value[1]["residual"])
        chosen_anchor = min(
            ((item[0]["anchor_residual_over_radius"][q], q, i)
             for i, item in enumerate(evaluated) for q in range(len(ANCHORS))),
            key=lambda value: value[0],
        )
        best_anchor_ratio, anchor_index, centre_index = chosen_anchor
        file = selected_dir / f"{oid}.npz"
        np.savez_compressed(
            file, material=result["material"], current=result["current"],
            prediction=result["predicted"], active_cells=cells, active_profile=profile,
            centre=np.array(chosen["centre"]), amplitude=result["amplitude"],
            background=result["background"],
        )
        records[oid] = dict(
            side=g.side, provided_template=meta["shape"], observations=measured.shape,
            noise_radius=radius, rank_source=str(rank_source.relative_to(HERE)),
            nominal_top5=shortlist.tolist(), selected_material=str(file.relative_to(OUT)),
            selected=chosen, nominal_best_residual_over_radius=float(errors[shortlist[0]]/radius),
            anchors_only_selected=dict(centre=evaluated[centre_index][0]["centre"],
                                       residual_over_radius=best_anchor_ratio,
                                       background=ANCHORS[anchor_index]),
            candidates=[v[0] for v in evaluated], total_forward_evaluations=sum(
                v[0]["full_forward_evaluations"] for v in evaluated
            ), **elapsed(started),
        )
        print(oid, chosen["centre"], f"fit {chosen['selected_residual_over_radius']:.5f}r",
              f"anchors {best_anchor_ratio:.5f}r", f"raw evals {records[oid]['total_forward_evaluations']}",
              flush=True)
    if len(records) != 18:
        raise ArithmeticError("Expected all eighteen originally frozen weak-background observations")
    save_json(OUT / "sealed_selection.json", dict(
        eligible_hypotheses="Five centres selected by original old-only nominal full-wave score",
        homogeneous_background_anchors=list(ANCHORS),
        bounded_uniform_core=[LOW, HIGH], bounded_uniform_exterior=[0., B],
        geometry_load=setup, count=len(records), records=records, timing=elapsed(began),
    ))


def replay():
    started = clock()
    records = json.loads((OUT / "sealed_selection.json").read_text())["records"]
    receipts = {}
    setup = {}
    replay_geometries = {}
    for side in (32, 96):
        one = clock()
        replay_geometries[side] = SectorGeometry(Geometry(geom_path(side)), B)
        setup[str(side)] = elapsed(one)
    for oid, row in records.items():
        one = clock()
        sector = replay_geometries[row["side"]]
        g = sector.g
        with np.load(OUT / row["selected_material"], allow_pickle=False) as frozen:
            x = frozen["material"]
            cells, profile = frozen["active_cells"], frozen["active_profile"]
            core, background = float(frozen["amplitude"]), float(frozen["background"])
            frozen_prediction = frozen["prediction"]
        with np.load(HERE / "public" / f"{oid}.npz", allow_pickle=False) as obs:
            measured, radius = obs["Y"], float(obs["r"])
        if not (LOW <= core <= HIGH and 0 <= background <= B):
            raise ArithmeticError(f"Invalid fitted point parameters for {oid}")
        if np.any(x[cells] != core * profile) or np.any(np.delete(x, cells) != background):
            raise ArithmeticError(f"Fitted point leaves the declared homogeneous template for {oid}")
        physical, current, code = g.forward(x)
        if code:
            raise ArithmeticError(f"Forward replay failed for {oid}")
        raw = float(np.linalg.norm(physical-measured))
        drift = float(np.linalg.norm(physical-frozen_prediction))
        state = float(np.linalg.norm(current-x[:,None]*(g.B+g.apply(current))))
        if drift > 1e-10 or state > 1e-8:
            raise ArithmeticError(f"Healthy original-operator replay failed for {oid}")
        accepted = None
        if raw < radius:
            accepted = certify_witness(sector, x, cells, current, measured, radius, background)
        receipts[oid] = dict(
            raw_residual_over_radius=raw/radius, state_equation_residual_frobenius=state,
            port_delta_to_frozen_prediction=drift, original_operator_witness=accepted,
            raw_inside_noise_ball=raw<radius, **elapsed(one),
        )
        print("replay", oid, f"{raw/radius:.5f}r",
              None if accepted is None else accepted["residual_ratio_upper"], flush=True)
    save_json(OUT / "replay.json", dict(records=receipts, setup=setup, timing=elapsed(started)))


def score():
    fitted = json.loads((OUT / "sealed_selection.json").read_text())
    replayed = json.loads((OUT / "replay.json").read_text())
    rows = []
    for oid, candidate in fitted["records"].items():
        row = replayed["records"][oid]
        with np.load(HERE / "offline" / f"{oid}.npz", allow_pickle=False) as label:
            truth = label["centre"]
        chosen = np.asarray(candidate["selected"]["centre"])
        offset = chosen-truth
        verified = row["original_operator_witness"]
        rows.append(dict(
            observation_id=oid, side=candidate["side"],
            template=candidate["provided_template"], true_xy=truth.tolist(),
            selected_xy=chosen.tolist(), exact_centre=bool(np.array_equal(chosen, truth)),
            position_error_cells=float(np.linalg.norm(offset)),
            generating_centre_in_nominal_top5=bool(
                truth[1]*candidate["side"] + truth[0] in candidate["nominal_top5"]
            ),
            nominal_raw_residual_over_radius=candidate["nominal_best_residual_over_radius"],
            anchors_only_exact_centre=(
                candidate["anchors_only_selected"]["centre"] == truth.tolist()
            ),
            anchors_only_residual_over_radius=candidate["anchors_only_selected"]["residual_over_radius"],
            fitted_core=candidate["selected"]["amplitude"],
            fitted_uniform_bg=candidate["selected"]["background"],
            fitted_raw_residual_over_radius=row["raw_residual_over_radius"],
            verified_witness=bool(verified is not None),
            certified_residual_upper_over_radius=None if verified is None else verified["residual_ratio_upper"],
            forward_evaluations=candidate["total_forward_evaluations"],
            fitting_cpu_s=candidate["cpu_s"], fitting_wall_s=candidate["wall_s"],
            replay_cpu_s=row["cpu_s"], replay_wall_s=row["wall_s"],
            state_equation_residual_frobenius=row["state_equation_residual_frobenius"],
        ))
    summary = {}
    for group in (32, 96, "all"):
        selected = [r for r in rows if group=="all" or r["side"] == group]
        residuals = [r["fitted_raw_residual_over_radius"] for r in selected]
        anchor_residuals = [r["anchors_only_residual_over_radius"] for r in selected]
        strict = [r["certified_residual_upper_over_radius"] for r in selected if r["verified_witness"]]
        summary[str(group)] = dict(
            objects=len(selected), truth_in_nominal_top5=sum(r["generating_centre_in_nominal_top5"] for r in selected),
            exact_centres=sum(r["exact_centre"] for r in selected),
            anchors_only_exact_centres=sum(r["anchors_only_exact_centre"] for r in selected),
            anchors_only_raw_inside=sum(r<1 for r in anchor_residuals),
            lower_than_original_nominal=sum(r["fitted_raw_residual_over_radius"] <
                                            r["nominal_raw_residual_over_radius"] for r in selected),
            fitted_raw_inside=sum(r<1 for r in residuals), verified_witnesses=len(strict),
            fitted_ratio_min_median_max=[min(residuals), float(np.median(residuals)), max(residuals)],
            verified_ratio_range=[min(strict), max(strict)] if strict else None,
            anchors_only_ratio_range=[min(anchor_residuals), max(anchor_residuals)],
            maximum_state_eq_residual=max(r["state_equation_residual_frobenius"] for r in selected),
            evaluated_full_wave_points=sum(r["forward_evaluations"] for r in selected),
            fitting_cpu_s=sum(r["fitting_cpu_s"] for r in selected),
            fitting_wall_s=sum(r["fitting_wall_s"] for r in selected),
            replay_cpu_s=sum(r["replay_cpu_s"] for r in selected),
            replay_wall_s=sum(r["replay_wall_s"] for r in selected),
        )
    save_json(OUT / "summary.json", dict(
        scope="Point-ranking comparator on a supplied outline, same old observations and radius. "
              "The homogeneous exterior and common core amplitude are two fitted parameters. "
              "Existence of a verified point does not certify any other centre infeasible.",
        summary=summary, rows=rows, full_fit_stage=fitted["timing"],
        full_replay_stage=replayed["timing"], geometry_load=fitted["geometry_load"],
        replay_setup=replayed["setup"],
    ))
    with (OUT / "per_object.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == "__main__":
    selected = argparse.ArgumentParser()
    selected.add_argument("stage", choices=("fit", "replay", "score"))
    stage = selected.parse_args().stage
    with threadpool_limits(limits=1,user_api="blas"):
        dict(fit=fit,replay=replay,score=score)[stage]()
