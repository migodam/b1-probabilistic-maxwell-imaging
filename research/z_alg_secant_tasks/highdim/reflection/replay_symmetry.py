"""Reindex archived full-wave observations through verified lattice symmetries."""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb, arb, ctx
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
TOL = 1e-10
ctx.prec = 128
ctx.threads = 1
MODELS = {
    "mirror_x": np.diag([-1., 1.]),
    "mirror_y": np.diag([1., -1.]),
    "turn_90": np.array([[0., -1.], [1., 0.]]),
    "turn_180": -np.eye(2),
    "turn_270": np.array([[0., 1.], [-1., 0.]]),
}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def points_from_angles(angles):
    angles = np.deg2rad(angles)
    return np.column_stack((np.cos(angles), np.sin(angles)))


def grid(side, width):
    axis = (np.arange(side) + .5) * width / side - width / 2
    y, x = np.meshgrid(axis, axis, indexing="ij")
    return np.column_stack((x.ravel(), y.ravel()))


def matched(points, linear_map):
    mapped = points @ linear_map.T
    distance, nearest = cKDTree(points).query(mapped)
    limit = TOL * max(1., np.max(np.abs(points)))
    trial = np.where(distance <= limit, nearest, -1)
    closed = trial >= 0
    while np.any(closed) and not np.all(closed[trial[closed]]):
        closed[closed] = closed[trial[closed]]
    trial[~closed] = -1
    return trial, float(max(distance[closed], default=0.))


def relative_error(a, b):
    return float(np.linalg.norm(a-b) / max(np.linalg.norm(a), np.finfo(float).tiny))


def angular_choices(tx, rx, source_perm, recv_perm, kind):
    fix_t = np.flatnonzero(source_perm == np.arange(len(tx)))
    fix_r = np.flatnonzero(recv_perm == np.arange(len(rx)))
    if not len(fix_t):
        return {"status": "no_axis_invariant_saved_source"}
    if not len(fix_r):
        return {"status": "no_axis_receiver"}
    for source in fix_t:
        direction = tx[source] / np.linalg.norm(tx[source])
        if kind == "line_source":
            direction = -direction
        for axis in sorted(fix_r, key=lambda j: (-np.dot(rx[j]/np.linalg.norm(rx[j]), direction), j)):
            baseline = rx[axis] / np.linalg.norm(rx[axis])
            det = baseline[0]*rx[:, 1]-baseline[1]*rx[:, 0]
            angle = np.rad2deg(np.arctan2(det, rx @ baseline))
            positive = [j for j in np.argsort(angle, kind="stable") if TOL < angle[j] < 180.-TOL]
            mirrored = [j for j in positive if recv_perm[j] >= 0 and recv_perm[j] != j]
            if not mirrored:
                continue
            closest = positive[0] if positive and positive[0] in mirrored else None
            valid_far = [j for j in mirrored if angle[j] <= 30. + TOL]
            far = valid_far[-1] if valid_far else None
            def pair(j):
                return None if j is None else {
                    "receiver": int(j), "reflected_receiver": int(recv_perm[j]), "offset_degrees": float(angle[j]),
                }
            return {
                "status": "axial_source_and_mirrored_receivers", "source": int(source),
                "axis": {"receiver": int(axis), "reflected_receiver": int(axis), "offset_degrees": 0.},
                "nearest_slot": pair(closest), "outermost_pair_le_30deg": pair(far),
                "same_physical_port_for_nearest_and_farthest": bool(closest is not None and closest == far),
            }
    return {"status": "axis_receiver_has_no_paired_off_axis_channel"}


def geometry_protocol(coords, fine, width, tx, rx, E, S, kind, include_turns):
    center = .5*(coords.min(axis=0)+coords.max(axis=0))
    local = coords-center
    tx_local = tx if kind == "plane_wave" else tx-center
    rx_local = rx-center
    fine_coords = grid(fine, width)
    result = {}
    for name, U in MODELS.items():
        if name.startswith("turn") and not include_turns:
            continue
        cell, cell_error = matched(local, U)
        fine_cell, fine_error = matched(fine_coords, U)
        source, tx_error = matched(tx_local, U)
        receive, rx_error = matched(rx_local, U)
        src = np.flatnonzero(source >= 0)
        rec = np.flatnonzero(receive >= 0)
        details = {
            "input_grid_matches": int(np.sum(cell >= 0)), "original_grid_matches": int(np.sum(fine_cell >= 0)),
            "matched_sources": int(len(src)), "matched_receivers": int(len(rec)),
            "matched_ports": int(len(src)*len(rec)), "all_ports": len(rx)*len(tx),
            "max_coordinate_errors": [cell_error, fine_error, tx_error, rx_error],
            "source_order": src.tolist(), "receiver_order": rec.tolist(),
            "source_image_order": source[src].tolist(), "receiver_image_order": receive[rec].tolist(),
            "source_field_error": None, "receiver_operator_error": None,
        }
        if not len(src) or not len(rec):
            details["status"] = "no_closed_source_receiver_product"
        elif np.any(cell < 0) or np.any(fine_cell < 0):
            details["status"] = "domain_grid_not_closed"
        else:
            de = relative_error(E[:, src], E[np.ix_(cell, source[src])])
            ds = relative_error(S[rec], S[np.ix_(receive[rec], cell)])
            details["source_field_error"], details["receiver_operator_error"] = de, ds
            if max(de, ds) > TOL:
                details["status"] = "incident_or_receiver_phase_law_not_verified"
            elif np.all(source[src] == src) and np.all(receive[rec] == rec):
                details["status"] = "only_fixed_channels_no_nontrivial_action"
            else:
                details["status"] = "verified_closed_permutation"
        if name.startswith("mirror"):
            details["single_source"] = angular_choices(tx_local, rx_local, source, receive, kind)
        result[name] = details
    allowed = [name for name in ("mirror_x", "mirror_y") if result[name]["status"] == "verified_closed_permutation"]
    return {
        "grid_center": center.tolist(), "tx_kind": kind,
        "source_coordinates_or_directions": tx.tolist(), "receivers": rx.tolist(), "symmetries": result,
        "largest_verified_mirror": max(allowed, key=lambda name: result[name]["matched_ports"]) if allowed else None,
        "primary_single_source_mirror": next((name for name in ("mirror_x", "mirror_y")
            if result[name]["status"] == "verified_closed_permutation"
            and result[name]["single_source"]["status"] == "axial_source_and_mirrored_receivers"), None),
    }


def plan(out):
    from archives import AUS, AUS_RESULTS, R8, R10, SRAP, acquisition, rel

    start, cpu = perf_counter(), process_time()
    geometries, views, excluded = {}, [], []
    inputs = {}

    def remember_geometry(geom, centers, E, S):
        index = {}
        for label, value in (("centers", centers), ("E", E), ("S", S)):
            match = next((k for k, array in inputs.items()
                          if k.startswith(label+"__") and np.array_equal(array, value)), None)
            if match is None:
                match = f"{label}__{len(inputs)}"
                inputs[match] = value
            index[label] = match
        geom["geometry_audit_arrays"] = index

    ring = R10/"data_public/twod"
    with np.load(ring/"operators.npz", allow_pickle=False) as a:
        g = geometry_protocol(a["centers"], 24, .96, a["sources"], a["receivers"], a["E0"], a["S"],
                              "line_source", include_turns=True)
        remember_geometry(g, a["centers"], a["E0"], a["S"])
        g.update({"frequency_hz": None, "wavenumber": float(a["k"]),
                  "operator_reference": rel(ring/"operators.npz")})
    geometries["r10_same_ring"] = g
    for file in sorted(ring.glob("r10_twod_*.npz")):
        views.append({"id": file.stem, "object": file.stem, "cohort": "R10", "layout": "r10_same_ring",
                      "raw_noisy": rel(file), "score_input": rel(R10/"sealed_truth/twod"/file.name),
                      "clean_saved": True, "original_side": 24, "inverse_side": 24, "domain_width": .96})

    aus_file = AUS_RESULTS/"PUBLIC_DATA.npz"
    # The angle rule comes from the physical R8 source; no Y or labels are loaded here.
    with np.load(aus_file, allow_pickle=False) as a:
        original, inverse, freq = (int(a[k]) for k in ("forward_side", "inverse_side", "frequency"))
    a = acquisition(inverse, freq, dense=False)
    g = geometry_protocol(a["centers"], original, 2., points_from_angles(18*np.arange(20)),
                          3.75*points_from_angles(9*np.arange(40)), a["E0"], a["S"],
                          "plane_wave", include_turns=True)
    remember_geometry(g, a["centers"], a["E0"], a["S"])
    g.update({"frequency_hz": freq, "wavenumber": float(a["k"]),
              "operator_reference": rel(R8/"src/paper_physics.py")+"::acquisition"})
    geometries["austria_20x40"] = g
    views.append({"id": AUS, "object": "austria_eps4", "cohort": "Austria", "layout": "austria_20x40",
                  "raw_noisy": rel(aus_file), "score_input": rel(R8/"data/paper_benchmarks_private"/f"{AUS}.npz"),
                  "clean_saved": False, "original_side": original, "inverse_side": inverse, "domain_width": 2.})

    for folder in ("data", "data_followup"):
        group = SRAP/folder
        for source in read_json(group/"cases.json"):
            if source.get("suite") != "LA_SYNTH" or source.get("fine") != 96 or source["side"] != 32:
                excluded.append({"id": source["case"], "reason": "outside_frozen_96_to_32_synthetic_contract"})
                continue
            key = "srap_"+source["array"]
            path = group/"public"/f"{source['case']}.npz"
            # Layout key contains the geometry; verify every source manifest agrees.
            source_setup = (source["tx_kind"], source["tx_plane_wave_propagation_deg"],
                            source["tx_line_source_m"], source["rx_m"])
            if key not in geometries:
                tx = (points_from_angles(source_setup[1]) if source_setup[0] == "plane_wave"
                      else np.array(source_setup[2]))
                rx = np.array(source_setup[3])
                with np.load(path, allow_pickle=False) as a:
                    g = geometry_protocol(a["centers"], 96, source["width"], tx, rx, a["E"], a["S"],
                                          source_setup[0], include_turns=False)
                    remember_geometry(g, a["centers"], a["E"], a["S"])
                    g.update({"frequency_hz": source["frequency"], "wavenumber": float(a["k"]),
                              "representative_geometry_archive": rel(path)})
                g["published_geometry_inputs"] = list(source_setup)
                geometries[key] = g
            elif geometries[key]["published_geometry_inputs"] != list(source_setup):
                raise ValueError("Layout name reused for a different physical geometry")
            views.append({"id": source["case"], "object": source["object_id"],
                          "cohort": "SRAP/"+source["cohort"], "layout": key, "raw_noisy": rel(path),
                          "score_input": rel(group/"truth"/f"{source['case']}.npz"), "clean_saved": False,
                          "original_side": 96, "inverse_side": 32, "domain_width": source["width"]})

    np.savez_compressed(out/"public_geometry.npz", **inputs)
    write(out/"SELECTION.json", {
        "scope": "geometry and saved source/receiver operators only; no response or material array read",
        "axis_convention": "mirror_x negates x; mirror_y negates y; both through the exact grid center",
        "selection": "x then y on coverage ties; first fixed source and one fixed receiver axis favoring forward direction; positive angular neighbor and last reflected pair within 30 degrees",
        "solver": "archived observations only, no new forward solve",
        "separate_instrument_hypothesis": "Each archived inspected complex sample and any future single output may have deterministic error at most 10% of the observed source row RMS. This condition is additional; the archived noise law is joint.",
        "view_manifest": views, "out_of_scope_views": excluded, "geometries": geometries,
        "cost": {"wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
                 "material_arrays_read": 0, "observation_arrays_read": 0, "new_full_wave_rhs": 0},
    })
    print("Planned", len(views), "views from", len(geometries), "geometries", flush=True)


def austria_exact_noise_radius(y):
    from archives import AUS, R8

    # The frozen R8 generator e/||e|| is used only to recover the scalar budget,
    # by solving ||e|| = .1 ||Y-e||. No clean complex values are reconstructed.
    table = read_json(R8/"delegated/paper_baselines/PAPER_CONTRACTS.json")["NIE2016"]["cases"]
    i = next(j for j, row in enumerate(table) if row["id"] == AUS)
    rng = np.random.default_rng(8600+i)
    u = rng.normal(size=y.shape)+1j*rng.normal(size=y.shape)
    u /= np.linalg.norm(u)
    f = .1
    b = f*f*np.vdot(u, y).real
    c = -f*f*np.linalg.norm(y)**2
    return float((-b+np.sqrt(b*b-(1-f*f)*c))/(1-f*f))


def norm(x):
    return float(np.linalg.norm(x))


def channels(data, meta, partner=False):
    what = "_image_order" if partner else "_order"
    return data[np.ix_(meta["receiver"+what], meta["source"+what])]


def lift_complex(values):
    return [acb(float(z.real), float(z.imag)) for z in np.ravel(values)]


def squared_length(numbers):
    return sum(z.real*z.real+z.imag*z.imag for z in numbers)


def rounded_ball(x):
    return [float(np.nextafter(float(x.lower()), -np.inf)),
            float(np.nextafter(float(x.upper()), np.inf))]


def strict_orbit(cy, qy, delta, involution):
    a, b = lift_complex(cy), lift_complex(qy)
    diffs = [(x-z)/2 for x, z in zip(a, b, strict=True)]
    odd = squared_length(diffs).sqrt()
    floor = arb(float(delta))
    if involution:
        evens = [(x+z)/2 for x, z in zip(a, b, strict=True)]
        even = squared_length(evens).sqrt()
        gain = even if even < odd else odd if odd < even else (
            min(even.lower(), odd.lower()).union(min(even.upper(), odd.upper())))
        extra = {"strict_E_interval": rounded_ball(even), "strict_O_interval": rounded_ball(odd)}
    else:
        inner = sum(x.conjugate()*z for x, z in zip(a, b, strict=True))
        defect = (squared_length(a)-abs(inner))/2
        gain = defect.nonnegative_part().sqrt()
        extra = {}
    return {
        **extra,
        "strict_calibrated_interval": rounded_ball(odd),
        "strict_gain_free_interval": rounded_ball(gain),
        "strict_gain_free_margin_over_delta": rounded_ball(gain/floor-1),
        "strict_calibrated_margin_over_delta": rounded_ball(odd/floor-1),
        "calibrated_true_under_conservative_radius": bool(odd > floor),
        "gain_free_true_under_conservative_radius": bool(gain > floor),
    }


def metrics(y, delta, change, name, clean, delta_upper):
    cy, qy = channels(y, change), channels(y, change, partner=True)
    even, odd = norm((cy+qy)/2), norm((cy-qy)/2)
    involution = name.startswith("mirror") or name == "turn_180"
    proof = strict_orbit(cy, qy, delta_upper, involution)
    values = {**proof}
    if involution:
        values.update({"E": even, "O": odd, "E_over_delta": even/delta, "O_over_delta": odd/delta,
                       "min_EO_over_delta": min(even, odd)/delta,
                       "no_orbit_pair_at_fixed_calibration": proof["calibrated_true_under_conservative_radius"],
                       "no_orbit_pair_even_with_free_complex_gain": proof["gain_free_true_under_conservative_radius"]})
    # Stable Gram-SVD form; equivalent to sqrt((||y||^2-|<y,Qy>|)/2).
    a = np.column_stack((cy.ravel(), qy.ravel()))
    d = 0. if a.shape[0] < 2 else float(np.linalg.svd(a, compute_uv=False)[-1]/np.sqrt(2))
    distance = norm(cy-qy)/2
    values.update({"free_gain_orbit_distance": d, "free_gain_orbit_distance_over_delta": d/delta,
                   "free_gain_orbit_excluded": proof["gain_free_true_under_conservative_radius"],
                   "calibrated_orbit_distance": distance,
                   "calibrated_orbit_excluded": proof["calibrated_true_under_conservative_radius"]})
    if clean is not None:
        original, image = channels(clean, change), channels(clean, change, partner=True)
        values["clean_change_distance"] = norm(original-image)
    return values


def gap_instrument(y, clean, delta, plan, meta):
    tx = plan["source"]
    ds, nr = y[:, tx], len(y)
    ready = np.array([j for j in meta["receiver_order"]])
    reflected = np.array(meta["receiver_image_order"])
    odd = norm((ds[ready]-ds[reflected])/2)
    output = {
        "source": tx, "matched_receiver_count": len(ready), "entire_receiver_row_closed": len(ready)==nr,
        "noisy_odd_norm": odd, "clean_odd_bound_from_joint_radius": [max(0., odd-delta), odd+delta],
        "row_rmS_signal": norm(ds)/np.sqrt(nr), "independent_hypothetical_per_sample_radius": .1*norm(ds)/np.sqrt(nr),
    }
    if clean is not None:
        output["clean_odd_norm"] = norm((clean[ready, tx]-clean[reflected, tx])/2)
    for arm in ("axis", "nearest_slot", "outermost_pair_le_30deg"):
        pair = plan[arm]
        if pair is None:
            output[arm] = None
            continue
        i, j = pair["receiver"], pair["reflected_receiver"]
        observed = abs(ds[i]-ds[j])
        discrepancy = np.sqrt(2)*delta if i != j else 0.
        lower, upper = max(0., observed-discrepancy), observed+discrepancy
        rho = output["independent_hypothetical_per_sample_radius"]
        if clean is not None:
            lower = upper = abs(clean[i, tx] - clean[j, tx])
        output[arm] = {
            **pair, "observed_complex_distance": float(observed),
            "true_complex_distance_or_joint_noise_bound": [float(lower), float(upper)],
            "distance_type": "saved_clean_field" if clean is not None else "noisy_triangle_bound",
            "original_joint_delta_as_absolute_single_channel": delta,
            "pair_excluded_using_archived_joint_budget": bool(lower > 2*delta),
            "pair_has_common_observation_using_archived_joint_budget": bool(upper <= 2*delta),
            "hypothetical_per_sample_model_excludes_pair": bool(
                lower > 2*rho if clean is not None else observed > 4*rho),
        }
    return output


def replay(out):
    from archives import WORKSPACE

    plan = read_json(out/"SELECTION.json")
    start, cpu = perf_counter(), process_time()
    rows, arrays = [], {}
    for view in plan["view_manifest"]:
        case, kind = view["id"], view["cohort"]
        geo = plan["geometries"][view["layout"]]
        with np.load(WORKSPACE/view["raw_noisy"], allow_pickle=False) as arc:
            y = arc["Y"]
            if kind == "R10":
                delta = float(arc["sigma"]*arc["noise_radius"])
                exact_radius = arb(float(arc["sigma"]))*arb(float(arc["noise_radius"]))
                delta_definition = "original CN99 whole matrix: saved sigma times chi-square radius"
            elif kind == "Austria":
                delta = austria_exact_noise_radius(y)
                # The original 0.1 relative noise law and triangle inequality give
                # ||noise|| <= ||Y||/9; the tiny slack covers binary64 generation rounding.
                exact_radius = squared_length(lift_complex(y)).sqrt()/(arb(9)-arb("1e-11"))
                delta_definition = "original R8 fixed 10% of clean matrix norm; recover scalar only from archived generator and Y"
            else:
                delta = float(arc["delta"])
                exact_radius = arb(delta)
                delta_definition = "original SRAP fixed 10% complex-noise joint matrix norm"
        delta_upper = rounded_ball(exact_radius)[1]
        y_clean = None
        if view["clean_saved"]:
            with np.load(WORKSPACE/view["score_input"], allow_pickle=False) as archive:
                y_clean = archive["Y_clean"]
        rowsy, colsy = (len(geo["receivers"]), len(geo["source_coordinates_or_directions"]))
        if y.shape != (rowsy, colsy):
            raise ValueError(f"Raw channel map differs from selected geometry for {case}")
        noise = None if y_clean is None else norm(y-y_clean)
        arrays[f"{case}__Y_noisy"] = y
        arrays[f"{case}__joint_delta"] = np.array(delta)
        arrays[f"{case}__proof_delta_upper"] = np.array(delta_upper)
        row = {
            **view, "archive_radius": delta, "archive_radius_definition": delta_definition,
            "strict_certificate_threshold": delta_upper,
            "noisy_matrix_frobenius": norm(y), "largest_verified_mirror": geo["largest_verified_mirror"],
            "primary_single_source_mirror": geo["primary_single_source_mirror"],
            "archive_realized_noise_norm": noise,
            "archive_realized_noise_inside_budget": None if noise is None else noise <= delta,
            "results": {},
        }
        if kind == "Austria":
            row["same_saved_matrix_historical_r12_radius"] = .12*norm(y)
        for name, design in geo["symmetries"].items():
            result = {"geometry_status": design["status"], "retained_ports": design["matched_ports"],
                      "all_ports": design["all_ports"]}
            if design["status"] in ("verified_closed_permutation", "only_fixed_channels_no_nontrivial_action"):
                result["all_data"] = metrics(y, delta, design, name, y_clean, delta_upper)
                if name.startswith("mirror"):
                    edge = design["single_source"]
                    if edge["status"] == "axial_source_and_mirrored_receivers":
                        result["one_source"] = gap_instrument(y,y_clean,delta,edge,design)
            row["results"][name] = result
        if y_clean is not None:
            arrays[f"{case}__saved_clean"] = y_clean
        rows.append(row)
    np.savez_compressed(out/"complex_readouts.npz", **arrays)
    write(out/"field_results.json", {
        "rows": rows, "material_GT_read": False,
        "radius_recovery": "cross-channel readout comparisons use ∣D_y−D_f∣≤√2 δ; virtual single output inherits δ, not δ/√count",
        "turn_conditions": "source and receiver both permuted; independent gain model uses smallest singular value of [Y,QY]/sqrt(2)",
        "cost": {"wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
                 "raw_noisy_matrices": len(rows), "stored_clean_matrices": sum(r["clean_saved"] for r in rows),
                 "new_full_wave_rhs": 0, "material_arrays_read": 0},
    })
    print("Replayed", len(rows), "saved views", flush=True)


def score(out):
    from archives import WORKSPACE

    start, cpu = perf_counter(), process_time()
    selected = read_json(out/"SELECTION.json")
    frozen = read_json(out/"field_results.json")
    if [x["id"] for x in selected["view_manifest"]] != [x["id"] for x in frozen["rows"]]:
        raise ValueError("The observation replay and exposure-frozen plan do not match")
    scores, aliases = {}, {}
    for row in frozen["rows"]:
        key = row["cohort"]+"/"+row["object"]
        aliases[row["id"]] = key
        if key in scores:
            continue
        side, width = row["inverse_side"], row["domain_width"]
        with np.load(WORKSPACE/row["score_input"], allow_pickle=False) as data:
            chi = np.asarray(data["m"], dtype=float).ravel()
        coords, result = grid(side, width), {}
        if len(chi) != len(coords):
            raise ValueError("Archived material score dimensions changed.")
        for i, axis in enumerate(("mirror_x", "mirror_y")):
            x = coords[:, i]
            result[axis] = {"negative_half_mean": float(np.mean(chi[x<0])),
                            "positive_half_mean": float(np.mean(chi[x>0])),
                            "half_mean_material_task_gap": float(abs(np.mean(chi[x<0])-np.mean(chi[x>0]))),
                            "domain_mean_first_moment_gap": float(abs(2*np.mean(x*chi)))}
        scores[key] = result
    write(out/"offline_task_scores.json", {
        "scope": "coarse archived material arrays opened only after geometry selection and complete data replay",
        "aliases": aliases, "scores": scores,
        "cost": {"wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
                 "distinct_materials_scored": len(scores), "new_full_wave_rhs": 0},
    })
    print("Scored", len(scores), "distinct materials offline", flush=True)


def frozen(out):
    receipt = read_json(out/"SELECTION.json")
    start, cpu = perf_counter(), process_time()
    layouts = receipt["geometries"]
    samples = {v["layout"]: v for v in receipt["view_manifest"]}
    rechecked = {}
    with np.load(out/"public_geometry.npz", allow_pickle=False) as geometry:
        for name, old in layouts.items():
            inputs = old["geometry_audit_arrays"]
            coordinates, E, S = [geometry[inputs[k]] for k in ("centers", "E", "S")]
            sample = samples[name]
            new = geometry_protocol(coordinates, sample["original_side"], sample["domain_width"],
                np.asarray(old["source_coordinates_or_directions"]), np.asarray(old["receivers"]), E, S,
                old["tx_kind"], any(key.startswith("turn") for key in old["symmetries"]))
            for key in ("largest_verified_mirror", "primary_single_source_mirror"):
                if new[key] != old[key]:
                    raise ValueError(f"Frozen geometry and original plan disagree for {name}/{key}")
            for key, definition in new["symmetries"].items():
                published = old["symmetries"][key]
                plan_fields = ("status", "source_order", "receiver_order",
                               "source_image_order", "receiver_image_order")
                if any(definition[field] != published[field] for field in plan_fields):
                    raise ValueError(f"Frozen port ordering and original plan disagree for {name}/{key}")
                if "single_source" in definition:
                    then, now = (item["single_source"] for item in (published, definition))
                    signature = lambda x: (x["status"], x.get("source"), *(
                        None if not x.get(group) else (x[group]["receiver"], x[group]["reflected_receiver"])
                        for group in ("axis", "nearest_slot", "outermost_pair_le_30deg")))
                    if signature(then) != signature(now):
                        raise ValueError(f"Frozen single-source geometry disagree for {name}/{key}")
            rechecked[name] = new
    rows = []
    with np.load(out/"complex_readouts.npz", allow_pickle=False) as archive:
        for selection in receipt["view_manifest"]:
            case = selection["id"]
            geo = rechecked[selection["layout"]]
            obs = archive[f"{case}__Y_noisy"]
            clean_key = f"{case}__saved_clean"
            clean = archive[clean_key] if clean_key in archive else None
            delta = float(archive[f"{case}__joint_delta"])
            delta_upper = float(archive[f"{case}__proof_delta_upper"])
            results = {}
            for name, sym in geo["symmetries"].items():
                if sym["status"] in ("verified_closed_permutation", "only_fixed_channels_no_nontrivial_action"):
                    domain = metrics(obs, delta, sym, name, clean, delta_upper)
                    if "single_source" in sym and sym["single_source"]["status"] == "axial_source_and_mirrored_receivers":
                        domain["one_source"] = gap_instrument(obs, clean, delta, sym["single_source"], sym)
                    results[name] = domain
            rows.append({"id": case, "primary": geo["largest_verified_mirror"], "results": results})
    write(out/"PUBLIC_RECHECK.json", {
        "rechecked_layouts": len(rechecked), "rechecked_views": len(rows), "new_full_wave_rhs": 0,
        "material_arrays_read": 0, "rows": rows,
        "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
    })
    print("Standalone public replay verified", len(rechecked), "layouts and", len(rows), "views", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=HERE/"results")
    p.add_argument("--stage", choices=("plan", "replay", "score", "frozen", "all"), default="all")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    stages = (plan, replay, score) if a.stage == "all" else (
        {"plan": plan, "replay": replay, "score": score, "frozen": frozen}[a.stage],)
    for stage in stages:
        stage(a.output)


if __name__ == "__main__":
    main()
