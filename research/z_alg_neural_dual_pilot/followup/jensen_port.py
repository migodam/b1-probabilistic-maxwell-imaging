"""Outward port-compressed Jensen bounds for the archived positive .75 anchors."""

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

sys.dont_write_bytecode = True

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh, solve
from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
WORKSPACE = PILOT.parent.parent
DEFAULT_CACHE = HERE / "results/port_enclosures.json"
CAP = .75
PORTS = 20
SOURCES = 10
SIZE = PORTS + SOURCES
GRID = tuple(i / 8 for i in range(9))
PRECISION = 112
THREADS = 2

sys.path.insert(0, str(PILOT / "src"))
import fixed_baselines as archived

original = archived._reference
SavedProfile = archived.SavedProfile
ball, star, mid = original.ball, original.star, original.mid
gershgorin, high, low = original.gershgorin, original.high, original.low
frobenius_squared = original.frobenius_squared


def started():
    return perf_counter(), process_time()


def elapsed(begin):
    return {
        "wall_seconds": perf_counter() - begin[0],
        "cpu_seconds": process_time() - begin[1],
    }


def atomic_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def relative_receipt(path):
    resolved = Path(path).resolve()
    return str(resolved.relative_to(WORKSPACE)) if resolved.is_relative_to(WORKSPACE) else resolved.name


def hermitian(a):
    return (a + a.conj().T) / 2


def block(matrix, i, j, rows, cols):
    return acb_mat(rows, cols, [matrix[r, c] for r in range(i, i + rows) for c in range(j, j + cols)])


def trace_real(matrix):
    return sum((matrix[i, i].real for i in range(matrix.nrows())), arb(0))


def point_pair(value):
    mantissa, exponent = value.man_exp()
    return str(mantissa), int(exponent)


def encode_real(value):
    return point_pair(value.mid()), point_pair(value.rad())


def decode_real(packed):
    midpoint, radius = ((int(mantissa), exponent) for mantissa, exponent in packed)
    return arb(midpoint, radius)


def encode_complex_matrix(matrix):
    return [(encode_real(value.real), encode_real(value.imag)) for value in matrix.entries()]


def decode_complex_matrix(packed):
    return acb_mat(SIZE, SIZE, [acb(decode_real(real), decode_real(imag)) for real, imag in packed])


@dataclass
class PortModel:
    label: str
    response: acb_mat
    receipt: dict


@dataclass
class Setup:
    anchors: tuple
    midpoint: PortModel
    receipt: dict


class MidpointProfile(SavedProfile):
    """Binary64 proposes the solves; normalized midpoint products stay in Arb."""

    def __init__(self, geometry, left, right):
        if left["mu"] != right["mu"]:
            raise ValueError("These archived anchors must share the exact physical scale.")
        real = tuple(zip(np.asarray(left["sector"], dtype=float).flat,
                         np.asarray(right["sector"], dtype=float).flat, strict=True))
        opt = tuple(zip(np.asarray(left["optical"], dtype=float).flat,
                        np.asarray(right["optical"], dtype=float).flat, strict=True))
        self.exact_real = [(arb(float(a)) + arb(float(b))) / 2 for a, b in real]
        self.exact_optical = [(arb(float(a)) + arb(float(b))) / 2 for a, b in opt]
        spec = {
            "mu": left["mu"],
            "sector": [float(v) for v in self.exact_real],
            "optical": [float(v) for v in self.exact_optical],
        }
        super().__init__(geometry, spec, CAP)

    def exact_scaled_coefficients(self):
        z = [acb(arb(CAP) * a, b) / 2 for a, b in zip(self.exact_real, self.exact_optical, strict=True)]
        return [self.exact_real[i] for i in self.g.blocks], [z[i] for i in self.g.blocks]


def enclose_response(action, normalized_gap, label, profile_receipt):
    begin = started()
    costs = action.g.cost.copy()
    m, defect, work = original.exact_common_response(action, nu=0)
    # A = K / alpha and l = T / alpha. For R = AX - l,
    # T^* K^{-1} T = -alpha m + alpha R^* A^{-1} R.
    scale = arb(float(action.w.scale))
    correction = scale * defect / arb(normalized_gap)
    m = -scale * m
    caps = [max(arb(0), correction[i, i].real.upper()) for i in range(SIZE)]
    response = acb_mat(SIZE, SIZE)
    for i, bound in enumerate(caps):
        response[i, i] = acb(m[i, i].real + arb(0).union(bound))
        for j in range(i + 1, SIZE):
            eps = (bound * caps[j]).sqrt().upper()
            noise = (-eps).union(eps)
            entry = (m[i, j] + m[j, i].conjugate()) / 2 + acb(noise, noise)
            response[i, j], response[j, i] = entry, entry.conjugate()
    receipt = {
        "profile": profile_receipt,
        "normalized_state_positive_gap_lower": normalized_gap,
        "native_state_positive_gap_lower": low(scale * arb(normalized_gap)),
        "stationary": work,
        "correction_operator_norm_upper": high(sum(caps, arb(0))),
        "largest_diagonal_correction_upper": high(max(caps)),
        "operations": dict(action.g.cost - costs),
        **elapsed(begin),
    }
    return PortModel(label, response, receipt)


def archived_anchors():
    profiles = {p.index: p for p in archived.fixed_profiles() if p.cap == CAP and p.index in (0, 2)}
    left, right = profiles[0], profiles[2]
    if not all(p.archived_inertia["nu"] == 0 and p.archived_inertia["positive_gap"] > 0.
               for p in (left, right)):
        raise ValueError("The two profiles must have exact archived positive native gaps.")
    return left, right


def load_setup(path, profiles):
    raw = json.loads(path.read_text())
    if raw["format"] != "dyadic-arb-balls-v1" or raw["state_cells"] != 9216 or raw["cap"] != CAP:
        raise ValueError("Archived reduced response has a different format or physical contract.")
    if raw["precision_bits"] != PRECISION or raw["portfolio_order"] != [p.index for p in profiles]:
        raise ValueError("Archived reduced response has different anchor settings.")
    models = tuple(
        PortModel(record["label"], decode_complex_matrix(record["response"]), record["receipt"])
        for record in raw["records"]
    )
    for model, p in zip(models[:2], profiles, strict=True):
        if model.receipt["profile"]["specification"] != p.specification:
            raise ValueError("The persisted anchor differs from the current archived specification.")
    if len(models) != 3 or tuple(v.label for v in models) != ("anchor_0", "anchor_2", "native_midpoint"):
        raise ValueError("Expected two anchor responses and the exact native midpoint response.")
    receipt = raw["setup"]
    receipt.setdefault("reusable_anchor_setup", {
        key: receipt["geometry_setup"][key] + sum(model.receipt[key] for model in models[:2])
        for key in ("wall_seconds", "cpu_seconds")
    })
    receipt.setdefault("additional_midpoint_control_setup", {
        key: models[2].receipt[key] for key in ("wall_seconds", "cpu_seconds")
    })
    return Setup(models[:2], models[2], receipt)


def prepare(cache=DEFAULT_CACHE, *, rebuild=False):
    """Enclose the 30×30 native responses and save a reusable representation."""
    ctx.prec, ctx.threads = PRECISION, THREADS
    cache = Path(cache)
    profiles = archived_anchors()
    if cache.is_file() and not rebuild:
        return load_setup(cache, profiles)
    beginning = started()
    with threadpool_limits(limits=1):
        geometry_time = started()
        geometry = archived.Geometry(*archived.native_geometry(), block_side=4)
        if (geometry.n, len(geometry.sensor), geometry.s) != (9216, PORTS, SOURCES):
            raise ValueError("Unexpected native operator or port order.")
        geometry_work = elapsed(geometry_time)
        models = []
        for p in profiles:
            spec = {
                "index": p.index,
                "specification": p.specification,
                "archived_positive_gap": p.archived_inertia,
                "coefficient_convention": "exact product of represented scale and normalized coefficient",
            }
            model = enclose_response(p.action(geometry), p.archived_inertia["positive_gap"],
                                     f"anchor_{p.index}", spec)
            models.append(model)
            print(model.label, model.receipt["wall_seconds"], flush=True)
        profile = MidpointProfile(geometry, profiles[0].specification, profiles[1].specification)
        gap = min(p.archived_inertia["positive_gap"] for p in profiles)
        spec = {
            "anchor_indices": [p.index for p in profiles],
            "t": .5,
            "gap_method": "positive convex combination with common exact scale",
            "coefficient_convention": "exact Arb half sums of the normalized anchor coefficients",
        }
        model = enclose_response(profile, gap, "native_midpoint", spec)
        print(model.label, model.receipt["wall_seconds"], flush=True)
        models.append(model)
    setup = {
        "geometry_setup": geometry_work,
        "profile_work": {model.label: model.receipt for model in models},
        "reusable_anchor_setup": {
            key: geometry_work[key] + sum(model.receipt[key] for model in models[:2])
            for key in geometry_work
        },
        "additional_midpoint_control_setup": {
            key: models[2].receipt[key] for key in geometry_work
        },
        "reusable_anchor_native_rhs": 60,
        "control_only_native_rhs": 30,
        "physical_forward_rhs": 0,
        "linear_solves_online": "20 by 20 port systems only",
        **elapsed(beginning),
    }
    raw = {
        "format": "dyadic-arb-balls-v1",
        "cap": CAP, "state_cells": geometry.n, "ports": PORTS, "sources": SOURCES,
        "precision_bits": PRECISION, "threads": THREADS,
        "portfolio_order": [p.index for p in profiles],
        "geometry_source": str(archived.FINE.relative_to(PILOT.parent.parent)),
        "setup": setup,
        "records": [{"label": model.label, "receipt": model.receipt,
                     "response": encode_complex_matrix(model.response)} for model in models],
    }
    atomic_json(cache, raw)
    return Setup(tuple(models[:2]), models[2], setup)


def source_bounds(matrix, *, upper):
    approximate, modes = eigh(hermitian(mid(matrix)), check_finite=False, driver="evd")
    approximate, modes = approximate[::-1], modes[:, ::-1]
    positive = int(np.count_nonzero(approximate > 0.))
    if positive:
        directions = ball(modes[:, :positive])
        gram_upper = max(arb(1), gershgorin(star(directions) * directions, upper=True))
        contraction_trace = trace_real(star(directions) * matrix * directions)
        lower = max(arb(0), (contraction_trace / gram_upper).lower())
    else:
        gram_upper, contraction_trace, lower = arb(1), arb(0), arb(0)
    record = {
        "squared_lower": max(0., low(lower)) if lower > 0. else 0.,
        "distance_lower": max(0., low(lower.sqrt())) if lower > 0. else 0.,
        "numeric_positive_trace": float(np.maximum(0., approximate).sum()),
        "trial_projector_rank": positive,
        "trial_gram_norm_upper": high(gram_upper),
        "trial_congruence_trace_lower": low(contraction_trace) if positive else 0.,
    }
    if upper:
        directions = ball(modes)
        gram_lower = gershgorin(star(directions) * directions, upper=False)
        if not gram_lower > 0.:
            raise ArithmeticError("The source-space proposal is not certified invertible.")
        congruence = star(directions) * matrix * directions
        for i, a in enumerate(approximate):
            congruence[i, i] -= arb(float(a))
        bump = frobenius_squared(congruence).sqrt().upper()
        ceilings = [max(arb(0), (arb(float(a)) + bump).upper()) / gram_lower for a in approximate]
        ceiling = sum(ceilings, arb(0))
        record.update(
            squared_upper=high(ceiling) if ceiling > 0. else 0.,
            distance_upper=high(ceiling.sqrt()) if ceiling > 0. else 0.,
            source_congruence_error_upper=high(bump),
        )
    return record


def pd_floor(matrix):
    floor = gershgorin(matrix, upper=False)
    if floor > 0.:
        return floor
    modes = ball(eigh(hermitian(mid(matrix)), check_finite=False, driver="evd")[1])
    floor = gershgorin(star(modes) * matrix * modes, upper=False)
    if not floor > 0.:
        return None
    return (floor / gershgorin(star(modes) * modes, upper=True)).lower()


def scalar_bounds(value):
    return {"lower": low(value), "upper": high(value)}


def source_forms(setup, y):
    yb = ball(y)
    eye = acb_mat(SOURCES, SOURCES)
    for i in range(SOURCES):
        eye[i, i] = 1
    xi = original.vertical(yb, eye)
    data_gram = star(yb) * yb
    q0, q1 = (a.response for a in setup.anchors)
    endpoints = [data_gram - star(xi) * a * xi for a in (q0, q1)]
    e = block(q1 - q0, 0, 0, PORTS, SIZE) * xi
    r0, r1 = (block(q, 0, 0, PORTS, PORTS) for q in (q0, q1))
    return endpoints, e, r0, r1, xi, data_gram


def numerical_port_candidates(hh, e, r0, r1):
    numeric_h = [hermitian(mid(h)) for h in hh]
    numeric_e = mid(e)
    numeric_r0, numeric_r1 = (hermitian(mid(q)) for q in (r0, r1))
    candidates = []
    for t in GRID[1:-1]:
        a = t * numeric_r0 + (1. - t) * numeric_r1
        x = solve(a, numeric_e, assume_a="pos", check_finite=False)
        trace_matrix = (1. - t) * numeric_h[0] + t * numeric_h[1] + t * (1. - t) * numeric_e.conj().T @ x
        eigenvalues = eigh(hermitian(trace_matrix), eigvals_only=True, check_finite=False, driver="evd")
        rank = int(np.count_nonzero(eigenvalues > 0.))
        candidates.append((t, float(np.maximum(eigenvalues, 0.).sum()), rank, x))
    return candidates


def port_form(hh, e, r0, r1, t_value, x, *, trace_upper=False):
    t = arb(t_value)
    a = t * r0 + (1 - t) * r1
    x = ball(x)
    ax = a * x
    e_x = star(e) * x
    affine = (1 - t) * hh[0] + t * hh[1]
    gap = (t * (1 - t)) * (e_x + star(e_x) - star(x) * ax)
    lower_form = affine + gap
    trace = trace_real(gap)
    trace_result = {"lower": max(0., low(trace)), "upper": None}
    if trace_upper:
        port_gap = pd_floor(a)
        defect = frobenius_squared(ax - e)
        if port_gap is not None:
            remainder = t * (1 - t) * defect / port_gap
            trace_result.update(upper=high(trace + remainder), port_gap_lower=low(port_gap),
                                port_solve_defect_upper=high(defect.sqrt()))
    return lower_form, affine, gap, trace_result


def evaluate(setup, y, radius, *, object_id=None, noise=None):
    """Evaluate both rigorous endpoints, the best interior proposal, and the exact midpoint."""
    beginning = started()
    y = np.asarray(y, dtype=np.complex128)
    radius = float(radius)
    if y.shape != (PORTS, SOURCES) or not np.isfinite(y).all() or not np.isfinite(radius) or radius < 0.:
        raise ValueError("Expected one finite 20×10 native observation and nonnegative radius.")
    with threadpool_limits(limits=1):
        form_begin = started()
        hh, e, r0, r1, xi, data_gram = source_forms(setup, y)
        timings = {"reduced_forms": elapsed(form_begin)}
        t0 = started()
        candidates = numerical_port_candidates(hh, e, r0, r1)
        best = max(candidates, key=lambda point: (point[1], -point[0]))
        middle = next(point for point in candidates if point[0] == .5)
        timings["interior_binary64_selection"] = elapsed(t0)
        t0 = started()
        ends = [source_bounds(matrix, upper=True) for matrix in hh]
        timings["strict_endpoints"] = elapsed(t0)
        t0 = started()
        selected_form, _, _, _ = port_form(hh, e, r0, r1, best[0], best[3])
        selected = source_bounds(selected_form, upper=False)
        timings["strict_port_selected"] = elapsed(t0)
        t0 = started()
        half_form, affine, gap_half, port_half_trace = port_form(hh, e, r0, r1, .5, middle[3], trace_upper=True)
        half = source_bounds(half_form, upper=True)
        timings["strict_port_midpoint_diagnostic"] = elapsed(t0)
        t0 = started()
        midpoint = data_gram - star(xi) * setup.midpoint.response * xi
        native_half = source_bounds(midpoint, upper=True)
        affine_half = source_bounds(affine, upper=True)
        exact_gain = scalar_bounds(trace_real(midpoint - affine))
        exact_gain["lower"] = max(0., exact_gain["lower"])
        remainder = scalar_bounds(trace_real(midpoint - affine - gap_half))
        remainder["lower"] = max(0., remainder["lower"])
        timings["strict_native_midpoint_control"] = elapsed(t0)
        t0 = started()
        # Concavity gives H(t) <= 2(1-t)Hmid - (1-2t)H1 on the left
        # half and H(t) <= 2t Hmid - (2t-1)H0 on the right. The convex
        # trace of the positive part of each affine roof is bounded by
        # its two endpoint values. These statements cover all real t.
        left_roof = source_bounds(2 * midpoint - hh[1], upper=True)
        right_roof = source_bounds(2 * midpoint - hh[0], upper=True)
        segment_upper = max(r["distance_upper"] for r in (left_roof, native_half, right_roof))
        timings["strict_continuous_segment_upper"] = elapsed(t0)
        t0 = started()
        grid_records = []
        for point in candidates:
            weight, score, rank, x = point
            if weight == best[0]:
                proof = selected
            elif weight == .5:
                proof = half
            else:
                grid_form, _, _, _ = port_form(hh, e, r0, r1, weight, x)
                proof = source_bounds(grid_form, upper=False)
            grid_records.append(
                {"t": weight, "numeric_squared_distance": score, "proposed_rank": rank,
                 "distance_lower": proof["distance_lower"], "squared_lower": proof["squared_lower"]}
            )
        timings["strict_grid_diagnostic"] = elapsed(t0)
    if half["distance_lower"] > native_half["distance_upper"]:
        raise ArithmeticError("The projected midpoint lower exceeds the exact native midpoint upper.")
    if port_half_trace["lower"] > exact_gain["upper"]:
        raise ArithmeticError("The projected Jensen trace lower exceeds the exact Jensen trace upper.")
    if max([native_half["distance_lower"]] + [r["distance_lower"] for r in ends + grid_records]) > segment_upper:
        raise ArithmeticError("A verified curve point exceeds the upper for the whole segment.")
    if exact_gain["lower"] > 0. and port_half_trace["upper"] is not None:
        fraction = {
            "lower": max(0., low(arb(port_half_trace["lower"]) / arb(exact_gain["upper"]))),
            "upper": min(1., high(arb(port_half_trace["upper"]) / arb(exact_gain["lower"]))),
        }
    else:
        fraction = None
    best_end = max(range(2), key=lambda index: ends[index]["distance_lower"])
    chosen_lower = selected["distance_lower"]
    endpoint_lower = ends[best_end]["distance_lower"]
    endpoint_upper = max(a["distance_upper"] for a in ends)
    best_portfolio = max(endpoint_lower, chosen_lower)
    selected.update(t=best[0], passes=chosen_lower > radius)
    for t, record in zip((0., 1.), ends, strict=True):
        record.update(t=t, passes=record["distance_lower"] > radius)
    for record in (native_half, half, affine_half):
        record.update(t=.5, passes=record["distance_lower"] > radius)
    online_keys = ("reduced_forms", "interior_binary64_selection", "strict_endpoints", "strict_port_selected")
    timings["port_portfolio_online"] = {
        key: sum(timings[stage][key] for stage in online_keys) for key in ("wall_seconds", "cpu_seconds")
    }
    timings["total"] = elapsed(beginning)
    return {
        "object_id": object_id, "noise": noise, "cap": CAP, "radius": radius,
        "endpoint_0": ends[0], "endpoint_2": ends[1],
        "best_endpoint": {"index": (0, 2)[best_end], "distance_lower": endpoint_lower,
                          "distance_upper": endpoint_upper, "passes": endpoint_lower > radius},
        "projected_best_interior": selected,
        "best_projected_portfolio": {
            "t": best[0] if chosen_lower > endpoint_lower else (0., 1.)[best_end],
            "distance_lower": best_portfolio, "passes": best_portfolio > radius,
        },
        "projected_midpoint": half, "native_midpoint": native_half, "affine_midpoint": affine_half,
        "continuous_segment_upper": segment_upper,
        "continuous_segment_upper_squared": max(r["squared_upper"] for r in (left_roof, native_half, right_roof)),
        "continuous_segment_at_or_below_radius": segment_upper <= radius,
        "continuous_segment_proof": {
            "left_extrapolation": left_roof, "midpoint": native_half, "right_extrapolation": right_roof,
        },
        "midpoint_energy": {
            "projected_jensen_trace": port_half_trace,
            "exact_jensen_trace": exact_gain,
            "missing_trace_from_chosen_port_trial": remainder,
            "captured_fraction": fraction,
        },
        "comparisons": {
            "projected_midpoint_lower_le_native_midpoint_upper": True,
            "certified_grid_lowers_le_continuous_segment_upper": True,
            "projected_strictly_better_than_both_endpoints": chosen_lower > endpoint_upper,
            "native_midpoint_strictly_better_than_both_endpoints": native_half["distance_lower"] > endpoint_upper,
            "projected_new_crossing_over_endpoints": best_portfolio > radius and endpoint_upper <= radius,
            "native_midpoint_new_crossing_over_endpoints": native_half["distance_lower"] > radius and endpoint_upper <= radius,
        },
        "interior_grid": grid_records,
        "timing": timings,
    }


def read_observations(path, *, validation=False):
    with np.load(path, allow_pickle=False) as data:
        y = data["Y"]
        radius = data["radius"]
        ids = data["object_id"]
        if validation:
            good = np.flatnonzero(data["split"] == "validation")
            noise = data["noise_fraction"].tolist().index(.02)
            if len(good) != 12:
                raise ValueError("Expected twelve archived independent validation objects.")
            return [(str(ids[i]), 2., y[i, noise], radius[i, noise]) for i in good]
        noise = data["noise_pct"] if "noise_pct" in data.files else 100 * data["noise_fraction"]
        if y.ndim == 3:
            if len(ids) != len(y) or radius.shape != (len(y),):
                raise ValueError("Flat observations require one ID and one radius per row.")
            noise = np.broadcast_to(noise, (len(y),))
            return [(str(a), float(b), c, d) for a, b, c, d in zip(ids, noise, y, radius, strict=True)]
        if y.ndim == 4 and radius.shape == y.shape[:2] and len(ids) == len(y):
            noise = np.broadcast_to(noise, radius.shape)
            return [(str(ids[i]), float(noise[i, j]), y[i, j], radius[i, j])
                    for i, j in np.ndindex(radius.shape)]
    raise ValueError("Unknown frozen observation schema.")


def historical_observations():
    from holdout import observations

    rows = [
        (object_id, noise, y, radius)
        for supplemental in (False, True)
        for object_id, noise, radius, y, caps in observations(supplement=supplemental)
        if CAP in caps
    ]
    if len(rows) != 84 or len({(oid, noise) for oid, noise, _, _ in rows}) != 84:
        raise ValueError("Exposure stress data must contain all 84 historical views.")
    return rows


def run(path, output, *, validation=False, historical=False, cache=DEFAULT_CACHE, rebuild=False):
    begin = started()
    setup = prepare(cache, rebuild=rebuild)
    setup_use = elapsed(begin)
    rows = []
    observations = historical_observations() if historical else read_observations(path, validation=validation)
    for oid, noise, y, radius in observations:
        row = evaluate(setup, y, radius, object_id=oid, noise=noise)
        rows.append(row)
        print(oid, noise, "endpoint", row["best_endpoint"]["distance_lower"], "port",
              row["best_projected_portfolio"]["distance_lower"], "mid",
              row["native_midpoint"]["distance_lower"], "radius", radius, flush=True)
    timings = {key: sum(row["timing"][key]["wall_seconds"] for row in rows)
               for key in rows[0]["timing"]} if rows else {}
    scope = "frozen fresh follow-up"
    if validation:
        scope = "validation only, existing 2 percent view per material"
    if historical:
        scope = "historical exposure stress, includes archived anchor donors"
    result = {
        "scope": scope,
        "input": "src/holdout.py primary and supplemental historical views" if historical else relative_receipt(path),
        "grid": GRID,
        "setup_cache": relative_receipt(cache),
        "setup_use": setup_use,
        "archived_setup": setup.receipt,
        "distinct_objects": len({row["object_id"] for row in rows}),
        "condition_count": len(rows),
        "evaluation_totals_wall_seconds": timings,
        "evaluation_totals_cpu_seconds": {
            key: sum(row["timing"][key]["cpu_seconds"] for row in rows) for key in timings
        },
        "rows": rows,
        **elapsed(begin),
    }
    atomic_json(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--validation", action="store_true")
    source.add_argument("--observations", type=Path)
    source.add_argument("--historical", action="store_true")
    source.add_argument("--prepare", action="store_true")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rebuild", action="store_true")
    opts = parser.parse_args()
    if opts.prepare:
        prepare(opts.cache, rebuild=opts.rebuild)
        return
    observations = PILOT / "data/synthetic_observations.npz" if opts.validation else opts.observations
    kind = "validation" if opts.validation else "historical" if opts.historical else "fresh"
    output = opts.output or HERE / "results" / f"{kind}_jensen.json"
    run(observations, output, validation=opts.validation, historical=opts.historical,
        cache=opts.cache, rebuild=opts.rebuild)


if __name__ == "__main__":
    main()
