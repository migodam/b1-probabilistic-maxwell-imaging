"""Data-only multimodal proposal branches on an exact active restriction of a 24² VIE."""

from collections import Counter
import argparse
import json
import sys

import numpy as np
from flint import ctx
from scipy.ndimage import gaussian_filter
from scipy.linalg import svd
from threadpoolctl import threadpool_limits

from pool_common import RESULTS, FLOW, FullWave, clock, elapsed, write_arrays, write_json
from box_certificates import Bounds

OUT = RESULTS / "structured_results"
PUBLIC = OUT / "public"
PRIVATE = OUT / "offline_input"
SEEDS = (118911, 230735, 229973, 368533)
RX = np.arange(16)
TX = np.array([0, 3, 6, 9])
SIDE = 24
PATCH = (slice(9, 16), slice(9, 16))
CAP = 5.5
N_STARTS = 4
TOP_K = (1, 3, 6)
MAX_CALLS = 42
BRANCHES = (
    "coverage_control", "anchor_local", "low_peak", "lower_peak", "high_peak", "higher_peak",
    "shift_right", "shift_down", "shift_left", "shift_up",
    "ring", "two_peaks", "ridge", "ramp", "plateau_near_3.5", "broad_near_3.5",
)


def filenames():
    for seed in SEEDS:
        yield seed, f"guided_tm24_{seed}__pw90x4__rx180x16__n0125"


def create_observations():
    if (PUBLIC / "geometry.npz").exists():
        raise FileExistsError(PUBLIC / "geometry.npz")
    PUBLIC.mkdir(parents=True, exist_ok=True)
    PRIVATE.mkdir(exist_ok=True)
    started = clock()
    sys.path.insert(0, str(FLOW / "src"))
    import physics

    geometry = physics.Solver(n=SIDE, k=4*np.pi, n_rx=32, n_illum=12)
    active = np.arange(SIDE**2).reshape(SIDE, SIDE)[PATCH].ravel()
    g, s, b = geometry.GD[np.ix_(active, active)], geometry.GS[:, active], geometry.incident[active]
    write_arrays(PUBLIC / "geometry.npz", green=g, sensor=s, incident=b, active=active,
                 points=geometry.points, receivers=geometry.receiver_pos, directions=geometry.illum_dir,
                 train_rx=RX, train_tx=TX, k=geometry.k, side=geometry.n)
    pos = (np.arange(7) - 3.) / 3
    yy, xx = np.meshgrid(pos, pos, indexing="ij")
    costs = {}
    for number, (seed, case) in enumerate(filenames()):
        begin = clock()
        rng = np.random.default_rng(seed)
        floor = rng.uniform(2.05, 2.45)
        theta = rng.uniform(-.45, .45)
        x_rot = np.cos(theta) * xx + np.sin(theta) * yy - rng.uniform(-.14, .14)
        y_rot = -np.sin(theta) * xx + np.cos(theta) * yy - rng.uniform(-.14, .14)
        if number == 0:
            shape = np.exp(-.5*((x_rot/.60)**2+(y_rot/.52)**2))
        elif number == 1:
            shape = np.exp(-((np.sqrt((x_rot/.89)**2+y_rot**2)-.65)/.27)**2)
        elif number == 2:
            shape = np.maximum(np.exp(-((x_rot-.47)**2+(y_rot-.24)**2)/.19),
                               np.exp(-((x_rot+.45)**2+(y_rot+.29)**2)/.15))
        else:
            shape = .5 + .5 * np.tanh(3*(x_rot + .38*y_rot))
        peak = rng.uniform(3.85, 4.65)
        material = np.asarray(floor + (peak - floor) * shape / shape.max()).ravel()
        model = FullWave(g, s, b, np.ones((32, 12)))
        clean = model.evaluate(material)[0]
        clean_train = clean[np.ix_(RX, TX)]
        perturbations = rng.standard_normal(clean_train.shape) + 1j*rng.standard_normal(clean_train.shape)
        raw_radius = .0125*np.linalg.norm(clean_train)
        observed = clean_train + perturbations*(raw_radius/np.linalg.norm(perturbations))
        radius = float(raw_radius + 1e-9*np.linalg.norm(clean_train))
        write_arrays(PUBLIC / f"{case}.npz", train=observed, radius=radius)
        write_arrays(PRIVATE / f"{case}.npz", material=material, noiseless_ports=clean,
                     raw_radius=raw_radius)
        costs[case] = dict(generation_operations=dict(model.counts), **elapsed(begin))
    write_json(OUT / "generation_costs.json", dict(cases=costs, elapsed=elapsed(started)))


def proposals(g, s, b, y, model, seed):
    start = clock()
    born = (s[:, None, :] * b.T[None, :, :]).reshape(y.size, g.shape[0])
    design = np.concatenate((born.real, born.imag))
    target = np.concatenate((y.real.ravel(), y.imag.ravel()))
    u, q, vh = svd(design, full_matrices=False, check_finite=False)
    use = q > q[0] * .04
    baseline = np.clip(vh[use].T @ ((u[:, use].T @ target) / q[use]), 0., CAP)
    anchor, _, anchor_cost = model.refine(baseline, np.zeros_like(baseline),
                                         np.full_like(baseline, CAP), max_calls=32)
    template = gaussian_filter(anchor.reshape(7, 7), .6, mode="nearest")
    pos = (np.arange(7) - 3.) / 3
    yy, xx = np.meshgrid(pos, pos, indexing="ij")
    m = np.median(template)
    shifted = [np.roll(template, step, axis=axis) for step, axis in
               [(2, 1), (2, 0), (-2, 1), (-2, 0)]]
    ring = m + .95 * (np.exp(-((np.hypot(xx, yy)-.70)/.24)**2)-.35)
    peaks = m + 1.15*(np.maximum(np.exp(-((xx-.50)**2+(yy+.18)**2)/.18),
                                 np.exp(-((xx+.50)**2+(yy-.18)**2)/.18))-.25)
    ridge = m + .8*(np.exp(-(yy/.35)**2)-.44)
    ramp = m + (xx - .55*yy) * .75
    centers = [template, template, .36*template+.35, .65*template+.30,
               1.23*template+.20, 1.48*template+.20, *shifted,
               ring, peaks, ridge, ramp, np.full_like(template, 3.5), np.full_like(template, 3.5)]
    radii = (CAP, .15, .36, .57, .63, .82, .42, .42, .42, .42, .43, .45, .45, .45, .70, 1.50)
    rng = np.random.default_rng(seed)
    samples, lows, highs = [], [], []
    for center, width in zip(centers, radii, strict=True):
        center = np.clip(center, .03, CAP-.03)
        lo = np.maximum(center-width, 0.)
        hi = np.minimum(center+width, CAP)
        for j in range(N_STARTS):
            direction = gaussian_filter(rng.normal(size=center.shape), .85)
            direction /= np.max(np.abs(direction))
            value = np.clip(center + float(j > 0) * min(width*.72, .75)*direction, lo, hi)
            samples.append(value.ravel())
        lows.append(lo.ravel())
        highs.append(hi.ravel())
    return np.array(samples), np.array(lows), np.array(highs), dict(
        center_full_wave_fit=anchor_cost, effective_born_svd_rank=int(use.sum()), **elapsed(start))


def model_inputs():
    with np.load(PUBLIC / "geometry.npz", allow_pickle=False) as ge:
        return ge["green"], ge["sensor"], ge["incident"]


def method_cost(prep, screening, rows, stage_counts):
    count = Counter(prep["center_full_wave_fit"]["counts"])
    count.update(stage_counts)
    for row in rows:
        count.update(row["fit"]["counts"])
    return dict(count=count,
                wall_seconds=prep["wall_seconds"]+screening["wall_seconds"]
                + sum(r["fit"]["wall_seconds"] for r in rows),
                cpu_seconds=prep["cpu_seconds"]+screening["cpu_seconds"]
                + sum(r["fit"]["cpu_seconds"] for r in rows))


def online():
    if (OUT / "online.json").exists():
        raise FileExistsError("Online predictions already exist.")
    begun = clock()
    setup_started = clock()
    g, all_s, all_b = model_inputs()
    s, b = all_s[RX], all_b[:, TX]
    ctx.prec, ctx.threads = 112, 1
    # Full row/column rank implies both unrestricted support port ranges are surjective.
    ranks = [int(np.linalg.matrix_rank(s)), int(np.linalg.matrix_rank(b))]
    setup = dict(**elapsed(setup_started), receiver_rank=ranks[0], incident_rank=ranks[1])
    result = {}
    for seed, case in filenames():
        with np.load(PUBLIC / f"{case}.npz", allow_pickle=False) as data:
            y, radius = data["train"], float(data["radius"])
        result[case] = online_case(g, s, b, y, radius, seed+1007, case)
    write_json(OUT / "online.json", dict(cases=result, geometry_setup=setup,
                                        total_run=elapsed(begun), verified_model="exact stored binary64"))


def online_case(g, s, b, y, radius, proposal_seed, case):
    model = FullWave(g, s, b, y)
    samples, lo, hi, prep = proposals(g, s, b, y, model, proposal_seed)
    starts = samples[::N_STARTS]
    if not np.all((samples >= np.repeat(lo, N_STARTS, axis=0)) & (samples <= np.repeat(hi, N_STARTS, axis=0))):
        raise ValueError("Sample generation escaped a proposal box.")

    point_start = clock()
    counts_before = model.counts.copy()
    point = np.stack([model.evaluate(m)[0] for m in samples])
    point_cost = elapsed(point_start)
    point_counts = model.counts - counts_before
    point_residuals = np.linalg.norm((point-y).reshape(len(samples), -1), axis=1)
    best_points = np.min(point_residuals.reshape(len(BRANCHES), N_STARTS), axis=1)
    ranking = np.argsort(best_points, kind="stable").tolist()
    bounds = Bounds(g, s, b, y)
    rows, refined, predictions = [], [], []
    for branch, center, lower_bound, upper_bound, best_point in zip(BRANCHES, starts, lo, hi, best_points, strict=True):
        neumann = bounds.neumann(lower_bound, upper_bound, radius)
        shifted = bounds.shifted(lower_bound, upper_bound, radius)
        m, predicted, fit = model.refine(center, lower_bound, upper_bound, MAX_CALLS)
        referee = clock()
        strict = bounds.point_upper(m)
        fit["independent_endpoint_check"] = dict(upper=strict, within_radius=strict<=radius, **elapsed(referee))
        if strict <= radius and max(neumann["directional_lower"], shifted["strict_lower"]) > radius:
            raise ArithmeticError("A feasible physical endpoint was incorrectly excluded.")
        rows.append(dict(branch=branch, best_sample_residual=float(best_point),
                         support_lower=0., neumann=neumann, shifted=shifted, fit=fit))
        refined.append(m)
        predictions.append(predicted)
    survivors = [c for c, row in enumerate(rows) if row["shifted"]["strict_lower"]<=radius]
    neum_rows = [r["neumann"] for r in rows]
    shift_rows = [r["shifted"] for r in rows]
    branches = list(range(len(BRANCHES)))
    complete = dict(wall_seconds=0., cpu_seconds=0.)
    branch_ids = {
        "all_cluster": (branches, complete, {}),
        "all_support_then_refine": (branches, complete, {}),
        "neumann_standard": ([c for c, r in enumerate(neum_rows) if r["standard_lower"]<=radius],
                             sum_costs(neum_rows), sum_counts(neum_rows, neumann_counts)),
        "neumann_directional": ([c for c, r in enumerate(neum_rows) if r["directional_lower"]<=radius],
                                sum_costs(neum_rows), sum_counts(neum_rows, neumann_counts)),
        "shifted": (survivors, sum_costs(shift_rows), {}),
    }
    branch_ids.update({f"top_{k}": (ranking[:k], point_cost, point_counts) for k in TOP_K})
    branch_ids["top_matched_count"] = (ranking[:len(survivors)], point_cost, point_counts)
    pipelines = {}
    for name, (visited, screening, stage_counts) in branch_ids.items():
        visited_rows = [rows[i] for i in visited]
        winner = min(visited, key=lambda c: rows[c]["fit"]["residual"], default=None)
        costs = method_cost(prep, screening, visited_rows, stage_counts)
        pipelines[name] = dict(winner=winner, branches=len(visited), visit=visited, **costs,
                               best_endpoint=rows[winner]["fit"]["residual"] if winner is not None else None,
                               strict_feasible=(rows[winner]["fit"]["independent_endpoint_check"]["within_radius"]
                                                if winner is not None else False))
    write_arrays(OUT / f"{case}.npz", train=y, radius=radius, samples=samples, lower=lo, upper=hi,
                 sample_predictions=point, refined=np.asarray(refined), train_predictions=np.asarray(predictions))
    print(f"{case}: shift {len(survivors):2d}/16, Neumann {pipelines['neumann_directional']['branches']:2d}/16, "
          f"shift actual min/r {pipelines['shifted']['best_endpoint']/radius:.3f}", flush=True)
    return dict(radius=radius, proposal_seed=proposal_seed, proposal_setup=prep, branch_records=rows,
                physical_point_ranking=ranking, pipelines=pipelines, point_screen_cost=point_cost,
                point_samples=len(samples), executed_counts=dict(model.counts))


def sum_costs(rows):
    return {key: sum(row[key] for row in rows) for key in ("cpu_seconds", "wall_seconds")}


def neumann_counts(row):
    return dict(lu=row["center_factorizations"], forward_rhs=row["center_incident_rhs"],
                extra_center_transfer_rhs=row["center_transfer_rhs"],
                outward_center_inverse_factors=row["verified_inverse_factors"])


def sum_counts(rows, key):
    result = Counter()
    for row in rows:
        result.update(key(row))
    return result


def offline():
    online_result = json.loads((OUT / "online.json").read_text())
    g, s, b = model_inputs()
    rx_test = np.setdiff1d(np.arange(s.shape[0]), RX)
    tx_test = np.setdiff1d(np.arange(b.shape[1]), TX)
    report, arrays = {}, {}
    for _, case in filenames():
        row = online_result["cases"][case]
        with np.load(PRIVATE / f"{case}.npz", allow_pickle=False) as truth_data:
            truth, clean = truth_data["material"], truth_data["noiseless_ports"]
        with np.load(OUT / f"{case}.npz", allow_pickle=False) as predictions:
            refined = predictions["refined"]
            contains = np.all((predictions["lower"]<=truth) & (truth<=predictions["upper"]), axis=1)
        r = row["radius"]
        with np.load(PUBLIC / f"{case}.npz", allow_pickle=False) as train:
            witnessed = Bounds(g, s[RX], b[:, TX], train["train"]).point_upper(truth)
        if witnessed > r:
            raise ArithmeticError("Generation failed the frozen physical feasibility check.")
        model = FullWave(g, s, b, clean)
        private_scores = []
        for i, (m, own) in enumerate(zip(refined, row["branch_records"], strict=True)):
            estimate = model.evaluate(m)[0]
            arrays[f"{case}__prediction_{i:02d}"] = estimate
            private_scores.append(dict(
                branch=own["branch"], truth_in_box=bool(contains[i]),
                full_image_relative_l2=float(np.linalg.norm(m-truth)/np.linalg.norm(truth)),
                inside_relative_l2=float(np.linalg.norm(m-truth)/np.linalg.norm(truth)),
                background_rms=0.,
                heldout_receiver_relative_l2=residual(estimate, clean, rx_test, TX),
                heldout_source_relative_l2=residual(estimate, clean, RX, tx_test)))
            if contains[i] and own["shifted"]["strict_lower"] > r:
                raise ArithmeticError("Shifted excluded a strict feasible truth-contained control.")
        report[case] = dict(truth_peak=float(max(truth)), truth_verified_data_upper=witnessed,
                            true_containment_count=int(contains.sum()),
                            count_background_cells=SIDE**2-len(truth),
                            metrics=private_scores, pipelines={
                                name: private_scores[values["winner"]] if values["winner"] is not None else None
                                for name, values in row["pipelines"].items()
                            }, offline_forward_count=dict(model.counts))
    write_arrays(OUT / "offline_full_ports.npz", **arrays)
    write_json(OUT / "offline.json", report)


def residual(predicted, truth, receivers, illuminations):
    ids = np.ix_(receivers, illuminations)
    return float(np.linalg.norm(predicted[ids]-truth[ids])/np.linalg.norm(truth[ids]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["generate", "online", "offline"])
    phase = parser.parse_args().phase
    with threadpool_limits(limits=1):
        {"generate": create_observations, "online": online, "offline": offline}[phase]()
