"""Data-selected restriction with independent coarse and fine full-wave solves."""

import argparse
import json
from pathlib import Path
import resource
from time import perf_counter, process_time

import numpy as np
from scipy.ndimage import binary_dilation
from scipy.sparse import csr_matrix, kron

from run_inversion import PENALTY, pair_centers, score
from src.designs import Design, all_designs, read_arm
from src.inversion import _adjacent, _bound_step, energy, fullwave
from src.ndt import BOUNDS, Host, grid, host

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SIDE = 24
ARMS = ("f1_t8_r8", "f2_t4_r8", "f4_t2_r8", "f8_t1_r8", "high_t8_r8",
        "loc_nominal", "f4_t8_r16", "oneside_f2_t4_r8")
DEVELOPMENT_CHECKPOINTS = (6, 12, 24)


def region(initial):
    absolute = abs(initial)
    threshold = max(.12, .15 * float(absolute.max()))
    mask = binary_dilation(absolute >= threshold, structure=np.ones((3, 3)))
    return np.flatnonzero(mask.ravel()), threshold


def hosts_for(n, design, cache_costs):
    for f in design.freqs:
        key = (n, f)
        if key not in cache_costs:
            t0, c0 = perf_counter(), process_time()
            model = host(n, f)
            cache_costs[key] = dict(
                resolution=n, frequency_ghz=f, wall_seconds=perf_counter() - t0,
                cpu_seconds=process_time() - c0, retained_host_bytes=sum(a.nbytes for a in model),
                plate_lu=1, complex_rhs=n * n + model.Eh.shape[1] + model.Sh.shape[0],
            )
    return [host(n, f) for f in design.freqs]


def restricted_models(hosts, n, design, active):
    stride = n // SIDE
    # Each material variable owns these distinct fine polarization currents.
    groups = np.arange(n*n).reshape(SIDE, stride, SIDE, stride)
    groups = groups.transpose(0, 2, 1, 3).reshape(SIDE * SIDE, stride**2)[active]
    cells = groups.ravel()
    selected = []
    for h in hosts:
        selected.append(Host(h.Eh[np.ix_(cells, design.tx)], h.Sh[np.ix_(design.rx, cells)],
                             h.Gh[np.ix_(cells, cells)], h.xy[cells]))
    return selected


def overlapped_models(hosts, n, design, active):
    # Area averaging supplies a fixed, non-nested comparison mesh.
    coarse_edges = np.arange(SIDE + 1) / SIDE
    fine_edges = np.arange(n + 1) / n
    left = np.maximum(fine_edges[:-1, None], coarse_edges[None, :-1])
    right = np.minimum(fine_edges[1:, None], coarse_edges[None, 1:])
    fraction = csr_matrix(np.maximum(right - left, 0.) * n)
    tie = kron(fraction, fraction, format="csr")[:, active]
    cells = np.flatnonzero(np.diff(tie.indptr))
    selected = [
        Host(h.Eh[np.ix_(cells, design.tx)], h.Sh[np.ix_(design.rx, cells)],
             h.Gh[np.ix_(cells, cells)], h.xy[cells])
        for h in hosts
    ]
    return selected, tie[cells]


def expand_material(x, n, tie=None):
    return np.repeat(x, (n // SIDE)**2) if tie is None else tie @ x


class RestrictedPenalty:
    def __init__(self, active):
        self.active = active
        inverse = np.full(SIDE**2, -1)
        inverse[active] = np.arange(len(active))
        self.a, self.b = _adjacent(SIDE)
        self.ja, self.jb = inverse[self.a], inverse[self.b]
        self.edge_a = self.ja >= 0
        self.edge_b = self.jb >= 0
        self.internal = self.edge_a & self.edge_b
        self.full = np.zeros(SIDE**2)
        self.diag = np.diag_indices(len(active))

    def evaluate(self, x, hessian=None):
        self.full[self.active] = x
        value, grad = PENALTY.evaluate(self.full)
        if hessian is not None:
            hessian[self.diag] += (PENALTY.amplitude / self.full.size
                                   / np.hypot(x, PENALTY.smooth))
            edge_weights = (PENALTY.tv / self.a.size /
                            np.hypot(self.full[self.a] - self.full[self.b], PENALTY.smooth))
            da, db = self.ja[self.edge_a], self.jb[self.edge_b]
            np.add.at(hessian, (da, da), edge_weights[self.edge_a])
            np.add.at(hessian, (db, db), edge_weights[self.edge_b])
            i, j = self.ja[self.internal], self.jb[self.internal]
            hessian[i, j] -= edge_weights[self.internal]
            hessian[j, i] -= edge_weights[self.internal]
        return value, grad[self.active]


def refine(models, x, active, y, reference, n, iterations, tie=None):
    start = perf_counter()
    stride_sq = (n // SIDE)**2
    penalty = RestrictedPenalty(active)
    weight = 1 / (reference * np.sqrt(y.size))
    counters = dict(factorizations=0, actual_complex_rhs=0, forward_evaluations=0, jacobian_evaluations=0)
    history, snapshots = [], {}

    def forward(material, jacobian=False):
        counters["forward_evaluations"] += 1
        counters["jacobian_evaluations"] += int(jacobian)
        if np.any(material):
            counters["factorizations"] += len(models)
            counters["actual_complex_rhs"] += sum(h.Eh.shape[1] + jacobian * h.Sh.shape[0] for h in models)
        return fullwave(models, expand_material(material, n, tie), jacobian=jacobian)

    if not len(active):
        return x, history, snapshots, counters, perf_counter() - start, "empty_region"

    stop_reason = "iteration_cap"
    for iteration in range(iterations):
        predicted, subcell_jac = forward(x, jacobian=True)
        # Forward remains a full block state; only the material derivative sums subcells.
        jac = (subcell_jac.reshape(y.size, len(active), stride_sq).sum(axis=2)
               if tie is None else subcell_jac @ tie)
        residual = (predicted - y).ravel() * weight
        jac_real = np.ascontiguousarray(np.r_[jac.real, jac.imag] * weight)
        residual_real = np.r_[residual.real, residual.imag]
        hessian = jac_real.T @ jac_real
        reg, gradient_reg = penalty.evaluate(x, hessian)
        objective = 0.5 * energy(residual_real) + reg
        gradient = jac_real.T @ residual_real + gradient_reg
        step = _bound_step(hessian, gradient, x)
        if np.max(abs(step)) < 5e-5:
            stop_reason = "small_step"
            break
        directional = gradient @ step
        for power in range(9):
            alpha = 0.5 ** power
            candidate = np.clip(x + alpha * step, *BOUNDS)
            new_prediction = forward(candidate)
            new_data = 0.5 * energy((new_prediction - y) * weight)
            new_objective = new_data + penalty.evaluate(candidate)[0]
            if new_objective <= objective + 1e-4 * alpha * directional:
                break
        else:
            stop_reason = "line_search"
            break
        x = candidate
        it = iteration + 1
        history.append(dict(iteration=it, wall_seconds=perf_counter() - start, objective=new_objective,
                            data_loss=new_data, accepted_step=alpha))
        if it in DEVELOPMENT_CHECKPOINTS:
            snapshots[it] = x.copy()
        if objective - new_objective <= 1e-8 * max(1., objective):
            stop_reason = "small_decrease"
            break
    return x, history, snapshots, counters, perf_counter() - start, stop_reason


def load_observations(path):
    with np.load(path, allow_pickle=False) as incoming:
        return {key: incoming[key] for key in (
            "reference_amplitude", "frequencies_ghz", "case_ids", "splits", "clean", "noise"
        )}


def read_seed_variants(include_sparse, only_sparse, seed_tag, designs):
    images, costs = {}, {}
    if not only_sparse:
        for suffix in ("", "_task_design"):
            path = RESULTS / f"inversion_materials{suffix}.npz"
            with np.load(path, allow_pickle=False) as source:
                saved = source["estimates"]
                for i, arm in enumerate(source["arms"]):
                    for j, case in enumerate(source["case_ids"]):
                        for k, method in enumerate(source["methods"]):
                            estimate = saved[i, j, k]
                            if np.isfinite(estimate).all():
                                images[(str(arm), str(case), str(method))] = estimate
            with (RESULTS / f"inversion_metrics{suffix}.json").open() as f:
                for row in json.load(f)["rows"]:
                    costs[(row["arm"], row["case"], row["method"])] = row
    if include_sparse or only_sparse:
        with np.load(RESULTS / f"sparse_materials{seed_tag}.npz", allow_pickle=False) as source:
            saved = source["images"]
            for i, arm in enumerate(source["names"]):
                for j, case in enumerate(source["case_ids"]):
                    estimate = saved[i, j]
                    if np.isfinite(estimate).all():
                        images[(str(arm), str(case), "sparse_born")] = estimate
        with (RESULTS / "task_fields_cost.json").open() as f:
            field_cost = {row["f"]: row["wall_s"] for row in json.load(f)["timing"] if row["grid"] == 48}
        sparse_path = RESULTS / f"sparse_cost{seed_tag}.json"
        if sparse_path.is_file():
            with sparse_path.open() as f:
                cost_by_arm = {row["arm"]: row for row in json.load(f) if row["grid"] == 48}
            for key in images:
                arm, _, method = key
                if method != "sparse_born" or arm not in cost_by_arm:
                    continue
                seed_cost = cost_by_arm[arm]
                costs[key] = dict(
                    shared_host_seconds=seed_cost["preparation_wall_s"]
                                        + sum(field_cost[f] for f in designs[arm].freqs),
                    solver_seconds=seed_cost["fit_wall_s_per_material"],
                    timing_scope="Per-case mean batch sparse fit; archived 48-grid derivative-bank and per-design preparation added cold.",
                    arm_batch_seconds=seed_cost["batched_fit_wall_s"],
                )
    return images, costs


def design_options():
    options = all_designs()
    with (RESULTS / "task_design.json").open() as f:
        for raw in json.load(f)["designs"]:
            options[raw["name"]] = Design(raw["name"], tuple(raw["freqs"]), tuple(raw["tx"]), tuple(raw["rx"]))
    return options


def compare(args):
    inputs = load_observations(args.data)
    cases = [str(x) for x in inputs["case_ids"]]
    reference = float(inputs["reference_amplitude"])
    options = design_options()
    seeds, costs = read_seed_variants(args.include_sparse, args.only_sparse, args.seed_tag, options)
    positions = grid(SIDE) * 1000
    cache_costs, rows, estimates, masks, development = {}, [], [], [], []
    development_only = args.dev_only
    resolutions = (24, 36, 48) if args.grid_sensitivity else (24, 48)
    variant_names = ("seed", "masked_seed", *(f"roi{n}" for n in resolutions))
    total_start = perf_counter()

    for arm in args.arms:
        design = options[arm]
        available = [
            key for key in seeds
            if key[0] == arm and key[1] in cases and (not development_only or key[1].startswith("dev_"))
            and (args.cases is None or key[1] in args.cases)
        ]
        if not available:
            continue
        both_hosts = {n: hosts_for(n, design, cache_costs) for n in resolutions}
        cold = {n: sum(cache_costs[n, f]["wall_seconds"] for f in design.freqs) for n in both_hosts}
        print(f"{arm} cold host " + " / ".join(f"{n}: {cold[n]:.3f}s" for n in resolutions), flush=True)
        for key in available:
            _, case, method = key
            y, sigma = read_arm(inputs, cases.index(case), design)
            t0 = perf_counter()
            seed = seeds[key]
            active, threshold = region(seed)
            masked = np.zeros((SIDE, SIDE))
            masked.ravel()[active] = seed.ravel()[active]
            mask = np.zeros(SIDE**2, dtype=bool)
            mask[active] = True
            selector_seconds = perf_counter() - t0
            images = [seed, masked]
            variants = {}
            snapshots = {}
            first_costs = costs.get(key, {})
            initial_solver_time = first_costs.get("solver_seconds")
            initial_host_time = first_costs.get("shared_host_seconds")
            row = dict(
                arm=arm, case=case, initializer=method, split=str(inputs["splits"][cases.index(case)]),
                frequencies_ghz=list(design.freqs), read_count=design.samples, shots=design.shots,
                noise_sd=sigma, threshold_eps=threshold, coarse_roi_cells=len(active),
                roi_fraction=len(active)/(SIDE**2), selector_seconds=selector_seconds,
                seed_solver_seconds=initial_solver_time, seed_host_cold_seconds=initial_host_time,
                seed_original_timing=first_costs if method == "sparse_born" else None,
            )
            for n, models in both_hosts.items():
                t1 = perf_counter()
                if n == 36:
                    restricted, tie = overlapped_models(models, n, design, active)
                else:
                    restricted, tie = restricted_models(models, n, design, active), None
                extract_seconds = perf_counter() - t1
                fitted, trace, states, counted, fit_seconds, stop_reason = refine(
                    restricted, masked.ravel()[active], active, y, reference, n, args.iterations, tie,
                )

                image = np.zeros((SIDE, SIDE))
                image.ravel()[active] = fitted
                images.append(image)
                if development_only:
                    snapshots[n] = {i: states.get(i, fitted) for i in DEVELOPMENT_CHECKPOINTS}
                # Material quality is scored after all estimation has finished.
                predicted_init = fullwave(restricted, expand_material(masked.ravel()[active], n, tie))
                predicted_fit = fullwave(restricted, expand_material(fitted, n, tie))
                residual_start = float(np.sqrt(np.mean(abs(predicted_init - y)**2)))
                residual_end = float(np.sqrt(np.mean(abs(predicted_fit - y)**2)))
                online = selector_seconds + extract_seconds + fit_seconds
                cached = online + initial_solver_time if initial_solver_time is not None else None
                # Avoid charging the same coarse plate preparation both to the seed and 24-grid correction.
                prep_n = cold[n] if len(active) and (n != 24 or method == "sparse_born") else 0.
                first_case_total = cached + initial_host_time + prep_n if (
                    cached is not None and initial_host_time is not None
                ) else None
                variants[f"roi{n}"] = dict(
                    **counted, accepted_iterations=len(trace), stop_reason=stop_reason,
                    solved_cells=restricted[0].Eh.shape[0], extract_seconds=extract_seconds,
                    fit_seconds=fit_seconds, correction_seconds=online, cold_host_seconds=cold[n],
                    end_to_end_cached_seconds=cached, end_to_end_first_case_seconds=first_case_total,
                    initial_residual_over_sigma=residual_start/sigma, final_residual_over_sigma=residual_end/sigma,
                    history=trace,
                )
            row["solvers"] = variants
            rows.append(row)
            estimates.append(images)
            masks.append(mask.reshape(SIDE, SIDE))
            development.append(snapshots)
            print(f"{arm:>20s} {case:>16s} {method:>12s} cells {len(active):3d} "
                  f"24/48 {variants['roi24']['fit_seconds']:.3f}/{variants['roi48']['fit_seconds']:.3f}s "
                  f"stop {variants['roi48']['stop_reason']}", flush=True)

    label = args.suffix
    np.savez_compressed(RESULTS / f"roi_materials{label}.npz", estimates=np.array(estimates),
                        masks=np.array(masks), arms=[r["arm"] for r in rows], cases=[r["case"] for r in rows],
                        initializers=[r["initializer"] for r in rows], variants=variant_names)
    with np.load(args.truth_data or args.data, allow_pickle=False) as f:
        truth = f["truth"]
    if args.data == RESULTS / "acquisition.npz":
        scored_pairs = pair_centers()
    elif args.truth_metadata:
        with args.truth_metadata.open() as stream:
            scenes = json.load(stream)["case_definitions"]
        scored_pairs = {
            scene["id"]: np.array([[p["x_mm"], p["y_mm"]] for p in scene["primitives"]])
            for scene in scenes if scene.get("shape") == "two_air_holes" or scene["id"].startswith("test_pair")
        }
    else:
        scored_pairs = {}
    for row, images, mask, checkpoints in zip(rows, estimates, masks, development, strict=True):
        target = truth[cases.index(row["case"])].ravel()
        pair = scored_pairs.get(row["case"])
        row["scores"] = {label: score(image.ravel(), target, positions, pair)
                         for label, image in zip(variant_names, images, strict=True)}
        valid = abs(target) >= .20
        row["truth_support_coverage_for_scoring_only"] = float(np.mean(mask.ravel()[valid])) if valid.any() else None
        if development_only:
            row["development_checkpoints"] = {}
            for n, cp in checkpoints.items():
                row["development_checkpoints"][n] = {}
                for i, q in cp.items():
                    image = np.zeros(SIDE**2)
                    image[mask.ravel()] = q
                    row["development_checkpoints"][n][i] = score(image, target, positions, pair)
    payload = dict(
        rows=rows, host_builds=list(cache_costs.values()), acquisition=args.data.name, sparse_seed_tag=args.seed_tag,
        iteration_cap=args.iterations, penalty=vars(PENALTY), material_box=BOUNDS,
        rss_max_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        elapsed_seconds=perf_counter() - total_start,
    )
    (RESULTS / f"roi_metrics{label}.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    return payload


def check_operators():
    design = all_designs()["high_t8_r8"]
    model = host(48, 30.)
    rng = np.random.default_rng(151145)
    active = np.sort(rng.choice(SIDE**2, 30, replace=False))
    x = rng.uniform(-.9, .9, len(active))
    direction = rng.normal(size=len(active))
    direction /= np.linalg.norm(direction)
    restricted = restricted_models([model], 48, design, active)
    value, subjac = fullwave(restricted, np.repeat(x, 4), jacobian=True)
    derivative = subjac.reshape(value.size, len(x), 4).sum(axis=2) @ direction

    full = np.zeros(SIDE**2)
    full[active] = x
    expanded = np.repeat(np.repeat(full.reshape(SIDE, SIDE), 2, axis=0), 2, axis=1)
    full_model = Host(model.Eh[:, design.tx], model.Sh[design.rx, :], model.Gh, model.xy)
    reference = fullwave([full_model], expanded)
    errors = []
    for step in (.001, .0001, .00001, .000001):
        plus = fullwave(restricted, np.repeat(x + step * direction, 4))
        minus = fullwave(restricted, np.repeat(x - step * direction, 4))
        estimate = (plus - minus).ravel() / (2 * step)
        errors.append(dict(step=step, directional_derivative_relative_error=float(
            np.linalg.norm(estimate - derivative) / np.linalg.norm(derivative)
        )))
    full_hessian = np.zeros((SIDE**2, SIDE**2))
    full_penalty, full_gradient = PENALTY.evaluate(full, full_hessian)
    local_hessian = np.zeros((len(active), len(active)))
    local_penalty, local_gradient = RestrictedPenalty(active).evaluate(x, local_hessian)
    result = dict(
        frequency_ghz=30., selected_coarse_cells=len(active), selected_fine_cells=4 * len(active),
        forward_restriction_relative_error=float(np.linalg.norm(value - reference) / np.linalg.norm(reference)),
        gradient_centered_differences=errors,
        penalty_value_absolute_error=abs(local_penalty - full_penalty),
        penalty_gradient_max_error=float(np.max(abs(local_gradient - full_gradient[active]))),
        penalty_hessian_max_error=float(np.max(abs(local_hessian - full_hessian[np.ix_(active, active)]))),
    )
    (RESULTS / "roi_operators.json").write_text(json.dumps(result, indent=2) + "\n")


def make_figures():
    import matplotlib.pyplot as plt

    rows_by_source = []
    for suffix in ("", "_sparse"):
        with (RESULTS / f"roi_metrics{suffix}.json").open() as stream:
            rows_by_source.extend(json.load(stream)["rows"])
    cases = ("test_hole", "test_pair_close", "test_pair_wide", "test_mixed", "test_empty")
    names = ("Hole", "3 mm pair", "6.2 mm pair", "Mixed contrast", "Intact")
    main_arm = "high_t8_r8"
    with np.load(RESULTS / "roi_materials_sparse.npz", allow_pickle=False) as pack:
        estimates = {str(case): image for case, image, arm in
                     zip(pack["cases"], pack["estimates"], pack["arms"], strict=True) if arm == main_arm}
    with np.load(RESULTS / "acquisition.npz", allow_pickle=False) as pack:
        truth = dict(zip(pack["case_ids"], pack["truth"], strict=True))

    fig, axes = plt.subplots(len(cases), 5, figsize=(9, 9), layout="constrained")
    titles = ("Truth", "Sparse linear seed", "Masked seed", "24-grid fit", "48-grid fit")
    for i, (case, title) in enumerate(zip(cases, names, strict=True)):
        images = (truth[case], *estimates[case])
        row = next(r for r in rows_by_source if r["case"] == case and r["arm"] == main_arm
                   and r["initializer"] == "sparse_born")
        axes[i, 0].set_ylabel(title)
        for j, (ax, data) in enumerate(zip(axes[i], images, strict=True)):
            graphic = ax.imshow(data, extent=(-12, 12, -12, 12), origin="lower",
                                vmin=-1.5, vmax=1.5, cmap="RdBu_r", interpolation="nearest")
            if i == 0:
                ax.set_title(titles[j], fontsize=10)
            if j:
                metric = row["scores"][("seed", "masked_seed", "roi24", "roi48")[j - 1]]
                value = metric["support_iou"]
                label = f"IoU {value:.2f}" if value is not None else f"{metric['false_support_cells']} false cells"
                ax.text(.03, .03, label, transform=ax.transAxes, fontsize=8,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=1.5))
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle("30 GHz · 64 complex reads · 2% reference noise")
    fig.colorbar(graphic, ax=axes, shrink=.55, label="Relative permittivity deviation")
    save_figure(fig, "roi_reconstruction")
    plt.close(fig)

    arms = ("f1_t8_r8", "f2_t4_r8", "f4_t2_r8", "f8_t1_r8", "high_t8_r8",
            "loc_nominal", "f4_t8_r16", "oneside_f2_t4_r8")
    arms_display = ("16 GHz", "2 / 30", "4 tones", "8 tones", "30 GHz",
                    "28 / 30", "4 tones rich", "One side")
    fig, axes = plt.subplots(1, 3, figsize=(12.7, 3.85), layout="constrained", sharey=True)
    for ax, method, title in zip(axes, ("fft_tsom", "projected_gn", "sparse_born"),
                                ("FFT-TSOM seed", "Projected Gauss–Newton seed", "Sparse linear seed"), strict=True):
        rows = [r for r in rows_by_source if r["initializer"] == method and r["split"] == "independent"
                and r["case"] != "test_empty"]
        for name, label, color, marker in (("seed", "Original seed", "#797979", "o"),
                ("roi24", "24-grid correction", "#b25740", "s"),
                ("roi48", "48-grid correction", "#137c8a", "D")):
            values = [np.mean([r["scores"][name]["rmse_eps"] for r in rows if r["arm"] == arm]) for arm in arms]
            ax.plot(np.arange(len(arms)), values, color=color, marker=marker, markersize=4, label=label)
        ax.set_title(title, fontsize=11)
        ax.set_xticks(np.arange(len(arms)), arms_display, rotation=45, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Mean material RMSE · 8 nonempty test objects")
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    save_figure(fig, "roi_error")
    plt.close(fig)


def save_figure(fig, name):
    stem = ROOT / "figures" / name
    fig.savefig(stem.with_suffix(".png"), dpi=185)
    fig.savefig(stem.with_suffix(".pdf"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arms", nargs="+", default=ARMS)
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--data", type=Path, default=RESULTS / "acquisition.npz")
    parser.add_argument("--truth-metadata", type=Path)
    parser.add_argument("--truth-data", type=Path)
    parser.add_argument("--grid-sensitivity", action="store_true")
    parser.add_argument("--seed-tag", default="")
    parser.add_argument("--dev-only", action="store_true")
    parser.add_argument("--include-sparse", action="store_true")
    parser.add_argument("--only-sparse", action="store_true")
    parser.add_argument("--iterations", type=int, default=6)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--check-operators", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.check_operators:
        check_operators()
    elif args.plot:
        make_figures()
    else:
        compare(args)
