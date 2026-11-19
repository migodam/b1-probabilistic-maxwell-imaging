"""Reconstruct from the sealed old-only actions and their shared measured reads."""
import argparse
from dataclasses import dataclass
import json
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from data14 import (ACTIONS, HALF_SPLIT_FIRST, Fit, FrequencyTM, MENU, ROOT, encode_complex,
                    joint_fit, load_cases, measured_reads)
from measure_select import starts


class FrequencyReuse(FrequencyTM):
    """One Maxwell factorization per distinct frequency at a material iterate."""
    def __init__(self, geometry):
        super().__init__(geometry)
        self.material = None
        self.cache = {}

    def response(self, m, delta=0., hessian=False, derivative=False):
        if not np.array_equal(m, self.material):
            self.material = np.asarray(m).copy()
            self.cache.clear()
        key = (delta, hessian, derivative)
        if key not in self.cache:
            self.cache[key] = super().response(m, delta, hessian, derivative)
        return self.cache[key]


def stage_indices(action, mode):
    deltas = (0., *action.deltas)
    order = tuple(sorted(range(len(deltas)), key=lambda i: deltas[i]))
    if mode == "cumulative3":
        return [tuple(i for i in order if deltas[i] <= delta) for delta in sorted(set(deltas))]
    if mode.startswith("fbch"):
        middle = (len(order) + 1) // 2
        return [order[:middle], order]
    return [tuple(range(len(deltas)))]


def make_stage(case, model, action, new, indices, wavelength_power):
    deltas = np.asarray((0., *action.deltas))[list(indices)]
    energy = np.asarray((1., *action.energies))[list(indices)]
    data = np.asarray((case.old, *new))[list(indices)]
    # WW exponent applies to the squared cost. The statistical whitening is still sqrt(energy).
    weights = np.sqrt(energy / (1. + deltas) ** wavelength_power)
    return Fit(model, deltas, data, weights=weights,
               gain_calibration=((case.sigma / case.gain_sd) ** 2, 1.))


def solve(case, action, new, old_starts, mode, budget):
    clock = perf_counter()
    model = FrequencyReuse(case.geometry)
    fit = joint_fit(case, model, action, new)
    seeds = starts(3 if mode.endswith("8") else 1)
    if mode == "single":
        seeds = seeds[:1]
    elif mode.startswith("root"):
        old = sorted(old_starts, key=lambda m: np.linalg.norm(fit.evaluate(m)[0]))
        seeds = old + seeds
    if mode.endswith("8"):
        seeds = seeds[:8]
    else:
        seeds = seeds[:3]
    plan = stage_indices(action, mode)
    wavelength_power = int(mode[-1]) if mode.startswith("fbch") else 0
    prediction = None
    final_loss = np.inf
    final_gain = None
    ending_losses = []
    stages_run = 0
    per_start = (budget - model.factorizations) // len(seeds)
    deltas = (0., *action.deltas)
    for seed in seeds:
        estimate = seed
        # Every method gets the same factorization cap. Reserve a full final score per start.
        remaining = per_start - len(set(deltas))
        for k, ix in enumerate(plan):
            n_unique = len({deltas[i] for i in ix})
            max_evals = max(1, remaining // (len(plan) - k) // n_unique)
            stage = make_stage(case, model, action, new, ix, wavelength_power)
            start_count = model.factorizations
            estimate, _ = stage.solve(estimate, max_evals)
            remaining -= model.factorizations - start_count
            stages_run += 1
        residual, _, _, _, gain = fit.evaluate(estimate)
        loss = float(np.linalg.norm(residual / case.sigma) ** 2)
        ending_losses.append(loss)
        if loss < final_loss:
            prediction, final_loss, final_gain = estimate, loss, gain
    if model.factorizations > budget:
        raise RuntimeError(f"Factorization cap exceeded for {case.case_id}, {action.name}, {mode}")
    return dict(prediction=prediction.tolist(), gain=encode_complex(final_gain), loss=final_loss,
                starts=len(seeds), optimization_stages=stages_run,
                start_losses=ending_losses, fit_factors=model.factorizations,
                fit_rhs=6 * model.factorizations,
                fit_wall_s=perf_counter() - clock)


@dataclass(frozen=True)
class Comparison:
    label: str
    action_label: str
    mode: str


def comparisons(actions, external_label):
    fixed = ("old_only", "repeat", "split_0.02", "uniform_0.02", "split_0.04", "split_0.08", "uniform_0.08")
    adaptive = ("fisher", "candidate_gap", "random", "random_at_cap", "geometry")
    if external_label:
        adaptive += (external_label,)
    selected = [Comparison(label, actions[label], mode) for label in (*fixed, *adaptive)
                for mode in ("single", "ordinary3", "root3")]
    selected.extend(Comparison(label, actions[label], mode)
                    for label in ("split_0.02", "uniform_0.02", "uniform_0.08")
                    for mode in ("cumulative3", "fbch0", "fbch1", "fbch2"))
    selected.extend(Comparison(label, actions[label], mode)
                    for label in ("split_0.02", "uniform_0.02", "split_0.08", *adaptive)
                    for mode in ("ordinary8", "root8"))
    return selected


def design_costs(choice, label, mode, external):
    costs = choice["costs"]
    root_cost = mode.startswith("root") or label in ("fisher", "candidate_gap")
    count = costs["old_correction_factors"] if root_cost else 0
    wall = costs["root_scan_s"] + costs["old_correction_s"] if root_cost else 0.
    for key in ("fisher", "candidate_gap", "geometry"):
        if label == key:
            count += costs[f"{key}_factors"]
            wall += costs[f"{key}_s"]
    rhs = 6 * count
    if root_cost:
        count += costs["root_validation_factors"]
        rhs += costs["root_validation_rhs"]
    if external:
        count += external.get("factors", 0)
        rhs += external.get("rhs", 0)
        wall += external.get("wall_s", 0.)
    return dict(design_factors=count, design_rhs=rhs, design_wall_s=wall,
                polynomial_paths=costs["polynomial_paths"] if root_cost else 0,
                polynomial_path_jacobians=costs["path_jacobian_evaluations"] if root_cost else 0,
                design_fallback=(choice["fallback_candidate_gap"] if label == "candidate_gap"
                                 else external.get("fallback")),
                design_certified=external.get("certified"), design_diameter=external.get("D_upper"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external", type=lambda name: ROOT / name,
                        help="Case-to-action map or case-to-{action,factors,wall_s,fallback} map.")
    parser.add_argument("--budget", default=780, type=int)
    parser.add_argument("--external-only", action="store_true")
    parser.add_argument("--external-label", default="region")
    parser.add_argument("--modes", nargs="+", help="Limit the requested reconstruction modes.")
    parser.add_argument("--half-only", action="store_true")
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    external = json.loads(args.external.read_text()) if args.external else {}
    choices = json.loads((ROOT / "results/measure_choices.json").read_text())
    previous = {}
    if args.external_only or args.half_only:
        for filename in ("measure_predictions.json", "measure_predictions_half.json"):
            cached = ROOT / "results" / filename
            if cached.is_file():
                base = json.loads(cached.read_text())
                if base["budget"] == args.budget:
                    previous.update({(r["case_id"], r["action"], r["mode"]): r for r in base["rows"]})
    fit_keys = ("prediction", "gain", "loss", "starts", "optimization_stages",
                "start_losses", "fit_factors", "fit_rhs", "fit_wall_s")
    records = []
    with np.load(ROOT / "results/measure_new_reads.npz") as observations, threadpool_limits(limits=1):
        for case, choice in zip(load_cases(), choices["cases"], strict=True):
            assert case.case_id == choice["case_id"]
            item = external.get(case.case_id, {})
            if isinstance(item, str):
                item = {"action": item}
            selected = dict(choice["actions"])
            if item:
                selected[args.external_label] = item["action"]
            studies = comparisons(selected, args.external_label if item else None)
            if args.external_only:
                studies = [v for v in studies if v.label == args.external_label]
            elif args.half_only:
                studies = [Comparison(label, action, mode)
                           for label, action in (("half_first", HALF_SPLIT_FIRST.name),
                                                 ("full_split", "split_0.02"))
                           for mode in ("root3", "ordinary8")]
            else:
                studies.extend(Comparison("sweep", a.name, "root3") for a in MENU)
            if args.modes:
                studies = [v for v in studies if v.mode in args.modes]
            recovered = {}
            for entry in studies:
                action = ACTIONS[entry.action_label]
                unique = (action.name, entry.mode)
                if unique not in recovered:
                    new = measured_reads(observations, case.case_id, action)
                    cached_fit = previous.get((case.case_id, *unique))
                    recovered[unique] = ({key: cached_fit[key] for key in fit_keys} if cached_fit
                                         else solve(case, action, new, np.asarray(choice["old_starts"]),
                                                    entry.mode, args.budget))
                result = recovered[unique]
                records.append(dict(case_id=case.case_id, cohort=case.cohort, label=entry.label,
                                    action=action.name, mode=entry.mode, **result,
                                    **design_costs(choice, entry.label, entry.mode,
                                                   item if entry.label == args.external_label else {})))
            print(case.case_id, len(recovered), sum(r["fit_factors"] for r in recovered.values()),
                  flush=True)
    result = dict(budget=args.budget, rows=records)
    path = ROOT / "results" / f"measure_predictions{args.tag}.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
