"""Select controls from old data, then write the choices before acquisition."""
import json
from itertools import combinations
from time import perf_counter

import numpy as np
from scipy.stats import qmc

from data14 import ACTIONS, BOX, Fit, FrequencyTM, MENU, R10, ROOT, load_cases

import sys
sys.path.insert(0, str(R10 / "src"))
from homotopy_control import run as find_roots


def starts(power=1):
    samples = qmc.Sobol(6, scramble=True, seed=14_20261119).random_base2(power)
    return [np.full(6, 12.), *(BOX[0] + (BOX[1] - BOX[0]) * samples)]


def action_record(action):
    return dict(name=action.name, deltas=action.deltas, energies=action.energies,
                energy=action.energy, exposures=action.exposures,
                bandwidth=action.bandwidth, max_offset=action.max_offset)


def profiled_information(y, j, lam):
    cross = np.column_stack((j.conj().T @ y, 1j * (j.conj().T @ y))).real
    gram = (j.conj().T @ j).real - cross @ cross.T / (np.vdot(y, y).real + lam)
    return max(float(np.linalg.eigvalsh(gram)[0]), 0.)


def select(case, number):
    work = perf_counter()
    scan = find_roots(case.geometry, case.old)
    zeros = np.asarray(scan["validated_roots"]).reshape(-1, 6)
    zeros = sorted(zeros, key=lambda m: np.linalg.norm(m - np.full(6, 12.)))
    lam = (case.sigma / case.gain_sd) ** 2
    ordinary = starts()
    model = FrequencyTM(case.geometry)
    root_time = perf_counter() - work
    if not zeros:
        fit = Fit(model, [0.], [case.old], gain_calibration=(lam, 1.))
        endpoint = [fit.solve(m, 130)[0] for m in ordinary]
        zeros_for_fit = [min(endpoint, key=lambda m: np.linalg.norm(fit.evaluate(m)[0]))]
    else:
        zeros_for_fit = zeros
    roots_cost = model.factorizations
    anchor = zeros_for_fit[0]
    deltas = sorted({0., *(delta for action in MENU for delta in action.deltas)})
    anchor_time = perf_counter() - work - root_time

    def grid(materials):
        return [{delta: model.response(m, delta) for delta in deltas} for m in materials]

    def blocks(grid_row, action, component):
        return np.concatenate([np.sqrt(e) * grid_row[d][component]
                               for e, d in zip((1., *action.energies), (0., *action.deltas))])

    def e_opt(raw):
        return {a.name: float(min(profiled_information(blocks(elem, a, 0), blocks(elem, a, 1), lam)
                                  for elem in raw)) for a in MENU}

    def choose(scores):
        return max(MENU, key=lambda a: (scores[a.name], -a.bandwidth, -a.exposures)).name

    def mismatch(u, v):
        # Direct orthogonal residual avoids subtracting nearly equal power levels.
        g = np.vdot(v, u) / np.vdot(v, v)
        return float(np.linalg.norm(u - g * v) ** 2)

    t0 = perf_counter()
    b0 = model.factorizations
    fisher_grid = grid([anchor])
    fisher = e_opt(fisher_grid)
    choose_fisher = choose(fisher)
    fisher_time = perf_counter() - t0
    fisher_cost = model.factorizations - b0

    t0 = perf_counter()
    b0 = model.factorizations
    layout = e_opt(grid([np.full(6, x) for x in (8., 16., 24.)]))
    choose_layout = choose(layout)
    geometry_cost = model.factorizations - b0
    geometry_time = perf_counter() - t0

    separation = {}
    if len(zeros) >= 2:
        t0 = perf_counter()
        b0 = model.factorizations
        points = [*fisher_grid, *grid(zeros[1:])]
        old = [np.r_[p[0.][0], np.sqrt(lam)] for p in points]
        directed = [(i, j, mismatch(u, v)) for i, u in enumerate(old)
                    for j, v in enumerate(old) if i != j]
        for action in MENU:
            pred = [np.r_[blocks(row, action, 0), np.sqrt(lam)] for row in points]
            separation[action.name] = float(np.sqrt(max(0., min(
                mismatch(pred[i], pred[j]) - base for i, j, base in directed))))
        choose_separation = choose(separation)
        separation_cost = model.factorizations - b0 + fisher_cost
        separation_time = perf_counter() - t0 + fisher_time
    else:
        choose_separation, separation_cost, separation_time = choose_fisher, fisher_cost, fisher_time

    random = np.random.default_rng(200309 + number).choice(
        [a.name for a in MENU if a.name != "repeat"])
    random_cap = np.random.default_rng(802471 + number).choice(
        [a.name for a in MENU if a.max_offset == .02])
    choices = dict(fisher=choose_fisher, candidate_gap=choose_separation,
                   geometry=choose_layout, random=str(random), random_at_cap=str(random_cap))
    fixed = {label: label for label in ("old_only", "repeat", "split_0.02", "uniform_0.02",
                                        "split_0.04", "split_0.08", "uniform_0.08")}
    choices.update(fixed)
    weights = dict(fisher=fisher, geometry=layout, candidate_gap=separation)
    classes = scan["path_classifications"]
    accepted = classes.get("accepted_distinct", 0) + classes.get("accepted_duplicate", 0)
    full_checks = accepted + classes.get("finite_failed_full_VIE", 0)
    costs = dict(root_scan_s=root_time, old_correction_s=anchor_time,
                 root_validation_factors=full_checks + accepted,
                 root_validation_rhs=3 * full_checks + 4 * accepted,
                 old_correction_factors=roots_cost, fisher_factors=fisher_cost,
                 fisher_s=fisher_time, candidate_gap_factors=separation_cost,
                 candidate_gap_s=separation_time, geometry_factors=geometry_cost,
                 geometry_s=geometry_time, polynomial_paths=scan["paths_attempted"],
                 path_jacobian_evaluations=scan["path_jacobian_evaluations"])
    return dict(case_id=case.case_id, cohort=case.cohort,
                exact_central_roots=[m.tolist() for m in zeros],
                old_starts=[m.tolist() for m in zeros_for_fit], actions=choices,
                fallback_candidate_gap=len(zeros) < 2, scores=weights, costs=costs)


def main():
    records = []
    for i, case in enumerate(load_cases()):
        result = select(case, i)
        records.append(result)
        print(case.case_id, len(result["exact_central_roots"]), result["actions"]["fisher"],
              result["actions"]["candidate_gap"], flush=True)
    path = ROOT / "results/measure_choices.json"
    path.write_text(json.dumps(dict(menu={name: action_record(a) for name, a in ACTIONS.items()},
                                    cases=records), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
