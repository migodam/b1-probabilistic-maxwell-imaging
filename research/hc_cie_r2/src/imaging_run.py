"""Run the frozen finite-jump experiment without importing any scoring code."""
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time
import zlib

import numpy as np
from scipy.linalg import lu_solve
from scipy.optimize import minimize
from scipy.stats import chi2
from threadpoolctl import threadpool_info


LAB = Path(__file__).resolve().parents[1]


def data_root():
    override = os.environ.get("HC_CIE_R2_INPUT_ROOT")
    if override:
        return Path(override).expanduser().resolve()
    markers = ("research/delegated/max_struct_r1", "research/research_max_struct_r13",
               "r13_worktree/research/research_max_struct_r13")
    return next((folder for folder in LAB.parents
                 if any((folder / marker).is_dir() for marker in markers)), LAB.parents[1])


def history_tree(root, version):
    subtree = f"research/research_max_struct_{version}"
    nested = Path("research/delegated/max_struct_r1") / f"{version}_worktree"
    alternatives = (nested / subtree, Path(f"{version}_worktree") / subtree, Path(subtree))
    required = "data_public/fresh48" if version == "r13" else "configs/primary_manifest.json"
    for path in alternatives:
        if (root / path / required).exists():
            return (root / path).parents[1]
    return root


REPO = data_root()
R13_TREE = history_tree(REPO, "r13")
R6_TREE = history_tree(REPO, "r6")
R13 = R13_TREE / "research/research_max_struct_r13"
R6 = R6_TREE / "research/research_max_struct_r6"


def recorded_input_path(recorded):
    original = REPO / recorded
    if original.exists():
        return original
    parts = Path(recorded).parts
    for marker, current_tree in (("r13_worktree", R13_TREE), ("r6_worktree", R6_TREE)):
        if marker in parts:
            return current_tree.joinpath(*parts[parts.index(marker) + 1:])
    return original


def json_write(path, data):
    def convert(value):
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Cannot serialize {type(value).__name__}")

    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, allow_nan=False, default=convert) + "\n")
    temp.replace(path)


def public_input(case, endpoint):
    keys = ("G", "E0", "E", "S", "Y", "grid_shape", "delta", "sigma", "m_init", "lower", "upper")
    if case.startswith("r13_fresh_"):
        path = R13 / "data_public/fresh48" / f"{case}.npz"
        if endpoint == "a_plus_arc_r7":
            parent = R13 / "supplement_aplus48/stack_run/RESULTS/aplus_stack112"
        else:
            parent = R13 / "RESULTS/confirmation_core48"
        paths = [path]
        snapshot = parent / f"{case}__{endpoint}" / "primary.npz"
    else:
        manifest = json.loads((R6 / "configs/primary_manifest.json").read_text())
        task = next(t for t in manifest["tasks"] if t["case_id"] == case and t["arm"] == endpoint)
        paths = [R6_TREE / p for p in task["public_files"]]
        snapshot = R6_TREE / task["output_dir"] / "prediction.npz"
    public = {}
    for path in paths:
        with np.load(path, allow_pickle=False) as raw:
            public.update({name: raw[name] for name in keys if name in raw})
    if "E0" not in public:
        public["E0"] = public["E"]
    with np.load(snapshot, allow_pickle=False) as raw:
        start, saved_y = raw["m"], raw["Y_prediction"]
    return public, start, saved_y, [str(p.relative_to(REPO)) for p in [*paths, snapshot]]


def tiles(shape, width):
    indices = np.arange(np.prod(shape)).reshape(shape)
    return [indices[y:y + width, x:x + width].ravel()
            for y in range(0, shape[0], width) for x in range(0, shape[1], width)]


def make_actions(model, start, state, lower, upper, config, mode, ledger):
    # L = I - D(start) G is the retained full-current factorization.
    # T = S L^-1 is the dressed receiver; K = G L^-1 maps impressed
    # current changes to total-field changes. No fixed global Gamma is formed.
    with model._physics.operation("BLOCK_RECEIVER", rhs=model.nr, receiver_transpose_RHS=model.nr):
        transfer = lu_solve(state["lu"], model.S.T, trans=1, check_finite=False).T
    need_k = mode in {"single", "pair_secant", "pair_exact"}
    internal = None
    if need_k:
        with model._physics.operation("BLOCK_PROPAGATION", rhs=model.n, block_propagation_RHS=model.n):
            internal = lu_solve(state["lu"], model.G.T, trans=1, check_finite=False).T

    ref_prior = model.reg.evaluate(start)[0]
    actions = []
    for tile, idx in enumerate(tiles(model.shape, config["block_side"])):
        basis = state["Q"][idx, :]
        for alpha in config["amplitudes"]:
            for sign in (-1, 1):
                value = np.clip(start[idx] + sign * alpha, lower, upper)
                change = value - start[idx]
                if not np.any(change):
                    continue
                material = start.copy()
                material[idx] = value
                rhs = change[:, None] * basis
                dressed = rhs
                if need_k:
                    local = np.eye(len(idx)) - change[:, None] * internal[np.ix_(idx, idx)]
                    dressed = np.linalg.solve(local, rhs)
                    ledger["small_LU"] += 1
                    ledger["small_LU_dimension_sum"] += len(idx)
                    ledger["small_RHS"] += model.ns
                actions.append({"tile": tile, "idx": idx, "material_value": value,
                                "alpha": alpha, "sign": sign, "change": change,
                                "signature": (transfer[:, idx] @ dressed).ravel(),
                                "prior_diff": float(model.reg.evaluate(material)[0] - ref_prior)})
    ledger["action_count"] += len(actions)
    return actions, transfer, internal, ref_prior


def pick_actions(model, start, state, Y, lower, upper, mode, config, rng, ledger):
    actions, transfer, internal, ref_prior = make_actions(
        model, start, state, lower, upper, config, mode, ledger)
    residual = (state["Y"] - Y).ravel()
    scale = np.vdot(Y, Y).real
    base = .5 * np.vdot(residual, residual).real / scale + ref_prior
    signatures = np.asarray([a["signature"] for a in actions])
    residual_linear = (signatures.conj() @ residual).real / scale
    stand_alone = np.sum(np.abs(signatures)**2, axis=1) / (2 * scale)
    predicted = residual_linear + stand_alone + np.asarray([a["prior_diff"] for a in actions])
    if mode == "single":
        pairs = [(i,) for i in range(len(actions))]
        base_scores = base + predicted
    else:
        pairs = [(i, j) for i, j in combinations(range(len(actions)), 2)
                 if actions[i]["tile"] != actions[j]["tile"]]
        i, j = np.asarray(pairs).T
        gram = (signatures @ signatures.conj().T).real / scale
        base_scores = base + predicted[i] + predicted[j] + gram[i, j]
    ledger["proposal_deck_count"] += len(pairs)
    n = config["proposal_audits"]
    if mode == "pair_random":
        picked = rng.permutation(len(pairs))[:n]
    elif mode == "pair_exact":
        shortlist = np.argsort(base_scores, kind="stable")[:config["exact_shortlist"]]
        coupled = []
        for index in shortlist:
            first, second = [actions[k] for k in pairs[index]]
            li, lj = first["idx"], second["idx"]
            support = np.r_[li, lj]
            change = np.r_[first["change"], second["change"]]
            # Sum all scattering paths between both changed material blocks.
            system = np.eye(len(support)) - change[:, None] * internal[np.ix_(support, support)]
            rhs = change[:, None] * state["Q"][support]
            current = np.linalg.solve(system, rhs)
            ledger["small_LU"] += 1
            ledger["small_LU_dimension_sum"] += len(support)
            ledger["small_RHS"] += model.ns
            row = residual + (transfer[:, support] @ current).ravel()
            material = start.copy()
            material[li] = first["material_value"]
            material[lj] = second["material_value"]
            score = .5 * np.vdot(row, row).real / scale + model.reg.evaluate(material)[0]
            coupled.append((float(score), int(index)))
        ledger["joint_exact_scored_count"] += len(shortlist)
        coupled.sort()
        picked = [index for _, index in coupled[:n]]
        for score, index in coupled:
            base_scores[index] = score
    else:
        picked = np.argsort(base_scores, kind="stable")[:n]
    nominations = []
    for index in picked:
        material = start.copy()
        descriptions = []
        for k in pairs[index]:
            action = actions[k]
            idx = action["idx"]
            material[idx] = action["material_value"]
            descriptions.append({"tile": action["tile"], "alpha": action["alpha"],
                                 "sign": action["sign"], "changed": int(np.count_nonzero(action["change"]))})
        nominations.append((material, float(base_scores[index]), descriptions))
    return nominations


def run_one(task):
    case, endpoint, mode, config = task
    wall0, cpu0 = time.perf_counter(), time.process_time()
    out = Path(config["output"]) / f"{case}__{endpoint}__{mode}"
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(R13_TREE))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from collections import Counter

    path_config = {k: v for k, v in config.items() if k != "output"}
    public, start, saved_y, input_paths = public_input(case, endpoint)
    model = MaxwellModel(public)
    lower, upper = float(public.get("lower", .02)), float(public.get("upper", 20.))
    Y = public["Y"]
    delta = float(public.get("delta", float(public.get("sigma", 0.)) *
                    np.sqrt(chi2.ppf(.99, 2 * Y.size) / 2)))
    unit = np.linalg.norm(Y)
    cap = config["budget_cpu_seconds"]
    actions = []
    materials = []
    fields = []
    ledger = Counter()
    rng = np.random.default_rng(config["seed"] + zlib.crc32(f"{case}:{endpoint}".encode()))
    best = None
    termination = None
    correction_iterations = None

    class StopBudget(Exception):
        """The current evaluation is outside the common paid prefix."""

    def check():
        if time.process_time() - cpu0 >= cap:
            raise StopBudget

    def audit(material, state, kind, forecast=None, blocks=None):
        nonlocal best
        cpu = time.process_time() - cpu0
        wall = time.perf_counter() - wall0
        norm = float(np.linalg.norm(state["Y"] - Y))
        data = .5 * (norm / unit)**2
        prior = float(model.reg.evaluate(material)[0])
        if np.min(material) < lower - 1e-12 or np.max(material) > upper + 1e-12:
            raise ValueError("Illegal trial: material outside shared box")
        phi = data + prior
        index = len(actions)
        row = dict(i=index, stage=kind, cpu=cpu, wall=wall, D=data, R=prior, Phi=phi,
                   relative_residual=float(norm / unit), residual_over_delta=norm / delta if delta > 0 else None,
                   state_relative_residual=state["state_relative_residual"], paid=bool(cpu <= cap),
                   forecast=forecast, blocks=blocks, cost=dict(model.cost), proposal_cost=dict(ledger))
        actions.append(row)
        materials.append(material.copy())
        fields.append(state["Y"].copy())
        if row["paid"] and (best is None or phi < actions[best]["Phi"]):
            best = index

    state = model.state(start)
    saved_error = float(np.linalg.norm(saved_y - state["Y"]) / max(np.linalg.norm(saved_y), 1e-100))
    if saved_error > 1e-7:
        raise ValueError(f"Historical endpoint port disagreement {saved_error:.8g}")
    audit(start, state, "endpoint")
    before_polish = best
    proposal_cpu = 0.
    chosen_jump = None

    try:
        if mode not in {"continue", "restart"}:
            proposal_started = time.process_time()
            proposals = pick_actions(model, start, state, Y, lower, upper, mode, config, rng, ledger)
            for material, forecast, blocks in proposals:
                check()
                audit(material, model.state(material), "jump", forecast, blocks)
            proposal_cpu = time.process_time() - proposal_started
            before_polish = best
            paid_jumps = [r for r in actions if r["stage"] == "jump" and r["paid"]]
            if paid_jumps:
                chosen_jump = min(paid_jumps, key=lambda r: r["Phi"])["i"]
        start_index = best
        if config.get("allow_uphill_seed") and chosen_jump is not None:
            start_index = chosen_jump
        polish_start = materials[start_index]
        if mode == "restart":
            polish_start = np.asarray(public.get("m_init", np.full(model.n, .1)), dtype=float)

        def objective(material):
            check()
            phi, grad = model.data_value_gradient(material, Y)
            audit(material, model.state(material), "direct")
            check()
            return 1e6 * phi, 1e6 * grad

        res = minimize(objective, polish_start, method="L-BFGS-B", jac=True,
                       bounds=[(lower, upper)] * model.n,
                       options=dict(maxiter=6000, maxfun=24000, ftol=1e-12, gtol=.01, maxls=30))
        termination = str(res.message)
        correction_iterations = int(res.nit)
    except StopBudget:
        termination = "PROCESS_CPU_LIMIT"

    predicted = np.asarray(materials)
    ports = np.asarray(fields)
    np.savez_compressed(out / "PREDICTIONS.npz", material=predicted,
                        selected_port=ports[best], selected_index=np.asarray(best),
                        endpoint_port=ports[0])
    receipt = dict(case=case, endpoint=endpoint, arm=mode, status="PREDICTION_COMPLETE",
                   endpoint_relative_port_disagreement=saved_error, config=path_config,
                   input_paths=input_paths, selected=best, selected_before_polish=before_polish,
                   chosen_jump=chosen_jump, proposal_cpu=proposal_cpu,
                   total_cpu=time.process_time() - cpu0, total_wall=time.perf_counter() - wall0,
                   actual_cpu_seconds_limit=cap, termination=termination,
                   correction_iterations=correction_iterations, counts=dict(model.cost),
                   proposal_counts=dict(ledger), threadpools=threadpool_info(),
                   without_truth=True, geometry_shape=list(model.shape), physical_counts=Y.shape,
                   lower=lower, upper=upper, chain=actions)
    json_write(out / "TRACE.json", receipt)
    return f"{case} {endpoint} {mode} {len(actions)} {actions[0]['Phi']:.6g} {actions[best]['Phi']:.6g}"


def tasks_for(config, endpoint_index=None, cases=None):
    tasks = []
    for key in ("r13", "r6"):
        available_cases = config[f"{key}_cases"]
        endpoints = config[f"{key}_endpoints"]
        if endpoint_index is not None:
            endpoints = [endpoints[endpoint_index]]
        for case in available_cases:
            if cases and case not in cases:
                continue
            for endpoint in endpoints:
                # Rotate arm schedule across start materials and across objects.
                arms = config["arms"]
                offset = zlib.crc32(f"{case}:{endpoint}".encode()) % len(arms)
                order = arms[offset:] + arms[:offset]
                tasks += [(case, endpoint, arm, config) for arm in order]
    return tasks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    parser.add_argument("--endpoint-index", type=int, choices=(0, 1, 2, 3))
    parser.add_argument("--cases", nargs="+")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    cfg["output"] = str(args.out.resolve())
    args.out.mkdir(parents=True, exist_ok=False)
    tasks = tasks_for(cfg, args.endpoint_index, args.cases)
    started = time.perf_counter()
    with ProcessPoolExecutor(args.workers) as pool:
        for result in pool.map(run_one, tasks):
            print(result, flush=True)
    manifest = []
    for case, endpoint, mode, _ in tasks:
        ident = f"{case}__{endpoint}__{mode}"
        receipt = json.loads((args.out / ident / "TRACE.json").read_text())
        assert receipt["status"] == "PREDICTION_COMPLETE"
        manifest.append({"case": case, "endpoint": endpoint, "arm": mode,
                         "directory": ident, "selected": receipt["selected"],
                         "candidate_count": len(receipt["chain"])})
    json_write(args.out / "PREDICTION_SEAL.json",
               {"id": cfg["id"], "without_truth": True, "workers": args.workers,
                "batch_wall": time.perf_counter() - started, "cpu_info": platform.uname()._asdict(),
                "n": len(manifest), "jobs": manifest})


if __name__ == "__main__":
    main()
