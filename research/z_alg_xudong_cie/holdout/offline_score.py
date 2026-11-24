"""Score complete sealed objects without exposing partial held-out results."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys
from time import sleep

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
if not __package__:
    sys.path.insert(0, str(HERE.parents[2]))
from research.z_alg_xudong_cie.holdout.data import (
    DEFAULT_ROOT, OBJECT_IDS, elapsed, planned_cases, require_freeze, timers, write_json,
)
from research.z_alg_xudong_cie.holdout.scoring import score as score_synthetic

sys.path.insert(0, str(HERE.parent / "chambers"))
from training import CASES as CHAMBERS
from score import score_material as score_chamber

BENCH = HERE.parent / "benchmark"
CHECKPOINT = DEFAULT_ROOT.parent / "sealed_offline_scoring"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def plan_groups(base_record, portfolio_record):
    base, portfolio = read_json(base_record), read_json(portfolio_record)
    methods = tuple(dict.fromkeys([*base["methods"], *portfolio["methods"]]))
    methods_synthetic = (*methods, base["extra_synthetic_method"])
    plan = {}
    for oid in OBJECT_IDS:
        observations = tuple(row["case"] for row in planned_cases() if row["object_id"] == oid)
        plan[("synthetic", oid)] = tuple((name, arm) for name in observations for arm in methods_synthetic)
    for name, case in CHAMBERS.items():
        if case.phase == "confirmation":
            plan[("chamber", name)] = tuple((name, arm) for arm in methods)
    return plan


def ready_inputs(directory, required):
    targets = [(directory / f"{case}__{method}.npz", case, method) for case, method in required]
    if not all(path.is_file() and path.with_suffix(".json").is_file() for path, _, _ in targets):
        return None
    ready = []
    for path, case, method in targets:
        receipt = read_json(path.with_suffix(".json"))
        if (receipt["observation"] != case or receipt["method"] != method or
                not receipt.get("sealed_at") or not receipt.get("public_only")):
            raise ValueError(f"Online receipt failed the public-data seal: {path.with_suffix('.json')}")
        ready.append((path, receipt))
    return ready


def comparison_to_online(score, receipt, category):
    if category == "synthetic":
        independent = score["native_fitting"].get(receipt["observation"])
    else:
        independent = score.get("train")
    prior = receipt.get("relative_native_training_residual")
    return dict(online_relative_residual=prior,
                independent_relative_residual=independent["relative_to_measured"] if independent else None,
                absolute_disagreement=abs(prior-independent["relative_to_measured"])
                if prior is not None and independent else None,
                original_sealed_at=receipt["sealed_at"], online_wall_seconds=receipt["wall_seconds"],
                online_cpu_seconds=receipt["cpu_seconds"])


def score_group(group, inputs, freeze, data_root):
    started = timers()
    kind, oid = group
    rows = []
    cache = {}
    for path, prior in inputs:
        before = timers()
        name, method = prior["observation"], prior["method"]
        meta = dict(case=name, arm=method, path=str(path.resolve()), array="m")
        try:
            with np.load(path, allow_pickle=False) as saved:
                m = saved["m"]
            cached = next(((answer, first) for material, answer, first in cache.get(name, [])
                           if np.array_equal(material, m)), None)
            if cached is not None:
                original, first = cached
                value = deepcopy(original)
                value["reused_identical_map_from_arm"] = first
                value["cached_original_score_seconds"] = original["wall_seconds"]
                value["native_forward"] = dict(value["native_forward"], physical_work_reused=True)
                value.update(elapsed(before))
            else:
                value = (score_synthetic(name, m, freeze_record=freeze, root=data_root)
                         if kind == "synthetic" else score_chamber(name, m, freeze_record=freeze))
                value["reused_identical_map_from_arm"] = None
                cache.setdefault(name, []).append((m, value, method))
            if kind == "synthetic":
                status = "scored" if value["native_forward"]["status"] == "scored" else "material_only"
            else:
                status = "scored" if value.get("held_transmit") else "material_only"
            rows.append(dict(status=status, source=meta, source_resolved_path=str(path.resolve()),
                             **value, online_crosscheck=comparison_to_online(value, prior, kind)))
        except (ArithmeticError, RuntimeError, OSError, ValueError, KeyError) as exc:
            rows.append(dict(status="unresolved", source=meta, reason=str(exc), **elapsed(before)))
    return dict(object_id=oid, category=kind, rows=rows, expected_rows=len(inputs), **elapsed(started))


def assemble(groups, completed, output):
    paths = dict(synthetic=output/"synthetic_scores.json", chamber=output/"chamber_scores.json")
    for category, target in paths.items():
        names = [group for group in groups if group[0] == category]
        if all(group in completed for group in names):
            samples = [completed[group] for group in names]
            rows = [row for item in samples for row in item["rows"]]
            write_json(target, dict(
                scope="sealed native blind scoring and actual offline charges; duplicate maps reused only within a view",
                rows=rows, independent_materials=len(names), complete_online_pairs=len(rows),
                aggregate_scoring_wall_seconds=sum(report["wall_seconds"] for report in samples),
                aggregate_scoring_cpu_seconds=sum(report["cpu_seconds"] for report in samples),
                expected=len(rows), status_counts={tag: sum(row["status"] == tag for row in rows)
                         for tag in ("scored", "material_only", "unresolved")},
            ))


def run(base_freeze, portfolio_freeze, predictions, data_root, output, checkpoint, wait):
    freeze = require_freeze(portfolio_freeze)
    groups = plan_groups(require_freeze(base_freeze), freeze)
    completed = {}
    checkpoint.mkdir(parents=True, exist_ok=True)
    while len(completed) < len(groups):
        for group, plan in groups.items():
            if group in completed:
                continue
            path = checkpoint / f"{group[0]}_{group[1]}.json"
            if path.is_file():
                report = read_json(path)
                if {(row["source"]["case"], row["source"]["arm"]) for row in report["rows"]} != set(plan):
                    raise ValueError(f"Saved offline checkpoint does not match the frozen arm list: {path}")
                completed[group] = report
                continue
            pending = ready_inputs(predictions, plan)
            if pending is None:
                continue
            result = score_group(group, pending, freeze, data_root)
            write_json(path, result)
            completed[group] = result
            print(group[0], group[1], len(result["rows"]), "sealed maps independently assessed", flush=True)
        if len(completed) == len(groups):
            assemble(groups, completed, output)
        if len(completed) < len(groups):
            if not wait:
                break
            sleep(15)
    return len(completed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-freeze", type=Path, default=BENCH / "FREEZE.json")
    parser.add_argument("--portfolio-freeze", type=Path, default=BENCH / "PORTFOLIO_FREEZE.json")
    parser.add_argument("--predictions", type=Path, default=BENCH / "predictions")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=Path, default=BENCH)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT)
    parser.add_argument("--wait", action="store_true")
    args = parser.parse_args()
    with threadpool_limits(1):
        count = run(args.base_freeze, args.portfolio_freeze, args.predictions, args.root,
                    args.output, args.checkpoint, args.wait)
    print(count, "complete independent objects", flush=True)
