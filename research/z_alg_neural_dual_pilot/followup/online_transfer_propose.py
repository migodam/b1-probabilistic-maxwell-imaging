"""Generate the frozen ordinary one-second coarse proposals on the public Transfer16 views."""

import argparse
import gzip
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for name in ("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from threadpoolctl import threadpool_limits


PILOT = Path(__file__).resolve().parents[1]
RESULTS = PILOT / "followup/results"
FROZEN = RESULTS / "online_transfer_frozen.json"
DATA = PILOT / "data/transfer16.npz"
CAP, SIDE, ARMS = .75, 8, "gradient1"
NOISES = (1, 2, 10)
sys.path.insert(0, str(PILOT / "src"))
from score_batch import clock, elapsed, shard_pair


def observations():
    with np.load(DATA, allow_pickle=False) as raw:
        if set(raw.files) != {"Y", "radius", "noise_pct", "object_id"}:
            raise ValueError("The canonical observations may contain only the noisy inputs.")
        y, radius, noise, ids = (raw[name] for name in ("Y", "radius", "noise_pct", "object_id"))
    if (
        y.shape != (16, 3, 20, 10) or radius.shape != (16, 3) or ids.shape != (16,)
        or len(set(ids.tolist())) != 16 or tuple(noise) != NOISES
        or not np.isfinite(y).all() or not np.isfinite(radius).all() or not np.all(radius > 0)
    ):
        raise ValueError("Invalid frozen Transfer16 field shapes, IDs, or paired radii.")
    return {str(oid): {
        level: (float(radius[i, j]), y[i, j], (CAP,)) for j, level in enumerate(NOISES)
    } for i, oid in enumerate(ids)}


def identifier(oid, noise):
    return f"{oid}__n{noise}__u{CAP:g}__b{SIDE}__{ARMS}"


def gzip_json(path, value):
    temporary = path.with_suffix(".tmp")
    with gzip.open(temporary, "wt", compresslevel=3) as stream:
        json.dump(value, stream, allow_nan=False, separators=(",", ":"))
    temporary.replace(path)


def row(coarse, oid, noise, radius, y):
    result = coarse.optimize(y, CAP, blocks=SIDE**2, warm_start=None, time_limits=(1.,))
    prefix = result["prefixes"][1.]
    candidate = prefix["best"]
    if candidate is not None and candidate["elapsed_wall_s"] > 1.:
        raise ArithmeticError("Selected evaluation completed beyond the frozen prefix.")
    return {
        "id": identifier(oid, noise), "oid": oid, "noise": noise, "radius": radius,
        "u": CAP, "side": SIDE, "arm": ARMS,
        "mu": candidate["mu"] if candidate else None, "lam": candidate["lam"] if candidate else None,
        "coarse_d2": candidate["d2"] if candidate else None, "nu": candidate["nu"] if candidate else None,
        "proposal_wall": prefix["wall_budget_used_s"], "proposal_cpu": prefix["completed_cpu_s"],
        "paid_wall": result["paid_wall_s"], "paid_cpu": result["paid_cpu_s"],
        "evals": prefix["eval_count"], "overrun_evals": result["over_budget_evals"],
        "initial": result["trace"][0], "trace": result["trace"], "optimizer_stop": result["stop_reason"],
        "last_eligible_wall": prefix["completed_wall_s"], "best_completed_wall": (
            candidate["elapsed_wall_s"] if candidate else None
        ),
    }


def run(shard):
    from coarse import Coarse

    started = clock()
    protocol = json.loads(FROZEN.read_text())
    cases = observations()
    rank, workers = shard
    selected = sorted(cases)[rank::workers]
    path = RESULTS / f"online_transfer_proposals-{rank}-of-{workers}.json.gz"
    if path.is_file():
        with gzip.open(path, "rt") as stream:
            previous = json.load(stream)
        if (
            previous["freeze_utc"] != protocol["freeze_utc"] or previous["shard"] != [rank, workers]
            or [item["id"] for item in previous["rows"]] != [
                identifier(oid, noise) for oid in selected for noise in NOISES
            ]
        ):
            raise ValueError("The existing proposal shard differs from this frozen selection.")
        print(path.name, "sealed proposals reused", flush=True)
        return
    with threadpool_limits(limits=1):
        oracle = Coarse()
        setup = elapsed(started)
        rows = []
        for oid in selected:
            for noise in NOISES:
                radius, y, _ = cases[oid][noise]
                result = row(oracle, oid, noise, radius, y)
                rows.append(result)
                print(result["id"], "eligible", result["evals"], "d2", result["coarse_d2"], flush=True)
    gzip_json(path, {"freeze_utc": protocol["freeze_utc"], "shard": shard, "objects": selected,
                     "setup": setup, "coarse_setup": oracle.setup_cost, "rows": rows, **elapsed(started)})


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--shard", type=shard_pair, default=(0, 1))
    run(cli.parse_args().shard)
