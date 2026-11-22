"""Time three proposal strategies on held-out observations, prior to native replay."""

import argparse
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
from time import perf_counter, process_time

for key in (
    "VECLIB_MAXIMUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from coarse import Coarse
from holdout import observations
from network import load_network, predict

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"


def compact(x):
    return json.dumps(
        x,
        allow_nan=False,
        separators=(",", ":"),
        default=lambda a: a.item() if isinstance(a, np.generic) else a.tolist(),
    )


def save(stream, x):
    stream.write(compact(x) + "\n")
    stream.flush()


def run(shard, supplement=False):
    if not (R / "holdout_open.json").exists():
        raise RuntimeError("Freeze training before evaluating the held-out fields.")
    destination = R / "supplement" if supplement else R
    destination.mkdir(parents=True, exist_ok=True)
    k, n = map(int, shard.split("/"))
    start, cpu = perf_counter(), process_time()
    torch.set_num_threads(1)
    out = destination / "proposals"
    out.mkdir(exist_ok=True)
    models = {
        t: load_network(R / p)
        for t, p in [
            ("neural", "neural_dual_weights.pt"),
            ("no_feature", "no_feature_weights.pt"),
        ]
    }
    co = Coarse()
    cases = list(observations(supplement=supplement))
    objects = sorted({q[0] for q in cases})[k::n]
    details = destination / "optimization_traces"
    details.mkdir(exist_ok=True)
    init = {
        "worker": shard,
        "wall": perf_counter() - start,
        "cpu": process_time() - cpu,
        "coarse": co.setup_cost,
        "opened": datetime.now(timezone.utc).isoformat(),
    }
    (out / f"setup_{k}_{n}.json").write_text(compact(init))
    with threadpool_limits(1):
        for oid in objects:
            file = out / f"{oid}.jsonl"
            done = set()
            if file.exists():
                old = [json.loads(v) for v in file.read_text().splitlines() if v]
                done = {r["id"] for r in old if "id" in r}
                if any(r.get("_complete") for r in old):
                    continue
            with (
                file.open("a") as output,
                gzip.open(details / f"{oid}.jsonl.gz", "at") as trace,
            ):
                for obj, noise, r, y, caps in cases:
                    if oid != obj:
                        continue
                    for u in caps:
                        for side in (4, 8):

                            def record(
                                arm, mu, lam, pwall, pcpu, d2=None, nu=None, **kw
                            ):
                                ident = f"{oid}__n{noise}__u{u:g}__b{side}__{arm}"
                                if ident in done:
                                    return
                                save(
                                    output,
                                    dict(
                                        id=ident,
                                        oid=oid,
                                        noise=noise,
                                        radius=r,
                                        u=u,
                                        side=side,
                                        arm=arm,
                                        mu=np.asarray(mu).flatten().tolist(),
                                        lam=np.asarray(lam).flatten().tolist(),
                                        proposal_wall=pwall,
                                        proposal_cpu=pcpu,
                                        coarse_d2=d2,
                                        nu=nu,
                                        **kw,
                                    ),
                                )
                                done.add(ident)

                            weights = {}
                            for arm, m in models.items():
                                t, c = perf_counter(), process_time()
                                mu, lam = predict(m, y, u, side, noise_radius=r)
                                dt, dc = perf_counter() - t, process_time() - c
                                weights[arm] = (mu, lam, dt, dc)
                                # This evaluation is an offline diagnostic. The zero-step deployment uses only predict().
                                rt, rc = perf_counter(), process_time()
                                diag = co.eval(y, u, mu, lam, smooth=False)
                                record(
                                    arm,
                                    mu,
                                    lam,
                                    dt,
                                    dc,
                                    diag["d2_exact"],
                                    diag["nu"],
                                    diagnostic_wall=perf_counter() - rt,
                                    diagnostic_cpu=process_time() - rc,
                                )
                            mu, lam, dt, dc = weights["neural"]
                            t, c = perf_counter(), process_time()
                            hy = co.refine(y, u, (mu, lam), steps=3)
                            hy_w, hy_c = perf_counter() - t, process_time() - c
                            best = hy["best"]
                            record(
                                "hybrid",
                                best["mu"],
                                best["lam"],
                                dt + hy_w,
                                dc + hy_c,
                                best["d2"],
                                best["nu"],
                                infer_wall=dt,
                                refine_wall=hy_w,
                                evals=hy["eval_count_total"],
                            )
                            save(
                                trace,
                                dict(
                                    oid=oid,
                                    noise=noise,
                                    u=u,
                                    side=side,
                                    arm="hybrid",
                                    **hy,
                                ),
                            )
                            dd = co.optimize(
                                y, u, blocks=side * side, time_limits=(1.0, 5.0)
                            )
                            save(
                                trace,
                                dict(
                                    oid=oid,
                                    noise=noise,
                                    u=u,
                                    side=side,
                                    arm="gradient",
                                    **dd,
                                ),
                            )
                            for sec, arm in ((1.0, "gradient1"), (5.0, "gradient5")):
                                q = dd["prefixes"][sec]
                                b = q["best"]
                                wall = (
                                    dd["paid_wall_s"]
                                    if sec == 5
                                    else q["wall_budget_used_s"]
                                )
                                cpu1 = (
                                    dd["paid_cpu_s"]
                                    if sec == 5
                                    else q["completed_cpu_s"]
                                )
                                record(
                                    arm,
                                    b["mu"],
                                    b["lam"],
                                    wall,
                                    cpu1,
                                    b["d2"],
                                    b["nu"],
                                    evals=q["eval_count"],
                                    time_limit=sec,
                                    completion_cpu_lower=q["completed_cpu_s"],
                                    last_wall=q["completed_wall_s"],
                                    overrun=dd["over_budget_evals"]
                                    if sec == 5
                                    else None,
                                )
                save(output, dict(_complete=True, oid=oid, records=len(done)))
            print(oid, len(done), round(perf_counter() - start, 1), flush=True)
    (out / f"worker_{k}_{n}.json").write_text(
        compact(
            dict(wall=perf_counter() - start, cpu=process_time() - cpu, objects=objects)
        )
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--supplement", action="store_true")
    args = ap.parse_args()
    run(args.shard, args.supplement)
