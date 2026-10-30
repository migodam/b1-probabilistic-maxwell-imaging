"""Fixed-work real 576-cell Maxwell throughput, not a material accuracy sweep."""
import os
for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
import multiprocessing as mp
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
import torch
from src.device import SharedOperators, ObjectiveLane


def load_public(old_root, observation):
    with np.load(Path(old_root) / "public_data/operators.npz") as f:
        p = {k: f[k] for k in f.files}
    with np.load(Path(old_root) / "public_data" / observation) as f:
        p.update({k: f[k] for k in f.files})
    return p


def legacy_module(root):
    sys.path.insert(0, str(Path(root) / "vendor"))
    spec = importlib.util.spec_from_file_location("r4_legacy", Path(root) / "src/__init__.py", submodule_search_locations=[str(Path(root) / "src")])
    m = importlib.util.module_from_spec(spec)
    sys.modules["r4_legacy"] = m
    spec.loader.exec_module(m)
    from r4_legacy.cuda_physics import CudaModel
    from r4_legacy.physics import Regularizer
    return CudaModel, Regularizer


def legacy_work(args):
    root, observation, objective, count, index = args
    torch.set_num_threads(1)
    p = load_public(root, observation)
    Model, Reg = legacy_module(root)
    model, reg = Model(p), Reg((24, 24))
    target, scale2 = p["Y"], float(np.vdot(p["Y"], p["Y"]).real)
    ref = model.forward(np.full(576, .1))
    norm = model.field_norm_squared(ref) / 12
    started = time.perf_counter()
    for i in range(count):
        # Public deterministic legal materials, no truth, digit label or template.
        x = .1 + .02 * (1 + np.sin(np.arange(576) * .07 + (i + index) * .11))
        state = model.forward(x)
        if objective == "dual":
            v, g, _ = model.smallgram(state, target)
            v, g = .5 * norm / scale2 * v, .5 * norm / scale2 * g
        elif objective == "fullj":
            g = model.jacobian(state)
            v = np.linalg.norm(state["y"] - target)
        else:
            v = .5 * np.linalg.norm(state["y"] - target) ** 2 / scale2
            g = model.data_gradient(state, target, scale2)
        if objective != "fullj":
            rv, rg, _ = reg.evaluate(x)
            v, g = v + rv, g + rg
    torch.cuda.synchronize()
    return {"evaluations": count, "worker_wall_s": time.perf_counter() - started,
            "last_objective": float(v), "counts": dict(model.counts)}


def stream_work(shared, p, objective, count, index):
    lane = ObjectiveLane(shared, p, objective=objective if objective != "fullj" else "data")
    started = time.perf_counter()
    for i in range(count):
        x = .1 + .02 * (1 + np.sin(np.arange(576) * .07 + (i + index) * .11))
        if objective == "fullj":
            lane.residual(x)
            lane.jacobian(x)
        else:
            lane.evaluate(x)
    return {"evaluations": count, "worker_wall_s": time.perf_counter() - started,
            "counts": dict(lane.counts), "host_result_wait_s": lane.wait_seconds}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy-root", required=True)
    ap.add_argument("--observation", default="obs_398f98f98cfe438e.npz")
    ap.add_argument("--evaluations", type=int, default=64)
    ap.add_argument("--output", default="results/throughput.json")
    a = ap.parse_args()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    p = load_public(a.legacy_root, a.observation)
    output = Path(a.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    begin = time.perf_counter()
    shared = SharedOperators(p)
    shared.ready.synchronize()
    setup_s = time.perf_counter() - begin
    gpu = torch.cuda.get_device_properties(0)
    result = {"hardware": {"name": gpu.name, "SMs": gpu.multi_processor_count,
                           "memory_bytes": gpu.total_memory, "torch": torch.__version__,
                           "cuda_build": torch.version.cuda, "platform": sys.platform},
              "scope": "fixed-work public-input Maxwell throughput; no GT or recovery claim",
              "dtype": "complex128", "physical_dimension": 576, "sources": 12, "receivers": 48,
              "resident_setup_seconds": setup_s, "setup_counts": shared.setup_counts,
              "rows": []}
    def save(row):
        result["rows"].append(row)
        output.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(row), flush=True)
    # Warm each tested API once; not an accuracy/precision test.
    for obj in ("data", "dual", "fullj"):
        legacy_work((a.legacy_root, a.observation, obj, 1, 0))
        stream_work(shared, p, obj, 1, 0)
    for objective in ("data", "dual", "fullj"):
        count = a.evaluations if objective != "fullj" else a.evaluations // 2
        for processes in (1, 2):
            t = time.perf_counter()
            args = [(a.legacy_root, a.observation, objective, count // processes, i) for i in range(processes)]
            if processes == 1:
                work = [legacy_work(args[0])]
                cold = 0.
            else:
                with mp.get_context("spawn").Pool(processes) as pool:
                    work = pool.map(legacy_work, args)
                # Parent wall explicitly includes process start/import and setup.
                cold = None
            wall = time.perf_counter() - t
            save({"kind": "legacy_processes", "objective": objective, "concurrency": processes,
                  "evaluations": count, "wall_seconds": wall, "evaluations_per_second": count / wall,
                  "includes_worker_setup": True, "workers": work})
        for streams in (1, 2, 4):
            t = time.perf_counter()
            with ThreadPoolExecutor(max_workers=streams) as pool:
                work = list(pool.map(lambda i: stream_work(shared, p, objective, count // streams, i), range(streams)))
            wall = time.perf_counter() - t
            save({"kind": "resident_streams", "objective": objective, "concurrency": streams,
                  "evaluations": count, "wall_seconds": wall, "evaluations_per_second": count / wall,
                  "includes_worker_setup": True, "shared_setup_separately": setup_s, "workers": work})
    # Batched candidate scoring, a different workload from objective+gradient.
    target = torch.as_tensor(p["Y"], device=shared.device)
    for batch in (1, 2, 4, 8, 16):
        materials = torch.as_tensor(np.stack([.1 + .02 * (1 + np.sin(np.arange(576) * .07 + i * .11)) for i in range(batch)]), device=shared.device)
        stream = torch.cuda.Stream()
        stream.wait_event(shared.ready)
        with torch.cuda.stream(stream):
            value, info = shared.score_batch(materials, target)
        stream.synchronize()
        t = time.perf_counter()
        with torch.cuda.stream(stream):
            for i in range(max(1, a.evaluations // batch)):
                value, info = shared.score_batch(materials, target)
        stream.synchronize()
        wall = time.perf_counter() - t
        n = max(1, a.evaluations // batch) * batch
        save({"kind": "batched_candidate_score", "objective": "forward_only", "concurrency": batch,
              "evaluations": n, "wall_seconds": wall, "evaluations_per_second": n / wall,
              "final_lu_info": info.cpu().tolist(), "scope": "GPU resident inputs and outputs; no optimizer host return"})
    # One targeted CUDA Graph attempt. Preserve an actual capture failure.
    try:
        batch = 4
        materials = materials[:batch].clone()
        stream = torch.cuda.Stream()
        stream.wait_stream(torch.cuda.current_stream())
        with torch.cuda.stream(stream):
            for _ in range(3):
                shared.score_batch(materials, target)
        stream.synchronize()
        graph = torch.cuda.CUDAGraph()
        with torch.cuda.graph(graph, stream=stream):
            graph_value, graph_info = shared.score_batch(materials, target)
        t = time.perf_counter()
        for _ in range(a.evaluations // batch):
            graph.replay()
        torch.cuda.current_stream().synchronize()
        wall = time.perf_counter() - t
        save({"kind": "cuda_graph_candidate_score", "objective": "forward_only", "concurrency": batch,
              "evaluations": a.evaluations, "wall_seconds": wall,
              "evaluations_per_second": a.evaluations / wall,
              "final_lu_info": graph_info.cpu().tolist(), "status": "OK"})
    except RuntimeError as exc:
        save({"kind": "cuda_graph_candidate_score", "status": "CAPTURE_FAILED", "error": str(exc)})
    result["memory_allocated_bytes"] = torch.cuda.max_memory_allocated()
    result["memory_reserved_bytes"] = torch.cuda.max_memory_reserved()
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
