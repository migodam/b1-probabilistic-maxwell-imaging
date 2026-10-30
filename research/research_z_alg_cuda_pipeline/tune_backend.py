"""One focused library-choice experiment, including an actual graph attempt."""
import os
for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[name] = "1"
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import time
import numpy as np
import torch
from benchmark import load_public, stream_work
from src.device import SharedOperators


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", required=True)
    p.add_argument("--output", default="results/backend_choice.json")
    a = p.parse_args()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    public = load_public(a.legacy_root, "obs_398f98f98cfe438e.npz")
    rows = []
    out = Path(a.output)
    def save(row):
        rows.append(row)
        out.write_text(json.dumps({"rows": rows, "dtype": "complex128", "profiler_supported_activities": [str(x) for x in torch.profiler.supported_activities()]}, indent=2))
        print(json.dumps(row), flush=True)
    for backend in ("default", "cusolver"):
        torch.backends.cuda.preferred_linalg_library(backend)
        shared = SharedOperators(public)
        shared.ready.synchronize()
        stream_work(shared, public, "data", 1, 0)
        t = time.perf_counter()
        with ThreadPoolExecutor(max_workers=4) as pool:
            work = list(pool.map(lambda i: stream_work(shared, public, "data", 32, i), range(4)))
        wall = time.perf_counter() - t
        save({"backend": backend, "workload": "data_gradient_4streams", "evaluations": 128, "wall_s": wall, "evaluations_per_second": 128 / wall})
        materials = torch.full((16, 576), .12, dtype=torch.float64, device=shared.device)
        target = torch.as_tensor(public["Y"], device=shared.device)
        stream = torch.cuda.Stream()
        stream.wait_event(shared.ready)
        with torch.cuda.stream(stream):
            shared.score_batch(materials, target)
        stream.synchronize()
        t = time.perf_counter()
        with torch.cuda.stream(stream):
            for _ in range(8):
                values, info = shared.score_batch(materials, target)
        stream.synchronize()
        wall = time.perf_counter() - t
        save({"backend": backend, "workload": "candidate_score_batch16", "evaluations": 128, "wall_s": wall, "evaluations_per_second": 128 / wall})
        if backend == "cusolver":
            try:
                static = materials[:4].clone()
                stream.wait_stream(torch.cuda.current_stream())
                with torch.cuda.stream(stream):
                    for _ in range(3):
                        shared.score_batch(static, target)
                stream.synchronize()
                graph = torch.cuda.CUDAGraph()
                with torch.cuda.graph(graph, stream=stream):
                    graph_values, graph_info = shared.score_batch(static, target)
                t = time.perf_counter()
                for _ in range(32):
                    graph.replay()
                torch.cuda.current_stream().synchronize()
                wall = time.perf_counter() - t
                save({"backend": backend, "workload": "candidate_graph_batch4", "status": "OK", "evaluations": 128, "wall_s": wall, "evaluations_per_second": 128 / wall, "lu_info": graph_info.cpu().tolist()})
            except RuntimeError as error:
                save({"backend": backend, "workload": "candidate_graph_batch4", "status": "CAPTURE_FAILED", "error": str(error)})


if __name__ == "__main__":
    main()
