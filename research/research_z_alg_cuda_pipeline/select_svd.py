"""Select an existing mature SVD operator on actual public-data TRF matrices."""
import os
for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[k] = "1"
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
import torch
from src.device import SharedOperators, ObjectiveLane
from benchmark import load_public


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy-root", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.legacy_root) / "vendor"))
    from scipy.linalg import svd
    from threadpoolctl import threadpool_limits
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    public = load_public(a.legacy_root, "obs_398f98f98cfe438e.npz")
    shared = SharedOperators(public)
    rows = []
    with threadpool_limits(limits=1):
        for stage in (8, 12, 24):
            lane = ObjectiveLane(shared, public, stage=stage)
            jac = lane.jacobian(np.full(stage * stage, .1))
            # Trust-region augmented shape, with a public positive diagonal.
            matrix = np.vstack((jac, np.eye(stage * stage) * .001))
            t = time.perf_counter()
            u, s, v = svd(matrix, full_matrices=False, check_finite=False)
            rows.append({"stage": stage, "shape": list(matrix.shape), "backend": "scipy_gesdd_cpu_1thread", "wall_s": time.perf_counter() - t})
            for driver in ("gesvd", "gesvdj"):
                t = time.perf_counter()
                x = torch.as_tensor(matrix, device=shared.device)
                u, s, v = torch.linalg.svd(x, full_matrices=False, driver=driver)
                # A CPU trust-region solver needs all three results, not just s.
                u.cpu(); s.cpu(); v.cpu()
                rows.append({"stage": stage, "shape": list(matrix.shape), "backend": "cusolver_" + driver, "wall_s": time.perf_counter() - t,
                             "timing_scope": "input upload + SVD + all U/s/Vh host returns; no precision change"})
            print(json.dumps(rows[-3:]), flush=True)
    Path("results/svd_choice.json").write_text(json.dumps({"rows": rows, "scope": "one public full nonlinear Jacobian per stage; kernel selection only; cold API cost included"}, indent=2))


if __name__ == "__main__":
    main()
