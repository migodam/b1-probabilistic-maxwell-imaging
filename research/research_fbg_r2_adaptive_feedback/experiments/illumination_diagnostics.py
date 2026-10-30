"""Run the fixed public Stage 3 illumination compatibility diagnostic on CPU."""
import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
from src.source_compatibility import analyze


def main():
    started = time.perf_counter()
    operator_dir = PROJECT / "data_contract" / "stage3_dev" / "operators"
    paths = sorted(operator_dir.glob("*.npz"))
    expected = {f"{layout}_p{p}" for layout in ("uniform", "random") for p in (3, 6, 12, 24)}
    if {path.stem for path in paths} != expected:
        raise RuntimeError("public Stage 3 operator set is not the expected eight acquisitions")
    side = 24
    yy, xx = np.meshgrid(np.linspace(0., 1., side), np.linspace(0., 1., side), indexing="ij")
    h = (1 + .25 * np.cos(np.pi * xx) * np.cos(np.pi * yy)).reshape(-1)
    h = h / np.max(np.abs(h))
    rows = []
    for path in paths:
        with np.load(path, allow_pickle=False) as z:
            row = analyze(z["G"], z["E0"], h, anchor=.1, step=.5)
        row["acquisition"] = path.stem
        row["operator_path"] = str(path.relative_to(PROJECT))
        rows.append(row)
    total = {key: sum(row["costs"][key] for row in rows) for key in rows[0]["costs"]}
    result = {
        "status": "COMPLETE", "utc": datetime.now(timezone.utc).isoformat(),
        "scope": "eight public Stage 3 scalar-TM acquisition operators; CPU mechanism diagnostic only",
        "direction": "h=(1+0.25 cos(pi x)cos(pi y))/max_abs on row-major 24x24 unit grid",
        "finite_current_proposal": "w=0.5 T, T=L^-1 diag(h)Q at m=0.1",
        "source_whitening": "Q and T right-multiplied by chol(Q*Q)^-*; no jitter",
        "compatibility_equation": "Lw=diag(u)(Q+Gw), shared real u rowwise least squares",
        "apparent_real_pass_rule": "global real compatibility residual <=1e-6 and no zero Q or finite V row",
        "projector_tangent": "exact dP_Q norm squared=2 tr(U*U H^-1), U=(I-P_Q)GT",
        "rational_increment_scope": "u=t h/(1+t gamma) only if GT_i=gamma_i Q_i and resulting u is real",
        "rows": rows, "total_costs": total, "total_wall_seconds": time.perf_counter() - started,
        "scientific_limit": "Exact infinitesimal current tangency does not establish finite shared-real-material compatibility. No imaging performance or A17 implementation gain was measured."
    }
    output = PROJECT / "results" / "illumination_current_compatibility.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(output)
    print("rows", len(rows), "wall_s", round(result["total_wall_seconds"], 3))


if __name__ == "__main__":
    main()
