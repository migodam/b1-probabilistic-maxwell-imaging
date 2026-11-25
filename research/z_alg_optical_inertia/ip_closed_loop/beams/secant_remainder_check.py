"""Recompute the native whole-box contraction factors used in the two-port remainder."""
import json

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, elapsed, geometry
from interval_certify import ProofGrid, hi
from native_paired import green_kernel_roundoff
from native_residual import ResidualBoundBox


def main():
    start = clocks()
    ctx.prec = 128
    ctx.threads = 3
    public = HERE.parent / "pair" / "power_96"
    command_rows = json.loads((public / "proposal.json").read_text())["rows"]
    grid = ProofGrid(geometry(96))
    kernel_roundoff = green_kernel_roundoff(grid)
    rows = {}
    with np.load(public / "public_classes.npz") as models:
        for key in command_rows:
            cells = models[key + "__cells"]
            widths = models[key + "__widths"]
            numeric, exact = grid.restricted(cells)
            factors = []
            for centre in models[key + "__centers"]:
                box = ResidualBoundBox(exact, numeric, centre, widths,
                                       kernel_roundoff * len(cells))
                factors.append({"T_times_width_upper": hi(box.eta),
                                "width_times_T_upper": hi(box.eta_reverse)})
            rows[key] = {"material_cells": len(cells), "classes": factors}
    out = {"source": "Original native Fourier operator and public old-data classes",
           "norm": "Frobenius upper bounds on the two operator norms, including numerical solve and kernel errors",
           "rows": rows, **elapsed(start)}
    atomic_write(HERE / "secant_remainder_96.json", out)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
