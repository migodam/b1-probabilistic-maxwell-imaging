"""Transfer the two coarse-selected, equal-cost MIMO material pairs to a finer VIE."""
import json
from time import perf_counter, process_time

import numpy as np
from scipy.special import ndtr
from threadpoolctl import threadpool_limits

from acquire import fft_field
from designs import all_designs, SIGMA_FRACTION
from material_witness import OUT
from ndt import Air, BACKGROUND


def main():
    t0, c0 = perf_counter(), process_time()
    proposals = json.loads((OUT / "material_weak_modes.json").read_text())["cases"]
    arms = all_designs()
    names = {"f2_t4_r8", "f4_t2_r8"}
    proposals = [p for p in proposals if {p["hard_design"], p["easy_design"]} == names]
    with np.load(OUT / "acquisition.npz") as pack:
        sigma = float(pack["reference_amplitude"]) * SIGMA_FRACTION
    n = 72
    frequencies = sorted(set(f for name in names for f in arms[name].freqs))
    physics = {f: Air(n, f) for f in frequencies}
    operators = {f: model.convolution() for f, model in physics.items()}
    cases = []
    for p in proposals:
        signals = {}
        max_residual = 0.
        ky = []
        for f, air in physics.items():
            outputs = []
            for member in ("minus", "plus"):
                pixels = np.zeros(24 * 24)
                pixels[p["pixel_ids"]] = p[member]
                q = pixels.reshape(24, 24).repeat(n // 24, 0).repeat(n // 24, 1)
                material = q.ravel() + BACKGROUND
                field, log = fft_field(air, operators[f], material)
                outputs.append(air.S @ (material[:, None] * field))
                max_residual = max(max_residual, log["max_relative_residual"])
                ky.append(log["krylov_iterations"])
            signals[f] = outputs[1] - outputs[0]

        def measure(name):
            d = arms[name]
            norm_sq = sum(np.linalg.norm(signals[f][np.ix_(d.rx, d.tx)]) ** 2 for f in d.freqs)
            delta = float(np.sqrt(d.shot_energy * norm_sq) / sigma)
            return {"gap_over_sigma": delta, "equal_prior_error": float(ndtr(-delta/np.sqrt(2.)))}

        g72 = {"hard": measure(p["hard_design"]), "easy": measure(p["easy_design"])}
        cases.append(dict(hard_design=p["hard_design"], easy_design=p["easy_design"],
                          patch_rms=p["patch_rms"], grid={"48": {"hard": p["fine_hard"], "easy": p["fine_easy"]}, "72": g72},
                          max_relative_linear_residual=max_residual, krylov_iterations=sum(ky)))
        print(p["hard_design"], g72, flush=True)
    payload = dict(cases=cases, grid_points=n*n, equal_energy=1., shots=8, complex_readouts=64,
                   evaluation="Unchanged two material contrasts selected at 24 cells per side. Direct FFT VIE without pivot or tangent approximation at 72 per side.",
                   wall_s=perf_counter()-t0, cpu_s=process_time()-c0)
    (OUT / "witness_transfer.json").write_text(json.dumps(payload, indent=2)+"\n")


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
