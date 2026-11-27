"""Independent original-grid zero/patch B refit using only its archived noisy observation."""
from pathlib import Path
from time import perf_counter, process_time

from continuation import (HERE, FROZEN, Model, VerifiedWave, BOX, np, realify,
                          solve, threadpool_limits, write_arrays, write_json)
from aperture_control import dark_apertures, layouts
from native_replay import one


def fit(model, x, lb, ub, max_updates=40):
    start = perf_counter(), process_time()
    current = np.clip(x, lb, ub)
    rows = []
    damping = .012
    for i in range(max_updates):
        state = model.forward(current)
        old = float(state["residual_ratio"])
        if old < .85:
            break
        diffuse = model.solve(current, model.SH, adjoint=True, rtol=2e-12).conj().T
        gradient = realify((diffuse[:, None, :] * state["field"].T[None, :, :])
                           .reshape(model.Y.size, model.n))
        target = realify((model.Y - state["pred"]).ravel())
        slope = gradient.T @ target
        free = (ub > lb) & (
            ((current > lb+1e-10) & (current < ub-1e-10)) |
            ((current <= lb+1e-10) & (slope > 0.)) |
            ((current >= ub-1e-10) & (slope < 0.)))
        if not free.any():
            break
        y = gradient[:, free]
        a = y @ y.T
        regularizer = damping * max(np.trace(a)/len(a), 1e-18)
        d = np.zeros(model.n)
        d[free] = y.T @ solve(a + regularizer*np.eye(len(a)), target, assume_a="pos")
        d *= min(1., .8 / max(abs(d).max(), 1e-30))
        best = old
        for frac in (1., .5, .25, .125, .0625):
            candidate = np.clip(current + frac*d, lb, ub)
            ratio = float(model.forward(candidate)["residual_ratio"])
            if ratio < old:
                best = ratio
                current = candidate
                damping = max(2e-6, damping * .59)
                break
        else:
            damping *= 3.
        rows.append(dict(update=i+1, old_ratio=old, accepted_ratio=best, free=int(free.sum()), damping=damping))
        if damping >= 10:
            break
    state = model.forward(current)
    return current, dict(numeric_over_r=state["residual_ratio"], trace=rows,
                         wall_s=perf_counter()-start[0], cpu_s=process_time()-start[1])


def main():
    t0, cpu0 = perf_counter(), process_time()
    p = FROZEN / "results" / "aperture_control" / "n32_c01.npz"
    with np.load(p, allow_pickle=False) as a:
        y, radius = a["B_interleaved_Y"], float(a["B_interleaved_r"])
        roi, patch = a["roi"].ravel(), a["patch"]
        seeds = {"A_vacuum":a["zero"], "B_previous_vacuum":a["B_interleaved_refit_zero"],
                 "B_previous_positive":a["B_interleaved_refit_patch"]}
    op_id = "n32_f400_a08"
    with np.load(FROZEN / "ports" / "public" / "operators" / f"{op_id}.npz",
                 allow_pickle=False) as a:
        public = dict(a)
    _, exciting, receiving = dark_apertures(public)
    tx, rx, scale, layout = layouts(public, exciting, receiving)["B_interleaved"]
    model = Model(layout, y, radius, "n32_f400_a08 / frozen B_interleaved", 32)
    wave = VerifiedWave(layout)
    output = {}
    arrays = dict(Y_B=y, radius_B=radius, patch=patch, roi=roi, source_ports=tx, receive_ports=rx)
    for name, start in seeds.items():
        for branch in ("vacuum", "high_patch"):
            if branch == "high_patch" and name != "B_previous_positive":
                continue
            lb, ub = np.zeros(model.n), np.full(model.n, BOX)
            if branch == "vacuum":
                ub[roi] = 0.
            else:
                lb[patch] = 2.
            idx = name + "_" + branch
            candidate, timing = fit(model, start, lb, ub)
            state = model.forward(candidate, fresh=True)
            strict, replay_j = one(wave, model, candidate, state["j"])
            output[idx] = dict(**timing, verified=strict, stop="discrepancy or bounded fit budget")
            arrays[idx+"_material"] = candidate
            arrays[idx+"_current"] = replay_j
            arrays[idx+"_prediction"] = model.S @ replay_j
            print(idx, timing["numeric_over_r"], strict["residual_upper_over_r"], flush=True)
    write_json(HERE / "B_independent_control.json", dict(case="n32_c01", radius=radius, responses=output,
            raw_input=str(p.relative_to(FROZEN)), shots=8, reads=64, frequency_MHz=400,
            roles=dict(tx=tx.tolist(),rx=rx.tolist(),source_scaling=scale),
            total=dict(wall_s=perf_counter()-t0, cpu_s=process_time()-cpu0, native_counts=dict(model.cost))))
    write_arrays(HERE / "B_independent_control.npz", **arrays)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
