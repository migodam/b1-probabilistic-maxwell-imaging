"""Refit original archived A and B acquisitions independently of their historical fit."""
import argparse
from collections import Counter
from time import perf_counter, process_time

from continuation import (HERE, FROZEN, Model, VerifiedWave, BOX, compensation, np,
                          solve, threadpool_limits, write_arrays, write_json, rd)
from b_control import dark_apertures, layouts, one
from aperture_joint import Joint, display

CASES = ("vacuum", "high_patch")
METHODS = ("new_B_endpoint", "prior_joint_endpoint", "old_A_endpoint")


def elapsed(start):
    return dict(wall_s=perf_counter()-start[0], cpu_s=process_time()-start[1])


def bounded_direction(jac, residual, low, high, trust, damping):
    """Damped data-space direction with trial bounds enforced inside the linear fit."""
    lo, hi = np.maximum(-trust, low), np.minimum(trust, high)
    gradient = jac.T @ residual
    free = ~((lo >= -1e-11) & (gradient >= 0)) & ~((hi <= 1e-11) & (gradient <= 0))
    free &= hi > lo
    result = np.zeros_like(lo)
    for update in range(5):
        if not free.any():
            break
        jf = jac[:, free]
        remaining = residual + jac[:, ~free] @ result[~free]
        gram = jf @ jf.T
        curvature = max(np.trace(gram)/len(gram), 1e-18)
        raw = -jf.T @ solve(gram + damping*curvature*np.eye(len(gram)),
                            remaining, assume_a="pos")
        now = np.flatnonzero(free)
        result[free] = np.clip(raw, lo[free], hi[free])
        saturated = (raw < lo[free]) | (raw > hi[free])
        if not saturated.any():
            break
        free[now[saturated]] = False
    return result


def fit(wave, x, low, high, max_major):
    start = perf_counter(), process_time()
    cm = Counter(wave.counts)
    x = np.clip(x, low, high)
    orig = wave.forward(x)
    state = orig
    trace = []
    damping = .005
    trust = .65
    stop = "major-iteration cap"
    best = state
    for it in range(max_major):
        worst = max(state["ratios"])
        if worst <= .89:
            stop = "both archived balls have numerical slack"
            break
        jac = wave.jacobian(state)
        # Keep the target acquisition weights fixed throughout each line search.
        ratios = np.asarray(state["ratios"])
        weights = np.power(ratios / max(ratios), .75).clip(.4, 1.)
        per_row = np.tile(np.repeat(weights, 64), 2)
        trace_row = dict(iteration=it+1, base_ratios=list(state["ratios"]),
                         damping=damping, trust_cell=trust, acquisition_weight=weights.tolist(), trials=[])
        residual = per_row * state["residual"]
        old = float(np.linalg.norm(residual))
        direction = bounded_direction(per_row[:, None]*jac, residual, low-x, high-x, trust, damping)
        for a in (1., .5, .25, .125, .0625):
            z = np.clip(x + a*direction, low, high)
            trial = wave.forward(z)
            value = float(np.linalg.norm(per_row * trial["residual"]))
            trace_row["trials"].append(dict(step=a, A_over_r=trial["ratios"][0],
                                           B_over_r=trial["ratios"][1], weighted_residual_l2=value))
            if max(trial["ratios"]) < max(best["ratios"]):
                best = trial
            if value < old*(1-1e-9):
                x, state = z, trial
                trust = min(1.3, trust*1.12)
                damping = max(7e-6, damping*.73)
                break
        else:
            trust *= .72
            damping = min(3., damping*3.)
        trace.append(trace_row)
        if trust < .004:
            stop = "linear trust-region floor"
            break
    return best, dict(initial=display(orig), best=display(best), traces=trace,
                      stop=stop, counts=dict(wave.counts-cm), **elapsed(start))


def physical_receipt(native, wave, x):
    first = native.forward(x, fresh=True)
    receipt, current = one(wave, native, x, first["j"])
    error = receipt["data_rounding_solve_upper"]/native.radius
    receipt["residual_lower_over_r"] = max(0., (receipt["numeric"]-error)
                             * (1-rd(native.Y.size+2)))
    return receipt, current, native.S@current


def main(args):
    start = perf_counter(), process_time()
    archived = FROZEN/"results"/"aperture_control"
    with np.load(archived/"geometry.npz", allow_pickle=False) as a:
        geom = dict(a)
    with np.load(archived/"n32_c01.npz", allow_pickle=False) as a:
        data = {k:a[k] for k in ("A_one_way_Y","A_one_way_r","B_interleaved_Y","B_interleaved_r")}
        old_a_vac = a["zero"]
        old_a_high = a["patch_candidate"]
    with np.load(HERE/"B_independent_control.npz", allow_pickle=False) as a:
        new_b_vac = a["B_previous_vacuum_vacuum_material"]
        new_b_high = a["B_previous_positive_high_patch_material"]
        roi, patch = a["roi"].ravel().astype(bool), np.sort(a["patch"])
    with np.load(FROZEN/"results"/"aperture_joint"/"n32_c01.npz", allow_pickle=False) as a:
        prior_joint_vac = a["roi_zero_x"].ravel()
        prior_joint_high = a["connected_patch_floor_2_x"].ravel()
    with np.load(FROZEN/"ports"/"public"/"operators"/"n32_f400_a08.npz", allow_pickle=False) as a:
        native = dict(a)
    _, exc, rec = dark_apertures(native)
    _, _, _, b_coeff = layouts(native, exc, rec)["B_interleaved"]
    a_model = Model.read("n32_c01")
    if not np.array_equal(a_model.Y, data["A_one_way_Y"]) or a_model.radius != float(data["A_one_way_r"]):
        raise ValueError("Archived A differs from the original observation")
    noisy_b = data["B_interleaved_Y"]
    r_b = float(data["B_interleaved_r"])
    b_model = Model(b_coeff, noisy_b, r_b, "n32_f400_a08 / B", 32)
    vw_a, vw_b = VerifiedWave(native), VerifiedWave(b_coeff)
    engine = Joint(geom, data)
    choices = {"vacuum":dict(zip(METHODS,(new_b_vac,prior_joint_vac,old_a_vac))),
               "high_patch":dict(zip(METHODS,(new_b_high,prior_joint_high,old_a_high)))}
    report = dict(case="n32_c01", A_observation="bmtm26_n32_c01_f400_a08_c48",
                  B_observation="original archived aperture_control/n32_c01 B_interleaved_Y",
                  frequency_MHz=400, A_radius=a_model.radius, B_radius=b_model.radius,
                  original_A_and_B=True, observation_contract="two separate deterministic Frobenius balls",
                  pulses=16, reads_complex=128, source_column_power_sum_W=2., pulse_dwell_s=None,
                  patch_cells=len(patch), roi_cells=int(roi.sum()), area_cell_m2=(2/32)**2,
                  trials={}, selected={}, heldout=None)
    results = dict(Y_A=a_model.Y,Y_B=b_model.Y,roi=roi,patch=patch)
    selected_x = {}
    area = (2/32)**2
    out = HERE/(f"joint_original_AB{args.label}.json")
    arrays_out = out.with_suffix(".npz")

    def checkpoint():
        report["cost"] = dict(**elapsed(start), optimizer=dict(engine.counts),
                              original_A_audit=dict(a_model.cost), original_B_audit=dict(b_model.cost))
        write_json(out, report)
        write_arrays(arrays_out, **results)

    for branch in args.families:
        low = np.zeros(engine.n)
        high = np.full(engine.n, BOX)
        if branch == "vacuum":
            high[roi] = 0.
        else:
            low[patch] = 2.
        report["trials"][branch] = {}
        ranked = []
        for name in args.seeds:
            trial, info = fit(engine, choices[branch][name], low, high, args.max_major)
            x = trial["x"]
            if np.any(x<low) or np.any(x>high):
                raise ValueError("Optimizer did not respect the frozen class")
            rec_a, ja, fa = physical_receipt(a_model, vw_a, x)
            rec_b, jb, fb = physical_receipt(b_model, vw_b, x)
            worst = max(rec_a["residual_upper_over_r"],rec_b["residual_upper_over_r"])
            strict = bool(worst<1.)
            report["trials"][branch][name] = dict(**info, native_A=rec_a,
                                                 native_B=rec_b, jointly_verified=strict)
            key = branch+"_"+name
            results[key+"_material"] = x
            results[key+"_A_current"] = ja
            results[key+"_B_current"] = jb
            results[key+"_A_prediction"] = fa
            results[key+"_B_prediction"] = fb
            ranked.append((worst, name, x))
            hit = min(ranked,key=lambda v:v[0])
            selected_x[branch] = hit[2]
            report["selected"][branch] = dict(seed=hit[1], maximum_A_B_residual_upper_over_r=hit[0],
                jointly_verified=hit[0]<1., patch_min=float(hit[2][patch].min()),
                patch_mean=float(hit[2][patch].mean()), roi_mean=float(hit[2][roi].mean()))
            print(branch, name, rec_a["residual_upper_over_r"],rec_b["residual_upper_over_r"],
                  strict, elapsed(start), flush=True)
            checkpoint()
            if strict:
                break
    if set(selected_x)==set(CASES):
        first, second = (selected_x[k] for k in CASES)
        report["exterior_cost_of_selected_pair"] = compensation(second[~roi]-first[~roi],area)
        if all(q["jointly_verified"] for q in report["selected"].values()):
            # The two material arrays and their native certification are written before
            # opening the held measurement. These values never select or modify either map.
            with np.load(FROZEN/"ports"/"public"/"cases"/"n32_c01.npz",allow_pickle=False) as a:
                yh, rh = a["Y_hold"], float(a["radius_hold"])
            h_coeff = dict(native, S=native["S_hold"])
            h_model = Model(h_coeff,yh,rh,"n32_f400_a08 / independent held receiver arc",32)
            h_wave = VerifiedWave(h_coeff)
            report["heldout"] = {}
            results["held_observation"] = yh
            for name, x in selected_x.items():
                check, curr, pred = physical_receipt(h_model,h_wave,x)
                report["heldout"][name] = check
                results["held_"+name+"_curr"] = curr
                results["held_"+name+"_prediction"] = pred
    checkpoint()


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--families",nargs="+",choices=CASES,default=list(CASES))
    parser.add_argument("--seeds",nargs="+",choices=METHODS,default=list(METHODS))
    parser.add_argument("--max-major",type=int,default=120)
    parser.add_argument("--label",default="")
    args=parser.parse_args()
    with threadpool_limits(1):
        main(args)
