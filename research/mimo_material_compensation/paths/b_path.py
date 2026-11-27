"""From a B-feasible vacuum, verify a B-only patch trajectory and independently audit A."""
from time import perf_counter, process_time

from continuation import (HERE, FROZEN, Model, VerifiedWave, Reduction, select_exterior,
                          run_path, compensation, exact_chart, np, threadpool_limits,
                          write_json, write_arrays, up, rd)
from b_control import dark_apertures, layouts, one


def main():
    start = perf_counter(), process_time()
    with np.load(HERE / "B_independent_control.npz", allow_pickle=False) as src:
        vacuum = src["B_previous_vacuum_vacuum_material"]
        high = src["B_previous_positive_high_patch_material"]
        roi, patch = src["roi"].ravel().astype(bool), np.sort(src["patch"])
        y = src["Y_B"]
        radius = float(src["radius_B"])
    original = Model.read("n32_c01")
    with np.load(FROZEN / "ports" / "public" / "operators" / "n32_f400_a08.npz",
                 allow_pickle=False) as src:
        native = dict(src)
    _, exc, rec = dark_apertures(native)
    tx, rx, scale, coeff = layouts(native, exc, rec)["B_interleaved"]
    b = Model(coeff, y, radius, "n32_f400_a08 / frozen B_interleaved", 32)
    b_wave = VerifiedWave(coeff)
    vp = b.forward(vacuum, fresh=True)
    receiver = b.solve(vacuum, b.SH, adjoint=True, rtol=2e-12).conj().T
    outside, spectrum = select_exterior(b, vacuum, vp["field"], receiver, roi, 256)
    support = np.r_[patch, outside]
    is_patch = np.isin(support, patch)
    reduced = Reduction.build(b, b_wave, vacuum, support, vp["j"], receiver)
    area = (2/32)**2
    trajectory, arr, knots = run_path(reduced, is_patch, roi, "ball_min_margin", 2., 2.**-10, area)
    exact, exact_positions = exact_chart(reduced, is_patch)
    a_wave = VerifiedWave(native)
    ab = {}
    native_outputs = {}
    independent_npz = {}
    for name, x in (("vacuum", vacuum), ("high_patch", high), ("trajectory_terminal", arr[-1])):
        ab[name] = {}
        for title, model, wav in (("A", original, a_wave), ("B", b, b_wave)):
            first = model.forward(x, fresh=True)
            receipt, current = one(wav, model, x, first["j"])
            numeric = receipt["numeric"]
            outward = receipt["data_rounding_solve_upper"]/model.radius
            receipt["residual_lower_over_r"] = max(0., (numeric-outward)
                                                      * (1-rd(model.Y.size+2)))
            ab[name][title] = receipt
            native_outputs[name+"_"+title+"_prediction"] = model.S@current
            independent_npz[name+"_"+title+"_current"] = current
    write_json(HERE/"B_independent_trajectory.json", dict(
        case="n32_c01", source="B_independent_control.npz", frequency_MHz=400,
        power_available_per_pulse_W=.125, source_column_power_sum_W=1.,
        pulse_dwell_s=None, reads_complex=64, pulses=8, radius_B=radius, radius_A=original.radius,
        B_layout=dict(transmit=tx.tolist(), receive=rx.tolist(), source_scaling=scale),
        material_box=[0., 4.5], physical_patch_cells=len(patch), patch_area_m2=area*len(patch),
        independent_endpoint_audits=ab, endpoint_pair_exterior_cost=compensation(high[~roi]-vacuum[~roi], area),
        endpoint_pair_region_means=[float(vacuum[roi].mean()), float(high[roi].mean())],
        endpoint_pair_patch_minima=[float(vacuum[patch].min()), float(high[patch].min())],
        fitted_pair_B_distance_over_two_r=float(np.linalg.norm(
            native_outputs["high_patch_B_prediction"] - native_outputs["vacuum_B_prediction"])/(2*radius)),
        reduction=reduced.receipt, exact_chart=exact, continuous=trajectory,
        reference_exterior_real_jacobian_singular_values=spectrum.tolist(),
        total=dict(wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1],
                   native_B_counts=dict(b.cost),native_A_counts=dict(original.cost))))
    write_arrays(HERE/"B_independent_trajectory.npz",
        B_observation=y, A_observation=original.Y, vacuum=vacuum, high_patch=high,
        verified_polyline=arr, Y_reduced=np.stack([k.prediction for k in knots]),
        patch=patch, roi=roi, support=support, K=reduced.kernel, D=reduced.transfer,
        e=reduced.field, F0=reduced.initial, **exact_positions, **native_outputs, **independent_npz)
    print("B-only path completed", trajectory["certified_patch_mean_interval"], flush=True)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
