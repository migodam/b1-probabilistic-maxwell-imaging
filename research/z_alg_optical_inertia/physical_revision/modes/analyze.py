#!/usr/bin/env python3
"""Native 96² current modes and six-object, fixed-weight source diagnostics."""
from __future__ import annotations

from contextlib import contextmanager
import csv
import json
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

from flint import ctx
import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, orth, subspace_angles, svd
from scipy.sparse.linalg import LinearOperator, eigsh, gmres
from scipy.special import hankel1, j1
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ROOT = STUDY.parents[1]
CURVE = STUDY / "submission/curve"
EVALUATION = STUDY / "synthetic/evaluation"
TSOM = ROOT / (
    "research/delegated/max_struct_r1/r14_worktree/research/"
    "research_fbg_r3_experiments/src"
)
sys.path[:0] = [str(CURVE), str(TSOM)]
from cie_tsom_baseline import fourier_apcs
from curve_inputs import full_ring_inputs
from envelope import FourierEnvelope
from native import Geometry, Native, Weights

NAMES = ("Slotted block", "Four lobes", "Double chevron",
         "Crescent and plate", "Staggered grooves", "Graded pinwheel")
WIDTH = 2.0
K0 = 2 * np.pi * 400e6 / 3e8
RNG = np.random.default_rng(54273)


@contextmanager
def measure(receipt, key):
    wall, cpu = perf_counter(), process_time()
    yield
    receipt[key] = dict(wall_seconds=perf_counter() - wall,
                        cpu_seconds=process_time() - cpu)


def fractions(x):
    return (np.abs(x) ** 2).sum(axis=0)


def rayleigh(v, x):
    return (v.conj() * x).sum(axis=0).real


def hermitian(x):
    return (x + x.conj().T) / 2


def fix_phase(x):
    maximum = x[np.argmax(np.abs(x), axis=0), np.arange(x.shape[1])]
    return x * np.exp(-1j * np.angle(maximum))


def spectrum(operator, count, which):
    guess = RNG.standard_normal(operator.shape[0]) + 1j * RNG.standard_normal(operator.shape[0])
    values, vectors = eigsh(operator, k=count, which=which, ncv=max(32, 2 * count + 4),
                             tol=1e-10, maxiter=900, v0=guess)
    order = np.argsort(values)
    if which == "LA":
        order = order[::-1]
    vectors = fix_phase(vectors[:, order])
    values = values[order]
    return values, vectors, np.linalg.norm(operator @ vectors - vectors * values, axis=0)


def native_weights(record):
    return Weights(record["cap"], record["scale"], tuple(record["weights"]),
                   tuple(record["optical"]))


def green_audit(g):
    step = WIDTH / g.side
    radius = step / np.sqrt(np.pi)
    kr = K0 * radius
    c = np.pi * kr * j1(kr) / 2
    kernel = ifft2(g.spec)
    analytic = np.array([1j * np.pi * kr / 2 * hankel1(1, kr) - 1,
                         1j * c * hankel1(0, K0 * step),
                         1j * c * hankel1(0, K0 * step * np.sqrt(2))])
    stored = kernel[[0, 0, 1], [0, 1, 1]]
    bin_axis = (np.arange(g.side) + 0.5) * step - WIDTH / 2
    xx, yy = np.meshgrid(bin_axis, bin_axis)
    xy = np.column_stack((xx.ravel(), yy.ravel()))
    return xy, c, stored[0], dict(
        stored_self_plus_offsets=np.column_stack((stored.real, stored.imag)).tolist(),
        reproduced_outgoing=np.column_stack((analytic.real, analytic.imag)).tolist(),
        max_complex_error=float(np.max(np.abs(analytic - stored))),
        max_error_if_conjugated=float(np.max(np.abs(analytic.conj() - stored))),
        disk_radiation_coefficient=float(c),
    )


def angular_power(xy, c, currents, count=720):
    theta = np.arange(count) * (2 * np.pi / count)
    dhat = np.column_stack((np.cos(theta), np.sin(theta)))
    power = np.zeros((count, currents.shape[1]))
    for lo in range(0, count, 48):
        hi = min(count, lo + 48)
        phases = np.exp(-1j * K0 * (dhat[lo:hi] @ xy.T))
        power[lo:hi] = np.abs(phases @ currents)**2
    return c * power.mean(axis=0), power


def minimum_in_span(x):
    return float(svd(x, compute_uv=False)[-1])


def current_diagnostics(g, xy, c, receiver_v, x):
    projected = x - receiver_v @ (receiver_v.conj().T @ x)
    p = projected / np.sqrt(fractions(projected))
    gx = g.numerical_green(x)
    ghx = (gx + g.numerical_green(x, adjoint=True)) / 2
    gp = g.numerical_green(p)
    power, angular = angular_power(xy, c, x)
    projected_power, _ = angular_power(xy, c, p)
    incident = g.incident / np.linalg.norm(g.incident, axis=0)
    width = np.max(np.abs(xy), axis=1)
    outside_central_square = width >= .75
    return dict(
        ker_S_energy=fractions(projected).tolist(),
        receiver_norm=np.linalg.norm(g.sensor @ x, axis=0).tolist(),
        radiation_v_star_R_v=(x.conj() * gx).sum(axis=0).imag.tolist(),
        radiation_angular_integral=power.tolist(),
        ker_S_normalized_radiation_v_star_R_v=(p.conj() * gp).sum(axis=0).imag.tolist(),
        ker_S_normalized_radiation_angular_integral=projected_power.tolist(),
        reactive_v_star_H_v=rayleigh(x, ghx).tolist(),
        incident_energy_per_unit_plane_wave=(np.abs(x.conj().T @ incident) ** 2).tolist(),
        incident_energy_mean=np.mean(np.abs(x.conj().T @ incident)**2, axis=1).tolist(),
        outside_central_1p5m_square_energy=fractions(x[outside_central_square]).tolist(),
        inverse_participation_support_fraction=(1 / (g.n * (np.abs(x)**4).sum(axis=0))).tolist(),
        homogeneous_constitutive_defect=np.linalg.norm(x - .75 * gx, axis=0).tolist(),
        span_min_homogeneous_constitutive_defect=minimum_in_span(x - .75 * gx),
    ), angular


def fourier_diagnostics(g, x, receiver_v, receiver_v_cut):
    p = x - receiver_v @ (receiver_v.conj().T @ x)
    qdark = orth(p)
    dark_normalized = p / np.sqrt(fractions(p))
    fx = fft2(x.reshape(g.side, g.side, -1), axes=(0, 1), norm="ortho")
    bands = []
    for mf in (1, 2, 3, 4, 6, 10):
        zapcs, metadata = fourier_apcs((g.side, g.side), mf, receiver_v_cut)
        zker = zapcs - receiver_v @ (receiver_v.conj().T @ zapcs)
        qapcs, qker = orth(zapcs), orth(zker)
        bins = np.asarray(metadata["frequencies"])
        bands.append(dict(
            mf=mf, columns=metadata["fourier_columns"],
            truncated_dpcs_apcs_rank=qapcs.shape[1], exact_ker_fourier_rank=qker.shape[1],
            raw_fourier_energy=fractions(fx[bins[:, 0], bins[:, 1]]).tolist(),
            truncated_dpcs_apcs_energy=fractions(qapcs.conj().T @ x).tolist(),
            exact_ker_fourier_energy=fractions(qker.conj().T @ x).tolist(),
            exact_ker_fourier_of_normalized_dark_energy=fractions(qker.conj().T @ dark_normalized).tolist(),
            angles_mode_space_to_tsom_deg=np.sort(np.degrees(subspace_angles(x, qapcs))).tolist(),
            angles_dark_mode_space_to_projected_fourier_deg=np.sort(
                np.degrees(subspace_angles(qdark, qker))).tolist(),
        ))
    shifted = np.fft.fftshift(fx, axes=(0, 1))
    radial_q = np.hypot(*np.meshgrid(np.arange(-g.side // 2, g.side // 2),
                                        np.arange(-g.side // 2, g.side // 2)))
    return dict(bands=bands, circular_dft_bands=[
        dict(radius_cycles_per_domain=k, energy=fractions(shifted[radial_q <= k]).tolist())
        for k in (1, 2, 3, 4, 6, 10)
    ]), shifted


def solve_columns(operator, rhs, preconditioner=None, rtol=2e-12, restart=65, maxiter=10):
    result = np.empty_like(rhs)
    iterations = []
    for j, col in enumerate(rhs.T):
        trace = []
        result[:, j], info = gmres(operator, col, M=preconditioner,
                                   rtol=rtol, atol=0., restart=restart,
                                   maxiter=maxiter, callback=trace.append,
                                   callback_type="pr_norm")
        iterations.append(len(trace))
        if info:
            raise ArithmeticError(f"Native GMRES failed on column {j}: {info}")
    defects = np.linalg.norm(operator @ result - rhs, axis=0) / np.linalg.norm(rhs, axis=0)
    return result, iterations, defects


def fixed_source_profiles(native, obs_rows, spectrum_values, modes):
    g, mu = native.g, native.w.scale
    basis = np.column_stack((g.sensor_star / mu, native.z[:, None] * g.incident))
    pre = LinearOperator(native.operator.shape, native.precondition, dtype=complex)
    stationary, iterations, defects = solve_columns(native.operator, basis, pre)
    compressed = hermitian(basis.conj().T @ stationary)
    source_h, source_v, source_l, source_x = {}, {}, {}, {}
    source_stats = {}
    for obs in obs_rows:
        transform = np.vstack((obs.y, np.eye(g.s)))
        h = hermitian(obs.y.conj().T @ obs.y - mu * transform.conj().T @ compressed @ transform)
        values, vectors = eigh(h, check_finite=False)
        values, vectors = values[::-1], vectors[:, ::-1]
        coupling = modes.conj().T @ (basis @ transform)
        positive_dual_term = mu * (coupling.conj().T / (-spectrum_values)) @ coupling
        source_v[obs.id], source_l[obs.id], source_h[obs.id] = vectors, coupling, h
        source_x[obs.id] = stationary @ transform
        coupling_eigenbasis = fractions(coupling @ vectors)
        dist = float(np.sqrt(np.maximum(0, values[len(modes.T):]).sum()))
        source_stats[obs.id] = dict(
            padded_radius=obs.padded_radius,
            source_eigenvalues=values.tolist(),
            positive_source_count=int(np.count_nonzero(values > 0)),
            numerical_fixed_weight_distance=dist,
            numerical_fixed_weight_distance_over_radius=dist / obs.padded_radius,
            neg_state_mode_rhs_energy_in_descending_source_eigenbasis=coupling_eigenbasis.tolist(),
            negative_state_schur_term_in_source_eigenbasis=rayleigh(vectors, positive_dual_term @ vectors).tolist(),
            remaining_positive_state_schur_term_in_source_eigenbasis=rayleigh(
                vectors, (obs.y.conj().T @ obs.y + positive_dual_term - h) @ vectors
            ).tolist(),
            measured_energy_in_source_eigenbasis=fractions(obs.y @ vectors).tolist(),
            incident_component_rhs_energy_in_source_eigenbasis=fractions(
                modes.conj().T @ (native.z[:, None] * g.incident) @ vectors).tolist(),
        )
    return source_stats, source_h, source_v, source_l, source_x, dict(
        common_rhs_iterations=iterations, max_relative_stationarity_defect=float(defects.max()))


def physical_currents(g, chi, diagonal):
    active = np.flatnonzero(chi)
    local = chi[active]

    def extend(v):
        full = np.zeros((g.n, v.reshape(len(active), -1).shape[1]), complex)
        full[active] = v.reshape(len(active), -1)
        return full

    def apply(v):
        values = v.reshape(len(active), -1)
        result = values - local[:, None] * g.numerical_green(extend(values))[active]
        return result.reshape(v.shape)

    op = LinearOperator((len(active), len(active)), apply, matmat=apply, dtype=complex)
    pre = LinearOperator(op.shape, lambda v: v / (1 - local * diagonal), dtype=complex)
    rhs = local[:, None] * g.incident[active]
    x, iterations, defects = solve_columns(op, rhs, pre, rtol=2e-11, restart=120, maxiter=15)
    return extend(x), iterations, defects


def array_key(obs):
    return obs.id.split("__")[0][-2:] + f"_pct{obs.noise_pct:g}"


def verified_union(row):
    return max(float(row[key] or 0) for key in (
        "strict_global_distance_lower", "strict_local_distance_lower",
        "fixed_0_distance_lower", "fixed_2_distance_lower",
    ))


def main():
    timing = {}
    total_start = perf_counter(), process_time()
    with measure(timing, "setup"):
        spec, sensor, incident, observations = full_ring_inputs()
        g = Geometry(spec, sensor, incident)
        _, sv_sensor, vh_sensor = svd(sensor, full_matrices=False, check_finite=False)
        receiver_rank = int(np.count_nonzero(sv_sensor > max(sensor.shape) * np.finfo(float).eps * sv_sensor[0]))
        cutoff_rank = int(np.count_nonzero(sv_sensor > 1e-3 * sv_sensor[0]))
        receiver_v = vh_sensor[:receiver_rank].conj().T
        receiver_v_cut = vh_sensor[:cutoff_rank].conj().T
        local_ledger = json.loads((CURVE / "local.75.json").read_text())
        chosen = next(row for row in local_ledger["rows"]
                      if row["id"] == "optfresh_04__pw360_10x20__pair1pct")
        w = native_weights(chosen["profile_key"])
        historic = next(p for p in local_ledger["profiles"]
                        if native_weights(p) == w and p["status"] == "native_gap_and_inertia_verified")
        native = Native(g, w)
        xy, c, diagonal, outgoing = green_audit(g)
    with measure(timing, "native_modes"):
        eigenvalues, modes, residual_norms = spectrum(native.operator, 6, "SA")
        if not np.allclose(eigenvalues, historic["numerical_smallest_scaled_K"],
                           rtol=0, atol=1e-8):
            raise ArithmeticError("Pencil eigenvalues differ from the original fixed-profile receipt")
        nneg = int(np.count_nonzero(eigenvalues < 0))
        x = modes[:, :nneg]
        print("K/μ0 eigenvalues", eigenvalues, flush=True)
    with measure(timing, "mode_controls_and_fourier"):
        mode_stats, angular = current_diagnostics(g, xy, c, receiver_v, x)
        mode_stats["pencil_projection_K_over_scale"] = rayleigh(x, native.operator @ x).tolist()
        print("K/μ0 ker energy, power, homogeneous defect:",
              mode_stats["ker_S_energy"], mode_stats["radiation_v_star_R_v"],
              mode_stats["homogeneous_constitutive_defect"], flush=True)
        fourier, shifted_fft = fourier_diagnostics(g, x, receiver_v, receiver_v_cut)
        print("TSOM projected mode-space angles:", [(row["mf"], row["angles_mode_space_to_tsom_deg"])
                                                    for row in fourier["bands"]], flush=True)
        mean = float(np.mean(w.optical))
        control = Native(g, Weights(w.cap, w.scale, w.real, (mean,) * len(w.real)))
        global_eigs, global_modes, global_residuals = spectrum(control.operator, 6, "SA")
        default_global_modes = global_modes[:, global_eigs < 0]
        reactive = LinearOperator(native.operator.shape, lambda v:
            ((g.numerical_green(v[:, None]) + g.numerical_green(v[:, None], adjoint=True)) / 2)[:, 0],
            dtype=complex)
        reactive_eigs, reactive_modes, reactive_residuals = spectrum(reactive, 6, "LA")
        noise = np.random.default_rng(46730)
        random = noise.normal(size=(g.n, 4)) + 1j * noise.normal(size=(g.n, 4))
        random = orth(random - receiver_v @ (receiver_v.conj().T @ random))
        controls = {}
        for name, q in (("next_two_K_positive", modes[:, nneg:]),
                        ("largest_four_eigenvalues_of_H", reactive_modes[:, :4]),
                        ("global_equal_mean_optical_negative_space", default_global_modes),
                        ("random_rank_four_exact_dark", random)):
            measurements, _ = current_diagnostics(g, xy, c, receiver_v, q)
            controls[name] = dict(measurements, basis_overlap_angles_deg=np.sort(
                np.degrees(subspace_angles(x, q))).tolist())
        controls["global_equal_mean_optical_negative_space"].update(
            eigenvalues=global_eigs.tolist(), eigen_residuals=global_residuals.tolist(),
            uniform_optical_weight=mean,
        )
        controls["largest_four_eigenvalues_of_H"].update(
            largest_six_eigenvalues=reactive_eigs.tolist(), eigen_residuals=reactive_residuals.tolist())
        dark = x - receiver_v @ (receiver_v.conj().T @ x)
        dark_gram_eigenvalues = eigh(hermitian(dark.conj().T @ dark), eigvals_only=True)
    with measure(timing, "fixed_weight_source_spectra"):
        chosen_rows = [o for o in observations if o.noise_pct in (1, 2)]
        source_stats, source_h, source_v, source_l, source_x, solves = fixed_source_profiles(
            native, chosen_rows, eigenvalues[:nneg], x)
        selected = source_stats[chosen["id"]]
        if not np.allclose(selected["source_eigenvalues"], chosen["source_form_numeric_eigenvalues"],
                           rtol=0, atol=1e-9):
            raise ArithmeticError("Source spectrum differs from the original fixed-profile receipt")
        print("1% numerical fixed profile distance/noise:",
              [(o.id, source_stats[o.id]["numerical_fixed_weight_distance_over_radius"])
               for o in chosen_rows if o.noise_pct == 1], flush=True)
    with measure(timing, "outward_inertia_recheck"):
        # This preserves the exact archived FFT operator and archived binary64 profile.
        ctx.prec, ctx.threads = 112, 1
        native.eigs, native.eig_basis = eigenvalues, modes
        interval_negative, negative_gap = native.ensure_negative(nneg, return_gap=True)
        interval_positive = FourierEnvelope(g).verify(w, nneg)
        if interval_positive["upper_count"] != nneg:
            raise ArithmeticError("Projected upper and lower inertia counts differ")
    with measure(timing, "outward_fixed_weight_source_recheck"):
        for obs in chosen_rows:
            row = source_stats[obs.id]
            source_values = np.asarray(row["source_eigenvalues"])
            strict = native.ensure_source(
                obs.y, nneg, interval_positive["positive_gap_lower"], source_x[obs.id],
                source_v[obs.id], source_values, negative_gap=negative_gap,
            )
            if strict["distance_upper"] < strict["distance_lower"]:
                raise ArithmeticError("Fixed profile distance bounds are inconsistent")
            row["fixed_weight_distance_interval"] = [strict["distance_lower"], strict["distance_upper"]]
            row["fixed_weight_distance_interval_over_radius"] = [
                strict["distance_lower"] / obs.padded_radius, strict["distance_upper"] / obs.padded_radius
            ]
            row["fixed_weight_enclosed_positive_source_eigenvalues_lower"] = strict[
                "positive_source_eigenvalues_lower"
            ]
            row["fixed_weight_source_form_error_upper"] = strict["upper_form_perturbation"]
            print("Fixed weight interval", obs.id, row["fixed_weight_distance_interval_over_radius"], flush=True)
    with measure(timing, "uniform_medium_outgoing_control"):
        chi_uniform = np.full(g.n, w.cap)
        j_uniform, its, defects = physical_currents(g, chi_uniform, diagonal)
        rhs = w.cap * g.incident
        gains = np.sqrt(eigh(hermitian(j_uniform.conj().T @ j_uniform),
                             hermitian(rhs.conj().T @ rhs), eigvals_only=True))
        homogeneous = dict(
            max_relative_constitutive_residual=float(defects.max()), iterations=its,
            coherent_relative_gain=gains.tolist(),
            mean_plane_wave_gain=float(np.linalg.norm(j_uniform)/np.linalg.norm(w.cap*g.incident)),
            negative_K_projection_current_energy=float(np.linalg.norm(x.conj().T @ j_uniform)**2
                                                        /np.linalg.norm(j_uniform)**2),
        )
    raw, object_stats = {}, []
    with measure(timing, "offline_truth_after_source_spectra"):
        curve = {r["id"]: r for r in csv.DictReader((CURVE / "noise_curve.csv").open())
                 if float(r["cap"]) == w.cap}
        for index, name in enumerate(NAMES, 1):
            key = f"{index:02d}"
            object_id = f"optfresh_{key}"
            with np.load(EVALUATION / f"{object_id}.npz", allow_pickle=False) as f:
                chi, family = f["m_fine"], str(f["family"])
            clean_id = f"{object_id}__pw360_10x20__n02"
            with np.load(EVALUATION / f"{clean_id}.npz", allow_pickle=False) as f:
                clean = f["Y_fine_clean"]
            born = sensor @ (chi[:, None] * incident)
            current, its, defects = physical_currents(g, chi, diagonal)
            support = chi > 0
            clean_norm = np.linalg.norm(clean)
            replay_error = np.linalg.norm(sensor @ current - clean) / clean_norm
            if replay_error > 1e-8:
                raise ArithmeticError(f"Clean reference full-wave mismatch at object {key}")
            source = []
            for obs in chosen_rows:
                if obs.object_id != object_id:
                    continue
                row = curve[obs.id]
                verified_d = verified_union(row)
                d = np.linalg.norm(obs.y - clean)
                if float(row["radius"]) != obs.padded_radius:
                    raise ArithmeticError("Different noise radii in public input and strict ledger")
                if obs.noise_pct == 2:
                    with np.load(STUDY / "synthetic/public" / f"{clean_id}.npz", allow_pickle=False) as f:
                        if not np.array_equal(obs.y, f["Y"]):
                            raise ArithmeticError("At 2%, the paired observation must be the original input")
                source.append(dict(
                    id=obs.id, noise_pct=obs.noise_pct, realized_noise_over_clean=float(d / clean_norm),
                    union_verified_distance_lower=verified_d,
                    union_verified_distance_over_radius=verified_d / obs.padded_radius,
                    strictly_excludes=verified_d > obs.padded_radius,
                    verified_global_distance_lower=float(row["strict_global_distance_lower"] or 0),
                    verified_local_distance_lower=float(row["strict_local_distance_lower"] or 0),
                    verified_16_region_distance_lower=max(
                        float(row["fixed_0_distance_lower"] or 0), float(row["fixed_2_distance_lower"] or 0)),
                    **source_stats[obs.id],
                ))
                raw.update({
                    "source_matrix_" + array_key(obs): source_h[obs.id],
                    "source_eigenvectors_" + array_key(obs): source_v[obs.id],
                    "source_mode_rhs_" + array_key(obs): source_l[obs.id],
                })
            object_stats.append(dict(
                index=key, id=object_id, name=name, source_family=family,
                peak_contrast=float(chi.max()),
                dielectric_mass_m2=float(chi.sum() * (WIDTH/g.side)**2),
                support_area_m2=float(support.sum() * (WIDTH/g.side)**2),
                clean_scattered_field_frobenius=float(clean_norm),
                born_relative_error=float(np.linalg.norm(clean - born) / clean_norm),
                born_scattered_field_frobenius=float(np.linalg.norm(born)),
                mode_energy_inside_actual_support=fractions(x[support]).tolist(),
                true_current_energy_in_mode_span=float(
                    np.linalg.norm(x.conj().T @ current)**2 / np.linalg.norm(current)**2),
                true_current_max_relative_state_defect=float(defects.max()),
                clean_replay_relative_error=float(replay_error),
                true_current_gmres_iterations=its, source=source,
            ))
            raw.update({"chi_" + key: chi.reshape(g.side, g.side), "clean_field_" + key: clean,
                        "born_field_" + key: born, "physical_currents_" + key: current})
    with measure(timing, "write_raw_arrays"):
        np.savez_compressed(HERE/"plot_data.npz", mode_eigenvectors=modes, mode_eigenvalues=eigenvalues,
                            mode_dft=shifted_fft, mode_angular_power=angular,
                            dark_projected_mode_eigenvectors=dark, global_K_modes=global_modes,
                            upper_H_modes=reactive_modes, uniform_cap_currents=j_uniform, **raw)
    output = dict(
        contract=dict(geometry="PW360_10x20", grid_side=g.side, input_observations=12,
                      width_m=WIDTH, frequency_hz=400e6, wave_speed_m_s=3e8,
                      receivers=g.sensor.shape[0], illuminations=g.s,
                      receiver_singular_values=sv_sensor.tolist(),
                      receiver_numerically_resolved_rank=receiver_rank,
                      old_tsom_relative_dpcs_cutoff=1e-3, old_tsom_dpcs_rank=cutoff_rank,
                      local_profile_source=chosen["id"], local_profile=w.as_dict(),
                      exposure="Six fixed, exposed synthetic materials, 1% and 2% paired views",
                      source_work_before_offline_regeneration=True),
        outgoing_green_audit=outgoing, numerical_eigenvalues=eigenvalues.tolist(),
        numerical_eigen_defects=residual_norms.tolist(), numerical_negative_count=nneg,
        eigenvector_orthogonality_defect=float(np.linalg.norm(modes.conj().T @ modes - np.eye(len(eigenvalues)))),
        outward_inertia_recheck=dict(exact_negative_count=nneg, negative_block_upper=interval_negative,
                                    **interval_positive),
        original_native_interval_profile=historic,
        exact_dark_projection_subspace_energy_extrema=[
            float(dark_gram_eigenvalues[0]), float(dark_gram_eigenvalues[-1])],
        mode_metrics=mode_stats, mode_controls=controls,
        fourier=fourier, uniform_cap_physical_response=homogeneous,
        source_linear_system=solves, objects=object_stats, native_operations=dict(g.cost),
        timing=timing, wall_seconds=perf_counter()-total_start[0],
        cpu_seconds=process_time()-total_start[1], host=platform.platform(),
        numpy_version=np.__version__, python_version=platform.python_version(),
    )
    (HERE/"results.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print("Finished", output["wall_seconds"], output["cpu_seconds"], flush=True)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
