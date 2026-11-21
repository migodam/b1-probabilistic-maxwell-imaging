"""Check negative hidden directions in the one-scalar sector relaxation."""

import argparse
import json
import os
from pathlib import Path
from time import perf_counter, process_time

for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[var] = "1"

import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import LinearOperator, eigsh
from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent
U_GRID = (.1, .5, 1., 2., 3., 4.)
SENSITIVITY_RELATIVE_CUTOFF = 1e-8
EPS = np.finfo(float).eps


def norm(x):
    return float(np.linalg.norm(x))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def rayleigh(hermitian, row_vectors):
    """Dominant eigenpair of the orthogonally projected Hermitian matrix."""
    n = len(hermitian)

    def project(vector):
        return vector - row_vectors @ (row_vectors.conj().T @ vector)

    def matvec(vector):
        return project(hermitian @ project(vector))

    rng = np.random.default_rng(911_000+n)
    trial = rng.standard_normal(n)+1j*rng.standard_normal(n)
    start = project(trial)
    eigenvalue, basis = eigsh(LinearOperator((n, n), matvec=matvec, dtype=complex),
                             k=1, which="LA", ncv=32, tol=1e-12, v0=start)
    v = project(basis[:, 0])
    v /= norm(v)
    value = float(np.vdot(v, hermitian@v).real)
    return v, {"eigensolver_largest": float(eigenvalue[0]), "vG_Hv_real": value,
               "projected_eigen_residual_l2": norm(matvec(v)-value*v),
               "critical_material_u_if_positive": 1/value if value > 0 else None}


def inspect_receivers(g, receiver):
    U, s, Z = svd(receiver, full_matrices=False)
    numerical = max(receiver.shape)*EPS*s[0]
    cutoff = SENSITIVITY_RELATIVE_CUTOFF*s[0]
    rank = int(np.sum(s > numerical))
    rank_cut = int(np.sum(s > cutoff))
    left_rows = Z.conj().T
    v, original = rayleigh(g, left_rows)
    original["receiver_residual_l2"] = norm(receiver@v)
    original["receiver_residual_over_sigma_max"] = original["receiver_residual_l2"]/float(s[0])
    if rank_cut < len(s):
        vc, relaxed = rayleigh(g, left_rows[:, :rank_cut])
        relaxed["receiver_residual_l2"] = norm(receiver@vc)
        relaxed["receiver_residual_over_sigma_max"] = relaxed["receiver_residual_l2"]/float(s[0])
    else:
        relaxed = original.copy()
    return (U, s, Z), v, {
        "receivers": len(receiver), "known_vector_dimension": int(receiver.shape[1]),
        "sigma_max": float(s[0]), "sigma_min": float(s[-1]), "sigma_min_over_sigma_max": float(s[-1]/s[0]),
        "rank_at_maxdim_epsilon": rank, "absolute_maxdim_epsilon_cutoff": float(numerical),
        "rank_at_relative_1e-8": rank_cut, "rank_sensitive_to_cutoff": rank != rank_cut,
        "all_measured_rows_projected": original,
        "projected_against_tsvd_retained_rows": relaxed,
    }


def sector_q(green, incident, current, u):
    total = incident+green@current
    per_cell = np.sum(abs(current)**2-u*np.real(current.conj()*total), axis=1)
    return float(np.sum(per_cell)), total, per_cell


def exact_box_fit(current, field, u):
    level = np.sum(np.real(current.conj()*field), axis=1)
    normalizer = np.sum(abs(field)**2, axis=1)
    m = np.divide(level, normalizer, where=normalizer > 0, out=np.zeros_like(level))
    attainable_m = np.clip(m, 0, u)
    deviation = current-attainable_m[:, None]*field
    return {
        "one_real_material_per_cell": True,
        "best_constitutive_residual_frobenius": norm(deviation),
        "relative_to_current_norm": norm(deviation)/norm(current),
        "max_constitutive_cell_row_residual": float(np.max(np.linalg.norm(deviation, axis=1))),
    }


def exact_null_diagnostic(S, G, v):
    """Conservative double precision enclosure for the represented operators.

    Project v onto the exact kernel using the Gram inverse. The reported
    directional loss uses a Frobenius bound rather than equating a near-null
    vector with an exact null vector.
    """
    p, n = S.shape
    product = S@S.conj().T
    approximate = np.linalg.inv(product)
    # Standard dot-product rounding inflation, deliberately larger than the
    # complex BLAS factors, with padded Frobenius norms.
    gamma = 32*(n+p+1)*EPS
    S_bound, P_bound, A_bound = (2*norm(x) for x in (S, product, approximate))
    defect = (2*norm(np.eye(p)-approximate@product)
              + gamma*A_bound*(P_bound+S_bound**2))
    if not defect < 1:
        return {"status": "row_gram_inverse_not_certified", "gram_inverse_defect_upper": defect}
    sigma_lower = np.sqrt((1-defect)/A_bound)
    vector_upper = norm(v)*(1+gamma)
    near_null_upper = 2*norm(S@v)+gamma*S_bound*vector_upper
    projection_error = near_null_upper/sigma_lower
    G_bound = 2*norm(G)
    center = float(np.vdot(v, G@v).real)
    lower_before_unit = (center-gamma*G_bound*vector_upper**2
                         - G_bound*projection_error*(2*vector_upper+projection_error))
    # Orthogonal projection cannot increase ||v||, so this denominator suffices.
    normalized_lower = lower_before_unit/vector_upper**2 if lower_before_unit > 0 else None
    return {
        "status": "positive_rayleigh_lower_for_exact_null_projection" if normalized_lower else "inconclusive",
        "gram_inverse_defect_upper": defect, "verified_sigma_min_lower": float(sigma_lower),
        "computed_vector_receiver_norm_upper": near_null_upper,
        "exact_kernel_projection_displacement_upper": float(projection_error),
        "unit_exact_kernel_vG_Hv_lower": normalized_lower,
        "u_threshold_ensuring_negative_direction": None if normalized_lower is None else 1/normalized_lower,
    }


def witnesses(S, G, H, E, Y, delta, factor, v):
    u_svd, s, z = factor
    # Check a direction for every requested u above the same frozen eigenvalue.
    # The well conditioned receiver chosen before looking at Y needs no truncation.
    r = int(np.sum(s > s[0]*max(S.shape)*EPS))
    J0 = z[:r].conj().T @ ((u_svd[:, :r].conj().T @ Y)/s[:r, None])
    d0, v_s = norm(S@J0-Y), norm(S@v)
    lambda_v = float(np.vdot(v,H@v).real)
    cert = exact_null_diagnostic(S, G, v)
    gamma = 32*(S.shape[1]+1)*EPS
    S_bound, Y_bound, J0_bound = (2*norm(x) for x in (S, Y, J0))
    initial_fit_upper = 2*d0+gamma*(S_bound*J0_bound+Y_bound)
    records = []
    for u in U_GRID:
        lead = u*lambda_v-1
        q0, _, _ = sector_q(G,E,J0,u)
        rec = {"u": u, "u_lambda_minus_one": lead, "q_J0": q0}
        if lead <= 0:
            records.append({**rec, "status": "below_observed_hidden_threshold"})
            continue
        j = J0[:, 0]
        b = 2*np.vdot(v, j-u*(H@j))-u*np.vdot(v, E[:, 0])
        phase = -b/abs(b) if b else 1.
        h = phase*v
        t = np.sqrt(max(q0, 0.)/lead)
        current = J0.copy()
        current[:, 0] += t*h
        q, fields, qcells = sector_q(G,E,current,u)
        residual = norm(S@current-Y)
        h_bound = 2*norm(h)
        h_receive_upper = 2*norm(S@h)+gamma*S_bound*h_bound
        fit_upper = (initial_fit_upper+t*h_receive_upper
                     + gamma*S_bound*(J0_bound+t*h_bound))
        J_upper, G_upper, E_upper = (2*norm(x) for x in (current, G, E))
        q_upper = q+32*(current.size+len(G)+1)*EPS*(
            J_upper**2+u*J_upper*(E_upper+G_upper*J_upper)
        )
        exact_lower = cert.get("unit_exact_kernel_vG_Hv_lower")
        records.append({
            **rec, "status": "sector_and_observation_fit" if q <= 0 and residual <= delta else "numerical_witness_failed",
            "source_index_zero_based": 0, "phase_rotation_real_imag": [float(phase.real), float(phase.imag)],
            "quadratic_t": float(t), "unfavorable_linear_term_after_phase": float(np.real(np.conj(phase)*b)),
            "q_final_direct": q, "q_final_expansion": float(q0-t*abs(b)-lead*t*t),
            "l2_current_after_update": norm(current), "l2_current_change": float(t),
            "data_fit_residuals": residual, "data_fit_triangle_upper_float": float(d0+t*v_s),
            "rounding_guarded_q_upper": float(q_upper),
            "rounding_guarded_data_fit_triangle_upper": float(fit_upper),
            "archived_data_ball_margin": delta-residual, "archived_data_ball_margin_by_triangle_float": delta-d0-t*v_s,
            "rounding_guarded_archived_margin_lower": float(delta-fit_upper),
            "violating_pixel_sectors": int(np.sum(qcells > 0)),
            "maximum_pixel_sector_excess": float(np.max(qcells)),
            "nearest_material_box_constitutive_test": exact_box_fit(current, fields, u),
            "certified_negative_exact_kernel_direction": bool(exact_lower is not None and u*exact_lower > 1),
        })
    return {
        "material_array_opened": False, "initial_fit_kind": "minimum norm realisation of complex LS with ordinary SVD cutoff",
        "initial_rank_retained": r, "initial_current_l2": norm(J0), "initial_data_residual": d0,
        "original_radius": delta, "receiver_v_l2": v_s, "strong_array_verified_null_diagnostic": cert,
        "tests": records,
    }


def run(source, output):
    start, cpu = perf_counter(), process_time()
    rows, physical, cached_green, used_receiver = [], {}, {}, []
    eigen_count = 0
    with np.load(source/"public_vie_inputs.npz", allow_pickle=False) as archive:
        identifiers = [k.removesuffix("__G_key") for k in archive.files if k.endswith("__G_key")]
        for case in identifiers:
            G_tag = str(archive[case+"__G_key"].item())
            if G_tag not in cached_green:
                green = archive[G_tag]
                H = (green+green.conj().T)/2
                _, unobstructed = rayleigh(H, np.empty((len(green), 0),dtype=complex))
                eigen_count += 1
                cached_green[G_tag] = (green,H,unobstructed)
            green,H,unobstructed = cached_green[G_tag]
            S = archive[case+"__S"]
            physical[case] = (green,H,S,archive[case+"__E0"])
            report = {"id": case, "green_square_size": len(green), "sources": physical[case][3].shape[1],
                      "unconstrained_receiver_bound_u0": unobstructed, "receiver_layouts": []}
            for label, count in (("first_1",1),("first_half",len(S)//2),("all",len(S))):
                selection = S[:count]
                match = next((v for v in used_receiver
                              if v[0] == G_tag and np.array_equal(v[1], selection)), None)
                if match is None:
                    factor, v, meta = inspect_receivers(H, selection)
                    used_receiver.append((G_tag,selection,factor,v,meta))
                    eigen_count += 1+int(meta["rank_sensitive_to_cutoff"])
                else:
                    _,_,factor,v,meta = match
                stage = {"view": label, **meta}
                report["receiver_layouts"].append(stage)
                if label == "all":
                    data = physical[case]
                    physical[case] = (*data, factor, v, stage)
                print(case,label,stage["all_measured_rows_projected"]["critical_material_u_if_positive"],
                      stage["sigma_min_over_sigma_max"], flush=True)
            rows.append(report)

    # The array is chosen using public operators only. This rule favors all
    # measured rows being well constrained and requires a sampled u witness.
    qualified = [(case, v) for case,v in physical.items()
                 if v[-1]["all_measured_rows_projected"]["critical_material_u_if_positive"] is not None and
                    v[-1]["all_measured_rows_projected"]["critical_material_u_if_positive"] < max(U_GRID)]
    if qualified:
        case, (G,H,S,E,factor,v,summary) = max(qualified, key=lambda item: item[1][-1]["sigma_min_over_sigma_max"])
        with np.load(source/"results"/f"{case}.npz", allow_pickle=False) as fields:
            Y = fields["Y_archive"]
        radius = read_json(source / "observation_radii.json")[case]["frozen_observation_ball"]
        witness = {"chosen_id": case, **witnesses(S,G,H,E,Y,radius,factor,v)}
    else:
        witness = {"status": "no_hidden_threshold_below_tested_material_caps"}
    write = {
        "model": "one global summed quadratic sector; pointwise material law is a different constraint",
        "receiver_order": "exact archived row order, first row, first floor(m/2), complete",
        "constraint_scope": "arithmetic on frozen operators, no new data or full-wave solve",
        "u_grid": U_GRID, "rows": rows, "example": witness,
        "example_selection_rule": "Among full receiver arrays with an eligible fixed u, use largest sigma_min/sigma_max; no Y before selection.",
        "geometry_eigenproblems": eigen_count, "unique_receiver_svd_problems": len(used_receiver),
        "material_arrays_read": 0, "maxwell_state_solves": 0,
        "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
    }
    output.write_text(json.dumps(write, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("witness",witness.get("chosen_id"),[(x["u"],x["status"]) for x in witness.get("tests",[])],flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=SOURCE)
    p.add_argument("--output", type=Path, default=HERE/"results.json")
    opts = p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        run(opts.source, opts.output)
