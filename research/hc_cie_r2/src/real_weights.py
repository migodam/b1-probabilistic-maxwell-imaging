#!/usr/bin/env python3
"""Audit the real-material weighted optical equality on four frozen observations."""
import argparse
import json
import os
from pathlib import Path
import platform
import resource
from time import perf_counter, process_time

for option in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[option] = "1"

import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve, solve, svd, svdvals

from spectra_inputs import load_observations


CASES = {
    "r13_r13_fresh_0000": "24² half ring",
    "r13_r13_fresh_0001": "24² full ring",
    "srap_s2d_010__pw360_10x20": "32² full ring",
    "srap_s2d_010__ls_line_same_10x20": "32² line",
}
SEED = 402909
CUTS = [1e-2, 1e-4, 1e-6, 1e-10]
EPS = np.finfo(float).eps


def norm(a):
    return float(np.linalg.norm(a))


def imh(a):
    return (a-a.conj().T)/(2j)


def arithmetic_gamma(operations):
    a = (8*operations+8)*EPS
    return a/(1-a)


def matrix_h(w, n, gn):
    return imh(n.conj().T@(w[:, None]*gn))


def gram_h(n, gn):
    p, k = n@n.conj().T, gn@gn.conj().T
    c = gn@n.conj().T
    g = .5*(k*p.conj()-c*c.T).real
    return (g+g.T)/2, p, k, c


def gram_formation_radius(g, n, gn, p, k, c):
    dim, dark_dim = n.shape
    rn, rv = np.linalg.norm(n, axis=1), np.linalg.norm(gn, axis=1)
    ev = arithmetic_gamma(dim)*np.linalg.norm(g, axis=1)*np.linalg.norm(n)
    ep = arithmetic_gamma(dark_dim)*np.outer(rn, rn)
    ek = (arithmetic_gamma(dark_dim)*np.outer(rv, rv)
          +np.outer(ev,rv)+np.outer(rv,ev)+np.outer(ev,ev))
    ec = arithmetic_gamma(dark_dim)*np.outer(rv,rn)+np.outer(ev,rn)
    err = (.5*(abs(p)*ek+abs(k)*ep+ep*ek+abs(c)*ec.T+abs(c.T)*ec+ec*ec.T)
           +arithmetic_gamma(1)*(abs(p*k)+abs(c*c.T)))
    return norm(err)


def diagonal_optical(j, gj, w):
    return imh(j.conj().T@(w[:, None]*gj))


def linear_signature(obs, n, gn, j, w):
    # IH(Z* A + C* Z) = IH(Z* (A-C)).
    e = obs.b+obs.g@j
    t = n.conj().T@(w[:, None]*e)-gn.conj().T@(w[:, None]*j)
    r = svdvals(t, check_finite=False)
    left, right = np.triu_indices(t.shape[1], k=1)
    values = np.sort(np.r_[r, np.repeat(np.hypot(r[left], r[right])/np.sqrt(2), 2)])[::-1]
    return t, dict(
        singular_values=values.tolist(), s_squared=t.shape[1]**2,
        minimum=float(values[-1]), maximum=float(values[0]),
        largest_over_source_fro=float(values[0]/norm(obs.b)),
        full_column_nullity_at_least=2*t.size-t.shape[1]**2,
        relative_rank={str(cut):int(np.sum(values>cut*values[0])) for cut in CUTS},
    )


def physical_current(obs, chi):
    a = -chi[:, None]*obs.g
    a.flat[::len(a)+1] += 1
    lu = lu_factor(a, overwrite_a=True, check_finite=False)
    return lu_solve(lu, chi[:, None]*obs.b, check_finite=False)


def weak_weights(g, n, gn):
    gram, p, k, c = gram_h(n, gn)
    size = len(n)
    one = np.ones(size)/np.sqrt(size)
    ihg = imh(g)
    h_const = n.conj().T@(ihg@n)/np.sqrt(size)
    const_norm = norm(h_const)
    h_const_scaled = h_const/const_norm
    # Direct trace cross-column; constant block is formed without subtracting O(1) Gram terms.
    const_cross = np.sum((gn@h_const_scaled)*n.conj(), axis=1).imag
    householder = np.eye(size)
    direction = householder[:, -1]-one
    householder -= 2*np.outer(direction,direction)/(direction@direction)
    perp = householder[:, :-1]
    restricted = perp.T@gram@perp
    restricted = (restricted+restricted.T)/2
    maximum_norm_scale = np.sqrt(max(np.max(np.diag(restricted)), EPS))
    complement_scaled = restricted/maximum_norm_scale**2
    cross = (perp.T@const_cross)/maximum_norm_scale
    factored = np.empty_like(gram)
    factored[:-1, :-1] = complement_scaled
    factored[:-1, -1] = cross
    factored[-1, :-1] = cross
    factored[-1, -1] = 1.
    scaled_values = eigh(factored, eigvals_only=True, check_finite=False)

    # Reconstruct the full Gram from stable constant-row traces before diagonalizing.
    scaled = np.r_[np.full(size-1, maximum_norm_scale), const_norm]
    stable_gram = scaled[:, None]*factored*scaled
    values, vectors = eigh(stable_gram, check_finite=False)
    weights = householder@vectors[:, 0]
    weights *= np.sign(one@weights)
    smallest = matrix_h(weights, n, gn)
    direct_min = norm(smallest)
    direct_min_norm_operator = float(np.max(abs(eigh(smallest, eigvals_only=True, check_finite=False))))
    gram_scale = gram_formation_radius(g, n, gn, p, k, c)
    gram_values = eigh(gram, eigvals_only=True, check_finite=False)
    rng = np.random.default_rng(SEED)
    w1, w2 = rng.standard_normal((2,size))
    w1 /= norm(w1)
    w2 /= norm(w2)
    h1, h2 = matrix_h(w1,n,gn), matrix_h(w2,n,gn)
    cross_trace = float(np.vdot(h1,h2).real)
    strongest = householder@vectors[:, -1]
    stable_sig = np.sqrt(np.maximum(values,0))
    eigen_residual = norm(stable_gram@vectors-vectors*values)
    gram_eigen_radius = gram_scale+norm(gram)*arithmetic_gamma(size)*4+eigen_residual
    eigen_positive = int(np.sum(gram_values>gram_eigen_radius))
    scaled_solution = solve(complement_scaled, cross, assume_a="pos", check_finite=False)
    correction = (const_norm/maximum_norm_scale)*scaled_solution
    schur_relative = 1-cross@scaled_solution
    schur_sigma = const_norm*np.sqrt(max(schur_relative,0.))
    schur_weight = one-perp@correction
    schur_direct = (h_const+matrix_h(-perp@correction,n,gn))/norm(schur_weight)
    # The reconstructed small squared eigenvalue is untrusted at the Gram error floor.
    # Publish the independently evaluated candidate as an upper witness, not a positive certificate.
    if eigen_positive == size-1:
        weights = schur_weight/norm(schur_weight)
        smallest = schur_direct
        direct_min = norm(smallest)
        direct_min_norm_operator = float(np.max(abs(eigh(smallest, eigvals_only=True,
                                                          check_finite=False))))
        stable_sig[0] = direct_min
    diagnostics = dict(
        normal_operator_rows=size, hermitian_output_real_dims=n.shape[1]**2,
        gram_formation_roundoff_fro=gram_scale,
        conservative_dense_gram_error_scale_fro=gram_eigen_radius,
        dense_gram_resolved_positive_count=eigen_positive,
        minimum_reported_is_only_direct_upper_witness=eigen_positive != size,
        numerically_unresolved_weight_count=size-eigen_positive,
        schur_congruence_complement_correction_norm=norm(correction),
        schur_complement_over_constant_norm_squared=float(schur_relative),
        schur_smallest_sigma_estimate=float(schur_sigma),
        schur_candidate_direct_norm=norm(schur_direct),
        direct_min_without_constant_split=norm(matrix_h(weights,n,gn)),
        stable_gram_eigendecomposition_residual_fro=eigen_residual,
        raw_gram_smallest_signed_eig=float(gram_values[0]),
        stabilized_gram_smallest_signed_eig=float(values[0]),
        stable_gram_eigenvectors_orthogonality_fro=norm(vectors.T@vectors-np.eye(size)),
        independently_scaled_gram_min=float(scaled_values[0]),
        independently_scaled_gram_max=float(scaled_values[-1]),
        constant_action_normalized=const_norm,
        constant_action_unscaled=const_norm*np.sqrt(size),
        constant_direct_scaled_vs_weighted_matrix_error=norm(h_const-matrix_h(one,n,gn)),
        direct_minimum_weight_norm=direct_min,
        direct_minimum_weight_spectral_norm=direct_min_norm_operator,
        direct_largest_weight_norm=norm(matrix_h(strongest,n,gn)),
        weight_minimum_cosine_with_constant=float(one@weights),
        weight_nonconstant_norm=norm(weights-one*(one@weights)),
        gram_trace_test=dict(direct=cross_trace, hadamard=float(w1@gram@w2),
                             absolute_error=float(abs(cross_trace-w1@gram@w2))),
        relative_effective_directions={str(cut): int(np.sum(stable_sig<=cut*stable_sig[-1]))
                                       for cut in CUTS},
        first_five_singular_values=stable_sig[:5].tolist(),
        last_five_singular_values=stable_sig[-5:].tolist(),
    )
    return weights, one, stable_sig, values, gram_values, scaled_values, diagnostics


def audit_state(obs,n,gn,w,material,j,port_anchor,delta_j):
    j_0 = j-n@(n.conj().T@j)
    dark_j = j-j_0
    dark_trial = n@(n.conj().T@delta_j)
    trial_green = obs.g@dark_trial
    d = n.conj().T@dark_j
    h = n.conj().T@dark_trial
    signature_affine, spectrum_affine = linear_signature(obs,n,gn,j_0,w)
    signature_material, spectrum_material = linear_signature(obs,n,gn,j,w)
    q = diagonal_optical(dark_j,obs.g@dark_j,w)
    l = imh(d.conj().T@signature_affine)
    r_0 = diagonal_optical(j_0,obs.b+obs.g@j_0,w)
    l_probe = imh(h.conj().T@signature_material)
    q_probe = diagonal_optical(dark_trial,trial_green,w)
    r_new = diagonal_optical(j+dark_trial,obs.b+obs.g@(j+dark_trial),w)
    r_phys = diagonal_optical(j,obs.b+obs.g@j,w)
    return dict(
        data_port_over_delta=float(port_anchor/obs.delta),
        material_range=[float(material.min()),float(material.max())],
        state_current_norm=norm(j), physically_data_feasible=bool(port_anchor<=obs.delta),
        baseline_optical_norm=norm(r_phys),
        linear_spectrum_on_same_physical_data_affine_representative=spectrum_affine,
        linear_spectrum_at_physical_material=spectrum_material,
        hidden_current_norm=norm(dark_j), affine_offset_optical=norm(r_0),
        finite_hidden_linear=norm(l), finite_hidden_quadratic=norm(q),
        quadratic_over_linear=float(norm(q)/norm(l)),
        affine_decomposition_absolute_roundoff=norm(r_0+l+q-r_phys),
        material_plus_two_percent_diff_norm=norm(delta_j),
        material_plus_two_percent_dark_diff_norm=norm(dark_trial),
        plus_two_percent_dark_linear=norm(l_probe),
        plus_two_percent_dark_quadratic=norm(q_probe),
        plus_two_percent_q_over_l=float(norm(q_probe)/norm(l_probe)),
        plus_two_percent_decomposition_absolute_roundoff=norm(r_new-r_phys-l_probe-q_probe),
    )


def run_observation(obs, refined_seed):
    start, cpu_start = perf_counter(), process_time()
    _, receive, vh = svd(obs.s, full_matrices=True, check_finite=False)
    rel_cut = max(obs.s.shape)*EPS
    rank = int(np.sum(receive>rel_cut*receive[0]))
    n = vh[rank:].conj().T
    gn = obs.g@n
    weights, one, stable_sig, gram_vals, raw_gram_vals, scaled_vals, report = weak_weights(obs.g,n,gn)
    comp_cost = [perf_counter()-start, process_time()-cpu_start]
    if refined_seed is None:
        base, source = obs.base, obs.provenance["base"]
    else:
        base, source = refined_seed
    j = physical_current(obs,base)
    trial_material = np.clip(1.02*base,obs.lower,obs.upper)
    jp = physical_current(obs,trial_material)
    port = norm(obs.s@j-obs.y)
    source_norm = norm(obs.b)
    constant_dark_source = norm(n.conj().T@obs.b)/np.sqrt(len(n))
    constant_dark_radiation = norm(n.conj().T@imh(obs.g))/np.sqrt(len(n))
    physical = {name:audit_state(obs,n,gn,w,base,j,port,jp-j)
                for name,w in (("minimum_weight",weights),("constant_weight",one))}
    delta_data = norm(obs.s@jp-obs.y)/obs.delta
    report.update(
        observation=obs.key, array=CASES[obs.key], material_id=obs.physical_id,
        source_count=obs.b.shape[1], receiver_count=obs.s.shape[0],
        dark_basis_dims=list(n.shape), full_s_receiver_singular_values=receive.tolist(),
        numerical_receive_relative_cut=rel_cut, receive_rank=rank,
        relative_last_receive=float(receive[-1]/receive[0]),
        receive_null_fro=norm(obs.s@n), dark_iso_error=norm(n.conj().T@n-np.eye(n.shape[1])),
        green_diagonal_real_min=float(obs.g.diagonal().real.min()),
        green_diagonal_imag_min=float(obs.g.diagonal().imag.min()),
        dark_normalized_incident_over_B_fro=constant_dark_source/source_norm,
        dark_normalized_radiation_fro=constant_dark_radiation,
        constant_linear_signal_upper_on_public_probe=constant_dark_source+
            2*constant_dark_radiation*norm(j),
        provenance=obs.provenance, physical_state_source=source,
        physical_candidate_plus_2percent_port_over_delta=delta_data, weights=physical,
        spectrum_wall_seconds=comp_cost[0], spectrum_cpu_seconds=comp_cost[1],
        forward_solves=2, forward_rhs=2*obs.b.shape[1],
        observation_wall_seconds=perf_counter()-start, observation_cpu_seconds=process_time()-cpu_start,
    )
    arrays = dict(curvature_singular_values_ascending=stable_sig, minimum_real_weight=weights,
                  stabilized_signed_gram=gram_vals, ordinary_signed_gram=raw_gram_vals,
                  scaled_signed_gram=scaled_vals, constant_weight=one, anchor_material=base,
                  perturbed_material=trial_material)
    return report,arrays


def plot(cases, arrays, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    colors = ("#316581", "#407764", "#b3563e", "#674c8b")
    fig, axes = plt.subplots(1,2,figsize=(11,4.2),layout="constrained")
    for i,(record, bundle) in enumerate(zip(cases,arrays)):
        sing = bundle["curvature_singular_values_ascending"]
        rank = np.arange(len(sing))/len(sing)
        axes[0].semilogy(rank,sing/sing[-1],color=colors[i],label=record["array"],lw=1.7)
        ctr = record["constant_action_normalized"]/sing[-1]
        axes[0].scatter(0,ctr,color=colors[i],marker="x")
        weights = record["weights"]["minimum_weight"]
        axes[1].scatter(i-.08,weights["quadratic_over_linear"],color=colors[i],
                        label="Finite physical dark current" if i == 0 else None,marker="o")
        axes[1].scatter(i+.08,weights["plus_two_percent_q_over_l"],color=colors[i],
                        label="+2% material perturbation, dark projection" if i == 0 else None,marker="s")
    axes[0].set(xlabel="Fraction of ascending weight spectrum",
                ylabel="Curvature singular value / maximum",title="Complete real-weight curvature spectrum")
    axes[0].xaxis.set_major_formatter(PercentFormatter(1))
    axes[0].legend(frameon=False,fontsize=9)
    axes[1].set_xticks(np.arange(len(cases)),[x["array"] for x in cases],rotation=17,ha="right")
    axes[1].set_yscale("log")
    axes[1].set(ylabel="Quadratic norm / linear norm", title="Weakest real weight · physical state probes")
    axes[1].axhline(1,ls=":",color="#777")
    axes[1].legend(frameon=False,fontsize=8)
    for ax in axes:
        ax.grid(alpha=.22)
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(output,dpi=160)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace",type=Path,required=True)
    p.add_argument("--out",type=Path)
    p.add_argument("--figure",type=Path)
    p.add_argument("--saved-physical",type=Path,help="path to public feasible_followup.npz")
    args = p.parse_args()
    project = Path(__file__).resolve().parents[1]
    out = args.out or project/"results/real_weights"
    out.mkdir(parents=True,exist_ok=True)
    figure = args.figure or project/"figures/real_weights.png"
    figure.parent.mkdir(parents=True,exist_ok=True)
    clock, cpu = perf_counter(), process_time()
    spectrum = json.loads((project/"results/spectra/feasible_followup.json").read_text())
    feasible_rows = spectrum["contrast"]
    chosen = min(feasible_rows,key=lambda row:row["physical_port_over_delta"])
    assert chosen["scale"] == 1 and chosen["physical_port_over_delta"] <= 1
    saved = args.saved_physical or project/"results/spectra/feasible_followup.npz"
    with np.load(saved,allow_pickle=False) as physical:
        physical_material = physical["base_material"].copy()

    design = {"r13_cases":["r13_fresh_0000","r13_fresh_0001"],"r15_cases":[],
              "r16_cases":{},"srap_objects":["s2d_010"],"srap_layouts":["pw360","ls_line_same"]}
    reports,arrays = [],[]
    for obs in load_observations(args.workspace.resolve(),design):
        alternative = None
        if obs.key == spectrum["observation"]:
            alternative = physical_material,"public spectra feasible physical follow-up at scale=1"
        elif obs.cohort == "r13":
            study = Path(obs.provenance["observation"]).parents[2]
            corrected = study/"RESULTS/confirmation_core48"/(obs.physical_id+"__cie")/"corrected.npz"
            with np.load(args.workspace/corrected,allow_pickle=False) as physical:
                alternative = physical["m"].ravel().copy(),str(corrected)
        report,bundle = run_observation(obs,alternative)
        np.savez_compressed(out/(obs.key+".npz"),**bundle)
        reports.append(report)
        arrays.append(bundle)
        print(obs.key,report["direct_minimum_weight_norm"],report["constant_action_normalized"],
              report["weights"]["minimum_weight"]["data_port_over_delta"],flush=True)
    plot(reports,arrays,figure)
    meta = dict(python=platform.python_version(),numpy=np.__version__,platform=platform.platform(),
                blas_threads=1,rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                wall_seconds=perf_counter()-clock,cpu_seconds=process_time()-cpu,
                reports=reports)
    (out/"receipt.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
