"""Full-frequency reflection and high-dimensional archived-material controls."""
import json
from time import perf_counter, process_time

from geometry_experiment import Block, Fit, FullTM, ROOT, R12, SIGMA, complex_pairs, material_rms

import numpy as np
from scipy.stats import chi2, norm

from mirror_spectrum import CENTER, V, Mirror


def coherent_mirror_block(model, offset, energy):
    green,b,receiver,*_ = model.wave.kernels(offset)
    weights = model.weights(offset)[0]
    incident = np.einsum("igs,gs->ig",b.reshape(len(green),2,2),weights)*np.sqrt(energy)
    return Block(green,incident,receiver,np.r_[CENTER,4.5,4.5,6.2,6.2],np.arange(2))


def mirror_control():
    t0,cpu0 = perf_counter(),process_time()
    same = Mirror([0,0])
    phase = Mirror([1,2])
    old = coherent_mirror_block(same,0.,1.)
    # One unit of energy in the new three-exposure band.
    offsets = (-.08,.012,.08)
    energy = 1/3
    blocks = [old,*[coherent_mirror_block(same,d,energy) for d in offsets]]
    swap = [1,0,3,2,5,4]
    identities = [dict(green=float(np.max(abs(b.green[np.ix_(swap,swap)]-b.green))),
                       incident=float(np.max(abs(b.incident[swap]-b.incident))),
                       receiver=float(np.max(abs(b.receiver[:,swap]-b.receiver))))
                  for b in blocks]
    model = FullTM(blocks)
    material = CENTER+.6*V
    partner = material[::-1]
    centers = [model.evaluate(m)[0] for m in (material,partner)]
    rng = np.random.default_rng(14_250819)
    noise = SIGMA/np.sqrt(2)*(rng.normal(size=len(centers[0]))+1j*rng.normal(size=len(centers[0])))
    eps = SIGMA*np.sqrt(chi2.ppf(.99,2*len(noise))/2)
    fixed = centers[0]+noise
    outcomes = []
    fit = Fit(model,fixed,eps)
    for start in (material,partner):
        found, work = fit.solve(start)
        outcomes.append(dict(fit.diagnostics(found,material),cost=work))
    invariant = dict(
        old=0., added_offsets=offsets, energies=[energy]*3,
        permutation_identities=identities, rms_pair=material_rms(material,partner),
        pair_mean_gap=float(np.linalg.norm(centers[0]-centers[1])),
        signal_classification_risk=.5, scores_of_terms=[fit.score(x)[0] for x in (material,partner)],
        fixed_observed=complex_pairs(fixed), outcomes=outcomes, wave_cost=dict(model.counts))

    # Keep the same old design. One delayed new feed pair breaks symmetry but
    # the budget permits a perturbation along the mean difference that reverses
    # the ordering of the two physically distinct materials.
    offset = .00025
    weak = FullTM([old,coherent_mirror_block(phase,offset,1.)])
    material = CENTER+.05*V
    partner = material[::-1]
    signal, competing = [weak.evaluate(m)[0] for m in (material,partner)]
    old_channels = old.receiver.shape[0]*old.incident.shape[1]
    added_difference = competing[old_channels:]-signal[old_channels:]
    delta = np.r_[np.zeros(old_channels),added_difference]
    gap = np.linalg.norm(added_difference)
    noise_radius = SIGMA*np.sqrt(chi2.ppf(.99,2*len(added_difference))/2)
    rows = []
    for part in (0., .5, .6):
        y = signal+part*delta
        objective = Fit(weak,y,noise_radius)
        scores = [objective.score(x)[0] for x in (material,partner)]
        start_results = []
        for m in (material,partner):
            found, work = objective.solve(m,cap=180)
            start_results.append(dict(objective.diagnostics(found,material),cost=work))
        rows.append(dict(t=part,old_noise_over_epsilon=0., new_noise_over_epsilon=float(part*gap/noise_radius),
                         scores=scores,winner=int(np.argmin(scores)), observations=complex_pairs(y),
                         fits=start_results))
    pure_noise = dict(offset=offset, energy=1., rms_pair=material_rms(material,partner),
                      new_gap=float(gap), noise_radius=float(noise_radius), sigma_complex_rms=SIGMA,
                      conditional_two_hypothesis_gaussian_error=float(norm.cdf(-gap/(np.sqrt(2)*SIGMA))),
                      tests=rows,wave_cost=dict(weak.counts))
    return dict(reflection=invariant, additive_only=pure_noise,
                wall_seconds=perf_counter()-t0,cpu_seconds=process_time()-cpu0)


def high_archive():
    from freq_high import GAIN_SIGMA, HALF_POWER_AMPLITUDE, Model, OUT, PLANS, WIDE, forward_complex, old_candidates, power_operators
    from high_endpoints import noise_radius

    def score(y, observed, eps, gamma):
        lam = (eps/gamma)**2
        gain = (np.vdot(y,observed)+lam)/(np.vdot(y,y).real+lam)
        additive = np.linalg.norm(gain*y-observed)/eps
        nuisance = abs(gain-1)/gamma
        return dict(score=float(np.hypot(additive,nuisance)),additive=float(additive),
                    gain_delta_over_gamma=float(nuisance),gain=[gain.real,gain.imag])

    hessians = {x["case"]:x for x in json.loads((R12/"results/high_hessian.json").read_text())}
    table = []
    for name in ("ring","shell","bridge"):
        t0,cpu0 = perf_counter(),process_time()
        _,public,candidates = old_candidates(name)
        with np.load(OUT/f"{name}_wide_acquisition.npz",allow_pickle=False) as f:
            old,wide,new_noise = f["old"],f["wide"],f["new_noise"]
        with np.load(OUT/f"{name}_acquisition.npz",allow_pickle=False) as f:
            selected = json.loads((OUT/f"{name}_design.json").read_text())["selected"]
            selected_reads = f[selected]
        with np.load(OUT/f"{name}_wide_7200_fits.npz",allow_pickle=False) as f:
            direct,cumulative = f["wide_direct_bad"],f["wide_cumulative_bad"]
        scores = json.loads((OUT/f"{name}_wide_7200_scores.json").read_text())
        solves = json.loads((OUT/f"{name}_wide_7200_solves.json").read_text())
        sigma = float(public["sigma"])
        new_radius = sigma*np.sqrt(chi2.ppf(.99,2*new_noise.size)/2)
        add_radius = np.hypot((1+GAIN_SIGMA)*noise_radius(public),new_radius)
        epsilon = np.sqrt(2)*add_radius
        gamma = np.sqrt(2)*GAIN_SIGMA
        offsets = PLANS[selected]
        models = {d:Model(public if d==0. else power_operators(public,d)) for d in set((0.,*WIDE,*offsets))}
        def joint(material, ds, reads):
            output = {d:forward_complex(models[d],material) for d in {0.,*ds}}
            y = np.r_[output[0.],*(HALF_POWER_AMPLITUDE*output[d] for d in ds)]
            return score(y,np.r_[old,*reads],epsilon,gamma)
        selection = [joint(m,offsets,selected_reads) for m in candidates]
        path = {}
        for label, m in (("wide_direct_bad",direct),("wide_cumulative_bad",cumulative)):
            path[label] = dict(
                **joint(m,WIDE,wide), historical_material_NRMSE=scores[label],
                historical_RHS=solves[label]["counts"]["actual_complex_RHS"])
        old_h = hessians.get(name, {})
        table.append(dict(case=name, selected=selected, candidate_selected_read_compatibility=selection,
                          wide_same_data=path, epsilon=epsilon, gamma=gamma, pre_sqrt2_additive_radius=float(add_radius),
                          archived_old_noise_radius=noise_radius(public),
                          new_actual_noise=float(np.linalg.norm(new_noise)), new_nominal_99_noise=float(new_radius),
                          old_profiled_candidate_gap=json.loads((OUT/f"{name}_design.json").read_text())["old_profile_gap"],
                          old_regularized_bad_endpoint=dict(
                              free=old_h.get("free"),
                              projected_grad_inf=old_h.get("projected_gradient_inf"),
                              face_hessian=old_h.get("phi_hessian_on_free_box_face",{}).get("min"),
                              noise_multiple=old_h.get("model_data_residual_over_noise")),
                          new_analysis_factorizations=sum(m.counts["LU"] for m in models.values()),
                          new_analysis_rhs=sum(m.counts["actual_complex_RHS"] for m in models.values()),
                          wall_seconds=perf_counter()-t0,cpu_seconds=process_time()-cpu0))
    return table


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("control",choices=("mirror","high"))
    args = parser.parse_args()
    result = mirror_control() if args.control=="mirror" else high_archive()
    (ROOT/"results"/f"geometry-{args.control}.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if args.control=="mirror":
        print(result["reflection"]["rms_pair"],result["reflection"]["pair_mean_gap"],
              result["additive_only"]["conditional_two_hypothesis_gaussian_error"],
              [x["scores"] for x in result["additive_only"]["tests"]])
    else:
        print([(r["case"],[x["score"] for x in r["candidate_selected_read_compatibility"]],
                {k: v["score"] for k,v in r["wide_same_data"].items()}) for r in result])
