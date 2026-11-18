"""Conditional paired frequency tests, local jets and true stationary continuation."""

from argparse import ArgumentParser
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.linalg import null_space
from scipy.optimize import brentq, root, root_scalar
from scipy.stats import qmc

from freq_model import BOX, Fit, FrequencyTM, ri

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "freq_r11_input.json"
R10 = ROOT.parent / "research_z_alg_r10" / "results"

# Magnitudes are declared before creating any new read.
MENU = np.array([.00025, .0005, .001, .002, .003, .004, .006, .008, .012, .02, .04, .08])
RNG_SEED = 20261119
ROBUST_SIGMA = 3e-7
GAIN_SD = 4e-4


def material_error(material, truth):
    return float(np.linalg.norm(material - truth) / np.sqrt(truth.size))


def corrected_pair(model, old, first, second):
    fit = Fit(model, [0.], [old])
    return np.array([fit.solve(x)[0] for x in [first, second]])


def candidate_basis(model, pair):
    middle = pair.mean(axis=0)
    a = np.linalg.norm(pair[1] - pair[0]) / 2
    v = (pair[1] - pair[0]) / (2 * a)
    w = null_space(v[None])
    j = ri(model.response(middle)[1])
    u = np.linalg.svd(j @ w, full_matrices=True)[0][:, :w.shape[1]]
    reference = model.response(pair[0])[0]

    def center_residual(z):
        m = middle + v*z[0] + w @ z[1:]
        y, j = model.response(m)
        jac = ri(j)
        return np.r_[u.T @ ri(y-reference), np.linalg.det(jac) /
                     np.prod(np.linalg.svd(jac @ w, compute_uv=False))]

    correction = root(center_residual, np.zeros(middle.size), tol=1e-8)
    center = middle + v*correction.x[0] + w @ correction.x[1:]
    if np.linalg.norm(correction.fun) < 1e-10 and np.linalg.norm(correction.x) < a/2:
        u, _, vh = np.linalg.svd(ri(model.response(center)[1]))
        v_new = vh[-1] * np.sign(vh[-1] @ v)
        return center, v_new, null_space(v_new[None]), a, True
    return middle, v, w, a, False


def predicted_profile(model, old, candidate, basis, gain_calibration):
    center, v, w, a, _ = basis
    y, j, h = model.response(center, hessian=True)
    yp, jd = model.response(center, derivative=True)[2:]
    ypc = model.response(candidate, derivative=True)[2]
    b = np.einsum("kij,i,j->k", h, a*v, a*v) / 2
    linear = j @ (a*v)
    extra_linear = jd @ (a*v)
    n = len(y)
    if gain_calibration:
        lam, g0 = gain_calibration
        # The conditional reference assumes data at nominal calibrated gain.
        nuisance = np.r_[y, y, np.sqrt(lam)]
        derivatives = np.column_stack([np.vstack([j@w, j@w, np.zeros(w.shape[1])]),
                                       nuisance, 1j*nuisance])
        residual = np.r_[y-old, y-old, 0.]
        quad = np.r_[b, b, 0.]
        tangent = np.r_[linear, linear, 0.]
        extra_constant = np.r_[np.zeros(n), yp-ypc, 0.]
        extra_tangent = np.r_[np.zeros(n), extra_linear, 0.]
    else:
        derivatives = np.vstack([j@w, j@w])
        residual = np.r_[y-old, y-old]
        quad = np.r_[b, b]
        tangent = np.r_[linear, linear]
        extra_constant = np.r_[np.zeros(n), yp-ypc]
        extra_tangent = np.r_[np.zeros(n), extra_linear]
    u = np.linalg.svd(ri(derivatives), full_matrices=True)[0][:, derivatives.shape[1]:]
    constant, tangent, quad, cd, td = [u.T @ ri(x) for x in (
        residual, tangent, quad, extra_constant, extra_tangent)]
    scale = np.linalg.norm(quad)
    constant, tangent, quad, cd, td = [x/scale for x in (constant, tangent, quad, cd, td)]

    def polynomial(delta):
        coeffs = np.vstack([constant + delta*cd, tangent + delta*td, quad])
        return sum((Polynomial(c)**2 for c in coeffs.T), Polynomial([0.])) / 2

    def discriminant(delta):
        a0, a1, a2, a3 = polynomial(delta).deriv().coef
        return 18*a3*a2*a1*a0 - 4*a2**3*a0 + a2*a2*a1*a1 - 4*a3*a1**3 - 27*a3*a3*a0*a0

    thresholds = {}
    for sign in (-1, 1):
        ds = np.r_[0., MENU] * sign
        brackets = [(di, dj) for di, dj in zip(ds[:-1], ds[1:]) if
                    discriminant(di) > 0 and discriminant(dj) <= 0]
        thresholds[str(sign)] = float(brentq(discriminant, *sorted(brackets[0]), xtol=1e-13)) if brackets else None
    return polynomial, thresholds


def exact_profile(fit, basis, t):
    center, v, w, a, _ = basis
    offset = center + (a*t)*v

    def stability(u):
        grad, hess = fit.evaluate(offset + w @ u, True)[2:4]
        return w.T @ grad, w.T @ hess @ w

    sol = root(stability, np.zeros(w.shape[1]), jac=True, tol=2e-10)
    m = offset + w @ sol.x
    residual, _, grad, hess, _ = fit.evaluate(m, True)
    h_strong = w.T @ hess @ w
    cross = v @ hess @ w
    g = a * v @ grad
    h = a*a * (v @ hess @ v - cross @ np.linalg.solve(h_strong, cross))
    return m, g, h, float(np.linalg.norm(w.T@grad)), float(np.linalg.eigvalsh(h_strong)[0])


def real_roots_of_polynomial(poly):
    return sorted(float(z.real) for z in poly.deriv().roots() if abs(z.imag) < 1e-8)


def point(fit, basis, t, truth):
    m, g, h, strong_gradient, strong_min = exact_profile(fit, basis, t)
    record = fit.describe(m)
    record.update(t=float(t), material=m.tolist(), material_RMS=material_error(m, truth),
                  profile_gradient=float(g), profile_hessian=float(h),
                  strong_gradient=strong_gradient, strong_hessian_min=strong_min)
    return record


def compare_folds(model, basis, pair, old, read, truth, calibration, polynomial, predicted):
    center, v, _, a, _ = basis
    target_x = (pair[1]-center)@v/a
    wrong_x = (pair[0]-center)@v/a
    zero = Fit(model, [0., 0.], [old, read(0.)], gain_calibration=calibration)
    gradient_scale = max(abs(exact_profile(zero, basis, wrong_x+.25)[1]), 1e-18)
    descriptions = {}
    for sign in (-1, 1):
        pred = predicted[str(sign)]
        grid = np.r_[0., MENU[:8]] * sign
        if pred is not None:
            grid = np.unique(np.r_[grid, np.linspace(0, pred, 12, endpoint=False), pred * .999])
        grid = sorted(grid, key=abs)
        previous = [wrong_x, (wrong_x+target_x)/2]
        trajectory = []
        for delta in grid:
            fit = Fit(model, [0., delta], [old, read(delta)], gain_calibration=calibration)
            pts = []
            for t in previous:
                def f(t):
                    record = exact_profile(fit, basis, t)
                    return record[1] / gradient_scale, record[2] / gradient_scale

                sol = root_scalar(f, x0=t, method="newton", fprime=True, maxiter=20, xtol=1e-9)
                state = point(fit, basis, sol.root, truth)
                if (abs(state["profile_gradient"]) > gradient_scale*1e-5 or
                        state["strong_gradient"] > 1e-11 or state["strong_hessian_min"] <= 0 or
                        state["active"] or np.any(np.asarray(state["material"]) <= BOX[0]) or
                        np.any(np.asarray(state["material"]) >= BOX[1])):
                    break
                pts.append(state)
            if len(pts) != 2 or pts[0]["profile_hessian"] <= 0 or pts[1]["profile_hessian"] >= 0:
                break
            previous = [s["t"] for s in pts]
            trajectory.append(dict(delta=float(delta), minimum=pts[0], saddle=pts[1],
                                   predicted_t=real_roots_of_polynomial(polynomial(delta))))
        fold = None
        if trajectory:
            last = trajectory[-1]
            t0 = np.mean(previous)
            d0 = pred if pred is not None else last["delta"]*1.4
            df_scale = max(abs(d0), .001)

            def residual(z):
                t, delta = z[0], df_scale*z[1]
                fit = Fit(model, [0., delta], [old, read(delta)], gain_calibration=calibration)
                _, g, h, _, _ = exact_profile(fit, basis, t)
                return np.array([g, h])/gradient_scale

            if np.sign(d0) == sign:
                sol = root(residual, [t0, d0/df_scale], options={"eps": 1e-9}, tol=1e-7)
                t, d = sol.x[0], df_scale*sol.x[1]
                if np.sign(d) == sign and abs(d) < .09 and np.linalg.norm(sol.fun) < 1e-4:
                    fit = Fit(model, [0., d], [old, read(d)], gain_calibration=calibration)
                    fold = point(fit, basis, t, truth)
                    fold["delta"] = float(d)
                    fold["profile_third"] = float((residual(sol.x+[1e-4,0])[1] -
                                                   residual(sol.x-[1e-4,0])[1])*gradient_scale/(2e-4))
                    fold["profile_gradient_delta"] = float((residual(sol.x+[0,1e-4])[0] -
                                                            residual(sol.x-[0,1e-4])[0])*gradient_scale/(2e-4*df_scale))
                    fold["last_tracked_delta"] = last["delta"]
                    fold["tracked_distance_in_t"] = float(abs(t-t0))
        descriptions[str(sign)] = dict(trajectory=trajectory, fold=fold)
    return descriptions


def fold_case(amplitude, sign):
    payload = json.loads(DATA.read_text())
    model = FrequencyTM(payload["geometry"])
    critical = np.array(payload["numerical_critical_point"])
    direction = np.array(payload["right_unit_null"])
    # New amplitudes and side are specified independently of the R11 outcome.
    truth = critical + sign*amplitude*direction
    initial = critical - sign*amplitude*direction
    old = model.response(truth)[0]
    fit = Fit(model, [0.], [old])
    current, _ = fit.solve(initial)
    y, jac, h = model.response(current, hessian=True)
    u, s, vh = np.linalg.svd(ri(jac))
    v = vh[-1]
    bend = ri(np.einsum("kij,i,j->k", h, v, v))
    jump = -2*s[-1]/(u[:,-1] @ bend)
    proposal = current + jump*v - .5*jump*jump*(vh[:-1].T/s[:-1]) @ (u[:,:-1].T @ bend)
    return model, truth, np.array([current, fit.solve(proposal)[0]])


def visible_gap(model, pair, delta, calibration=None):
    f0, f1 = [model.response(m, 0.)[0] for m in pair]
    g0, g1 = [model.response(m, delta)[0] for m in pair]
    raw = np.linalg.norm(g1-g0)
    if calibration is None:
        return float(raw), float(np.linalg.norm(np.r_[f1-f0, g1-g0]))
    lam, _ = calibration
    # Symmetrize the profiled distance between the two nominal hypotheses.
    p0, p1 = np.r_[f0,g0,np.sqrt(lam)], np.r_[f1,g1,np.sqrt(lam)]
    gap_01 = np.linalg.norm(p1 - p0*np.vdot(p0,p1)/np.vdot(p0,p0))
    gap_10 = np.linalg.norm(p0 - p1*np.vdot(p1,p0)/np.vdot(p1,p1))
    return float(raw), float(min(gap_01,gap_10))


def new_read_generator(model, truth, noise, gain, short_noise):
    cache = {}
    def read(d, old=False, short=None):
        key = (d, old, short)
        if key not in cache:
            additive = noise[0 if old else 1] if short is None else 2 * short_noise[short]
            cache[key] = gain*model.response(truth,d)[0]+additive
        return cache[key]
    return read


def score_fit(fit, m, truth, cost):
    return dict(material_RMS=material_error(m, truth), factorizations=cost,
                material=m.tolist(), **fit.describe(m))


def baselines(model, pair, truth, old, read, calibration, delta, baseline_seed):
    fit = Fit(model, [0., delta], [old, read(delta)], gain_calibration=calibration)
    cap = 100
    start = pair[0]
    direct, n0 = fit.solve(start, cap)
    second, n1 = fit.solve(pair[1], cap)
    draws = qmc.Sobol(truth.size, scramble=True, seed=baseline_seed).random_base2(m=1)
    perturbation = 2*draws[0]-1
    perturbation *= np.linalg.norm(pair[1]-pair[0])/np.linalg.norm(perturbation)
    generic_start = np.clip(start+perturbation, BOX[0]+1e-8, BOX[1]-1e-8)
    other, n2 = fit.solve(generic_start, cap)
    pick = lambda x,y: min((x,y), key=lambda m: np.vdot(fit.evaluate(m)[0],fit.evaluate(m)[0]).real)

    # Useful no-additional-acquisition continuation on the same final samples.
    continuation = start
    cost = 0
    for w in (.25, .5, .75, 1.):
        weighted = Fit(model, [0.,delta], [old,read(delta)], weights=(1.,w),
                       gain_calibration=calibration)
        continuation, work = weighted.solve(continuation, cap//2)
        cost += work

    # A fixed old exposure and four new quarter-duration exposures.
    # Weight the quarter-duration responses by sqrt(time)=1/2.
    # Thus the new endpoint has reduced SNR and observations differ from direct TRF.
    times = np.linspace(delta/4, delta, 4)
    values = [read(d, short=i) for i, d in enumerate(times)]
    cumulative = start
    cumulative_cost = 0
    for j in range(1,5):
        cumulative_fit = Fit(model, [0., *times[:j]], [old,*values[:j]],
                             weights=(1., *(.5*np.ones(j))), gain_calibration=calibration)
        cumulative, used = cumulative_fit.solve(cumulative, cap//2)
        cumulative_cost += used
    return dict(delta=float(delta), direct=score_fit(fit,direct,truth,n0),
                pair_multistart=score_fit(fit,pick(direct,second),truth,n0+n1),
                generic_multistart=score_fit(fit,pick(direct,other),truth,n0+n2),
                same_read_continuation=score_fit(fit,continuation,truth,cost),
                cumulative=score_fit(cumulative_fit,cumulative,truth,cumulative_cost),
                cumulative_offsets=times.tolist(), pair_member=score_fit(fit,second,truth,n1),
                generic_member=score_fit(fit,other,truth,n2))


def gain_calibration_rule(noisy):
    return (ROBUST_SIGMA**2/GAIN_SD**2, 1.) if noisy == "profiled" else None


def run_one(model, truth, designs, tag, noisy=False, seed=RNG_SEED, natural=False):
    rng = np.random.default_rng(seed)
    shape = (6, len(model.response(truth)[0]))
    noise = np.zeros(shape, dtype=complex)
    gain = 1.
    calibration = gain_calibration_rule(noisy)
    if noisy:
        noise = (rng.normal(size=noise.shape)+1j*rng.normal(size=noise.shape))*ROBUST_SIGMA/np.sqrt(2)
        if calibration:
            gain = 1 + GAIN_SD*.65*np.exp(.63j)
    read = new_read_generator(model, truth, noise[:2], gain, noise[2:])
    old = read(0., old=True)
    pair = corrected_pair(model, old, *designs)
    prefit = np.asarray([np.linalg.norm(model.response(m)[0] - old) for m in pair])
    if np.linalg.norm(pair[1]-pair[0]) < 1e-3 or np.max(prefit) > 1e-8:
        return dict(tag=tag, skipped="old-data hypotheses failed validation",
                    old_residuals=prefit.tolist(), candidates=pair.tolist())
    f0, _, yp0, jp0 = model.response(pair[0],derivative=True)
    f1, _, yp1, _ = model.response(pair[1],derivative=True)
    gap0, gamma = f1-f0, yp1-yp0
    basis = candidate_basis(model, pair)
    polynomial, predicted = predicted_profile(model, old, pair[1], basis, calibration)
    reversed_basis = basis[0], -basis[1], basis[2], basis[3], basis[4]
    _, reverse = predicted_profile(model, old, pair[0], reversed_basis, calibration)
    exact = compare_folds(model, basis, pair, old, read, truth, calibration, polynomial, predicted) if basis[4] else {}
    grid = np.r_[-MENU[::-1], MENU]
    selection = []
    for d in grid:
        raw, gap = visible_gap(model, pair, d, calibration)
        actual = model.response(pair[1],d)[0]-model.response(pair[0],d)[0]
        selection.append(dict(delta=float(d), new_pair_gap=raw, joint_profile_gap=gap,
                              first_jet_gap=float(np.linalg.norm(gap0+d*gamma)),
                              first_jet_relative_error=float(np.linalg.norm(actual-gap0-d*gamma)/raw)))
    policy = max((s for s in selection if abs(s["delta"])<=.003),
                 key=lambda s:s["joint_profile_gap"])
    magnitudes = [abs(policy["delta"])]
    if all(d is not None for d in predicted.values()) and not natural:
        magnitudes.append(float(np.mean(np.abs(list(predicted.values())))))
    tests = [baselines(model, pair, truth, old, read, calibration, magnitude*sign, seed)
             for magnitude in magnitudes for sign in (-1,1)]
    h = 1e-5
    plus, jp = model.response(pair[0],h)
    minus, jm = model.response(pair[0],-h)
    jet_check = float(np.linalg.norm(yp0-(plus-minus)/(2*h))/np.linalg.norm(yp0))
    mixed_check = float(np.linalg.norm(jp0-(jp-jm)/(2*h))/np.linalg.norm(jp0))
    hf = Fit(model, [0., .002], [old,read(.002)], gain_calibration=calibration)
    vector = np.arange(truth.size)**.5/5
    hgrad = hf.evaluate(pair[0], True)[3] @ vector
    hp = (hf.evaluate(pair[0]+1e-3*vector)[2]-hf.evaluate(pair[0]-1e-3*vector)[2])/2e-3
    hess_check = float(np.linalg.norm(hgrad-hp)/np.linalg.norm(hp))
    results = dict(tag=tag, k0=model.k0, noisy=noisy, material=truth.tolist(), candidates=pair.tolist(),
                   seed=seed, natural=natural, old_pair_gap=float(np.linalg.norm(gap0)),
                   pair_frequency_jet_norm=float(np.linalg.norm(gamma)),
                   response_jet_relative_FD=jet_check, mixed_jet_relative_FD=mixed_check,
                   hessian_directional_relative_FD=hess_check, selection=selection,
                   selected_offset=policy["delta"], comparisons=tests,
                   observed_old=ri(old).tolist(), applied_complex_additive_error=ri(noise[:2].reshape(-1)).tolist(),
                   candidate_material_RMS=[material_error(m,truth) for m in pair],
                   a=basis[3], estimated_fold_center=basis[0].tolist(),
                   old_weak_direction=basis[1].tolist(), estimated_fold=basis[4],
                   predicted=predicted, predicted_reversed=reverse, exact=exact,
                   factorizations=model.factorizations,
                   noise_rms_per_complex_channel=ROBUST_SIGMA if noisy else 0.,
                   gain=[float(np.real(gain)), float(np.imag(gain))])
    return results


def main():
    out = ROOT/"results"/"freq_signed.json"
    records = []
    for amp, sign, noise, repeat in [
        (.23,1,False,0), (.47,-1,False,0), (.83,1,False,0),
        (.83,1,"known_gain",0), (.83,1,"profiled",0), (.83,1,"profiled",1), (.83,1,"profiled",2),
    ]:
        tag = f"r11_a{amp}_side{sign:+d}" + (f"_{noise}_{repeat}" if noise else "")
        result = run_one(*fold_case(amp,sign), tag,noise,seed=RNG_SEED+repeat)
        records.append(result)
        record(result)
        out.write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')

    natural_results = []
    cases = json.loads((R10/"branch_holdout.json").read_text())["cases"]
    with np.load(R10/"natural_inputs.npz") as archive:
        for tag in archive["case_ids"]:
            geometry = cases[tag]["geometry"]
            # Archive order is fixed by discovery, and includes no truth-based ranking.
            old_hypotheses = archive[f"{tag}/clean_witnesses"]
            if len(old_hypotheses) < 2:
                result = dict(tag=tag, archive_candidate_count=len(old_hypotheses),
                              skipped="fewer than two old-data roots in the historical archive")
            else:
                model = FrequencyTM(geometry)
                truth = archive[f"{tag}/truth"]
                result = run_one(model, truth, old_hypotheses[:2], tag, natural=True)
                result["archive_candidate_count"] = len(old_hypotheses)
                result["cohort_provenance"] = "R10 natural_inputs.npz clean witnesses, first two in discovery order"
                result["operator_abs_difference"] = max(
                    float(np.linalg.norm(model.kernels(0.)[0]-archive[f"{tag}/H"])),
                    float(np.linalg.norm(model.kernels(0.)[1]-archive[f"{tag}/B0"])),
                    float(np.linalg.norm(model.kernels(0.)[2][0]-archive[f"{tag}/R"])))
            natural_results.append(result)
            record(result)
            (ROOT/"results"/"freq_natural.json").write_text(
                json.dumps(natural_results, indent=2, allow_nan=False)+'\n')


def record(result):
    if "skipped" in result:
        print(result["tag"], result["skipped"], flush=True)
        return
    actual = {s:(len(p["trajectory"]), p["fold"]["delta"] if p["fold"] else None)
              for s,p in result["exact"].items()}
    print(result["tag"], "old", result["candidate_material_RMS"], "pred",
          result["predicted"], "actual", actual, flush=True)


def swapped_case():
    model, _, pair = fold_case(.83,1)
    # The two exact old roots generate the same old measurement. Exchange only
    # which one is used to produce new reads, and force the formerly true start.
    result = run_one(model, pair[0], pair[::-1], "r11_a0.83_same_old_swapped")
    record(result)
    (ROOT/"results"/"freq_same_old_swapped.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--same-old-swap", action="store_true")
    args = parser.parse_args()
    if args.same_old_swap:
        swapped_case()
    else:
        main()
