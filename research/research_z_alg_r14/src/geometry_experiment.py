"""Full scalar-TM data, matched optimization, and directional witnesses for R14."""
import json
import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import sys
from time import perf_counter, process_time

for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares
from scipy.stats import chi2, qmc

ROOT = Path(__file__).resolve().parents[1]
R10 = ROOT.parent / "research_z_alg_r10"
R12 = ROOT.parent / "research_z_alg_r12"
sys.path.insert(0, str(R12 / "src"))
from freq_model import FrequencyTM


BOX = (.5, 30.)
SIGMA = 3e-7


def realify(z):
    return np.stack((z.real, z.imag), axis=1).reshape((2 * len(z), *z.shape[1:]))


def material_rms(x, truth):
    return float(np.linalg.norm(x - truth) / np.sqrt(truth.size))


def complex_pairs(x):
    return np.column_stack((x.real, x.imag)).tolist()


@dataclass
class Block:
    green: np.ndarray
    incident: np.ndarray
    receiver: np.ndarray
    host: np.ndarray
    active: np.ndarray


@dataclass
class Read:
    name: str
    blocks: list
    data: np.ndarray
    old_count: int
    seeds: np.ndarray
    epsilon: float
    gamma: float
    truth: np.ndarray
    physical_gain: complex
    known_data: np.ndarray | None = None


class FullTM:
    """A common LU gives the full state, exact material Jacobian and optional Hessian."""

    def __init__(self, blocks):
        self.blocks = blocks
        self.counts = Counter()

    def evaluate(self, material, second=False):
        ys, js, hs = [], [], []
        for b in self.blocks:
            m = b.host.copy()
            m[b.active] = material
            n = len(m)
            step = lu_factor(np.eye(n) - m[:, None] * b.green)
            source_count = b.incident.shape[1]
            curr = lu_solve(step, m[:, None] * b.incident)
            e = (b.incident + b.green @ curr)[b.active]
            a = lu_solve(step, b.receiver.T, trans=1).T[:, b.active]
            ys.append((b.receiver @ curr).T.reshape(-1))
            js.append(np.einsum("ri,is->sri", a, e).reshape(-1, len(material)))
            self.counts.update(factorizations=1, rhs=source_count + b.receiver.shape[0])
            if second:
                gq = (b.green[b.active] @ lu_solve(step, np.eye(n)[:, b.active]))
                gq = gq[:, b.active]
                mixed = np.einsum("ri,ij,js->srij", a, gq, e).reshape(-1, len(material), len(material))
                hs.append(mixed + mixed.swapaxes(1, 2))
                self.counts.update(rhs=len(material))
        return np.concatenate(ys), np.concatenate(js), np.concatenate(hs) if second else None


class Fit:
    """Variable projection of one shared complex gain under the stated ellipsoid."""

    def __init__(self, wave, observed, epsilon, gamma=0.):
        self.wave = wave
        self.observed = observed
        self.epsilon = epsilon
        self.penalty = epsilon / gamma if gamma else None
        self.scale = max(np.linalg.norm(observed), np.finfo(float).tiny)
        self.last = None

    def state(self, m, second=False):
        if self.last is not None and np.array_equal(self.last[0], m) and not second:
            return self.last[1]
        y, j, h = self.wave.evaluate(m, second)
        if self.penalty is None:
            g, r, jr = 1. + 0.j, y - self.observed, j
        else:
            lam = self.penalty**2
            denom = np.vdot(y, y).real + lam
            g = (np.vdot(y, self.observed) + lam) / denom
            gd = (j.conj().T @ self.observed - 2*g*np.real(y.conj() @ j)) / denom
            r = np.r_[g*y - self.observed, self.penalty*(g - 1)]
            jr = np.vstack((g*j + y[:, None]*gd, self.penalty*gd))
        hh = None
        if second:
            fixed_j = g*j
            hh = np.real(fixed_j.conj().T @ fixed_j +
                         np.einsum("k,kij->ij", r[:len(y)].conj(), g*h))
            if self.penalty is not None:
                cross = np.real(fixed_j.conj().T @ np.column_stack((y, 1j*y)) +
                                np.column_stack((j.conj().T @ r[:len(y)],
                                                 -1j*j.conj().T @ r[:len(y)])))
                hh -= cross @ cross.T / denom
        result = r, jr, g, hh
        self.last = m.copy(), result
        return result

    def score(self, m):
        r, _, g, _ = self.state(m)
        return float(np.linalg.norm(r)/self.epsilon), g

    def solve(self, m, bounds=BOX, cap=140):
        cost, t0, cpu0 = self.wave.counts.copy(), perf_counter(), process_time()
        lo, hi = bounds
        interior = min(1e-8, np.min(np.asarray(hi)-np.asarray(lo))/1e3)
        fit = least_squares(lambda x: realify(self.state(x)[0])/self.scale,
                            np.clip(m, lo + interior, hi - interior),
                            jac=lambda x: realify(self.state(x)[1])/self.scale,
                            bounds=bounds, method="trf",
                            ftol=1e-13, xtol=1e-13, gtol=2e-15, max_nfev=cap)
        return fit.x, dict(self.wave.counts-cost, nfev=fit.nfev,
                           status=fit.status, wall_seconds=perf_counter()-t0,
                           cpu_seconds=process_time()-cpu0)

    def diagnostics(self, m, truth, bounds=BOX, second=True):
        r, j, g, h = self.state(m, second=second)
        grad = np.real(j.conj().T @ r)/self.scale**2
        lo, hi = bounds
        free = (m > lo+1e-7) & (m < hi-1e-7)
        hess_min = float(np.linalg.eigvalsh(h[np.ix_(free,free)]/self.scale**2)[0]) if free.any() and second else None
        projected = m-np.clip(m-grad, lo, hi)
        return dict(material=m.tolist(), material_rms=material_rms(m,truth),
                    additive_over_epsilon=float(np.linalg.norm(r[:len(self.observed)])/self.epsilon),
                    score=float(np.linalg.norm(r)/self.epsilon), gain=[float(g.real),float(g.imag)],
                    projected_grad_inf=float(np.linalg.norm(projected,np.inf)),
                    grad_inf=float(np.linalg.norm(grad,np.inf)),
                    free_hessian_min=hess_min, active_boxes=int((~free).sum()))


def r10_data(group, arm="active"):
    with np.load(R10 / "results" / f"{group}_inputs.npz", allow_pickle=False) as d:
        for name in d["case_ids"]:
            n = len(d[f"{name}/H"])
            block = Block(d[f"{name}/H"], d[f"{name}/{arm}/B"],
                          d[f"{name}/R"].reshape(1,-1), d[f"{name}/host"], np.arange(n))
            gain = 1 + complex(d[f"{name}/gain_error"])
            y = d[f"{name}/{arm}/y"]
            yield Read(str(name), [block], y, len(d[f"{name}/y_old"]),
                       d[f"{name}/central_roots"], float(d["eps"]), float(d["gamma"]),
                       d[f"{name}/truth"], gain, y/gain)


def frequency_block(physical, offset):
    g, b, s, *_ = physical.kernels(offset)
    n = len(g)
    return Block(g, b, s, np.zeros(n), np.arange(n))


def r12_frequency_data():
    """Replay fixed archived old noise and the paired new noise vector; no redraw per fit."""
    cfg = json.loads((R12 / "src" / "freq_r11_input.json").read_text())
    mates = json.loads((R12 / "results" / "freq_signed.json").read_text())
    wave = FrequencyTM(cfg["geometry"])
    eps = SIGMA*np.sqrt(chi2.ppf(.99,12)/2)
    for idx in (3, 4):
        row = mates[idx]
        target = np.array(row["material"])
        candidates = np.array(row["candidates"])
        gain = complex(*row["gain"])
        old = np.array(row["observed_old"])
        old = old[::2] + 1j*old[1::2]
        noise = np.array(row["applied_complex_additive_error"])
        noise = (noise[::2] + 1j*noise[1::2]).reshape(2,-1)
        # Same two realized additive vectors as in R12. The two signed offsets
        # are a fixed development pair and are analyzed separately.
        for offset in (-.001545,.001545):
            blocks = [frequency_block(wave,0.),frequency_block(wave,offset)]
            fnew = wave.response(target,offset)[0]
            y = np.r_[old, gain*fnew+noise[1]]
            gamma = 4e-4 if idx==4 else 0.
            # Outer ellipsoid of the independently bounded noise and gain disks.
            epsilon = np.sqrt(2)*eps if gamma else eps
            if gamma:
                gamma *= np.sqrt(2)
            yield Read(f"fold-{idx}-{offset:+g}", blocks, y, len(old), candidates, epsilon, gamma,
                       target, gain, y/gain)


def stationary_witness(wave, start, reference, g_m=1., g_n=1., epsilon=1., gamma=0.):
    """A specific chord may disprove positive direction even when outputs differ."""
    a, _, _ = wave.evaluate(reference)
    records = []
    for t in (0., .001, .01, .04, .1, .2, .4):
        m = start + t*(reference-start)
        g = g_m + t*(g_n-g_m)
        h = m-reference
        dg = g-g_n
        response, jac, _ = wave.evaluate(m)
        derivative = (g*(jac@h) + dg*response)/epsilon
        chord = (g*response-g_n*a)/epsilon
        direction = np.vdot(derivative,chord).real
        if gamma:
            direction += abs(dg/gamma)**2
        h_norm2 = np.dot(h,h) + (abs(dg/gamma)**2 if gamma else 0.)
        records.append(dict(t=t, direction_per_norm2=float(direction/h_norm2),
                            chord_over_epsilon=float(np.linalg.norm(chord)),
                            material=m.tolist(), gain=[float(g.real),float(g.imag)]))
    return records


def run(read, seed_id, use_known_gain=False, certify=True):
    wave = FullTM(read.blocks)
    known = use_known_gain or not read.gamma
    gamma = 0. if known else read.gamma
    data = read.known_data if use_known_gain else read.data
    eps = read.epsilon/abs(read.physical_gain) if use_known_gain else read.epsilon
    fit = Fit(wave, data, eps, gamma)
    n = len(read.blocks[0].active)
    draws = qmc.Sobol(n, scramble=True, seed=49019+seed_id).random_base2(2)
    cold = np.full(n, np.mean(BOX))
    starts = [cold, *(.5+29.5*draws[:3])]
    costs = {}
    begin, t0 = wave.counts.copy(), perf_counter()
    central_scores = []
    old_scores = []
    new_scores = []
    for x in read.seeds:
        y, _, _ = wave.evaluate(x)
        yo, fo = data[:read.old_count], y[:read.old_count]
        if gamma:
            lam = (eps/gamma)**2
            old_gain = (np.vdot(fo,yo)+lam)/(np.vdot(fo,fo).real+lam)
            old_residual = np.r_[old_gain*fo-yo, np.sqrt(lam)*(old_gain-1)]
        else:
            old_gain, old_residual = 1., fo-yo
        old_scores.append(float(np.linalg.norm(old_residual)/eps))
        new_scores.append(float(np.linalg.norm(old_gain*y[read.old_count:]-data[read.old_count:])/eps))
        central_scores.append(fit.score(x)[0])
    # All endpoint choices are data-only. Truth enters diagnostics below.
    if central_scores:
        pick = int(np.argmin(central_scores))
        starts.insert(0, read.seeds[pick])
    else:
        pick = None
    costs["candidate_scoring"] = dict(wave.counts-begin, wall_seconds=perf_counter()-t0)

    solutions, costs["TRF"] = [], []
    for start in starts:
        found, cost = fit.solve(start)
        costs["TRF"].append(cost)
        solutions.append(found)
    candidates = [fit.score(x)[0] for x in solutions]
    winner_index = int(np.argmin(candidates))
    selected = solutions[winner_index]
    selected_score, selected_gain = fit.score(selected)

    restart = []
    for idx, x in enumerate(read.seeds):
        if idx != pick:
            found, cost = fit.solve(x)
            costs["TRF"].append(cost)
            solutions.append(found)
            restart.append(idx)
    search = dict(direct_index=0, generic_count=len(starts), selected_index=winner_index,
                  selected_score=selected_score, alternative_seed_indices=restart)

    certificates = {}
    restrictions = {}
    secants = []
    interior_minima = []
    if certify:
        from geometry_bounds import directions
        tic = perf_counter()
        certificates = directions(read.blocks, selected, selected_gain, eps, gamma, data)
        costs["boxes"] = dict(wall_seconds=perf_counter()-tic)
        width = certificates["relative_width"]
        if width > 0:
            limits = (np.maximum(BOX[0], selected-width),
                      np.minimum(BOX[1], selected+width))
            origin = starts[1] if pick is not None else starts[0]
            regional, cost = fit.solve(origin, limits)
            restrictions = dict(
                box_half_width=width, TRF_cost=cost,
                point=fit.diagnostics(regional,read.truth,bounds=limits),
                truth_inside=bool(np.all(read.truth>=limits[0]) and np.all(read.truth<=limits[1])),
                truth_gain_inside=bool(not gamma or
                    max(abs((read.physical_gain-selected_gain).real),abs((read.physical_gain-selected_gain).imag))
                    <= gamma*certificates["gain_half_width"]),
            )

    # No optimizer, certificate center, or starting box uses these truth diagnostics.
    initial = []
    for m, old, new, score in zip(read.seeds,old_scores,new_scores,central_scores,strict=True):
        initial.append(dict(material=m.tolist(), score=score, old_score=old,
                            new_score_at_old_gain=new, material_rms=material_rms(m,read.truth)))
    begin, t0 = wave.counts.copy(), perf_counter()
    diagnostics = [fit.diagnostics(m, read.truth) for m in solutions]
    true_score = fit.score(read.truth)[0]
    costs["diagnostics"] = dict(wave.counts-begin, wall_seconds=perf_counter()-t0)
    if read.name.startswith("fold") or read.name in {"fresh6-00","fresh6-09","confirm6-00"}:
        from geometry_bounds import rigorous_witness
        begin, t0 = wave.counts.copy(), perf_counter()
        for idx, point in enumerate(diagnostics):
            if point["material_rms"] <= .1 or idx==winner_index:
                continue
            truth_gain = read.physical_gain if gamma else 1.
            witness = dict(fit_index=idx, to_selected=stationary_witness(
                wave, solutions[idx], selected, complex(*point["gain"]), selected_gain, eps, gamma),
                to_truth=stationary_witness(wave, solutions[idx], read.truth,
                    complex(*point["gain"]), truth_gain, eps, gamma))
            for target, reference, gain0 in (("to_selected",selected,selected_gain),
                                              ("to_truth",read.truth,truth_gain)):
                strongest = min(witness[target],key=lambda z:z["direction_per_norm2"])
                if strongest["direction_per_norm2"] < 0:
                    witness[target+"_certified_negative"] = dict(t=strongest["t"],**rigorous_witness(
                        read.blocks,np.array(strongest["material"]),reference,complex(*strongest["gain"]),
                        gain0,eps,gamma))
            secants.append(witness)
        costs["chords"] = dict(wave.counts-begin,wall_seconds=perf_counter()-t0)
    if read.name == "fold-3-+0.001545" and not gamma and certify:
        from geometry_bounds import local_minimum
        begin, t0 = wave.counts.copy(), perf_counter()
        for idx, point in enumerate(diagnostics):
            if not (.1 < point["material_rms"] < 1.) or point["active_boxes"]:
                continue
            m = solutions[idx].copy()
            for _ in range(3):
                r, j, _, h = fit.state(m,second=True)
                m -= np.linalg.solve(h,np.real(j.conj().T@r))
            interior_minima.append(dict(fit_index=idx,proof=local_minimum(
                read.blocks,m,data,eps,read.truth,BOX)))
        costs["wrong_minimum_verification"] = dict(
            wave.counts-begin,wall_seconds=perf_counter()-t0)

    return dict(case=read.name, known_gain=known, epsilon=eps, gamma=gamma, truth=read.truth.tolist(),
                physical_gain=[read.physical_gain.real,read.physical_gain.imag], true_score=true_score,
                observed=complex_pairs(data), start_candidates=initial, selected_candidate=pick,
                multistart=search, fits=diagnostics, certificates=certificates, regional=restrictions,
                chord_witnesses=secants, certified_wrong_minima=interior_minima, costs=costs)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("group", choices=("natural", "confirm", "fold"))
    p.add_argument("--no-certificates", action="store_true")
    p.add_argument("--known-gain", action="store_true")
    p.add_argument("--arm", choices=("active", "random", "fisher"), default="active")
    p.add_argument("--cases", nargs="+")
    args = p.parse_args()
    dataset = r12_frequency_data() if args.group=="fold" else r10_data(args.group,args.arm)
    result = []
    for i, read in enumerate(dataset):
        if args.cases and read.name not in args.cases:
            continue
        row = run(read,i,args.known_gain,not args.no_certificates)
        result.append(row)
        direct = row["fits"][0]
        selected = row["fits"][row["multistart"]["selected_index"]]
        reg = row["regional"]
        print(read.name, "candidates=", len(read.seeds), "truth_score=",round(row["true_score"],4),
              "direct=",round(direct["material_rms"],5),round(direct["score"],4),
              "selected=",round(selected["material_rms"],5),round(selected["score"],4),
              "certified=",row["certificates"].get("relative_width"),
              "truth_in_box=",reg.get("truth_inside"), flush=True)
    suffix = ("-known" if args.known_gain else "") + ("" if args.arm=="active" else f"-{args.arm}")
    target = ROOT / "results" / f"geometry-{args.group}{suffix}.json"
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
