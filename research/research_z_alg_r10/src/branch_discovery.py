"""Conditional branch discovery with one common numerically accepted Maxwell solution."""

import argparse
import json
import os
from pathlib import Path
from statistics import median
import sys
from time import perf_counter, process_time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
R2 = ROOT.parent / "research_z_alg_r2"
R5 = ROOT.parent / "research_z_alg_r5_measurement"
sys.path.insert(0, str(R2 / "src"))
from physics import Model

LO, HI = 0.5, 30.0
SPAN = HI - LO
MATERIAL_TOL = 1e-3  # RMS contrast units for branch deduplication
RESPONSE_TOL = 1e-10
METHODS = ("TRF", "guided TRF", "shifted deflation", "ordinary secant",
           "joint reciprocal", "projected reciprocal")
WEIGHTED = "radius-weighted reciprocal"
ALG_METHODS = (*METHODS[3:], WEIGHTED)


def pack(y):
    return np.stack((y.real, y.imag), axis=-1).reshape(-1)


def rms(m):
    return np.linalg.norm(m) / np.sqrt(m.size)


class Experiment:
    def __init__(self, geometry, observation=None):
        self.model = Model(**geometry)
        self.h, self.b = self.model.G, self.model.E0
        self.r = self.model.S[0]
        self.n = self.model.n
        self.eye = np.eye(self.n, dtype=complex)
        if observation is not None:
            self.obs = observation
            self.scale = np.linalg.norm(observation)

    def response(self, m):
        lu = lu_factor(self.eye - m[:, None] * self.h, check_finite=False)
        current = lu_solve(lu, m[:, None] * self.b, check_finite=False)
        field = self.b + self.h @ current
        transfer = lu_solve(lu, self.r, trans=1, check_finite=False)
        derivative = (transfer[:, None] * field).T
        jac = np.stack((derivative.real, derivative.imag), axis=1).reshape(-1, self.n)
        return self.r @ current, jac

    def dual(self, m):
        return np.linalg.solve((self.eye - m[:, None] * self.h).T, self.r) * m

    def validate(self, m):
        if not np.all(np.isfinite(m)) or np.any(m < LO - 1e-10) or np.any(m > HI + 1e-10):
            return None
        state = self.model.forward(m, independent=True)
        err = np.linalg.norm(state["y"][0] - self.obs) / self.scale
        valid = err <= RESPONSE_TOL and state["state_relative_residual"] <= 1e-10
        valid &= state["resolvent_sigma_min"] >= 1e-6
        if not valid:
            return None
        _, jac = self.response(m)
        return {
            "m": m.tolist(), "response_relative": float(err),
            "state_relative": state["state_relative_residual"],
            "resolvent_sigma_min": state["resolvent_sigma_min"],
            "old_material_jac_sigma_min": float(np.linalg.svd(jac, compute_uv=False)[-1]),
        }


def trf(p, m0, check, excluded=(), max_evals=120, whiten=None, tolerance=1e-14):
    previous = None
    residual, jacobian = None, None

    def evaluate(x):
        nonlocal previous, residual, jacobian
        check()
        if previous is not None and np.array_equal(previous, x):
            return
        previous = x.copy()
        y, jac = p.response(LO + SPAN * x)
        if whiten is None:
            residual, jacobian = pack(y - p.obs) / p.scale, jac * (SPAN / p.scale)
        else:
            derivative = jac[::2] + 1j * jac[1::2]
            transformed = whiten.T @ derivative
            residual = pack((y - p.obs) @ whiten)
            jacobian = SPAN * np.stack((transformed.real, transformed.imag), axis=1).reshape(jac.shape)
        if excluded:
            # Shifted norm deflation: D = product(1 + ||x-root||^(-2)).
            factor = 1.
            grad_log_d = np.zeros_like(x)
            for root in excluded:
                diff = x - root
                distance_sq = max(float(diff @ diff), 1e-20)
                factor *= 1. + 1. / distance_sq
                grad_log_d -= 2. * diff / (distance_sq * (1. + distance_sq))
            jacobian = factor * (jacobian + np.outer(residual, grad_log_d))
            residual = factor * residual

    def fun(x):
        evaluate(x)
        return residual

    def jac(x):
        evaluate(x)
        return jacobian

    sol = least_squares(
        fun, (np.clip(m0, LO, HI) - LO) / SPAN, jac=jac, bounds=(0., 1.),
        ftol=1e-13, xtol=1e-13, gtol=tolerance, max_nfev=max_evals,
    )
    m = LO + SPAN * sol.x
    err = np.linalg.norm(p.response(m)[0] - p.obs) / p.scale
    return m, float(err), sol.nfev


class TimeExpired(Exception):
    pass


def no_deadline():
    pass


def find_anchor(p, sobol, orthogonalize):
    start = process_time()
    min_error = None
    count = 0
    _, sigma, vh = np.linalg.svd(p.b, full_matrices=False)
    whitening = (vh.conj().T / np.maximum(sigma, sigma[0] * 1e-12)) * sigma[0] / p.scale
    settings = [(None, "raw")]
    if orthogonalize and sigma[-1] < sigma[0] / 100:
        settings.append((whitening, "orthogonalized old channels"))
    for transform, name in settings:
        for m0 in [np.full(p.n, (HI + LO) / 2), *sobol[:64]]:
            count += 1
            m, error, _ = trf(p, m0, no_deadline, whiten=transform, max_evals=240 if orthogonalize else 120)
            min_error = error if min_error is None else min(error, min_error)
            if error <= RESPONSE_TOL:
                accepted = p.validate(m)
                if accepted is not None:
                    return m, {"starts": count, "fit_metric": name,
                               "cpu_seconds": process_time() - start, **accepted}
    return None, {"starts": count, "cpu_seconds": process_time() - start, "best_residual": min_error}


def guided_seeds(p, m, sobol):
    _, jac = p.response(m)
    _, _, vh = np.linalg.svd(jac)
    seeds = []
    for i, draw in enumerate(sobol):
        if i < 8:
            d = (-1 if i % 2 else 1) * vh[-1 - (i // 4)]
            distance = np.min(np.where(d > 0, (HI - m) / np.maximum(d, 1e-30),
                                       (m - LO) / np.maximum(-d, 1e-30)))
            fraction = 0.38 if i % 4 < 2 else 0.82
            seeds.append(np.clip(m + fraction * distance * d, LO, HI))
        seeds.append(draw)
    return seeds


class Secant:
    def __init__(self, p, anchor, check):
        self.p, self.check = p, check
        self.w = p.dual(anchor)
        self.v = p.r + self.w @ p.h
        self.f0 = np.imag(self.w * self.v.conj())
        left_vectors, singular, _ = np.linalg.svd(p.b, full_matrices=True)
        rank = np.sum(singular > 1e-12 * singular[0])
        self.q = left_vectors[:, rank:].conj().T
        self.qh = self.q @ p.h
        self.d = self.q.shape[0]
        self.scale = np.linalg.norm(self.w * self.v.conj())
        self.step_scale = max(np.linalg.norm(self.w), 1e-12)
        self.reciprocal_scale = max(np.linalg.norm(self.q * self.qh.conj()), 1e-12)

    def start(self, m):
        shift = self.p.dual(m) - self.w
        return shift @ self.q.conj().T

    def quadratic_seeds(self):
        a, b = self.q * self.v.conj(), self.w * self.qh.conj()
        jac = np.concatenate(((a + b).imag, (a - b).real)).T * self.step_scale / self.scale
        u, sigma, vh = np.linalg.svd(jac, full_matrices=False)
        seeds = []
        for component in [len(sigma) - 1, len(sigma) - 2]:
            weak = vh[component]
            direction = self.step_scale * (weak[:self.d] + 1j * weak[self.d:])
            quadratic = ((direction @ self.q) * (direction @ self.qh).conj()).imag / self.scale
            curvature = u[:, component] @ quadratic
            if abs(curvature) < 1e-12:
                continue
            step = -sigma[component] / curvature
            if not 1e-6 < abs(step) < 2.:
                continue
            strong = sigma > max(10. * sigma[component], 1e-10)
            correction = -((u[:, strong].T @ quadratic) / sigma[strong]) @ vh[strong]
            corrected = self.step_scale * (correction[:self.d] + 1j * correction[self.d:])
            for candidate in (step * direction + step**2 * corrected, step * direction):
                w = self.w + candidate @ self.q
                v = self.p.r + w @ self.p.h
                if np.min(np.abs(v)) > 1e-15 and np.linalg.norm(candidate) <= 2. * self.step_scale:
                    seeds.append((np.clip((w / v).real, LO, HI), candidate))
        return seeds

    def recover(self, delta):
        wl = self.w + delta
        vl = self.p.r + wl @ self.p.h
        if np.any(np.abs(vl) < 1e-15):
            return None, None, None
        material = wl / vl
        m = material.real
        if not np.all(np.isfinite(m)) or np.any(m < LO - 1e-7) or np.any(m > HI + 1e-7):
            return None, None, None
        err = np.linalg.norm(self.p.response(m)[0] - self.p.obs) / self.p.scale
        return m, float(err), float(np.max(np.abs(material.imag)))

    def ordinary(self, initial):
        def value(x, jac):
            self.check()
            t = self.step_scale * (x[:self.d] + 1j * x[self.d:])
            w, v = self.w + t @ self.q, self.v + t @ self.qh
            if not jac:
                return ((w * v.conj()).imag - self.f0) / self.scale
            a, b = self.q * v.conj(), w * self.qh.conj()
            return self.step_scale * np.concatenate(((a + b).imag, (a - b).real)).T / self.scale

        x0 = np.r_[initial.real, initial.imag] / self.step_scale
        sol = least_squares(lambda x: value(x, False), x0, jac=lambda x: value(x, True),
                            ftol=1e-12, xtol=1e-12, gtol=1e-13, max_nfev=120)
        t = self.step_scale * (sol.x[:self.d] + 1j * sol.x[self.d:])
        return t @ self.q, sol.nfev

    def reciprocal(self, initial, project, weighted=False):
        j = int(np.argmax(np.abs(initial)))
        if abs(initial[j]) < 1e-15:
            initial = np.eye(self.d, dtype=complex)[j]
        indices = np.arange(self.d) != j
        d = self.d - 1
        z0 = initial[indices] / initial[j]
        dn = np.concatenate((self.q[indices], 1j * self.q[indices]))
        dh = np.concatenate((self.qh[indices], 1j * self.qh[indices]))
        da = (dn * self.v.conj() + self.w * dh.conj()).imag.T
        db = (dn * self.v.conj() - self.w * dh.conj()).real.T
        a_derivative = np.stack((da, db), axis=1)
        previous = None
        cached = None

        def terms(x, use_projection):
            nonlocal previous, cached
            self.check()
            if previous is not None and np.array_equal(previous, x):
                return cached
            previous = x.copy()
            coef = np.ones(self.d, dtype=complex)
            coef[indices] = x[:d] + 1j * x[d:2*d]
            n, nh = coef @ self.q, coef @ self.qh
            a = np.column_stack(((n * self.v.conj() + self.w * nh.conj()).imag,
                                 (n * self.v.conj() - self.w * nh.conj()).real))
            c = (n * nh.conj()).imag
            dc = (dn * nh.conj() + n * dh.conj()).imag.T
            u, sigma, vh = np.linalg.svd(a, full_matrices=False)
            valid = sigma > 1e-12 * sigma[0]
            inv_s = np.divide(1, sigma, where=valid, out=np.zeros_like(sigma))
            s = -((vh.T * inv_s) @ (u.T @ c)) if use_projection else x[-2:] / self.step_scale
            raw = a @ s + c
            draws = np.einsum("ijk,j->ik", a_derivative, s) + dc
            if use_projection:
                # Differentiate A.T r = 0. The extra residual term is needed away from a root.
                pinv_a = (vh.T * inv_s) @ u.T
                normal_pinv = (vh.T * inv_s**2) @ vh
                ds = -pinv_a @ draws - normal_pinv @ np.einsum("ijk,i->jk", a_derivative, raw)
                draws += a @ ds
            norm_sq = float(np.vdot(n, n).real)
            factor = norm_sq * self.reciprocal_scale
            d_norm = 2 * (dn @ n.conj()).real
            jac_z = (draws - np.outer(raw, d_norm / norm_sq)) / factor
            jac = jac_z if use_projection else np.column_stack((jac_z, a / (factor * self.step_scale)))
            residual = raw / factor
            if weighted and use_projection:
                # Restore some amplitude sensitivity without a vanishing weight at t -> 0.
                radius_ratio_sq = (self.step_scale**2 / norm_sq) * (s @ s)
                weight = .05 + .95 / (1. + radius_ratio_sq)
                d_ratio = (self.step_scale**2 / norm_sq) * (2. * s @ ds - (s @ s) * d_norm / norm_sq)
                d_weight = -.95 * d_ratio / (1. + radius_ratio_sq)**2
                jac = weight * jac + np.outer(residual, d_weight)
                residual = weight * residual
            cached = (residual, jac, n, s)
            return cached

        x0 = np.r_[z0.real, z0.imag]
        if project:
            if not x0.size:
                _, _, n, s = terms(x0, True)
                ts = complex(*s)
                return (None if abs(ts) < 1e-14 else n / ts.conjugate()), 1
            bound = np.full(x0.size, 2.)
        else:
            s0 = terms(x0, True)[3]
            x0 = np.r_[x0, s0 * self.step_scale]
            bound = np.r_[np.full(2*d, 2.), np.inf, np.inf]
        # The joint arm's least-squares initial amplitude is counted in this attempt.
        previous = None
        sol = least_squares(lambda x: terms(x, project)[0], x0,
                            jac=lambda x: terms(x, project)[1], bounds=(-bound, bound),
                            ftol=1e-12, xtol=1e-12, gtol=1e-13, max_nfev=120)
        _, _, n, s = terms(sol.x, project)
        ts = complex(*s)
        return (None if abs(ts) < 1e-14 else n / ts.conjugate()), sol.nfev


def search(p, anchor, name, plain_seeds, budget, quadratic_seeds):
    wall = perf_counter()
    start = process_time()
    attempts, accepted = [], []
    best_residual = None

    def check():
        if process_time() - start >= budget:
            raise TimeExpired

    seeds = plain_seeds if name == "TRF" else guided_seeds(p, anchor, plain_seeds)
    starts = [(m, None) for m in seeds]
    setup = name in ALG_METHODS or (quadratic_seeds and name != "TRF")
    secant = Secant(p, anchor, check) if setup else None
    if quadratic_seeds and secant is not None and secant.d:
        starts = secant.quadratic_seeds() + starts
    for i, (m0, projected_start) in enumerate(starts):
        began = process_time()
        try:
            check()
            pre_polish = None
            imag_defect = None
            if name not in ALG_METHODS:
                excluded = [(anchor - LO) / SPAN, *((np.asarray(r["m"]) - LO) / SPAN for r in accepted)]
                m, error, evaluations = trf(
                    p, m0, check, excluded=excluded if name == "shifted deflation" else (),
                )
            elif secant.d:
                initial = secant.start(m0) if projected_start is None else projected_start
                if name == "ordinary secant":
                    delta, evaluations = secant.ordinary(initial)
                else:
                    delta, evaluations = secant.reciprocal(initial, project=name != "joint reciprocal",
                                                           weighted=name == WEIGHTED)
                m, error, imag_defect = secant.recover(delta) if delta is not None else (None, None, None)
                if (m is not None and rms(m - anchor) > MATERIAL_TOL
                        and RESPONSE_TOL < error < 1e-6 and imag_defect < 1e-5):
                    pre_polish = error
                    m, error, extra = trf(p, m, check, max_evals=12)
                    evaluations += extra
            else:
                break
            record = {"start": i + 1, "nfev": evaluations, "data_residual": error,
                      "quadratic_seed": projected_start is not None}
            if imag_defect is not None:
                record["imaginary_material_defect"] = imag_defect
            if pre_polish is not None:
                record["pre_polish_residual"] = pre_polish
            if error is not None:
                best_residual = error if best_residual is None else min(error, best_residual)
                if (error <= RESPONSE_TOL and rms(m - anchor) > MATERIAL_TOL
                        and all(rms(m - root["m"]) > MATERIAL_TOL for root in accepted)):
                    validated = p.validate(m)
                    if validated is not None:
                        accepted.append({
                            **validated, "start": i + 1, "cpu_seconds": process_time() - start,
                            "distance_to_anchor_RMS": float(rms(m - anchor)),
                        })
                        record["new_branch"] = len(accepted)
            record["cpu_seconds"] = process_time() - began
            attempts.append(record)
        except TimeExpired:
            attempts.append({"start": i + 1, "time_limited": True, "cpu_seconds": process_time() - began})
            break
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as error:
            attempts.append({"start": i + 1, "error": str(error), "cpu_seconds": process_time() - began})
    return {
        "cpu_seconds": process_time() - start, "wall_seconds": perf_counter() - wall,
        "best_search_residual": best_residual, "new_branches": accepted, "attempts": attempts,
    }


def cases(fresh_count, stress_count, cohort, seed_offset):
    if cohort in ("all", "archive"):
        meta = json.loads((R5 / "inputs/cases.json").read_text())
        for name, case in sorted(meta.items()):
            with np.load(R5 / "inputs" / (name + ".npz")) as data:
                observation = data["old_clean"].copy()
            yield "R5-" + name, case["geometry"], observation, None
    if cohort == "archive":
        return
    g = json.loads((R2 / "experiments/GENERAL6_CONTRACT.json").read_text())["geometry"]
    geometry_rng = np.random.default_rng(421006 + seed_offset)
    truth_rng = np.random.default_rng(421029 + seed_offset)
    for family in range(3):
        base = np.asarray(g["centers"])
        centers = base + geometry_rng.uniform(-.045, .045, base.shape)
        angles = np.array([-2.5, -.9, 2.42]) + .17 * (family + 1)
        sources = (1.04 * np.column_stack((np.cos(angles), np.sin(angles)))).tolist()
        geometry = dict(centers=centers.tolist(), radius=.045 + .002 * family,
                        sources=sources, k=6.5 + .8 * family,
                        receivers=[[.21 + .063 * family, -.96 + .025 * family]])
        for i in range(family, fresh_count, 3):
            truth = truth_rng.uniform(10., 30., size=6)
            yield f"fresh6-{i:02d}", geometry, None, truth
    angles = np.array([-2.95, -1.88, -.60, .39, 1.40, 2.58])
    sources = (1.16 * np.column_stack((np.cos(angles), np.sin(angles)))).tolist()
    xy = np.array([(x, y) for y in [-.31, 0., .33] for x in [-.47, -.14, .17, .49]])
    centers = xy + geometry_rng.uniform(-.035, .035, xy.shape)
    geometry = dict(centers=centers.tolist(), radius=.045, sources=sources, k=8.1, receivers=[[.31, .92]])
    for i in range(stress_count):
        truth = truth_rng.uniform(10., 30., size=12)
        yield f"fresh12-{i:02d}", geometry, None, truth


def archive_references(case):
    # Historical references only opened after completing all discovery searches.
    found = []
    for method in ("Q", "TRF64", "reciprocal64"):
        data = json.loads((R5 / "results" / case / "clean" / f"{method}.json").read_text())
        for row in data["accepted_roots"]:
            material = np.asarray(row)
            if all(rms(material - root) > MATERIAL_TOL for root in found):
                found.append(material)
    return found


def evaluate(budget, fresh_count, stress_count, cohort, tag, seed_offset, quadratic_seeds):
    results, arrays, scoring_labels = {}, {}, {}
    methods_to_search = (*METHODS, WEIGHTED) if quadratic_seeds else METHODS
    for j, (name, geometry, archived_obs, truth) in enumerate(cases(fresh_count, stress_count, cohort, seed_offset)):
        if truth is None:
            p = Experiment(geometry, archived_obs)
        else:
            generator = Model(**geometry)
            obs = generator.forward(truth, independent=True)["y"][0]
            p = Experiment(geometry, obs)
        # Only geometry, observations and the common starting materials enter discovery.
        sobol = LO + SPAN * qmc.Sobol(p.n, scramble=True, seed=421121 + j + seed_offset).random_base2(11)
        anchor, anchor_result = find_anchor(p, sobol, orthogonalize=quadratic_seeds)
        methods = {method: search(p, anchor, method, sobol, budget, quadratic_seeds)
                   for method in methods_to_search} if anchor is not None else {}
        results[name] = {"geometry": geometry, "common_anchor": anchor_result, "methods": methods}
        scoring_labels[name] = truth
        for key, value in (("B", p.b), ("H", p.h), ("R", p.r), ("y", p.obs)):
            arrays[f"{name}/{key}"] = value
        if anchor is not None:
            arrays[f"{name}/anchor"] = anchor
            for method, record in methods.items():
                arrays[f"{name}/{method}"] = np.array([branch["m"] for branch in record["new_branches"]]).reshape(-1, p.n)
        print(name, "anchor", anchor_result["starts"], " / branches",
              {m: len(v["new_branches"]) for m, v in methods.items()}, flush=True)

    for name, result in results.items():
        truth = scoring_labels[name]
        if truth is None:
            case = name.split("-", 1)[1]
            archived_roots = archive_references(case)
            with np.load(R5 / "private/truth.npz") as data:
                truth = data[case]
        else:
            archived_roots = []
        arrays[f"{name}/truth"] = truth
        if archived_roots:
            arrays[f"{name}/archive_references"] = np.array(archived_roots)
        result["archive_reference_count"] = len(archived_roots)
        if not result["methods"]:
            continue
        anchor = np.asarray(result["common_anchor"]["m"])
        result["common_anchor"]["truth_RMS"] = float(rms(anchor - truth))
        for record in result["methods"].values():
            roots = [anchor] + [np.asarray(root["m"]) for root in record["new_branches"]]
            record["best_truth_RMS"] = min(float(rms(m - truth)) for m in roots)
            record["archive_endpoints_recovered"] = sum(any(rms(m - root) <= MATERIAL_TOL for m in roots)
                                                       for root in archived_roots)
    base = ROOT / "results" / f"branch{tag}"
    base.with_suffix(".json").write_text(json.dumps({
        "per_method_cpu_seconds": budget, "box": [LO, HI], "response_relative_tolerance": RESPONSE_TOL,
        "root_dedup_RMS": MATERIAL_TOL, "case_seed_offset": seed_offset,
        "quadratic_seeds_in_guided_methods": quadratic_seeds, "cases": results,
    }, indent=2, allow_nan=False))
    np.savez_compressed(base.with_suffix(".npz"), **arrays)


def coverage(rows):
    eligible = [row for row in rows if row["methods"]]
    report = {
        "targets": len(rows), "anchors": len(eligible),
        "median_anchor_cpu_ms": 1000. * median(row["common_anchor"]["cpu_seconds"] for row in eligible),
        "truth_matches_at_anchor": sum(row["common_anchor"]["truth_RMS"] <= MATERIAL_TOL for row in eligible),
        "methods": {},
    }
    for method in eligible[0]["methods"]:
        records = [row["methods"][method] for row in eligible]
        first = [r["new_branches"][0]["cpu_seconds"] for r in records if r["new_branches"]]
        branches = [branch for r in records for branch in r["new_branches"]]
        report["methods"][method] = {
            "targets_with_extra_root": len(first), "extra_roots": len(branches),
            "median_cpu_ms_to_first_extra_root": 1000. * median(first) if first else None,
            "median_starts_completed": median(sum(not a.get("time_limited") for a in r["attempts"]) for r in records),
            "post_search_truth_matches": sum(r["best_truth_RMS"] <= MATERIAL_TOL for r in records),
            "archive_endpoints_matched": sum(r["archive_endpoints_recovered"] for r in records),
            "max_relative_data_residual_of_extra_roots": max((b["response_relative"] for b in branches), default=None),
            "extra_roots_using_material_polish": sum("pre_polish_residual" in r["attempts"][b["start"] - 1]
                                                   for r in records for b in r["new_branches"]),
        }
    return report


def figures():
    import matplotlib.pyplot as plt

    development = json.loads((ROOT / "results/branch.json").read_text())["cases"]
    archive = json.loads((ROOT / "results/branch_archive.json").read_text())["cases"]
    holdout = json.loads((ROOT / "results/branch_holdout.json").read_text())["cases"]
    diagnostic = json.loads((ROOT / "results/branch_diagnostics.json").read_text())["cases"]
    panels = {
        "R5 archive": list(archive.values()),
        "Held-out 6-cell": [v for k, v in holdout.items() if k.startswith("fresh6-")],
        "Held-out 12-cell": [v for k, v in holdout.items() if k.startswith("fresh12-")],
    }
    statistics = {title: coverage(rows) for title, rows in panels.items()}
    for size in (6, 12):
        statistics[f"Development {size}-cell"] = coverage([v for k, v in development.items()
                                                           if k.startswith(f"fresh{size}-")])
    (ROOT / "results/branch_summary.json").write_text(json.dumps(statistics, indent=2, allow_nan=False))

    methods = (*METHODS, WEIGHTED)
    colors = dict(zip(methods, ["#777b82", "#51914a", "#a24f48", "#15344e", "#927ab7", "#007ca9", "#dd811b"]))
    styles = dict(zip(methods, [":", "--", "--", "-", ":", "-", "-."]))
    names = dict(zip(methods, ["Multistart TRF", "Guided TRF", "Shifted deflation", "Ordinary secant",
                              "Joint reciprocal", "Projected reciprocal", "Radius-weighted reciprocal"]))
    with plt.rc_context({"font.size": 10, "axes.spines.right": False, "axes.spines.top": False}):
        fig, axes = plt.subplots(2, 3, figsize=(12, 6.8))
        for ax, (title, rows) in zip(axes[0], panels.items()):
            maximum = 0
            for method in methods:
                times = sorted(1000. * row["methods"][method]["new_branches"][0]["cpu_seconds"] for row in rows
                               if row["methods"][method]["new_branches"])
                maximum = max(maximum, len(times))
                ax.step([.2, *times, 320.], [0, *range(1, len(times) + 1), len(times)], where="post",
                        color=colors[method], ls=styles[method], lw=1.8, label=names[method])
            ax.set(xscale="log", xlim=(.2, 340), ylim=(-.12, maximum + .5), xticks=(.3, 1, 3, 10, 30, 100, 320),
                   xticklabels=("0.3", "1", "3", "10", "30", "100", "320"),
                   xlabel="Conditional search CPU (ms)", title=f"{title} ({len(rows)} targets)")
            ax.set_yticks(np.arange(0, maximum + 1, max(1, maximum // 5)))
            ax.grid(axis="y", color="#dddfe2", lw=.55)
        axes[0, 0].set_ylabel("Targets with an extra root")
        for ax, (case, title) in zip(axes[1], [("R5-T0", "R5 T0"), ("R5-T2", "R5 T2"),
                                             ("fresh12-00", "12-cell dev. 00")]):
            row = diagnostic[case][0]
            angles = np.array(list(row["basin_hits"]), dtype=float)
            for method, key, mark in [("ordinary secant", "ordinary", "o"), ("joint reciprocal", "joint", "s"),
                                      ("projected reciprocal", "projected", "^"), (WEIGHTED, "weighted", "d")]:
                rate = np.array([v[key] for v in row["basin_hits"].values()]) / row["basin_trials_per_angle"]
                ax.plot(angles, rate, color=colors[method], lw=1.5, ls=styles[method], marker=mark, ms=4)
            ax.set(xlim=(0, 31), ylim=(-.05, 1.06), xticks=angles,
                   title=f"{title}  |  cond(L) = {row['radial_L_condition']:.0f}",
                   xlabel="Perturbation of branch direction (°)")
            ax.grid(axis="y", color="#dddfe2", lw=.55)
        axes[1, 0].set_ylabel("Fraction returning to branch")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncols=4, frameon=False, fontsize=9.5)
        fig.subplots_adjust(left=.065, right=.985, bottom=.17, top=.94, hspace=.67, wspace=.28)
        fig.savefig(ROOT / "figures/branch_comparison.png", dpi=210)
        fig.savefig(ROOT / "figures/branch_comparison.pdf")
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=float, default=.32)
    parser.add_argument("--fresh", type=int, default=18)
    parser.add_argument("--stress", type=int, default=4)
    parser.add_argument("--cohort", choices=("all", "archive", "fresh"), default="all")
    parser.add_argument("--tag", default="")
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--quadratic-seeds", action="store_true")
    parser.add_argument("--figures", action="store_true")
    args = parser.parse_args()
    if args.figures:
        figures()
    else:
        evaluate(args.budget, args.fresh, args.stress, args.cohort, args.tag, args.seed_offset, args.quadratic_seeds)
