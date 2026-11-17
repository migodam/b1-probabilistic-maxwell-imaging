"""Independent 2D TM experiment for two-shot recovery of a prescribed cell."""
import os
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import argparse
import csv
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares, minimize
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
R2 = ROOT.parent / "research_z_alg_r2"
if not (R2 / "src/physics.py").is_file():
    R2 = ROOT.parents[3] / "z_alg_r2_publication/repository/research/research_z_alg_r2"
sys.path.insert(0, str(R2 / "src"))
from physics import Model, ReducedKernel

LOWER, UPPER = 0.5, 30.0
TARGET = 4
N_CASES = 24
NOISE_PER_SHOT = (0.0, 1e-6, 1e-4, 5e-4, 2e-3)
CLEAN_EPS = 3e-12
APERTURES = ("ring", "arc", "wide_ring", "wide_arc")
SEED = 9172626


def geometry(aperture):
    """Cell accounting and aperture prescription independent of unknown materials."""
    j, i = np.indices((7, 7))
    pitch = .30 if aperture.startswith("wide") else .145
    cells = np.column_stack((pitch * (i.ravel() - 3), pitch * (j.ravel() - 3)))
    active = (7 * j[1::2, 1::2] + i[1::2, 1::2]).ravel()
    x, y = cells.T * (.145 / pitch)
    host = 2.2 + .8 * np.cos(4 * x + y) + 1.6 * np.exp(-((x + .23)**2 + (y - .18)**2) / .11)
    host[active] = 0.
    if aperture.endswith("ring"):
        angles = np.linspace(-np.pi, np.pi, 20, endpoint=False) + .071
    else:
        angles = np.pi + np.deg2rad(np.linspace(-75, 75, 20))
    sources = 1.55 * np.column_stack((np.cos(angles), np.sin(angles)))
    receiver = np.array([[.38, 1.10]])
    return cells, active, host, sources, receiver


def design(aperture):
    cells, active, host, sources, receiver = geometry(aperture)
    start = time.process_time()
    model = Model(centers=cells, sources=sources, receivers=receiver, radius=.041, k=2*np.pi)
    kernel = ReducedKernel(model, host, active=active)
    host_seconds = time.process_time() - start

    start = time.process_time()
    H, B = kernel.K, kernel.q
    Q = j0(model.k * np.linalg.norm(sources[:, None] - sources[None, :], axis=2))
    eigenvalues, eigenvectors = np.linalg.eigh(Q)
    # On this physical subspace, any unit-radiated-power feed satisfies ||a||_2 <= 4.
    keep = eigenvalues >= 1 / 16
    P = eigenvectors[:, keep] / np.sqrt(eigenvalues[keep])
    U, D, Vh = np.linalg.svd(B @ P, full_matrices=False)
    rank = int(np.count_nonzero(D > 1e-12 * D[0]))
    directions = np.column_stack((np.eye(len(active))[:, TARGET], H[:, TARGET]))
    U = U[:, :rank]
    V = Vh[:rank].conj().T
    # Whitened minimum-power preimages in the allowed physical source subspace.
    preimages = V @ ((U.conj().T @ directions) / D[:rank, None])
    basis, _, _ = np.linalg.svd(preimages, full_matrices=False)
    focus = P @ basis
    coordinates = basis.conj().T @ preimages
    synthesis = B @ focus @ coordinates - directions
    relative_synthesis = np.linalg.norm(synthesis, axis=0) / np.linalg.norm(directions, axis=0)
    available = bool(np.max(relative_synthesis) < 1e-10)
    full = P @ V if rank == len(active) else None
    source_seconds = time.process_time() - start

    meta = dict(
        aperture=aperture, cells=cells, active=active, host=host, sources=sources, receiver=receiver,
        radius=model.radius, wavenumber=model.k, target_active_index=TARGET,
        target_cell_index=active[TARGET], Q_eigenvalues=eigenvalues, Q_eigenvalue_floor=1/16,
        Q_retained_modes=int(np.count_nonzero(keep)), incident_rank=rank, incident_singular_values=D,
        design_has_exact_span_to_roundoff=available, relative_synthesis_error=relative_synthesis,
        max_Q_orthonormality_error=np.linalg.norm(focus.conj().T @ Q @ focus - np.eye(2)),
        focus_component_minimum_powers=np.sum(abs(preimages)**2, axis=0),
        focus_feed_l2_norms=np.linalg.norm(focus, axis=0),
        host_setup_cpu_seconds=host_seconds, host_setup_counts=dict(model.counts),
        source_design_cpu_seconds=source_seconds, focus=focus, focus_coordinate_matrix=coordinates,
        full_completion_condition=(D[0] / D[-1] if full is not None else None)
    )
    return model, kernel, Q, P, focus, coordinates, full, available, meta


def pack(y):
    return np.column_stack((y.real, y.imag)).reshape(-1)


class Response:
    """Reduced single-receiver model with cached analytic material Jacobian."""
    def __init__(self, kernel, feeds):
        self.H = kernel.K
        self.R = kernel.L.ravel()
        self.F = kernel.q @ feeds
        self.identity = np.eye(kernel.r)
        self.z = None
        self.nlu, self.nrhs = 0, 0
        self.jacobian = None

    def predict(self, z, derivative=False):
        if self.z is None or not np.array_equal(self.z, z):
            self.z = z.copy()
            self.factor = lu_factor(self.identity - z[:, None] * self.H, check_finite=False)
            self.t = lu_solve(self.factor, self.R, trans=1, check_finite=False)
            self.w = self.t * z
            self.y = self.w @ self.F
            self.jacobian = None
            self.nlu += 1
            self.nrhs += 1
        if derivative and self.jacobian is None:
            current = lu_solve(self.factor, z[:, None] * self.F, check_finite=False)
            field = self.F + self.H @ current
            j = (self.t[:, None] * field).T
            self.jacobian = np.stack((j.real, j.imag), axis=1).reshape(-1, len(z))
            self.nrhs += self.F.shape[1]
        return self.jacobian if derivative else self.y


def bounded_noise(rng, shots):
    # Uniform in the ball in 2s real dimensions, subsequently scaled by sqrt(s)*sigma.
    n = rng.standard_normal(shots) + 1j * rng.standard_normal(shots)
    return n * rng.uniform()**(1 / (2 * shots)) / np.linalg.norm(n)


def confidence(d, R, coords, epsilon):
    """Real q satisfying |u-qv| <= epsilon*||c_e-q*c_h|| on the material box."""
    start = time.process_time()
    ce, ch = coords.T
    u, v = d @ ce, R + d @ ch
    norm_v_sq = abs(v)**2
    q0 = np.clip(np.real(np.conj(v) * u) / norm_v_sq, LOWER, UPPER) if norm_v_sq else (LOWER+UPPER)/2
    residual0, feed0 = u - q0*v, ce - q0*ch
    a = norm_v_sq - epsilon**2 * np.vdot(ch, ch).real
    b = -2 * np.real(np.conj(v) * residual0 - epsilon**2 * np.vdot(ch, feed0))
    c = abs(residual0)**2 - epsilon**2 * np.vdot(feed0, feed0).real
    # The shift avoids subtracting O(q^2) coefficients to measure a tight interval.
    if a == 0:
        roots = [] if b == 0 else [-c / b]
    else:
        discriminant = b*b - 4*a*c
        if discriminant < 0:
            roots = []
        else:
            stable = -.5 * (b + np.copysign(np.sqrt(discriminant), b))
            roots = [-b / (2*a)] if stable == 0 else [stable/a, c/stable]
    internal = sorted({float(q0 + r) for r in roots if LOWER <= q0 + r <= UPPER})
    endpoints = [LOWER, *internal, UPPER]
    segments = []
    for lo, hi in zip(endpoints, endpoints[1:]):
        x = .5 * (lo+hi) - q0
        if (a*x+b)*x+c <= 0:
            segments.append([lo, hi])
    for root in internal:
        if not any(lo <= root <= hi for lo, hi in segments):
            segments.append([root, root])
    segments.sort()
    merged = []
    for lo, hi in segments:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = hi
        else:
            merged.append([lo, hi])
    interval = [merged[0][0], merged[-1][1]] if merged else None
    return dict(
        interval=interval, components=merged,
        width=(interval[1]-interval[0] if interval else None),
        estimate=(float(np.mean(interval)) if interval else None),
        real_ls=q0, observed_denominator_abs=abs(v),
        numerically_uninformative=(interval == [LOWER, UPPER]),
        cpu_seconds=time.process_time()-start
    )


def trf(response, d, initials):
    begin = time.process_time()
    lu_before, rhs_before = response.nlu, response.nrhs
    scale = max(np.linalg.norm(d), 1e-3)
    rows = []
    for m in initials:
        fit = least_squares(
            lambda z: pack(response.predict(z)-d)/scale, m,
            jac=lambda z: response.predict(z, derivative=True)/scale,
            bounds=(LOWER, UPPER), max_nfev=100,
            ftol=3e-12, xtol=3e-12, gtol=3e-12,
        )
        rows.append(dict(z=fit.x, residual=np.linalg.norm(fit.fun)*scale, status=fit.status,
                         nfev=fit.nfev, njev=fit.njev))
    best = min(rows, key=lambda row: row["residual"])
    return dict(
        estimate=best["z"][TARGET], best=best, starts=rows,
        target_spread=float(np.ptp([row["z"][TARGET] for row in rows])),
        nlu=response.nlu-lu_before, nrhs=response.nrhs-rhs_before,
        cpu_seconds=time.process_time()-begin,
    )


def profiles(response, d, fit, epsilon, independent, feeds):
    """Search two feasible target extrema and verify selected witnesses with full VIE."""
    begin = time.process_time()
    lu_before, rhs_before = response.nlu, response.nrhs
    scale = max(np.linalg.norm(d), 1e-3)
    epsilon_scaled = epsilon / scale

    def constraint(z):
        r = pack(response.predict(z)-d) / scale
        return epsilon_scaled - np.linalg.norm(r)

    def constraint_jac(z):
        r = pack(response.predict(z)-d) / scale
        J = response.predict(z, derivative=True) / scale
        norm = np.linalg.norm(r)
        return -((r / norm) @ J) if norm else np.zeros(len(z))

    best = fit["best"]["z"]
    feasible = [row["z"] for row in fit["starts"] if row["residual"] <= epsilon]
    trials = []
    target = np.eye(len(best))[TARGET]
    if feasible and epsilon > 10 * CLEAN_EPS:
        for sign in (1, -1):
            warm = min(feasible, key=lambda z: sign * z[TARGET])
            trial = minimize(
                lambda z: sign*z[TARGET], warm, jac=lambda _: sign * target,
                bounds=[(LOWER, UPPER)] * len(best), method="SLSQP",
                constraints=[dict(type="ineq", fun=constraint, jac=constraint_jac)],
                options=dict(maxiter=90, ftol=2e-11),
            )
            residual = np.linalg.norm(response.predict(trial.x)-d)
            # A small retreat toward an already feasible point tolerates SLSQP
            # returning the just-outside side of the constraint surface.
            if residual > epsilon:
                for frac in (1e-6, 1e-4, 1e-2):
                    candidate = (1-frac)*trial.x + frac*best
                    candidate_residual = np.linalg.norm(response.predict(candidate)-d)
                    if candidate_residual <= epsilon:
                        feasible.append(candidate)
                        break
            else:
                feasible.append(trial.x)
            trials.append(dict(sign=sign, z=trial.x, residual=residual,
                               status=trial.status, iterations=trial.nit))
    extrema = []
    for sign in (1, -1):
        for z in sorted(feasible, key=lambda a: sign*a[TARGET]):
            residual = np.linalg.norm(independent(z) @ feeds-d)
            if residual <= epsilon:
                extrema.append(dict(z=z, residual_l2=residual, target=z[TARGET]))
                break
    return dict(
        independently_validated_witnesses=extrema,
        target_spread=(abs(extrema[1]["target"]-extrema[0]["target"]) if len(extrema) == 2 else None),
        search_trials=trials, scope="demonstrated feasible spread, lower bound on diameter",
        nlu=response.nlu-lu_before, nrhs=response.nrhs-rhs_before,
        cpu_seconds=time.process_time()-begin,
    )


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist() if not np.iscomplexobj(value) else {"real":value.real.tolist(), "imag":value.imag.tolist()}
    if isinstance(value, complex):
        return dict(real=value.real, imag=value.imag)
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def run():
    start = time.process_time()
    models = {key: design(key) for key in APERTURES}
    for key, data in models.items():
        meta = data[-1]
        print(key, "ranks", meta["Q_retained_modes"], meta["incident_rank"], "synthesis",
              meta["relative_synthesis_error"], "spectrum", meta["incident_singular_values"], flush=True)
    rng_material = np.random.default_rng(SEED)
    truth = rng_material.uniform(LOWER, UPPER, (N_CASES, len(models["ring"][1].active)))
    starts_rng = np.random.default_rng(SEED+1)
    starts = np.vstack((np.full(truth.shape[1], (LOWER+UPPER)/2),
                        starts_rng.uniform(LOWER, UPPER, (11, truth.shape[1]))))
    observations, rows, cases, realized = [], [], [], {}
    for aperture, data in models.items():
        model, kernel, Q, P, focus, coords, full, available, meta = data
        R = kernel.L.ravel()
        material_rows = []
        host_y = kernel.y0.ravel()
        active = kernel.active
        start_data = time.process_time()

        def independent(z):
            m = kernel.m0.copy()
            m[active] += z
            return model.forward(m, independent=True)["y"].ravel() - host_y

        for case, z in enumerate(truth):
            rng_random = np.random.default_rng(SEED+100+case)
            random_q = np.linalg.qr(
                rng_random.standard_normal((P.shape[1], 2)) +
                1j*rng_random.standard_normal((P.shape[1], 2))
            )[0]
            feeds = dict(focused=focus, random=P @ random_q)
            if full is not None:
                feeds["completion_equal_energy"] = full * np.sqrt(2 / len(active))
                feeds["completion_equal_shot"] = full
                feeds["completion_equal_noise_energy"] = feeds["completion_equal_energy"]
            designs = {arm: Response(kernel, A) for arm, A in feeds.items()}
            clean_dictionary = independent(z)
            full_state = model.state
            first_born = (R*z) @ kernel.q
            w_engine = Response(kernel, np.eye(model.ns))
            clean_reduced = w_engine.predict(z)
            reciprocal = (R + w_engine.w @ kernel.K)[TARGET]
            target_identity_error = abs(w_engine.w[TARGET] - z[TARGET]*reciprocal)
            material_rows.append(dict(
                case=case, z=z, clean_library_l2=np.linalg.norm(clean_dictionary),
                relative_Born_response_disagreement=np.linalg.norm(clean_dictionary-first_born)/np.linalg.norm(clean_dictionary),
                max_full_reduced_library_error=np.max(abs(clean_reduced-clean_dictionary)),
                target_clean_reciprocal_abs=abs(reciprocal), target_identity_error=target_identity_error,
                full_VIE_relative_residual=full_state["state_relative_residual"],
                full_VIE_sigma_min=full_state["resolvent_sigma_min"]
            ))
            pairs = {}
            completion_coords = {}
            for arm, A in feeds.items():
                realized[f"{aperture}_{case:02d}_{arm}_feeds"] = A
                pairs[arm] = (clean_dictionary @ A, np.diag(A.conj().T @ Q @ A).real)
                if arm.startswith("completion"):
                    completion_coords[arm] = np.linalg.solve(
                        kernel.q @ A, np.column_stack((np.eye(len(active))[:, TARGET], kernel.K[:, TARGET]))
                    )
            for level, sigma in enumerate(NOISE_PER_SHOT):
                noise_rng = np.random.default_rng(SEED+4000+case*len(NOISE_PER_SHOT)+level)
                noise_two = bounded_noise(noise_rng, 2)
                noise_full = bounded_noise(noise_rng, len(active))
                for arm, A in feeds.items():
                    clean, powers = pairs[arm]
                    noise_dimension = 2 if arm == "completion_equal_noise_energy" else len(clean)
                    noise_radius = np.sqrt(noise_dimension) * sigma
                    epsilon = noise_radius + CLEAN_EPS
                    noise = noise_radius * (noise_two if len(clean)==2 else noise_full)
                    d = clean + noise
                    fit = trf(designs[arm], d, starts) if arm != "completion_equal_noise_energy" else None
                    witnesses = None
                    if arm in ("focused", "random"):
                        witnesses = profiles(designs[arm], d, fit, epsilon, independent, A)
                    inversion = None
                    direct = None
                    if arm == "focused" and available:
                        inversion = confidence(d, R[TARGET], coords, epsilon)
                    elif arm.startswith("completion"):
                        dstart = time.process_time()
                        w = np.linalg.solve((kernel.q @ A).T, d)
                        quotient = w / (R + w @ kernel.K)
                        direct = dict(estimate=np.clip(quotient[TARGET].real, LOWER, UPPER),
                                      quotient=quotient, cpu_seconds=time.process_time()-dstart)
                        inversion = confidence(d, R[TARGET], completion_coords[arm], epsilon)
                    if inversion is not None:
                        inversion["covers_truth"] = any(lo <= z[TARGET] <= hi for lo, hi in inversion["components"])
                        inversion["error"] = abs(inversion["estimate"] - z[TARGET]) if inversion["estimate"] is not None else None
                    row = dict(
                        aperture=aperture, case=case, sigma=sigma, arm=arm,
                        truth=z[TARGET], shots=len(clean), shot_powers=powers,
                        total_emitted_energy=float(np.sum(powers)), radius=epsilon, actual_noise_l2=np.linalg.norm(noise),
                        clean_response_l2=np.linalg.norm(clean), interval_inversion=inversion,
                        trf=fit, feasible_profile=witnesses, direct_quotient=direct
                    )
                    observations.append(row)
                    rows.append(dict(
                        aperture=aperture, case=case, sigma=sigma, arm=arm, truth=z[TARGET],
                        shots=len(clean), total_energy=np.sum(powers),
                        bounded_SNR_dB=(20*np.log10(np.linalg.norm(clean)/noise_radius) if sigma else None),
                        interval_lo=(inversion["interval"][0] if inversion and inversion["interval"] else None),
                        interval_hi=(inversion["interval"][1] if inversion and inversion["interval"] else None),
                        interval_width=(inversion["width"] if inversion else None),
                        covers_truth=(inversion["covers_truth"] if inversion else None),
                        uninformative=(inversion["numerically_uninformative"] if inversion else None),
                        interval_estimate=(inversion["estimate"] if inversion else None),
                        interval_error=(inversion["error"] if inversion else None),
                        interval_cpu_s=(inversion["cpu_seconds"] if inversion else None),
                        trf_estimate=(fit["estimate"] if fit else None),
                        trf_error=(abs(fit["estimate"]-z[TARGET]) if fit else None),
                        trf_best_residual=(fit["best"]["residual"] if fit else None),
                        trf_cpu_s=(fit["cpu_seconds"] if fit else None),
                        trf_reduced_lu=(fit["nlu"] if fit else None),
                        trf_reduced_rhs=(fit["nrhs"] if fit else None),
                        witness_spread=(witnesses["target_spread"] if witnesses else None),
                        direct_error=(abs(direct["estimate"]-z[TARGET]) if direct else None),
                    ))
                    realized[f"{aperture}_{case:02d}_{arm}_noise{level}_observed"] = d
            if (case+1)%6 == 0:
                print(aperture, case+1, "of", len(truth), "CPU elapsed", round(time.process_time()-start_data, 1), flush=True)
        cases.append(dict(aperture=aperture, physics_rows=material_rows, model_operations=dict(model.counts),
                          model_wall_seconds=dict(model.timings)))
    out = ROOT / "results"
    with (out/"physical_individual_rows.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (out/"physical_raw.json").write_text(json.dumps(
        dict(seed=SEED, contrast_bounds=[LOWER, UPPER], numeric_margin=CLEAN_EPS, starts=starts,
             design={key: val[-1] for key,val in models.items()}, physics_checks=cases,
             observations=observations, cpu_seconds=time.process_time()-start),
        default=encode, separators=(",", ":"))+"\n")
    np.savez_compressed(out/"physical_observations.npz", truths=truth, **realized)
    print("Completed", len(rows), "observations in", round(time.process_time()-start, 1), "CPU seconds", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-only", action="store_true")
    args = parser.parse_args()
    if args.design_only:
        for case in APERTURES:
            meta = design(case)[-1]
            print(json.dumps({key: value for key, value in meta.items() if key not in
                              ("cells", "sources", "focus", "focus_coordinate_matrix")}, default=encode))
    else:
        run()
