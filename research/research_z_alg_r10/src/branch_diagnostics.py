"""Reproduce the retrospective direction-basin and local Jacobian diagnostics."""

import json

import branch_discovery as b
import numpy as np
from scipy.optimize._numdiff import approx_derivative


CASES = (
    "R5-T0", "R5-T2", "R5-T4", "fresh6-00", "fresh6-15",
    "fresh12-00", "fresh12-01", "fresh12-02",
)
ANGLES = (1, 3, 7, 15, 30)
TRIALS_PER_ANGLE = 12


def diagnose():
    results_dir = b.ROOT / "results"
    summary = json.loads((results_dir / "branch.json").read_text())
    rng = np.random.default_rng(421136)
    results = {
        "note": "Retrospective perturbations of independently accepted ordinary-secant roots; "
                "these roots were never given to the discovery methods.",
        "cases": {},
    }
    with np.load(results_dir / "branch.npz") as data:
        for name in CASES:
            row = summary["cases"][name]
            p = b.Experiment(row["geometry"], data[name + "/y"])
            anchor = data[name + "/anchor"]
            secant = b.Secant(p, anchor, b.no_deadline)
            seed_index = list(summary["cases"]).index(name)
            sobol = b.qmc.Sobol(p.n, scramble=True, seed=421121 + seed_index).random_base2(7)
            plain = b.LO + b.SPAN * sobol
            inits = np.array([secant.start(m) for m in b.guided_seeds(p, anchor, plain)[:16]])
            fold = secant.quadratic_seeds()
            a = secant.q * secant.v.conj()
            v = secant.w * secant.qh.conj()
            jac = np.concatenate(((a + v).imag, (a - v).real)).T
            _, jac_singular, jac_vh = np.linalg.svd(jac)
            rows = []
            for root in data[name + "/ordinary secant"]:
                direction = secant.start(root)
                norm_direction = np.linalg.norm(direction)
                q = direction / norm_direction
                n = q @ secant.q
                h = q @ secant.qh
                a = np.column_stack(
                    ((n * secant.v.conj() + secant.w * h.conj()).imag,
                     (n * secant.v.conj() - secant.w * h.conj()).real)
                )
                c = (n * h.conj()).imag
                singular = np.linalg.svd(a, compute_uv=False)
                s = -np.linalg.lstsq(a, c, rcond=1e-12)[0]
                seed_angles = np.degrees(
                    np.arccos(np.clip(abs(inits @ q.conj()) / np.linalg.norm(inits, axis=1), 0, 1))
                )
                fold_angles = [
                    float(np.degrees(np.arccos(np.clip(abs(fq @ q.conj()) / np.linalg.norm(fq), 0, 1))))
                    for _, fq in fold
                ]
                weak_direction_coefficients = np.abs(jac_vh @ np.r_[q.real, q.imag])
                relative_radius = norm_direction / secant.step_scale
                rows.append({
                    "material_RMS": float(b.rms(root - anchor)),
                    "relative_secant_size": float(relative_radius),
                    "radial_L_condition": float(singular[0] / singular[-1]),
                    "radial_L_singulars": singular.tolist(),
                    "projective_zero_residual": float(np.linalg.norm(a @ s + c) / secant.reciprocal_scale),
                    "exact_angular_coefficients_on_linearized_singular_basis": weak_direction_coefficients.tolist(),
                    "closest_of_first_16_seed_angles_deg": float(seed_angles.min()),
                    "fold_seed_angles_deg": fold_angles,
                    "radial_reciprocal_weight": float(.05 + .95 / (1 + 1 / relative_radius**2)),
                    "linearized_ordinary_singulars": jac_singular.tolist(),
                    "jacobian_checks_at_root": check_jacobians(secant, direction),
                    "basin_trials_per_angle": TRIALS_PER_ANGLE,
                    "basin_hits": angular_basin(secant, root, q, norm_direction, rng),
                })
            results["cases"][name] = rows
            print(name, "nearest seed angle", rows[0]["closest_of_first_16_seed_angles_deg"],
                  "fold seed angles", rows[0]["fold_seed_angles_deg"],
                  "basin hits", rows[0]["basin_hits"], flush=True)
    (results_dir / "branch_diagnostics.json").write_text(json.dumps(results, indent=2))


def check_jacobians(secant, direction):
    original = b.least_squares
    captures = []

    def capture(fun, x, *, jac, **kwargs):
        value = fun(x)
        analytic = jac(x)
        finite_diff = approx_derivative(fun, x, method="3-point", abs_step=1e-6)
        captures.append({
            "residual_norm": float(np.linalg.norm(value)),
            "jac_relative_error": float(np.linalg.norm(analytic - finite_diff) / np.linalg.norm(finite_diff)),
            "jac_singular_values": np.linalg.svd(analytic, compute_uv=False).tolist(),
        })
        return original(fun, x, jac=jac, **kwargs)

    methods = (
        ("ordinary", lambda: secant.ordinary(direction)),
        ("joint", lambda: secant.reciprocal(direction, False)),
        ("projected", lambda: secant.reciprocal(direction, True)),
        ("radius-weighted", lambda: secant.reciprocal(direction, True, True)),
    )
    derivatives = {}
    b.least_squares = capture
    try:
        for title, solve in methods:
            try:
                solve()
            except (ValueError, np.linalg.LinAlgError):
                continue
            derivatives[title] = captures[-1]
    finally:
        b.least_squares = original
    return derivatives


def angular_basin(secant, root, q, norm_direction, rng):
    basin = {}
    for angle in ANGLES:
        hits = dict(ordinary=0, joint=0, projected=0, weighted=0)
        theta = np.radians(angle)
        for _ in range(TRIALS_PER_ANGLE):
            dq = rng.normal(size=secant.d) + 1j * rng.normal(size=secant.d)
            dq -= q * (dq @ q.conj())
            dq /= np.linalg.norm(dq)
            init = norm_direction * (np.cos(theta) * q + np.sin(theta) * dq)
            solutions = (
                ("ordinary", secant.ordinary(init)),
                ("joint", secant.reciprocal(init, False)),
                ("projected", secant.reciprocal(init, True)),
                ("weighted", secant.reciprocal(init, True, True)),
            )
            for key, (delta, _) in solutions:
                material = secant.recover(delta)[0] if delta is not None else None
                if material is not None and b.rms(material - root) < b.MATERIAL_TOL:
                    hits[key] += 1
        basin[str(angle)] = hits
    return basin


if __name__ == "__main__":
    diagnose()
