"""Bound the capped four-region SDP from below using only its receivers and trace cap."""
from pathlib import Path
import json
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY / "submission/curve"))
from curve_inputs import full_ring_inputs
sys.path.insert(0, str(STUDY))
from core.certificate import ball, gershgorin, safe_float, star

GRID = (0., 10., 100., 250., 350., 400., 1e3, 1e4, 1e5, 1e6, 1e8)


def eigen_interval(form):
    mid = np.array(form.tolist(), dtype=complex)
    _, vectors = eigh((mid + mid.conj().T) / 2)
    trial = ball(vectors[:, -1:])
    rayleigh = ((star(trial) * form * trial)[0, 0] /
                (star(trial) * trial)[0, 0]).real
    basis = ball(vectors)
    floor = gershgorin(star(basis) * basis, upper=False)
    if not floor > 0:
        raise ArithmeticError("The source basis could not be verified.")
    upper = gershgorin(star(basis) * form * basis, upper=True) / floor
    return rayleigh, upper


def solve(geometry, ports, energy):
    y = ball(ports)
    norm_lo, norm_hi = eigen_interval(star(y) * y)
    if energy == 0:
        return dict(cap=0., lower=safe_float(norm_lo.sqrt(), upper=False),
                    upper=safe_float(norm_hi.sqrt(), upper=True), exact_law="operator norm of Y")
    gram, ev, left = geometry
    rotated = left.conj().T @ ports

    def neg(log_eta):
        eta = np.exp(log_eta)
        inner = rotated.conj().T @ (rotated / (ev + eta)[:, None])
        return -eta * (eigh(inner, eigvals_only=True)[-1] - energy)

    proposal = minimize_scalar(neg, bounds=(-34., 10.), method="bounded")
    eta = arb(float(np.exp(proposal.x)))
    weighted = gram + acb_mat(np.eye(gram.nrows(), dtype=int).tolist()) * eta
    h = star(y) * weighted.solve(y)
    rayleigh, _ = eigen_interval(h)
    lower_sq = eta * (rayleigh - arb(energy))
    lower = safe_float(lower_sq.sqrt(), upper=False) if lower_sq > 0 else 0.
    return dict(cap=energy, lower=max(0., lower), upper=safe_float(norm_hi.sqrt(), upper=True),
                receiver_eta=float(eta.mid()), upper_witness="zero current")


def main():
    wall, cpu = perf_counter(), process_time()
    _, sensor, _, observations = full_ring_inputs()
    ctx.prec = 128
    s = ball(sensor)
    gram = s * star(s)
    ev, left = eigh(sensor @ sensor.conj().T)
    geometry = gram, ev, left
    result = {}
    for obs in observations:
        if obs.object_id not in ("optfresh_01", "optfresh_05") or obs.noise_pct != 2:
            continue
        needed = star(ball(obs.y)) * gram.solve(ball(obs.y))
        lower, upper = eigen_interval(needed)
        result[obs.object_id] = {
            "observation": obs.id, "radius": obs.padded_radius,
            "necessary_energy_for_exact_fit_all_directions": [
                safe_float(lower, upper=False), safe_float(upper, upper=True)],
            "rows": [solve(geometry, obs.y, energy) for energy in GRID],
        }
    (HERE / "energy_port_control.json").write_text(json.dumps({
        "scope": "Uniform coherent maximum of the per-direction four-region minima. Lower bounds drop both local material constraints. Upper bounds use physical zero current.",
        "energy": "trace of current-current block; Euclidean cell convention",
        "precision_bits": ctx.prec, "objects": result, "wall_seconds": perf_counter() - wall,
        "cpu_seconds": process_time() - cpu}, indent=2) + "\n")


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
