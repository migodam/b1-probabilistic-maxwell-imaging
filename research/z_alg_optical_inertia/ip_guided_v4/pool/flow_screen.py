"""Numerical proposal search for signed native 1024-cell L4 box screens."""

import json
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh, solve
from threadpoolctl import threadpool_limits

from pool_common import FLOW, clock, elapsed, write_json
from flow_audit import OUTPUT


def numeric_search(g, s, b, y, lo, hi):
    started = clock()
    m, r = (hi + lo)/2, (hi - lo)/2
    p = -m[:,None] * g
    p.flat[::len(m)+1] += 1.
    b_mid, b_rad = m[:,None] * b, r[:,None] * b
    rad_g = r[:,None] * g
    base_k = s.conj().T @ s
    base_l = s.conj().T @ y
    base_c = y.conj().T @ y
    shifted_k = p.conj().T @ p - rad_g.conj().T @ rad_g
    shifted_k = .5*(shifted_k+shifted_k.conj().T)
    shifted_l = p.conj().T @ b_mid + rad_g.conj().T @ b_rad
    shifted_c = b_mid.conj().T @ b_mid - b_rad.conj().T @ b_rad
    lossy = (g - g.conj().T)/(2j)
    scale = np.linalg.norm(s,2)**2
    result = []
    # Optical 0 probes the sector alone, both signs of 1 probe exact optical equalities.
    # Lower edge of the exact optical PSD interval was not used.
    for optical in (0., -.5, .5):
        for ratio in (2.**-10, 2.**-7, 2.**-4):
            tau = float(scale*ratio)
            k = base_k + tau*(shifted_k+optical*lossy)
            lin = base_l + tau*(shifted_l+(.5j*optical)*b)
            c = base_c + tau*shifted_c
            try:
                factor = cho_factor(k, lower=True, check_finite=False)
                x = cho_solve(factor, lin, check_finite=False)
                nu, edges = 0, []
            except np.linalg.LinAlgError:
                ev, vv = eigh(k, subset_by_index=(0,1), check_finite=False)
                edges = ev.tolist()
                if ev[1] < 0:
                    result.append(dict(tau=tau, optical=optical, lambda_low=edges,
                                       numeric_negative_inertia_at_least=2, numeric_lower=0.))
                    continue
                x = solve(k, lin, assume_a="her", check_finite=False)
                nu = 1
            h = c-lin.conj().T@x
            ev = np.maximum(0., np.linalg.eigvalsh(.5*(h+h.conj().T))[:2-nu])
            result.append(dict(tau=tau, optical=optical, lambda_low=edges,
                               numeric_negative_inertia=nu, numeric_lower=float(np.sqrt(ev.sum()))))
    return dict(weights=result, best=max(r["numeric_lower"] for r in result), **elapsed(started))


def main():
    begin = clock()
    data = json.loads((OUTPUT/"online.json").read_text())
    continuation = json.loads((OUTPUT/"continuation.json").read_text())["current_pass"]
    sys.path.insert(0, str(FLOW/"src"))
    import physics

    setup = clock()
    geometry = physics.Solver()
    setup_cost = elapsed(setup)
    receipts = {}
    for key in data["cases"]:
        row = continuation[key]
        with np.load(OUTPUT/f"{key}.npz", allow_pickle=False) as arrays:
            y, rxs, lo, hi = arrays["train"], arrays["used_rx"], arrays["lows"], arrays["highs"]
        gathered = {}
        for box in range(len(lo)):
            r = numeric_search(geometry.GD, geometry.GS[rxs], geometry.incident, y, lo[box], hi[box])
            r["online_continued_endpoint_residual"] = row["refined_box_results"][str(box)]["residual"]
            gathered[str(box)] = r
            print(f"{key} box {box:2d} numeric L4/r={r['best']/row['radius']:.4f}; "
                  f"endpoint/r={r['online_continued_endpoint_residual']/row['radius']:.3f}",flush=True)
            receipts[key] = dict(radius=row["radius"],boxes=gathered)
            write_json(OUTPUT/"screen_proposals.json",dict(cases=receipts, setup=setup_cost,
                total=elapsed(begin), status="numeric proposals only"))
    write_json(OUTPUT/"screen_proposals.json",dict(cases=receipts, setup=setup_cost,
        total=elapsed(begin), status="complete numeric proposals; exact proof required for any exclusion"))


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
