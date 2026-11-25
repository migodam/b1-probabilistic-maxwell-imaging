"""Rigorous selected endpoint and perturbation checks for signed native A13 boxes."""

import argparse
import json
import sys

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits

from pool_common import FLOW, FullWave, clock, elapsed, write_arrays, write_json
from box_certificates import ball, identity, lower, midpoint, norm, spectral_upper, star, upper
from flow_audit import OUTPUT
from replay import exact_surjectivity

PROPOSAL_RANK = 120
SEED = 624
TARGET_PRECONDITIONER_DEFECT = .98


def row_scale(matrix, values):
    cols = matrix.ncols()
    entries = matrix.entries()
    return acb_mat(matrix.nrows(), cols,
                   [p*value for i, value in enumerate(values) for p in entries[i*cols:(i+1)*cols]])


def column_scale(matrix, values):
    cols = matrix.ncols()
    return acb_mat(matrix.nrows(), cols,
                   [p*values[i%cols] for i, p in enumerate(matrix.entries())])


def bound_norm(a):
    return norm(a).upper()


class SignedNative:
    """Ball proof via a low-rank preconditioner of the physical matrix I-XG."""

    def __init__(self, green, sensor, incident):
        start = clock()
        self.g, self.s, self.b = map(ball, (green, sensor, incident))
        rng = np.random.default_rng(SEED)
        projection = green @ rng.standard_normal(size=(green.shape[1], PROPOSAL_RANK))
        projection = green @ (green.conj().T @ projection)
        u = np.linalg.qr(projection, mode="reduced")[0]
        v = u.conj().T @ green
        self.u, self.v = ball(u), ball(v)
        self.i = identity(PROPOSAL_RANK)
        self.e = self.g - self.u*self.v
        self.u_norm = spectral_upper(self.u, u)
        self.v_norm = spectral_upper(star(self.v), v.conj().T)
        self.s_norm = spectral_upper(star(self.s), sensor.conj().T)
        self.e_norm = norm(self.e).upper()
        self.g_norm = (self.u_norm*self.v_norm + self.e_norm).upper()
        width = self.e.ncols()
        values = self.e.entries()
        self.square_by_row = [sum((x.real*x.real+x.imag*x.imag
                                   for x in values[i*width:(i+1)*width]), arb(0)).upper()
                              for i in range(self.e.nrows())]
        self.setup = dict(approximation_rank=PROPOSAL_RANK, seed=SEED,
                          norm_green_upper=upper(self.g_norm),
                          approximation_defect_frob_upper=upper(self.e_norm), **elapsed(start))

    def inverse(self, x):
        mu = row_scale(self.u, x)
        k = (self.i-self.v*mu).solve(self.i)
        left = mu*k
        lnorm = bound_norm(left)
        t = (1+lnorm*self.v_norm).upper()
        me = sum((v*v*sq for v,sq in zip(x, self.square_by_row, strict=True)), arb(0)).sqrt().upper()
        e = (t*me).upper()
        used_targeted_defect = False
        if e >= TARGET_PRECONDITIONER_DEFECT:
            # Exploit that the large inverse modes see only V X (G-U V).
            overlap = bound_norm(column_scale(self.v, x)*self.e)
            e = min(e, (me+lnorm*overlap).upper())
            used_targeted_defect = True
        if not e < 1:
            raise ArithmeticError(f"Unresolved physical inverse; defect upper is {upper(e):.6g}.")
        return k, left, t, e, (t/(1-e)).upper(), used_targeted_defect

    def point(self, material, state, observed):
        start = clock()
        x = [arb(float(value)) for value in material]
        _, _, _, e, inv, targeted = self.inverse(x)
        j = ball(state)
        res = bound_norm(j-row_scale(self.b+self.g*j, x))
        port_diff = bound_norm(self.s*j-ball(observed))
        err = (port_diff + self.s_norm*inv*res).upper()
        return dict(strict_residual_upper=upper(err), approximate_observed_residual_upper=upper(port_diff),
                    state_defect_upper=upper(res), inverse_defect_upper=upper(e),
                    inverse_norm_upper=upper(inv), used_targeted_defect=targeted, **elapsed(start))

    def contraction(self, lo, hi, qvec, state, y):
        start = clock()
        endpoints = [(arb(float(a)), arb(float(b))) for a,b in zip(lo, hi, strict=True)]
        center = [(a+b)/2 for a,b in endpoints]
        halfwidth = [(b-a)/2 for a,b in endpoints]
        d_norm = max(v.upper() for v in halfwidth)
        k, left, t, e, inv, targeted = self.inverse(center)
        # First column(s) are X0 B forward currents; last column is the singular witness.
        j, v, ports = ball(state[:, :-1]), ball(state[:, -1:]), ball(y)
        witness = ball(qvec[:, None])
        field_v = self.g*v
        check = v-row_scale(field_v, center)-row_scale(witness, halfwidth)
        v_defect = bound_norm(check)
        kappa_low = (norm(field_v).lower()-self.g_norm*inv*v_defect)/norm(witness).upper()

        inverse_change = (t*e/(1-e)).upper()
        kvr = k * column_scale(self.v, halfwidth)
        kv_norm = spectral_upper(star(kvr), midpoint(kvr).conj().T)
        kappa_hi = (self.u_norm*kv_norm + self.e_norm*t*d_norm +
                    self.g_norm*inverse_change*d_norm).upper()
        standard_lower = 0.
        if kappa_hi < 1:
            approx_center_field = self.b+self.g*j
            state_diff = inv*bound_norm(j-row_scale(approx_center_field, center))
            ec = bound_norm(approx_center_field) + self.g_norm*state_diff
            center_residual_lower = norm(self.s*j-ports).lower()-self.s_norm*state_diff
            st = self.s + (self.s*left)*self.v
            strad = column_scale(st, halfwidth)
            optical_radius = (spectral_upper(star(strad), midpoint(strad).conj().T)
                              + self.s_norm*inverse_change*d_norm)*ec/(1-kappa_hi)
            standard_lower = max(0., lower(center_residual_lower-optical_radius))
        return dict(strict_kappa_lower=max(0.,lower(kappa_low)), strict_kappa_upper=upper(kappa_hi),
                    strict_standard_lower=standard_lower, inverse_defect_upper=upper(e),
                    used_targeted_defect=targeted, **elapsed(start))


def point_runner(cases, verifier, solver):
    continuation = json.loads((OUTPUT/"continuation.json").read_text())["current_pass"]
    receipts = {}
    with np.load(OUTPUT/"continuation_arrays.npz", allow_pickle=False) as arrays:
        for name, row in cases.items():
            native = continuation[name]
            with np.load(OUTPUT/f"{name}.npz", allow_pickle=False) as case:
                sensor, observed = solver.GS[case["used_rx"]], case["train"]
                assert np.array_equal(case["used_rx"], cases[name]["used_rx"])
                lo, hi = case["lows"], case["highs"]
            model = FullWave(solver.GD, sensor, solver.incident, observed)
            result = {}
            for label, v in native["refined_box_results"].items():
                if v["residual"] > row["radius"]:
                    continue
                material = arrays[f"{name}__box{int(label):02d}__material"]
                if not np.all((lo[int(label)]<=material)&(material<=hi[int(label)])):
                    raise ArithmeticError("The point witness left the box.")
                point_time = clock()
                _, j, _ = model.evaluate(material)
                try:
                    r = verifier.point(material, j, observed)
                    r["strictly_feasible"] = r["strict_residual_upper"] <= row["radius"]
                except (ArithmeticError, ZeroDivisionError, ValueError) as exc:
                    r = {"strictly_feasible": False, "reason": str(exc)}
                r["total_full_wave_and_audit"] = elapsed(point_time)
                result[label] = r
                print(f"{name} box {label}: upper/r={r.get('strict_residual_upper', np.inf)/row['radius']:.6f}", flush=True)
            receipts[name] = dict(boxes=result, full_wave_refactor_counts=dict(model.counts))
            write_json(OUTPUT/"strict_point_witnesses.json", dict(cases=receipts, geometry_setup=verifier.setup,
                       scope="Exact stored finite matrices and exact saved branch endpoints, both permitted signs."))
    return receipts


def neumann_runner(cases, verifier, solver):
    receipts = {}
    for name, row in cases.items():
        with np.load(OUTPUT/f"{name}.npz", allow_pickle=False) as case:
            lo, hi, y = case["lows"], case["highs"], case["train"]
        result = {}
        for label in range(len(lo)):
            start = clock()
            g, b = solver.GD, solver.incident
            m, d = .5*(lo[label]+hi[label]), .5*(hi[label]-lo[label])
            p = -m[:,None]*g
            p.flat[::len(m)+1] += 1
            factor = lu_factor(p, check_finite=False)
            v = np.random.default_rng(75011).standard_normal(len(m)).astype(complex)
            for _ in range(16):
                v = d*lu_solve(factor, g.conj().T@(g@lu_solve(factor,d*v,check_finite=False)),
                               trans=2,check_finite=False)
                v /= np.linalg.norm(v)
            state = lu_solve(factor, np.column_stack((m[:,None]*b, d*v)), check_finite=False)
            q = g@state[:,-1]
            try:
                r = verifier.contraction(lo[label], hi[label], v, state, y)
            except (ArithmeticError, ZeroDivisionError, ValueError) as exc:
                r = dict(reason=str(exc), strict_standard_lower=0.)
            r["numeric_kappa_lower_proxy"] = float(np.linalg.norm(q))
            r["count"] = dict(center_lu=1, power_witness_rhs=32, center_state_rhs=3)
            r["total"] = elapsed(start)
            result[str(label)] = r
            print(f"{name} box {label}: kappa lower={r.get('strict_kappa_lower', np.nan):.5f}, "
                  f"upper={r.get('strict_kappa_upper', np.nan):.5f}; "
                  f"Neumann/r={r['strict_standard_lower']/row['radius']:.5f}",flush=True)
            receipts[name] = dict(radius=row["radius"], branches=result)
            write_json(OUTPUT/"native_neumann.json",dict(cases=receipts,geometry_setup=verifier.setup))
    return receipts


def main(mode):
    ctx.prec, ctx.threads = 96, 1
    begin = clock()
    sys.path.insert(0,str(FLOW/"src"))
    import physics

    solver = physics.Solver()
    cases = json.loads((OUTPUT/"online.json").read_text())["cases"]
    for name, row in cases.items():
        with np.load(OUTPUT/f"{name}.npz", allow_pickle=False) as saved:
            row["used_rx"] = saved["used_rx"]
    receivers = next(iter(cases.values()))["used_rx"]
    if not all(np.array_equal(row["used_rx"],receivers) for row in cases.values()):
        raise ValueError("Different archived receiver masks.")
    sensors = solver.GS[receivers]
    if not (OUTPUT/"native_geometry.npz").exists():
        write_arrays(OUTPUT/"native_geometry.npz", green=solver.GD, sensor=solver.GS,
                     incident=solver.incident, points=solver.points, used_rx=receivers)
    write_json(OUTPUT/"support_only.json", dict(port_range=exact_surjectivity(sensors, solver.incident),
                scope="Full-grid signed branches, unrestricted supported currents and symmetric transfer."))
    verifier = SignedNative(solver.GD, sensors, solver.incident)
    print("exact native setup",verifier.setup,flush=True)
    if mode in ("point","all"):
        point_runner(cases, verifier, solver)
    if mode in ("neumann","all"):
        neumann_runner(cases, verifier, solver)
    print(elapsed(begin),flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase",choices=("point","neumann","all"))
    with threadpool_limits(limits=1):
        main(parser.parse_args().phase)
