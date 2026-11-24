"""Outward analytic-TM forward audits and receiver dual checks."""
import argparse
from pathlib import Path
import json
import sys

from flint import acb, acb_mat, arb, ctx
import numpy as np
from scipy.linalg import cholesky, eigh, solve_triangular, svd
from threadpoolctl import threadpool_limits

from run import (HERE, STUDY, adj, clock, elapsed, finite_geometry,
                 native_geometry, relax_receiver, square_root, write_json)

sys.path[:0] = [str(STUDY / "submission/curve"), str(STUDY)]
from core.certificate import ball, frobenius_squared, gershgorin, safe_float, star
from native import Geometry as NativeExactFFT


def num(mat):
    return np.asarray([complex(x) for x in mat.entries()]).reshape(mat.nrows(), mat.ncols())


def norm(mat):
    return frobenius_squared(mat).sqrt()


def hi(value):
    return safe_float(value, upper=True)


def lo(value):
    return safe_float(value, upper=False)


def eye(n):
    a = acb_mat(n,n)
    for i in range(n):
        a[i,i] = 1
    return a


def minor(a, rows, cols):
    return acb_mat(len(rows), len(cols), [a[int(i),int(j)] for i in rows for j in cols])


def join(a, b):
    return acb_mat(a.nrows(), a.ncols() + b.ncols(),
                   [v[i,j] for i in range(a.nrows()) for v in (a,b) for j in range(v.ncols())])


def trace_real(a):
    return sum((a[i,i].real for i in range(a.nrows())), arb(0))


def exact_number(x):
    return arb(str(x))


def unit_circle(n):
    return [(arg.cos(), arg.sin()) for arg in (2 * arb.pi() * i / n for i in range(n))]


def h1(order, x):
    return acb(x.bessel_j(order), x.bessel_y(order))


class FiniteExact:
    """Analytic disk collocation with exact coordinates, ports and incidence."""
    def __init__(self, geo, width):
        width = exact_number(width)
        side, n = geo.side, geo.n
        k = 8 * arb.pi() / 3
        q = k * width / side / arb.pi().sqrt()
        coeff = arb.pi() * q * q.bessel_j(1) / 2
        kernel = {}
        for dy in range(side):
            for dx in range(side):
                if dx + dy:
                    kernel[dx,dy] = acb(0, coeff) * h1(0, k * width / side * arb(dx*dx + dy*dy).sqrt())
        kernel[0,0] = acb(0, arb.pi() * q / 2) * h1(1, q) - 1
        self.g = acb_mat(n,n,[kernel[abs(i%side-j%side),abs(i//side-j//side)] for i in range(n) for j in range(n)])
        points = [(width*(2*(i%side) + 1 - side)/(2*side),
                   width*(2*(i//side) + 1 - side)/(2*side)) for i in range(n)]
        illuminations = unit_circle(geo.inputs)
        receivers = [(9*u/5,9*v/5) for u,v in unit_circle(geo.p)]
        self.b = acb_mat(n,geo.inputs,[(arb(10) / geo.inputs).sqrt() * acb(0,k*(x*u + y*v)).exp()
                                      for x,y in points for u,v in illuminations])
        self.s = acb_mat(geo.p,n,[(arb(20) / geo.p).sqrt() * acb(0,coeff) *
                                 h1(0,k*((x-u)**2 + (y-v)**2).sqrt())
                                 for u,v in receivers for x,y in points])
        self.r = (self.g - star(self.g)) / acb(0,2)
        self.r_negative_bound = arb(0)

    def green(self, a):
        return self.g * a

    def radiation(self, a):
        return self.r * a


class NativeExact:
    """Exact represented v3 coefficients, with an analytic radiation guard."""
    def __init__(self, geo):
        self.fft = NativeExactFFT(geo.spectrum, geo.s, geo.b)
        self.b, self.s = self.fft.incident_ball, self.fft.sensor_ball
        full = self.fft.spec_ball
        self.r_spectrum = [[acb(z.imag) for z in row] for row in full]
        self.r_negative_bound = self.radiation_guard()

    def radiation_guard(self):
        # Exact inverse DFT of the represented Hermitian multiplier. Compare
        # every offset that occurs in the restricted native matrix to its
        # analytic disk kernel. Entrywise discrepancy gives an n*eps bound.
        p, n = self.fft.pad, self.fft.side
        horizontal = [acb.dft(row,inverse=True) for row in self.r_spectrum]
        vertical = [acb.dft([row[j] for row in horizontal],inverse=True) for j in range(p)]
        k = 8*arb.pi()/3
        q = k*2/n/arb.pi().sqrt()
        coeff = arb.pi()*q*q.bessel_j(1)/2
        analytic = {d2: coeff*(k*2/n*arb(d2).sqrt()).bessel_j(0)
                    for d2 in {i*i+j*j for i in range(n) for j in range(n)}}
        differences = [abs(vertical[j%p][i%p]-analytic[i*i+j*j]).upper()
                       for i in range(1-n,n) for j in range(1-n,n)]
        return (n*n*max(differences)).upper()

    def green(self, a):
        return self.fft.exact_green(a)

    def radiation(self, a):
        save = self.fft.spec_ball
        self.fft.spec_ball = self.r_spectrum
        try:
            return self.fft.exact_green(a)
        finally:
            self.fft.spec_ball = save


def verify_observation(geo, exact, x, tangent, y, rho):
    start = clock()
    nominal, active, current, its = geo.forward(x, float(tangent))
    x_cap = max(float(v) for v in x)
    # A=D(1/chi)-G on the active set, with -Im(A) >= alpha - guard.
    reciprocal = [acb(1,-tangent) / ((1 + tangent*tangent) * arb(float(x[i]))) for i in active]
    alpha_lower = tangent / ((1 + tangent*tangent) * arb(x_cap)) - exact.r_negative_bound
    if not alpha_lower > 0:
        raise ArithmeticError("Forward coercivity was not certified.")
    s = minor(exact.s, range(geo.p), active)
    b = minor(exact.b, active, range(geo.inputs))
    current_b = ball(current)
    if isinstance(exact, FiniteExact):
        action = minor(exact.g, active, active) * current_b
    else:
        full = np.zeros_like(geo.b)
        full[active] = current
        action = minor(exact.green(ball(full)), active, range(geo.inputs))
    error = acb_mat(len(active),geo.inputs,[b[i,j] - reciprocal[i] * current_b[i,j] + action[i,j]
                   for i in range(len(active)) for j in range(geo.inputs)])
    outward_error = norm(s) * norm(error) / alpha_lower
    forward = s * current_b
    upper = norm(ball(y) - forward) + outward_error
    margin = arb(rho) - upper
    if not margin > 0:
        raise ArithmeticError("Generating material's exact response was not enclosed in the data ball.")
    return dict(rho=rho, radius_margin_lower=lo(margin), forward_residual_upper=hi(upper),
                nominal_clean_error_upper=hi(norm(ball(nominal)-forward) + outward_error),
                coercivity_lower=lo(alpha_lower), active_cells=len(active),
                state_defect_upper=hi(norm(error)), **elapsed(start))


def gram_interval(geo, exact, cap, tangent):
    start = clock()
    alpha = tangent / ((1 + tangent*tangent) * cap)
    q_b = join(exact.b, star(exact.s))
    gap = alpha - exact.r_negative_bound
    if not gap > 0:
        raise ArithmeticError("Dissipative inverse was not certified.")
    alpha_num = float(alpha)
    proposal, its = geo.radiation_inverse(np.column_stack((geo.b, adj(geo.s))), alpha_num)
    x_b = ball(proposal)
    ax = exact.radiation(x_b) + x_b * alpha
    e = q_b - ax
    # Q*M^-1 Q = trial_form + E*M^-1 E.
    trial_form = star(x_b) * q_b + star(q_b) * x_b - star(x_b) * ax
    defect_squared = frobenius_squared(e)
    error = (defect_squared / gap).upper()
    return trial_form, error, proposal, dict(alpha_lower=lo(alpha), loss_pd_gap_lower=lo(gap),
              gram_shift_upper=hi(error), defect_frob_upper=hi(defect_squared.sqrt()),
              inverse_iterations=its, **elapsed(start))


def blocks(geo, trial, error):
    s, p = geo.inputs, geo.p
    left, right = range(s), range(s,s+p)
    cb = minor(trial, left,left)/4
    rs = minor(trial, right,right)
    fc = minor(trial, right,left)*acb(0,arb(1)/2)
    return rs+eye(p)*error, cb+eye(s)*(error/4), fc


def dual_evaluate(geo, y, z, trial, error, *, rank_cutoff=1.e-9):
    start = clock()
    rs, cb, fc = blocks(geo, trial, error)
    lr, lc = square_root(num(rs)), square_root(num(cb))
    r_eig = eigh(lr)
    c_eig = eigh(lc)
    # Small eigenvalues can be truncated. The exact trial Z is the product
    # of the two dyadic matrices subsequently verified.
    rhalf = r_eig[1][:,r_eig[0] > rank_cutoff * r_eig[0][-1]]
    chalf = c_eig[1][:,c_eig[0] > rank_cutoff * c_eig[0][-1]]
    ri = (rhalf / r_eig[0][-rhalf.shape[1]:]) @ adj(rhalf)
    ci = (chalf / c_eig[0][-chalf.shape[1]:]) @ adj(chalf)
    a, singular, bh = svd(lr @ z @ lc, full_matrices=False)
    u = ball((ri @ a) * np.sqrt(singular))
    v = ball((ci @ adj(bh)) * np.sqrt(singular))
    zb = u * star(v)
    nrm = norm(zb)
    if not nrm > 0:
        return dict(dual_lower=0., measured_norm_upper=hi(nrm), **elapsed(start))
    # Nuclear norm <= ||R^(1/2) U||_F ||C^(1/2) V||_F.
    left = trace_real(star(u) * rs * u)
    right = trace_real(star(v) * cb * v)
    support_upper = (left.upper() * right.upper()).sqrt()
    center_error = (min(geo.inputs,geo.p)*arb(1)).sqrt() * error / 2
    normalized_bound = ((trace_real(star(zb) * (ball(y) - fc)) - support_upper -
                        nrm * center_error) / nrm)
    return dict(dual_lower=max(0.,lo(normalized_bound)), measured_norm_upper=hi(nrm),
                support_upper=hi(support_upper), center_allowance_upper=hi(center_error), **elapsed(start))


def primal_evaluate(geo, y, trial, error):
    start = clock()
    rs, cb, fc = blocks(geo, trial, arb(0))
    _, _, point = relax_receiver(y, num(rs), num(cb), num(fc), return_point=True)
    center = num(fc)
    for shrink in (1.e-5, 1.e-4, 1.e-3, .01):
        proposal = center + (1.-shrink)*(point-center)
        point_b = ball(proposal)
        offset = point_b - fc
        block = acb_mat(geo.inputs+geo.p, geo.inputs+geo.p, [
            item for slab in (join(cb, star(offset)), join(offset, rs)) for item in slab.entries()
        ])
        try:
            chol = cholesky(num(block), lower=True)
        except np.linalg.LinAlgError:
            continue
        inv = ball(solve_triangular(chol, np.eye(len(chol)), lower=True))
        gap = gershgorin(inv * block * star(inv), upper=False)
        if gap > 0:
            # The exact block differs from the trial block by a PSD
            # congruence of E*M^-1*E. The same point is therefore feasible.
            return dict(primal_upper=hi(norm(ball(y)-point_b)), block_congruence_gap_lower=lo(gap),
                        contraction_shrink=shrink, **elapsed(start))
    return dict(primal_status="positive_block_unverified", **elapsed(start))


def scalar_primal(geo, exact, y, trial, error, proposal, alpha):
    start = clock()
    rs, cb, fc = blocks(geo, trial, error)
    # Certify a scalar cap at least this large by lower bounds on the two
    # terms of the same triangle inequality used by the scalar baseline.
    # Every upper cap chosen there contains this smaller ball too.
    cb_mid = num(cb)
    _, cb_vec = eigh(cb_mid)
    direction = ball(cb_vec[:,[-1]])
    cb_lower = trace_real(star(direction) * (cb - eye(geo.inputs)*(error/4)) * direction) / frobenius_squared(direction)
    j_num = proposal[:, :geo.inputs] / 2
    _, _, vh = svd(j_num, full_matrices=False)
    j_direction = ball(adj(vh[:, :])[:,[0]])
    # A convenient lower bound needs only M status and the certified Gram
    # defect upper, so no duplicate M applications are paid here.
    gap = alpha - exact.r_negative_bound
    inverse_error = (error/gap).sqrt()
    j_lower = norm(ball(j_num)*j_direction) / norm(j_direction) - inverse_error/2
    current_cap = (max(arb(0),j_lower.lower()) + (max(arb(0),cb_lower.lower())/alpha).sqrt())**2
    if current_cap <= 0:
        return dict(scalar_status="nonpositive_coherent_cap", **elapsed(start))
    ss = exact.s * star(exact.s)
    _, _, point = relax_receiver(y, num(ss), float(current_cap)*np.eye(geo.inputs), np.zeros_like(y),
                                  maxiter=1200, return_point=True)
    for shrink in (1.e-5, 1.e-4, 1.e-3, .01):
        point_b = ball((1-shrink)*point)
        source = eye(geo.inputs)*current_cap.lower()
        block = acb_mat(geo.inputs+geo.p,geo.inputs+geo.p, [
            item for slab in (join(source,star(point_b)),join(point_b,ss)) for item in slab.entries()
        ])
        try:
            chol = cholesky(num(block), lower=True)
        except np.linalg.LinAlgError:
            continue
        inv = ball(solve_triangular(chol, np.eye(len(chol)), lower=True))
        gap = gershgorin(inv * block * star(inv), upper=False)
        if gap > 0:
            return dict(scalar_primal_upper=hi(norm(ball(y)-point_b)), cap_lower=lo(current_cap),
                        block_congruence_gap_lower=lo(gap), **elapsed(start))
    return dict(scalar_status="positive_block_unverified", **elapsed(start))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=["finite_0", "finite_1", "native"])
    parser.add_argument("--forward", default="all", choices=["all","low","none"])
    parser.add_argument("--audit-input", type=Path)
    args = parser.parse_args()
    start = clock()
    data = json.loads((HERE / f"{args.study}.json").read_text())
    native = args.study == "native"
    width = {"finite_0":.6,"finite_1":1.25}.get(args.study)
    geo = native_geometry() if native else finite_geometry(width)
    exact = NativeExact(geo) if native else FiniteExact(geo,width)
    out = dict(study=args.study, precision_bits=ctx.prec,
               model=("exact v3 binary64 spectrum/S/B and analytic radiation comparison"
                      if native else "analytic disk-collocation TM with exact dimensions, phases and Hankel coefficients"),
               guard_on_negative_radiation=hi(exact.r_negative_bound), inputs=[], projected=[])
    target = HERE / f"{args.study}_verified.json"
    if args.audit_input:
        old = json.loads(args.audit_input.read_text())
        out["inputs"] = old["inputs"]
        out["projected"] = old["projected"]
    already = {o["id"] for o in out["inputs"]}
    complete = {(o["tangent"], o["cap"]) for o in out["projected"]}
    with np.load(HERE / f"{args.study}_inputs.npz",allow_pickle=False) as saved:
        for row in data["observations"]:
            ident, name = row["id"], row["object_id"]
            if ident in already or args.forward == "none" or (args.forward == "low" and not name.endswith("_low")):
                continue
            tangent = exact_number(row["tangent"])
            result = verify_observation(geo, exact, saved[name], tangent, saved[ident], row["rho"])
            print("forward", ident, result["radius_margin_lower"], flush=True)
            out["inputs"].append(dict(id=ident, **result))
            write_json(target, {**out, **elapsed(start)})
        for entry in data["receivers"]:
            t, u = entry["tangent"], entry["cap"]
            if (t,u) in complete:
                continue
            # Review every finite-grid numerical crossing. At native resolution
            # verify all proposed positives, plus the cap-1 point of each loss.
            comparisons = entry["rows"]
            selected = [row for row in comparisons if row["ellipsoid"]["dual"] > row["rho"] or
                        (native and u in (.75,1.)) or (not native and u in (.75,1.) and t >= .05)]
            if not selected:
                continue
            trial, error, proposal, setup = gram_interval(geo, exact, exact_number(u), exact_number(t))
            rows = []
            for row in selected:
                ident = row["id"]
                z = saved[f"{ident}__cap{u:g}"]
                result = dual_evaluate(geo, saved[ident], z, trial, error)
                primal = primal_evaluate(geo, saved[ident], trial, error)
                scalar = None
                if not native and result["dual_lower"] > row["rho"]:
                    alpha = exact_number(t) / ((1 + exact_number(t)**2) * exact_number(u))
                    scalar = scalar_primal(geo, exact, saved[ident], trial, error, proposal, alpha)
                print("dual", ident, u, result["dual_lower"]/row["rho"],
                      "primal", primal.get("primal_upper",0.)/row["rho"], flush=True)
                rows.append(dict(id=ident, rho=row["rho"], excludes=(result["dual_lower"]>row["rho"]),
                                 **result, primal=primal, scalar=scalar))
            out["projected"].append(dict(tangent=t, cap=u, setup=setup, rows=rows))
            write_json(target,{**out, **elapsed(start)})
    write_json(target,{**out, **elapsed(start)})


if __name__ == "__main__":
    ctx.prec = 128
    ctx.threads = 1
    with threadpool_limits(1):
        main()
