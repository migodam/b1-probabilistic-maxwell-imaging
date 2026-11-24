"""Finite-grid lossy TM forward data and matrix bounds."""
import json
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, solve, svd
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import LinearOperator, cg, gmres
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]
TANGENTS = (0.01, 0.05, 0.2)
CAPS = (0.75, 1., 1.25, 1.5, 2., 3.)
PERCENT = 0.02
SEED = 20261125014
MATERIAL_SEED = 20261125013


def elapsed(start):
    return dict(wall_seconds=perf_counter() - start[0], cpu_seconds=process_time() - start[1])


def clock():
    return perf_counter(), process_time()


def herm(x):
    return (x + x.conj().T) / 2


def adj(x):
    return x.conj().T


def square_root(mat):
    values, vectors = eigh(herm(mat))
    return (vectors * np.sqrt(np.maximum(values, 0.))) @ adj(vectors)


def descending(mat):
    return eigh(herm(mat), eigvals_only=True)[::-1]


def tail(mat, skip=0):
    return float(np.sqrt(np.maximum(descending(mat)[skip:], 0).sum()))


def cap_coefficients(cap, tangent):
    beta = 1 / ((1 + tangent**2) * cap)
    return beta, tangent * beta


def make_noise(clean, object_index):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, object_index]))
    direction = rng.normal(size=clean.shape) + 1j * rng.normal(size=clean.shape)
    noise = direction * (PERCENT * np.linalg.norm(clean) / np.linalg.norm(direction))
    y = clean + noise
    return y, float(np.linalg.norm(y - clean) + 1.e-6 * np.linalg.norm(y))


def write_json(path, value):
    tmp = path.with_suffix(".partial")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


class Geometry:
    def __init__(self, side, s, b, *, dense=None, spectrum=None):
        self.side, self.s, self.b = side, s, b
        self.n, self.p, self.inputs = len(b), len(s), b.shape[1]
        self.dense, self.spectrum = dense, spectrum
        self.dk = None if dense is None else (dense - adj(dense)) / (2j)
        self.hk = None if dense is None else herm(dense)
        if spectrum is not None:
            self.rad_spectrum = spectrum.imag.copy()
        self.columns = 0

    def green(self, matrix, *, radiation=False):
        if self.dense is not None:
            return (self.dk if radiation else self.dense) @ matrix
        vector = matrix.reshape(self.side, self.side, -1)
        spectrum = self.rad_spectrum if radiation else self.spectrum
        pad = len(spectrum)
        self.columns += vector.shape[-1]
        transformed = fft2(vector, s=(pad, pad), axes=(0, 1))
        values = ifft2(spectrum[:, :, None] * transformed, axes=(0, 1))[:self.side, :self.side]
        return values.reshape(matrix.shape)

    def radiation_inverse(self, values, alpha):
        if self.dense is not None:
            return solve(self.dk + alpha * np.eye(self.n), values, assume_a="pos"), None

        def apply(vector):
            return self.green(vector, radiation=True) + alpha * vector

        operator = LinearOperator((self.n, self.n), matvec=apply, dtype=np.complex128)
        output, iterations = [], []
        for rhs in values.T:
            trace = []
            x, info = cg(operator, rhs, rtol=3e-12, atol=0, maxiter=450,
                         callback=lambda _: trace.append(None))
            if info:
                raise ArithmeticError(f"Radiation solve {info}")
            output.append(x)
            iterations.append(len(trace))
        return np.column_stack(output), iterations

    def forward(self, x, tangent):
        active = np.flatnonzero(x)
        rate = ((1 + 1j * tangent) * x[active])
        reciprocal = 1 / rate
        if self.dense is not None:
            system = np.diag(reciprocal) - self.dense[np.ix_(active, active)]
            trial = solve(system, self.b[active], assume_a="gen")
            iterations = None
        else:
            def apply(vector):
                current = np.zeros(self.n, complex)
                current[active] = vector
                return reciprocal * vector - self.green(current)[active]

            system = LinearOperator((len(active), len(active)), matvec=apply, dtype=np.complex128)
            pre = LinearOperator(system.shape, matvec=lambda v: v / reciprocal, dtype=np.complex128)
            trial, iterations = [], []
            for rhs in self.b[active].T:
                trace = []
                sol, info = gmres(system, rhs, M=pre, rtol=2e-12, atol=0, restart=70, maxiter=8,
                                  callback=lambda _: trace.append(None), callback_type="pr_norm")
                if info:
                    raise ArithmeticError(f"Full wave solve {info}")
                trial.append(sol)
                iterations.append(len(trace))
            trial = np.column_stack(trial)
        return self.s[:, active] @ trial, active, trial, iterations


def positions(side, width):
    line = (np.arange(side) + 0.5) * width / side - width / 2
    return np.column_stack([v.ravel() for v in np.meshgrid(line, line)])


def directions(n):
    angle = 2 * np.pi * np.arange(n) / n
    return np.column_stack((np.cos(angle), np.sin(angle)))


def finite_geometry(width, side=12, p=20, sources=10):
    points = positions(side, width)
    wavenumber = 8 * np.pi / 3
    q = wavenumber * width / side / np.sqrt(np.pi)
    factor = np.pi * q * jv(1, q) / 2
    r = np.linalg.norm(points[:, None] - points[None], axis=2)
    np.fill_diagonal(r, 1)
    dense = 1j * factor * hankel1(0, wavenumber * r)
    np.fill_diagonal(dense, 1j*np.pi*q*hankel1(1, q)/2 - 1)
    rx = 1.8 * directions(p)
    rs = np.linalg.norm(rx[:, None] - points[None], axis=2)
    s = np.sqrt(20 / p) * 1j * factor * hankel1(0, wavenumber * rs)
    b = np.sqrt(10 / sources) * np.exp(1j * wavenumber * (points @ directions(sources).T))
    return Geometry(side, s, b, dense=dense)


def finite_material(width, index, side=12):
    rng = np.random.default_rng(np.random.SeedSequence([MATERIAL_SEED, index]))
    points = positions(side, width)
    shift = rng.uniform(-0.06, 0.06, size=2) * width
    angle = rng.uniform(-0.8, 0.8)
    rot = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
    x, y = ((points - shift) @ rot / width).T
    a, b = rng.uniform([.20, .16], [.27, .21])
    strength = rng.integers(72, 91) / 25.
    ellipse = ((x - 0.04 * (index % 2)) / a)**2 + (y / b)**2 <= 1
    coef = ellipse * np.where(y > 0, 1., .82)
    if index % 2:
        disk = ((x + .27)**2 + (y - .10)**2) <= .13**2
        coef = np.maximum(coef, .70 * disk)
    return strength * coef


def native_geometry():
    with np.load(STUDY / "synthetic/public/fine_geometry.npz") as f:
        spec = f["spectrum"]
    with np.load(STUDY / "physical_revision/array/observations/arrays.npz") as f:
        s, b = f["S32"], f["B32"]
    return Geometry(96, s, b, spectrum=spec)


def receiver_setup(geo, alpha):
    started = clock()
    q = np.column_stack((geo.b, adj(geo.s)))
    x, its = geo.radiation_inverse(q, alpha)
    c = herm(adj(q) @ x)
    s = geo.inputs
    rs, cb, fc = c[s:, s:], c[:s, :s] / 4, .5j * c[s:, :s]
    # Uniform coherent-current cap from the same, richer ellipsoid.
    jcenter = x[:, :s] / 2
    scalar_cap = (svd(jcenter, compute_uv=False)[0] + np.sqrt(descending(cb)[0] / alpha))**2
    return (rs, cb, fc, scalar_cap), {"inverse_iterations": its, **elapsed(started)}, x


def relax_receiver(y, rs, cb, fc, *, maxiter=2500, return_point=False):
    """FISTA on the primal contraction, with a dual candidate at every step."""
    start = clock()
    lr, lc = square_root(rs), square_root(cb)
    residual_c = y - fc
    scale = np.linalg.norm(residual_c)
    lipschitz = descending(rs)[0] * descending(cb)[0]
    v = np.zeros_like(y)
    extrap = v.copy()
    a = 1.
    best_dual, best_primal, z_best = 0., scale, v
    best_point = fc
    stop = "iteration_limit"
    for it in range(1, maxiter + 1):
        gradient = lr @ (lr @ extrap @ lc - residual_c) @ lc
        u, sig, vh = svd(extrap - gradient / lipschitz, full_matrices=False)
        now = (u * np.minimum(sig, 1.)) @ vh
        residual = residual_c - lr @ now @ lc
        norm = np.linalg.norm(residual)
        if norm < best_primal:
            best_primal = float(norm)
            best_point = y - residual
        if norm > 0:
            z = residual / norm
            dual = np.vdot(z, residual_c).real - svd(lr @ z @ lc, compute_uv=False).sum()
            if dual > best_dual:
                best_dual, z_best = float(dual), z.copy()
        if best_primal <= 2.e-6 * scale:
            stop = "small_primal_upper"
            break
        if best_primal - best_dual <= 2.e-6 * scale:
            stop = "primal_dual_gap"
            break
        new_a = (1 + np.sqrt(1 + 4 * a*a)) / 2
        extrap = now + ((a - 1) / new_a) * (now - v)
        a, v = new_a, now
    stats = dict(dual=best_dual, primal=best_primal, gap=best_primal - best_dual,
                 iterations=it, stop=stop, **elapsed(start))
    return (stats, z_best, best_point) if return_point else (stats, z_best)


def inspect_profile(geo, y, beta, alpha, mu, optical, mode, cb, x_recip, best):
    """A legal weighted current quadratic. Dividing by mu improves scaling."""
    if mode == "old":
        # beta ||J||^2 - Re(J*E) <= 0 and Im(J*E) <= 0.
        a_re, a_im, diag = 1., optical, beta
    elif mode == "fixed_phase":
        # Add optical times [Im(J*E) + t Re(J*E) = 0].
        a_re, a_im, diag = 1. - optical * alpha / beta, optical, beta
    elif mode == "shift":
        # Add optical times [J*(R+alpha I)J - Im(B*J) <= 0].
        a_re, a_im, diag = 1., optical, beta + optical * alpha
    else:
        raise ValueError(mode)
    k = adj(geo.s) @ geo.s / mu - a_re * geo.hk + a_im * geo.dk + diag * np.eye(geo.n)
    z = (a_re + 1j * a_im) / 2
    l = adj(geo.s) @ y / mu + z * geo.b
    kval, kv = eigh(herm(k))
    if min(abs(kval)) < 1e-11:
        return
    ji = kv @ ((adj(kv) @ l) / kval[:, None])
    hy = herm(adj(y) @ y - mu * adj(l) @ ji)
    nu = int((kval < 0).sum())
    d = tail(hy, nu)
    candidate = dict(distance=d, mu=mu, weight=optical, nu=nu,
                     min_absolute_eigenvalue=float(np.min(abs(kval))))
    if best[mode]["distance"] < d or "nu" not in best[mode]:
        best[mode] = candidate

    if mode != "old" or not nu:
        return
    # Project all negative directions, then split off the most expensive ones
    # in the induced ellipsoid metric. Preserve their rank contribution.
    w = np.sqrt(-mu * kval[:nu, None]) * adj(kv[:, :nu])
    center = w @ (.5j * x_recip[:, :geo.inputs] - ji)
    inverse = solve(geo.dk + alpha * np.eye(geo.n), adj(w), assume_a="pos")
    vals, q = eigh(herm(w @ inverse))
    center = adj(q) @ center
    for free in range(min(nu, geo.inputs - 1) + 1):
        if free == nu:
            continue
        keep = nu - free
        grams = center[:keep]
        eigen = vals[:keep]
        oracle = float(eigen[-1])

        def objective(parameter):
            t = oracle * (1 + np.exp(parameter))
            matcap = t * cb + t * adj(grams) @ (grams / (t - eigen[:, None]))
            return tail(hy - matcap, free)

        # Tail under a matrix cap is not generally unimodal.
        parameters = np.linspace(-7, 6, 14)
        candidates = [objective(v) for v in parameters]
        index = int(np.argmax(candidates))
        d = candidates[index]
        parameter = float(parameters[index])
        if 0 < index < len(parameters) - 1 and d > 0:
            opt = minimize_scalar(lambda v: -objective(v), bounds=(parameters[index-1], parameters[index+1]),
                                  method="bounded", options={"maxiter": 22})
            if -opt.fun > d:
                d = -float(opt.fun)
                parameter = float(opt.x)
        label = "full_projected_negative" if free == 0 else "partial_projected_negative"
        if best[label]["distance"] < d or "nu" not in best[label]:
            best[label] = dict(distance=d, mu=mu, weight=optical, nu=nu, free=free,
                              projected_multiplier=oracle * (1 + np.exp(parameter)),
                              gamma_max=oracle)


def finite_baselines(geo, y, beta, alpha, cb, x_recip):
    start = clock()
    labels = ("old", "fixed_phase", "shift", "full_projected_negative", "partial_projected_negative")
    best = {mode: {"distance": 0.} for mode in labels}
    scale = descending(geo.s @ adj(geo.s))[0] / max(beta, .2)
    for mu in scale * np.logspace(-3, .5, 8):
        for weight in (0., .1, .5, 2., 10., 40.):
            inspect_profile(geo, y, beta, alpha, mu, weight, "old", cb, x_recip, best)
            inspect_profile(geo, y, beta, alpha, mu, weight, "shift", cb, x_recip, best)
        for weight in (-2., -.5, 0., .2, .8, 1., 1.2, 2.):
            inspect_profile(geo, y, beta, alpha, mu, weight * beta / alpha, "fixed_phase", cb, x_recip, best)
    return {"best": best, "profiles_by_arm": [48, 64, 48], **elapsed(start)}


def run_block(geo, names, materials, caps, outname, *, compare):
    totals = clock()
    receivers, observations, measurements, contrasts, zc = [], [], {}, {}, {}
    for index, (name, x) in enumerate(zip(names, materials, strict=True)):
        contrasts[name] = x
        for tangent in TANGENTS:
            start = clock()
            clean, active, current, its = geo.forward(x, tangent)
            y, rho = make_noise(clean, index)
            ident = f"{name}__fixed_t{tangent:g}__pair2pct"
            measurements[ident] = y
            observations.append(dict(id=ident, object_id=name, tangent=tangent,
                                     rho=rho, peak=float(np.max(x)), active=len(active),
                                     forward_iterations=its, **elapsed(start)))
            print(ident, rho, flush=True)
    for tangent in TANGENTS:
        for u in caps:
            bias = geo.columns
            beta, alpha = cap_coefficients(u, tangent)
            (rs, cb, fc, scalar_cap), time, x = receiver_setup(geo, alpha)
            rows = []
            for obs in observations:
                if obs["tangent"] != tangent:
                    continue
                ident = obs["id"]
                y = measurements[ident]
                base, z = relax_receiver(y, rs, cb, fc)
                # A scalar ellipsoid in current space, with an operator-norm
                # constraint and zero center. The receiving shape is retained.
                energy, _ = relax_receiver(y, geo.s @ adj(geo.s), scalar_cap * np.eye(geo.inputs),
                                           np.zeros_like(y), maxiter=800)
                record = dict(id=ident, rho=obs["rho"], ellipsoid=base, scalar_energy=energy)
                if compare:
                    record["state_screen"] = finite_baselines(geo, y, beta, alpha, cb, x)
                zc[f"{ident}__cap{u:g}"] = z
                rows.append(record)
            receivers.append(dict(tangent=tangent, cap=u, alpha=alpha, beta=beta,
                scalar_coherent_energy_cap=float(scalar_cap), rs_singular_extrema=descending(rs)[[0,-1]].tolist(),
                cb_singular_extrema=descending(cb)[[0,-1]].tolist(), setup=time, rows=rows,
                green_application_columns=geo.columns-bias))
            print(outname, tangent, u, [(a["id"].split("__")[0],
                round(a["ellipsoid"]["dual"] / a["rho"], 5),
                round(a["ellipsoid"]["primal"] / a["rho"], 5)) for a in rows], flush=True)
            write_json(HERE / f"{outname}.json", dict(observations=observations, receivers=receivers,
                       green_application_columns=geo.columns, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /
                       (2**20 if sys.platform == "darwin" else 1024), **elapsed(totals)))
            np.savez_compressed(HERE / f"{outname}_inputs.npz", **measurements, **contrasts, **zc)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "finite"
    if mode == "finite":
        for index, width in enumerate((.6, 1.25)):
            names = [f"loss_tm_{2*index+i+1:02}" for i in range(2)]
            x = [finite_material(width, 2*index+i) for i in range(2)]
            run_block(finite_geometry(width), names + [names[0] + "_low"], x + [x[0]/8], CAPS,
                      f"finite_{index}", compare=True)
    elif mode == "native":
        names = [f"optfresh_{i:02}" for i in range(1, 7)]
        fields = []
        for name in names:
            with np.load(STUDY / f"synthetic/evaluation/{name}.npz") as f:
                fields.append(f["m_fine"])
        run_block(native_geometry(), names + [names[0] + "_low"], fields + [fields[0]/8], CAPS,
                  "native", compare=False)
    elif mode == "finite_quick":
        for index, width in enumerate((.6, 1.25)):
            names = [f"loss_tm_{2*index+i+1:02}" for i in range(2)]
            x = [finite_material(width, 2*index+i) for i in range(2)]
            run_block(finite_geometry(width), names + [names[0] + "_low"], x + [x[0]/8], CAPS,
                      f"finite_{index}", compare=False)
    else:
        raise ValueError(mode)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
