"""Replay two frozen coarse optical block profiles with the native fine FFT."""
import argparse
import json
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import LinearOperator, eigsh, gmres
from threadpoolctl import threadpool_limits

from fourier_bound import DATA, PUBLIC, STUDY, h, read

CASES = ("optfresh_03__pw180_10x20__n10", "optfresh_04__pw360_10x20__n02")
SIDE = 96
N = SIDE ** 2
FACTOR = (SIDE // 32) ** 2


class NativeProfile:
    """Hermitian K/mu, formed only through native-grid FFT actions."""

    def __init__(self, spec, sensor, incident, data, cap, mu, lam):
        self.spec = spec
        self.pad = spec.shape
        self.sensor = sensor
        self.sensor_star = sensor.conj().T
        self.cap = cap
        self.mu = mu
        z = .5 * (cap + 1j * np.broadcast_to(lam, (N,)) / mu)
        self.z = z[:, None]
        self.rhs = self.sensor_star @ data / mu + self.z * incident
        self.observation_gram = h(data.conj().T @ data)
        self.cost = dict(operator_applications=0, complex_columns=0, green_fft_convolutions=0)
        self.operator = LinearOperator((N, N), matvec=self.action, matmat=self.action, dtype=np.complex128)

    def green(self, x, *, adjoint=False):
        ns = x.shape[1]
        fk = fft2(x.reshape(SIDE, SIDE, ns), s=self.pad, axes=(0, 1))
        kernel = self.spec.conj() if adjoint else self.spec
        result = ifft2(kernel[:, :, None] * fk, axes=(0, 1))[:SIDE, :SIDE]
        return result.reshape(N, ns)

    def action(self, x):
        original_shape = x.shape
        x = x.reshape(N, -1)
        self.cost["operator_applications"] += 1
        self.cost["complex_columns"] += x.shape[1]
        self.cost["green_fft_convolutions"] += 2*x.shape[1]
        answer = (x + self.sensor_star @ (self.sensor @ x) / self.mu
                  - self.z * self.green(x) - self.green(self.z.conj() * x, adjoint=True))
        return answer.reshape(original_shape)

    def diagnostic(self, radius):
        start, cpu = perf_counter(), process_time()
        ns = self.rhs.shape[1]
        rng = np.random.default_rng(54273)
        trial = rng.standard_normal(N) + 1j*rng.standard_normal(N)
        eigs, modes = eigsh(self.operator, k=ns+3, which="SA", ncv=45, maxiter=500,
                             tol=2.e-11, v0=trial)
        order = np.argsort(eigs)
        eigs, modes = eigs[order], modes[:, order]
        residuals = np.linalg.norm(self.action(modes)-modes*eigs, axis=0)
        nu = int(np.count_nonzero(eigs < 0.))
        result = dict(native_side=SIDE, dimension=N, mu=self.mu, cap=self.cap,
                      lowest_scaled_K_eigenvalues=eigs.tolist(),
                      lowest_scaled_K_residuals=residuals.tolist(), numeric_inertia=nu,
                      radius=radius)
        if nu >= ns:
            result.update(status="numerical_illumination_budget_exhausted", distance=0., margin=-radius)
        elif np.min(abs(eigs)) < 10*np.max(residuals):
            result.update(status="numerical_inertia_unresolved", distance=None, margin=None)
        else:
            solutions = np.empty_like(self.rhs)
            solves = []
            for j, column in enumerate(self.rhs.T):
                trace = []
                solutions[:, j], info = gmres(self.operator, column, restart=100, maxiter=10,
                                              rtol=5.e-11, atol=0., callback=trace.append,
                                              callback_type="pr_norm")
                solves.append(dict(iterations=len(trace), solver_info=int(info)))
            backward = np.linalg.norm(self.action(solutions)-self.rhs)/np.linalg.norm(self.rhs)
            result.update(solves=solves, relative_stationary_residual=float(backward))
            if backward > 1.e-8 or any(r["solver_info"] for r in solves):
                result.update(status="numerical_stationarity_unresolved", distance=None, margin=None)
            else:
                form = h(self.observation_gram - self.mu * self.rhs.conj().T @ solutions)
                spectrum = eigvalsh(form, check_finite=False)[::-1]
                distance = float(np.sqrt(np.maximum(spectrum[nu:], 0.).sum()))
                self.solution, self.source_form = solutions, form
                result.update(status="numerical_crossing" if distance > radius else "numerical_below_radius",
                              distance=distance, margin=distance-radius, H_eigenvalues=spectrum.tolist())
        return dict(**result, operations=self.cost, wall_seconds=perf_counter()-start,
                    cpu_seconds=process_time()-cpu)


def repeated_blocks(proposal):
    blocks = np.asarray(proposal["block_lambda"]).reshape(4, 4)
    return np.repeat(np.repeat(blocks, SIDE//4, axis=0), SIDE//4, axis=1).ravel()


def main():
    observations = {m["id"]: m for m in json.loads((PUBLIC/"manifest.json").read_text())["observations"]}
    coarse = json.loads((STUDY/"results/core/gradient_certificates.json").read_text())["cases"]
    rows = {}
    for name in CASES:
        meta = observations[name]
        spec, sensor, incident, data = read(name, meta["geometry"])
        record = coarse[name]
        radius = meta["radius_noise_fro"]
        mapped = repeated_blocks(record["block_proposal"])
        case_record = dict(radius=radius, coarse_cap=record["cap"],
                           coarse_block_proposal=record["block_proposal"],
                           coarse_scalar_proposal=record["scalar_proposal"], native={})
        for label, source, lam in (
            ("frozen_spatial", record["block_proposal"], mapped),
            ("frozen_matched_scalar", record["scalar_proposal"], record["scalar_proposal"]["lambda"]),
        ):
            operator = NativeProfile(spec, sensor, incident, data, record["cap"], source["mu"]/FACTOR, lam/FACTOR)
            result = operator.diagnostic(radius)
            case_record["native"][label] = result
            print(name, label, result["numeric_inertia"], result["distance"], result["margin"],
                  result["wall_seconds"], flush=True)
        rows[name] = case_record
        (DATA/"fine_spatial_followup.json").write_text(json.dumps(
            dict(scope="Numerical diagnostics on all 9216 native generation-grid cells",
                 mapping="Repeat each frozen coarse optical block on a 24×24 fine block; μ and λ divided by 9",
                 cases=rows), indent=2, allow_nan=False)+"\n")


def strict_replay():
    from scipy.linalg import cholesky, eigh
    from flint import acb, acb_mat, arb, ctx

    sys.path.insert(0, str(STUDY/"core"))
    from certificate import ball, frobenius_squared, gershgorin, safe_float, star

    started, cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = 112, 2
    name = CASES[0]
    spec, s, b, y = read(name, "PW180_10x20")
    row = json.loads((DATA/"fine_spatial_followup.json").read_text())["cases"][name]
    proposal = row["coarse_block_proposal"]
    mu, cap = proposal["mu"]/FACTOR, proposal["u"]
    lam = repeated_blocks(proposal)/FACTOR
    native = NativeProfile(spec, s, b, y, cap, mu, lam)
    numeric = native.diagnostic(row["radius"])
    if numeric["status"] != "numerical_crossing":
        raise ArithmeticError("This replay did not suggest a physical crossing.")

    padding = len(spec)
    beta, shift, alpha, gap = map(arb, (.875, .234375, .6875, 1./256))
    scale, upper = arb(mu), arb(cap)
    lo, hi = arb(float(lam.min()))/scale, arb(float(lam.max()))/scale
    baseline, contrast = (lo+hi)/2, (hi-lo)/2
    signs = np.where(lam == lam.max(), 1, -1)
    if len(np.unique(lam)) != 2 or len(spec.shape) != 2 or spec.shape[0] != spec.shape[1]:
        raise ValueError("The Fourier enclosure expects exactly two real weight levels.")
    block_signs = signs.reshape(SIDE, SIDE)[::SIDE//4, ::SIDE//4]
    qmax = 10
    modes_1d = list(range(-qmax, qmax+1))
    qs = [(j % padding, k % padding) for j in modes_1d for k in modes_1d]
    keep = set(qs)
    radius_comp = alpha + beta*contrast**2
    raw = [[acb(float(z.real), float(z.imag)) for z in r] for r in spec]
    omitted_minimum = arb(1000)
    for j, spec_row in enumerate(raw):
        for k, z in enumerate(spec_row):
            if (j,k) in keep:
                continue
            symbol_slack = 1-upper*z.real+baseline*z.imag-radius_comp
            gamma2 = (z.real-shift)**2+z.imag**2
            slack = (symbol_slack-gamma2/(4*beta)).lower()
            omitted_minimum = min(omitted_minimum, slack)
    if not omitted_minimum > 0:
        raise ArithmeticError("An omitted Fourier symbol is not positive.")

    # The two fields in each Fourier symbol are E x and E (M x), M=±1.
    # Their exact full Fourier norms agree. The omitted 2×2 blocks are
    # positive by the checked Schur complement, giving A >= alpha I+F C F*.
    roots = [acb(arb(2)*k/padding).exp_pi_i() for k in range(padding)]
    diffs = {(j-i) % padding for i in modes_1d for j in modes_1d}
    block_len = SIDE//4
    sums = {delta: [sum((roots[delta*k % padding] for k in range(j*block_len,(j+1)*block_len)), acb(0))
                    for j in range(4)] for delta in diffs}
    simple = {delta: sum(values, acb(0)) for delta, values in sums.items()}
    cross = {(dy,dx): sum((int(block_signs[j,k])*sums[dy][j]*sums[dx][k]
                           for j in range(4) for k in range(4)), acb(0))/padding**2
             for dy in diffs for dx in diffs}
    outer = {(dy,dx): simple[dy]*simple[dx]/padding**2 for dy in diffs for dx in diffs}
    f = acb_mat(len(modes_1d), SIDE, [roots[-k*j % padding] for k in modes_1d for j in range(SIDE)])
    shat, signed_shat = [], []
    sh = s.conj().T
    shb = ball(sh)
    for column in sh.T:
        orig = (f*ball(column.reshape(SIDE,SIDE))*f.transpose())/padding
        signed = (f*ball((signs*column).reshape(SIDE,SIDE))*f.transpose())/padding
        shat.append(orig.entries())
        signed_shat.append(signed.entries())
    ss = star(shb)*shb
    r, rec = len(qs), len(s)
    n = 2*r + rec
    unity = alpha-gap
    small = acb_mat(n,n)
    for i,(yi,xi) in enumerate(qs):
        for j,(yj,xj) in enumerate(qs):
            dy, dx = (yj-yi)%padding,(xj-xi)%padding
            standard = outer[dy,dx]/unity
            weighted = cross[dy,dx]/unity
            small[i,j] = small[i+r,j+r] = standard
            small[i,j+r] = small[i+r,j] = weighted
        for j in range(rec):
            entry0 = shat[j][i]/unity
            entry1 = signed_shat[j][i]/unity
            small[i,2*r+j], small[2*r+j,i] = entry0, entry0.conjugate()
            small[i+r,2*r+j], small[2*r+j,i+r] = entry1, entry1.conjugate()
    for i in range(rec):
        for j in range(rec):
            small[2*r+i,2*r+j] = ss[i,j]/unity
        small[2*r+i,2*r+i] += scale
    positives = rec
    for i,(j,k) in enumerate(qs):
        z = raw[j][k]
        coeff0 = 1-upper*z.real+baseline*z.imag-radius_comp
        coeff1 = beta*contrast**2
        coupling = acb(0,contrast)*(z-shift).conjugate()/2
        det = coeff0*coeff1-coupling.real**2-coupling.imag**2
        if not (det < 0 or det > 0):
            raise ArithmeticError("Retained 2×2 Fourier block inversion is unresolved.")
        positives += 1+int(det>0)
        small[i,i] += coeff1/det
        small[i+r,i+r] += coeff0/det
        small[i,i+r] -= coupling/det
        small[i+r,i] -= coupling.conjugate()/det
    gram_seconds = perf_counter()-started

    mid = h(np.array([[complex(x) for x in r] for r in small.tolist()]))
    values, modes = eigh(mid, check_finite=False, driver="evd")
    negatives = n-positives
    if int(np.count_nonzero(values < 0)) != negatives:
        raise ArithmeticError("The Fourier lower envelope did not give numerical positive inertia.")
    step = float(np.min(np.abs(values))/2)
    negative = modes[:,:negatives]
    w = negative*np.sqrt(-2*values[:negatives])
    lifted = mid+w@w.conj().T
    lifted.flat[::n+1] -= step
    ff = ball(np.tril(cholesky(lifted,lower=True,check_finite=False)))
    wb, nb = ball(w), ball(negative)
    residual = small+wb*star(wb)-ff*star(ff)
    for i in range(n):
        residual[i,i] -= arb(step)
    defect = frobenius_squared(residual).sqrt().upper()
    neg_bound = gershgorin(star(nb)*(small*nb+arb(step)*nb),upper=True)
    if not (arb(step)>defect and neg_bound<0):
        raise ArithmeticError("Outward small-matrix inertia unresolved.")
    inertia_seconds = perf_counter()-started-gram_seconds

    # H needs only its two useful upper eigen-directions. Evaluate the
    # full, untruncated FFT using the represented source/current proposals.
    _, hv = eigh(native.source_form, check_finite=False)
    v = hv[:, -2:][:, ::-1]
    vb = ball(v)
    xb = ball(native.solution @ v)
    yv, bv = ball(y)*vb, ball(b)*vb
    mult = [acb(upper,baseline+int(x)*contrast)/2 for x in signs]

    def mult_rows(m, scalars):
        return acb_mat(N,m.ncols(),[u*scalars[i] for i in range(N) for u in
                                   (m[i,j] for j in range(m.ncols()))])

    def exact_fft_apply(x, *, adjoint=False):
        z = acb(0)
        cols = []
        for c in range(x.ncols()):
            rows = [acb.dft([x[SIDE*i+j,c] for j in range(SIDE)]+[z]*(padding-SIDE)) for i in range(SIDE)]
            transformed = []
            for j in range(padding):
                d = acb.dft([rows[i][j] for i in range(SIDE)]+[z]*(padding-SIDE))
                filtered = [d[i]*(raw[i][j].conjugate() if adjoint else raw[i][j]) for i in range(padding)]
                transformed.append(acb.dft(filtered, inverse=True)[:SIDE])
            native_rows = [acb.dft([transformed[j][i] for j in range(padding)],inverse=True)[:SIDE]
                           for i in range(SIDE)]
            cols.append([k for r in native_rows for k in r])
        return acb_mat(x.ncols(),N,sum(cols,[])).transpose()

    az = (xb+shb*(star(shb)*xb)/scale-mult_rows(exact_fft_apply(xb),mult)
          -exact_fft_apply(mult_rows(xb,[x.conjugate() for x in mult]),adjoint=True))
    rhs = shb*yv/scale+mult_rows(bv,mult)
    defect_source = frobenius_squared(az-rhs).upper()
    h0 = star(yv)*yv+scale*(star(xb)*az-star(rhs)*xb-star(xb)*rhs)
    correction = (scale*defect_source/gap).upper()
    gram = star(vb)*vb
    retained = []
    for k in range(1,3):
        projected = acb_mat(k,k,[h0[i,j] for i in range(k) for j in range(k)])
        gram1 = acb_mat(k,k,[gram[i,j] for i in range(k) for j in range(k)])
        den = gershgorin(gram1,upper=True)
        num = gershgorin(projected,upper=False)
        # The stationary error was bounded in the trial basis already.
        floor = ((num-correction)/den).lower()
        retained.append(max(arb(0),floor))
    lower = sum(retained,arb(0)).sqrt()
    target = arb(row["radius"])+arb("0.000001")*frobenius_squared(ball(y)).sqrt()
    out = dict(case=name, cap=cap, mu=mu, scope="9216 native cells, represented binary64 generating FFT symbol and ports",
        profile="frozen spatial coarse proposal divided by 9", numeric_distance=numeric["distance"],
        outward_distance_lower=safe_float(lower,upper=False), frozen_noise=row["radius"],
        noise_with_generation_allowance_upper=safe_float(target,upper=True),
        margin_with_generation_allowance_lower=safe_float(lower-target,upper=False),
        excluded_at_padded_radius=bool(lower>target),
        spectral_floor_scaled_K_verified=float(gap), omitted_symbol_slack_lower=safe_float(omitted_minimum,upper=False),
        Fourier_cutoff=qmax, alpha=float(alpha), beta=float(beta), real_symbol_shift=float(shift),
        proof_dimension=n, proof_negatives=negatives, verified_native_K_negatives=0,
        small_inertia_defect_upper=safe_float(defect,upper=True),
        small_inertia_shift=step, source_stationary_squared_defect_upper=safe_float(defect_source,upper=True),
        source_schur_correction_upper=safe_float(correction,upper=True),
        source_eigenvalue_lowers=[safe_float(x,upper=False) for x in retained],
        precision_bits=ctx.prec, gram_wall_seconds=gram_seconds, inertia_wall_seconds=inertia_seconds,
        wall_seconds=perf_counter()-started, cpu_seconds=process_time()-cpu)
    (DATA/"fine_spatial_verified.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(name,"outward",out["outward_distance_lower"],out["margin_with_generation_allowance_lower"],
          out["wall_seconds"],flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    args = parser.parse_args()
    with threadpool_limits(1):
        strict_replay() if args.verify else main()
