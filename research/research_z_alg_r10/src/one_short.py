"""Find and separate physical response branches with one unmeasured source direction left."""

import json
import time

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares, minimize
from scipy.special import j0

from lift_active import (
    BOUNDS, ROOT, Model, ReducedKernel, Response, affine_response, error, features,
    full_data, json_value, lifted_coefficients, material_from_row, pack, random_basis,
)
from noise_sets import EPS, GAIN, noisy_vector


def physical_row(H, R, z):
    return np.linalg.solve((np.eye(len(z)) - z[:, None] * H).T, R) * z


def one_line(H, R, B, y, noisy):
    """Inspect the exact affine lift, or use its closest rank-seven line for starts."""
    w0, N, bsing, _ = affine_response(B, y)
    if N.shape[1] != 2:
        raise ValueError("Expected exactly two missing complex source directions")
    c, M = lifted_coefficients(w0, N, H, R)
    tscale = max(np.linalg.norm(w0) / np.sqrt(len(w0)), 1e-10)
    fscale = np.r_[np.full(4, tscale), np.full(4, tscale**2)]
    Q = M * fscale
    U, s, V = np.linalg.svd(Q, full_matrices=False)
    rank = int(np.sum(s > s[0]*1e-11))
    origins = []
    roots = None
    if rank == 8 and (not noisy or s[7] > s[6]*.10):
        full = V.T @ ((U.T @ -c) / s)
        origins.append(("full_lift", full))
    if rank == 7 or noisy:
        origin = V[:7].T @ ((U[:, :7].T @ -c) / s[:7])
        direction = V[7]
        dt = direction[:2] + 1j*direction[2:4]
        j = int(np.argmax(abs(dt)))
        t = origin[j] + 1j*origin[j+2]
        poly = [abs(dt[j])**2, 2*np.real(t.conjugate()*dt[j])-direction[4+j],
                abs(t)**2-origin[4+j]]
        # A homogeneous change in Hermitian features with dt=0 gives a linear
        # intersection condition. A vanishing scalar constraint leaves this chart inconclusive.
        if poly[0] <= 1e-24:
            poly = poly[1:]
        roots = [float(root.real) for root in np.roots(poly) if abs(root.imag) < 1e-10]
        origins.extend(("lift_line", origin + rho*direction) for rho in roots)

    candidates = []
    for source, vector in origins:
        t = tscale*(vector[:2] + 1j*vector[2:4])
        w = w0 + N @ t
        denom = R+w @ H
        if np.min(abs(denom)) < 1e-13*np.max(abs(denom)):
            continue
        z = material_from_row(w, H, R)
        lie = features(t)/fscale - vector
        feature_defect = np.linalg.norm(lie)/max(np.linalg.norm(vector), 1e-10)
        real_defect = np.linalg.norm(z.imag)/max(np.linalg.norm(z), 1e-10)
        is_legal = np.all((z.real >= BOUNDS[0]) & (z.real <= BOUNDS[1]))
        mat = np.clip(z.real, *BOUNDS)
        actual = physical_row(H, R, mat)
        response_defect = np.linalg.norm(actual @ B-y)/max(np.linalg.norm(y), 1e-10)
        if noisy or (feature_defect < 1e-7 and real_defect < 1e-7 and response_defect < 1e-8 and is_legal):
            candidates.append(dict(source=source, material=mat, actual=actual, row=w,
                                   feature_defect=feature_defect, real_defect=real_defect,
                                   initial_response_defect=response_defect,
                                   minimum_resolvent_singular=np.linalg.svd(
                                       np.eye(len(mat))-mat[:, None]*H,compute_uv=False)[-1]))
    diag = dict(lift_rank=rank, lift_relative_singular=s/s[0], incident_relative_singular=bsing/bsing[0],
                line_roots=roots, algebraic_start_count=len(candidates))
    return candidates, N, diag


def gain_profile(predicted, observed):
    residual = observed-predicted
    g = np.vdot(predicted, residual)/(np.vdot(predicted, predicted).real+(EPS/GAIN)**2)
    value = np.linalg.norm((residual-g*predicted)/EPS)**2+abs(g/GAIN)**2
    return value, g


def fit_physical(kernel, feeds, host, observed, initial):
    """Physical data score with one receiver gain and an additive Euclidean noise ball."""
    response = Response(kernel, feeds)
    k = len(initial)

    def fun(x):
        predicted = response.predict(x[:k])+host
        g = GAIN*(x[k]+1j*x[k+1])
        return np.r_[pack(((1+g)*predicted-observed)/EPS), x[k:]]

    def jac(x):
        predicted = response.predict(x[:k])+host
        D = response.predict(x[:k], derivative=True)
        g = GAIN*(x[k]+1j*x[k+1])
        J = (1+g)*(D[::2]+1j*D[1::2])
        J = np.column_stack((J, GAIN*predicted, 1j*GAIN*predicted))/EPS
        return np.vstack((np.stack((J.real, J.imag), axis=1).reshape(-1,k+2),
                          np.c_[np.zeros((2,k)),np.eye(2)]))

    _, g = gain_profile(response.predict(initial)+host, observed)
    lower = np.r_[np.full(k, BOUNDS[0]), -np.inf, -np.inf]
    upper = np.r_[np.full(k, BOUNDS[1]), np.inf, np.inf]
    fit = least_squares(fun, np.r_[initial, g.real/GAIN, g.imag/GAIN], jac=jac, bounds=(lower, upper),
                        max_nfev=200, gtol=1e-11, ftol=1e-11, xtol=1e-11)
    return fit.x[:k], np.linalg.norm(fit.fun)**2, response.nlu


def old_candidates(kernel, feeds, host, observed, noisy, *, gain_sweep=False):
    begin = time.process_time()
    H, R = kernel.K, kernel.L.ravel()
    B = kernel.q @ feeds
    candidates, N, diagnostic = one_line(H, R, B, observed-host, noisy)
    recovered = []
    if noisy and gain_sweep:
        # The sweep uses only the known receiver-gain radius and the current data.
        trial_gains = GAIN*np.ravel([rho*np.exp(1j*np.arange(12)*np.pi/6)
                                   for rho in (.25,.50,.75,1.)])
        candidates = list(candidates)
        for gain in trial_gains:
            approximate_clean = observed/(1+gain)-host
            starts, _, _ = one_line(H,R,B,approximate_clean,True)
            for c in starts:
                if not any(error(c["material"], prior["material"]) < .025 for prior in candidates):
                    candidates.append(dict(c, assumed_gain_for_seed=[gain.real,gain.imag]))
        diagnostic["gain_sweep_starts"] = len(candidates)
    for c in candidates:
        z = c["material"]
        if noisy:
            z, score, nlu = fit_physical(kernel, feeds, host, observed, z)
        else:
            score, _ = gain_profile(c["actual"] @ B+host, observed)
            nlu = 0
        old = [r for r in recovered if error(r["material"], z) < 1e-4]
        diagnostic.setdefault("polish_scores", []).append(score)
        if score <= 1.0 and not old:
            recovered.append(dict(c, material=z, actual=physical_row(H,R,z), joint_score=score, nlu=nlu))
    diagnostic.update(accepted=len(recovered), cpu_s=time.process_time()-begin)
    return recovered, N, diagnostic


def normalize(x):
    return x/np.linalg.norm(x)


def geometry_e_opt(B, P):
    """Globally maximize sigma_min([B,Pc]) over ||c||=1 by Schur bisection."""
    u, s, _ = np.linalg.svd(B, full_matrices=False)
    proj = u.conj().T @ P
    leftover = P-u @ proj
    G = leftover.conj().T @ leftover
    hi = s[-1]**2
    lo = 0.
    c = np.linalg.eigh(G)[1][:,-1]

    def actual(c):
        return np.linalg.svd(np.column_stack((B,P@c)),compute_uv=False)[-1]

    best = actual(c)
    # At the cap the coupling into the weakest old singular space must vanish.
    # Restrict there before evaluating the finite Schur complement.
    weak = s**2-hi <= 1e-10*hi
    Z = null_space(proj[weak])
    if len(Z.T):
        regular = proj[~weak] @ Z
        capped = Z.conj().T @ G @ Z - regular.conj().T @ ((hi/(s[~weak]**2-hi))[:,None]*regular)
        at_cap = Z @ np.linalg.eigh(capped)[1][:,-1]
        actual_cap = actual(at_cap)
        if actual_cap >= (1-1e-10)*s[-1]:
            return at_cap
        if actual_cap > best:
            best, c = actual_cap, at_cap

    for _ in range(43):
        tau = (hi+lo)/2
        values, vectors = np.linalg.eigh(G-proj.conj().T @ ((tau/(s**2-tau))[:,None]*proj))
        proposal = vectors[:,-1]
        separation = actual(proposal)
        if separation > best:
            best, c = separation, proposal
        if values[-1] >= tau:
            lo = tau
        else:
            hi = tau
    return c


def a_opt_fisher(kernel, p_old, source_modes, cands, starts):
    """Local A-optimality, averaged over old-data physical candidates, including gain."""
    if not cands:
        return starts[0]
    k = kernel.r
    ns = source_modes.shape[1]
    all_feeds = np.column_stack((p_old,source_modes))
    f = Response(kernel, all_feeds)
    host = (kernel.y0 @ all_feeds).ravel()
    blocks = []
    for c in cands:
        z = c["material"]
        yy = f.predict(z)+host
        dj = f.predict(z, derivative=True)
        nuisance = np.column_stack((GAIN*yy, 1j*GAIN*yy))
        j = np.column_stack((dj[::2]+1j*dj[1::2],nuisance))/EPS
        jo = j[:p_old.shape[1]]
        prior = np.diag(np.r_[np.zeros(k),1.,1.])
        normal = jo.real.T @ jo.real+jo.imag.T @ jo.imag+prior
        covariance = np.linalg.inv(normal)
        blocks.append((j[p_old.shape[1]:], covariance, covariance[:,:k]@covariance[:k,:]))

    def reduction(c):
        center = normalize(c[:ns]+1j*c[ns:])
        ans = 0.0
        for J, C, Cz in blocks:
            x = center @ J
            A = np.vstack((x.real,x.imag))
            U = A @ C @ A.T
            val = A @ Cz @ A.T
            ans += np.trace(np.linalg.solve(np.eye(2)+U, val))
        return ans

    scale = max(reduction(np.r_[x.real,x.imag]) for x in starts)
    def objective(c):
        return -reduction(c)/max(scale,1e-12)

    fits = [minimize(objective,np.r_[c.real,c.imag],method="L-BFGS-B",
                     options={"maxiter":40,"maxfun":2500}) for c in starts]
    best = min(fits,key=lambda fit:fit.fun).x
    return normalize(best[:ns]+1j*best[ns:])


def protocols(kernel, modes, old, cands, N, random):
    P = kernel.q @ modes
    B = kernel.q @ old
    g = geometry_e_opt(B,P)
    if len(cands) >= 2:
        responses = [c["actual"] @ P for c in cands]
        gap = max((a-b for i,a in enumerate(responses) for b in responses[i+1:]),key=np.linalg.norm)
        c = normalize(gap.conj())
    else:
        # There is no discovered ambiguity to target.
        c = g
    f = a_opt_fisher(kernel,old,modes,cands,[g,random])
    projection = np.linalg.svd(N.T @ P,full_matrices=False)[2][0].conj()
    radios = dict(branch=modes @ c, geometry_e_opt=modes @ g, local_fisher=modes @ f,
                  random=modes @ random, repeat_old=old[:,0], geometry_projection=modes @ projection)
    return radios


def load_case(k):
    file = np.load(ROOT/"results"/f"lift_k{k}.npz")
    model = Model(centers=file["centers"], sources=file["sources"], receivers=file["receivers"],
                  radius=float(file["disk_radius"]), k=float(file["wavenumber"]))
    return file, ReducedKernel(model,file["host"],active=file["active"])


def observe(kernel, total, old, noise, g, a):
    feeds = np.column_stack((old,a))
    B = kernel.q @ feeds
    host = (kernel.y0 @ feeds).ravel()
    nominal = total @ feeds
    y = (1+g)*nominal + noise
    return feeds, B, host, nominal, y


def current_case(k, source, kernel, sample, noisy, *, save, gain_sweep=False):
    H,R = kernel.K,kernel.L.ravel()
    old = source["p_old"][:,:-1]
    modes = source["observable_modes"]
    truth = np.vstack((source["pair_materials"],source["heldout_materials"]))
    total = source["independent_all_source_data"][sample]+source["y_host"]
    rng = np.random.default_rng(394760+100*k+sample)
    noise_raw = noisy_vector(rng,k)*EPS*.55
    graw = .8*GAIN*np.exp(2j*np.pi*rng.uniform())
    g = graw if noisy else 0.
    noise = np.r_[noise_raw[:-2],noise_raw[-1]] if noisy else np.zeros(k-1)
    old_nominal = total @ old
    old_y = (1+g)*old_nominal+noise[:-1]
    old_host = (kernel.y0 @ old).ravel()
    cands,N,diag = old_candidates(kernel,old,old_host,old_y,noisy,gain_sweep=gain_sweep)
    rnd = np.random.default_rng(8871+k+1000*sample)
    rnd = noisy_vector(rnd,k)
    selected_at = time.process_time()
    beams = protocols(kernel,modes,old,cands,N,rnd)
    select_cpu = time.process_time()-selected_at
    results, records = {}, []
    for name,a in beams.items():
        p, B, hy, y0, y = observe(kernel,total,old,noise,g,a)
        score = [gain_profile(physical_row(H,R,z) @ B+hy,y)[0] for z in truth[:2]]
        old_candidate = [gain_profile(c["actual"] @ B+hy,y)[0] for c in cands]
        local_scores = [fit_physical(kernel,p,hy,y,c["material"]) for c in cands]
        best = min(local_scores,key=lambda v:v[1]) if local_scores else None
        exact_gap = abs((source["independent_all_source_data"][0]-source["independent_all_source_data"][1]) @ a)
        nominal_candidate_gap = abs((cands[0]["actual"]-cands[1]["actual"]) @ B[:,-1]) if len(cands)==2 else None
        sv = np.linalg.svd(B,compute_uv=False)
        results[name] = dict(
            prototype_exact_gap=exact_gap, nominal_discovered_gap=nominal_candidate_gap,
            incident_sigma_min=sv[-1], incident_condition=sv[0]/sv[-1],
            geometric_completion_projection=np.linalg.norm(N.T @ B[:,-1]),
            power_inner_product_vs_branch=abs(np.vdot(a,source["gram"] @ beams["branch"])),
            unit_radiated_power=np.vdot(a,source["gram"] @ a).real,
            prototype_joint_scores=score, old_candidates_joint_scores=old_candidate,
            refined_from_old_score=[x[1] for x in local_scores],
            rmse_refined=error(best[0],truth[sample]) if best else None)
        records.append((p,B,hy,y0,y))
    if save:
        series = "noisy" if noisy else "clean"
        if gain_sweep:
            series += "_repaired"
        np.savez_compressed(ROOT/"results"/f"short_instance_{k}_{series}.npz",
                            k=k,H=H,R=R,protocols=np.array(list(beams)),feeds=np.stack([r[0] for r in records]),
                            B=np.stack([r[1] for r in records]),host=np.stack([r[2] for r in records]),
                            nominal=np.stack([r[3] for r in records]),observed=np.stack([r[4] for r in records]),
                            z=truth[:2],first_sample_truth=truth[sample],gain=g,noise=noise,eps=EPS,gamma=GAIN,
                            old_candidates=np.array([c["material"] for c in cands]),
                            old_scores=np.array([c["joint_score"] for c in cands]),
                            old_N=N, q=kernel.q, modes=modes, gram=source["gram"])
    energy = None
    if sample in (0,1):
        energy = {}
        for name,a in beams.items():
            rows = []
            other = physical_row(H,R,truth[1-sample])
            for power in np.logspace(-8, 0, 33):
                _, B, hy, _, y = observe(kernel,total,old,noise,g,np.sqrt(power)*a)
                objective = gain_profile(other @ B+hy,y)[0]
                rows.append((power,objective))
            energy[name] = dict(rows=rows, first_excluded_unit_energy=next(
                (power for power,objective in rows if objective>1.),None))
    for c in cands:
        c.update(scored_pair_rmse=[error(c["material"],z) for z in truth[:2]],
                 scored_truth_rmse=error(c["material"],truth[sample]))
        del c["row"],c["actual"]
    return dict(k=k,sample=sample,noisy=noisy,old_shots=k-2,new_shots=k-1,
                diagnostics=diag,candidates=cands,results=results,selection_cpu_s=select_cpu,
                energy_schedule=energy,gain_sweep=gain_sweep)


def third_geometry():
    """One predetermined 18-cell experiment with a 237-degree transmit aperture."""
    k = 18
    rng = np.random.default_rng(17822)
    nx, ny = 6, 3
    xx, yy = np.meshgrid(.51*(np.arange(nx)-2.5), .48*(np.arange(ny)-1))
    cells = np.c_[xx.ravel(), yy.ravel()]+rng.uniform(-.049,.049,(k,2))
    xy = np.array([[a,b] for b in (-.24,.24) for a in np.linspace(-1.02,1.02,5)])
    centers = np.r_[cells,xy]
    host = np.r_[np.zeros(k),rng.uniform(1.3,4.0,len(xy))]
    theta = np.deg2rad(np.linspace(-132.,105.,42))
    sources = np.c_[2.6*np.sin(theta),2.45*np.cos(theta)]
    receivers = np.array([[.51,1.78]])
    model = Model(centers=centers,sources=sources,receivers=receivers,radius=.068,k=2*np.pi)
    kernel = ReducedKernel(model,host,active=np.arange(k))
    gram = j0(model.k*np.linalg.norm(sources[:,None]-sources[None,:],axis=-1))
    g, U = np.linalg.eigh(gram)
    keep = g >= 1/16
    P = U[:,keep]/np.sqrt(g[keep])
    _, s, V = np.linalg.svd(kernel.q @ P,full_matrices=False)
    if len(s) < k or s[-1]/s[0] < 1e-11:
        raise ValueError(f"Third layout cannot realize full source span: {s}")
    modes = P @ V.conj().T
    truth = rng.uniform(1.1,10.9,(14,k))
    delta = np.array([full_data(model,kernel,z) for z in truth])
    # The pair is used only to create the engineered old blind array.
    gap = (delta[0]-delta[1]) @ modes
    p_old = modes @ null_space(gap[None,:]) @ random_basis(rng,(k-1,k-1))
    D = dict(pair_materials=truth[:2],heldout_materials=truth[2:],independent_all_source_data=delta,
             y_host=kernel.y0.ravel(),p_old=p_old,gram=gram,observable_modes=modes,
             retained_q_modes=int(keep.sum()),full_incidence_singular=s,H=kernel.K,R=kernel.L.ravel(),
             q=kernel.q,active_support=kernel.active,host=host,centers=centers,sources=sources,
             receivers=receivers,aperture_degrees=237.,disk_radius=model.radius)
    np.savez_compressed(ROOT/"results"/"short_geometry_18.npz",**D)
    return D,kernel


def main():
    rows = []
    for k in (12,24,18):
        source,kernel = load_case(k) if k in (12,24) else third_geometry()
        for sample in range(14):
            for noisy in (True,False):
                row = current_case(k,source,kernel,sample,noisy,save=sample==0)
                rows.append(row)
                print(k,sample,noisy,row["diagnostics"]["accepted"],
                      [(n, round(v["prototype_exact_gap"],5),round(v["incident_sigma_min"],5))
                      for n,v in row["results"].items()],flush=True)
    repaired = current_case(18,source,kernel,0,True,save=True,gain_sweep=True)
    (ROOT/"results"/"short_gain_sweep.json").write_text(
        json.dumps(repaired,indent=2,default=json_value,allow_nan=False)+"\n")
    (ROOT/"results"/"short_comparison.json").write_text(
        json.dumps(rows,indent=2,default=json_value,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
