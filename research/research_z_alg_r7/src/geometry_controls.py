"""Equal-power nonlinear illumination controls and a receiver-darkness stress test."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import null_space, solve_triangular
from scipy.optimize import brentq, least_squares
from scipy.special import j0
from scipy.stats import qmc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
Q0 = 15.25
DELTA = 2.
CAP = 4.
N = 6
I = np.eye(N)
METHODS = ["matched", "focus", "local_sensitivity", "Jacobian_null", "random"]
APERTURES = ["ring", "arc120", "arc60"]
NAMES = ["Matched", "Host focus", "Local sensitivity", "Jacobian null", "Random"]
COLORS = ["#007181", "#bb651d", "#584ea3", "#4c8a45", "#888888"]


def physics_module(path):
    if path is None:
        path = ROOT.parent / "research_z_alg_r2/src/physics.py"
    spec = importlib.util.spec_from_file_location("pulse_tm_geometry", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_bank(aperture):
    span = {"arc120": 2*np.pi/3, "arc60": np.pi/3}.get(aperture)
    t = .137 + (np.arange(8)*np.pi/4 if span is None else np.linspace(-span/2, span/2, 8))
    return 1.05*np.column_stack([np.cos(t), np.sin(t)])


def jacobian_rows(model, state):
    receive = np.linalg.solve(state["L"].T, model.S.T).T
    return receive[:, :N, None] * state["E"][None, :N, :], receive[:, :N]


def power(a, gram):
    return np.vdot(a, gram@a).real


def whitening(gram):
    lower = np.linalg.cholesky(gram)
    return solve_triangular(lower.conj().T, np.eye(len(gram)), lower=False)


def target_drive(B, b, drive_weight):
    W = whitening(drive_weight)
    return W @ np.linalg.lstsq(B @ W, b, rcond=None)[0]


def matched_capped(B, b, gram):
    """Maximize the common scale under radiation power <=1 and drive norm <=4."""
    identity = np.eye(len(gram))
    def optimum(lam):
        a = target_drive(B, b, gram + lam * identity)
        return a / np.sqrt(power(a, gram))

    a = optimum(0.)
    if np.linalg.norm(a) <= CAP:
        return a, "radiation"
    a = target_drive(B, b, identity)
    a *= CAP / np.linalg.norm(a)
    if power(a, gram) <= 1:
        return a, "drive"
    lam = brentq(lambda t: np.linalg.norm(optimum(np.exp(t))) - CAP, -35., 35.)
    return optimum(np.exp(lam)), "both"


def designs(i, reduced, jac, gram, rng):
    B, H = reduced.q, reduced.K
    b = I[:, i] - Q0*H[:, i]
    W = whitening(gram)
    def normalized(a):
        return a / np.sqrt(power(a, gram))
    def align(a, reference):
        return a * np.exp(-1j*np.angle(np.vdot(reference, gram @ a)))

    matched = normalized(target_drive(B, b, gram))
    focus = normalized(target_drive(B, I[:, i], gram))
    sensitivity = normalized(W @ (jac[i] @ W).conj())
    nuisance = np.arange(N) != i
    null = null_space(jac[nuisance] @ W)
    aggressive = normalized(W @ null @ (jac[i] @ W @ null).conj())
    aggressive = align(aggressive, matched)
    random = normalized(rng.normal(size=B.shape[1])+1j*rng.normal(size=B.shape[1]))
    capped, constraint = matched_capped(B, b, gram)
    return np.column_stack([matched, focus, sensitivity, aggressive, random]), capped, constraint


def physical(model, z, host):
    m = host.copy()
    m[:N] = z
    return model.forward(m, independent=True)


def find_dark(model, reduced, host, rng):
    """Seek a single receiver zero with two varying nuisance cells and three fixed ones."""
    q = np.full(N, Q0)
    attempts = []
    def visibility(z):
        return np.linalg.solve((I-z[:, None]*reduced.K).T, reduced.L.T).T
    # Only two nuisance values vary per solve. Keep a reproducible bounded search log.
    witnesses = []
    for i in range(N):
        for pair in itertools.combinations(np.delete(np.arange(N), i), 2):
            indices = list(pair)
            def vector(x):
                z = q.copy()
                z[indices] = x
                v = visibility(z)[0, i]
                return np.array([v.real, v.imag])
            best = None
            for x0 in [np.array([Q0]*2), *rng.uniform(.5, 30., (2, 2))]:
                opt = least_squares(vector, x0, bounds=(.5,30.), gtol=1e-15, ftol=1e-15, xtol=1e-15,
                                    max_nfev=120)
                err = np.linalg.norm(vector(opt.x))
                if best is None or err < best[0]:
                    best = (err, opt.x, opt.jac)
                if err < 1e-12:
                    break
            err, x, deriv = best
            attempts.append([i, *indices, *x, err])
            if err < 1e-12:
                z = q.copy()
                z[indices] = x
                lo, hi = z.copy(), z.copy()
                lo[i], hi[i] = .5, 30.
                ylo, yhi = physical(model, lo, host)["y"], physical(model, hi, host)["y"]
                witnesses.append(dict(cell=int(i), varying=[int(j) for j in indices], z=z.tolist(),
                                      abs_visibility=np.abs(visibility(z)[:, i]).tolist(),
                                      nuisance_jacobian_singular_values=np.linalg.svd(deriv,compute_uv=False).tolist(),
                                      max_abs_feed_response_gap=np.max(abs(yhi-ylo),axis=-1).tolist(),
                                      ylo=np.stack([ylo.real,ylo.imag], axis=-1).tolist(),
                                      yhi=np.stack([yhi.real,yhi.imag], axis=-1).tolist()))
    return np.array(attempts), witnesses


def run(physics):
    started = time.perf_counter()
    cpu = time.process_time()
    geom = json.loads((ROOT.parent / "research_z_alg_r5_measurement/amendment_codex/results/H12_contract.json").read_text())
    centers = np.array(geom["centers"])
    host = np.array(geom["known_host"])
    receivers = np.array([[.17, .86], [-.17, -.86]])
    k = 2*np.pi
    zref = np.full(N, Q0)
    quasi = qmc.Sobol(N-1,scramble=True,seed=18112026).random_base2(8)*29.5 + .5
    rng = np.random.default_rng(2026111807)
    output = dict(nuisance=quasi, receivers=receivers, centers=centers, host=host,
                  q=Q0, delta=DELTA, methods=METHODS, apertures=APERTURES)
    summary = {"setup": "8 coherent external feeds, two receiver masks, disk radius .04, k=2pi",
               "reference": Q0, "positive_shift": DELTA, "nuisance_per_cell": len(quasi),
               "constraint_cap": CAP, "summaries": {}}
    for aperture in APERTURES:
        sources = source_bank(aperture)
        gram = j0(k*np.linalg.norm(sources[:,None]-sources[None,:], axis=-1))
        model = physics.Model(centers=centers, sources=sources, receivers=receivers, radius=.04, k=k)
        reduced = physics.ReducedKernel(model, host, active=np.arange(N))
        ref_state = physical(model, zref, host)
        jac, nominal_v = jacobian_rows(model, ref_state)
        ref = ref_state["y"]
        drives, caps, caps_status = [], [], []
        # Axis order: tested cell, sample, receiver, primitive feed. Explicit full current solves.
        y = np.empty((N, len(quasi), 2, 2, 8), complex)
        direct_identity = []
        for i in range(N):
            A, capped, constraint = designs(i, reduced, jac[0], gram, rng)
            drives.append(A)
            caps.append(capped)
            caps_status.append(constraint)
            other = np.arange(N) != i
            for s, nuisance in enumerate(quasi):
                z = zref.copy()
                z[other] = nuisance
                for shift, val in enumerate([Q0, Q0+DELTA]):
                    z[i] = val
                    state = physical(model, z, host)
                    y[i,s,shift] = state["y"]
                if s % 16 == 0:
                    v = np.linalg.solve((I-z[:,None]*reduced.K).T, reduced.L.T).T[:,i]
                    b = I[:,i] - Q0*reduced.K[:,i]
                    alpha = np.vdot(b,reduced.q@A[:,0])/np.vdot(b,b)
                    direct_identity.append(np.max(abs((state["y"]-ref)@A[:,0] - alpha*DELTA*v)))
        drives, caps = np.array(drives), np.array(caps)
        centered = np.einsum("iszrf,ifm->iszrm", y-ref, drives)
        paired = centered[:,:,1]-centered[:,:,0]
        capped_diff = np.einsum("isrf,if->isr", y[:,:,1]-y[:,:,0], caps)
        drift = np.abs(centered[:,:,0,0,:])
        sep = np.abs(paired[:,:,0,:])
        power_errors = [abs(power(a,gram)-1.) for bank in drives for a in bank.T]
        sample_summary = {}
        for j, name in enumerate(METHODS):
            sample_summary[name] = dict(null_drift_q50=float(np.median(drift[...,j])),
                                        null_drift_max=float(np.max(drift[...,j])),
                                        pair_gap_quantiles=np.quantile(sep[...,j],[0,.05,.5,.95,1]).tolist(),
                                        altered_centered_q50=float(np.median(abs(centered[:,:,1,0,j]))))
        alphas = []
        for i in range(N):
            b = I[:,i]-Q0*reduced.K[:,i]
            alphas.append(abs(np.vdot(b,reduced.q@drives[i,:,0])/np.vdot(b,b)))
        summary["summaries"][aperture] = dict(
            methods=sample_summary, max_power_error=max(power_errors),
            max_direct_identity_abs=max(direct_identity),
            max_relative_null_drive_gap=max(np.linalg.norm(drives[:,:,0]-drives[:,:,3],axis=1)/
                                            np.linalg.norm(drives[:,:,0],axis=1)),
            min_nominal_receiver_visibility=float(min(abs(nominal_v[0]))),
            target_field_scale=alphas, matched_feed_norm=np.linalg.norm(drives[:,:,0],axis=1).tolist(),
            capped_constraints=caps_status, capped_feed_norm=np.linalg.norm(caps,axis=1).tolist(),
            capped_power=[power(a,gram) for a in caps], capped_pair_median=float(np.median(abs(capped_diff[:,:,0]))),
            capped_to_power_only_median_by_cell=np.median(abs(capped_diff[:,:,0])/sep[:,:,0],axis=1).tolist(),
            forward_count=model.counts["INDEPENDENT_FORWARD"],
            condition_B=float(np.linalg.cond(reduced.q)), condition_gram=float(np.linalg.cond(gram)))
        output.update({aperture+"_"+key: value for key, value in dict(
            sources=sources, gram=gram, drives=drives, capped_drives=caps, y_ref=ref, y=y,
            centered=centered, paired=paired, capped_pair=capped_diff).items()})
        if aperture == "ring":
            attempts, dark = find_dark(model, reduced, host, np.random.default_rng(1962))
            output["dark_search_columns"] = ["cell","vary1","vary2","z1","z2","abs_v0"]
            output["dark_search"] = attempts
            summary["receiver_dark_witnesses"] = dark
            summary["dark_search_full_solves"] = 2*len(dark)

    np.savez_compressed(ROOT/"results/geometry_controls_raw.npz", **output)
    figures(output)
    summary["wall_seconds_with_plots"] = time.perf_counter()-started
    summary["cpu_seconds_with_plots"] = time.process_time()-cpu
    (ROOT/"results/geometry_controls.json").write_text(json.dumps(summary, indent=2)+"\n")


def figures(data):
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    fig, axes = plt.subplots(1,3,figsize=(13.8,4.3), layout="constrained")
    for j, name in enumerate(NAMES):
        for axis, values in [(axes[0],abs(data["ring_centered"][:,:,0,0,j])),
                             (axes[1],abs(data["ring_paired"][:,:,0,j]))]:
            v=np.sort(np.maximum(values.ravel(),1e-18))
            axis.plot(v, np.arange(1,len(v)+1)/len(v), color=COLORS[j], label=name,
                      ls="--" if j==3 else "-", lw=1.65)
    for axis, title, xlabel in zip(axes[:2],["A  Fixed target, varying surroundings","B  Target shifted by +2"],
                                   ["Deviation from nominal null reading","Paired change at same surroundings"]):
        axis.set(xscale="log",ylabel="Empirical cumulative fraction",xlabel=xlabel,title=title,ylim=(0,1.02))
        axis.grid(True,alpha=.18)
    axes[0].legend(fontsize=8, frameon=False, loc="lower right")
    for aperture,style in zip(APERTURES,["-","--",":"]):
        raw=abs(data[aperture+"_paired"][:,:,0,0])
        cap=abs(data[aperture+"_capped_pair"][:,:,0])
        for capname, curve, color in [("radiation",raw,"#007181"),("radiation + drive",cap,"#a44d2b")]:
            if capname == "radiation + drive" and np.allclose(cap, raw):
                continue
            axes[2].plot(np.arange(1,N+1),np.median(curve,axis=1),"o"+style,color=color,
                         label=f"{aperture}, {capname}")
    axes[2].set(yscale="log",xlabel="Target cell (1–6)",ylabel="Median paired change",
                title="C  Matched drive under aperture limits", xticks=range(1,N+1))
    axes[2].legend(fontsize=8,frameon=False)
    axes[2].grid(True,alpha=.18)
    for ext in ["pdf","png"]:
        fig.savefig(ROOT / f"figures/geometry_controls.{ext}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--physics", type=Path, help="R2 src/physics.py (default: sibling checkout)")
    run(physics_module(parser.parse_args().physics))
