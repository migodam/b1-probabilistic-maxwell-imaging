"""Physical comparison of coordinated patch illumination with local anchor designs."""
from physical_experiment import Model, ReducedKernel, Response, encode, pack

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial import polynomial as poly
from scipy.linalg import null_space
from scipy.optimize import least_squares, minimize
from scipy.special import j0

from patch_pencil import (
    chain_kernel, design_pencil, focused_polynomial, local_anchor_sources, noisy_components,
    profile_distance, whitening,
)

ROOT = Path(__file__).resolve().parents[1]
LOW, HIGH = 1., 15.
CENTER, SCALE = 8., 7.
SEED = 18061126
N_CASES = 16
NOISE = (0., 2e-5, 2e-4)
ROUNDING = 8e-13
LAYOUTS = ("primary", "confirmation")
REGIONS = {"L": [1, 5, 6], "square": [1, 2, 5, 6]}
ARMS = ("pencil", "anchor", "random", "completion")
GRID = np.linspace(LOW, HIGH, 113)


def host_model(layout):
    rr, cc = np.indices((6, 6))
    compact = np.column_stack((cc.ravel()-2.5, rr.ravel()-2.5))
    active = (6*rr[1:4, 1:5] + cc[1:4, 1:5]).ravel()
    x, y = compact.T
    if layout == "primary":
        cells = .31 * compact
        host = 2.7 + .65*np.cos(.35*x+.61*y) + 1.1*np.exp(-.22*((x+1)**2+(y-1.5)**2))
        theta = 2*np.pi*np.arange(10)/10 + .213
        source_radii = np.full(10, 1.76)
        receiver = np.array([[.50, 1.56]])
    else:
        theta_grid = -.17
        rotation = np.array([[np.cos(theta_grid), -np.sin(theta_grid)],
                             [np.sin(theta_grid), np.cos(theta_grid)]])
        cells = .34 * compact @ rotation.T + np.array([.04,-.02])
        host = 3.1 + .68*np.sin(.57*x-.43*y) + .80*np.exp(-.31*((x-1.)**2+(y+1.4)**2))
        theta = 2*np.pi*np.arange(10)/10 + .397 + .084*np.sin(1.7*np.arange(10))
        source_radii = 1.82 + .09*np.cos(.8*np.arange(10)+.25)
        receiver = np.array([[-.63, 1.57]])
    host[active] = 0.
    sources = source_radii[:,None] * np.column_stack((np.cos(theta), np.sin(theta)))
    model = Model(centers=cells, sources=sources, receivers=receiver, radius=.045, k=2*np.pi)
    ker = ReducedKernel(model, host, active=active)
    gram = j0(model.k * np.linalg.norm(sources[:, None] - sources[None, :], axis=2))
    W = whitening(gram)
    return model, ker, gram, W


def target_reciprocal(H, R, z):
    w = np.linalg.solve((np.eye(len(z)) - z[:, None] * H).T, R) * z
    return R + w @ H


def coefficient_span_diagnostics(H, whitened_sources, target, degree):
    """Inspect the union of incident coefficients of all minimum-degree chains."""
    E = np.eye(len(H))[:, target]
    left = null_space(whitened_sources.conj().T).conj().T
    b0, b1 = E - CENTER * H @ E, -SCALE * H @ E
    basis, _ = chain_kernel(left @ b0, left @ b1, degree)
    P = basis.reshape(degree + 1, len(target), -1)
    coeff = np.zeros((degree+2, len(H), basis.shape[1]), complex)
    for j, block in enumerate(P):
        coeff[j] += b0 @ block
        coeff[j+1] += b1 @ block
    columns = np.concatenate(tuple(coeff), axis=1)
    normalized = columns / np.linalg.norm(columns, axis=0)
    return dict(
        minimum_degree_chain_count=basis.shape[1],
        all_field_coefficient_columns=columns.shape[1],
        all_field_coefficient_rank=int(np.linalg.matrix_rank(columns)),
        all_field_coefficient_singular_values=np.linalg.svd(columns, compute_uv=False),
        column_normalized_field_singular_values=np.linalg.svd(normalized,compute_uv=False),
        center_evaluated_template_rank=int(np.linalg.matrix_rank(P[0])),
        center_evaluated_template_singular_values=np.linalg.svd(P[0],compute_uv=False),
    )


def subspaces(layout, model, ker, qgram, W):
    H, R, bw = ker.K, ker.L.ravel(), ker.q @ W
    rref = target_reciprocal(H,R,np.full(ker.r, CENTER))
    controls = {}
    for r, target in [(2,[1, 5]), (5,[1, 2, 5, 6, 10])]:
        control = design_pencil(H, bw, target, CENTER, SCALE, rref)
        controls[r] = (None if control is None else
                       dict(degree=control["degree"], shots=control["shots"].shape[1],
                            defect=control["relative_synthesis"]))
    region_data, diagnostics = {}, {}
    for region, target in REGIONS.items():
        start = time.process_time()
        pencil = design_pencil(H, bw, target, CENTER, SCALE, rref)
        design_cpu = time.process_time() - start
        if pencil is None:
            raise ValueError(f"No polynomial template for {layout} {region}")
        Qp = pencil["shots"]
        anchors = np.linspace(LOW,HIGH,Qp.shape[1])
        nominal_rows = []
        for aq in anchors:
            zref = np.full(ker.r,CENTER)
            zref[target] = aq
            nominal_rows.append(target_reciprocal(H,R,zref))
        start = time.process_time()
        Qa = local_anchor_sources(H, bw, target, anchors, nominal_rows)
        anchor_cpu = time.process_time() - start
        if Qa.shape != Qp.shape:
            raise ValueError(f"Incomparable measured ranks for {layout} {region}: {Qa.shape}, {Qp.shape}")
        angle_sines = np.linalg.svd(Qa - Qp @ (Qp.conj().T @ Qa), compute_uv=False)
        projector_diff = np.linalg.norm(Qp @ Qp.conj().T - Qa @ Qa.conj().T, ord=2)
        nuisance = np.setdiff1d(np.arange(ker.r),target)
        observability = np.vstack((R[target], H[np.ix_(nuisance,target)]))
        obs_singular = np.linalg.svd(observability, compute_uv=False)
        diagnostics[region] = dict(
            target=target, degree=pencil["degree"], ranks=pencil["ranks"],
            regional_observability_rank=int(np.linalg.matrix_rank(observability)),
            regional_observability_singular_values=obs_singular,
            smallest_zero_degree_chain_singular=pencil["ranks"][0]["singular"][-1],
            relative_synthesis=pencil["relative_synthesis"],
            shots=Qp.shape[1], anchors=anchors,
            principal_angles_deg=np.rad2deg(np.arcsin(np.clip(angle_sines,0,1))),
            projector_distance=projector_diff, design_cpu_s=design_cpu, anchor_cpu_s=anchor_cpu,
            coefficient_span_diagnostics=coefficient_span_diagnostics(H,bw,target,pencil["degree"]),
            template=pencil,
        )
        region_data[region] = (pencil, Qa, anchors)
    return region_data, dict(
        layout=layout, cells=model.centers, sources=model.sources,
        receiver=model.receivers, active=ker.active, host=ker.m0,
        disk_radius=model.radius, gram_eigenvalues=np.linalg.eigvalsh(qgram),
        active_incident_singular_values=np.linalg.svd(bw,compute_uv=False),
        reduced_dimension=ker.r, host_setup_counts=dict(model.counts),
        regions=diagnostics, diagnostics_other_region_sizes=controls,
    )


class PatchResponse:
    def __init__(self, kernel, feeds, region):
        self.response = Response(kernel, feeds)
        self.target = np.array(region)
        self.nuisance = np.setdiff1d(np.arange(kernel.r), self.target)
        self.dimension = 1+len(self.nuisance)

    def material(self, x):
        z = np.empty(self.response.H.shape[0])
        z[self.target] = x[0]
        z[self.nuisance] = x[1:]
        return z

    def predict(self, x, derivative=False):
        z = self.material(x)
        if not derivative:
            return self.response.predict(z)
        J = self.response.predict(z,derivative=True)
        return np.column_stack((np.sum(J[:,self.target],axis=1), J[:,self.nuisance]))


def fit_trf(forward, data, initials):
    start = time.process_time()
    nlu0, nrhs0 = forward.response.nlu, forward.response.nrhs
    scale = max(np.linalg.norm(data),1e-3)
    rows = []
    for x0 in initials:
        sol = least_squares(
            lambda x: pack(forward.predict(x)-data)/scale, x0,
            jac=lambda x: forward.predict(x,derivative=True)/scale,
            bounds=(LOW,HIGH), max_nfev=100, ftol=3e-12, xtol=3e-12, gtol=3e-12,
        )
        rows.append(dict(x=sol.x,residual_l2=np.linalg.norm(sol.fun)*scale,status=sol.status,
                         nfev=sol.nfev,njev=sol.njev))
    best = min(rows,key=lambda row: row["residual_l2"])
    return dict(best=best, starts=rows, cpu_s=time.process_time()-start,
                reduced_lu=forward.response.nlu-nlu0, reduced_rhs=forward.response.nrhs-nrhs0)


def witness(forward, data, fit, epsilon, indep_feed):
    """Attempt a checked feasible pair from a bounded real local profile search."""
    best = fit["best"]["x"]
    good = [row["x"] for row in fit["starts"] if row["residual_l2"]<=epsilon]
    start = time.process_time()
    scale = max(np.linalg.norm(data),1e-3)
    qdir = np.eye(forward.dimension)[0]
    trials = []

    def feasible(x):
        return (epsilon - np.linalg.norm(forward.predict(x)-data))/scale

    def jac(x):
        res = pack(forward.predict(x)-data)
        J = forward.predict(x,derivative=True)
        n = np.linalg.norm(res)
        return -(res / n) @ J / scale if n else np.zeros(len(x))

    if good and epsilon > 2 * ROUNDING:
        for sign in (1,-1):
            warm = min(good,key=lambda x:sign*x[0])
            sol = minimize(
                lambda x: sign*x[0], warm, jac=lambda _:sign*qdir,
                method="SLSQP", bounds=[(LOW,HIGH)]*len(best),
                constraints=[dict(type="ineq",fun=feasible,jac=jac)],
                options=dict(maxiter=120,ftol=2e-10),
            )
            trials.append(dict(sign=sign,x=sol.x,residual_l2=np.linalg.norm(forward.predict(sol.x)-data),status=sol.status))
            for frac in (0.,1e-5,1e-3,.05):
                x = (1-frac)*sol.x + frac*best
                if feasible(x)>=0:
                    good.append(x)
                    break
    extrema = []
    for sign in (1,-1):
        for x in sorted(good,key=lambda x:sign*x[0]):
            z = forward.material(x)
            residual = np.linalg.norm(indep_feed(z)-data)
            if residual <= epsilon:
                extrema.append(dict(q=x[0],z=z,residual_l2=residual))
                break
    return dict(extrema=extrema, trials=trials, cpu_s=time.process_time()-start,
                spread=extrema[1]["q"]-extrema[0]["q"] if len(extrema)==2 else None)


def real_roots(g, truth):
    roots = poly.polyroots(g)
    close_real = [CENTER + SCALE*root.real for root in roots
                  if abs(root.imag) <= 5e-8 and -1-1e-12 <= root.real <= 1+1e-12]
    distinct = [float(q) for q in close_real if abs(q-truth)>1e-5]
    return roots, distinct


def profile_grid(qgrid, d, H, R, fields, target):
    return np.asarray([profile_distance(q,d,H,R,fields,target) for q in qgrid])


def check_sign_cases():
    experiments = [
        ("endpoint root",np.array([-1.,1.]),np.array([[1.,0.]]),0.,[[1.,1.]]),
        ("isolated real root",np.array([.5,1.]),np.array([[1.,0.]]),0.,[[-.5,-.5]]),
        ("irrational isolated roots",np.array([-.5,0.,1.]),np.array([[1.,0.,0.]]),0.,
         [[-np.sqrt(.5),-np.sqrt(.5)],[np.sqrt(.5),np.sqrt(.5)]]),
        ("nonzero epsilon tangent",np.array([-.5+.25j,0.,1.]),np.array([[1.,0.,0.]]),.25,
         [[-np.sqrt(.5),-np.sqrt(.5)],[np.sqrt(.5),np.sqrt(.5)]]),
        ("one connected interval",np.array([0.,1.]),np.array([[1.,0.]]),.3,[[-.3,.3]]),
        ("two intervals",np.array([.3,0.]),np.array([[0.,1.]]),1.,[[-1.,-.3],[.3,1.]]),
        ("whole box",np.array([.3,0.]),np.array([[1.,0.]]),1.,[[-1.,1.]]),
        ("empty",np.array([1.3,0.]),np.array([[1.,0.]]),1.,[]),
    ]
    results = []
    for name,g,c,e,expect in experiments:
        got = noisy_components(g,c,e)
        good = len(got)==len(expect) and np.allclose(got,expect,atol=2e-13,rtol=0.)
        results.append(dict(name=name, actual=got, expected=expect, pass_check=bool(good)))
    return results


def run():
    started = time.process_time()
    layouts = {}
    numerical = {}
    for layout in LAYOUTS:
        model,ker,Q,W = host_model(layout)
        regions, diag = subspaces(layout,model,ker,Q,W)
        for region,d in diag["regions"].items():
            print(layout,region,"degree",d["degree"],"shots",d["shots"],"angles",
                  d["principal_angles_deg"],"source mismatch",d["relative_synthesis"],flush=True)
        layouts[layout] = (model,ker,Q,W,regions)
        numerical[layout] = diag
    energy = min(2., min(
        16 * d["gram_eigenvalues"][0] * min(r["shots"] for r in d["regions"].values())
        for d in numerical.values()))
    print("Common total radiated energy",energy,flush=True)
    out, csv_rows, saved = [], [], {}
    check_cases = check_sign_cases()
    counts = dict(profile_identity_error_max=0., scalar_grid_disagreements=0,
                  scalar_grid_points=0, trf_derivative_relative_errors={})
    for design_id, (layout, layout_data) in enumerate(layouts.items()):
        model,ker,Q,W,regions = layout_data
        H, R = ker.K, ker.L.ravel()
        host_y = ker.y0.ravel()
        clock_start = time.process_time()

        def independent(z):
            return model.forward(ker.m0+ker.U@z,independent=True)["y"].ravel()-host_y

        for reg_id,(region,(design,Qa,anchors)) in enumerate(regions.items()):
            target = REGIONS[region]
            N = np.setdiff1d(np.arange(ker.r),target)
            prefix = f"{layout}_{region}"
            P = design["shots"]
            shot_count = P.shape[1]
            coords = design["coord"]/np.sqrt(energy/shot_count)
            observe_design = {**design,"coord":coords}
            material_rng = np.random.default_rng(SEED+1000*design_id+100*reg_id)
            full_parameters = material_rng.uniform(LOW,HIGH,(N_CASES,1+len(N)))
            start_rng = np.random.default_rng(SEED+1000*design_id+100*reg_id+1)
            start_points = np.vstack((np.full(1+len(N),CENTER),start_rng.uniform(LOW,HIGH,(7,1+len(N)))))
            saved[prefix+"_truth_parameters"] = full_parameters
            saved[prefix+"_initial_parameters"] = start_points
            for case,x_true in enumerate(full_parameters):
                noise_rng = np.random.default_rng(SEED+8000+1000*design_id+100*reg_id+case)
                random_z = noise_rng.standard_normal((model.ns,shot_count))+1j*noise_rng.standard_normal((model.ns,shot_count))
                random_q = np.linalg.qr(random_z)[0]
                Qall = dict(pencil=P,anchor=Qa,random=random_q,completion=np.eye(model.ns))
                feeds = {arm: W @ U * np.sqrt(energy/U.shape[1]) for arm,U in Qall.items()}
                forwards = {arm: PatchResponse(ker,A,target) for arm,A in feeds.items()}
                z = forwards["pencil"].material(x_true)
                clean = independent(z)
                reduced = Response(ker,np.eye(model.ns)).predict(z)
                clean_check = float(np.max(abs(clean-reduced)))
                true_born = (R*z) @ ker.q
                vtrue = target_reciprocal(H,R,z)
                if case == 0:
                    x0=start_points[0]
                    analytic = forwards["pencil"].predict(x0,derivative=True)
                    eye=np.eye(len(x0))
                    fd=np.column_stack([pack(forwards["pencil"].predict(x0+1e-4*unit)-forwards["pencil"].predict(x0-1e-4*unit))/2e-4 for unit in eye])
                    counts["trf_derivative_relative_errors"][prefix] = np.linalg.norm(analytic-fd)/np.linalg.norm(analytic)
                for arm, A in feeds.items():
                    saved[f"{prefix}_{case:02d}_{arm}_feed"] = A
                clean_g = focused_polynomial(clean @ feeds["pencil"],R[target],observe_design,CENTER,SCALE)
                clean_roots, nontruth = real_roots(clean_g,x_true[0])
                t_true=(x_true[0]-CENTER)/SCALE
                true_template = poly.polyval(t_true,design["templates"].T)
                for level, epsilon in enumerate(NOISE):
                    eta = noise_rng.standard_normal(model.ns)+1j*noise_rng.standard_normal(model.ns)
                    radial = noise_rng.uniform()
                    bound = epsilon+ROUNDING
                    level_rows = {}
                    for arm in ARMS:
                        U,A = Qall[arm],feeds[arm]
                        noise = eta @ U
                        noise *= (epsilon*radial)/np.linalg.norm(noise)
                        d = clean @ A + noise
                        saved[f"{prefix}_{case:02d}_{arm}_noise{level}"] = d
                        y = d
                        fit = fit_trf(forwards[arm],y,start_points)
                        profile = profile_grid(GRID,d,H,R,ker.q @ A,target)
                        exact_profile = profile_distance(x_true[0],d,H,R,ker.q@A,target)
                        away = abs(GRID[:,None]-anchors).min(axis=1)>1e-6
                        outside_truth = abs(GRID-x_true[0])>.25
                        off = away & outside_truth
                        certificate = None
                        if arm=="pencil":
                            begin=time.process_time()
                            g = focused_polynomial(d,R[target],observe_design,CENTER,SCALE)
                            components = noisy_components(g, coords, bound)
                            scalar_cpu = time.process_time()-begin
                            values = [[CENTER+SCALE*lo,CENTER+SCALE*hi] for lo,hi in components]
                            hull = [values[0][0],values[-1][1]] if values else None
                            width = hull[1]-hull[0] if hull else None
                            contains_truth = any(lo <= x_true[0] <= hi for lo,hi in values)
                            certificate = dict(
                                components=values, hull=hull, width=width, covers=contains_truth,
                                total_component_length=sum(hi-lo for lo,hi in values),
                                full_box_hull=hull==[LOW,HIGH], clean_roots=clean_roots,
                                clean_nontruth_real_roots=nontruth, scalar_cpu_s=scalar_cpu,
                                true_clean_visibility=abs(vtrue[target] @ true_template),
                                truth_g=poly.polyval(t_true,g),
                                estimate=sum(hull)/2 if hull else None, g=g,
                            )
                            tt=(GRID-CENTER)/SCALE
                            c=poly.polyval(tt,coords.T)
                            vv=abs(poly.polyval(tt,g))
                            ofeeds=np.linalg.norm(c,axis=0)
                            scalar_distance=vv/ofeeds
                            counts["profile_identity_error_max"] = max(counts["profile_identity_error_max"],
                                np.max(abs(scalar_distance[profile[:,1]==1]-profile[profile[:,1]==1,0]),initial=0.))
                            polynomial_test=(vv <= bound*ofeeds)
                            interval_test=np.array([any(l<=q<=u for l,u in values) for q in GRID])
                            counts["scalar_grid_disagreements"]+=np.count_nonzero(polynomial_test != interval_test)
                            counts["scalar_grid_points"]+=len(GRID)
                        w = None
                        if case in (0,1) and epsilon == NOISE[-1]:
                            w=witness(forwards[arm],y,fit,bound,lambda m:independent(m)@A)
                        level_rows[arm] = d
                        row = dict(
                            layout=layout, region=region, case=case, arm=arm, epsilon=epsilon,
                            truth=x_true[0], shots=len(d), interval=certificate,
                            trf=fit, profile_grid=profile, profile_at_truth=exact_profile,
                            profile_off_anchor_rejections=np.mean(profile[off,0]>bound),
                            profile_off_anchor_null_dim_nonzero=np.mean(profile[off,1]>0),
                            noise_l2=np.linalg.norm(noise), noise_radius_with_rounding=bound,
                            total_radiated_energy=np.trace(A.conj().T@Q@A).real,
                            maximum_feed_norm=np.linalg.norm(A,axis=0).max(),
                            bounded_SNR_dB=(20*np.log10(np.linalg.norm(clean@A)/epsilon) if epsilon else None),
                            full_reduced_diagnostic=clean_check,
                            relative_Born_response_disagreement=np.linalg.norm(clean-true_born)/np.linalg.norm(clean),
                            witness=w,
                        )
                        out.append(row)
                        iv = certificate or {}
                        est=iv.get("estimate")
                        csv_rows.append(dict(
                            layout=layout,region=region,case=case,arm=arm,epsilon=epsilon,truth=x_true[0],shots=len(d),
                            energy=row["total_radiated_energy"],max_feed_norm=row["maximum_feed_norm"],
                            interval_width=iv.get("width"),interval_covers=iv.get("covers"),
                            interval_component_count=len(iv["components"]) if iv else None,
                            nontruth_clean_roots=len(iv["clean_nontruth_real_roots"]) if iv else None,
                            interval_error=abs(est-x_true[0]) if est is not None else None,
                            trf_error=abs(fit["best"]["x"][0]-x_true[0]),trf_residual=fit["best"]["residual_l2"],
                            profile_off_anchor_rejections=row["profile_off_anchor_rejections"],
                            profile_off_anchor_null_dim_nonzero=row["profile_off_anchor_null_dim_nonzero"],
                            profile_at_truth=exact_profile[0],
                            witness_spread=w["spread"] if w else None,
                        ))
                    transform=P.conj().T@Qa
                    if region=="L":
                        equivariance_error=np.linalg.norm(level_rows["pencil"]@transform-level_rows["anchor"])
                        maximum=numerical[layout].get("L_anchor_pencil_max_observation_equivariance_error",0.)
                        numerical[layout]["L_anchor_pencil_max_observation_equivariance_error"]=max(maximum,equivariance_error)
                if (case+1)%8==0:
                    print(layout,region,"done",case+1,"of",N_CASES,"CPU",time.process_time()-clock_start,flush=True)
    directory=ROOT/"results"
    with (directory/"pencil_individual.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=csv_rows[0].keys())
        writer.writeheader()
        writer.writerows(csv_rows)
    raw=dict(seed=SEED,noise_norms=NOISE,common_total_energy=energy, q_test_grid=GRID,
             q_box=(LOW,HIGH),numeric_margin=ROUNDING,design=numerical,
             simple_interval_checks=check_cases,numerical_checks=counts,
             observations=out,cpu_seconds=time.process_time()-started)
    (directory/"pencil_raw.json").write_text(json.dumps(raw,default=encode,separators=(",",":"))+"\n")
    np.savez_compressed(directory/"pencil_observations.npz",**saved)
    print("Completed",len(out),"observation rows with checks",counts,"in",time.process_time()-started,"CPU s",flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-only",action="store_true")
    args=parser.parse_args()
    if args.design_only:
        for layout in LAYOUTS:
            print(layout)
            model,ker,Q,W=host_model(layout)
            _,d=subspaces(layout,model,ker,Q,W)
            print(json.dumps({**d,"regions":{
                key:{k:v for k,v in x.items() if k!="template"} for key,x in d["regions"].items()
            }},default=encode))
        print(check_sign_cases())
    else:
        run()
