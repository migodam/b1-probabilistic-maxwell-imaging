"""Outward checks using shared full-resolution stationary responses."""

import argparse
from collections import defaultdict
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh, svd
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh, gcrotmk, gmres
from threadpoolctl import threadpool_limits

from screen import (CAPS, CURVE, HERE, RECEIVERS, SOURCES, frozen_profiles, grid_profile,
                    inputs, load_json, h, save_json, times)

sys.path.insert(0, str(CURVE))
from envelope import FourierEnvelope
from native import (Geometry, Native, Weights, ball, frobenius_squared, gershgorin,
                    high, low, mid, mult_rows, star)


def digest_profile(profile, cap):
    scale, wr, op = profile["mu"], profile["sector"], profile["optical"]
    if profile["normalized"]:
        # Treat the binary64 products as the exact native weights. The normalized
        # design in the prior proof agrees with these ratios to binary roundoff.
        wr, op = [x*scale for x in wr], [x*scale for x in op]
    return Weights(cap,scale,tuple(wr),tuple(op))


class NativeState(Native):
    def __init__(self, geo, weights):
        self.isotropic = weights.is_global
        super().__init__(geo,weights)
        if self.isotropic:
            ratio = weights.optical[0]/weights.scale
            symbol = -weights.cap*geo.spec.real + ratio*geo.spec.imag
            self.numeric_symbol = symbol.astype(complex)
            u, t = arb(weights.cap), arb(weights.optical[0])/arb(weights.scale)
            self.exact_symbol = [[acb(-u*v.real+t*v.imag) for v in row] for row in geo.spec_ball]

    def numeric_action(self, x):
        if not self.isotropic:
            return super().numeric_action(x)
        g = self.g
        shape = x.shape
        x = x.reshape(g.n,-1)
        spec = g.spec
        try:
            g.spec = self.numeric_symbol
            z = g.numerical_green(x)
        finally:
            g.spec = spec
        return (x + g.sensor_star@(g.sensor@x)/self.w.scale + z).reshape(shape)

    def exact_action(self, x):
        if not self.isotropic:
            return super().exact_action(x)
        g = self.g
        spec = g.spec_ball
        try:
            g.spec_ball = self.exact_symbol
            z = g.exact_green(x)
        finally:
            g.spec_ball = spec
        return x + star(g.sensor_ball)*(g.sensor_ball*x)/arb(self.w.scale) + z

    def modes(self, expected, *, equilibrate=False):
        start,cpu = perf_counter(),process_time()
        rng = np.random.default_rng(54273)
        initial = rng.standard_normal(self.g.n)+1j*rng.standard_normal(self.g.n)
        k = max(6,expected+4)
        self.eigen_operator = self.operator
        if equilibrate:
            _, sig, vstar = svd(self.g.sensor,full_matrices=False,check_finite=False)
            attenuation = 1-np.sqrt(self.w.scale/(self.w.scale+sig**2))

            def half(x):
                a = x.reshape(self.g.n,-1)
                b = a-vstar.conj().T@(attenuation[:,None]*(vstar@a))
                return b.reshape(x.shape)

            def scaled(x):
                return half(self.numeric_action(half(x)))

            self.eigen_operator=LinearOperator((self.g.n,self.g.n),matvec=scaled,
                                                matmat=scaled,dtype=np.complex128)
        while True:
            try:
                val,vec = eigsh(self.eigen_operator,k=k,ncv=max(36,3*k+4) if equilibrate else max(30,2*k+4),
                                 which="SA",tol=5e-11,maxiter=1400,v0=initial)
            except ArpackNoConvergence:
                if equilibrate:
                    raise
                result=self.modes(expected,equilibrate=True)
                result.update(retried_after_unscaled_arpack=True,**times(start,cpu))
                return result
            if val.max() >= 0:
                break
            k += 6
        order=np.argsort(val)
        self.eigs=val[order]
        v=vec[:,order]
        self.eig_basis=half(v) if equilibrate else v
        self.nu=int(np.count_nonzero(self.eigs<0.))
        defect=self.eigen_operator@v-v*self.eigs
        return dict(numerical_smallest_scaled_K=self.eigs.tolist(),
                    numerical_eigen_defect=float(np.linalg.norm(defect)),
                    numerical_space="receiver-equilibrated congruence" if equilibrate else "native K/mu",
                    **times(start,cpu))


def check_state(native, envelope, expected, *, equilibrate=False):
    start,cpu = perf_counter(),process_time()
    cost0=native.g.cost.copy()
    numerical = native.modes(expected,equilibrate=equilibrate)
    nu=native.nu
    neg = native.ensure_negative(nu,return_gap=True)
    neg_upper, neg_gap = neg if nu else (None,None)
    rec=dict(numeric=numerical,native_nu_at_least=nu,negative_congruence_upper=neg_upper,
             negative_gap_lower=neg_gap)
    try:
        proof=envelope.verify(native.w,nu)
        rec.update(exact_nu=nu,positive_gap_lower=proof["positive_gap_lower"],
                   stable_state_proof=proof)
    except ArithmeticError as err:
        rec["positive_gap_error"]=str(err)
    return dict(**rec,operations=dict(native.g.cost-cost0),**times(start,cpu))


def bin_rhs(native, incident):
    g=native.g
    return np.column_stack((g.sensor_star/native.w.scale,native.z[:,None]*incident))


def outward_rhs(native, incident):
    g=native.g
    _,z=native.exact_scaled_coefficients()
    port=star(g.sensor_ball)/arb(native.w.scale)
    launch=mult_rows(ball(incident),z)
    return acb_mat(g.n,port.ncols()+launch.ncols(),[
        a[i,j] for i in range(g.n) for a in (port,launch) for j in range(a.ncols())])


def reusable_stationarity(native, source_bank):
    """Apply the native operator once to each common source or receiver response."""
    t,c=perf_counter(),process_time()
    g=native.g
    before=g.cost.copy()
    incident=np.column_stack(list(source_bank.values()))
    rhs=bin_rhs(native,incident)
    pre=LinearOperator((g.n,g.n),matvec=native.precondition,dtype=np.complex128)
    prop = np.empty_like(rhs)
    recycled, iterations, solver_fallback = [], [], []
    for j, b in enumerate(rhs.T):
        steps=[]
        x,info=gcrotmk(native.operator,b,M=pre,CU=recycled,m=25,k=12,
                        maxiter=16,rtol=8.e-12,atol=0.,callback=lambda q:steps.append(1))
        if info:
            x,info=gmres(native.operator,b,x0=x,M=pre,restart=100,maxiter=12,
                         rtol=8.e-12,atol=0.)
            solver_fallback.append(j)
        # An unconverged numerical trial may still give an interval result via its
        # actual defect. Retain its charge and distinguish it in the ledger.
        residual=float(np.linalg.norm(native.numeric_action(x)-b)/np.linalg.norm(b))
        if not np.isfinite(residual):
            raise ArithmeticError(f"Stationary trial has an invalid residual at {j}.")
        prop[:,j]=x
        iterations.append(dict(gcrotmk_cycles=len(steps),solver_status=int(info),relative_residual=residual))
    numeric=times(t,c)
    t1,c1=perf_counter(),process_time()
    xb=ball(prop)
    L=outward_rhs(native,incident)
    AX=native.exact_action(xb)
    defect=AX-L
    stat=star(xb)*AX-star(xb)*L-star(L)*xb
    gram=star(defect)*defect
    idx={}
    start=len(g.sensor)
    for s in source_bank:
        idx[s]=list(range(len(g.sensor)))+list(range(start,start+s))
        start += s
    return stat,gram,idx,dict(common_columns=rhs.shape[1],common_receiver_columns=len(g.sensor),
        common_incident_columns=incident.shape[1],solver_iterations=iterations,
        solver_fallback_columns=solver_fallback,numeric_stationarity=numeric,
        interval_stationarity=times(t1,c1),frobenius_defect_upper=high(frobenius_squared(defect).sqrt()),
        operations=dict(g.cost-before),**times(t,c))


def submatrix(a, indices):
    return acb_mat(len(indices),len(indices),[a[i,j] for i in indices for j in indices])


def analyze_row(y, state, mu, stat, gram, indices, radius):
    start,cpu=perf_counter(),process_time()
    p,s=y.shape
    nu=state["exact_nu"]
    yy=ball(y)
    eye=acb_mat(s,s)
    for i in range(s):
        eye[i,i]=1
    all_y=acb_mat(p+s,s,yy.entries()+eye.entries())
    ya=star(all_y)
    trial_form=star(yy)*yy+arb(mu)*(ya*submatrix(stat,indices)*all_y)
    error=ya*submatrix(gram,indices)*all_y
    eta=arb(mu)*sum((error[i,i].real for i in range(s)),arb(0)).upper()
    scale_lower=eta/arb(state["positive_gap_lower"])
    d0=min(state["positive_gap_lower"],state["negative_gap_lower"]) if nu else state["positive_gap_lower"]
    scale_upper=eta/arb(d0)
    eig,v=eigh(h(mid(trial_form)),driver="evd",check_finite=False)
    eig,v=eig[::-1],v[:,::-1]
    v=ball(v)
    vg=star(v)*v
    vf=star(v)*trial_form*v
    gmin=gershgorin(vg,upper=False)
    gmax=gershgorin(vg,upper=True)
    if not gmin > 0:
        raise ArithmeticError("Source eigenbasis Gram has unresolved positive lower bound.")
    lower=[]
    for k in range(nu+1,s+1):
        idx=list(range(k))
        left=gershgorin(submatrix(vf,idx),upper=False)
        if not left>0:
            break
        left=left/gershgorin(submatrix(vg,idx),upper=True)-scale_lower
        if not left>0:
            break
        lower.append(left.lower())
    for i,l in enumerate(eig):
        vf[i,i]-=arb(float(l))
    congruence_shift=frobenius_squared(vf).sqrt().upper()
    # The first correction controls the congruence of H0. The second is
    # already in the original source coordinates and follows via Weyl.
    upper=[]
    signed_upper=[]
    signed_lower=[]
    for value in eig:
        right=(arb(float(value))+congruence_shift).upper()
        left=(arb(float(value))-congruence_shift).lower()
        ehi=(right/(gmin if right>0 else gmax)+scale_upper).upper()
        elo=(left/(gmax if left>0 else gmin)-scale_upper).lower()
        signed_upper.append(ehi)
        signed_lower.append(elo)
        upper.append(max(arb(0),ehi))
    upper=upper[nu:]
    lower_dist=max(0.,low(sum(lower,arb(0)).sqrt()))
    upper_dist=high(sum(upper,arb(0)).sqrt()) if any(v>0 for v in upper) else 0.
    count_lower=len(lower)
    count_upper=sum(bool(val>0) for val in upper)
    margin=low(arb(lower_dist)-arb(radius))
    return dict(status="native_source_interval",negative_state_count=nu,
                omitted_leading_source_eigenvalues=nu,
                approximate_source_eigenvalues=eig.tolist(),
                source_eigenvalues_lower=[low(x) for x in signed_lower],
                source_eigenvalues_upper=[high(x) for x in signed_upper],
                first_retained_source_eigenvalue=float(eig[nu]),
                first_retained_source_eigenvalue_upper=high(signed_upper[nu]),
                positive_retained_eigenvalue_count_lower=count_lower,
                positive_retained_eigenvalue_count_upper=count_upper,
                lower_retained_eigenvalues=[low(x) for x in lower],
                upper_positive_part_retained_eigenvalues=[high(x) if x>0 else 0. for x in upper],
                distance_lower=lower_dist,distance_upper=upper_dist,
                margin_lower=margin,strictly_excludes=margin>0.,
                lower_defect_correction_upper=high(scale_lower),
                absolute_gap_defect_correction_upper=high(scale_upper),**times(start,cpu))


def frontier(rows, checked, attempted, verified):
    choices=defaultdict(set)
    fixed={p["name"]:p for p in frozen_profiles()}
    for row in rows:
        ident=row["id"]
        if ident in verified or attempted[ident]>=3:
            continue
        mark=row["padded_radius"]*.985
        for f in row["fixed_profiles"]:
            if f["name"] not in checked and f["distance"]>mark:
                choices[f["name"]].add(ident)
        for f in row["scalar_grid_best"]:
            if (f"scalar_grid_{f['index']}" not in checked and f["distance"]>mark):
                choices[f"scalar_grid_{f['index']}"].add(ident)
    if not choices:
        return None
    name=max(choices,key=lambda k:(len(choices[k]),k in fixed,k))
    it=choices[name]
    return fixed[name] if name in fixed else grid_profile(int(name.split("_")[-1])),it


def projected_expected(profile, priors):
    if not priors:
        return None
    if profile["name"].startswith("scalar_grid"):
        index=int(profile["name"].split("_")[-1])
        vals=[d["nu"] for r in priors for d in r["scalar_grid_best"] if d["index"]==index]
    else:
        vals=[d["nu"] for d in priors[0]["fixed_profiles"] if d["name"]==profile["name"]]
    return vals[0] if vals else None


def check_profile(profile,u,geos,envelopes,sensors,source_bank,observations,selection,prior, *,
                  equilibrate=False):
    start,cpu=perf_counter(),process_time()
    st = dict(name=profile["name"],profile=profile,cap=u,observations=[])
    weights=digest_profile(profile,u)
    side=profile["partition"]
    if side not in geos:
        # For the state/distance common RHS method the geometry only needs one B.
        # The source bank is assembled independently for each selected port count.
        geos[side]=Geometry(geos[2].spec,sensors,source_bank[SOURCES[0]],block_side=side)
        envelopes[side]=FourierEnvelope(geos[side])
    kw=dict(geos[side].cost)
    native=NativeState(geos[side],weights)
    state=check_state(native,envelopes[side],prior if prior is not None else int(38*u-20),
                      equilibrate=equilibrate)
    st["state"]=state
    print("state", len(sensors), u, profile["name"], state.get("exact_nu"),
          state.get("positive_gap_error"), round(state["wall_seconds"],2), flush=True)
    subset={s:v for s,v in observations.items() if selection is None or any(x["id"] in selection for x,y in v)}
    can_use=[s for s in subset if s>state["native_nu_at_least"] and "exact_nu" in state]
    if can_use:
        try:
            stat,gram,indices,sol=reusable_stationarity(native,{s:source_bank[s] for s in can_use})
            st["shared_source_response"]=sol
        except ArithmeticError as err:
            st["shared_source_response_error"]=str(err)
    for s,ops in subset.items():
        for x,y in ops:
            r=dict(id=x["id"],object_id=x["object_id"],sources=s,receivers=x["receivers"],
                   noise_pct=x["noise_pct"],radius=x["padded_radius"])
            if s<=state["native_nu_at_least"]:
                rr=dict(status="native_inertia_exhaustion",negative_state_count_at_least=state["native_nu_at_least"],
                        distance_lower=0.,distance_upper=0.,margin_lower=-r["radius"],strictly_excludes=False)
                if "exact_nu" in state:
                    rr["negative_state_count"]=state["exact_nu"]
            elif s in can_use and "shared_source_response" in st:
                rr=analyze_row(y,state,weights.scale,stat,gram,indices[s],r["radius"])
            else:
                rr=dict(status="native_source_unresolved",strictly_excludes=False)
            st["observations"].append(dict(**r,**rr))
    st["total_operations"]={
        k:value-kw.get(k,0) for k,value in geos[side].cost.items()}
    st.update(times(start,cpu))
    return st


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receivers",type=int,choices=RECEIVERS,required=True)
    parser.add_argument("--fixed-only",action="store_true")
    parser.add_argument("--high-grid-controls",action="store_true")
    parser.add_argument("--equilibrate",action="store_true")
    parser.add_argument("--cap",type=float,choices=CAPS,action="append")
    args=parser.parse_args()
    ctx.prec,ctx.threads=112,2
    start,cpu=perf_counter(),process_time()
    spec,sensors,sources,observations=inputs()
    p=args.receivers
    prior_file=HERE/"screen.json"
    priors=load_json(prior_file)["records"] if prior_file.is_file() else []
    geos={2:Geometry(spec,sensors[p],sources[SOURCES[0]])}
    envelopes={2:FourierEnvelope(geos[2])}
    results=dict(scope="Exact represented binary64 discrete coefficients, full 9216-cell FFT",
                receiver_count=p,precision_bits=ctx.prec,arb_threads=ctx.threads,
                profiles=[],fixed_only=args.fixed_only,geometry_setup=times(start,cpu))
    suffix="_recovery" if args.equilibrate else ""
    out=HERE/(f"high_grid_p{p}{suffix}.json" if args.high_grid_controls else f"native_p{p}{suffix}.json")
    if args.high_grid_controls:
        for u in args.cap or CAPS[1:]:
            records=[r for r in priors if r["receivers"]==p and r["cap"]==u and
                     r["sources"]==32 and r["scalar_grid_best"]]
            if not records:
                continue
            best=max(records,key=lambda r:(r["scalar_grid_best"][0]["distance"]/r["padded_radius"],
                                            r["scalar_grid_best"][0]["first_retained"]/r["padded_radius"]**2))
            grid=best["scalar_grid_best"][0]
            profile=grid_profile(grid["index"])
            try:
                st=check_profile(profile,u,geos,envelopes,sensors[p],sources,
                                 observations[p],{best["id"]},grid["nu"],
                                 equilibrate=args.equilibrate)
            except (ArithmeticError,ValueError,RuntimeError) as err:
                st=dict(name=profile["name"],profile=profile,cap=u,observations=[],
                        unresolved_error=str(err))
            st.update(selection_id=best["id"], projected_winner=grid)
            results["profiles"].append(st)
            results.update(times(start,cpu))
            save_json(out,results)
            print("high-cap control",p,u,profile["name"],st.get("state",{}).get("exact_nu"),
                  st.get("unresolved_error"),flush=True)
        return
    for u in args.cap or CAPS:
        pri=[r for r in priors if r["receivers"]==p and r["cap"]==u]
        attempted=defaultdict(int)
        verified=set()
        checked=set()
        profile=next(frozen_profiles())
        selection=None
        while True:
            name=profile["name"]
            st=check_profile(profile,u,geos,envelopes,sensors[p],sources,observations[p],selection,
                             projected_expected(profile,pri))
            for rr in st["observations"]:
                if rr["strictly_excludes"]:
                    verified.add(rr["id"])
                if selection and rr["id"] in selection:
                    attempted[rr["id"]]+=1
            checked.add(name)
            results["profiles"].append(st)
            results.update(times(start,cpu))
            save_json(out,results)
            print(p,u,name,st["state"].get("exact_nu"),st["state"]["native_nu_at_least"],
                  len(verified),round(st["wall_seconds"],2),st["state"].get("positive_gap_error",""),
                  flush=True)
            next_one=frontier(pri,checked,attempted,verified) if len(checked)<12 and not args.fixed_only else None
            if next_one is None:
                break
            profile,selection=next_one


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
