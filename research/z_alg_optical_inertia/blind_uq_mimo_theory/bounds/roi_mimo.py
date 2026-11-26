"""Blind physical endpoint fits on new observations for frozen stationary ROI 2."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve, solve
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/"bounds"/"roi_mimo_fits"
PUB = ROOT/"ports"/"public"
sys.path[:0] = [str(ROOT/p) for p in ("imaging","diagnostics","ports")]
from shared import Model, write_json, write_arrays
from cie import CIE, CurrentState
from verified_wave import VerifiedWave, up, fro, rd
from model import DiskGreen, network_terms, operator
from sweep import change_network, normalized_reactance_error

ETAS = (0,5,25,60)
COUNTS = (6,8,12,16)
U = 4.5


class ActiveModel(Model):
    """Exact elimination of zero material cells, independently audited by the Green action."""

    def solve(self, x, b, adjoint=False, rtol=None):
        if self.green.dense is None:
            return super().solve(x,b,adjoint,rtol)
        mode="physical_adjoint" if adjoint else "physical_forward"
        self.cost[mode+"_rhs"]+=b.shape[1]
        state=self.green.dense
        active=np.flatnonzero(x)
        inactive=np.flatnonzero(x==0)
        if not len(active):
            return b.copy()
        if not adjoint or getattr(self,"fact_material",None) is None or not np.array_equal(self.fact_material,x):
            system=np.eye(len(active),dtype=complex)-x[active,None]*state[np.ix_(active,active)]
            self.factor=lu_factor(system,check_finite=False)
            self.fact_material=x.copy()
            self.cost["physical_active_lu"]+=1
            self.cost["active_cells_factorized"]+=len(active)
        z=b.copy()
        local=b[active].copy()
        if len(inactive) and not adjoint and np.any(b[inactive]):
            local+=x[active,None]*(state[np.ix_(active,inactive)]@b[inactive])
        z[active]=lu_solve(self.factor,local,trans=2 if adjoint else 0,check_finite=False)
        if len(inactive) and adjoint:
            z[inactive]+=state[np.ix_(active,inactive)].conj().T@(x[active,None]*z[active])
        defect=(z-self.green.adjoint(x[:,None]*z)-b if adjoint else
                z-x[:,None]*self.green.apply(z)-b)
        residual=float(np.max(np.linalg.norm(defect,axis=0)/np.maximum(np.linalg.norm(b,axis=0),1e-300)))
        self.cost["max_"+mode+"_relative_defect"]=max(self.cost["max_"+mode+"_relative_defect"],residual)
        if not np.isfinite(residual) or residual>1e-8:
            raise ArithmeticError(f"Active {mode} original equation defect {residual}")
        return z


def public_operator(side,q):
    name=f"n{side}_f400_a{q:02}"+("" if q==8 else "_c048")
    with np.load(PUB/"operators"/f"{name}.npz") as values:
        return name,{k:values[k] for k in values.files}


def view_tag(q,e):
    return f"f400_a{q:02}_c048_e{e:03}"


def observation(case,q,e):
    with np.load(PUB/"cases"/f"{case}_mimo.npz") as arr:
        key=view_tag(q,e)
        return arr[key+"_Y"],float(arr[key+"_radius"])


def seed_files(case):
    base=ROOT/"imaging"/"results"/case/"stationary"
    with np.load(base/"baseline.npz") as arr:
        mask=arr["rois"][1].ravel()
        image=arr["chi"].ravel()
    seeds={"stationary_image":image}
    for path,key,tag in ((base/"b1_controls.npz","physical_chi","unrestricted_public"),
                         (ROOT/"results"/"gn_profiles"/case/"roi_02.npz","material","prior_zero_gn"),
                         (ROOT/"results"/"stable_profiles"/case/"roi_02.npz","material","prior_low"),
                         (ROOT/"results"/"stable_profiles"/case/"roi_02_pin.npz","material","prior_high_pin")):
        if not path.is_file():
            continue
        with np.load(path) as arr:
            if "roi" in arr and not np.array_equal(arr["roi"].ravel(),mask):
                continue
            if key in arr:
                seeds[tag]=arr[key].ravel()
    return mask,seeds


def box_start(model,seeds,lo,hi):
    starts=[]
    for tag,x in seeds.items():
        candidate=np.clip(x,lo,hi)
        if any(np.array_equal(candidate,a[2]) for a in starts):
            continue
        try:
            ratio=float(model.forward(candidate)["residual_ratio"])
        except (ArithmeticError,ValueError):
            continue
        starts.append((ratio,tag,candidate))
    if not starts:
        raise ArithmeticError("No qualified initial state solve for this material box")
    best=min(starts,key=lambda x:x[0])
    return best[2],best[1],{tag:ratio for ratio,tag,_ in starts}


def gn(model,seed,lo,hi,steps,target):
    x=seed.copy()
    damping=.002
    history=[]
    reason="iterations_exhausted"
    for iteration in range(steps):
        state=model.forward(x)
        old=float(state["residual_ratio"])
        if old<target:
            reason="qualified_numerical_target"
            break
        receiver=model.solve(x,model.SH,adjoint=True).conj().T
        jac=(receiver[:,None,:]*state["field"].T[None,:,:]).reshape(model.Y.size,model.n)
        a=np.vstack((jac.real,jac.imag))
        difference=(model.Y-state["pred"]).ravel()
        error=np.r_[difference.real,difference.imag]
        gradient=a.T@error
        movable=hi>lo
        near_l=x<=lo+1e-9
        near_h=x>=hi-1e-9
        active=movable&(~near_l | (gradient>0))&(~near_h | (gradient<0))
        if not np.any(active):
            reason="no_free_descent_direction"
            break
        scale=max(float(np.linalg.norm(a[:,active])**2/min(a.shape[0],active.sum())),1e-19)
        for j in range(4):
            ids=np.flatnonzero(active)
            if not len(ids):
                break
            local=a[:,ids]
            normal=local@local.T
            normal.flat[::len(error)+1]+=damping*scale
            step=local.T@solve(normal,error,assume_a="pos",check_finite=False)
            bad=(near_l[ids]&(step<0))|(near_h[ids]&(step>0))
            if not np.any(bad) or j==3:
                break
            active[ids[bad]]=False
        if not len(ids):
            reason="active_set_exhausted"
            break
        delta=np.zeros_like(x)
        delta[ids]=step*min(1.,1.2/max(abs(step).max(),1e-16))
        accepted=None
        tried=0
        for fraction in (1.,.5,.25,.125,.0625):
            proposal=np.clip(x+fraction*delta,lo,hi)
            if np.array_equal(x,proposal):
                continue
            tried+=1
            try:
                value=float(model.forward(proposal)["residual_ratio"])
            except (ArithmeticError,ValueError):
                continue
            if value<old*(1-1e-7):
                accepted=(value,proposal,fraction)
                break
        if accepted is None:
            damping=min(1e3,damping*8.)
        else:
            x=accepted[1]
            damping=max(1e-6,damping*.6) if accepted[2]>=.5 else min(.5,damping*1.5)
        history.append(dict(iteration=iteration+1,residual_before=old,
                            residual_accepted=accepted[0] if accepted else None,trial_evaluations=tried,
                            damping=damping,free=int(len(ids))))
        if damping>=1e3:
            reason="damping_restart_limit"
            break
    st=model.forward(x)
    return x,dict(residual_ratio=float(st["residual_ratio"]),state_relative_defect=float(st["state_relative_defect"]),
                  stop=reason,history=history)


def verify_wave(model,checker,x,accept_numeric=1.):
    elapsed=time.perf_counter()
    cpu=time.process_time()
    state=model.forward(x)
    ratio=float(state["residual_ratio"])
    row=dict(peak=None,verified_inside_noise=False,numerical_residual_ratio=ratio,
             original_relative_state_defect=float(state["state_relative_defect"]))
    if ratio<accept_numeric:
        try:
            cert,proof=checker.inverse(x)
            cur=state["j"]
            rhs=x[:,None]*checker.B
            error=checker.solve_error(x,cur,rhs,up(rd(1)*max(x)*fro(checker.B)),cert)
            pred=checker.S@cur
            residual=up(fro(pred-model.Y)+checker.snorm*error+
                        rd(model.n)*checker.snorm*fro(cur)+rd(2)*(fro(pred)+fro(model.Y)))
            row.update(verified_inside_noise=bool(residual<model.radius),hard_noise_radius=model.radius,
                       residual_bound_upper=residual,residual_ratio_upper=up(residual/model.radius),
                       receiver_current_correction_upper=error,inverse=proof)
        except ArithmeticError as exc:
            row["verification_error"]=str(exc)
    else:
        row["verification_error"]="no numerical crossing"
    return row|dict(wall_s=time.perf_counter()-elapsed,cpu_s=time.process_time()-cpu)


def try_box(model,checker,seeds,lo,hi,mask,steps,target):
    x,seed_tag,seed_rs=box_start(model,seeds,lo,hi)
    winner,replay=gn(model,x,lo,hi,steps,target)
    proof=verify_wave(model,checker,winner)
    proof["peak"]=float(winner[mask].max())
    return winner,dict(seed=seed_tag,initial_ratios=seed_rs,fit=replay,proof=proof)


def pick_pins(m,x,mask):
    state=m.forward(x)
    receiver=m.solve(x,m.SH,adjoint=True).conj().T
    cells=np.flatnonzero(mask)
    impact=np.linalg.norm(receiver[:,cells],axis=0)*np.linalg.norm(state["field"][cells],axis=1)
    increment=np.maximum(U-x[cells],0.)
    chosen=[int(cells[np.argmax(x[cells])])]
    for j in np.argsort(impact*increment,kind="stable"):
        cell=int(cells[j])
        if cell not in chosen:
            chosen.append(cell)
        if len(chosen)>=3:
            break
    return chosen


def nominal_one(case,q,op_id,op,checker,limit,matched_warm=False,followup=False,released=False):
    clock,cpu=time.perf_counter(),time.process_time()
    mask,seeds=seed_files(case)
    y,r=observation(case,q,0)
    solve_model=Model if released else ActiveModel
    m=solve_model(op,y,r,op_id,int(op["grid_shape"][0]))
    lo,hi=np.zeros(m.n),np.full(m.n,U)
    preferred={k:v for k,v in seeds.items() if k in ("stationary_image","unrestricted_public")}
    rootseed=preferred or seeds
    if followup:
        for smaller in (6,8):
            path=OUT/f"n32_q{smaller:02}"/f"{case}.npz"
            if path.exists():
                with np.load(path) as earlier:
                    if not np.array_equal(earlier["roi"],mask):
                        raise ArithmeticError("Mismatched ROI in blind follow-up seeds")
                    for key in ("high","low"):
                        if key in earlier.files:
                            seed=f"q{smaller:02}_{key}"
                            rootseed[seed]=earlier[key]
                            seeds[seed]=earlier[key]
    warm=None
    if matched_warm:
        warmclock=time.perf_counter(),time.process_time()
        c=CIE(m)
        st=CurrentState(seeds["stationary_image"],np.empty((0,q),complex),np.empty((0,2),int),0)
        warmed=[]
        for stage in range(3):
            st,rec,audits=c.advance(st,stage,40)
            warmed.append((float(m.forward(st.x)["residual_ratio"]),st.x.copy()))
        selected=min(warmed,key=lambda a:a[0])
        rootseed["matched_current_seed"]=selected[1]
        seeds["matched_current_seed"]=selected[1]
        warm=dict(projected_current_updates=120,ratios=[r for r,_ in warmed],
                  wall_s=time.perf_counter()-warmclock[0],cpu_s=time.process_time()-warmclock[1])
    point,free_report=try_box(m,checker,rootseed,lo,hi,mask,limit if followup else min(20,limit),.76)
    seeds_low={k:v for k,v in seeds.items() if k in ("stationary_image","unrestricted_public","prior_zero_gn","matched_current_seed")}
    if "prior_zero_gn" not in seeds_low and "prior_low" in seeds:
        seeds_low["prior_low"]=seeds["prior_low"]
    seeds_low["same_layout_fit"]=point
    if followup:
        for key in ("q06_low","q08_low"):
            if key in seeds:
                seeds_low[key]=seeds[key]
    attempts_low=[]
    late_seed=None
    best_attempt=np.inf
    low=None
    low_proof=None
    for cap in (0.,.25,.75,1.5,2.75):
        box_hi=hi.copy()
        box_hi[mask]=cap
        niter=limit if cap==0 else min(36 if followup else 18,limit)
        candidate,rec=try_box(m,checker,seeds_low,lo,box_hi,mask,niter,.78)
        attempts_low.append(dict(cap=cap,**rec))
        if rec["fit"]["residual_ratio"]<best_attempt:
            best_attempt=rec["fit"]["residual_ratio"]
            late_seed=candidate
        if rec["proof"]["verified_inside_noise"]:
            low,low_proof=candidate,rec["proof"]
            break
        seeds_low["previous_cap"]=candidate
    if low is None and free_report["proof"]["verified_inside_noise"]:
        low,low_proof=point,free_report["proof"]
    if matched_warm and late_seed is not None and not free_report["proof"]["verified_inside_noise"]:
        if best_attempt<free_report["fit"]["residual_ratio"]:
            revised,revised_report=try_box(m,checker,{"best_low_box":late_seed},lo,hi,mask,
                                           min(90 if followup else 26,limit),.78)
            if revised_report["fit"]["residual_ratio"]<free_report["fit"]["residual_ratio"]:
                revised_report["original_free_seed_residual"]=free_report["fit"]["residual_ratio"]
                point,free_report=revised,revised_report
                if low is None and free_report["proof"]["verified_inside_noise"]:
                    low,low_proof=point,free_report["proof"]
    hi_seeds={"same_layout_fit":point,"stationary_image":seeds["stationary_image"]}
    if "prior_high_pin" in seeds:
        hi_seeds["prior_high_pin"]=seeds["prior_high_pin"]
    elif "unrestricted_public" in seeds:
        hi_seeds["unrestricted_public"]=seeds["unrestricted_public"]
    attempts_high=[]
    high=None
    high_proof=None
    # A passing previously evaluated map at the prior cap needs no optimizer.
    if free_report["proof"]["verified_inside_noise"] and free_report["proof"]["peak"]==U:
        high,high_proof=point,free_report["proof"]
    else:
        for j,cell in enumerate(pick_pins(m,point,mask)):
            box_lo=lo.copy()
            box_lo[cell]=U
            niter=limit if j==0 else min(20,limit)
            candidate,rec=try_box(m,checker,hi_seeds,box_lo,hi,mask,niter,.78)
            attempts_high.append(dict(pixel=cell,**rec))
            if rec["proof"]["verified_inside_noise"]:
                high,high_proof=candidate,rec["proof"]
                break
    if high is None and free_report["proof"]["verified_inside_noise"]:
        high,high_proof=point,free_report["proof"]
    kept=[(v,p) for v,p in ((low,low_proof),(high,high_proof)) if p is not None]
    if kept:
        low,low_proof=min(kept,key=lambda a:a[1]["peak"])
        high,high_proof=max(kept,key=lambda a:a[1]["peak"])
    return dict(case=case,q=q,training_observation=view_tag(q,0),array_operator=op_id,roi="roi_02",
                region_cells=int(mask.sum()),noise_radius=r,
                incident_available_W=float(op["available_power"]),
                available_per_source_W=float(op["available_power"])/q,
                certified_outer=[0.,U],nominal_inner_peak_range=None if not kept else [low_proof["peak"],high_proof["peak"]],
                nominal_inner_width=None if not kept else high_proof["peak"]-low_proof["peak"],
                zero_and_ceiling_compatible=bool(kept and low_proof["peak"]==0 and high_proof["peak"]==U),
                nominal_endpoint_witnesses={"low":low_proof,"high":high_proof},
                free_arm=free_report,matched_current_warm=warm,attempts_low=attempts_low,attempts_high=attempts_high,
                endpoint_hull_only=True,cost=dict(m.cost),wall_s=time.perf_counter()-clock,
                cpu_s=time.process_time()-cpu),mask,low,high


def op2(disk,q,e,bare):
    net=network_terms(disk,q,48.)
    if e:
        net=change_network(net,e/1000*normalized_reactance_error(net["positions"]))
    return operator(disk,net,bare)


def calibration_replay(case,q,op_id,builders,maps):
    out=[]
    for e,(op,checker,delta_norm) in builders.items():
        y,r=observation(case,q,e)
        clock,cpu=time.perf_counter(),time.process_time()
        normalized=Model(op,y,r,op_id,int(op["grid_shape"][0]))
        sides={}
        for label,(x,peak) in maps.items():
            if x is None:
                continue
            proof=verify_wave(normalized,checker,x)
            proof["peak"]=peak
            sides[label]=proof
        hull=sorted(p["peak"] for p in sides.values() if p["verified_inside_noise"])
        out.append(dict(calibration_level=e/1000,calibration_error_norm_ohms=delta_norm,
                        observation=view_tag(q,e),measurement_noise_radius=r,
                        same_two_nominal_materials_only=True,matched_reconstructed_operator=True,
                        certified_outer=[0.,U],inner_peak_range=[hull[0],hull[-1]] if hull else None,
                        inner_width=hull[-1]-hull[0] if hull else None,
                        endpoints=sides,cost=dict(normalized.cost),wall_s=time.perf_counter()-clock,
                        cpu_s=time.process_time()-cpu))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--q",type=int,choices=COUNTS,action="append",required=True)
    ap.add_argument("--case",action="append")
    ap.add_argument("--limit",type=int,default=30)
    ap.add_argument("--replay-only",action="store_true")
    ap.add_argument("--skip-existing",action="store_true")
    ap.add_argument("--matched-warm",action="store_true")
    ap.add_argument("--followup",action="store_true")
    ap.add_argument("--released-solve-path",action="store_true")
    args=ap.parse_args()
    cases=args.case or [f"n32_c{k:02}" for k in range(1,9)]
    sides={int(case.split("_")[0][1:]) for case in cases}
    for side in sorted(sides):
        if side!=32:
            raise ValueError("This bounded GN endpoint runner is scoped to 32²; native trials use other proof budgets")
        disk=DiskGreen(side,400)
        bare=disk.dense()
        for q in args.q:
            group=(OUT/"extended_budget" if args.followup else OUT)/f"n{side}_q{q:02}"
            group.mkdir(parents=True,exist_ok=True)
            init=time.perf_counter(),time.process_time()
            op_id,op=public_operator(side,q)
            checker=VerifiedWave(op)
            builders={}
            generated=op2(disk,q,0,bare)
            if "G" not in op and not args.released_solve_path:
                # Same binary64 Green coefficients as the original FFT state.
                # The dense matrix is a factorization proposal on this 32² grid.
                op["G"]=generated["G"]
            disagreement={k:float(np.max(abs(generated[k]-op[k]))) for k in
                          ("S","E0","spectrum","low_rank_left","low_rank_right","M_nominal")}
            if max(disagreement.values())>1e-12:
                raise ArithmeticError(f"Nominal recomputation differs for q={q}")
            for e in ETAS[1:]:
                other=op2(disk,q,e,bare)
                delta=np.linalg.norm(other["M_nominal"]-op["M_nominal"],2)
                builders[e]=other,VerifiedWave(other),float(delta)
            write_json(group/"setup.json",dict(q=q,side=side,nominal_agreement=disagreement,
                      input_operator_id=op_id,known_synthetic_replay=True,
                      setup_wall_s=time.perf_counter()-init[0],setup_cpu_s=time.process_time()-init[1]))
            for case in [a for a in cases if a.startswith(f"n{side}_")]:
                path=group/f"{case}.json"
                if args.skip_existing and path.exists() and "calibration_replay" in json.loads(path.read_text()):
                    continue
                if args.replay_only:
                    report=json.loads(path.read_text())
                    with np.load(path.with_suffix(".npz")) as retained:
                        mask=retained["roi"]
                        endpoints=[retained[k] if k in retained else None for k in ("low","high")]
                    low,high=endpoints
                else:
                    report,mask,low,high=nominal_one(case,q,op_id,op,checker,args.limit,args.matched_warm,args.followup,args.released_solve_path)
                    retained={"roi":mask}
                    if low is not None:retained["low"]=low
                    if high is not None:retained["high"]=high
                    write_arrays(path.with_suffix(".npz"),**retained)
                    write_json(path,report)
                    print(case,q,"nominal",report["nominal_inner_peak_range"],round(report["wall_s"],1),flush=True)
                replay=calibration_replay(case,q,op_id,builders,{
                      tag:(x,float(x[mask].max()) if x is not None else None)
                      for tag,x in (("low",low),("high",high))})
                write_json(path,report|dict(calibration_replay=replay))
                print(case,q,"η replay",[a["inner_width"] for a in replay],flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
