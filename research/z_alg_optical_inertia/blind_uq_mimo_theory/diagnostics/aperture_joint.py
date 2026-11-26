"""Joint original and interleaved nominal full-wave data-space GN control."""
from collections import Counter
import argparse
import csv
import json

from aperture_control import HERE, ROOT, DEST, elapsed, now, write_arrays, write_json
import numpy as np
from scipy.linalg import cho_factor,cho_solve,lu_factor,lu_solve
from threadpoolctl import threadpool_limits

CFG=json.loads((HERE/"aperture_joint_freeze.json").read_text())
OUT=ROOT/"results"/"aperture_joint"
A,B="A_one_way","B_interleaved"


class Joint:
    def __init__(self,geom, data):
        self.G=geom["G"]
        self.B=geom["E0_both"][:,[0,1,2,3,4,5,6,7,8,10,12,14]]
        self.S=geom["S_both"][[1,3,5,7,8,9,10,11,12,13,14,15]]
        self.tx=((0,1,2,3,4,5,6,7),(0,2,4,6,8,9,10,11))
        self.rx=((4,5,6,7,8,9,10,11),(0,1,2,3,5,7,9,11))
        self.Y=[data[k+"_Y"] for k in (A,B)]
        self.r=[float(data[k+"_r"]) for k in (A,B)]
        self.counts=Counter()
        self.n=self.G.shape[0]

    def forward(self,x):
        self.counts["full_wave_LU"]+=1
        self.counts["forward_rhs"]+=self.B.shape[1]
        a=-x[:,None]*self.G
        a.flat[::self.n+1]+=1.
        factor=lu_factor(a,check_finite=False)
        block=x[:,None]*self.B
        j=lu_solve(factor,block,check_finite=False)
        e=self.B+self.G@j
        defect=float(np.linalg.norm(j-x[:,None]*e)/max(np.linalg.norm(block),1e-300))
        if defect>1e-9 or not np.all(np.isfinite(j)):
            raise ArithmeticError("Joint current failed original Green residual")
        sj=self.S@j
        blocks=[(sj[np.ix_(rx,tx)]-y)/r for rx,tx,y,r in zip(self.rx,self.tx,self.Y,self.r)]
        rs=np.concatenate([p.ravel() for p in blocks])
        ratios=[float(np.linalg.norm(p)) for p in blocks]
        return dict(x=x.copy(),state_j=j,field=e,LU=factor,ratios=ratios,
                    value=.5*float(np.vdot(rs,rs).real),residual=np.r_[rs.real,rs.imag],
                    state_defect=defect)

    def jacobian(self,st):
        self.counts["receiver_transpose_rhs"]+=self.S.shape[0]
        self.counts["jacobian_builds"]+=1
        transfer=lu_solve(st["LU"],self.S.T,trans=1,check_finite=False).T
        j=[]
        for rx,tx,r in zip(self.rx,self.tx,self.r):
            block=np.einsum("ri,is->rsi",transfer[list(rx)],st["field"][:,list(tx)]).reshape(-1,self.n)
            j.append(block/r)
        j=np.vstack(j)
        return np.vstack((j.real,j.imag))


def display(st):
    return dict(objective=st["value"],A_over_r=st["ratios"][0],B_over_r=st["ratios"][1],
                physical_state_relative_defect=st["state_defect"])


def minimize_joint(model, initial, lower, upper):
    start=now()
    prior=Counter(model.counts)
    n=model.n
    x=np.clip(initial,lower,upper)
    st=model.forward(x)
    initial_report=display(st)
    damping=CFG["initial_regularization_curvature_fraction"]
    trust=.85
    accepted=0
    last_pg=None
    trace=[]
    stop="MAX_MAJOR_ITERATIONS"
    for it in range(CFG["cap_major_iterations_each_start"]):
        if max(st["ratios"])<=CFG["reference_discrepancy_target_each_radius"]:
            stop="JOINT_DISCREPANCY"
            break
        jac=model.jacobian(st)
        r=st["residual"]
        grad=jac.T@r
        last_pg=float(np.max(np.abs(x-np.clip(x-grad,lower,upper))))
        entry=dict(old=display(st),major=it+1,accepted_updates=accepted,
                   projected_gradient_inf=last_pg,damping=damping,trust_cell=trust,trials=[])
        trace.append(entry)
        if last_pg<=CFG["projected_gradient_inf_tolerance"]:
            stop="PROJECTED_STATIONARITY"
            break
        free=(upper>lower)&~((x<=lower+1e-11)&(grad>=0))&~((x>=upper-1e-11)&(grad<=0))
        if not free.any():
            stop="NO_FREE_DESCENT_COORDINATES"
            break
        jf=jac[:,free]
        gram=jf@jf.T
        curvature=float(np.trace(gram)/len(gram))
        gram.flat[::len(gram)+1]+=max(1e-16,damping*curvature)
        direction=np.zeros(n)
        direction[free]=-jf.T@cho_solve(cho_factor(gram,check_finite=False),r,check_finite=False)
        direction*=min(1.,trust/max(np.max(np.abs(direction)),1e-300))
        improved=False
        for alpha in CFG["line_search"]:
            trial_x=np.clip(x+alpha*direction,lower,upper)
            trial=model.forward(trial_x)
            entry["trials"].append(dict(alpha=alpha,**display(trial)))
            if trial["value"] < st["value"]*(1-1e-10):
                x,st=trial_x,trial
                accepted+=1
                improved=True
                trust=min(1.7,trust*1.18)
                damping=max(2e-5,damping*.72)
                break
        if not improved:
            trust*=.64
            damping=min(10.,damping*4.)
            if trust<.008:
                stop="TRUST_REGION_STALLED"
                break
    if max(st["ratios"])<=CFG["reference_discrepancy_target_each_radius"]:
        stop="JOINT_DISCREPANCY"
    return st,dict(start=initial_report,final=display(st),accepted_updates=accepted,
                   major_iterations=len(trace),projected_gradient_inf=last_pg,stop=stop,
                   counts=dict(model.counts-prior),trace=trace,**elapsed(start))


def solve_case(oid,geom):
    t=now()
    inputs=("zero","patch_candidate","roi","patch",A+"_Y",A+"_r",B+"_Y",B+"_r",
            B+"_refit_zero",B+"_refit_patch")
    with np.load(DEST/f"{oid}.npz",allow_pickle=False) as d:
        data={key:d[key] for key in inputs}
    wave=Joint(geom,data)
    # Candidate inspection is shared and reported separately for equal budgets.
    before=now()
    names=("A_zero","A_patch","B_zero","B_patch")
    seeds=dict(zip(names,(data["zero"],data["patch_candidate"],data[B+"_refit_zero"],data[B+"_refit_patch"])))
    inspect={k:display(wave.forward(v)) for k,v in seeds.items()}
    paid_inspect=dict(counts=dict(wave.counts),**elapsed(before))
    zeros=np.zeros(wave.n)
    upper=np.full(wave.n,4.5)
    lower_patch=zeros.copy()
    lower_patch[data["patch"]]=2.
    seed_pair_a=min(names[:2],key=lambda name:inspect[name]["objective"])
    seed_pair_b=min(names[2:],key=lambda name:inspect[name]["objective"])
    branches=(
        ("roi_zero",("A_zero","B_zero"),zeros,np.where(data["roi"],0.,upper)),
        ("connected_patch_floor_2",("A_patch","B_patch"),lower_patch,upper),
        ("unrestricted",(seed_pair_a,seed_pair_b),zeros,upper),
    )
    results={}
    for key,initial_names,lo,hi in branches:
        options=[]
        costs=[]
        for name in sorted(initial_names,key=lambda name:inspect[name]["objective"]):
            state,info=minimize_joint(wave,seeds[name],lo,hi)
            info["seed"]=name
            options.append(state)
            costs.append(info)
            if info["stop"]=="JOINT_DISCREPANCY":
                break
        index=min(range(len(options)),key=lambda i:(not max(options[i]["ratios"])<1,options[i]["value"]))
        selected=options[index]
        x=selected["x"]
        if np.any(x<lo) or np.any(x>hi):
            raise ArithmeticError("Joint selection escaped its frozen class")
        results[key]=dict(chosen_start=costs[index]["seed"],**display(selected),
                          joint_in_ball_num=bool(max(selected["ratios"])<1),starts=costs,
                          final_current=selected["state_j"],material=x)
    arrays={}
    for name, log in results.items():
        arrays[name+"_J"]=log.pop("final_current")
        arrays[name+"_x"]=log.pop("material").reshape(32,32)
    write_arrays(OUT/f"{oid}.npz",**arrays)
    result=dict(case=oid,fit_families=results,paid_seed_inspection=paid_inspect,
                seed_metrics=inspect,total_counts=dict(wave.counts),**elapsed(t))
    write_json(OUT/f"{oid}.json",result)
    print(oid,{n:(round(v["A_over_r"],4),round(v["B_over_r"],4),v["joint_in_ball_num"]) for n,v in results.items()},flush=True)
    return result


def run(selected=None):
    t=now()
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(DEST/"geometry.npz",allow_pickle=False) as f:
        geom=dict(f)
    cases=CFG["cases"] if selected is None else selected
    for case in cases:
        solve_case(case,geom)
    if selected is None:
        aggregate()
    else:
        write_json(OUT/f"batch_{cases[0]}_{cases[-1]}.json",dict(cases=cases,**elapsed(t)))


def aggregate():
    results=[json.loads((OUT/f"{c}.json").read_text()) for c in CFG["cases"]]
    with (OUT/"summary.csv").open("w") as f:
        keys=("case","family","A_over_r","B_over_r","objective","joint_in_ball_num","chosen_start",
              "accepted_updates_total","major_iterations_total","fit_wall_s","fit_cpu_s","full_wave_LU","forward_rhs")
        tab=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore")
        tab.writeheader()
        for c in results:
            for family,branch in c["fit_families"].items():
                z=branch["starts"]
                tab.writerow(dict(case=c["case"],family=family,**branch,
                    accepted_updates_total=sum(r["accepted_updates"] for r in z),
                    major_iterations_total=sum(r["major_iterations"] for r in z),
                    fit_wall_s=sum(r["wall_s"] for r in z),fit_cpu_s=sum(r["cpu_s"] for r in z),
                    full_wave_LU=sum(r["counts"].get("full_wave_LU",0) for r in z),
                    forward_rhs=sum(r["counts"].get("forward_rhs",0) for r in z)))
    write_json(OUT/"summary.json",dict(cases=CFG["cases"],total_available_W=2,
                   acquisition_shots=16,shared_physical_rhs_per_eval=12,
                   object_wall_s_sum=sum(r["wall_s"] for r in results),
                   object_cpu_s_sum=sum(r["cpu_s"] for r in results),
                   costs_are_for_complete_actual_trials=True))


def verify():
    from verified_wave import VerifiedWave, fro, rd, up
    start=now()
    cases=[]
    for case in CFG["cases"]:
        doc=OUT/f"{case}.json"
        if not doc.exists():
            break
        cases.append(json.loads(doc.read_text()))
    with np.load(DEST/"geometry.npz",allow_pickle=False) as f:
        geom=dict(f)
    tx=((0,1,2,3,4,5,6,7),(0,2,4,6,8,9,10,11))
    output={}
    for family in CFG["families"]:
        first=next((r["case"] for r in cases if r["fit_families"][family]["joint_in_ball_num"]),None)
        if first is None:
            continue
        with np.load(OUT/f"{first}.npz",allow_pickle=False) as d:
            x,j=d[family+"_x"].ravel(),d[family+"_J"]
        with np.load(DEST/f"{first}.npz",allow_pickle=False) as d:
            data={name+key:d[name+key] for name in (A,B) for key in ("_Y","_r")}
        q={}
        inv=None
        for name,cols in zip((A,B),tx):
            op=dict(geom,S=geom[name+"_S"],E0=geom[name+"_E0"])
            w=VerifiedWave(op)
            if inv is None:
                inv,proof=w.inverse(x)
            js=j[:,list(cols)]
            y,r=data[name+"_Y"],float(data[name+"_r"])
            rhs=x[:,None]*w.B
            jerror=w.solve_error(x,js,rhs,up(rd(3)*x.max()*fro(w.B)),inv)
            wave=w.S@js
            err=up(rd(len(x))*fro(w.S)*fro(js)+fro(w.S)*jerror+
                   rd(2)*(fro(wave)+fro(y)))
            residual=up(fro(wave-y)+err)
            q[name]=dict(radius=r,residual_upper=residual,residual_upper_over_r=up(residual/r),
                         current_error_upper=jerror,numeric=float(np.linalg.norm(wave-y)/r))
        output[family]=dict(case=first,inverse=proof,layouts=q,verified=all(v["residual_upper_over_r"]<1 for v in q.values()))
    write_json(OUT/"verified.json",dict(candidates=output,
       earliest_completed_case_prefix=[r["case"] for r in cases],**elapsed(start)))
    print({k:(v["case"],v["verified"],[a["residual_upper_over_r"] for a in v["layouts"].values()])
           for k,v in output.items()},flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    group=ap.add_mutually_exclusive_group()
    group.add_argument("--verify-only",action="store_true")
    group.add_argument("--cases",nargs="+",choices=CFG["cases"])
    group.add_argument("--aggregate",action="store_true")
    arg=ap.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        if arg.verify_only:
            verify()
        elif arg.aggregate:
            aggregate()
        else:
            run(arg.cases)
