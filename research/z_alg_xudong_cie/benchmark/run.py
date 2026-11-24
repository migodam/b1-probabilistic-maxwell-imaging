"""Fit training voltages under a real elapsed-time budget. No material or held-port input."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

for name in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
from scipy.fft import ifft2
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
RESEARCH = STUDY.parent
for directory in (RESEARCH/"srap_cie_array_study_20261121", STUDY/"holdout", STUDY/"chambers",
                  STUDY/"cross_cie", STUDY/"pruning", STUDY/"minpeak_native", STUDY/"diagnosis"):
    sys.path.insert(0,str(directory))
import data as collection
import training as chambers
import engine
import cross_engine
from search import Model as NativeForward
from branches import propose_branches, expand, as_branch
from native_box import NativeBox
from diagnose import get_case as development_case, CASES as DEV

def now():
    return datetime.now(timezone.utc).isoformat()

def clocks():
    return perf_counter(), process_time()

def spent(start):
    return dict(wall=perf_counter()-start[0],cpu=process_time()-start[1])

def load(name):
    if name in chambers.CASES:
        p=chambers.load_training(name)
        case=dict(p["metadata"],suite="LA_FRESNEL",side=p["grid_shape"][0])
        if chambers.CASES[name].acquisition=="Manitoba":
            case.update(suite="MANITOBA",public_schedule=dict(modes=[6,10,16],betas=[6.,2.,.5],rank=5))
        return case,p,None
    if name in DEV:
        case,p=development_case(name)
        key=case.get("array")
        if name.startswith(("optfresh","s2d")):
            key="PW180_10x20"
            fine=RESEARCH/"z_alg_optical_inertia/synthetic/public/fine_geometry.npz"
            with np.load(fine,allow_pickle=False) as a:
                native=dict(E=a[f"{key}__B"],S=a[f"{key}__S"],spectrum=a["spectrum"])
            return case,p,native
        return case,p,None
    case,p=collection.get_case(name)
    return case,p,collection.get_native_case(name)[1]

def make_model(p,case,clock):
    model=engine.legacy.PaperModel(p,case,clock)
    if "spectrum" in p:
        model.green.spectrum=p["spectrum"]
        model.green.pad=p["spectrum"].shape
        model.green.diagonal=ifft2(p["spectrum"])[0,0]
    return model

class TimedNative(NativeForward):
    def green(self,z,adjoint=False):
        self.clock.check()
        return super().green(z,adjoint)

def lift(m,n):
    side=int(np.sqrt(len(m)))
    if n%side:
        raise ValueError("Material lift must preserve the box exactly.")
    return np.repeat(np.repeat(m.reshape(side,side),n//side,axis=0),n//side,axis=1).ravel()

def generators(model,p,case,method,initial,upper,settings,measured):
    if method=="dual" and np.isscalar(upper) and upper==20.:
        return engine.run_cie_arm(model,p["Y"],case,"dual",initial,float(p["delta"]))
    return cross_engine.iterate(model,p["Y"],case,method,initial,upper,passive=measured,
                                **settings.get(method,{}))

class Benchmark:
    def __init__(self,name,arm,config,*,shared=None,allocated_wall=0.,allocated_cpu=0.):
        self.name,self.arm,self.config=name,arm,config
        self.lab=name in chambers.CASES
        self.shared=shared
        self.allocated_wall,self.allocated_cpu=allocated_wall,allocated_cpu
        self.total_budget=config["lab_seconds" if self.lab else "synthetic_seconds"]
        self.deadline=self.total_budget-allocated_wall
        self.search_limit=self.deadline-config["lab_audit_seconds" if self.lab else "synthetic_audit_seconds"]
        self.clock=engine.legacy.Clock(self.deadline)
        self.trace=[]; self.runs=[]; self.boxes=[]; self.native_rows=[]
        self.case,self.p,self.fine=load(name)
        self.model=make_model(self.p,self.case,self.clock)
        self.setup=self.clock.wall()
        self.y=self.p["Y"];self.radius=float(self.p["delta"])
        self.air=np.zeros(self.model.n)
        self.air_phi=.5+float(self.model.reg.evaluate(self.air)[0])
        self.short=[];self.pool_size=0
        self.by_start={}
        self.last=self.air.copy()
        self.prefix_J=None;self.prefix_material=None;self.prefix_step=None
        self.fine_counts=Counter()
        self.flags=[]
        self.native_cost={"wall":0.,"cpu":0.}
        self.proof_cost={"wall":0.,"cpu":0.}
        self.winner=self.air.copy()
        self.final_phi,self.final_residual=self.air_phi,1.
        self.winner_label="air"
        self.candidates=[dict(id="air",coarse_phi=self.air_phi,coarse_residual=1.)]
        self.coarse_best=self.air
        self.coarse_phi=self.air_phi
        self.coarse_id="air"

    def audit_coarse(self,m,step,stage,label):
        price=clocks()
        val,res,_=self.model.objective(m,self.y)
        self.clock.check()
        if not np.isfinite(val):
            raise ArithmeticError("Nonfinite full-wave objective.")
        key=f"{label}/stage{stage}/step{step}"
        item=dict(id=key,coarse_phi=float(val),coarse_residual=float(res),
                  wall_at=self.clock.wall(),**spent(price))
        self.candidates.append(item);self.pool_size+=1
        if val<self.coarse_phi:
            self.coarse_phi,self.coarse_best,self.coarse_id=val,m.copy(),key
        if not self.fine:
            if val<self.final_phi:
                self.winner,self.final_phi,self.final_residual,self.winner_label=m.copy(),val,res,key
        else:
            if self.arm in ("dual_then_restart","dual_then_pixel"):
                if label not in self.by_start or val<self.by_start[label][0]:
                    self.by_start[label]=(float(val),m.copy(),key)
            if all(np.linalg.norm(m-v[1]) > 1e-5*np.sqrt(len(m)) for v in self.short):
                self.short.append((float(val),m.copy(),key))
                self.short.sort(key=lambda v:v[0])
                self.short=self.short[:self.config["shortlist"]]

    def trajectory(self,method,label,initial,upper,cut,*,prefix=False):
        self.clock.limit=min(self.search_limit,cut)
        if self.clock.wall()>=self.clock.limit:
            self.runs.append(dict(label=label,status="no_time"));return
        start=clocks(); before=self.clock.counts.copy(); stop="schedule_complete"; last=0
        settings=self.config["method_settings"]
        if self.arm=="dual_then_pixel" and method=="pixel":
            settings={"pixel":dict(strength=.25,sector_strength=0.,cadence=16,trust=.08,
                                    stationary_stages=False,unconditional_pr_reset=False)}
        gen=generators(self.model,self.p,self.case,method,initial.copy(),upper,settings,self.lab)
        try:
            for m,step,stage,event in gen:
                self.clock.check()
                self.last=m.copy();last=step
                frame=gen.gi_frame.f_locals
                if prefix:
                    if step>=self.config["prefix_steps"]:
                        self.prefix_material=m.copy();self.prefix_J=frame["J"].copy()
                        self.prefix_step=step
                        self.audit_coarse(m,step,stage,label)
                        stop="branch_prefix_complete";break
                    continue
                end=bool(frame.get("end",False))
                if step in (1,64) or step%192==0 or end:
                    self.audit_coarse(m,step,stage,label)
        except engine.legacy.TimeLimit:
            stop="budget"
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as err:
            stop="solver_failure";self.flags.append(f"{label}: {err}")
        finally:
            gen.close()
        self.runs.append(dict(label=label,method=method,steps=last,status=stop,
                              counts=dict(self.clock.counts-before),**spent(start)))

    def region_branches(self):
        cfg=self.config["branches"]
        return propose_branches(self.prefix_material,self.prefix_J,**cfg)["branches"]

    def chamber_checks(self,branches):
        from restricted import score_hard
        start=clocks()
        rows=self.model.ports.indices(self.model.ns)
        common=sorted(set.intersection(*(set(v.tolist()) for v in rows)))
        if len(common)>16:
            common=[common[i] for i in np.linspace(0,len(common)-1,16,dtype=int)]
        if common:
            packed=np.column_stack([self.y[[list(rr).index(j) for j in common],s]
                                   for s,rr in enumerate(rows)])
            geometry=NativeBox(self.model.green.spectrum,self.model.ports.S[common],self.model.E)
        for branch in branches:
            row=dict(branch=branch.name,tile_caps=branch.caps,decision="RETAIN",
                     reason="target/model error has no absolute laboratory bound",
                     common_receiver_ids=common)
            if common and branch.name in ("ranked_support","shifted_support") and self.clock.wall()<self.search_limit*.5:
                try:
                    conditional=score_hard(geometry,packed,self.radius,branch,t_grid=(8.,),
                                           alpha_grid=(2.**-14,),max_active=576,max_proofs=1)
                    row["conditional_finite_model"]=conditional
                except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
                    row["conditional_error"]=str(exc)
            self.boxes.append(row)
        self.proof_cost=spent(start)
        return branches

    def screen(self,branches):
        if self.lab:
            return self.chamber_checks(branches)
        if not self.fine:
            for branch in branches:
                self.boxes.append(dict(branch=branch.name,tile_caps=branch.caps,decision="RETAIN",
                    reason="physical target/model error bound unavailable" if self.lab else "native model unavailable"))
            return branches
        mark=clocks()
        import screen as pruning_screen
        table={}
        def remaining():
            return min(self.search_limit*self.config["proof_fraction"],self.search_limit)-self.clock.wall()
        if remaining()>0:
            cert=self.shared if self.shared is not None else pruning_screen.Checks(
                self.fine["spectrum"],self.fine["S"],self.fine["E"])
            for branch in branches:
                if remaining()<=0: break
                if branch.name=="wildcard":continue
                begun=clocks()
                try:
                    r=cert.score(self.y,self.radius,branch)
                    eligible=self.clock.wall()<=min(self.search_limit*self.config["proof_fraction"],self.search_limit)
                    table[branch.name]=dict(**r,within_proof_budget=eligible)
                    if not eligible:
                        table[branch.name]["decision"]="RETAIN"
                except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
                    table[branch.name]=dict(branch=branch.name,decision="RETAIN",reason=str(exc),**spent(begun))
        keep=[]
        for b in branches:
            r=table.get(b.name,dict(branch=b.name,decision="RETAIN",reason="open branch or proof time exhausted"))
            r["tile_caps"]=b.caps
            self.boxes.append(r)
            if r["decision"]!="PRUNE":
                keep.append(b)
        self.proof_cost=spent(mark)
        return keep

    def execute_search(self):
        cfg=self.config
        if self.arm in ("dual_then_restart","dual_then_pixel"):
            self.trajectory("dual","first_dual",self.air,20.,self.search_limit/2)
            is_pixel=self.arm=="dual_then_pixel"
            self.trajectory("pixel" if is_pixel else "dual","second_start",
                            self.air if is_pixel else self.air+.65,20.,self.search_limit)
            entries=list(self.by_start.values())+self.short
            self.short=[]
            for val,m,key in entries:
                if all(np.linalg.norm(m-v[1])>1e-5*np.sqrt(len(m)) for v in self.short):
                    self.short.append((val,m,key))
                if len(self.short)==cfg["shortlist"]:break
            return
        if self.arm not in ("multistart","matched_boxes","pruned","pruned_reuse"):
            self.trajectory(self.arm,self.arm,self.air,20.,self.search_limit)
            return
        if self.arm=="multistart":
            total=3
            for i,value in enumerate((0.,.65,1.5)):
                t=self.clock.wall()
                self.trajectory("dual",f"start_{value:g}",self.air+value,20.,t+(self.search_limit-t)/(total-i))
            return
        self.trajectory("dual","prefix",self.air,20.,self.search_limit,prefix=True)
        if self.prefix_material is None:
            self.flags.append("Prefix incomplete: open fallback.")
            self.trajectory("dual","fallback",self.air,20.,self.search_limit);return
        bs=self.region_branches()
        if self.arm in ("pruned","pruned_reuse"):
            bs=self.screen(bs)
        else:
            self.boxes=[dict(branch=b.name,tile_caps=b.caps,decision="RUN_CONTROL") for b in bs]
        for i,b in enumerate(bs):
            # The open branch always has an equal share of remaining search time.
            upper=20. if b.name=="wildcard" else expand(b,self.model.side)
            m=np.minimum(upper,self.prefix_material)
            t=self.clock.wall()
            self.trajectory("dual",b.name,m,upper,t+(self.search_limit-t)/(len(bs)-i))

    def finalize_native(self):
        if self.fine is None:return
        clock=self.clock;clock.limit=self.deadline
        start=clocks()
        solver=TimedNative(self.fine["spectrum"],self.fine["S"],self.fine["E"],self.y);solver.clock=clock
        yn=np.linalg.norm(self.y)
        for q,m,label in self.short:
            if clock.wall()>=self.deadline: break
            begin=clocks(); before=solver.count.copy()
            try:
                fr,_,_=solver.full(lift(m,96),gradient=False)
                clock.check()
                rel=float(fr/yn);phi=float(.5*rel**2+self.model.reg.evaluate(m)[0]);clock.check()
                self.native_rows.append(dict(id=label,native_residual=rel,native_phi=phi,
                                            coarse_phi=q,accepted_in_nominal_ball=bool(fr<=self.radius),
                                            counts=dict(solver.count-before),**spent(begin)))
                if phi<self.final_phi:
                    self.winner,self.final_phi,self.final_residual,self.winner_label=m.copy(),phi,rel,label
            except engine.legacy.TimeLimit:
                self.native_rows.append(dict(id=label,status="deadline",counts=dict(solver.count-before),**spent(begin)));break
            except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
                self.native_rows.append(dict(id=label,status="failed",message=str(exc),
                                            counts=dict(solver.count-before),**spent(begin)))
        self.native_cost=spent(start);self.fine_counts=solver.count.copy()

    def run(self,out):
        self.execute_search()
        self.finalize_native()
        rec=dict(observation=self.name,method=self.arm,setup_seconds=self.setup,
            wall_seconds=self.clock.wall()+self.allocated_wall,cpu_seconds=self.clock.cpu()+self.allocated_cpu,
            allocated_geometry_wall=self.allocated_wall,allocated_geometry_cpu=self.allocated_cpu,
            budget_seconds=self.total_budget,
            search_wall_cap=self.search_limit,native_rank_cost=self.native_cost,proof_cost=self.proof_cost,
            original_operator="lift32_to_native96" if self.fine else "laboratory_native_or_development_grid",
            selected=self.winner_label,relative_native_training_residual=self.final_residual,
            selected_phi=self.final_phi,coarse_best_id=self.coarse_id,coarse_best_phi=self.coarse_phi,
            radius=self.radius,training_norm=float(np.linalg.norm(self.y)),search_runs=self.runs,
            certificates=self.boxes,coarse_audits=self.candidates,native_audits=self.native_rows,
            coarse_counts=dict(self.clock.counts),native_counts=dict(self.fine_counts),
            execution_notes=self.flags,public_only=True)
        rec["radius_status"]=("numerical_in_ball" if self.final_residual*np.linalg.norm(self.y)<=self.radius else "outside_ball")
        out.mkdir(parents=True,exist_ok=True)
        target=out/f"{self.name}__{self.arm}"
        io=clocks()
        collection.write_npz(target.with_suffix(".npz"),m=self.winner,coarse_best=self.coarse_best,
                             last_completed=self.last)
        rec["output_array_io"]=spent(io)
        rec["sealed_at"]=now()
        collection.write_json(target.with_suffix(".json"),rec)
        return [self.name,self.arm,round(rec["wall_seconds"],2),self.winner_label,round(self.final_residual,4)]

def run_observation(args):
    name,arms,idx,settings,out=args
    np.seterr(under="ignore")
    offset=idx%len(arms)
    with threadpool_limits(1):
        for arm in arms[offset:]+arms[:offset]:
            target=out/f"{name}__{arm}"
            if target.with_suffix(".json").exists() and target.with_suffix(".npz").exists():
                continue
            result=Benchmark(name,arm,settings).run(out)
            print(*result,flush=True)
    return name

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config",type=Path,default=HERE/"FREEZE.json")
    ap.add_argument("--output",type=Path,default=HERE/"predictions")
    ap.add_argument("--case",action="append")
    ap.add_argument("--method",action="append")
    ap.add_argument("--workers",type=int,default=3)
    a=ap.parse_args()
    settings=json.loads(a.config.read_text())
    cases=a.case or [r["case"] for r in collection.planned_cases()]+[
                      k for k,v in chambers.CASES.items() if v.phase=="confirmation"]
    tasks=[(c,a.method or settings["methods"],i,settings,a.output) for i,c in enumerate(cases)]
    with ProcessPoolExecutor(a.workers) as pool:
        work=[pool.submit(run_observation,t) for t in tasks]
        for f in as_completed(work): f.result()
