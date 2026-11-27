"""Targeted native96 inner search using just one previously acquired complex scalar."""
from fractions import Fraction
import json
import os
from pathlib import Path
import sys

for v in ("OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[v]="1"
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
ENV=HERE.parents[1]
CONTROL=HERE.parent/"n96"
sys.path.insert(0,str(ENV/"postmeasurement"))
sys.path.insert(0,str(ENV))

from attack import BareOpen, elapsed, real, timer
from expand import AllColumnJacobian, additional_candidates
from coherent_modes import DEN
from fluence_loop import readout
from matched_amplitude import read_pair
from passive_replay import VerifiedWave, fro, load_cover, rd, stored_response, up
from probe import Circuit
from flint import arb
import numpy as np
from scipy.optimize import lsq_linear


SIDE=96
BANDS=(400,460)
RX=np.arange(8,16)
TX=np.arange(8)


def write_file(path,data):
    with path.open("x") as f:
        json.dump(data,f,indent=2,allow_nan=False,ensure_ascii=False,default=lambda s:s.item())
        f.write("\n")


def controller():
    prop=json.loads((CONTROL/"proposal.json").read_text())
    decision=json.loads((CONTROL/"decision.json").read_text())["readings"]
    with np.load(CONTROL/"controls.npz",allow_pickle=False) as a:
        integers=[a[k] for k in ("tx_ticks","rx_ticks","amp_ticks")]
        dt,dr,ds=integers
        tx=(dt[0]+1j*dt[1])/DEN
        rx=(dr[0]+1j*dr[1])/DEN
        ap=ds/DEN
        A=a["old_observed"]
    with np.load(CONTROL/"new_readings.npz",allow_pickle=False) as d:
        assert set(d.files)==set(decision)
        observed={k:(A.copy(),d[k]) for k in d.files}
    if prop["shots"]!=1 or tx.shape!=(8,1) or rx.shape!=(8,1) or ap.shape!=(1,) or A.shape!=(8,8):
        raise ValueError("Need the one-scalar frozen design")
    if any(s[1].shape!=(1,) for s in observed.values()):
        raise ValueError("The added acquired wave must be only one complex scalar")
    for ticks in (dt,dr,ds):
        if int(np.sum(ticks*ticks))>DEN**2:
            raise ArithmeticError("Frozen weights exceed the unit norm or shot power")
    return prop,decision,observed,tx,rx,ap,integers


def scalar_map(matrix,tx,rx,ap,amplitude):
    p=np.zeros((16,16),complex)
    p[8:16,:8]=matrix
    return readout(p,tx,rx,ap,amplitude)


def project_columns(j,tx,rx,ap,amp):
    k=(rx[:,0].conj()[:,None]*tx[:,0][None,:]).ravel()
    complex_columns=j[:64]+1j*j[64:]
    value=(ap[0]*amp)*(k@complex_columns)
    return np.vstack((value.real,value.imag))


def proof(material,y,r,tx,rx,ap,amp):
    start=timer()
    if not np.all((material>=0)&(material<=4.5)):
        raise ArithmeticError("Outside admitted physical material box")
    ids=np.flatnonzero(material)
    notes={}
    source={}
    for band,measure in zip(BANDS,y):
        work=timer()
        w=Circuit(SIDE,band,ids)
        n=SIDE**2
        wave=VerifiedWave(dict(spectrum=w.symbol,grid_shape=(SIDE,SIDE),S=w.c,E0=w.b,
                               low_rank_left=np.zeros((n,1)),low_rank_right=np.zeros((1,n))))
        s,scerr,st=stored_response(w,wave,material[ids])
        transfer,e,_=load_cover(s,scerr,1.,0.,0.)
        if band==400:
            c=8**-.5
            ec=up(float(abs(1/arb(8).sqrt()-arb(c)).upper()))
            b=transfer[8:16,:8]
            f=b*c
            err=up((c+ec)*e+ec*fro(b)+rd(4)*c*fro(b))
        else:
            f=readout(transfer,tx,rx,ap,amp)
            err=up((abs(amp)*abs(ap[0]))*e+rd(80)*(fro(transfer)+fro(f)))
        subtraction_guard=up(rd(4)*(fro(f)+fro(measure)))
        gap=up((fro(f-measure)+err+subtraction_guard)/r)
        notes[str(band)]=dict(residual_upper_over_r=gap,numerical_residual_over_r=float(np.linalg.norm(f-measure)/r),
                              determinant_guard=st,response_error_fro_upper=err,
                              material_lu=1,bare_forward_rhs=25,inverse_rhs=len(ids),**elapsed(work))
        source[f"strict_{band}_prediction"]=f
        source[f"strict_{band}_scattering_25ports"]=s
    return notes,source,elapsed(start)


class Fit:
    def __init__(self,obs,r,tx,rx,ap,amp):
        self.obs=obs
        self.r=r
        self.tx,self.rx,self.ap,self.amp=tx,rx,ap,amp
        self.jac=[AllColumnJacobian(SIDE,400,1.),AllColumnJacobian(SIDE,460,8.)]
        self.trial_cost={"400":[],"460":[]}

    def project(self,g):
        return scalar_map(g,self.tx,self.rx,self.ap,self.amp)

    def values(self,x):
        a,b=[k.evaluate(x) for k in self.jac]
        rst=[real(a[0]-self.obs[0])/self.r,real(self.project(b[0])-self.obs[1])/self.r]
        js=[a[1]/self.r,project_columns(b[1],self.tx,self.rx,self.ap,self.amp)/self.r]
        return rst,js

    def trial_wave(self,support):
        models=[BareOpen(SIDE,400,support,1.),BareOpen(SIDE,460,support,8.)]
        for b,w in zip(BANDS,models):
            self.trial_cost[str(b)].append(w.cost)
        return models

    def distance(self,fields):
        a,b=fields
        return [float(np.linalg.norm(a-self.obs[0])/self.r),
                float(np.linalg.norm(self.project(b)-self.obs[1])/self.r)]


def stage(tag,label,obs,r,tx,rx,ap,amp,refs,src,roi):
    stamp=timer()
    challenger=1 if label=="high" else 0
    x=refs[challenger].copy()
    family=float(x[roi][0])
    if not np.all(x[roi]==family) or family!=(0. if label=="high" else 2.):
        raise ValueError("Wrong forced core")
    fit=Fit(obs,r,tx,rx,ap,amp)
    dist,js=fit.values(x)
    cells=additional_candidates(SIDE,x,js,dist,src,80)
    scope=np.unique(np.r_[src,cells])
    free=scope[~np.isin(scope,roi)]
    local=np.searchsorted(scope,free)
    fn=fit.trial_wave(scope)
    trace=[]
    initial=np.array([np.linalg.norm(q) for q in dist])
    best=x.copy()
    best_gap=max(initial)
    trust=.65
    certificate=None
    scalars={}
    for i in range(12):
        step=timer()
        a0,a1=[float(np.linalg.norm(z)) for z in dist]
        if max(a0,a1)<.998:
            o,fields,t=proof(x,obs,r,tx,rx,ap,amp)
            trace.append(dict(iteration=i,start=[a0,a1],strict_replay=o,attempted_seconds=t))
            if max(z["residual_upper_over_r"] for z in o.values())<1:
                certificate=o
                scalars=fields
                best=x.copy()
                best_gap=max(a0,a1)
                break
        # Solve a real 128+2 row subproblem. The second frequency contributes
        # solely its one acquired complex scalar and no other port response.
        bias=float(np.clip((a0/max(a1,1e-8))**1.5,.5,2.25))
        code=np.vstack((bias*js[0][:,free],js[1][:,free]))
        y=-np.r_[bias*dist[0],dist[1]]
        d=lsq_linear(code,y,bounds=(np.maximum(-trust,-x[free]),np.minimum(trust,4.5-x[free])),
                     method="trf",lsq_solver="lsmr",lsmr_tol=4e-5,lsmr_maxiter=130,max_iter=20,tol=1e-5)
        current=x[scope].copy()
        trials=[]
        exact_stop=False
        for alpha in (1.,.5,.25):
            local_x=current.copy()
            local_x[local]=np.clip(current[local]+alpha*d.x,0,4.5)
            ys=[w.fields(local_x,jac=False)[0] for w in fn]
            gaps=fit.distance(ys)
            tx_x=np.zeros_like(x)
            tx_x[scope]=local_x
            trials.append((max(gaps),tx_x,gaps,alpha))
            if max(gaps)<.998:
                o,fields,t=proof(tx_x,obs,r,tx,rx,ap,amp)
                trace.append(dict(iteration=i,trial=gaps,strict_replay=o,attempted_seconds=t))
                if max(z["residual_upper_over_r"] for z in o.values())<1:
                    certificate=o
                    scalars=fields
                    best=tx_x
                    best_gap=max(gaps)
                    exact_stop=True
                    break
        q=min(trials,key=lambda z:z[0])
        print(tag,"core",family,"iter",i,"from",(a0,a1),"to",q[2],flush=True)
        trace.append(dict(iteration=i,old_ball_weight=bias,trust=trust,start=[a0,a1],tries=[z[2] for z in trials],
                          linear_iterations=int(d.nit),**elapsed(step)))
        if exact_stop:
            break
        if q[0]<best_gap:
            best=q[1]
            best_gap=q[0]
        if q[0]>=max(a0,a1)-1e-5:
            trust*=.45
            if trust<.075:
                break
        else:
            x=q[1]
            dist,js=fit.values(x)
            if q[0]<max(a0,a1)*.80:
                trust=min(trust*1.2,1.1)
    final_targeted_certificate=certificate
    final_replay_extra_cpu=0.
    if certificate is None:
        certificate,scalars,t=proof(best,obs,r,tx,rx,ap,amp)
        final_replay_extra_cpu=t["cpu_s"]
    if not np.all(best[roi]==family):
        raise ArithmeticError("The certificate candidate lost the core")
    is_fit=max(k["residual_upper_over_r"] for k in certificate.values())<1.
    # Reuse only the public other endpoint. It is separately audited against
    # the same acquired fields and can complete a certified diameter test.
    opposite=refs[1-challenger]
    other,other_fields,secondary_cost=proof(opposite,obs,r,tx,rx,ap,amp)
    other_is_fit=max(k["residual_upper_over_r"] for k in other.values())<1.
    with (HERE/f"{tag}_candidate.npz").open("xb") as f:
        np.savez_compressed(f,material_opposite_class=best,material_public_same_class=opposite,roi_indices=roi,
                            observed_old_400=obs[0],observed_added_scalar_460=obs[1],searched_territory=scope,
                            chosen_outer_sites=cells,**scalars,**{"same_"+k:v for k,v in other_fields.items()})
    max_edit=max(abs(Fraction(a)-Fraction(b)) for a,b in zip(refs[challenger],best))
    val=dict(sample=tag,initial_sealed_decision=label,core_opposite_class=family,
             core_same_class=float(opposite[roi][0]),core_cells=len(roi),physical_core_area_m2=len(roi)*(2/SIDE)**2,
             independent_strict_opposite=certificate,independent_strict_same=other,
             opposite_status="both balls verified" if is_fit else "unresolved",
             same_status="both balls verified" if other_is_fit else "unresolved",
             certified_possible_mean_diameter_lower=2 if is_fit and other_is_fit else None,
             tolerance_for_calibrated_reference=.01,max_edit_from_public_opposite_spec=up(float(max_edit)),
             active_scope_count=len(scope),added_exterior_coordinates=len(cells),trace=trace,
             trial_support_construction="Frozen public reference footprint plus one data-gradient selection of extra exterior sites",
             fullgrid_jacobian_cost={str(b):m.cost for b,m in zip(BANDS,fit.jac)},
             trial_physics_cost=fit.trial_cost,
             final_extra_replay_cpu=final_replay_extra_cpu,public_control_audit_cost=secondary_cost,
             early_certificate=(final_targeted_certificate is not None),**elapsed(stamp))
    write_file(HERE/f"{tag}.json",val)
    print("STRICT",tag,family,[v["residual_upper_over_r"] for v in certificate.values()],
          "same",[v["residual_upper_over_r"] for v in other.values()],is_fit,flush=True)
    return val


def main():
    start=timer()
    p,decision,obs,tx,rx,ap,ticks=controller()
    refs,supp,roi,r,oid=read_pair()
    if p["detector_radius"]!=r or np.sum(roi)!=164:
        raise ArithmeticError("Wrong frozen nominal radius or ROI")
    roi=np.flatnonzero(roi)
    outcomes={}
    for name,num in decision.items():
        label=num["label"]
        if label not in ("low","high"):
            raise ValueError("No sealed opposite class")
        outcomes[name]=stage(name,label,obs[name],r,tx,rx,ap,p["amplitude"],refs,supp,roi)
    write_file(HERE/"results.json",dict(side=SIDE,observation_id=oid,
              new_action=p["selected"],compensated_material="Shared real nondispersive box [0,4.5]",
              acquisition_inputs=["controls.npz::tx_ticks,rx_ticks,amp_ticks,old_observed",
                                  "proposal.json::amplitude,detector_radius","new_readings.npz","decision.json"],
              old_output_shape=[8,8],added_output_shape=[1],residual_real_rows_by_band=[128,2],
              old_normalization="one Watt over eight legacy shot columns, each readout scaled by sqrt(1/8)",
              new_normalization="unscaled A transfer projected by the exact frozen transmit and receive weights, then scaled by amp_ticks/4096 and recorded amplitude",
              recorded_amplitude=p["amplitude"],energy=p["available_energy_upper"],radius=r,
              outcomes=outcomes,**elapsed(start)))


if __name__=="__main__":
    main()
