"""Fit actual frozen two-frequency data with the core held at either 0 or 2."""
import argparse
import json
from pathlib import Path
import os
import sys
from time import perf_counter, process_time

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
    os.environ[name]="1"
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
ENV=HERE.parent
sys.path.insert(0,str(ENV))
from probe import Circuit, R, at_load
from passive_replay import VerifiedWave, fro, up, rd, stored_response, load_cover
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares, minimize
import numpy as np
from flint import arb

TX=np.arange(8)
RX=np.arange(8,16)
BANDS=(400,460)
DISCREPANCY=0.975
GEOMETRY_ORIGIN=(.125,.0625)


def timer():
    return perf_counter(),process_time()


def elapsed(mark):
    return dict(wall_s=perf_counter()-mark[0],cpu_s=process_time()-mark[1])


def write_json(file,body):
    dest=HERE/file
    tmp=dest.with_suffix(".tmp")
    tmp.write_text(json.dumps(body,indent=2,allow_nan=False,default=lambda x:x.item())+"\n")
    tmp.replace(dest)


def real(a):
    return np.r_[a.real.ravel(),a.imag.ravel()]


def read(side):
    prior=ENV/(f"phase_{side}_two_shells.npz" if side==32 else "wave_96_2.0_1.npz")
    with np.load(prior,allow_pickle=False) as q:
        ix=q["ix"]
        origgroup=q["group"]
        keys=("plus_2.0","minus_2.0") if side==32 else ("plus","minus")
        x0,x1=[q[k] for k in keys]
        r=float(q["r"])
    branch=ENV/("low_energy_32" if side==32 else "low_energy_96_amplitude2")
    old_path=ENV/"acquisition_32"/"old_and_predictions.npz" if side==32 else branch/"selector.npz"
    old_key="old" if side==32 else "old_reading_total"
    with np.load(old_path,allow_pickle=False) as f:
        old=f[old_key]
    proposal=json.loads((branch/"proposal.json").read_text())
    decisions=json.loads((branch/"decision.json").read_text())["results"]
    energy=(1.,float(proposal["operating_energy"]))
    if r!=float(proposal["radius"]):
        raise ArithmeticError("Inconsistent radius")
    with np.load(branch/"readings.npz",allow_pickle=False) as f:
        if set(f.files)!=set(decisions):
            raise ValueError("Unexpected keys in actual low-power measurement")
        measured={name:(old.copy(),f[name]) for name in decisions}
    if any(a.shape!=(8,8) for pair in measured.values() for a in pair):
        raise ValueError("Unexpected recorded readout shape")
    core=ix[origgroup==0]
    centers={0.:np.zeros(side**2),2.:np.zeros(side**2)}
    centers[0.][ix]=x1
    centers[2.][ix]=x0
    assert np.all(x0[origgroup==0]==2.) and np.all(x1[origgroup==0]==0.)
    in_paths=dict(public_initial_material=prior.name,old=f"{old_path.relative_to(ENV)}::{old_key}",
                  actual_new=str((branch/"readings.npz").relative_to(ENV)),
                  new_energy=str((branch/"proposal.json").relative_to(ENV))+"::operating_energy",
                  prior_classification=str((branch/"decision.json").relative_to(ENV)))
    return ix,core,centers,measured,{k:v["label"] for k,v in decisions.items()},r,energy,in_paths


class BareOpen:
    """Bare volume solve with the exact circuit reduction for a zero-current probe."""
    def __init__(self,side,band,ix,energy):
        obj=Circuit(side,band,ix)
        self.G=obj.g
        self.B=obj.b[:,:24]
        self.C=obj.c[:24]
        self.Z0=obj.z0[:24,:24]+R*np.eye(24)
        self.amp=(energy/8)**.5
        self.cost=dict(material_lu=0,material_forward_rhs=0,material_adjoint_rhs=0,circuit_inverse=0)

    def fields(self,x,jac):
        a=-x[:,None]*self.G
        a.flat[::len(x)+1]+=1.
        lu=lu_factor(a,check_finite=False)
        j=lu_solve(lu,x[:,None]*self.B,check_finite=False)
        z=self.Z0-1j*(self.C@j)
        q=np.linalg.solve(z,np.eye(24))
        out=-2*R*self.amp*q[np.ix_(RX,TX)]
        self.cost["material_lu"]+=1
        self.cost["material_forward_rhs"]+=24
        self.cost["circuit_inverse"]+=1
        if not jac:
            return out,None
        u=lu_solve(lu,(q[RX]@self.C).T,trans=1,check_finite=False).T
        e=(self.B+self.G@j)@q[:,TX]
        grad=-2j*R*self.amp*(u[:,None,:]*e.T[None,:,:]).reshape(64,-1)
        self.cost["material_adjoint_rhs"]+=8
        return out,np.concatenate((grad.real,grad.imag),axis=0)


class CandidateReached(Exception):
    pass


class Fit:
    def __init__(self,side,ix,core_ids,t,obs,r,energy):
        self.side=side
        self.ix=ix
        self.roi=np.isin(ix,core_ids)
        self.ext=~self.roi
        self.y=np.array([real(a) for a in obs])
        self.r=r
        self.t=t
        self.model=[BareOpen(side,b,ix,w) for b,w in zip(BANDS,energy)]
        self.x=np.zeros(len(ix))
        self.x[self.roi]=t

    def start(self,full,degree):
        y,x=np.divmod(self.ix[self.ext],self.side)
        u=(x+.5)*2/self.side-1-GEOMETRY_ORIGIN[0]
        v=(y+.5)*2/self.side-1-GEOMETRY_ORIGIN[1]
        rho=np.hypot(u,v)
        theta=(np.arctan2(v,u)+2*np.pi)%(2*np.pi)
        if degree:
            radii=np.r_[np.linspace(.15,.50,13),.60,.75,1.0,1.5,2.0]
            ring=np.searchsorted(radii,rho,side="left")
            sector=np.floor(degree*theta/(2*np.pi)).astype(int)
            ids=ring*degree+np.minimum(sector,degree-1)
            _,self.design=np.unique(ids,return_inverse=True)
            eye=np.eye(int(self.design.max())+1)
            self.mask=eye[self.design]
            den=np.bincount(self.design)
            guess=np.bincount(self.design,weights=full[self.ix[self.ext]])/den
        else:
            self.design=np.arange(np.sum(self.ext))
            self.mask=None
            guess=full[self.ix[self.ext]]
        self.last={}
        self.best=None
        self.calls=0
        return np.clip(guess,0,4.5)

    def evaluate(self,params,jac=True):
        current=self.last
        if current and np.array_equal(current["params"],params) and ("j" in current or not jac):
            return current
        self.x[self.ext]=params[self.design]
        pred=[]
        act=[]
        for m in self.model:
            z,g=m.fields(self.x,jac=jac)
            pred.append(z)
            if jac:
                cropped=g[:,self.ext]
                act.append(cropped if self.mask is None else cropped@self.mask)
        dy=(np.array([real(x) for x in pred])-self.y)/self.r
        residuals=np.linalg.norm(dy,axis=1)
        out=dict(params=params.copy(),x=self.x.copy(),p=pred,res=dy,ratios=residuals)
        if jac:
            out["j"]=np.array(act)/self.r
        self.calls+=1
        self.last=out
        if self.best is None or max(residuals)<max(self.best["ratios"]):
            self.best=out
        return out

    def ls(self,z):
        info=self.evaluate(z)
        if max(info["ratios"])<DISCREPANCY:
            raise CandidateReached()
        return info["res"].ravel()

    def jac(self,z):
        return self.evaluate(z)["j"].reshape(256,-1)

    def balanced(self,z):
        info=self.evaluate(z)
        ratios2=np.sum(info["res"]**2,axis=1)
        if max(info["ratios"])<DISCREPANCY:
            raise CandidateReached()
        positive=np.maximum(ratios2-.96**2,0.)
        d=np.einsum("ai,aij->aj",info["res"],info["j"])
        return .25*(positive@positive),positive@d

    def run(self,full,degree,iters):
        start=timer()
        p0=self.start(full,degree)
        initial=self.evaluate(p0,jac=True)["ratios"].tolist()
        status="initial"
        nit=0
        if max(initial)>=DISCREPANCY:
            try:
                if degree and len(p0)<160:
                    out=least_squares(self.ls,p0,jac=self.jac,bounds=(0.,4.5),
                                      method="trf",max_nfev=iters,gtol=1e-8,xtol=1e-10,ftol=1e-10)
                    nit=int(out.nfev)
                else:
                    out=minimize(self.balanced,p0,jac=True,method="L-BFGS-B",
                                 bounds=[(0,4.5)]*len(p0),options=dict(maxiter=iters,maxls=20,ftol=1e-10,gtol=1e-8))
                    nit=int(out.nit)
                status=str(out.message)
            except CandidateReached:
                status="strict replay requested after numerical discrepancy"
        elif max(initial)<DISCREPANCY:
            status="public start is already inside numerical discrepancy"
        st=self.best
        x=np.zeros(self.side**2)
        x[self.ix]=st["x"]
        summary=dict(sectors=degree,num_free_coordinates=len(p0),initial=initial,best=st["ratios"].tolist(),
                     calls=self.calls,accepted_iteration_counter=nit,stop=status,**elapsed(start))
        return x,st["p"],summary

    def ledger(self):
        return {str(b):m.cost for b,m in zip(BANDS,self.model)}


def strict_replay(side,x,obs,r,energy,cores):
    out={}
    waves={}
    if np.any((x<0.)|(x>4.5)) or not np.all(x[cores]==x[cores][0]):
        raise ArithmeticError("Illegal candidate")
    ix=np.flatnonzero(x)
    for band,e,y in zip(BANDS,energy,obs):
        now=timer()
        w=Circuit(side,band,ix)
        dummy=dict(spectrum=w.symbol,grid_shape=(side,side),S=w.c,E0=w.b,
                   low_rank_left=np.zeros((side*side,1)),low_rank_right=np.zeros((1,side*side)))
        wave=VerifiedWave(dummy)
        s,serr,q=stored_response(w,wave,x[ix])
        open,oe,_=load_cover(s,serr,1.,0.,0.)
        a=open[np.ix_(RX,TX)]
        amp=(e/8)**.5
        delta_amp=up(float(abs(arb(e/8).sqrt()-arb(amp)).upper()))
        f=a*amp
        err=up((amp+delta_amp)*oe+delta_amp*fro(a)+rd(4)*amp*fro(a))
        residual=up((fro(f-y)+err+rd(3)*(fro(f)+fro(y)))/r)
        out[str(band)]=dict(strict_residual_over_r_upper=residual,numerical_residual_over_r=float(np.linalg.norm(f-y)/r),
                            response_error_upper=err,strict_material_lu=1,strict_inverse_rhs=len(ix),
                            strict_physical_rhs=25,certificate=q,**elapsed(now))
        waves[f"cert_open_{band}"]=f
        waves[f"cert_s_{band}"]=s
    return out,waves


def candidate(side,tag,family,start,cases):
    ix,core,centers,data,labels,r,energy,paths=cases
    clock=timer()
    profile=Fit(side,ix,core,family,data[tag],r,energy)
    if side==32:
        path=((4,60),(8,90),(0,120))
    else:
        path=((4,38),(8,45),(0,50))
    if start:
        path=path[:1]
    records=[]
    audits=None
    waves={}
    initial=centers[family]
    numeric_best=np.inf
    working=initial
    material=None
    for degree,iters in path:
        proposal,raw,numbers=profile.run(working,degree,iters)
        records.append(numbers)
        value=max(numbers["best"])
        print(side,tag,family,degree,numbers["initial"],numbers["best"],flush=True)
        if value<numeric_best:
            numeric_best=value
            material=proposal
            waves={"numeric_400":raw[0],"numeric_460":raw[1]}
        if value<1.002:
            receipt,raw_waves=strict_replay(side,proposal,data[tag],r,energy,core)
            numbers["verified"]=receipt
            if max(k["strict_residual_over_r_upper"] for k in receipt.values())<1.:
                audits=receipt
                material=proposal
                waves=dict(waves,**raw_waves)
                print("VERIFIED",side,tag,family,labels[tag],[(j,v["strict_residual_over_r_upper"]) for j,v in audits.items()],flush=True)
                break
        working=material
    key=f"candidate_{side}_{tag}_{int(family)}"
    x=material
    data_dump=dict(material=x,roi_indices=core,original_support=ix,
                   original_reference=centers[family],observed_400=data[tag][0],observed_460=data[tag][1],**waves)
    with (HERE/(key+".npz")).open("wb") as f:
        np.savez_compressed(f,**data_dump)
    ans=dict(side=side,specimen=tag,core_fixed=family,declared_earlier_class=labels[tag],
             is_opposite_to_earlier_class=(family==0.)==(labels[tag]=="high"),
             actual_observation_sources=paths,strict_status="both balls verified" if audits else "unresolved",
             independent_audit=audits,stages=records,
             allowed_exterior="All exterior cells within [0,4.5]. This search starts in the public support, explores shells, sectors then individual exterior cells there.",
             nx=side,core_cells=len(core),area_m2=len(core)*(2/side)**2,
             rho_abs=r,old_total_available_W=energy[0],new_total_available_W=energy[1],
             observed_layout="8x8 total calibrated power waves, 25th probe open, sum shot power divided by eight",
             optimizer_cost=profile.ledger(),**elapsed(clock))
    write_json(key+".json",ans)
    return ans


def main(side,scheduling):
    case=read(side)
    names=case[4]
    choices=[(tag,0. if label=="high" else 2.) for tag,label in names.items()]
    if scheduling!="opposite":
        choices+= [(tag,2. if label=="high" else 0.) for tag,label in names.items()]
    start=timer()
    result=[candidate(side,tag,family,False,case) for tag,family in choices]
    write_json(f"summary_{side}_{scheduling}.json",dict(side=side,arms=result,**elapsed(start)))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("side",type=int,choices=[32,96])
    p.add_argument("--schedule",choices=["opposite","all"],default="all")
    args=p.parse_args()
    main(args.side,args.schedule)
