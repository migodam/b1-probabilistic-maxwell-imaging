"""Policy-fixed MIMO pairings, separate frequency states and objective phase diagram."""
import argparse
import csv
import json
import re
import time

import numpy as np
from scipy.linalg import eigh, norm
from threadpoolctl import threadpool_limits

from run import OPTICAL, profile
from sector import HERE, PUBLIC, UPPER, Green, adj, dark_trial, herm
from verify import ExactUniform, ctx, verify_dark


def geometry(green,s,b,upper,active,ratio):
    trial,scattered=dark_trial(green,s,active,rank=104)
    n=min(len(trial.T), len(s))
    mass=herm(adj(trial)@trial)
    prop=adj(trial)@scattered
    kr=mass-upper*herm(prop)
    ki=(prop-adj(prop))/(2j)
    rb=adj(s@trial)@(s@trial)
    modes=np.linalg.svd(s,compute_uv=False)
    source=np.linalg.svd(b,compute_uv=False)
    profiles=[]
    best=None
    for t in ratio:
        eig,v=eigh(herm(rb+kr+t*ki),driver="evd",check_finite=False)
        d=int((eig<-1e-8).sum())
        score=(d,-eig[n-1])
        profiles.append(dict(ratio=t,null_negative_lower_numeric=d,s_th_quadratic=float(eig[n-1])))
        if best is None or score < best[0]:
            best=score,trial@v[:,:n],t
    return profiles,best,dict(receiver_singular_values=modes.tolist(),incident_singular_values=source.tolist(),
                    receiver_10pct=int((modes>0.1*modes[0]).sum()),incident_10pct=int((source>0.1*source[0]).sum()))


class Transfer:
    def __init__(self,g,s,b,cap,ratio,w):
        mix=w*(cap+1j*ratio)/2
        symbol=w*(1-cap*g.spectrum.real+ratio*g.spectrum.imag)
        alpha=w/2
        keep=np.flatnonzero(symbol.ravel()<alpha)
        fy,fx=np.unravel_index(keep,g.spectrum.shape)
        y,x=np.indices((g.side,g.side))
        phase=(y.reshape(-1,1)*fy+x.reshape(-1,1)*fx)*(2j*np.pi/g.pad)
        f=np.column_stack((np.exp(phase)/g.pad,g.left,adj(g.right),adj(s)))
        r=g.left.shape[1]
        m=len(keep)
        p=len(s)
        d=np.r_[1/(symbol.ravel()[keep]-alpha),np.zeros(2*r),np.ones(p)]
        o=np.zeros((len(d),len(d)),complex)
        o.flat[::len(d)+1]=d
        o[m:m+r,m+r:m+2*r]=-np.eye(r)/mix.conjugate()
        o[m+r:m+2*r,m:m+r]=-np.eye(r)/mix
        small=herm(o+adj(f)@f/alpha)
        eig,v=eigh(small,driver="evd",check_finite=False)
        self.nu=int(np.sum(eig>0)-r-p)
        g0=np.column_stack((adj(s),mix*b))
        projection=adj(v)@(adj(f)@g0)
        self.schur=herm(adj(g0)@g0/alpha-adj(projection)@(projection/eig[:,None])/alpha**2)
        self.w=w
        self.ratio=ratio
        self.cap=cap

    def distance(self,y):
        if self.nu >= y.shape[1]:
            return 0.
        z=np.vstack((y,np.eye(y.shape[1])))
        e=eigh(herm(adj(y)@y-adj(z)@self.schur@z),eigvals_only=True,check_finite=False)[::-1]
        return float(np.maximum(e[max(0,self.nu):],0).sum()**0.5)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side",type=int,required=True)
    parser.add_argument("--strict",action="store_true")
    parser.add_argument("--operator-only",action="store_true")
    parser.add_argument("--stationary-mask",action="store_true")
    args=parser.parse_args()
    ctx.prec,ctx.threads=112,2
    clock=time.perf_counter()
    out={"side":args.side,"rows":[],"geometry":[],"frequency_products":[]}
    led=HERE.parent/"ports"/"results"/f"array_grid_n{args.side}.csv"
    source={}
    for z in csv.DictReader(led.open()) if led.exists() else ():
        key=(z["case_id"],int(z["frequency_MHz"]),int(z["q"]),int(float(z["radiation_scale_ohms"])),
             int(round(float(z["normalized_calibration_error"])*1000)))
        source[key]=z
    obs={p.stem.removesuffix("_mimo"):np.load(p) for p in sorted((PUBLIC/"cases").glob(f"n{args.side}_*_mimo.npz"))}
    if not obs and not args.operator_only:
        raise FileNotFoundError("Paired MIMO observations not ready")
    case=next(iter(obs),f"n{args.side}_c01")
    mask_path=HERE.parent/"imaging"/"results"/case
    mask_path=mask_path/("stationary/baseline.npz" if args.stationary_mask else "baseline.npz")
    out["mask_version"]="stationary" if args.stationary_mask else "published-stop"
    suffix=("_geometry" if args.operator_only else "")+("_stationary" if args.stationary_mask else "")
    target=HERE/"results"/f"mimo_{args.side}{suffix}.json"
    mask=np.load(mask_path)["roi"].ravel()
    for op in sorted((PUBLIC/"operators").glob(f"n{args.side}_f*_a*.npz")):
        if not re.fullmatch(r"n\d+_f\d+_a\d+(_c\d+)?",op.stem):
            continue
        key0=op.stem.split("_",1)[1]
        if not "_c" in key0:
            key0+="_c048"
        f,q,c=(int(x[1:]) for x in key0.split("_"))
        primary=f"{key0}_e000_Y"
        if not args.operator_only and primary not in obs[case]:
            continue
        p=np.load(op)
        green=Green.read(op,args.side)
        s,b=p["S"],p["E0"]
        results=[geometry(green,s,b,UPPER,np.arange(green.n),OPTICAL),
                 geometry(green,s,b,UPPER,np.flatnonzero(~mask),OPTICAL)]
        data=dict(frequency=f,sources=q,receivers=q,coupling=c,operator=op.stem,
                  roi_origin=case,roi_cells=int(mask.sum()),**results[0][2])
        for (spec,best,_),region in zip(results,("whole_prior","primary_roi_zero"),strict=True):
            data[region]=spec
            data[region+"_largest_surviving_tail_numeric"]=max(0,q-best[0][0])
            if args.strict and c==48:
                active=np.arange(green.n) if region=="whole_prior" else np.flatnonzero(~mask)
                pcase=profile(np.full(green.n,UPPER),s,optical=best[2])
                data[region+"_strict_receiver_null"]=verify_dark(green,s,best[1],active,
                              pcase.high,pcase.weight,pcase.optical)
        out["geometry"].append(data)
        if args.operator_only:
            target.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
            print(op.stem, "geometry only", flush=True)
            continue
        scale=float(norm(s,2)**2)
        transfers=[Transfer(green,s,b,0.25,t,w*scale) for t in (0.,1.,4.) for w in (0.012,0.06,0.3)]
        for cid, arrays in obs.items():
            for e in (0,5,25,60):
                view=f"{key0}_e{e:03d}"
                if f"{view}_Y" not in arrays:
                    continue
                y=arrays[view+"_Y"]
                r=float(arrays[view+"_radius"])
                cal=source[(cid,f,q,c,e)]
                enlargement=float(cal["Y_train_calibration_bound_absolute"])
                cal_full=float(cal["Y_train_calibration_bound_full"])
                trial=max(((t.distance(y),t) for t in transfers),key=lambda a:a[0])
                row=dict(case=cid,frequency=f,channels=q,coupling=c,calibration_error=e/1000,radius=r,
                         calibration_enlargement_absolute=enlargement,
                         calibration_enlargement_if_full_transfer_independently_observed=cal_full,
                         global_cap_test=0.25,global_relaxed_distance=trial[0],surrogate_nu=trial[1].nu,
                         margin_nominal_radius_only=trial[0]-r,margin_absolute_envelope=trial[0]-r-enlargement,
                         margin_full_transfer_sensitivity=trial[0]-r-cal_full,
                         width_ordinary_primary_roi_prior=UPPER,
                         roi_tail_numeric_upper=data["primary_roi_zero_largest_surviving_tail_numeric"])
                if args.strict and c==48 and cid==case and (e==0 or trial[0] > r+enlargement):
                    proof=ExactUniform(green,s,b,y).proof(0.25,trial[1].ratio,trial[1].w,r+enlargement)
                    row["arb_global"]=proof
                out["rows"].append(row)
        print(op.stem, "global/ROI compressed min",
              *(min(a["null_negative_lower_numeric"] for a in t[0]) for t in results), flush=True)
        out["wall_seconds"]=time.perf_counter()-clock
        target.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    by={}
    for z in out["rows"]:
        by.setdefault((z["case"],z["channels"],z["coupling"],z["calibration_error"]),[]).append(z)
    for key, rows in by.items():
        if len(rows)<2:
            continue
        joint_d2=sum(i["global_relaxed_distance"]**2 for i in rows)
        joint_r2=sum((i["radius"]+i["calibration_enlargement_absolute"])**2 for i in rows)
        row=dict(case=key[0],channels=key[1],coupling=key[2],calibration_error=key[3],
                 independent_frequencies=sorted(i["frequency"] for i in rows),
                 total_basis_columns=sum(i["channels"] for i in rows),
                 joint_squared_distance=joint_d2,joint_squared_radius=joint_r2,joint_margin=joint_d2-joint_r2,
                 any_frequency_has_positive_margin=any(i["margin_absolute_envelope"]>0 for i in rows),
                 zero_roi_tail_columns_upper=sum(i["roi_tail_numeric_upper"] for i in rows))
        out["frequency_products"].append(row)
    out["wall_seconds"]=time.perf_counter()-clock
    target.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    with threadpool_limits(1):
        main()
