"""Frozen role-swap observation control on the original loaded q8 array."""
from pathlib import Path
import argparse
import csv
import json
import os
import sys

for var in ("OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","OMP_NUM_THREADS"):
    os.environ[var] = "1"
sys.dont_write_bytecode = True

import numpy as np
from scipy.ndimage import label
from scipy.optimize import Bounds,minimize
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PORT=ROOT/"ports"
DEST=ROOT/"results"/"aperture_control"
FREEZE=json.loads((HERE/"aperture_control_freeze.json").read_text())
sys.path.insert(0,str(ROOT/"imaging"))
sys.path.insert(0,str(PORT))
from shared import Model, elapsed, energy, now, write_json, write_arrays
from model import DiskGreen, Green, physical_solve


def dark_apertures(base):
    first_count=8
    drive=2*np.sqrt(FREEZE["termination_ohm"]/first_count)*np.eye(24,16)
    inv=np.linalg.solve(base["M_nominal"],np.eye(24))
    excitation=base["E0_bare"] @ (inv @ drive)
    receiver=-1j*np.sqrt(FREEZE["termination_ohm"])*(inv @ base["S_bare"])
    # Reference orientations use the exact earlier released coefficient arrays.
    excitation[:,:8]=base["E0"]
    receiver[8:16]=base["S"]
    return drive,excitation,receiver


def layouts(base, exc, rec):
    out={}
    for name, spec in FREEZE["layouts"].items():
        t=np.array(spec["tx"])
        s=np.array(spec["rx"])
        reference_power=float(FREEZE["available_source_W_total"])
        norm=np.sqrt(reference_power * 8 / len(t))
        out[name]=(t,s,norm,dict(base,E0=exc[:,t]*norm,S=rec[s]))
    return out


def load_pair(oid):
    p=ROOT/"results"/"extended"/f"{oid}.json"
    r=json.loads(p.read_text())
    if not 0 <= r["upper_residual_over_r"] < 1:
        raise ValueError("The connected-patch source must already be certified")
    with np.load(p.with_suffix(".npz"),allow_pickle=False) as a:
        patch=a["patch"]
        roi=a["roi"].ravel().astype(bool)
        positive=a["material"].ravel()
    with np.load((ROOT/r["zero_file"]).with_suffix(".npz"),allow_pickle=False) as a:
        zero=a["material"].ravel()
    cover=np.zeros(roi.size,bool)
    cover[patch]=True
    side=int(np.sqrt(roi.size))
    if label(cover.reshape(side,side))[1] != 1 or not roi[patch].all() or positive[patch].min()<2 or np.any(zero[roi]):
        raise ValueError("Incorrect frozen ROI material contract")
    return zero,positive,roi,patch,r


def independent_clean(base, drive, truth, bare):
    started=now()
    disk=DiskGreen(32,400)
    current, cost=physical_solve(truth,disk.kernel.apply,base["E0_bare"],
                                  disk.kernel.diagonal,bare)
    loaded=base["M_nominal"]-1j*(base["S_bare"] @ current)
    loaded_current=np.linalg.solve(loaded,drive)
    empty_current=np.linalg.solve(base["M_nominal"],drive)
    response=np.sqrt(float(base["reference_ohms"]))*(empty_current-loaded_current)
    return response,loaded,loaded_current,cost,elapsed(started)


def native_fixed(m, base, excitation, rec):
    start=now()
    g=Green(base)
    j,cost=physical_solve(m,g.apply,excitation,g.diagonal,base["G"])
    field=rec @ j
    return field,j,cost,elapsed(start)


class StopFit(Exception):
    pass


def fit(m, start, lower, upper):
    t=now()
    best=m.audit(start)
    x_best=start.copy()
    opt_iter=0
    trace=[]
    def objective(x):
        nonlocal best,x_best
        st=m.forward(x)
        info=m.audit(x)
        trace.append(dict(residual_over_r=st["residual_ratio"],**elapsed(t)))
        if info["residual_over_radius"] < best["residual_over_radius"]:
            best,x_best=info,x.copy()
        if best["residual_over_radius"] <= FREEZE["refit"]["discrepancy_target"]:
            raise StopFit("DISCREPANCY_TARGET")
        value, grad=m.data_grad(x)
        return value,grad
    def callback(x):
        nonlocal opt_iter
        opt_iter+=1
    try:
        ans=minimize(objective,start,jac=True,method="L-BFGS-B",bounds=Bounds(lower,upper),
                     callback=callback,
                     options=dict(maxiter=FREEZE["refit"]["accepted_update_cap_per_arm"],
                                  maxls=25,ftol=1e-13,gtol=1e-9))
        status=str(ans.message)
    except StopFit as e:
        status=str(e)
    state=m.forward(x_best)
    return x_best,state["j"].copy(),dict(best,status=status,updates=opt_iter,
        evaluations=len(trace),trace=trace,counts=dict(m.cost),**elapsed(t))


def compare(oid,base,definitions,drive,exc,rec,bare):
    started=now()
    zero, positive, roi, patch, proof = load_pair(oid)
    with np.load(PORT/"public"/"cases"/f"{oid}.npz",allow_pickle=False) as z:
        original_y=z["Y"]
        original_r=float(z["radius"])
    with np.load(PORT/"private"/f"{oid}.npz",allow_pickle=False) as z:
        truth=z["material"].ravel()
        original_clean=z["Y_clean"]
    qtrue,loaded,portcurrent,true_cost,true_time=independent_clean(base,drive,truth,bare)
    fz,jz,cz,tz=native_fixed(zero,base,exc,rec)
    fp,jp,cp,tp=native_fixed(positive,base,exc,rec)
    yparity=float(np.linalg.norm(qtrue[8:16,:8]-original_clean)/np.linalg.norm(original_clean))
    if yparity>1e-8:
        raise ArithmeticError("Independent loaded acquisition disagrees with original clean A")
    output=dict(case=oid,side=32,roi_index=proof["roi_index"],roi_cells=int(roi.sum()),
                patch_cells=int(len(patch)),max_original_clean_A_relative_difference=yparity,
                truth_generation_cost=true_cost,truth_generation_timing=true_time,
                fixed_costs=dict(zero=cz,patch=cp,zero_time=tz,patch_time=tp),by_layout={})
    arrays=dict(zero=zero,patch_candidate=positive,roi=roi,patch=patch,actual_full_unit_a=qtrue[:16])
    for name,(tx,rx,scale,operator) in definitions.items():
        t=now()
        ix=np.ix_(rx,tx)
        nchan=len(tx)*len(rx)
        r=original_r*np.sqrt(nchan/64)
        clean=qtrue[ix]*scale
        if name=="A_one_way":
            y=original_y
        else:
            layout_seed=1 if name=="B_interleaved" else 2
            rng=np.random.default_rng(np.random.SeedSequence([FREEZE["noise_seed_root"],
                                           int(oid[-2:]),layout_seed]))
            direction=rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape)
            y=clean+(2*r/3/np.linalg.norm(direction))*direction
        a,b=fz[ix]*scale,fp[ix]*scale
        zero_ratio=float(np.linalg.norm(a-y)/r)
        patch_ratio=float(np.linalg.norm(b-y)/r)
        if name=="A_one_way" and max(zero_ratio,patch_ratio)>=1:
            raise RuntimeError("Reference measurements no longer admit both fixed candidates")
        delta=float(np.linalg.norm(a-b)/(2*r))
        ii=portcurrent[:,tx]*scale
        aa=base["Z_AA"]+loaded-base["M_nominal"]
        accepted=np.real(np.sum(ii[tx,np.arange(len(tx))].conj()*(aa@ii)[tx,np.arange(len(tx))]))
        d=dict(radius=r,shape=list(y.shape),channels=nchan,shots=len(tx),tx=tx.tolist(),rx=rx.tolist(),
               s_per_column=float(scale),per_channel_radius=r/np.sqrt(nchan),
               actual_noise_over_r=float(np.linalg.norm(y-clean)/r),
               total_available_W=float(energy(drive[:,tx]*scale)/(4*FREEZE["termination_ohm"])),
               total_real_driven_accepted_W=float(accepted),
               fixed_separation_over_two_r=delta,fixed_zero_residual_over_r=zero_ratio,
               fixed_patch_residual_over_r=patch_ratio,
               fixed_clean_zero_error_over_r=float(np.linalg.norm(a-clean)/r),
               fixed_clean_patch_error_over_r=float(np.linalg.norm(b-clean)/r),
               fixed_balls_disjoint_numeric=bool(delta>1),
               fixed_compatible_count=int(zero_ratio<1)+int(patch_ratio<1))
        for key, arr in (("Y",y),("Y_clean",clean),("r",r),("F_zero",a),("F_patch",b),
                         ("J_zero",jz[:,tx]*scale),("J_patch",jp[:,tx]*scale)):
            arrays[name+"_"+key]=arr
        if name in FREEZE["refit"]["layouts"]:
            models=[Model(operator,y,r,name,32) for _ in range(2)]
            zero_limit=np.where(roi,0.,4.5)
            positive_floor=np.zeros_like(positive)
            positive_floor[patch]=2.
            x0,c0,m0=fit(models[0],zero,np.zeros_like(zero),zero_limit)
            x1,c1,m1=fit(models[1],positive,positive_floor,np.full_like(positive,4.5))
            if np.any(x0[roi]) or np.min(x1[patch])<2:
                raise ArithmeticError("Refit escaped its frozen material family")
            d.update(zero_refitted=m0,patch_refitted=m1,
                     refitted_compatible_count=int(m0["residual_over_radius"]<1)+int(m1["residual_over_radius"]<1),
                     refitted_separation_over_two_r=float(np.linalg.norm(operator["S"]@(c0-c1))/(2*r)))
            for key,arr in (("refit_zero",x0),("refit_patch",x1),("refit_J_zero",c0),("refit_J_patch",c1)):
                arrays[name+"_"+key]=arr
        d.update(**elapsed(t))
        output["by_layout"][name]=d
    output.update(**elapsed(started))
    write_arrays(DEST/f"{oid}.npz",**arrays)
    write_json(DEST/f"{oid}.json",output)
    print(oid,{k:(round(v["fixed_separation_over_two_r"],3),round(v["fixed_zero_residual_over_r"],3),
                  round(v["fixed_patch_residual_over_r"],3),
                  None if "zero_refitted" not in v else (round(v["zero_refitted"]["residual_over_radius"],3),
                  round(v["patch_refitted"]["residual_over_radius"],3))) for k,v in output["by_layout"].items()},flush=True)
    return output


def execute():
    started=now()
    DEST.mkdir(parents=True,exist_ok=True)
    with np.load(ROOT/FREEZE["reference"],allow_pickle=False) as f:
        base={k:f[k] for k in f.files}
    drive,exc,rec=dark_apertures(base)
    definitions=layouts(base,exc,rec)
    bare=DiskGreen(32,400).dense()
    shared_timing=elapsed(started)
    write_arrays(DEST/"geometry.npz",**{k:base[k] for k in ("G","spectrum","grid_shape","low_rank_left","low_rank_right",
                  "diagonal_dressed","positions","M_nominal","reference_ohms")},
                 drive=drive,E0_both=exc,S_both=rec,
                 **{f"{name}_{key}": val for name, v in definitions.items()
                    for key,val in (("E0",v[3]["E0"]),("S",v[3]["S"]))})
    results=[compare(oid,base,definitions,drive,exc,rec,bare) for oid in FREEZE["cases"]]
    fields=("case","layout","channels","shots","total_available_W","total_real_driven_accepted_W","radius",
            "fixed_separation_over_two_r","fixed_zero_residual_over_r","fixed_patch_residual_over_r","fixed_compatible_count",
            "refitted_separation_over_two_r","refitted_zero_r","refitted_patch_r","refitted_compatible_count","wall_s","cpu_s")
    with (DEST/"summary.csv").open("w") as fp:
        writer=csv.DictWriter(fp,fieldnames=fields,extrasaction="ignore")
        writer.writeheader()
        for case in results:
            for name,value in case["by_layout"].items():
                writer.writerow(dict(case=case["case"],layout=name,**value,
                    refitted_zero_r=value.get("zero_refitted",{}).get("residual_over_radius"),
                    refitted_patch_r=value.get("patch_refitted",{}).get("residual_over_radius")))
    write_json(DEST/"summary.json",dict(shared_setup=shared_timing,case_ids=FREEZE["cases"],
               strict_scope="fixed candidate and refitted numerical point evaluation, separate exact forecast audit",
               case_receipts=[r["case"]+".json" for r in results],**elapsed(started)))


def verify_case(oid):
    from verified_wave import VerifiedWave, fro, rd, up

    t=now()
    with np.load(DEST/"geometry.npz",allow_pickle=False) as z:
        geom={key:z[key] for key in z.files}
    with np.load(DEST/f"{oid}.npz",allow_pickle=False) as z:
        data={key:z[key] for key in z.files}
    names=FREEZE["strict_replay"]["layouts"]
    out={}
    retained_inverse={}
    for layout in names:
        s,b=geom[layout+"_S"],geom[layout+"_E0"]
        w=VerifiedWave(dict(geom,S=s,E0=b))
        y=data[layout+"_Y"]
        r=float(data[layout+"_r"])
        def audit(x, j, key):
            if key not in retained_inverse:
                retained_inverse[key]=w.inverse(x)
            inv, invlog=retained_inverse[key]
            rhs=x[:,None]*b
            err_rhs=up(rd(3)*float(x.max())*fro(b))
            err_j=w.solve_error(x,j,rhs,err_rhs,inv)
            f=s@j
            err_f=up(rd(s.shape[1])*fro(s)*fro(j) + fro(s)*err_j)
            return f,err_f,dict(inverse=invlog,operator_response_error=err_f,
                                current_error_bound=err_j)
        def ratio_bounds(v, err, divisor, *subtracted):
            gamma=rd(2)
            err=up(err + gamma*sum(fro(x) for x in subtracted))
            low=np.linalg.norm(v)/(1+rd(2*v.size))
            return [float(np.nextafter(max(0.,low-err*(1+rd(5)))/divisor,-np.inf)),
                    up((fro(v)+err)/divisor)]
        rec={}
        for option in ("fixed","refitted"):
            if option=="fixed":
                x0,x1=data["zero"],data["patch_candidate"]
                j0,j1=data[layout+"_J_zero"],data[layout+"_J_patch"]
                suffix=""
            else:
                x0,x1=data[layout+"_refit_zero"],data[layout+"_refit_patch"]
                j0,j1=data[layout+"_refit_J_zero"],data[layout+"_refit_J_patch"]
                suffix=layout
            p0,e0,c0=audit(x0,j0,"zero"+suffix)
            p1,e1,c1=audit(x1,j1,"patch"+suffix)
            rec[option]=dict(zero_check=c0,patch_check=c1,
                zero_residual_over_r=ratio_bounds(p0-y,e0,r,p0,y),
                patch_residual_over_r=ratio_bounds(p1-y,e1,r,p1,y),
                separation_over_two_r=ratio_bounds(p0-p1,e0+e1,2*r,p0,p1))
        out[layout]=rec
    write_json(DEST/f"{oid}_verified.json",dict(case=oid,by_layout=out,
       scope="One-sided exact stored binary64 nominal coefficients with Arb-enclosed DFT roots",
       **elapsed(t)))


def cross_audit():
    started=now()
    with np.load(ROOT/FREEZE["reference"],allow_pickle=False) as op:
        base={key:op[key] for key in op.files}
    stacked_s=np.vstack((base["S"],base["S_hold"]))
    rows=[]
    for case in FREEZE["cases"]:
        with np.load(PORT/"public"/"cases"/f"{case}.npz",allow_pickle=False) as obs:
            y,r,y_hold,r_hold=obs["Y"],float(obs["radius"]),obs["Y_hold"],float(obs["radius_hold"])
        record=json.loads((DEST/f"{case}.json").read_text())
        with np.load(DEST/f"{case}.npz",allow_pickle=False) as obj:
            for layout in FREEZE["refit"]["layouts"]:
                for branch in ("zero","patch"):
                    x=obj[layout+"_refit_"+branch]
                    p,_,cost,timing=native_fixed(x,base,base["E0"],stacked_s)
                    a_ratio=float(np.linalg.norm(p[:8]-y)/r)
                    h_ratio=float(np.linalg.norm(p[8:]-y_hold)/r_hold)
                    new=float(record["by_layout"][layout][branch+"_refitted"]["residual_over_radius"])
                    rows.append(dict(case=case,layout=layout,material_family=branch,
                         original_A_over_r=a_ratio,original_hold_over_r=h_ratio,new_actual_over_r=new,
                         fits_both_A_and_new=bool(max(new,a_ratio)<1),original_solve_LU=cost["LU"],
                         original_rhs=cost["physical_rhs"],**timing))
    with (DEST/"cross_audit.csv").open("w") as fp:
        writer=csv.DictWriter(fp,fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    write_json(DEST/"cross_audit_cost.json",dict(rows=len(rows),
      audit_scope="No further fit. Original A and fixed held rows evaluated from one original dressed state.",
      wall_s_total=elapsed(started)["wall_s"],cpu_s_total=elapsed(started)["cpu_s"]))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    phase=ap.add_mutually_exclusive_group()
    phase.add_argument("--verify-only",action="store_true")
    phase.add_argument("--cross-audit",action="store_true")
    cli=ap.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        if cli.verify_only:
            verify_case(FREEZE["strict_replay"]["case"])
        elif cli.cross_audit:
            cross_audit()
        else:
            execute()
