"""Certify transferred blind ROI endpoints on the actual swept data balls."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import lu_factor, lu_solve

from model import COUNTS, R, STRENGTHS, DiskGreen, network_terms, port_point_matrices
from native_roi_path_verify import checked_solve, enclose_matrix, product_enclosure, scaled_rhs, times
from native_verify import ExplicitDots, gamma, tiny
from sweep import tag, union_positions
from verify_witness import as_acb, ball, frobenius, outward

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FIT=ROOT/"bounds/roi_mimo_fits"
CASES=[f"n32_c{i:02}" for i in range(1,9)]


def candidates(case):
    accepted={}
    mask=None
    for path in sorted(FIT.glob(f"**/{case}*.npz")):
        with np.load(path,allow_pickle=False) as data:
            trial_mask=data["roi"].ravel()
            if mask is None:
                mask=trial_mask
            if not np.array_equal(mask,trial_mask):
                raise ArithmeticError(f"Inconsistent frozen second ROI: {path}")
            receipt=json.loads(path.with_suffix(".json").read_text())
            if receipt.get("nominal_inner_peak_range") is None:
                continue
            evidence=receipt["nominal_endpoint_witnesses"]
            for key in ("low","high"):
                proof=evidence.get(key)
                if proof is None or not proof["verified_inside_noise"]:
                    continue
                x=np.asarray(data[key],dtype=np.float64).ravel()
                if x.shape != (1024,) or not np.all((0<=x)&(x<=4.5)):
                    raise ArithmeticError(f"Invalid material in {path}")
                peak=float(np.max(x[mask]))
                if peak != proof["peak"]:
                    raise ArithmeticError(f"Stored ROI peak disagrees with proof in {path}")
                origin=dict(path=str(path.relative_to(ROOT)),key=key,
                            nominal_source_array=int(receipt["q"]),proof_peak=peak)
                token=x.tobytes()
                if token in accepted:
                    accepted[token]["sources"].append(origin)
                else:
                    accepted[token]=dict(x=x,peak=peak,sources=[origin])
    if mask is None or not np.any(mask) or not accepted:
        raise ArithmeticError(f"Missing fixed ROI or legal maps for {case}")
    return mask,list(accepted.values())


def bare_certificate(x,G,E,S,dots):
    """One 72-column bare VIE solve with an explicit inverse-defect enclosure."""
    start=perf_counter(),process_time()
    ids=np.flatnonzero(x)
    x=x[ids]
    n=len(ids)
    Gr=G[np.ix_(ids,ids)]
    A=-x[:,None]*Gr
    A.flat[::n+1]+=1
    eA=gamma(6)*(arb(n).sqrt()+arb(float(np.max(x)))*dots.norm(Gr))+8*n*tiny
    lu=lu_factor(np.asfortranarray(A),overwrite_a=True,check_finite=False)
    inverse=lu_solve(lu,np.eye(n,dtype=complex),check_finite=False)
    norm_inv=dots.norm(inverse)
    AX,product_error=dots.product(A,inverse)
    defect=(dots.norm(np.eye(n)-AX)+gamma(5)*(arb(n).sqrt()+dots.norm(AX))
            +product_error+eA*norm_inv+8*n*tiny)
    del inverse,AX
    if not defect<1:
        raise ArithmeticError(f"Bare inverse defect too large ({defect})")
    bound=norm_inv/(1-defect)
    b,eb=scaled_rhs(x,E[ids],dots)
    j,ej,forward=checked_solve(A,lu,eA,bound,b,eb,dots)
    F_raw,eF=product_enclosure(S[:,ids],arb(0),j,ej,dots)
    return enclose_matrix(F_raw,eF),dict(
        active_cells=n,bare_LU=1,verified_inverse_columns=n,bare_physical_rhs=E.shape[1],
        inverse_defect_upper=outward(defect)[1],inverse_Fro_bound=outward(bound)[1],
        bare_transfer_error_Fro_upper=outward(eF)[1],forward=forward,**times(start)),(ids,j,bound)


def tighter_bare_transfer(x,G,E,S,state):
    """Use the same verified inverse bound and compute the bare defect in Arb."""
    started=perf_counter(),process_time()
    ids,j,bound=state
    mat=[arb(float(v)) for v in x[ids]]
    sub=G[np.ix_(ids,ids)]
    n=len(ids)
    A=acb_mat(n,n,[int(i==k)-mat[i]*as_acb(v)
                   for i,row in enumerate(sub) for k,v in enumerate(row)])
    src=E[ids]
    rhs=acb_mat(n,src.shape[1],[mat[i]*as_acb(v) for i,row in enumerate(src) for v in row])
    J=ball(j)
    V=ball(S[:,ids])
    defect=frobenius(rhs-A*J).upper()
    error=(frobenius(V)*bound*defect).upper()
    approx=V*J
    transfer=acb_mat(approx.nrows(),approx.ncols(),
                     [v+acb(arb(0,error),arb(0,error)) for v in approx.entries()])
    return transfer,dict(exact_Arb_residual_Fro_upper=outward(defect)[1],
                         refined_bare_transfer_error_Fro_upper=outward(error)[1],**times(started))


def circuit_models(disk):
    out={}
    for q in COUNTS:
        for omega in STRENGTHS:
            net=network_terms(disk,q,omega)
            M0=ball(net["matrix"])
            I=ball(np.eye(M0.nrows()))
            out[q,int(omega)]=dict(M0=M0,I=I,empty=M0.solve(I),g=[arb(float(g)) for g in net["gain"]],
                                    rx=net["rx"],hold=net["hold"],tx=np.arange(q))
    return out


def circuit_response(transfer,indices,network):
    """Exact loaded impedance solution and reference subtraction in port units."""
    n=len(indices)
    gg=network["g"]
    M=network["M0"]-acb_mat(n,n,[acb(0,gg[i]*gg[j])*transfer[int(a),int(b)]
                                 for i,a in enumerate(indices) for j,b in enumerate(indices)])
    F=2*R*(network["empty"]-M.solve(network["I"]))
    q=len(network["tx"])
    amplitude=arb(q).sqrt()
    def masked(rx):
        return acb_mat(len(rx),q,[F[int(r),int(t)]/amplitude for r in rx for t in network["tx"]])
    return masked(network["rx"]), masked(network["hold"])


def observations(archive,q,omega):
    key=tag(400,q,float(omega),0)
    return dict(identifier=key,Y=ball(archive[key+"_Y"]),
                Y_hold=ball(archive[key+"_Y_hold"]),
                radius=arb(float(archive[key+"_radius"])),
                radius_hold=arb(float(archive[key+"_radius_hold"])))


def trial(transfer,indices,network,observed):
    start=perf_counter(),process_time()
    fit,hold=circuit_response(transfer,indices,network)
    r=frobenius(fit-observed["Y"])/observed["radius"]
    rh=frobenius(hold-observed["Y_hold"])/observed["radius_hold"]
    return dict(train_over_r=outward(r),hold_over_r=outward(rh),train_strict_inside=bool(r<1),
                train_strict_outside=bool(r>1),hold_strict_inside=bool(rh<1),**times(start))


def reference_crosswalk():
    """Establish the identical baseline observation contract for the continuation."""
    base=FIT/"extended_budget/n32_q12/n32_c02.json"
    continued=base.with_name("n32_c02_continuation.json")
    primary=json.loads(base.read_text())
    continuation=json.loads(continued.read_text())
    path=HERE/"public/cases/n32_c02_mimo.npz"
    key=tag(400,12,48,0)
    with np.load(path) as a:
        radius=float(a[key+"_radius"])
        shape=list(a[key+"_Y"].shape)
    proofs=continuation["nominal_endpoint_witnesses"]
    exact_match=(primary["training_observation"]==key and primary["noise_radius"]==radius
                 and all(proof["hard_noise_radius"]==radius for proof in proofs.values()))
    if not exact_match or shape!=[12,12]:
        raise ArithmeticError("The direct continuation and coupling replay do not share raw Y and r")
    return dict(raw_archive=str(path.relative_to(ROOT)),Y_key=key+"_Y",r_key=key+"_radius",
                Y_shape=shape,radius=radius,nominal_source_archive=str(continued.with_suffix(".npz").relative_to(ROOT)),
                nominal_source_keys=["low","high"],same_observation_and_noise_contract=exact_match,
                direct_dressed_proof_upper={k:p["residual_ratio_upper"] for k,p in proofs.items()})


def summarize_grid(rows):
    result={}
    for q in COUNTS:
        for omega in STRENGTHS:
            group=[r for r in rows if r["q"]==q and r["radiation_ohms"]==int(omega)]
            result[f"q{q:02}_c{int(omega):03}"]=dict(count=len(group),q=q,radiation_ohms=int(omega),
                shares_both_extremes=sum(r["exact_full_hull"] for r in group),
                positive_inner_width=sum(r["known_inner_width"] is not None for r in group),
                one_peak_only=sum(r["status"]=="single_peak" for r in group),
                no_verified_peak=sum(r["status"]=="no_verified_peak" for r in group),
                training_and_hold_pass=sum(r["train_and_hold_common_maps"]>0 for r in group),
                inner_widths_by_case={r["case"]:r["known_inner_width"] for r in group})
    return result


def chart(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig,ax=plt.subplots(figsize=(13.0,4.6),layout="constrained")
    layout=[(q,int(s)) for s in STRENGTHS for q in COUNTS]
    look={(r["case"],r["q"],r["radiation_ohms"]):r for r in rows}
    matrix=np.full((len(CASES),len(layout)),np.nan)
    for j,(q,s) in enumerate(layout):
        for i,case in enumerate(CASES):
            row=look[case,q,s]
            if row["known_inner_width"] is not None:
                matrix[i,j]=row["known_inner_width"]
    cmap=plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#e7ebee")
    heat=ax.imshow(matrix,aspect="auto",cmap=cmap,vmin=0,vmax=4.5)
    for j,(q,s) in enumerate(layout):
        for i,case in enumerate(CASES):
            row=look[case,q,s]
            width=row["known_inner_width"]
            value=f"{width:.2f}" if width is not None else ("•" if row["status"]=="single_peak" else "—")
            ax.text(j,i,value,va="center",ha="center",color="white" if width is not None and width>3 else "#13243a",fontsize=10)
    for j in (3.5,7.5):
        ax.axvline(j,color="white",linewidth=3)
    axis=ax.secondary_xaxis("top")
    axis.set_xticks([1.5,5.5,9.5],["16 Ω","48 Ω","144 Ω"],fontweight="bold")
    axis.tick_params(length=0)
    ax.set_xticks(range(len(layout)),[f"{q}×{q}" for q,_ in layout])
    ax.set_yticks(range(len(CASES)),[c[-3:] for c in CASES])
    ax.set_xlabel("400 MHz, 1 W available source power, actual observation at each setting")
    ax.set_ylabel("Frozen blind object and second ROI")
    fig.suptitle("ROI peak width from certified transferred material maps",fontweight="bold",fontsize=15)
    bar=fig.colorbar(heat,ax=ax,pad=.016)
    bar.set_label("Certified lower bound on peak interval width")
    ax.text(.0,-.19,"• One certified peak    — No transferred map verified",transform=ax.transAxes)
    fig.savefig(out,dpi=180)
    plt.close(fig)


def csv_write(path,rows):
    def scalar_row(r):
        out={}
        for k,v in r.items():
            if k in ("low_sources","high_sources","source_paths"):
                out[k]=json.dumps(v,ensure_ascii=False)
            elif isinstance(v,list) and len(v)==2:
                out[k+"_lower"],out[k+"_upper"]=v
            else:
                out[k]=v
        return out
    scaled=[scalar_row(r) for r in rows]
    keys=dict.fromkeys(k for r in scaled for k in r)
    with path.open("w") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(keys))
        writer.writeheader()
        writer.writerows(scaled)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads",type=int,default=4)
    parser.add_argument("--no-plot",action="store_true")
    args=parser.parse_args()
    ctx.prec,ctx.threads=128,1
    start=perf_counter(),process_time()
    disk=DiskGreen(32,400)
    positions,indices=union_positions()
    E,S,_=port_point_matrices(disk,positions)
    G=disk.dense()
    models=circuit_models(disk)
    setup=times(start)
    rows,details,charges=[],[],[]
    with TemporaryDirectory() as tmp:
        dy=Path(tmp)/"explicit_dots.dylib"
        subprocess.run(["clang","-O3","-fno-fast-math","-ffp-contract=off","-dynamiclib",
                        "-o",str(dy),str(HERE/"certify_dot.c")],check=True)
        dots=ExplicitDots(dy,args.threads)
        for case in CASES:
            mask,materials=candidates(case)
            with np.load(HERE/"public/cases"/f"{case}_mimo.npz") as a:
                measured={(q,s):observations(a,q,s) for q,s in models}
            trials=defaultdict(list)
            for i,entry in enumerate(materials):
                transfer,cost,state=bare_certificate(entry["x"],G,E,S,dots)
                identifier=f"{case}_m{i+1:02}"
                refined=None
                refine_cost=None
                for (q,s),model in models.items():
                    acquisition=measured[q,s]
                    initial=trial(transfer,indices[q],model,acquisition)
                    needs_refinement=not (initial["train_strict_inside"] or initial["train_strict_outside"])
                    if needs_refinement and refined is None:
                        refined,refine_cost=tighter_bare_transfer(entry["x"],G,E,S,state)
                        print(case,identifier,"recomputed original bare residual in Arb",refine_cost,flush=True)
                    evidence=trial(refined,indices[q],model,acquisition) if needs_refinement else initial
                    item=dict(case=case,q=q,radiation_ohms=s,observation=acquisition["identifier"],
                            map_id=identifier,roi_peak=entry["peak"],source_paths=entry["sources"],
                            initial_train_over_r=initial["train_over_r"],used_exact_Arb_defect=needs_refinement,
                            initial_small_solve_wall_seconds=initial["wall_seconds"],**evidence)
                    details.append(item)
                    trials[q,s].append(item)
                charges.append(dict(case=case,map_id=identifier,sources=entry["sources"],roi_peak=entry["peak"],
                              exact_Arb_defect_refinement=refine_cost,**cost))
            for (q,s),points in trials.items():
                eligible=[p for p in points if p["train_strict_inside"]]
                peaks=sorted({p["roi_peak"] for p in eligible})
                low=next((p for p in eligible if p["roi_peak"]==peaks[0]),None) if peaks else None
                high=next((p for p in eligible if p["roi_peak"]==peaks[-1]),None) if peaks else None
                row=dict(case=case,roi="roi_02",region_cells=int(mask.sum()),q=q,radiation_ohms=s,
                         observation=measured[q,s]["identifier"],noise_radius=outward(measured[q,s]["radius"]),
                         checked_maps=len(points),strict_feasible_maps=len(eligible),
                         strict_outside_maps=sum(p["train_strict_outside"] for p in points),
                         train_and_hold_common_maps=sum(p["hold_strict_inside"] for p in eligible),
                         known_inner_width=peaks[-1]-peaks[0] if len(peaks)>=2 else None,
                         inner_peak_endpoints=[peaks[0],peaks[-1]] if peaks else None,
                         exact_full_hull=bool(peaks and peaks[0]==0 and peaks[-1]==4.5),
                         status="two_or_more_peaks" if len(peaks)>=2 else "single_peak" if peaks else "no_verified_peak",
                         certified_outer_hull=[0.,4.5],
                         low_witness_id=low["map_id"] if low else None,high_witness_id=high["map_id"] if high else None,
                         low_sources=low["source_paths"] if low else None,high_sources=high["source_paths"] if high else None,
                         lower_of_map_distance_lowers=min(p["train_over_r"][0] for p in points))
                rows.append(row)
            print(case,len(materials),"frozen maps; verified widths",
                  [(q,s,next(r["known_inner_width"] for r in rows if (r["case"],r["q"],r["radiation_ohms"])==(case,q,s)))
                   for (q,s) in models],flush=True)
    output=HERE/"results"
    grid=summarize_grid(rows)
    csv_write(output/"roi_coupling_cells.csv",rows)
    csv_write(output/"roi_coupling_endpoints.csv",details)
    if not args.no_plot:
        chart(rows,HERE/"figures"/"roi_coupling.png")
    (output/"roi_coupling.json").write_text(json.dumps(dict(
        contract="One shared frozen feasible material pool per blind case, second ROI. No truth, held inputs, or refits for testing feasibility.",
        coefficient_scope="Regenerated binary64 bare disk-pulse matrices/point ports and empty impedance, interpreted as exact represented acquisition coefficients",
        included_perturbation_eta=0,frequency_MHz=400,total_available_power_W=1,arb_precision_bits=ctx.prec,
        fitted_case_count=len(CASES),paired_grid_cells=len(rows),endpoint_replays=len(details),
        near_threshold_refinements=sum(p["used_exact_Arb_defect"] for p in details),
        matched_q12_c02_48=reference_crosswalk(),
        bare_LU=sum(c["bare_LU"] for c in charges),
        verified_inverse_columns=sum(c["verified_inverse_columns"] for c in charges),
        bare_physical_rhs=sum(c["bare_physical_rhs"] for c in charges),grid=grid,rows=rows,
        bare_passes=charges,setup=setup,**times(start)),indent=2)+"\n")
    print(grid,flush=True)


if __name__=="__main__":
    main()
