"""Apply trace/interlacing to the existing all-pixel sparse witness and score its center."""
from time import perf_counter,process_time

from flint import arb,ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from outward import clipped_best,f_norm,rowscale,sq,top_sum
from replay import HERE,ball,gershgorin,herm,high,load_json,low,save_json,star,times

CASE="u100_best"
BUDGETS=(26,29,52,59,60,64)


def source_center_upper(j,e,h,proof,weight,mu,cap,greenmax,m,budgets):
    start=perf_counter(),process_time()
    eig,rot=eigh(herm(h),check_finite=False)
    dim=len(eig)-m
    q=rot[:,:dim]
    qb=ball(q)
    qg=star(qb)*qb
    qlo,qhi=gershgorin(qg,upper=False),gershgorin(qg,upper=True)
    if not qlo>0:
        raise ArithmeticError("Source test frame is singular.")
    form=star(qb)*ball(h)*qb
    for i,value in enumerate(eig[:dim]):
        form[i,i]-=arb(float(value))
    source_error=arb(proof["source_form_saved_to_trial_F_error_upper"])+arb(
                 proof["absolute_gap_defect_correction_upper"])
    shift=(f_norm(form)+source_error*qhi).upper()
    base=sum((max(arb(0),(arb(float(value))+shift).upper()) for value in eig[:dim]),arb(0))
    jp,ep=j@q,e@q
    material=clipped_best(jp,ep,float(cap))
    row_a=rowscale(ball(e)*qb,[arb(float(x)) for x in material])-ball(j)*qb
    row_power=np.array([high(mu*arb(float(w))*sum((sq(row_a[t,l]) for l in range(dim)),arb(0)))
                      for t,w in enumerate(weight)])
    jerr0=arb(proof["center_error_frobenius_upper"])
    jerr=jerr0+arb(proof["center_saved_to_trial_F_error_upper"])
    eerr=greenmax*jerr0+arb(proof["total_center_field_saved_to_trial_F_error_upper"])
    drift=(mu.sqrt()*arb(float(max(weight))).sqrt()*qhi.sqrt()*(jerr+cap*eerr)).upper()
    radius=arb(proof["radius"])
    out=[]
    for k in budgets:
        center_bound=(top_sum(row_power,k).sqrt()+drift)**2
        distance=((base+center_bound)/qlo).sqrt()
        out.append(dict(id=proof["id"],row_budget=k,free_negative_directions=m,
                        source_spectral_term_upper=high(base/qlo),signed_center_term_upper=high(center_bound/qlo),
                        ratio_upper=high(distance/radius),
                        cannot_exclude=bool(distance<radius),**times(*start)))
    return out,material,row_power


def main():
    start=perf_counter(),process_time()
    raw=np.load(HERE/f"{CASE}.npz",allow_pickle=False)
    strict=load_json(HERE/f"{CASE}_strict.json")
    with np.load(HERE/f"{CASE}_strict.npz",allow_pickle=False) as previous:
        weights=previous["frame_weight_upper"]
        full_witness=previous["contraction_row_power_upper_13"]
    state=strict["verified_state"]
    nu=state["exact_nu"]
    test=next(r for r in strict["contractions"] if r["contraction_dimension"]==nu)
    dn=arb(test["negative_energy_lower"])
    drift=arb(test["norm_projection_error_upper"])
    counts=[]
    for k in BUDGETS:
        bound=high((top_sum(full_witness,k).sqrt()+drift)**2/dn)
        counts.append(dict(rows=k,trace_upper=bound,maximum_compensated_dimension=min(nu,int(np.floor(bound)))))
    meta=load_json(HERE/f"{CASE}.json")
    cap=arb(meta["profile"]["cap"])
    mu=arb(meta["profile"]["scale"])
    green=arb(strict["projector"]["green_norm_upper"])
    c0=strict.get("outward_centers")
    if c0 is None:
        c0=[dict(**r["outward_baseline"],id=r["id"],radius=r["padded_radius"]) for r in meta["rows"]]
    ans,archive=[],{}
    for c in sorted({r["maximum_compensated_dimension"] for r in counts}):
        budgets=[r["rows"] for r in counts if r["maximum_compensated_dimension"]==c]
        for i,(j,e,h,proof) in enumerate(zip(raw["center"],raw["center_field"],raw["hs"],c0,strict=True)):
            vals,x,d=source_center_upper(j,e,h,proof,weights,mu,cap,green,nu-c,budgets)
            ans.extend(vals)
            archive[f"center_x_c_{c}_object_{i+1}"]=x
            archive[f"center_cost_c_{c}_object_{i+1}"]=d
        print(c,budgets, [(k,max(x["ratio_upper"] for x in ans if x["row_budget"]==k)) for k in budgets],
              flush=True)
    save_json(HERE/"deflated_sparse_tightening.json",dict(case=CASE,
        scope="Any diagonal nonnegative frame with the stated support cardinality, no positive deflation, arbitrary signed center matrix under the same source containment",
        all_pixel_sorted_witness=True, all_pixel_count=len(full_witness),contraction_dimension=nu,
        method="At least c eigenvalues of the compressed full negative Gram would be at least 1; its upper trace gives c <= floor(trace_upper)",
        counts=counts,source_upper_bounds=ans,operations=dict(no_additional_Green_columns=True),total=times(*start)))
    np.savez_compressed(HERE/"deflated_sparse_tightening.npz",**archive)


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
