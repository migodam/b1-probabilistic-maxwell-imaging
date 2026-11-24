"""Verify an all-pixel common negative direction from a positive-frame LP dual."""
from time import perf_counter,process_time

from flint import acb_mat, arb, ctx
import numpy as np
from scipy.optimize import linprog
from threadpoolctl import threadpool_limits

from outward import clipped_best,f_norm,rowscale
from replay import HERE, ball, configure, high, load_json, low, save_json, star, times


def main():
    start=perf_counter(),process_time()
    case="u100_best"
    geo,native,_,_=configure(case)
    strict=load_json(HERE/f"{case}_strict.json")
    archive=np.load(HERE/f"{case}_strict.npz",allow_pickle=False)
    weights=archive["frame_weight_upper"]
    state=strict["projector"]
    leak=arb(state["ritz_positive_leak_upper"])
    umin=arb(state["negative_projected_gram_min_lower"])
    kneg=arb(state["negative_K_norm_upper"])
    cap=arb(native.w.cap)
    op_row_cap=1+cap*arb(state["green_norm_upper"])
    src=np.load(HERE/f"{case}_conditioned.npz",allow_pickle=False)
    trial=np.load(HERE/f"{case}_all_pixel_probe.npz",allow_pickle=False)
    eig=src["negative_eig"]
    v=src["negative_vec"]
    z=v[:,0:1]/np.sqrt(-eig[0])
    gz=src["negative_gn"][:,0:1]
    x=clipped_best(z,gz,native.w.cap)
    obj=abs(z[:,0]-x*gz[:,0])**2
    d=trial["frame_densities"]
    lp=linprog(-obj,A_ub=d,b_ub=np.ones(len(d)),
               bounds=np.column_stack((np.zeros(len(weights)),weights)),method="highs")
    if not lp.success:
        raise ArithmeticError(lp.message)
    dual=np.maximum(0.,-lp.ineqlin.marginals)
    keep=dual>1.e-12
    dual=dual[keep]
    p,gp=trial["positive_trial"][:,keep],trial["positive_trial_green"][:,keep]

    zb,vb=ball(z),ball(v)
    az=zb-rowscale(geo.exact_green(zb),[arb(float(i)) for i in x])
    # Use the represented z in the exact split, and ensure the tiny root
    # representation error is also covered by computing its actual V scalar.
    root=arb(float(1/np.sqrt(-eig[0])))
    vz=acb_mat(len(z),1,[vb[j,0] for j in range(len(z))])
    z_leak=(leak*root+f_norm(zb-vz*root)).upper()
    negative_energy=-(star(zb)*native.exact_action(zb))[0,0].real
    negative_density=np.array([high((abs(az[i,0])+op_row_cap*z_leak)**2) for i in range(len(z))])

    pb=ball(p)
    pg=geo.exact_green(pb)
    ap=native.exact_action(pb)
    out_positives=[]
    positive_density=np.empty((len(dual),len(z)))
    positive_norm=[]
    for i,beta in enumerate(dual):
        this=p[:,i:i+1]
        this_b=ball(this)
        nerr=((f_norm(star(vb)*this_b)+leak*f_norm(this_b))/umin.sqrt()).upper()
        h=(star(this_b)*acb_mat(len(z),1,[ap[j,i] for j in range(len(z))]))[0,0].real
        energy=(h+kneg*nerr**2).upper()
        positive_norm.append(high(energy))
        endpoint=abs(p[:,i]-native.w.cap*gp[:,i])>abs(p[:,i])
        shift=op_row_cap*nerr
        for j,last in enumerate(endpoint):
            value=pb[j,i]-cap*pg[j,i] if last else pb[j,i]
            amplitude=max(arb(0),(abs(value)-shift).lower())
            positive_density[i,j]=max(0.,low(amplitude**2))
        out_positives.append(dict(weight=float(beta),positive_energy_upper=high(energy),
                                  trial_negative_leak_upper=high(nerr)))
    money=sum((arb(float(beta))*arb(float(e)) for beta,e in zip(dual,positive_norm,strict=True)),arb(0))
    rest=arb(0)
    for j,f in enumerate(weights):
        dominated=sum((arb(float(a))*arb(float(b)) for a,b in zip(
            dual,positive_density[:,j],strict=True)),arb(0))
        slack=max(arb(0),(arb(float(negative_density[j]))-dominated).upper())
        rest+=arb(float(f))*slack
    upper=(money+rest).upper()
    margin=(negative_energy-upper).lower()
    result=dict(case=case,verified_state=strict["verified_state"],scope_pixels=len(weights),
        full_compensation_impossible=bool(margin>0),strict_margin_lower=low(margin),
        energy_lower=low(negative_energy),weighted_negative_trace_upper=high(upper),
        maximum_compensated_dimension=len(eig)-1 if margin>0 else len(eig),
        frame_dual_term_upper=high(money),positive_row_price_mismatch_upper=high(rest),
        positive_witnesses=out_positives,
        numerical_lp_upper_proposal=-float(lp.fun),operations=dict(geo.cost),**times(*start))
    out=HERE/f"{case}_all_pixel_strict"
    save_json(out.with_suffix(".json"),result)
    np.savez_compressed(out.with_suffix(".npz"),negative_current=z,negative_witness_x=x,
        negative_density_upper=negative_density,positive_density_lower=positive_density,
        positive_witness_ids=np.flatnonzero(keep),dual=dual)
    print(result,flush=True)


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
