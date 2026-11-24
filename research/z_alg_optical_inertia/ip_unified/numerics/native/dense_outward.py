"""Verify the selected dense positive-frame trace duals with native balls."""
from time import perf_counter,process_time

from flint import acb_mat,arb,ctx
import numpy as np
from threadpoolctl import threadpool_limits

from outward import clipped_best,f_norm,rowscale,sq
from replay import (HERE,ball,configure,gershgorin,high,load_json,low,save_json,star,times)


def main():
    start=perf_counter(),process_time()
    case="u100_best"
    geo,native,_,_=configure(case)
    strict=load_json(HERE/f"{case}_strict.json")
    projector=strict["projector"]
    trace=load_json(HERE/f"{case}_dense_frame.json")
    cache=np.load(HERE/f"{case}_dense_frame.npz",allow_pickle=False)
    spec=np.load(HERE/f"{case}_conditioned.npz",allow_pickle=False)
    positive=np.load(HERE/f"{case}_all_pixel_probe.npz",allow_pickle=False)
    w=np.load(HERE/f"{case}_strict.npz",allow_pickle=False)["frame_weight_upper"]
    v,eig=spec["negative_vec"],spec["negative_eig"]
    v0=ball(v)
    green_v=geo.exact_green(v0)
    av=native.exact_action(v0)
    vk=star(v0)*av
    n,g=v/np.sqrt(-eig),spec["negative_gn"]
    cap=arb(native.w.cap)
    leak=arb(projector["ritz_positive_leak_upper"])
    umin=arb(projector["negative_projected_gram_min_lower"])
    maxneg=arb(projector["negative_K_norm_upper"])
    opnorm=1+cap*arb(projector["green_norm_upper"])
    pos_trial,pos_green=positive["positive_trial"],positive["positive_trial_green"]
    traced,archive={},{}
    for budget,rr in trace["data"].items():
        timer=perf_counter(),process_time()
        best=min(rr["trace_bounds"], key=lambda z:(len(eig)-z["dimension"]+z["sampled_lp_trace_upper"]))
        dim,mode=best["dimension"],best["basis"]
        basis=(np.eye(len(eig)) if mode==0 else spec["homogeneous_subspace"])[:,:dim]
        ncoord=ball(basis/np.sqrt(-eig)[:,None])
        j,gj=v0*ncoord,green_v*ncoord
        material=clipped_best(n@basis,g@basis,native.w.cap)
        x=j-rowscale(gj,[arb(float(x)) for x in material])
        jdrift=(opnorm*leak*f_norm(ncoord)).upper()
        energy=-gershgorin(star(ncoord)*vk*ncoord,upper=True)
        rows=cache[f"pixels_{budget}"]
        desired=np.array([high((sum((sq(x[q,k]) for k in range(dim)),arb(0)).sqrt()+jdrift)**2)
                          for q in rows])
        y=cache[f"trace_dual_{budget}_{mode}_{dim}"]
        active=np.flatnonzero(y>1.e-12)
        if np.any(active>=pos_trial.shape[1]):
            raise ArithmeticError("A chosen trace dual uses a sampled witness with no archived native vector.")
        p=pos_trial[:,active]
        p0=ball(p)
        gp=geo.exact_green(p0)
        kp=native.exact_action(p0)
        all_usage=np.empty((len(active),len(rows)))
        right=[]
        for k,wi in enumerate(active):
            r=acb_mat(len(v),1,[p0[x,k] for x in range(len(v))])
            resid=((f_norm(star(v0)*r)+leak*f_norm(r))/umin.sqrt()).upper()
            qa=(star(p0)*kp)[k,k].real
            right.append(high(qa+maxneg*resid**2))
            strong=np.abs(p[:,k]-native.w.cap*pos_green[:,wi])>np.abs(p[:,k])
            for l,i in enumerate(rows):
                q=p0[i,k]-cap*gp[i,k] if strong[i] else p0[i,k]
                d=max(arb(0),(abs(q)-opnorm*resid).lower())
                all_usage[k,l]=max(0.,low(d**2))
        money=sum((arb(float(y[i]))*arb(float(rhs)) for i,rhs in zip(active,right,strict=True)),arb(0))
        remainder=arb(0)
        for idx,i in enumerate(rows):
            service=sum((arb(float(y[j]))*arb(float(pp))
                         for j,pp in zip(active,all_usage[:,idx],strict=True)),arb(0))
            missing=max(arb(0),(arb(float(desired[idx]))-service).upper())
            remainder+=arb(float(w[i]))*missing
        strict_trace=high((money+remainder)/energy)
        c=min(len(eig),len(eig)-dim+int(np.floor(strict_trace)))
        out=dict(row_count=len(rows),basis_number=mode,dimension=dim,energy_lower=low(energy),
                 normalized_trace_upper=strict_trace,positive_dual_usage_upper=high(money),
                 rejected_price_upper=high(remainder),strict_compensated_dimension_upper=c,
                 candidate_best_necessary_capacity_upper_numeric=rr["best_necessary_capacity_upper"],
                 active_positive_witnesses=active.tolist(),**times(*timer))
        traced[budget]=out
        archive[f"budget_{budget}_negative_test_x"]=material
        archive[f"budget_{budget}_negative_density_upper"]=desired
        archive[f"budget_{budget}_positive_density_lower"]=all_usage
        archive[f"budget_{budget}_positive_energy_upper"]=right
        archive[f"budget_{budget}_negative_basis"]=basis
        print(case,out,flush=True)
    out=HERE/f"{case}_dense_outward"
    save_json(out.with_suffix(".json"),dict(case=case,data=traced,operations=dict(geo.cost),**times(*start)))
    np.savez_compressed(out.with_suffix(".npz"),**archive)


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
