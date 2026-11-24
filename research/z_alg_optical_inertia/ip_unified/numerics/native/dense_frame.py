"""Numerical partial-compensation screen on 208, 416 and 9,216 pixels."""
from time import perf_counter, process_time

import numpy as np
from scipy.optimize import linprog
from threadpoolctl import threadpool_limits

from outward import clipped_best
from replay import HERE, load_json, save_json, times


def run_upper(score,density,weight):
    opt=linprog(-score,A_ub=density,b_ub=np.ones(len(density)),
                bounds=np.column_stack((np.zeros(len(weight)),weight)),method="highs")
    if not opt.success:
        raise ArithmeticError(opt.message)
    return -float(opt.fun),opt.x,-opt.ineqlin.marginals


def main():
    start=perf_counter(),process_time()
    label="u100_best"
    vals=load_json(HERE/f"{label}.json")
    cap=vals["profile"]["cap"]
    source=np.load(HERE/f"{label}_conditioned.npz",allow_pickle=False)
    pos=np.load(HERE/f"{label}_all_pixel_probe.npz",allow_pickle=False)
    wc=np.load(HERE/f"{label}_strict.npz",allow_pickle=False)["frame_weight_upper"]
    nu=len(source["negative_eig"])
    n=source["negative_vec"]/np.sqrt(-source["negative_eig"])
    g=source["negative_gn"]
    more=12
    rng=np.random.default_rng(26011250)
    rot=rng.standard_normal((pos["positive_trial"].shape[1],more))+1j*rng.standard_normal(
                                  (pos["positive_trial"].shape[1],more))
    rot/=np.linalg.norm(rot,axis=0)
    pp,pg=pos["positive_trial"]@rot,pos["positive_trial_green"]@rot
    density=np.concatenate((pos["frame_densities"],np.maximum(abs(pp)**2,abs(pp-cap*pg)**2).T))
    directions=np.column_stack((np.eye(nu),source["homogeneous_subspace"]))
    a,b=n@directions,g@directions
    material=np.clip(np.real(a.conj()*b)/np.maximum(abs(b)**2,np.finfo(float).tiny),0.,cap)
    score=abs(a-material*b)**2
    priority=np.max(score,axis=1)*wc
    idx=source["pixel_index"].tolist()
    initial=set(idx)
    idx.extend(int(k) for k in np.argsort(priority)[::-1] if k not in initial)
    idx=np.array(idx)
    proof={}
    archive={}
    for budget in (16*nu,32*nu,len(n)):
        t,c=perf_counter(),process_time()
        pixel=idx[:budget]
        d=density[:,pixel]
        w=wc[pixel]
        attempted=[run_upper(q[pixel],d,w) for q in score.T]
        i=int(np.argmax([r[0] for r in attempted]))
        best,best_weight,ddual=attempted[i]
        traces=[]
        for basis_index,basis in enumerate([np.eye(nu),source["homogeneous_subspace"]]):
            for dim in range(1,nu+1):
                a,b=n@basis[:,:dim],g@basis[:,:dim]
                x=clipped_best(a,b,cap)
                trace=np.sum(abs(a-x[:,None]*b)**2,axis=1)
                upper,weight,lpdual=run_upper(trace[pixel],d,w)
                traces.append(dict(dimension=dim,basis=basis_index,
                    sampled_lp_trace_upper=upper,
                    compensated_dimension_upper_numeric=min(nu,nu-dim+int(np.floor(upper)))))
                archive[f"trace_dual_{budget}_{basis_index}_{dim}"]=lpdual
        proof[str(budget)]=dict(best_direction=i,best_necessary_capacity_upper=best,
            all_direction_upper=[float(a[0]) for a in attempted],trace_bounds=traces,
            favorable_comp_dimension_upper_numeric=min(r["compensated_dimension_upper_numeric"] for r in traces),
            **times(t,c))
        archive[f"pixels_{budget}"]=pixel
        archive[f"candidate_best_weight_{budget}"]=best_weight
        archive[f"candidate_dual_{budget}"]=ddual
        print(label,budget,proof[str(budget)]["best_necessary_capacity_upper"],
              proof[str(budget)]["favorable_comp_dimension_upper_numeric"],flush=True)
    archive["trial_frame_densities"]=density
    archive["directions"]=directions
    archive["row_priority"]=priority
    out=HERE/f"{label}_dense_frame"
    save_json(out.with_suffix(".json"),dict(case=label,scope="Necessary screen using sampled positive frame inequalities; candidate LP weights need not satisfy the full frame. All results herein numerical.",
        data=proof,positive_test_columns=len(density),**times(*start)))
    np.savez_compressed(out.with_suffix(".npz"),**archive)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
