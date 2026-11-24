"""Numerically test source deep-negative absorption, without using it as a proof."""
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from outward import clipped_best
from replay import HERE,herm,load_json,save_json,times

CASE="u100_best"
BUDGETS=(60,64)
PENALTIES=(0.,1.e-5,1.e-4,.001,.01,.1,1.,10.)


def sq_rows(a):
    return np.sum(abs(a)**2,axis=1)


def highest_rows(r,k):
    return float(np.partition(r,-k)[-k:].sum())


def main():
    setup=perf_counter(),process_time()
    raw=np.load(HERE/f"{CASE}.npz",allow_pickle=False)
    w=np.load(HERE/f"{CASE}_strict.npz",allow_pickle=False)["frame_weight_upper"]
    run=load_json(HERE/f"{CASE}.json")
    check=load_json(HERE/"deflated_sparse_tightening.json")
    mu,cap=run["profile"]["scale"],run["profile"]["cap"]
    nu=run["verified_state"]["exact_nu"]
    result=[]
    for meta,j,e,h in zip(run["rows"],raw["center"],raw["center_field"],raw["hs"],strict=True):
        eigen,v=eigh(herm(h),check_finite=False)
        rows=[]
        for k in BUDGETS:
            c=next(r["maximum_compensated_dimension"] for r in check["counts"] if r["rows"]==k)
            dim=len(eigen)-nu+c
            eig,q=eigen[:dim],v[:,:dim]
            r=meta["padded_radius"]
            best=None
            ntry=0
            for d in range(1,dim):
                if eig[d-1]>=0:
                    break
                deep=q[:,:d]/np.sqrt(-eig[:d])
                weak=q[:,d:]
                jd,ed=j@deep,e@deep
                jw,ew=j@weak,e@weak
                hweak=np.maximum(eig[d:],0.).sum()
                for a in PENALTIES:
                    xx=clipped_best(np.column_stack((jw,np.sqrt(a)*jd)),
                                    np.column_stack((ew,np.sqrt(a)*ed)),cap)
                    deep_cost=highest_rows(mu*w*sq_rows(jd-xx[:,None]*ed),k)
                    weak_cost=highest_rows(mu*w*sq_rows(jw-xx[:,None]*ew),k)
                    ntry+=1
                    if deep_cost>=1:
                        continue
                    upper=float(np.sqrt(hweak+weak_cost/(1-deep_cost))/r)
                    if best is None or upper<best["ratio_upper_numeric"]:
                        best=dict(ratio_upper_numeric=upper,theta_upper_numeric=deep_cost,
                                  weak_center_term_numeric=weak_cost,deep_modes=d,proposal_tradeoff=a)
            baseline=next(t["ratio_upper"] for t in check["source_upper_bounds"]
                          if t["row_budget"]==k and t["id"]==meta["id"])
            rows.append(dict(row_budget=k,source_compression_dimension=dim,checked_proposals=ntry,
                             best_numeric_split=best,previous_outward_unsplit_upper=baseline))
        result.append(dict(id=meta["id"],rows=rows))
        print(result[-1],flush=True)
    save_json(HERE/"deflated_deep_probe.json",dict(scope="Numerical diagnostic only; no outward split-source bound used",
        source_rule="When Ddeep*Adeep^(-1/2) has squared norm <= theta<1, Young's inequality absorbs the deep block and charges the weak center by 1/(1-theta)",
        results=result,operations=dict(no_additional_Green_columns=True),total=times(*setup)))


if __name__=="__main__":
    with threadpool_limits(1):
        main()
