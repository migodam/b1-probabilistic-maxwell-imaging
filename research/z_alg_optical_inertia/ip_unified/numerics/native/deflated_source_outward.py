"""Outward source-center caps for the uniform-frame deflation experiment."""
from time import perf_counter,process_time

from flint import arb,ctx
import numpy as np
from threadpoolctl import threadpool_limits

from deflated_sparse_tightening import source_center_upper
from replay import HERE,high,load_json,save_json,times

CASE="u100_best"


def main():
    t,c0=perf_counter(),process_time()
    raw=np.load(HERE/f"{CASE}.npz",allow_pickle=False)
    strict=load_json(HERE/f"{CASE}_strict.json")
    frame=load_json(HERE/"deflated.json")
    state=strict["verified_state"]
    nu=state["exact_nu"]
    meta=load_json(HERE/f"{CASE}.json")
    cap=arb(meta["profile"]["cap"])
    mu=arb(meta["profile"]["scale"])
    g_norm=arb(strict["projector"]["green_norm_upper"])
    centers=strict.get("outward_centers")
    if centers is None:
        centers=[dict(**r["outward_baseline"],id=r["id"],radius=r["padded_radius"]) for r in meta["rows"]]
    c_list=sorted({r[f]["maximum_compensated_dimension"] for r in frame["results"] for f in ("all_box","young")})
    tests={}
    for c in c_list:
        tests[c]=[]
        for j,e,h,proof in zip(raw["center"],raw["center_field"],raw["hs"],centers,strict=True):
            seq,_,_=source_center_upper(j,e,h,proof,np.ones(len(j)),mu,cap,g_norm,nu-c,[len(j)])
            tests[c].append(seq[0])
    results=[]
    for r in frame["results"]:
        q=r["weak_positive_modes"]
        for scope in ("all_box","young"):
            proof=r[scope]
            c=proof["maximum_compensated_dimension"]
            w=arb(proof["frame"]["weight_upper"])
            six=[]
            for m,s in zip(meta["rows"],tests[c],strict=True):
                d=(arb(s["source_spectral_term_upper"])+w*arb(s["signed_center_term_upper"])).sqrt()/arb(
                        m["padded_radius"])
                six.append(dict(id=m["id"],ratio_upper=high(d),strict_no_crossing=bool(d<1)))
            results.append(dict(weak_positive_modes=q,scope=scope,
                compensated_dimension_upper=c,signed_center_ratio_upper=six))
            print(q,scope,c,max(r["ratio_upper"] for r in six),sum(r["strict_no_crossing"] for r in six),
                  flush=True)
    save_json(HERE/"deflated_source_outward.json",dict(scope="Requires the augmented block compression to include all pure-source coordinates. Uniform nonnegative frame weights, arbitrary signed center, native data.",
                      effective_center_source_tests=tests,results=results,
                      operations=dict(no_additional_Green_columns=True),total=times(t,c0)))


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
