"""One native gap refinement by lifting a small Ritz span.

The added receiving rows are a mathematical device. They are never used as
synthetic observations or passed into the final source response.
"""

import argparse
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from screen import CAPS, HERE, RECEIVERS, SOURCES, inputs, load_json, save_json, times
from verify_native import (FourierEnvelope, Geometry, NativeState, Weights, analyze_row,
                           ball, digest_profile, frobenius_squared, gershgorin, high, low,
                           reusable_stationarity, star)


def lifted_gap(native, weights, neg):
    t,c=perf_counter(),process_time()
    decomposition=native.modes(neg+9)
    if native.nu!=neg:
        raise ArithmeticError("The refined Ritz proposal disagrees with the proven negative count.")
    geo=native.g
    added=np.sqrt(weights.scale*8.)*native.eig_basis.conj().T
    lifted=Geometry(geo.spec,np.vstack((geo.sensor,added)),geo.incident)
    env=FourierEnvelope(lifted)
    proposal=env.verify(weights,0)
    v=ball(added[neg:].conj().T)
    av=native.exact_action(v)
    d=v*0
    for j in range(v.ncols()):
        value=arb(float(native.eigs[neg+j]))
        for i in range(v.nrows()):
            d[i,j]=v[i,j]*value
    defect=frobenius_squared(av-d).upper()
    gram=star(v)*v
    tform=star(v)*av
    g0=arb(proposal["positive_gap_lower"])
    if not gershgorin(gram,upper=False)>0:
        raise ArithmeticError("The represented positive Ritz columns were dependent.")
    for gtrial in (2**-6,2**-8,2**-10,2**-12,2**-15,2**-18):
        g=arb(gtrial)
        if g>=g0:
            continue
        margin=gershgorin(tform-g*gram,upper=False)-defect/(g0-g)
        if margin>0:
            return dict(exact_nu=neg,positive_gap_lower=gtrial,
                        lifted_fourier_proof=proposal,verified_lift_width=len(added),
                        positive_ritz_schur_margin_lower=low(margin),
                        positive_ritz_defect_upper=high(defect.sqrt()),
                        refined_numeric_ritz=decomposition,**times(t,c))
    raise ArithmeticError("The native lifted complement and positive Ritz span did not separate.")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receivers",type=int,choices=(0,*RECEIVERS),required=True)
    args=parser.parse_args()
    p=args.receivers
    ctx.prec,ctx.threads=112,2
    started,cpu=perf_counter(),process_time()
    spec,sensors,sources,observations=inputs()
    if p:
        operators=sensors[p]
        old=load_json(HERE/f"native_p{p}.json")
        choose=[r for r in old["profiles"] if r["name"]=="stable_scalar_optical" and
                "exact_nu" not in r["state"] and r["state"]["native_nu_at_least"]<32]
        cases=[(r["cap"],digest_profile(r["profile"],r["cap"]),r["state"],r["name"],r["profile"])
                for r in choose]
    else:
        operators=np.zeros((1,len(sources[SOURCES[0]])),complex)
        old=load_json(HERE/"bare_inertia.json")
        cases=[(r["cap"],Weights.global_(r["cap"],1.,0.),r,"bare_real_sector",{})
               for r in old["rows"] if "exact_nu" not in r]
    g=Geometry(spec,operators,sources[SOURCES[0]])
    out=dict(scope="Rank lift used only to verify the native state complement; physical measurements unchanged",
             receiver_count=p,precision_bits=ctx.prec,profiles=[])
    for u,w,prior,name,profile in cases:
        t,c=perf_counter(),process_time()
        pre=g.cost.copy()
        new=dict(name=name,profile=profile,cap=u,state=dict(prior),observations=[])
        try:
            native=NativeState(g,w)
            new["state"].update(lifted_gap(native,w,prior["native_nu_at_least"]))
            s_use=[s for s in SOURCES if s>new["state"]["exact_nu"]]
            if s_use and p:
                stat,gram,indices,sol=reusable_stationarity(native,{s:sources[s] for s in s_use})
                new["shared_source_response"]=sol
                for s in s_use:
                    for x,y in observations[p][s]:
                        r=dict(id=x["id"],object_id=x["object_id"],sources=s,receivers=p,noise_pct=x["noise_pct"],
                               radius=x["padded_radius"])
                        ytest=analyze_row(y,new["state"],w.scale,stat,gram,indices[s],r["radius"])
                        new["observations"].append(dict(**r,**ytest))
        except (ArithmeticError,RuntimeError,ValueError) as err:
            new["refinement_error"]=str(err)
        new.update(**times(t,c),operations=dict(g.cost-pre))
        out["profiles"].append(new)
        out.update(times(started,cpu))
        save_json(HERE/f"gap_refinement_p{p}.json",out)
        print(p,u,new["state"].get("exact_nu"),new.get("refinement_error"),
              new["wall_seconds"],flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
