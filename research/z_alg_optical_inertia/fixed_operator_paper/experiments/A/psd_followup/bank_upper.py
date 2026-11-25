"""Check whether any of the fixed bank's positive-state tails can clear the 90 remaining gaps."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from common import PUBLIC, RESULTS, load, write_json
from bounds import (Acquisition, arb, assemble, ball, bt, ctx, frob, high, herm, low,
                    propose_frame, select, validate_frame, weight_vector)
from replay import original
from run import PARAMS
import numpy as np
from scipy.linalg import eigh

ctx.prec=192
OUT=HERE/"positive_tail_upper_audit"


def compare(acq,ell,upp):
    results=[]
    n,s=acq.n,acq.s
    state,source=range(n),range(n,n+s)
    for key in PARAMS:
        exact=original(acq,ell,upp,*weight_vector(acq.n,*key))
        K=select(exact,state,state)
        numerical=assemble(acq,ell,upp,*key)
        eigs,U=eigh(numerical[:n,:n],check_finite=False)
        min_eig=float(eigs[0])
        entry={"weights":key,"smallest_K_eigenvalue_proposal":min_eig}
        if min_eig<0.:
            u=ball(U[:,:1])
            value=(bt(u)*K*u)[0,0].real
            result=dict(K_has_verified_negative_direction=bool(value<0),
                        negative_direction_upper=high(value),**entry)
            results.append(result)
            continue
        F,edge=propose_frame(numerical[:n,:n],n)
        if F is None:
            results.append(dict(error="positive_K_frame_missing",**entry))
            continue
        pd,krecord=validate_frame(K,F)
        if not pd>0.:
            results.append(dict(error="positive_K_unverified",**entry))
            continue
        L=select(exact,state,source)
        H=select(exact,source,source)-bt(L)*K.solve(L)
        midpoint=np.array([complex(v) for v in H.entries()]).reshape(s,s)
        d,V=eigh(herm(midpoint),check_finite=False)
        V=ball(V)
        E_g=frob(bt(V)*V-ball(np.eye(s))).upper()
        E_h=frob(bt(V)*H*V-ball(np.diag(d))).upper()
        if not E_g<1.:
            results.append(dict(error="source_frame_rank_unverified",**entry))
            continue
        upper=sum((max(arb(0),(arb(float(eig))+E_h).upper())/(1-E_g)
                   for eig in d),arb(0)).upper()
        results.append(dict(K_PSD=True,K_frame_gap=krecord["stretched_positive_gap_lower"],
                            H_Gram_error_upper=high(E_g),H_diagonalization_error_upper=high(E_h),
                            strictly_below_noise=bool(upper<arb(acq.radius)**2),
                            positive_tail_upper_over_r2=high(upper/arb(acq.radius)**2),**entry))
    return results


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--begin",type=int,default=1)
    parser.add_argument("--end",type=int,default=24)
    parser.add_argument("--summary",action="store_true")
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    if args.summary:
        all_items=[]
        total_time=0.
        for path in sorted(OUT.glob("pair_*.json")):
            receipt=json.loads(path.read_text())
            all_items.extend(receipt["rows"])
            total_time+=receipt["seconds"]
        assert len(all_items)==90
        n_total=0
        counts=Counter()
        positive_upp=[]
        unresolved=[]
        for x in all_items:
            for r in x["weights"]:
                n_total+=1
                if r.get("K_has_verified_negative_direction",False):
                    counts["verified_indefinite_or_negative"]+=1
                elif r.get("strictly_below_noise",False):
                    counts["verified_PSD_positive_trace_below_noise"]+=1
                    positive_upp.append(r["positive_tail_upper_over_r2"])
                else:
                    unresolved.append((x["pair_id"],x["node"],r))
        write_json(HERE/"positive_tail_upper_audit.json",dict(
            scope="All 72 same-bank weights on each of the 90 residual primary C-explicit and strengthened B-unresolved boxes.",
            boxes=len(all_items),weights=n_total,proved=dict(counts),
            verified_PSD_tail_upper_over_r2_max=max(positive_upp) if positive_upp else None,
            unresolved=unresolved,precision_bits=ctx.prec,wall_seconds_sum=total_time,
            interpretation="Every H positive-tail bound that missed the threshold is only for the frozen finite bank. Globally feasible continuous PSD multipliers are outside this conclusion."))
        return
    data=load(PUBLIC)
    for k in range(args.begin,args.end+1):
        path=OUT/f"pair_{k:02d}.json"
        if path.exists():
            continue
        rows=json.loads((HERE/"positive_tail"/path.name).read_text())["rows"]
        started=perf_counter()
        chosen=[]
        for row in rows:
            if row["C"]!="EXCLUDED" or row["status"]!="UNRESOLVED":
                continue
            t,b,w,box=row["node"]
            ell,upp=data["boxes"][k-1,b,w,box]
            acq=Acquisition(data["G"],data["S"],data["B"],data["Y"][k-1,t,b],
                            float(data["radii"][k-1,t,b]))
            chosen.append(dict(pair_id=row["pair_id"],node=row["node"],weights=compare(acq,ell,upp)))
        write_json(path,dict(pair=k,rows=chosen,seconds=perf_counter()-started))
        print(k,len(chosen),flush=True)


if __name__=="__main__":
    main()
