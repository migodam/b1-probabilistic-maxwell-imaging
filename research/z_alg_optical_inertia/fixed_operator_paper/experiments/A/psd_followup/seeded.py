"""Rescreen Phase-A B's original seeded search after exact vacuum elimination."""
from collections import Counter
from time import perf_counter

from solve import (Acquisition, ArrayStore, HERE, arb, eigh, forms, hadamard, np, project,
                   propose_frame, save, select_public, strict_schur, write_json, low, independent)
from baseline_psd import potential_plans, summarize
from bounds import weight_vector
from run import PARAMS

STAGE_MAP={"matrix_qcqp_hadamard":"H","source_diagonal_qcqp":"I"}


def one(acq,ell,upp,keep,frames):
    start=perf_counter()
    a0,W,L=forms(acq.G,acq.S,acq.B,acq.Y,ell,upp)
    modes=(np.eye(acq.s),hadamard(acq.s))
    trials={}
    for params in PARAMS:
        w,o=[v[keep] for v in weight_vector(len(keep),*params)]
        key=tuple(w)+tuple(o)
        if key in trials:
            continue
        value=summarize(acq,a0+np.einsum("k,kij->ij",w,W)+np.einsum("k,kij->ij",o,L),w,o,modes)
        if value is not None:
            trials[key]=value
    initial=list(trials.values())
    initial_count=len(initial)
    for bank in (initial,[r for r in initial if not np.any(r["optical"])]):
        if not bank:
            continue
        best=max(bank,key=lambda r:r["trace"])
        w,o=best["weight"],best["optical"]
        variants=[]
        for i in range(acq.n):
            for factor in (.2,5.):
                w_new=w.copy()
                w_new[i]*=factor
                variants.append((w_new,o))
            if np.any(o):
                oo=o.copy()
                oo[i]*=-1
                variants.append((w,oo))
        for w,o in variants:
            key=tuple(w)+tuple(o)
            if key in trials:
                continue
            value=summarize(acq,a0+np.einsum("k,kij->ij",w,W)+np.einsum("k,kij->ij",o,L),w,o,modes)
            if value is not None:
                trials[key]=value
    choices=list(trials.values())
    plans=potential_plans(modes,choices)
    elapsed=perf_counter()-start
    verification=perf_counter()
    bounds=[]
    winner=None
    cached={}
    for numeric,variant,group in plans:
        lower=arb(0)
        accepted=[]
        for i,q in group:
            row=choices[i]
            if i not in cached:
                frame,edge=propose_frame(row["border"][:acq.n,:acq.n],acq.n)
                if frame is None:
                    continue
                H,rec=strict_schur(acq,ell,upp,row["weight"],row["optical"],frame)
                cached[i]=H,rec,frame
            H,rec,frame=cached[i]
            if H is None:
                continue
            typ=STAGE_MAP.get(variant,"I")
            value=project(H,q,acq.s,typ)
            lower+=max(arb(0),value.lower())
            accepted.append(dict(weights=row["weight"].tolist(),optical=row["optical"].tolist(),
                                 source_mode=q,basis=typ,
                                 frame=frames.add(frame),strict_mode_lower=low(value),verified_state=rec))
        excludes=bool(lower>arb(acq.radius)**2)
        bound=dict(variant=variant,strict_squared_lower=max(0.,low(lower)),excluded=excludes,
                   terms=accepted)
        if excludes and winner is None:
            winner=bound
        bounds.append(bound)
    return dict(excluded=winner is not None,initial_legal_PSD_seeds=initial_count,
                tested_legal_weights=len(trials),method_attempts=bounds,
                proposal_seconds=elapsed,arb_seconds=perf_counter()-verification)


def main():
    source,nodes=select_public()
    store=ArrayStore()
    outputs=[]
    for idx,(pid,target,slot,pair,row) in enumerate(nodes):
        t,b,w,_=row["box_positions"]
        ell,upp=source["boxes"][pair,b,w,slot]
        keep=upp>0.
        ar=Acquisition(source["G"][np.ix_(keep,keep)],source["S"][:,keep],source["B"][keep],
                       source["Y"][pair,t,b],float(source["radii"][pair,t,b]))
        outcome=one(ar,ell[keep],upp[keep],keep,store)
        variant_payload={}
        for j,attempt in enumerate(outcome["method_attempts"]):
            terms=[term | dict(S_space=term["basis"]) for term in attempt["terms"]]
            variant_payload[str(j)]=dict(excluded=attempt["excluded"],proof_terms=terms)
        independent_checks=independent(ar,ell[keep],upp[keep],variant_payload,store.arrays)
        outcome["independent_replay_seconds"]=sum(x["elapsed_seconds"] for x in independent_checks.values())
        outcome["independent_replay_count"]=len(independent_checks)
        physical=row["witness"]["audit"]["proved"]
        assert not (outcome["excluded"] and physical)
        status="EXCLUDED" if outcome["excluded"] else "WITNESS_FOUND" if physical else "UNRESOLVED"
        outputs.append(dict(pair_id=pid,target_id=target,condition_id=row["condition_id"],
                            status=status,proof=outcome))
    write_json(HERE/"seeded_same_vacuum.json",dict(
        scope="Original Phase-A seeded B proposals and legal local perturbations, with the same exact vacuum elimination as both follow-up methods",
        states=Counter(outcome["status"] for outcome in outputs),
        proposal_seconds=sum(r["proof"]["proposal_seconds"] for r in outputs),
        strict_seconds=sum(r["proof"]["arb_seconds"] for r in outputs),
        nodes=outputs))
    save(HERE/"seeded_same_vacuum_frames.npz",**store.arrays)


if __name__=="__main__":
    main()
