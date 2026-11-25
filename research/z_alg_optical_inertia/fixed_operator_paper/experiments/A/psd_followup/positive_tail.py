"""All-box PSD strengthening using the common bank's full positive Schur trace."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from common import PUBLIC, RESULTS, digest, herm, load, save, star, write_json
from bounds import (Acquisition, arb, assemble, ball, bt, ctx, low, propose_frame,
                    select, strict_pd, validate_frame, weight_vector)
from replay import original
from run import ArrayStore, PARAMS
import numpy as np
from scipy.linalg import eigh


ctx.prec=192
ctx.threads=1
OUT=HERE/"positive_tail"
DECISIONS=RESULTS/"DECISION_SEAL.json"


def numeric_bank(acq,ell,upp):
    eligible=[]
    near_boundary=0
    for key in PARAMS:
        A=assemble(acq,ell,upp,*key)
        n=acq.n
        K=A[:n,:n]
        values,basis=eigh(K,check_finite=False)
        if values[0]<=1.e-12*max(1.,values[-1]):
            near_boundary+=bool(values[0]>0.)
            continue
        coupling=star(basis)@A[:n,n:]
        H=herm(A[n:,n:]-star(coupling)@(coupling/values[:,None]))
        h,U=eigh(H,check_finite=False)
        pos=h>0.
        score=float(np.sum(h[pos]))
        eligible.append((score,key,K,U[:,pos],float(values[0])))
    eligible.sort(key=lambda a:a[0],reverse=True)
    return eligible,near_boundary


def positive_projector_trace(acq,A,Kframe,F):
    """Return tr(P_F H), with P_F the exact projector onto range(F)."""
    n,s=acq.n,acq.s
    if Kframe.shape!=(n,n) or F.shape[0]!=s or F.shape[1]==0:
        return None
    state,src=range(n),range(n,n+s)
    K=select(A,state,state)
    pd,krec=validate_frame(K,Kframe)
    if not pd>0:
        return None
    L=select(A,state,src)
    H=select(A,src,src)-bt(L)*K.solve(L)
    Q=ball(F)
    gram=bt(Q)*Q
    is_pd,grec=strict_pd(gram)
    if not is_pd:
        return None
    chosen=bt(Q)*H*Q
    value=gram.solve(chosen)
    trace=sum((value[j,j].real for j in range(value.nrows())),arb(0))
    return trace,krec,grec


def compare_one(acq,ell,upp,row,frames):
    start=perf_counter()
    bank,near = numeric_bank(acq,ell,upp)
    after=perf_counter()
    failures=Counter()
    cert=None
    replay_time=0.
    for score,key,K,F,smallest in bank:
        if score<=acq.radius**2:
            break
        try:
            Kframe,edge=propose_frame(K,acq.n)
            if Kframe is None:
                failures["K_positive_frame_not_proposed"]+=1
                continue
            proof=positive_projector_trace(acq,assemble(acq,ell,upp,*key,strict=True),Kframe,F)
            if proof is None:
                failures["positive_state_or_projector_Gram_not_proved"]+=1
                continue
            lower,kpd,fpd=proof
            if not lower>arb(acq.radius)**2:
                failures["strict_projection_below_noise"]+=1
                continue
            start_replay=perf_counter()
            w,o=weight_vector(acq.n,*key)
            replay=positive_projector_trace(acq,original(acq,ell,upp,w,o),Kframe,F)
            replay_time+=perf_counter()-start_replay
            if replay is None or not replay[0]>arb(acq.radius)**2:
                raise ArithmeticError("independent projected PSD certificate replay failed")
            cert=dict(weights=list(key),numeric_smallest_K=smallest,strict_K_inertia=(acq.n,0,0),
                      positive_K=kpd,projector_rank=F.shape[1],positive_projector_Gram=fpd,
                      state_frame=frames.add(Kframe),source_frame=frames.add(F),
                      numeric_trace_positive_part=score,
                      strict_squared_residual_lower=low(lower),
                      strict_residual_over_radius_lower=low(lower.sqrt()/arb(acq.radius)),
                      replay_strict_squared_residual_lower=low(replay[0]))
            break
        except (ValueError,ZeroDivisionError) as exc:
            failures[type(exc).__name__]+=1
    if cert and row["status"]["B"]=="WITNESS_FOUND":
        raise ArithmeticError("PSD positive spectral part contradicts strict original witness")
    return dict(common_bank_trials=len(PARAMS),numerical_positive_K_trials=len(bank),
                numerical_near_rank_positive_K_trials=int(near),
                numerical_max_positive_projector_lower_squared=bank[0][0] if bank else 0.,
                numerically_above_noise_proposals=sum(p[0]>acq.radius**2 for p in bank),
                certificate=cert,
                proposal_seconds=after-start,
                online_proof_seconds=perf_counter()-after-replay_time,
                independent_proof_seconds=replay_time,failures=failures,
                status=("EXCLUDED" if cert else row["status"]["B"]))


def shard(begin,end):
    a=load(PUBLIC)
    primary_seal=json.loads(DECISIONS.read_text())
    OUT.mkdir(exist_ok=True)
    for pair in range(begin-1,end):
        p=OUT/f"pair_{pair+1:02d}.json"
        if p.exists():
            continue
        raw=RESULTS/"raw"/p.name
        rep=RESULTS/"replay_inertia"/p.name
        assert primary_seal["files"][f"raw/{p.name}"]==digest(raw)
        assert json.loads(rep.read_text())["pass_all"]
        rows=json.loads(raw.read_text())["rows"]
        t0=perf_counter()
        records=[]
        store=ArrayStore()
        current=None
        for row in rows:
            t,b,w,slot=row["box_positions"]
            if (t,b)!=current:
                current=t,b
                acq=Acquisition(a["G"],a["S"],a["B"],a["Y"][pair,t,b],float(a["radii"][pair,t,b]))
            ell,upp=a["boxes"][pair,b,w,slot]
            decision=compare_one(acq,ell,upp,row,store)
            pC=row["methods"]["C"]["proof"]
            records.append(dict(node=row["box_positions"],pair_id=row["pair_id"],target_id=row["target_id"],
                                condition_id=row["condition_id"],old_B=row["status"]["B"],
                                C=row["status"]["C"],C_proof_negative_K=None if pC is None else
                                pC["numeric_K_negative_count"],**decision))
        save(OUT/f"pair_{pair+1:02d}_frames.npz",**store.arrays)
        write_json(p,dict(pair=pair+1,rows=records,wall_seconds=perf_counter()-t0,
                          primary_B_replay_path=str(rep.relative_to(RESULTS))))
        print(pair+1,Counter(r["status"] for r in records),flush=True)


def summary():
    paths=sorted(OUT.glob("pair_*.json"))
    assert len(paths)==24
    rows=[r for p in paths for r in json.loads(p.read_text())["rows"]]
    assert len(rows)==2688
    cost=Counter()
    old,new,c=Counter(),Counter(),Counter()
    comparison=Counter()
    failures=Counter()
    for r in rows:
        old[r["old_B"]]+=1
        new[r["status"]]+=1
        c[r["C"]]+=1
        for k in ("proposal_seconds","online_proof_seconds","independent_proof_seconds"):
            cost[k]+=r[k]
        failures.update(r["failures"])
        comparison[("primary_C","old_B",r["C"],r["old_B"])]+=1
        comparison[("primary_C","strong_B",r["C"],r["status"])]+=1
    primary_nu0_gap=[r for r in rows if r["C"]=="EXCLUDED" and r["old_B"]!="EXCLUDED"
                     and r["C_proof_negative_K"]==0]
    no_nu0_gap=[r for r in primary_nu0_gap if r["status"]!="EXCLUDED"]
    added=[r for r in rows if r["status"]=="EXCLUDED" and r["old_B"]!="EXCLUDED"]
    receipt=dict(scope="Post-hoc strengthening of the primary B finite-bank diagnostic. Uses no truth or outside data.",
                 node_count=len(rows),common_weights_per_node=72,
                 all_numerical_positive_K_bank_proposals=sum(r["numerical_positive_K_trials"] for r in rows),
                 near_rank_numerical_positive_K_proposals=sum(r["numerical_near_rank_positive_K_trials"] for r in rows),
                 strict_projected_positive_part_proofs=sum(r["certificate"] is not None for r in rows),
                 primary_B=old,strengthened_B=new,primary_C=c,
                 new_primary_B_exclusions=len(added),new_primary_B_only_exclusions=sum(r["C"]!="EXCLUDED" for r in added),
                 primary_C_nu0_gaps=len(primary_nu0_gap),strictly_recovered_primary_C_nu0_gaps=len(primary_nu0_gap)-len(no_nu0_gap),
                 unresolved_primary_C_nu0_gaps=[(r["pair_id"],r["node"]) for r in no_nu0_gap],
                 combined_table={"|".join(k):v for k,v in comparison.items()},
                 timings_seconds=cost,certificate_failures=failures,
                 numeric_search_claim="All 72 proposed weights scanned. Positive K and projector rigor checked for selected qualifying certificates, not all weights globally optimized.",
                 B_original_independent_replay="Reuses frozen Phase A replay receipts, each verified per pair.",
                 replay="Arb 192 positive K/Gram, stationary Schur solve, exact source projector trace. Reassembled original row block form.",
                 accepted_frames=sorted(str(p.relative_to(HERE)) for p in OUT.glob("*.npz")))
    write_json(HERE/"positive_tail_summary.json",receipt)
    with (HERE/"positive_tail_nodes.csv").open("w",newline="") as handle:
        cols=["pair_id","target_id","condition_id","old_B","C","C_proof_negative_K","status",
              "numerical_positive_K_trials","numerically_above_noise_proposals","strict_ratio",
              "proposal_seconds","online_proof_seconds","independent_proof_seconds"]
        writer=csv.DictWriter(handle,fieldnames=cols,extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            proof=row["certificate"]
            writer.writerow(row|dict(strict_ratio=proof["strict_residual_over_radius_lower"] if proof else None))
    print(json.dumps(receipt,indent=2))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--begin",type=int,default=1)
    parser.add_argument("--end",type=int,default=24)
    parser.add_argument("--summary",action="store_true")
    args=parser.parse_args()
    if args.summary:
        summary()
    else:
        shard(args.begin,args.end)


if __name__=="__main__":
    main()
