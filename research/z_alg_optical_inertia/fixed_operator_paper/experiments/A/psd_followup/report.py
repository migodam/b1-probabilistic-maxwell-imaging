"""Seal the 32 public decisions, score frozen target flags and publish complete costs."""
from collections import Counter,defaultdict
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from common import RESULTS, digest, write_json

VARIANTS=("spatial_psd","signed_local_qcqp","source_diagonal_qcqp","source_hadamard_matrix_qcqp")


def main():
    out=HERE/"nodes"
    nodes=sorted(out.glob("*.json"))
    assert len(nodes)==32,len(nodes)
    data=[json.loads(p.read_text()) for p in nodes]
    assert all(len(d["B_variants"])==4 and len(d["B_independent_replay"])==4 for d in data)
    initial=[json.loads(p.read_text()) for p in sorted((HERE/"diagnostic_initial_clarabel"/"nodes").glob("*.json"))]
    seeded=json.loads((HERE/"seeded_same_vacuum.json").read_text())["nodes"]
    assert len(seeded)==32
    paths=nodes+sorted(out.glob("*.npz"))+[HERE/"seeded_same_vacuum.json",HERE/"seeded_same_vacuum_frames.npz"]
    receipt=HERE/"DECISION_SEAL.json"
    frozen={str(p.relative_to(HERE)):digest(p) for p in paths}
    if receipt.exists():
        assert json.loads(receipt.read_text())["files"]==frozen
    else:
        write_json(receipt,dict(utc=datetime.now(timezone.utc).isoformat(),files=frozen,
                                selected_before_score=True))
    with (RESULTS/"all_nodes.csv").open() as h:
        labels={(r["pair_id"],r["target_id"],r["condition_id"]):(
            r["generating"]=="True",r["contains_generating_material"]=="True")
                for r in csv.DictReader(h)}
    statuses=defaultdict(Counter)
    opposite=defaultdict(Counter)
    disagreement=[]
    solvers=Counter()
    variant_totals=defaultdict(Counter)
    time=Counter()
    worst_B=0.
    worst_C=0.
    worst_B_union=0.
    union_diff=[]
    for idx,row in enumerate(data):
        key=tuple(row[name] for name in ("pair_id","target_id","condition_id"))
        generating,contained=labels[key]
        b=row["B_independent_state"]
        c=row["C_known_zero"]["C_known_zero"]
        truth_condition="truth" if contained else "opposite_noncover"
        seed=seeded[idx]
        assert tuple(seed[name] for name in ("pair_id","target_id","condition_id"))==key
        bs=seed["status"]
        combined="EXCLUDED" if "EXCLUDED" in (bs,b) else b
        if combined!=c:
            union_diff.append(dict(pair_id=key[0],target_id=key[1],condition_id=key[2],combined=combined,C=c))
        for m,s in (("new_B",b),("C_known_zero",c),("seeded_B_same_vacuum",bs),
                    ("B_union_same_vacuum",combined),
                    ("original_B_no_exact_zero_elimination",row["old_B_status"])):
            statuses[m][s]+=1
            if not contained:
                opposite[m][s]+=1
            if s=="EXCLUDED" and contained:
                raise ArithmeticError("a truth-containing node was removed")
        if c!=b:
            disagreement.append(dict(pair_id=key[0],target_id=key[1],condition_id=key[2],
                                     C=c,B=b,contained=contained))
        for variant in VARIANTS:
            v=row["B_variants"][variant]
            proposal=v.get("proposal",{})
            for solver in proposal.get("solver_receipts",[]):
                solvers[solver["solver"],solver["status"]]+=1
            variant_totals[variant]["excluded"]+=v["excluded"]
            variant_totals[variant]["proof_terms"]+=len(v.get("proof_terms",[]))
            variant_totals[variant]["solver_wall_seconds"]+=proposal.get("solver_wall_seconds",0.)
            variant_totals[variant]["proposal_build_seconds"]+=proposal.get("build_seconds",0.)
            variant_totals[variant]["strict_seconds"]+=v.get("certificate_seconds",0.)
            variant_totals[variant]["independent_replay_seconds"]+=row["B_independent_replay"][variant]["elapsed_seconds"]
            variant_totals[variant]["budget_failures"]+=v.get("status") in ("wall_cap","unattempted_wall_cap")
            time["B_independent_replay"]+=row["B_independent_replay"][variant]["elapsed_seconds"]
        time["B_online"]+=row["B_end_to_end_online_seconds"]
        time["shared_setup"]+=row["shared_setup_seconds"]
        time["C_proposal"]+=row["C_known_zero"]["proposal_seconds"]
        time["C_strict_including_independent_replay"]+=row["C_known_zero"]["strict_including_replay_seconds"]
        time["B_seeded_proposal"]+=seed["proof"]["proposal_seconds"]
        time["B_seeded_strict"]+=seed["proof"]["arb_seconds"]
        time["B_seeded_replay"]+=seed["proof"]["independent_replay_seconds"]
        b_node=row["B_end_to_end_online_seconds"]+row["shared_setup_seconds"]+sum(
            result["elapsed_seconds"] for result in row["B_independent_replay"].values())
        b_union_node=b_node+sum(seed["proof"][k] for k in
                                ("proposal_seconds","arb_seconds","independent_replay_seconds"))
        worst_B=max(worst_B,b_node)
        worst_B_union=max(worst_B_union,b_union_node)
        worst_C=max(worst_C,row["C_total_proposal_and_replay_seconds"]+row["shared_setup_seconds"])
    write_json(HERE/"summary.json",dict(
        design="Exact vacuum b=0, width ratio 3, all 8 independent structure pairs, both measurements and both public boxes",
        independent_pairs=8,targets=16,nodes=32,opposite_noncover=sum(opposite["new_B"].values()),
        methods=statuses,opposite_noncover_outcomes=opposite,variant_counts_and_cost=variant_totals,
        solver_outcomes={"|".join(k):n for k,n in solvers.items()},
        initial_clarabel_only_counts=Counter(x["B_independent_state"] for x in initial),
        C_B_disagreement=disagreement,B_union_C_disagreement=union_diff,zero_truth_exclusion=True,
        aggregate_wall_seconds=time,maximum_node_seconds=dict(B=worst_B,B_union_sum=worst_B_union,C=worst_C),
        permitted_B_C_per_node_seconds=10.,
        strict_replay="Arb 192-bit independent original-coordinate assembly verified every saved bound and all decisions",
        implication="Claims compare the finite tested cone and inherited inertia bank on this exposed structure diagnostic. The result does not bound the best possible matrix-valued multiplier."))


if __name__=="__main__":
    main()
