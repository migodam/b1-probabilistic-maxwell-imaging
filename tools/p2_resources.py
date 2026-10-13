"""Append-only attribution supplements, never rewrite old ledger receipts."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.accounting import Ledger,dump
def main():
    out=ROOT/'research/joint_p2/results';ledger=Ledger(out);rows=[json.loads(a) for a in (out/'JOB_LEDGER.jsonl').read_text().splitlines()];overlap={'p2-fresh-continuous-reference':0,'p2-fresh-continuous-reference-batched':0,'p2-context-partition-refinement-resume':0}
    dump(out/'JOB_ATTRIBUTION_CORRECTION.json',dict(reason='historical job global start/end count deltas included unrelated concurrent reservations; unique BEGIN rows have always enforced the sole4000 budget',immutable_original_rows=True,overrides=[dict(job_id=j,raw_global_delta=next(a['new_physics'] for a in rows if a['id']==j),actual_own_new_physics=v) for j,v in overlap.items()],owner_tags_added_to_future_BEGIN=True,sole_total_new_physics=ledger.count(),not_a_budget_reset=True))
    summary=ledger.summary();receipts=[]
    for p in (out/'jobs').glob('*/CHILD_CPU.json'):
        receipts.append(dict(path=str(p.relative_to(out)),**json.loads(p.read_text())))
    summary.update(controlled_wall_measure='sum of metered numerical job wall including failures/repeats; overlapping CPU jobs charged conservatively',job_physics_attribution='JOB_ATTRIBUTION_CORRECTION.json; never sum historical global deltas',child_cpu_receipts=receipts,remaining_stage_cap={s:cap-summary['by_stage'][s] for s,cap in {'mechanism':400,'reference':2000,'fresh':800,'reserve':800}.items()})
    dump(out/'RESOURCE_SUMMARY.json',summary)
if __name__=='__main__':main()
