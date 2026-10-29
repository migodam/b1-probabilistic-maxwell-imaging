import json,time
from collections import Counter,defaultdict
from accounting import P,Book,write

def main():
    book=Book('FINAL_ACCOUNTING_AND_PROVENANCE_ADDENDUM')
    rows=[json.loads(s) for s in (P/'COST_LEDGER.jsonl').read_text().splitlines()]
    # Preserve six original combined-snapshot IDs; add referential links, never rewrite old rows.
    src=[r for r in rows if r['run']=='INDEPENDENT_70DPS_FULL_MAXWELL'];paid=[];links=[]
    for row in src:
        if row['event']=='BEGIN' and row['kind']=='HIGH_PRECISION_FULL_SOURCE_SOLVE':paid.append(row['action'])
        if row['event']=='SAVE':
            links.append(dict(original_snapshot=row['action'],path=row['path'],paid_source_actions=paid[-2:]))
    book.append(dict(event='COMBINED_HIGH_PRECISION_SNAPSHOT_LINK_ADDENDUM',links=links,new_complex_RHS=0,note='Original snapshot IDs were preceding END sequence numbers; original data/ledger unchanged'))
    write(P/'audit/HIGH_PRECISION_SNAPSHOT_LINK_ADDENDUM.json',dict(status='PROVENANCE_LINKS_EXPLICIT',links=links,new_RHS=0))
    book.end()
    rows=[json.loads(s) for s in (P/'COST_LEDGER.jsonl').read_text().splitlines()]
    begins=[r for r in rows if r['event']=='BEGIN'];ends=[r for r in rows if r['event']=='RUN_END']
    total=sum(r.get('actual_complex_RHS',0) for r in begins);formula=sum(r.get('closed_form_source_RHS_equivalent',0) for r in begins)
    jobs=[]
    for e in ends:
        br=[r for r in begins if r['run']==e['run']]
        jobs.append(dict(run=e['run'],status=e['status'],actual_RHS=sum(r.get('actual_complex_RHS',0) for r in br),formula_equivalent=sum(r.get('closed_form_source_RHS_equivalent',0) for r in br),wall=e['controlled_wall_seconds'],CPU=e['controlled_CPU_seconds'],peak_RSS=e['peak_RSS_bytes']))
    active={r['run'] for r in rows if r['event']=='RUN_START'}-{r['run'] for r in ends}
    assert not active and total+formula<=30000 and all(j['actual_RHS']+j['formula_equivalent']<=10000 for j in jobs)
    LU=sum(r.get('LU',0) for r in begins)
    out=dict(status='CLOSED',actual_complex_RHS=total,closed_form_source_RHS_equivalent=formula,conservative_combined=total+formula,
        hard_cap_actual=30000,stronger_combined_cap=30000,unused_actual=30000-total,unused_combined=30000-total-formula,
        LU_factorizations=LU,closed_form_material_states=sum(r.get('states',0) for r in begins),
        forward_new_states_explicit=sum(r.get('new_states',0) for r in begins),legacy_forward_states_extra=8,
        physical_and_algebraic_job_wall_seconds=sum(j['wall'] for j in jobs),job_CPU_seconds=sum(j['CPU'] for j in jobs),peak_RSS=max(j['peak_RSS'] for j in jobs),
        neural_training=0,GPU_hours=0,physical_workers=1,BLAS_threads=1,live_physical_LU_max=1,jobs=jobs,
        accounting_qualification='Full controlled job clocks and actions are measured; drafting/literature/plotting/delivery have separate receipts, not fabricated into the physics clock',
        historical_cost='No R6/R7 cost transfer or double count',raw_dtype='float64/complex128, np.savez_compressed lossless',no_active_jobs=True)
    write(P/'RESOURCE_SUMMARY.json',out)
    write(P/'CLOSURE.json',dict(status='CLOSED_COMPLETE_E0_E1_E2_E3_E4',actual_RHS=total,formula_equivalent=formula,combined=total+formula,cap=30000,unused_combined=30000-total-formula,
        reason='Same physical coverage; optimized equal-trace total saving below20%; no isolated second branch discovered; local probability numerical pass with global mass unqualified',no_active_jobs=True,utc=time.time()))
    print(json.dumps({k:v for k,v in out.items() if k!='jobs'},indent=2))
if __name__=='__main__':main()
