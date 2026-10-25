"""Saved evidence/code/link consistency only; no PDE or new scientific gate."""
from pathlib import Path
import json,re
import numpy as np
from wslam.accounting import BudgetBook,json_default

ROOT=Path(__file__).resolve().parents[1]

def main():
    required=['START_HERE.md','EXPERIMENT_CONTRACT.json','PHYSICS_PROTOCOL.md',
        'PRIOR_AND_LIKELIHOOD.md','DYNAMIC_COMPENSATION_THEORY.md',
        'E1_MAXWELL_OBSERVABILITY.md','A0_ONLINE_RECONSTRUCTION.md','E2_BELIEF_MAP.md',
        'E3_REVISITATION_ABLATION.md','E4_POSTERIOR_VALIDATION.md','E5_ACTIVE_SENSING.md',
        'SCIENTIFIC_DECISION.md','GATE_SUMMARY.json','RESEARCH_W_SLAM_REPORT_ZH.md',
        'PHYSICS_LEDGER.jsonl','FAILURE_LEDGER.jsonl','REPRODUCE.md']
    for name in required:
        assert (ROOT/name).is_file(),name
    compiled=[]
    for folder in ['src','tools','tests']:
        for p in sorted((ROOT/folder).rglob('*.py')):
            compile(p.read_text(),str(p),'exec');compiled.append(str(p.relative_to(ROOT)))
    links=[]
    for name in required:
        if not name.endswith('.md'):continue
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',(ROOT/name).read_text()):
            if target.startswith(('http:','https:','#')):continue
            path=target.split('#')[0]
            assert (ROOT/path).exists(),(name,target)
            links.append(dict(source=name,target=target))
    with np.load(ROOT/'results/e2e3/inference/POSTERIORS.npz',allow_pickle=False) as z:
        assert bool(z['complete']) and z['completed_mask'].all()
        ref=z['weights'];mw=z['method_weights'];methods=z['methods'].astype(str).tolist()
        assert ref.shape==(12,2,3,2950) and mw.shape==(12,2,3,11,2950)
        assert np.isfinite(ref).all() and np.isfinite(mw).all()
        assert np.allclose(ref.sum(-1),1,rtol=0,atol=1e-10)
        assert np.allclose(mw.sum(-1),1,rtol=0,atol=1e-10)
        difference=float(np.max(np.abs(mw[:,:,:,methods.index('INCREMENTAL_JOINT_SMOOTHER')]-ref)))
        assert difference<=1e-10
    logs=json.loads((ROOT/'results/e2e3/inference/METHOD_RESULTS.json').read_text())
    assert len(logs)==792
    guard_rows=[r for r in logs if r.get('criterion') and not r['criterion'].get('defined',False)]
    assert all(r['action']=='REOPEN' for r in guard_rows)
    proxy=[r for r in logs if r['method']=='ONLINE_TASK_WEIGHTED_PROXY' and r['time']>1]
    assert all(r['criterion']['current_reference_access_for_decision'] is False for r in proxy)
    for rel,count in [('results/e4/comparison/RESULTS.json',384)]:
        report=json.loads((ROOT/rel).read_text());assert len(report['rows'])==count
    e5=json.loads((ROOT/'results/e5/finite_sensing/MANIFEST.json').read_text())
    assert e5['status']=='COMPLETE' and e5['completed_cases']==24
    assert len(json.loads((ROOT/'results/e5/finite_sensing/ONLINE_RESULTS.json').read_text()))==168
    budgets=BudgetBook().snapshot()
    for stage,values in budgets.items():
        assert all(v>=0 for v in values['remaining'].values()),stage
        assert all(values['used'][k]<=v for k,v in values['cap'].items()),stage
    out=dict(status='PASS',scope='DELIVERY_CONSISTENCY_ONLY_NOT_SCIENTIFIC_APPROVAL',
        compiled_sources=compiled,local_report_links=links,reference_weight_error=difference,
        conservative_gaussian_reopen_records=len(guard_rows),complete_E4_rows=384,complete_E5_rows=168,
        budget_snapshot=budgets,new_PDE=0,new_gates=False)
    (ROOT/'FINAL_DELIVERY_CHECK.json').write_text(json.dumps(out,indent=2,default=json_default)+'\n')
    print(json.dumps(dict(status=out['status'],compiled=len(compiled),local_links=len(links),new_PDE=0)))

if __name__=='__main__':main()
