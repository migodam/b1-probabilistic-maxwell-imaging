"""Describe frozen scored endpoints and paid CPU slices; no new solves or GT reads."""
import csv,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def main():
    s=ROOT/'RESULTS/confirmation_scoring';rows=json.loads((s/'SCORES.json').read_text())['rows']
    budgets=json.loads((s/'BUDGET_SCORES.json').read_text())['rows']
    out=ROOT/'RESULTS/summary';out.mkdir(exist_ok=True)
    groups=list(dict.fromkeys((r['method'],r.get('dpcs_rule')) for r in rows));summary=[];slices=[]
    for method,rule in groups:
        rr=[r for r in rows if r['method']==method and r.get('dpcs_rule')==rule]
        summary.append(dict(method=method,rule=rule,objects=len(rr),successes=sum(r['material_image_success'] for r in rr),
            mean_NRMSE=float(np.mean([r['material_nrmse'] for r in rr])),median_NRMSE=float(np.median([r['material_nrmse'] for r in rr])),
            mean_CPU_seconds=float(np.mean([r.get('standalone_CPU_seconds',r['cpu_seconds']) for r in rr])),
            scope='six independent finite-model objects; descriptive not population significance'))
        if rule=='two_frozen_candidates':continue
        for bound in (30.,60.):
            selected=[]
            for r in rr:
                b=[x for x in budgets if x['task_id']==r['task_id'] and x.get('standalone_CPU_seconds',x['CPU_seconds'])<=bound]
                chosen=max(b,key=lambda x:x.get('standalone_CPU_seconds',x['CPU_seconds'])) if b else None
                if r.get('standalone_CPU_seconds',r['cpu_seconds'])<=bound:chosen=r
                if chosen is not None:selected.append(chosen)
            slices.append(dict(method=method,rule=rule,CPU_budget=bound,covered_objects=len(selected),
                successes=sum(r['material_image_success'] for r in selected),
                mean_NRMSE=float(np.mean([r['material_nrmse'] for r in selected])) if selected else None,
                selection='latest actually paid original-Phi incumbent by CPU; final exact endpoint if task already completed; no interpolation or GT selection',
                status='POSTHOC_BUDGET_DESCRIPTION; not algorithm calibration',
                expensive_profile_cost='all current Gram/Cholesky included in CPU; not free physics'))
    for name,data in [('CONFIRMATION_SUMMARY',summary),('PAID_CPU_SLICES',slices)]:
        (out/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
        with (out/(name+'.csv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    print(json.dumps(dict(summary=summary,CPU_slices=slices),indent=2))

if __name__=='__main__':main()
