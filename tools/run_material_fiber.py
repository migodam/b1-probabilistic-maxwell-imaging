from pathlib import Path
import sys,json,argparse,traceback,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,Book,Target,BudgetStop,BlockStop,write,save
from material_fiber.optimize import reconstruct

def stage_A(t):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());records=[];points={};methods=['fixed_LM','adaptive_TR','material_TR','coarse_TR']
    # Post-closure hardening: generation must finish in a different OS process.
    # The executed A version e5bd9d3 used same-process function isolation; disclosed.
    if not (OUT/'A/OBSERVATIONS_GENERATED.json').exists():
        raise RuntimeError('Run the separately budgeted observation process first; no online truth reads')
    dev=[r['parent'] for r in parents if r['split']=='train'][:4]+['P1_obs_B']
    for parent in dev:
        for method in methods:
            init=t.load('OLD_B') if parent=='P1_obs_B' else None
            p,r=reconstruct(t,parent,method,8,init,run=f'{parent}_{method}_checkpoint8');points[parent,method]=p;records.append(r)
    ranking=[]
    for method in methods:
        rr=[r for r in records if r['method']==method]
        ratios=[]
        for r in rr:
            hist=json.loads((OUT/r['history']).read_text());ratios.append(r['objective']/max(hist[0]['objective'],1e-30))
        ranking.append(dict(method=method,qualified=sum(r['data_feasible'] for r in rr),median_objective_ratio=float(np.median(ratios)),RHS=sum(r['incremental_RHS'] for r in rr)))
    ranking.sort(key=lambda r:(-r['qualified'],r['median_objective_ratio'],r['RHS'],r['method']))
    top=[r['method'] for r in ranking[:2]];write(OUT/'A/OPTIMIZER_SELECTION_8.json',dict(ranking=ranking,continued=top,GT_used=False))
    finalists=[]
    for parent in dev:
        for method in top:
            p,r=reconstruct(t,parent,method,8,points[parent,method],run=f'{parent}_{method}_checkpoint16',iteration_offset=8)
            points[parent,method]=p;finalists.append(r)
    scores=[]
    for method in top:
        rr=[r for r in finalists if r['method']==method];rat=[]
        for r in rr:
            first=json.loads((OUT/f'A/histories/{r["parent"]}_{method}_checkpoint8.json').read_text())[0]['objective'];rat.append(r['objective']/first)
        scores.append(dict(method=method,qualified=sum(r['data_feasible'] for r in rr),median_objective_ratio=float(np.median(rat))))
    scores.sort(key=lambda r:(-r['qualified'],r['median_objective_ratio'],r['method']));winner=scores[0]['method']
    write(OUT/'A/OPTIMIZER_FROZEN.json',dict(winner=winner,scores=scores,selection='train/exposed B only; no GT',heldout_GT_opened=False))
    anchors={}
    for parent in dev[:-1]:
        pp=[points[parent,m] for m in top];y=t.observation(parent);p=min(pp,key=lambda p:(not t.qualified(p,y),t.objective(p,y)));anchors[parent]=p.key
    for r in parents:
        parent=r['parent']
        if parent in anchors:continue
        p,rec=reconstruct(t,parent,winner,16,run=f'{parent}_{winner}_frozen');anchors[parent]=p.key;finalists.append(rec)
    report=[]
    for r in parents:
        p=t.load(anchors[r['parent']]);report.append(t.record(p,r['parent'],winner,phase='frozen_anchor',split=r['split']))
    write(OUT/'A/ANCHORS.json',anchors);write(OUT/'A/SUMMARY.json',dict(winner=winner,parents=report,qualified=sum(r['data_feasible'] for r in report),count=len(report),GT_online=False,old_P1_ledger_modified=False))
    print(json.dumps(dict(stage='A',qualified=sum(r['data_feasible'] for r in report),parents=len(report),winner=winner,consumed=dict(t.book.used['A']))),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=list('ABCDE'));ap.add_argument('--recover');args=ap.parse_args()
    b=Book(args.stage,args.recover);t=Target(b);status='COMPLETE';error=None
    try:
        if args.stage=='A':stage_A(t)
        elif args.stage=='B':
            if args.recover:
                from material_fiber.baseline_audit import run
                run(t)
            else:
                from material_fiber.explore import stage_B
                stage_B(t)
        elif args.stage=='C':
            from material_fiber.global_search import stage_C
            stage_C(t)
        elif args.stage=='D':
            from material_fiber.learn import stage_D
            stage_D(t)
        else:
            from material_fiber.acceleration import stage_E
            stage_E(t)
    except BudgetStop as e:status='BUDGET_LIMIT';error=repr(e);print(error,flush=True)
    except BaseException as e:
        status='ERROR';error=repr(e);write(OUT/f'{args.stage}/FAILURE.json',dict(error=error,traceback=traceback.format_exc()));raise
    finally:b.close(status,error)

if __name__=='__main__':main()
