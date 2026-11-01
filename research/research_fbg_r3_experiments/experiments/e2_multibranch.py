#!/usr/bin/env python3
"""Frozen public E2 mechanism runs. This program never opens private labels.

Predictions are sealed in NPZ before any downstream scorer may read them.
Run from repository root using the project's existing NumPy/SciPy environment.
"""
from pathlib import Path
import argparse
from collections import Counter
import gzip
import json
import sys
import time
import numpy as np
from threadpoolctl import threadpool_limits

STUDY=Path(__file__).resolve().parents[1]
REPO=STUDY.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import load_model
from research.research_fbg_r3_experiments.src.algebraic_multiroot import run_multibranch


def clean(x):
    if isinstance(x,np.ndarray):
        return clean(x.tolist())
    if isinstance(x,(np.integer,np.floating)):
        return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):
        return None
    if isinstance(x,complex):
        return dict(real=x.real,imag=x.imag)
    if isinstance(x,dict):
        return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [clean(v) for v in x]
    return x


def plan(wall_seconds=1.):
    case_ids=json.loads((STUDY/'data/e1_case_ids.json').read_text())
    groups=sorted({(r['side'],r['layout']) for r in case_ids})
    chosen=[min((r for r in case_ids if (r['side'],r['layout'])==g),key=lambda r:r['case_id']) for g in groups]
    return dict(status='FROZEN_BEFORE_PREDICTIONS',label_access=False,
                case_selection='Lexicographic anonymous case ID per public side/layout; no morphology metadata',
                cases=chosen,starts=list(range(5)),sizes=[4,6],
                methods=['algebraic','block_trf','block_multistart','random','full_multistart','block_multistart_matched'],
                config=dict(max_supports=3,channels=1,max_candidates=3,corrector_nfev=40,
                            block_nfev=30,block_starts=6,proposal_evals=180,
                            seed=20261101,wall_seconds=wall_seconds,support_order='gradient'),
                support_selection='Per k: top2 paid all-data gradient blocks and one first unused cyclic block',
                common_prior='Frozen model smooth=1e-4 and TV=1e-4; common full-material TRF',
                stopping='One 1-second wall cap per run, checked between paid operations; max3 corrector seeds, 40 objective evaluations per seed',
                opportunity_scope='Equal upper seed quota and wall cap; actual valid/selected seed counts reported, no duplicate padding',
                outputs_closed_before_scoring=True)


def summarize(prefix):
    """Produce label-free summaries from closed artifacts, without new solves."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    closed=json.loads(Path(str(prefix)+'_results.json').read_text())
    rows=closed['rows']
    algebraic=[r for r in rows if r['method']=='algebraic']
    legal=[]
    gated_proposal_endpoints=0
    gated_fallback_runs=0
    with gzip.open(str(prefix)+'_ledger.jsonl.gz','rt') as f:
        for line in f:
            run=json.loads(line)
            if run.get('method')=='algebraic':
                legal.extend(r for r in run['roots'] if r.get('status')=='LEGAL_SUBSET_ROOT')
                gated_proposal_endpoints+=sum(c['relative_residual']<=.001 and c['provenance'].get('kind')!='NO_PROPOSAL_ANCHOR_FALLBACK' for c in run['corrected'])
                gated_fallback_runs+=int(run['relative_residual']<=.001 and any(c['provenance'].get('kind')=='NO_PROPOSAL_ANCHOR_FALLBACK' for c in run['corrected']))
    totals=[]
    for method in closed['plan']['methods']:
        selected=[r for r in rows if r['method']==method]
        totals.append(dict(method=method,runs=len(selected),data_gate_passes=sum(r['relative_residual']<=.001 for r in selected),
                           median_relative_residual=float(np.median([r['relative_residual'] for r in selected])),
                           total_wall_seconds=sum(r['wall_seconds'] for r in selected),
                           total_model_cost=dict(sum((Counter(r['model_cost']) for r in selected),Counter())),
                           total_kernel_algebraic_counts=dict(sum((Counter(r['kernel_algebraic_counts']) for r in selected),Counter())),
                           corrected_seed_distribution=dict(Counter(r['corrected_count'] for r in selected)),
                           recorded_phase_wall={p:sum(r.get('phase_wall',{}).get(p,0.) for r in selected)
                                                for p in ('reference','proposal','corrector')},
                           phase_wall_recorded_runs=sum('phase_wall' in r for r in selected)))
    matched=[]
    for a in algebraic:
        b=next(r for r in rows if (r['case_id'],r['k'],r['start'],r['method'])==
               (a['case_id'],a['k'],a['start'],'block_multistart_matched'))
        matched.append((a,b))
    pairs=[d for r in algebraic for d in r.get('same_chart_material_root_pairwise_rms',[])]
    setup={}
    for case in closed['plan']['cases']:
        _,public=load_model(STUDY/'data/e1_public'/(case['case_id']+'.npz'))
        setup[str((case['side'],case['layout']))]=float(public.get('green_setup_wall',np.nan))
    result=dict(status='CLOSED_PREDICTIONS_LABEL_FREE_SUMMARY',methods=totals,
                endpoint_filter_counts=dict(sum((Counter(r['filter_counts']) for r in algebraic),Counter())),
                charts=sum(len(r['charts']) for r in algebraic),
                chart_status_counts=dict(sum((Counter(c['status'] for c in r['charts']) for r in algebraic),Counter())),
                legal_roots=len(legal),legal_all_data_gate_passes=sum(r['all_data_relative_residual']<=.001 for r in legal),
                maximum_subset_relative_residual=max(r['subset_relative_residual'] for r in legal),
                legal_all_data_relative_residual_range=[min(r['all_data_relative_residual'] for r in legal),max(r['all_data_relative_residual'] for r in legal)],
                actual_material_root_pairwise_rms=dict(count=len(pairs),minimum=min(pairs),maximum=max(pairs),median=float(np.median(pairs))),
                algae_vs_matched_score=dict(better=sum(a['score']<b['score']-1e-8 for a,b in matched),
                                           tie=sum(abs(a['score']-b['score'])<=1e-8 for a,b in matched),
                                           worse=sum(a['score']>b['score']+1e-8 for a,b in matched)),
                raw_legal_root_reaching_gate_after_corrector=gated_proposal_endpoints,
                algebraic_gate_runs_with_empty_root_anchor_fallback=gated_fallback_runs,
                shared_preassembled_operator_green_setup_seconds=setup,
                construction_caveat='Per-run public NPZ loading and model constructor wall was not recorded; common cached operator import count=1/run. Shared original Green assembly times above are not per-run charges. All finite-kernel setup RHS, full model solves, small solves and algebraic work are recorded.',
                phase_caveat='Initial five execution-first runs reused without repetition; phase split missing for these five, total wall and all LU/RHS counters present.',
                independence_caveat='Four anonymous observations, not 40 independent materials; k4/k6 share observations and starts. No GT scores in this file.')
    Path(str(prefix)+'_summary.json').write_text(json.dumps(clean(result),indent=2)+'\n')
    fig,axes=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    axes[0,0].scatter([max(r['subset_relative_residual'],1e-16) for r in legal],
                      [r['all_data_relative_residual'] for r in legal],c=['#247ba0' if len(r['support'])==4 else '#e07830' for r in legal],s=30)
    axes[0,0].set(xscale='log',yscale='log',xlabel='Selected-channel relative residual',ylabel='All-data relative residual',title='33 legal subset roots; 0 all-data roots')
    axes[0,0].axhline(.001,color='crimson',linestyle='--',label='All-data gate 0.001')
    axes[0,0].legend(fontsize=8)
    labels=['Q6/Q20','block TRF','block multi','random','full multi','matched multi']
    colors=['#247ba0']+['#a9b8c7']*5
    axes[0,1].bar(labels,[r['data_gate_passes'] for r in totals],color=colors)
    axes[0,1].set(ylabel='Data-gate passes / 40 runs',ylim=(0,5),title='Q passes: both from empty-root anchor fallback')
    axes[0,1].tick_params(axis='x',labelrotation=25)
    axes[1,0].bar(labels,[r['total_model_cost'].get('LU',0) for r in totals],color=colors)
    axes[1,0].set(ylabel='Full physical LU count (all 40 runs)',title='Actual full solves, including correctors')
    axes[1,0].tick_params(axis='x',labelrotation=25)
    axes[1,1].bar(labels,[r['total_wall_seconds'] for r in totals],color=colors)
    axes[1,1].set(ylabel='Measured total wall / seconds',title='Equal wall ceilings were not binding')
    axes[1,1].tick_params(axis='x',labelrotation=25)
    fig.suptitle('FBG-R3-E / N2: fixed-outside branch proposals, public observations only',fontsize=12)
    (STUDY/'figures').mkdir(exist_ok=True)
    target=STUDY/'figures'/('e2_'+prefix.name.removeprefix('e2_')+'_summary.png')
    fig.savefig(target,dpi=160)
    plt.close(fig)
    print(json.dumps(dict(summary=str(prefix)+'_summary.json',figure=str(target))))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--wall-seconds',type=float,default=1.)
    ap.add_argument('--tag',default='mechanism')
    ap.add_argument('--max-runs',type=int)
    ap.add_argument('--resume-prefix',type=Path)
    ap.add_argument('--summarize-prefix',type=Path)
    args=ap.parse_args()
    if args.summarize_prefix:
        summarize(args.summarize_prefix)
        return
    frozen=plan(args.wall_seconds)
    outdir=STUDY/'results'
    outdir.mkdir(exist_ok=True)
    prefix=outdir/('e2_'+args.tag)
    Path(str(prefix)+'_plan.json').write_text(json.dumps(frozen,indent=2)+'\n')
    results,materials,predictions,archives=[],{},{},{}
    proposal_materials={}
    if args.resume_prefix:
        old=args.resume_prefix
        results=json.loads(Path(str(old)+'_results.json').read_text())['rows']
        for suffix,destination in (('_materials.npz',materials),('_predictions.npz',predictions),('_candidate_materials.npz',archives)):
            with np.load(str(old)+suffix) as f:
                destination.update({k:f[k] for k in f.files})
    completed={r['run_id']:r for r in results}
    started=time.perf_counter()
    with threadpool_limits(limits=1), gzip.open(str(prefix)+'_ledger.jsonl.gz','wt') as detail:
        if args.resume_prefix:
            with gzip.open(str(args.resume_prefix)+'_ledger.jsonl.gz','rt') as old:
                for line in old:
                    detail.write(line)
                    saved=json.loads(line)
                    for i,c in enumerate(saved.get('proposals',[])):
                        proposal_materials[saved['run_id']+f'_proposal{i}']=np.asarray(c['material'])
        for case in frozen['cases']:
            public=STUDY/'data/e1_public'
            starts=np.load(public/f"starts_{case['side']}.npz")['starts']
            for size in frozen['sizes']:
                for start in frozen['starts']:
                    for method in frozen['methods']:
                        if args.max_runs is not None and len(results)>=args.max_runs:
                            break
                        model,data=load_model(public/(case['case_id']+'.npz'))
                        cfg=dict(frozen['config'],method=method,sizes=(size,))
                        run_id=f"{case['case_id']}_k{size}_s{start}_{method}"
                        if run_id in completed:
                            continue
                        if method=='block_multistart_matched':
                            alg=completed[f"{case['case_id']}_k{size}_s{start}_algebraic"]
                            cfg.update(method='block_multistart',max_candidates=alg.get('corrected_count',0))
                        try:
                            out=run_multibranch(model,data['Y'],starts[start],cfg)
                            materials[run_id]=out['m']
                            predictions[run_id]=out['prediction']
                            for i,c in enumerate(out['corrected']):
                                archives[run_id+f'_candidate{i}']=c['m']
                            for i,c in enumerate(out['proposals']):
                                proposal_materials[run_id+f'_proposal{i}']=np.asarray(c['material'])
                            detail.write(json.dumps(clean(dict(run_id=run_id,**out)),ensure_ascii=False)+'\n')
                            roots=out['roots']
                            distances=[float(np.linalg.norm(np.asarray(r['material'])-starts[start])/np.sqrt(model.n))
                                       for r in roots if r.get('status')=='LEGAL_SUBSET_ROOT']
                            legal=[r for r in roots if r.get('status')=='LEGAL_SUBSET_ROOT']
                            pair_distances=[float(np.linalg.norm(np.asarray(a['material'])-np.asarray(b['material']))/np.sqrt(model.n))
                                            for i,a in enumerate(legal) for b in legal[:i]
                                            if a['support']==b['support'] and a['receiver']==b['receiver']]
                            row=dict(run_id=run_id,case_id=case['case_id'],side=case['side'],layout=case['layout'],
                                     k=size,start=start,method=method,status=out['status'],
                                     relative_residual=out['relative_residual'],score=out['score'],
                                     initial_relative_residual=out['initial_relative_residual'],
                                     wall_seconds=out['wall_seconds'],proposal_count=out['proposal_count'],
                                     corrected_count=out['corrected_count'],model_cost=out['model_cost'],
                                     phase_wall=out['phase_wall'],phase_model_cost=out['phase_model_cost'],
                                     kernel_algebraic_counts=out['ledger']['counts'],filter_counts=out['filter_counts'],
                                     charts=out['charts'],legal_root_anchor_material_rms=distances,
                                     same_chart_material_root_pairwise_rms=pair_distances,
                                     qualified_subset_roots=sum(r.get('status')=='LEGAL_SUBSET_ROOT' and r.get('subset_relative_residual',1.)<=1e-6 for r in roots),
                                     qualified_all_data_roots=sum(r.get('status')=='LEGAL_SUBSET_ROOT' and r.get('all_data_relative_residual',1.)<=.001 for r in roots),
                                     closed_candidate_material_keys=[run_id+f'_candidate{i}' for i in range(len(out['corrected']))])
                        except Exception as exc:
                            row=dict(run_id=run_id,case_id=case['case_id'],side=case['side'],layout=case['layout'],
                                     k=size,start=start,method=method,status='FAILED_EXECUTION',error=repr(exc),
                                     model_cost=dict(model.cost))
                            import traceback
                            detail.write(json.dumps(dict(run_id=run_id,error=traceback.format_exc()))+'\n')
                        results.append(row)
                        completed[run_id]=row
                        print(json.dumps(clean({k:row.get(k) for k in ('run_id','status','relative_residual','proposal_count','corrected_count','wall_seconds','error')})),flush=True)
                    if args.max_runs is not None and len(results)>=args.max_runs:
                        break
                if args.max_runs is not None and len(results)>=args.max_runs:
                    break
            if args.max_runs is not None and len(results)>=args.max_runs:
                break
    np.savez_compressed(str(prefix)+'_materials.npz',**materials)
    np.savez_compressed(str(prefix)+'_predictions.npz',**predictions)
    np.savez_compressed(str(prefix)+'_candidate_materials.npz',**archives)
    np.savez_compressed(str(prefix)+'_proposal_materials.npz',**proposal_materials)
    final=dict(status='PREDICTIONS_CLOSED_LABELS_NOT_READ',plan=frozen,
               elapsed_wall_seconds=time.perf_counter()-started,run_count=len(results),
               rows=results,prediction_files=[str(prefix)+'_materials.npz',str(prefix)+'_predictions.npz',
                                            str(prefix)+'_candidate_materials.npz'],
               scoring='NOT_RUN_BY_N2_WORKER; parent may invoke scorer only after this closure record')
    Path(str(prefix)+'_results.json').write_text(json.dumps(clean(final),indent=2)+'\n')
    print(json.dumps(dict(status=final['status'],runs=len(results),results=str(prefix)+'_results.json')),flush=True)
    summarize(prefix)


if __name__=='__main__':
    main()
