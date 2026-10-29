"""Publication figures after frozen scoring. GT is used only for labeled plots."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parents[1]


def step_value(xs,ys,query):
    xs=np.asarray(xs);ys=np.asarray(ys)
    indices=np.searchsorted(xs,query,side='right')-1
    return np.where(indices>=0,ys[np.maximum(indices,0)],np.nan)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--split',choices=['dev','holdout','confirm'],default='holdout')
    p.add_argument('--output-prefix',default='holdout')
    p.add_argument('--presentation-note',default='V6.2 fixed reciprocal FullJ control for all cases; original native control panels archived; no new physics or candidate changes')
    p.add_argument('--baseline-token',default='FullJReciprocal1024')
    p.add_argument('--baseline-label',default='Reciprocal FullJ')
    p.add_argument('--r3-token',required=True)
    p.add_argument('--r3-label',required=True)
    a=p.parse_args()
    payload=json.loads(a.scores.read_text())
    if not payload.get('truth_read_after_prediction_freeze'):
        raise RuntimeError('Post-freeze scoring receipt required')
    rows=payload['rows']
    case_ids=sorted(set(r['case_id'] for r in rows))
    def one(case,token):
        matches=[r for r in rows if r['case_id']==case and token in Path(r['file']).name]
        if len(matches)!=1:raise ValueError((case,token,len(matches)))
        return matches[0]
    fig,axes=plt.subplots(len(case_ids),4,figsize=(13,2.4*len(case_ids)),squeeze=False,layout='constrained')
    norm=Normalize(.5,30)
    for rownum,case in enumerate(case_ids):
        baseline=one(case,a.baseline_token);r3=one(case,a.r3_token)
        with np.load(ROOT/'private_truth'/a.split/(case+'.npz'),allow_pickle=False) as f:
            truth=f['m']
        with np.load(ROOT/'public_data'/a.split/(case+'.npz'),allow_pickle=False) as f:
            centers=f['centers']
        bpred=json.loads((ROOT/baseline['file']).read_text());rpred=json.loads((ROOT/r3['file']).read_text())
        vectors=[truth,np.asarray(bpred['candidates'][baseline['top1_index']]['m']),
                 np.asarray(rpred['candidates'][r3['top1_index']]['m']),
                 np.asarray(rpred['candidates'][r3['oracle_best_index']]['m'])]
        titles=['Ground truth',f"{a.baseline_label} top1\nRMS {baseline['top1_RMS']:.3g} | rel. residual {baseline['top1_relative_residual']:.2g}",
                f"{a.r3_label} top1\nRMS {r3['top1_RMS']:.3g} | rel. residual {r3['top1_relative_residual']:.2g}",
                f"{a.r3_label} oracle best of {r3['candidate_count']}\nRMS {r3['oracle_best_RMS']:.3g} | GT selected"]
        for col,(vector,title) in enumerate(zip(vectors,titles)):
            ax=axes[rownum,col]
            im=ax.scatter(centers[:,0],centers[:,1],c=vector,cmap='viridis',norm=norm,s=400,marker='s',edgecolor='white',linewidth=.8)
            for xy,value in zip(centers,vector):
                ax.text(*xy,f'{value:.1f}',ha='center',va='center',fontsize=7,color='white' if value<14 else 'black')
            ax.set(xlim=(-.43,.44),ylim=(-.32,.32),aspect='equal',xticks=[],yticks=[],title=title)
            ax.title.set_fontsize(9)
            if col==0:ax.set_ylabel(case)
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Material value (public box 0.5–30)',shrink=.75)
    fig.suptitle(f'{len(case_ids)} fresh blind targets; all 12 material values unknown; maximum 3 s per method\nTop1 uses all-data residual. Oracle column uses GT only after predictions froze.',fontsize=12)
    destination=ROOT/'figures';destination.mkdir(exist_ok=True)
    for ext in ('png','pdf'):fig.savefig(destination/f'{a.output_prefix}_materials.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)

    names=sorted({Path(r['file']).stem.split(r['case_id']+'_',1)[1] for r in rows})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    t=np.linspace(0,3,301)
    summaries=[]
    for method_index,name in enumerate(names):
        subset=[one(case,name) for case in case_ids]
        time_values=[]
        for r in subset:
            trace=r['online_anytime_top1']
            if not trace:
                trace=[dict(seconds=x['seconds'],top1_RMS=x['online_archive_top1_RMS']) for x in r['archive_available_by_time']]
            time_values.append(step_value([x['seconds'] for x in trace],[x['top1_RMS'] for x in trace],t))
        vals=np.asarray(time_values)
        valid=np.all(np.isfinite(vals),axis=0)
        median=np.array([np.median(v[np.isfinite(v)]) if np.any(np.isfinite(v)) else np.nan for v in vals.T])
        axes[0,0].plot(t[valid],median[valid],label=name,color=f'C{method_index}')
        for casevals in vals:axes[0,0].plot(t,casevals,alpha=.15,linewidth=.7,color=f'C{method_index}')
        maxk=max(r['candidate_count'] for r in subset)
        ks=np.arange(1,maxk+1)
        oracle=[]
        for r in subset:
            curve=r['oracle_best_at_K'];ys=[x['best_RMS'] for x in curve]
            oracle.append(step_value(np.arange(1,len(ys)+1),ys,ks))
        axes[0,1].plot(ks,np.median(oracle,axis=0),label=name,color=f'C{method_index}')
        summaries.append(dict(name=name,full_rhs=np.mean([r['cost']['full_complex_rhs'] for r in subset]),
                              reduced_rhs=np.mean([r['cost']['reduced_complex_rhs'] for r in subset]),
                              algebraic_rhs=np.mean([r['cost']['counts'].get('algebraic_complex_rhs',0) for r in subset]),
                              top1=[r['top1_RMS'] for r in subset],oracle=[r['oracle_best_RMS'] for r in subset],
                              accepted=[r['data_accepted_count'] for r in subset]))
    axes[0,0].set(xlabel='Numerical wall seconds (maximum 3 s)',ylabel='Material RMS of online residual top1',title='Online image quality: median + individual cases',xlim=(0,3))
    axes[0,1].set(xlabel='K returned candidates, residual-ranked (log scale)',ylabel='Offline oracle best material RMS',title='Coverage diagnostic; GT never selects top1',xscale='log')
    x=np.arange(len(summaries));width=.25
    axes[1,0].bar(x-width,[s['full_rhs'] for s in summaries],width,label='Full RHS (dimension 12)')
    axes[1,0].bar(x,[s['reduced_rhs'] for s in summaries],width,label='Reduced RHS (varying dimension)')
    axes[1,0].bar(x+width,[s['algebraic_rhs'] for s in summaries],width,label='Algebraic LS RHS (12x9 QR)')
    axes[1,0].set(xticks=x,xticklabels=[s['name'] for s in summaries],ylabel='Mean actual complex RHS columns per case',title='Paid solver work; no full-equivalent conversion')
    axes[1,0].tick_params(axis='x',rotation=20)
    for i,s in enumerate(summaries):
        jitter=np.linspace(-.12,.12,len(case_ids))
        axes[1,1].scatter(i+jitter,s['top1'],marker='o',color=f'C{i}',label=s['name']+' top1')
        axes[1,1].scatter(i+jitter,s['oracle'],marker='x',color=f'C{i}')
        for offset,top,oracle in zip(jitter,s['top1'],s['oracle']):axes[1,1].plot([i+offset,i+offset],[oracle,top],color=f'C{i}',alpha=.4)
    axes[1,1].set(xticks=x,xticklabels=[s['name'] for s in summaries],ylabel='Material RMS for each fresh case',title='Circles: residual top1; crosses: offline oracle')
    axes[1,1].tick_params(axis='x',rotation=20)
    for ax in axes.ravel():ax.grid(alpha=.15)
    for ax in (axes[0,0],axes[0,1],axes[1,0]):ax.legend(fontsize=7)
    fig.suptitle('Equal maximum 3 s numerical budgets; four held-out targets; all methods share public initialization\nIndependent validation, generation and all exploratory paid costs are retained separately in COST_LEDGER.json.',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(destination/f'{a.output_prefix}_cost_coverage.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    provenance='FIGURE_PROVENANCE.json' if a.output_prefix=='holdout' else a.output_prefix+'_FIGURE_PROVENANCE.json'
    (destination/provenance).write_text(json.dumps(dict(scores=str(a.scores),freeze_manifest=payload['freeze_manifest'],baseline_token=a.baseline_token,baseline_label=a.baseline_label,r3_token=a.r3_token,r3_label=a.r3_label,cases=case_ids,truth_read_for_postfreeze_plot_only=True,summaries=summaries,presentation_version=a.presentation_note),indent=2)+'\n')
    print(json.dumps(dict(figures=[f'figures/{a.output_prefix}_materials.png',f'figures/{a.output_prefix}_cost_coverage.png'],case_count=len(case_ids))))


if __name__=='__main__':main()
