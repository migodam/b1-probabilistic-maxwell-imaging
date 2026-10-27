"""Source-backed CDEV figures, using corrected full-anchor phase masks only."""
import os
os.environ['MPLBACKEND']='Agg'
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
TRACK=Path(__file__).resolve().parent
METHODS=('physics','strongCNN','distilled_w128_n32768','deterministic_matched')
COLORS={'physics':'#555555','strongCNN':'#c87a22','distilled_w128_n32768':'#2776ad','deterministic_matched':'#40935c'}
LABELS={'physics':'Physics anchor','strongCNN':'Frozen strong CNN','distilled_w128_n32768':'Frozen distilled','deterministic_matched':'Frozen deterministic'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=TRACK/'figures_CDEV_revision2')
    args=parser.parse_args()
    if not args.output.resolve().is_relative_to(TRACK.resolve()) or args.output.exists():raise FileExistsError('new Track C figure destination required')
    correction=json.loads((TRACK/'PHASE_CORRECTION_RECEIPT.json').read_text())
    if not correction.get('all20_principal_risk_cost_gate_rows_identical'):raise ValueError('phase correction qualification required')
    sources=[TRACK/'evaluation_CDEV/CDEV_RESULTS.json',TRACK/'evaluation_CDEV_revision1/CDEV_RESULTS.json']
    results=[json.loads(path.read_text()) for path in sources]
    if any('anchor q0 standardized latent field' not in result.get('phase_definition','') for result in results):raise ValueError('corrected common baseline-anchor phase masks only')
    rows=results[0]['all_arms']+[row for row in results[1]['all_arms'] if row['config']!='native96_probit64']
    configs=list(dict.fromkeys(row['config'] for row in rows));args.output.mkdir(parents=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,ax=plt.subplots(figsize=(10.8,6.5))
    ax.axhline(1.01,color='#b32938',linestyle='--',label='1% risk noninferiority limit')
    ax.axvline(.2,color='#b32938',linestyle=':',label='20% wall-saving limit')
    for method in METHODS:
        selected=[r for r in rows if r['method']==method and r['config']!='native96_probit64']
        ax.scatter([r['wall_saving_lower_one_sided_98_333333'] for r in selected],
                   [r['MSE_ratio_upper_one_sided_98_333333'] for r in selected],s=65,c=COLORS[method],label=LABELS[method],alpha=.85)
        for row in selected:
            if row['CDEV_pass']:
                x,y=row['wall_saving_lower_one_sided_98_333333'],row['MSE_ratio_upper_one_sided_98_333333']
                ax.scatter([x],[y],s=145,facecolors='none',edgecolors=COLORS[method],linewidths=1.8)
                text=row['config'].replace('native','N').replace('_probit','/P')
                offset=((-34,-17) if row['config']=='native48_probit64' else (7,-17)) if method=='physics' else (5,15)
                ax.annotate(text,(x,y),xytext=offset,textcoords='offset points',fontsize=8,color=COLORS[method])
    ax.set(xlabel='Lower one-sided 98.333% whole-wall saving',ylabel='Upper one-sided 98.333% MSE ratio to SAME estimator at full96/64',
           title='CDEV64: stage allocation can pass where joint truncation fails')
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x,_:f'{100*x:.0f}%'));ax.legend(fontsize=8,ncol=2,loc='upper right');ax.grid(alpha=.15)
    fig.text(.01,.015,'Development selection only. Parent-paired bootstrap20,000; seed losses averaged first. Timing measured under then-current load.',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1));fig.savefig(args.output/'RISK_COST_GATES.png');fig.savefig(args.output/'RISK_COST_GATES.pdf');plt.close(fig)
    # Decompose actual recorded Original phase costs; no inferred free cache cost.
    parents=[row['parent'] for row in json.loads((TRACK.parent/'configs/parents_CDEV.json').read_text())['parents']]
    native_wall=[];other_wall=[];native_rhs=[];other_rhs=[]
    for config in configs:
        reports=[json.loads((TRACK/'contexts/CDEV'/p/config/'REPORT.json').read_text()) for p in parents]
        nwall=np.array([r['original_native_report']['wall_seconds'] for r in reports]);whole=np.array([r['end_to_end_anchor_wall_seconds'] for r in reports])
        nrhs=np.array([8*r['original_native_report']['fresh_forward_states']+32*r['original_native_report']['receiver_J32_calls'] for r in reports])
        total=np.array([r['path_resource_cost']['RHS'] for r in reports])
        native_wall.append(nwall.mean());other_wall.append((whole-nwall).mean());native_rhs.append(nrhs.mean());other_rhs.append((total-nrhs).mean())
    names=[name.replace('adaptive_','adaptive\n').replace('native','N').replace('_probit','/P') for name in configs]
    fig,axes=plt.subplots(2,1,figsize=(11.5,7.5),sharex=True);x=np.arange(len(configs))
    for ax,a,b,label in [(axes[0],native_wall,other_wall,'Recorded wall seconds / parent'),(axes[1],native_rhs,other_rhs,'Paid complex RHS / parent')]:
        ax.bar(x,a,color='#537f9e',label='Native phase');ax.bar(x,b,bottom=a,color='#d79545',label='Probit + remaining context/report work')
        ax.set_ylabel(label);ax.grid(axis='y',alpha=.15)
    axes[0].set_title('Actual CDEV64 paired path costs (cache reuse never removes cold phase cost)');axes[0].legend(fontsize=9)
    axes[1].set_xticks(x,names,fontsize=9);fig.tight_layout();fig.savefig(args.output/'PHASE_COSTS.png');fig.savefig(args.output/'PHASE_COSTS.pdf');plt.close(fig)
    # Task phases are fixed by the original full anchor, never by truth.
    phase_configs=['native96_probit64','native24_probit16','native48_probit32','native24_probit64','native48_probit64']
    phase_names=['host_plateau','transition','inclusion_plateau'];fig,axes=plt.subplots(1,3,figsize=(14,4.8),sharey=True)
    for ax,phase in zip(axes,phase_names):
        for k,method in enumerate(('physics','distilled_w128_n32768','deterministic_matched')):
            values=[next(r for r in rows if r['config']==config and r['method']==method)['phase_risks'][phase]['MSE'] for config in phase_configs]
            ax.plot(range(len(values)),values,'o-',label=LABELS[method],color=COLORS[method])
        ax.set_title(phase.replace('_',' '));ax.set_xticks(range(len(phase_configs)),[c.replace('native','N').replace('_probit','/P') for c in phase_configs],rotation=35,ha='right');ax.set_yscale('log');ax.grid(alpha=.15)
    axes[0].set_ylabel('Complex material MSE within common phase masks');axes[-1].legend(fontsize=8)
    fig.suptitle('Corrected phase risks: masks defined by full96/64 anchor q0, common to all candidates')
    fig.tight_layout();fig.savefig(args.output/'CORRECTED_ANCHOR_PHASE_RISKS.png');fig.savefig(args.output/'CORRECTED_ANCHOR_PHASE_RISKS.pdf');plt.close(fig)
    benchmark_path=TRACK/'clean_inference_benchmark30/BENCHMARK.json'
    extra=[]
    if benchmark_path.exists():
        benchmark=json.loads(benchmark_path.read_text());names=list(benchmark['arms'])
        labels=['Original full deterministic','Ordinary low48/P64 deterministic','Ordinary low24/P64 deterministic','qKKT low24/P64 deterministic','Original full distilled']
        means=[benchmark['arms'][name]['mean_wall_seconds']*1000 for name in names]
        p95=[benchmark['arms'][name]['p95_wall_seconds']*1000 for name in names]
        fig,ax=plt.subplots(figsize=(10,4.7));y=np.arange(len(names))
        ax.barh(y,means,color=['#40935c','#7bb08b','#7bb08b','#9363a8','#2776ad']);ax.plot(p95,y,'|',markersize=15,color='#333333',label='95th percentile call, not a confidence bound')
        for n,value in enumerate(means):ax.text(value+.04,n,f'{value:.3f} ms',va='center',fontsize=9)
        ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Provider + feature + head wall milliseconds per parent');ax.set_title('Clean CPU1 inference: 30 repeats x64 parents x3 seeds; warm batch1 calls')
        ax.grid(axis='x',alpha=.15);ax.legend(fontsize=8,loc='lower right');fig.tight_layout();fig.savefig(args.output/'CLEAN_HEAD_MICROCOST.png');fig.savefig(args.output/'CLEAN_HEAD_MICROCOST.pdf');plt.close(fig);extra=['CLEAN_HEAD_MICROCOST']
    # Preserve the exact plotted source numbers in a parallel static table.
    with (args.output/'GATE_SOURCE.csv').open('w') as stream:
        keys=['config','method','MSE_ratio','MSE_ratio_upper_one_sided_98_333333','whole_wall_saving','wall_saving_lower_one_sided_98_333333','CDEV_pass']
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader();writer.writerows({k:r[k] for k in keys} for r in rows)
    metadata=dict(status='COMPLETE',source_results=[str(p) for p in sources],source_phase_correction=str(TRACK/'PHASE_CORRECTION_RECEIPT.json'),
        phase_masks='Original full96/64 anchor q0; no GT phase partition',new_RHS=0,GT_opened_by_plotter=False,
        figures=['RISK_COST_GATES','PHASE_COSTS','CORRECTED_ANCHOR_PHASE_RISKS',*extra],scope='CDEV development and recorded costs only; no D confirmation',
        clean_head_benchmark=str(benchmark_path) if benchmark_path.exists() else None)
    (args.output/'FIGURE_SOURCES.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata),flush=True)


if __name__=='__main__':main()
