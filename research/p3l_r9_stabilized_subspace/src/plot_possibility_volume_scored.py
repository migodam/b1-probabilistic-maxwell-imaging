"""E5 figures from completed offline_statistics ONLY. Never imports physics or truth.
Drafted while E5 runs; execute only after root completes the sealed scoring step.
"""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parents[1]
DEFAULT=P/'results/possibility_volume/offline_statistics'
LABELS={'FIXED_Q8_RAY':'Fixed weak line','DYNAMIC_Q8_RAY':'Dynamic complement','RANDOM66_FULL':'Ordinary random66'}
COLORS={'FIXED_Q8_RAY':'#176c97','DYNAMIC_Q8_RAY':'#dc7b2b','RANDOM66_FULL':'#66845c'}
ALLOWED={'SCORING_COMPLETE.json','SUMMARY.json','WIDTH_COVERAGE_COST_POINTS.json','FULL_SYSTEM_COSTS.json'}

def read_scored(folder,name):
    if name not in ALLOWED:raise RuntimeError('Only already-scored summary files are allowed')
    if folder.name!='offline_statistics' or 'private' in folder.parts:raise RuntimeError('Expected public offline_statistics directory')
    path=folder/name
    if not path.exists():raise RuntimeError('Scored output not available: '+str(path))
    return json.loads(path.read_text())

def store(fig,path,overwrite):
    if path.exists() and not overwrite:raise RuntimeError('Figure exists; explicit --overwrite required')
    fig.tight_layout();fig.savefig(path,dpi=180);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)

def run(folder,out,overwrite=False):
    folder=Path(folder).absolute();out=Path(out).absolute()
    receipt=read_scored(folder,'SCORING_COMPLETE.json')
    if receipt.get('status')!='COMPLETE':raise RuntimeError('Root offline scoring is not complete')
    s=read_scored(folder,'SUMMARY.json')
    if s.get('status')!='OFFLINE_SCORED_NO_POSTERIOR_CLAIM' or s['parent_count']!=32:raise RuntimeError('Complete32-parent scored contract required')
    points=read_scored(folder,'WIDTH_COVERAGE_COST_POINTS.json');costs=read_scored(folder,'FULL_SYSTEM_COSTS.json')
    arms=list(s['cloud_width_summary']);out.mkdir(parents=True,exist_ok=True)

    # Pointwise binomial intervals do not become simultaneous bands.
    fig,ax=plt.subplots(figsize=(8,4.8))
    for arm in arms:
        rows=sorted([r for r in s['coverage'] if r['method']==arm and abs(r['epsilon']-.075)<1e-12],key=lambda r:r['K'])
        xx=np.array([r['K'] for r in rows]);yy=np.array([r['rate'] for r in rows])
        lo=np.array([r['Clopper_Pearson95'][0] for r in rows]);hi=np.array([r['Clopper_Pearson95'][1] for r in rows])
        ax.errorbar(xx,yy,yerr=[yy-lo,hi-yy],marker='o',capsize=3,color=COLORS[arm],label=LABELS[arm])
    ax.set(xlabel='Nonanchor proposal prefix K (shared qualified anchor included)',ylabel='Material-tube coverage / all32 parents',
           ylim=(-.03,1.03),xticks=[0,1,4,8,16],title='E5 actual-member tube coverage, epsilon=0.075')
    ax.text(.02,.02,'Pointwise95% exact binomial intervals; not posterior probability',transform=ax.transAxes,fontsize=8)
    ax.grid(alpha=.2);ax.legend(fontsize=9);store(fig,out/'coverage_prefix.png',overwrite)

    specs=[('coverage_fixed_minus_dynamic','Coverage: fixed − dynamic','parent hit-rate difference'),
           ('coverage_fixed_minus_random','Coverage: fixed − random66','parent hit-rate difference'),
           ('offline_nearest_MSE_fixed_minus_dynamic','Offline nearest MSE: fixed − dynamic','MSE difference'),
           ('all_attempt_standalone_wall_fixed_minus_dynamic','All-attempt cost: fixed − dynamic','seconds difference')]
    fig,axes=plt.subplots(4,1,figsize=(8,7.4))
    for ax,(key,title,unit) in zip(axes,specs):
        row=s['primary_four_contrasts'][key];ax.axvline(0,color='#777',lw=1)
        if row.get('status')=='ESTIMATED':
            value=row['difference'];lo,hi=row['percentile_interval']
            ax.errorbar(value,0,xerr=[[value-lo],[hi-value]],fmt='o',capsize=5,color='#176c97')
            ax.set_title(title+f" ({100*row['confidence']:.2f}% paired-parent bootstrap)",fontsize=10)
        else:ax.text(.05,.5,row.get('status','NOT_ESTIMABLE'),transform=ax.transAxes);ax.set_title(title,fontsize=10)
        ax.set(yticks=[],xlabel=unit);ax.grid(axis='x',alpha=.2)
    fig.suptitle('Four predeclared contrasts; bootstrap is approximate with32 parents',fontsize=11)
    store(fig,out/'primary_contrasts.png',overwrite)

    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    for arm in arms:
        rows=[p for p in points if p['method']==arm]
        finite=[r for r in rows if r['pair_diameter'] is not None and r['nearest_MSE'] is not None]
        x=[r['pair_diameter'] for r in finite];y=[np.sqrt(r['nearest_MSE']) for r in finite]
        axes[0].scatter(x,y,s=28,alpha=.7,color=COLORS[arm],label=LABELS[arm]+f" ({len(finite)}/32 nonempty)")
        finite_cost=[r for r in rows if r['standalone_wall_seconds'] is not None]
        axes[1].scatter([r['standalone_wall_seconds'] for r in finite_cost],
                        [r['qualified_unique'] for r in finite_cost],s=28,alpha=.7,color=COLORS[arm])
    axes[0].axhline(.075,color='#888',ls='--',lw=1);axes[0].set(xlabel='Qualified-cloud pair diameter (material RMS)',
        ylabel='Offline nearest-member material RMS',title='Width versus offline coverage capacity')
    axes[1].set(xlabel='Standalone all-attempt wall / s',ylabel='Qualified unique material members',
                title='All attempts, anchor and shared overhead included')
    axes[0].legend(fontsize=8);[a.grid(alpha=.2) for a in axes]
    fig.suptitle('Empty clouds stay failures in coverage; absent nearest errors are not imputed',fontsize=10)
    store(fig,out/'width_coverage_cost.png',overwrite)

    fig,ax=plt.subplots(figsize=(8,4.6));xx=np.arange(len(arms));width=.25
    action_names=[('anchor_baseline_MSE','Anchor'),('uniform_MSE','Uniform finite cloud'),('likelihood_weighted_MSE','Likelihood-weighted finite cloud')]
    for j,(field,label) in enumerate(action_names):
        means=[];errs=[[],[]]
        for arm in arms:
            row=s['deployable_approximate_actions'][arm][field];mean=row['mean'];lo,hi=row['bootstrap95_percentile']
            means.append(mean);errs[0].append(max(0.,mean-lo));errs[1].append(max(0.,hi-mean))
        ax.bar(xx+(j-1)*width,means,width,yerr=errs,capsize=3,label=label)
    ax.set_xticks(xx,[LABELS[a] for a in arms]);ax.set(ylabel='Mean material MSE / all32 parents',
        title='Deployable approximate ensemble actions; not posterior means')
    ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2);store(fig,out/'approximate_actions.png',overwrite)

    meta=dict(source_directory=str(folder),source_files=sorted(ALLOWED),new_Maxwell_RHS=0,truth_arrays_read=False,
              parent_count=32,scope='Figures summarize already-scored frozen output only; no endpoint or model selection',
              omissions='Nearest-error scatter excludes empty clouds only visually; all32 remain in coverage and risk denominators')
    (out/'PLOT_PROVENANCE.json').write_text(json.dumps(meta,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--statistics-dir',type=Path,default=DEFAULT)
    p.add_argument('--output',type=Path,default=P/'figures/e5');p.add_argument('--overwrite',action='store_true')
    a=p.parse_args();run(a.statistics_dir,a.output,a.overwrite)
