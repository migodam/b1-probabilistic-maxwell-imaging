"""Post-seal figures from frozen predictions and offline scoring only."""
import argparse, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from score_results import read, METHODS

def make(scores_path,truth_dir,out):
    scores=read(scores_path)
    if not scores.get('truth_disclosed_only_after_seal'):raise RuntimeError('Independent sealed scoring required')
    frozen=read(scores['seal']);receipts={}
    for file in frozen['receipts']:
        receipt=read(file['path']) if file.get('path') else {'task':file['task'],'status':'NOT_RUN'}
        receipts[receipt['task']['task_id']]=receipt
    out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True);truth_dir=pathlib.Path(truth_dir)
    rows=scores['rows'];lookup={(r['backend'],r['shape'],r['amplitude'],r['noise'],r['method']):r for r in rows}
    # One plate for every amplitude/noise/backend: CPU contains every digit, all three contrasts.
    for backend in ('cpu','cuda'):
        for noise in ('clean','noisy'):
            for amplitude in (.4,2.5,12.):
                for kind in ('digit','shape'):
                    labels=[f'digit{i}' for i in range(1,10)] if kind=='digit' else ['ring','two_islands','curved_bar']
                    labels=[s for s in labels if any((backend,s,amplitude,noise,m) in lookup for m in METHODS)]
                    if not labels:continue
                    fig,axes=plt.subplots(len(labels),5,figsize=(11.5,max(2.5,len(labels)*1.55)),squeeze=False)
                    for i,label in enumerate(labels):
                        reference=next(lookup[(backend,label,amplitude,noise,m)] for m in METHODS if (backend,label,amplitude,noise,m) in lookup)
                        with np.load(truth_dir/(reference['observation_id']+'.npz'),allow_pickle=False) as f:truth=f['m'].reshape(24,24)
                        for j,method in enumerate((None,)+METHODS):
                            ax=axes[i,j];row=lookup.get((backend,label,amplitude,noise,method)) if method else None
                            material=truth if method is None else (receipts[row['task_id']].get('result') or {}).get('m') if row else None
                            if material is not None:
                                im=ax.imshow(np.asarray(material).reshape(24,24),cmap='viridis',vmin=.02,vmax=amplitude+.15,origin='upper')
                                if row and row['finite_endpoint']:ax.set_xlabel(f"NRMSE {row['NRMSE']:.3g} | IoU {row['IoU_Otsu']:.2f}",fontsize=8)
                            else:ax.text(.5,.5,'FAILED / MISSING',ha='center',va='center',transform=ax.transAxes,fontsize=7)
                            ax.set_xticks([]);ax.set_yticks([])
                            if j==0:ax.set_ylabel(label,fontsize=10)
                            if i==0:ax.set_title('Truth (offline)' if method is None else method,fontsize=9)
                    fig.suptitle(f'{backend.upper()} · {noise} · contrast amplitude {amplitude:g}\nCommon material scale [.02,{amplitude+.15:g}]; Otsu segmentation uses reconstruction only',fontsize=11)
                    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(out/f'{backend}_{kind}_A{amplitude:g}_{noise}.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,backend in zip(axes,('cpu','cuda')):
        for method in METHODS:
            pairs=[r for r in scores['noise_pairs'] if r['backend']==backend and r['method']==method]
            if pairs:ax.scatter([r['clean_NRMSE'] for r in pairs],[r['noisy_NRMSE'] for r in pairs],label=method,s=22,alpha=.75)
        lim=max([max(r['clean_NRMSE'],r['noisy_NRMSE']) for r in scores['noise_pairs'] if r['backend']==backend]+[1.])
        ax.plot([0,lim],[0,lim],'k--',lw=.8);ax.set(xlabel='Clean NRMSE',ylabel='Paired noisy NRMSE',title=backend.upper());ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(out/'paired_noise.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for i,backend in enumerate(('cpu','cuda')):
        pairs=[r for r in scores['pairwise'] if r['backend']==backend];ax=axes[i,0]
        labels=[r['left']+' / '+r['right'] for r in pairs];pos=np.arange(len(pairs))
        wins=np.array([r['left_wins'] for r in pairs]);ties=np.array([r['ties'] for r in pairs]);losses=np.array([r['right_wins'] for r in pairs])
        ax.barh(pos,wins,label='Left wins');ax.barh(pos,ties,left=wins,label='Tie');ax.barh(pos,losses,left=wins+ties,label='Right wins')
        for j,r in enumerate(pairs):ax.text(wins[j]+ties[j]+losses[j]+.2,j,f"missing {r['missing_pairs']}",va='center',fontsize=7)
        ax.set_yticks(pos,labels,fontsize=7);ax.set_title(backend.upper()+' paired NRMSE outcomes');ax.legend(fontsize=7)
        ax=axes[i,1];summary=[r for r in scores['summary'] if r['backend']==backend]
        cells=[[r['method'],f"{r['n_finite']}/{r['n_planned']}",
                f"{r['LU_sum']:.0f}" if r['LU_sum'] is not None else '?',
                f"{r['actual_complex_RHS_sum']:.0f}" if r['actual_complex_RHS_sum'] is not None else '?',
                f"{r['wall_seconds_sum']:.1f}" if r['wall_seconds_sum'] is not None else '?'] for r in summary]
        ax.axis('off');table=ax.table(cellText=cells,colLabels=['Method','Finite','LU sum','RHS sum','Wall s'],loc='center',cellLoc='center');table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.8)
        ax.set_title(backend.upper()+' costs incl. failed tasks where known\nUnknown counters remain missing; no CPU/GPU speedup claim',fontsize=10)
    fig.tight_layout();fig.savefig(out/'paired_wins_and_costs.png',dpi=180);plt.close(fig)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--scores',required=True);ap.add_argument('--truth',default='private_truth');ap.add_argument('--out',default='figures')
    a=ap.parse_args();make(a.scores,a.truth,a.out)
