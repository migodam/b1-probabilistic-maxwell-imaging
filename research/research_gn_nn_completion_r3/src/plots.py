"""All-case physical-scale material panels and cost/defect charts."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def make_plots(scores_path, *, study=None, output=None):
    scores_path=Path(scores_path)
    scores=json.loads(scores_path.read_text())
    study=Path(study)if study else scores_path.parents[1]
    output=Path(output)if output else scores_path.parent/'figures'
    output.mkdir(parents=True,exist_ok=True)
    seal=json.loads(Path(scores['seal_path']).read_text())
    data=Path(seal['data_root'])
    indexed={}
    for row in scores['rows']:
        indexed.setdefault(row['obs_id'],[]).append(row)
    index=[]
    for obs_id,rows in indexed.items():
        source=rows[0]
        with np.load(data/source['truth_path'],allow_pickle=False)as z:
            truth=np.asarray(z['m']).reshape(24,24)
        anchor_path=study/'anchors'/f'{obs_id}.npz'
        anchor=None
        if anchor_path.exists():
            with np.load(anchor_path,allow_pickle=False)as z:
                anchor=np.asarray(z['m0']).reshape(24,24)
        panels=[]
        for row in rows:
            path=study/'evaluation'/row['run']/f'{obs_id}.npz'
            if not path.exists():
                panels.append((row,None,None))
                continue
            with np.load(path,allow_pickle=False)as z:
                raw=z['m_proposal_clipped'].reshape(24,24)if 'm_proposal_clipped'in z else anchor
                if 'top1'in z:
                    corrected=z['top1'].reshape(24,24)
                else:
                    corrected=None
                panels.append((row,raw,corrected))
        # Every displayed estimate shares one material color scale per case,
        # including failed/outlier proposals; no per-panel rescaling or cropping.
        finite_arrays=[a for a in (truth,anchor)if a is not None]+[a for _,raw,corrected in panels for a in (raw,corrected)if a is not None]
        vmin=min(.02,min(float(a.min())for a in finite_arrays))
        vmax=max(float(a.max())for a in finite_arrays)
        for page,start in enumerate(range(0,max(1,len(panels)),5)):
            group=panels[start:start+5]
            columns=2+len(group)
            fig,axes=plt.subplots(2,columns,figsize=(2.6*columns,5.5),squeeze=False,layout='constrained')
            image=None
            for line in range(2):
                for col,(label,array)in enumerate([('GT',truth),('GN anchor',anchor)]+[(p[0]['run'],p[1+line])for p in group]):
                    ax=axes[line,col]
                    if array is not None:
                        image=ax.imshow(array,cmap='viridis',vmin=vmin,vmax=vmax,interpolation='nearest')
                    else:
                        ax.text(.5,.5,'NO FEASIBLE\nSELECTION'if line else'NOT RUN',ha='center',va='center',transform=ax.transAxes)
                    ax.set_xticks([]);ax.set_yticks([])
                    ax.set_title(label,fontsize=9)
                    if col==0:
                        ax.set_ylabel('Raw full proposal'if line==0 else'Online selected after correction',fontsize=10)
            if image is not None:
                fig.colorbar(image,ax=axes.ravel().tolist(),label='Physical material contrast',shrink=.82)
            fig.suptitle(f'{obs_id} | {source["split"]} | noise={source["noise"]} | {source.get("distribution","ID")}\nUniform case scale [{vmin:.3g}, {vmax:.3g}]',fontsize=11)
            path=output/f'{obs_id}_page{page+1}.png'
            fig.savefig(path,dpi=150);plt.close(fig)
            index.append(dict(obs_id=obs_id,file=path.name,runs=[p[0]['run']for p in group],vmin=vmin,vmax=vmax,selection='all observations, all runs, manifest order'))
    runs=list(dict.fromkeys(r['run']for r in scores['rows']))
    palette=plt.get_cmap('tab10')
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for i,run in enumerate(runs):
        group=[r for r in scores['rows']if r['run']==run and 'wall_seconds'in r]
        color=palette(i%10)
        axes[0].scatter([r['wall_seconds']for r in group],[r['costcounts'].get('LU',0)for r in group],label=run,color=color,s=15,alpha=.65)
        axes[1].scatter([r['costcounts'].get('actual_complex_RHS',0)for r in group],[r['wall_seconds']for r in group],color=color,s=15,alpha=.65)
        selected=[r for r in group if r.get('top1')]
        axes[2].scatter([r['top1']['relative_residual']for r in selected],[r['top1']['NRMSE']for r in selected],color=color,s=15,alpha=.65)
    axes[0].axvline(10,color='k',ls='--',lw=.8);axes[0].axhline(100,color='k',ls='--',lw=.8)
    axes[0].set(xlabel='Incremental total wall seconds',ylabel='Charged LU count')
    axes[1].axvline(6000,color='k',ls='--',lw=.8);axes[1].axhline(10,color='k',ls='--',lw=.8)
    axes[1].set(xlabel='Charged complex RHS',ylabel='Incremental total wall seconds')
    axes[2].axhline(.15,color='k',ls='--',lw=.8);axes[2].axvline(.001,color='gray',ls='--',lw=.8);axes[2].axvline(.003,color='k',ls=':',lw=.8)
    axes[2].set(xlabel='Full-VIE physical data relative residual',ylabel='Material NRMSE / GT amplitude',xscale='log')
    if runs:axes[0].legend(fontsize=7)
    fig.savefig(output/'all_costs_and_material.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for i,run in enumerate(runs):
        group=[r for r in scores['rows']if r['run']==run and r.get('auxiliary')]
        if not group:continue
        color=palette(i%10)
        axes[0].scatter([max(r['auxiliary']['auxiliary_data_relative_residual'],1e-16)for r in group],[max(r['auxiliary']['constitutive_relative_defect'],1e-16)for r in group],label=run,color=color,s=16,alpha=.7)
        selected=[r for r in group if r.get('top1')]
        axes[1].scatter([max(r['auxiliary']['auxiliary_data_relative_residual'],1e-16)for r in selected],[max(r['top1']['relative_residual'],1e-16)for r in selected],color=color,s=16,alpha=.7)
    axes[0].set(xlabel='Auxiliary current data relative residual',ylabel='Auxiliary constitutive relative defect',xscale='log',yscale='log')
    axes[1].set(xlabel='Auxiliary current data relative residual',ylabel='Selected full-VIE data relative residual',xscale='log',yscale='log')
    handles,labels=axes[0].get_legend_handles_labels()
    if handles:axes[0].legend(fontsize=7)
    fig.savefig(output/'all_auxiliary_vs_physical_defects.png',dpi=160);plt.close(fig)
    (output/'FIGURE_INDEX.json').write_text(json.dumps(dict(scores=str(scores_path.resolve()),case_panels=index,no_cherrypicking=True,scale='one physical material scale per observation across all arms and raw/corrected panels'),indent=2)+'\n')
    return output
