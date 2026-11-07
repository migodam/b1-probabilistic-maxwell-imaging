"""All-case physical-scale material panels and cost/defect charts."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def make_plots(scores_path, *, study=None, output=None, replay_frozen=False):
    scores_path=Path(scores_path)
    scores=json.loads(scores_path.read_text())
    study=Path(study)if study else scores_path.parents[1]
    replay_frozen=bool(replay_frozen or scores.get('replay_frozen',False))
    output=Path(output)if output else scores_path.parent/'figures'
    output.mkdir(parents=True,exist_ok=True)
    from research.research_gn_nn_completion_r3.src.score import remap_path
    original_study=Path(scores.get('original_study',str(Path(scores['seal_path']).parent.parent)))
    seal_path=remap_path(scores['seal_path'],original_study,study)if replay_frozen else Path(scores['seal_path'])
    seal=json.loads(seal_path.read_text())
    if replay_frozen:
        source_study=Path(seal['study'])
        seal=dict(seal,manifest=dict(seal['manifest'],path=str(remap_path(seal['manifest']['path'],source_study,study))),
                  data_root=str(remap_path(seal['data_root'],source_study,study)))
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
        # Frozen physical scales shared by every case in each amplitude cohort.
        vmin,vmax=.02,(2.65 if source.get('level')=='stress' else .55)
        saturated=[]
        for row,raw,corrected in panels:
            for stage,array in [('raw',raw),('selected',corrected)]:
                if array is not None:
                    saturated.append(dict(run=row['run'],stage=stage,below=int(np.count_nonzero(array<vmin)),above=int(np.count_nonzero(array>vmax))))
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
                    sat=0 if array is None else int(np.count_nonzero((array<vmin)|(array>vmax)))
                    ax.set_title(label+(f' | saturated {sat}'if sat else''),fontsize=9)
                    if col==0:
                        ax.set_ylabel('Raw full proposal'if line==0 else'Online selected after correction',fontsize=10)
            if image is not None:
                fig.colorbar(image,ax=axes.ravel().tolist(),label='Physical material contrast (out-of-range colors saturated)',shrink=.82,extend='both')
            fig.suptitle(('Historical already-disclosed replay | 'if replay_frozen else'')+f'{obs_id} | {source["split"]} | noise={source["noise"]} | {source.get("distribution","ID")}\nFixed amplitude-cohort scale [{vmin:.3g}, {vmax:.3g}]',fontsize=11)
            path=output/f'{obs_id}_page{page+1}.png'
            fig.savefig(path,dpi=150);plt.close(fig)
            index.append(dict(obs_id=obs_id,file=path.name,runs=[p[0]['run']for p in group],vmin=vmin,vmax=vmax,selection='all observations, all runs, manifest order',saturation=saturated))
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
    axes[0].set(xlabel='Recorded seconds (weakreg45 includes its own anchor)',ylabel='Charged LU count')
    axes[1].axvline(6000,color='k',ls='--',lw=.8);axes[1].axhline(10,color='k',ls='--',lw=.8)
    axes[1].set(xlabel='Charged complex RHS',ylabel='Recorded method seconds; cost scopes differ')
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
    (output/'FIGURE_INDEX.json').write_text(json.dumps(dict(scores=str(scores_path.resolve()),case_panels=index,no_cherrypicking=True,replay_frozen=replay_frozen,truth_status='already disclosed historical replay'if replay_frozen else'original sealed scoring',scale='fixed across every low-amplitude case [.02,.55], every stress case [.02,2.65]; saturation declared'),indent=2)+'\n')
    _compact_manifest_gallery(scores,study,data,output,seal=seal,replay_frozen=replay_frozen)
    _current_mechanism_figure(scores,output)
    from html import escape
    gallery=['<!doctype html><meta charset="utf-8"><title>All case material gallery</title>', '<style>body{font:16px system-ui;margin:24px;background:#f5f6f8}img{width:100%;height:auto;background:white}section{margin-bottom:24px}h1{font-size:24px}</style>', '<h1>All case material gallery</h1><p>Manifest order, all cases and runs; fixed physical scales low [.02,.55], stress [.02,2.65]. Saturation is declared in each panel. No performance-based selection.</p>']
    for case in index:
        gallery.append(f'<section><h2>{escape(case["obs_id"])}</h2><img loading="lazy" src="{escape(case["file"])}" alt="All arm raw and selected material maps"></section>')
    (output/'GALLERY.html').write_text('\n'.join(gallery)+'\n')
    return output


def _compact_manifest_gallery(scores,study,data,output,*,seal,replay_frozen=False):
    """Choose first clean low case per digit/family in original manifest order."""
    by_obs={}
    for row in scores['rows']:
        by_obs.setdefault(row['obs_id'],[]).append(row)
    manifest=json.loads(Path(seal['manifest']['path']).read_text())
    chosen={}
    order=[f'digit{d}'for d in range(1,10)]+['primitive_composites','islands','smooth']
    for record in manifest['records']:
        if record['obs_id']not in by_obs or record.get('level','low')!='low' or float(record['noise'])!=0:
            continue
        key=f'digit{record["digit"]}'if record.get('family')=='digits'else record.get('family')
        if key in order and key not in chosen:
            chosen[key]=record
    selected=[(key,chosen[key])for key in order if key in chosen]
    runs=[]
    for row in scores['rows']:
        run=row['run']
        if run not in runs and any(a in run for a in ('v_image','v_physics','j_full','j_null')):
            runs.append(run)
    # First available run per declared arm, from the frozen seal's run order.
    arms=('v_image','v_physics','j_full','j_null','j_null_current')
    run_columns=[]
    for arm in arms:
        matches=[r for r in runs if (f'_{arm}_' in f'_{r}_' or r==arm)]
        if arm=='j_null':matches=[r for r in matches if 'j_null_current'not in r]
        if matches:run_columns.append(matches[0])
    if not selected:return
    fig,axes=plt.subplots(len(selected),2+len(run_columns),figsize=(2.0*(2+len(run_columns)),1.8*len(selected)),squeeze=False,layout='constrained')
    for line,(key,record)in enumerate(selected):
        obs=record['obs_id']
        with np.load(data/record['truth_path'],allow_pickle=False)as z:truth=z['m'].reshape(24,24)
        anchor_path=study/'anchors'/f'{obs}.npz'
        anchor=None
        if anchor_path.exists():
            with np.load(anchor_path,allow_pickle=False)as z:anchor=z['m0'].reshape(24,24)
        panels=[('GT',truth),('GN',anchor)]
        for run in run_columns:
            path=study/'evaluation'/run/f'{obs}.npz';array=None
            if path.exists():
                with np.load(path,allow_pickle=False)as z:
                    if 'top1'in z:array=z['top1'].reshape(24,24)
            panels.append((run,array))
        for col,(title,array)in enumerate(panels):
            ax=axes[line,col];ax.set_xticks([]);ax.set_yticks([])
            if array is not None:
                im=ax.imshow(array,cmap='viridis',vmin=.02,vmax=.55,interpolation='nearest')
                sat=int(np.count_nonzero((array<.02)|(array>.55)))
                if sat:ax.text(.02,.02,f'sat:{sat}',transform=ax.transAxes,fontsize=7,color='white',backgroundcolor='black')
            else:ax.text(.5,.5,'NOT RUN / FAIL',transform=ax.transAxes,ha='center',va='center',fontsize=7)
            if line==0:ax.set_title(title,fontsize=8)
            if col==0:ax.set_ylabel(key,fontsize=9)
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Physical material contrast; fixed low cohort [.02,.55]',shrink=.5,extend='both')
    fig.suptitle(('Historical already-disclosed replay | 'if replay_frozen else'')+'Preselected first clean low case per digit/family, manifest order; online selected materials',fontsize=12)
    fig.savefig(output/'preselected_digits1to9_and_nondigits.png',dpi=150);plt.close(fig)
    (output/'PRESELECTED_MANIFEST.json').write_text(json.dumps(dict(rule='first clean low observation per digit1..9 and each nondigit family, in original manifest order, no metric-based selection',records=[dict(selection_key=k,**r)for k,r in selected],runs=run_columns,range=[.02,.55],replay_frozen=replay_frozen),indent=2)+'\n')


def _current_mechanism_figure(scores,output):
    fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    runs=list(dict.fromkeys(r['run']for r in scores['rows']if r.get('current_metrics')))
    palette=plt.get_cmap('tab10')
    for index,run in enumerate(runs):
        rows=[r for r in scores['rows']if r['run']==run and r.get('current_metrics')]
        color=palette(index%10)
        x=[max(r['current_metrics']['auxiliary_current_relative_L2_error'],1e-16)for r in rows]
        qualified=[r for r in rows if r['current_metrics'].get('constitutive_relative_defect')is not None]
        axes[0].scatter([max(r['current_metrics']['auxiliary_current_relative_L2_error'],1e-16)for r in qualified],[max(r['current_metrics']['constitutive_relative_defect'],1e-16)for r in qualified],label=run,color=color,s=15,alpha=.6)
        raw=[r for r in rows if r['current_metrics'].get('raw_full_material_NRMSE')is not None]
        axes[1].scatter([max(r['current_metrics']['auxiliary_current_relative_L2_error'],1e-16)for r in raw],[max(r['current_metrics']['raw_full_material_NRMSE'],1e-16)for r in raw],color=color,s=15,alpha=.6)
        selected=[r for r in rows if r['current_metrics'].get('selected_material_NRMSE')is not None]
        axes[2].scatter([max(r['current_metrics']['projected_target_current_error_relative'],1e-16)for r in selected],[max(r['current_metrics']['selected_material_NRMSE'],1e-16)for r in selected],color=color,s=15,alpha=.6)
    axes[0].set(xlabel='Auxiliary current relative L2 error vs GT',ylabel='Constitutive relative defect',xscale='log',yscale='log')
    axes[1].set(xlabel='Auxiliary current relative L2 error vs GT',ylabel='Raw full material NRMSE / amplitude',xscale='log',yscale='log')
    axes[2].set(xlabel='Increment error relative to P0(jGT−j0)',ylabel='Selected material NRMSE / amplitude',xscale='log',yscale='log')
    if runs:axes[0].legend(fontsize=7)
    fig.savefig(output/'current_error_vs_constitutive_and_material.png',dpi=160);plt.close(fig)
