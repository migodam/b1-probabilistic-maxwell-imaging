"""Compact final risk figures and fixed-first-four offline material galleries.

Run under a stage D analysis receipt, after both fresh scorings close. No
physics, model selection, checkpoint mutation, or posterior reinterpretation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write


def read(path):return json.loads((OUT/path).read_text())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-directory',type=Path,default=OUT/'figures/main')
    args=parser.parse_args()
    for path in ('L/evaluation/COMPLETE.json','L/trajectory_evaluation/COMPLETE.json'):
        if read(path)['status']!='COMPLETE':raise RuntimeError('All fresh scoring must close first')
    dest=args.output_directory
    dest.mkdir(parents=True,exist_ok=False)
    result=read('L/evaluation/RESULTS.json')
    parents=result['parents'];assert len(parents)==256
    with np.load(OUT/'L/evaluation/RAW_RISKS.npz') as z:
        names=z['method'].astype(str).tolist();risk=z['fullgrid_MSE'].copy()
    index={n:i for i,n in enumerate(names)}
    selected=['native','probit','matched_point176','direct','fused']
    # The identities are prespecified. Include the raw unsupervised prior mean
    # if its stored field uses a different explicit name; never choose by risk.
    selected=[n for n in selected if n in index]
    labels={'exactprior_material_mean':'Prior material mean','native':'Native TRF (96 cap)',
        'probit':'Strong probit point','matched_point176':'Ordinary extra16 control',
        'direct':'Measurement-only supervised','fused':'Selected local geometry CNN'}
    fig,axes=plt.subplots(1,3,figsize=(16,4.5),constrained_layout=True)
    means=[risk[:,index[n]].mean() for n in selected]
    axes[0].barh(np.arange(len(selected)),means,color=['C3' if n=='fused' else 'C0' for n in selected])
    axes[0].set(yticks=np.arange(len(selected)),yticklabels=[labels[n] for n in selected],xlabel='Full-grid complex material MSE',title='256 fresh parents; no exclusions')
    axes[0].invert_yaxis()
    axes[0].text(.98,.02,'Prior-only MSE = '+format(result['summaries']['fullgrid_MSE']['exactprior_material_mean']['mean'],'.4f')+' (table)',
        ha='right',va='bottom',transform=axes[0].transAxes,fontsize=9)
    p=risk[:,index['probit']];f=risk[:,index['fused']]
    axes[1].scatter(p,f,s=13,alpha=.7);lim=max(p.max(),f.max())
    axes[1].plot([0,lim],[0,lim],'k--',lw=1)
    axes[1].set(xlabel='Strong probit MSE',ylabel='Selected risk action MSE',title='One dot = one independent parent')
    axes[2].hist(p-f,bins=32);axes[2].axvline(0,color='black',lw=1)
    axes[2].set(xlabel='Probit MSE minus learned MSE',ylabel='Parent count',title='Both gains and failures retained')
    for ext in ('png','pdf'):fig.savefig(dest/('PRIMARY_RISK.'+ext),dpi=170)
    plt.close(fig)
    manifest=read('METHODS_FROZEN_FOR_LOCKED.json')
    families=['local_anchor_cnn','local_geometry_cnn','local_capacity_cnn','local_permuted_geometry_cnn',
        'weak_decoder','weak_strong_decoder','random_decoder','all_decoder_matched','u2_decoder']
    values=[];family_labels=[]
    for family in families:
        entries=[e for e in manifest['ablation_checkpoints'] if e['name'].startswith(family+'_n') and '_n512_' in e['name']]
        if len(entries)!=3:raise ValueError('All three frozen seeds required: '+family)
        values.append([risk[:,index['ablation__'+e['name']]].mean() for e in entries]);family_labels.append(family)
    fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True)
    for k,v in enumerate(values):
        ax.scatter([k-.16,k,k+.16],v,s=35)
        ax.plot([k-.25,k+.25],[np.mean(v)]*2,color='black',lw=2)
    ax.axhline(p.mean(),color='C3',ls='--',label='Strong probit point')
    ax.set(xticks=np.arange(len(values)),xticklabels=[s.replace('_','\n') for s in family_labels],ylabel='Full-grid complex MSE',title='Prespecified secondary controls; 256 parents, three trained seeds')
    ax.legend()
    for ext in ('png','pdf'):fig.savefig(dest/('GEOMETRY_CONTROLS.'+ext),dpi=170)
    plt.close(fig)
    cfg=read('configs/TARGET.json')
    limits=((cfg['host_chi'][0],cfg['host_chi'][0]+cfg['inclusion_increment_real_interval'][1]),
            (cfg['host_chi'][1],cfg['host_chi'][1]+cfg['inclusion_increment_imag_interval'][1]))
    galleries=[]
    for parent in parents[:4]:
        with np.load(OUT/'inputs/offline_only'/f'{parent}.npz') as z:truth=z['chi'].copy()
        with np.load(OUT/'L/estimates'/f'{parent}.npz') as z:
            images=[truth]+[z[n].copy() for n in ('native','probit','matched_point176','fused')]
        with np.load(OUT/'L/trajectory_estimates'/f'{parent}.npz') as z:
            images.append(z['dynamic_curvature__arc_length_trapezoid'].copy())
        titles=['GT (offline only)','Native point','Strong probit point','Ordinary extra16','Learned risk action','Curvature arc risk action']
        fig,axes=plt.subplots(2,6,figsize=(18,6),constrained_layout=True)
        for row,component in enumerate(('real','imag')):
            for col,(image,title) in enumerate(zip(images,titles)):
                im=axes[row,col].imshow(getattr(image,component).reshape(24,24),origin='lower',cmap='viridis',vmin=limits[row][0],vmax=limits[row][1])
                axes[row,col].set_title(title+'\n'+component);axes[row,col].set_xticks([]);axes[row,col].set_yticks([])
            fig.colorbar(im,ax=axes[row,:],shrink=.8,label='Susceptibility chi')
        fig.suptitle(parent+' | fixed first four, shared color scales | risk actions are not posterior samples')
        for ext in ('png','pdf'):fig.savefig(dest/(parent+'.'+ext),dpi=160)
        plt.close(fig);galleries.append(parent)
    refpath=OUT/'B/jobs/B_field2_direct_CDF_DEV030_r0/RESULT.json'
    if refpath.exists():
        ref=json.loads(refpath.read_text())
        # The raw refinement summary carries the exact direct-integral CDFs.
        moments=OUT/'B/jobs/B_field2_direct_CDF_DEV030_r0'
        arrays=sorted(moments.glob('GL_*_MOMENTS.npz'),key=lambda p:int(p.stem.split('_')[1]))
        if arrays:
            with np.load(moments/'FROZEN_PARTITION.npz') as z:cuts=[z['CDF_a_cut_axis0'],z['CDF_a_cut_axis1']]
            fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
            for path in arrays:
                with np.load(path) as z:
                    for axis,ax in enumerate(axes):ax.plot(cuts[axis],z[f'directly_integrated_marginal_CDF_axis{axis}'],'o-',label=path.stem.split('_')[1])
            for axis,ax in enumerate(axes):ax.set(xlabel=f'Absolute Gaussian coordinate a{axis}',ylabel='Directly integrated CDF',ylim=(0,1),title='Conditional field2 SLICE, not FULL');ax.legend(title='GL order')
            for ext in ('png','pdf'):fig.savefig(dest/('FIELD_SLICE_DIRECT_CDF.'+ext),dpi=170)
            plt.close(fig)
    write(dest/'MANIFEST.json',dict(status='COMPLETE',new_physics_RHS=0,offline_gallery_parents=galleries,
        gallery_rule='First four frozen fresh IDs, never selected by error',
        risk_unit='Independent parent; training seeds not counted as additional parents',
        files=[p.name for p in sorted(dest.iterdir()) if p.suffix in ('.png','.pdf')],
        interpretation='Point risk, trajectory risk actions and conditional-SLICE probability remain separate'))


if __name__=='__main__':main()
