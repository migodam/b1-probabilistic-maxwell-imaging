"""Actual-array development figures; uses only already separately scored GT."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write

def main():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    risk=OUT/'results/development_risk_first';density=OUT/'results/development_density_first'
    if not (risk/'scores_id_recovery/RESULTS.json').exists() or not (density/'scores/RESULTS.json').exists():raise ValueError('independent scoring must complete before GT visualization')
    folder=OUT/'figures/development_first'
    if folder.exists():raise FileExistsError('preserve previous figures')
    b=Book('I','plot_development_first',analysis_only=True,max_wall=180)
    status='COMPLETE';error=None
    try:
        folder.mkdir(parents=True);parents=[f'R4DEV_{i:05d}' for i in range(4)]
        names=['OFFLINE truth','Physics point','Frozen R2 CNN','Reoptimized R2 class','Exact-bin moments','Full66 Spline mean']
        images=[]
        for parent in parents:
            with np.load(OUT/'inputs/offline_only'/f'{parent}.npz') as z:truth=z['chi']
            with np.load(OUT/'results/contexts/dev'/f'{parent}.npz') as z:anchor=z['chi0']
            values=[truth,anchor]
            for mid in ['r2_frozen','local_cnn_s2026102501','bin_phase8_s2026102501']:
                with np.load(risk/'online'/f'{mid}.npz') as z:
                    index=np.flatnonzero(z['parent'].astype(str)==parent)[0];values.append(z['prediction_clipped'][index])
            with np.load(density/'online/spline_s1'/f'{parent}.npz') as z:values.append(z['material_mean'])
            images.append(values)
        fig,axes=plt.subplots(8,6,figsize=(13.8,16.8),layout='constrained')
        for i,values in enumerate(images):
            for part,(lo,hi,cmap) in enumerate([(.1,1.9,'viridis'),(.04,.22,'magma')]):
                for j,image in enumerate(values):
                    ax=axes[2*i+part,j];a=image.real if part==0 else image.imag
                    shown=ax.imshow(a.reshape(24,24),origin='lower',interpolation='nearest',vmin=lo,vmax=hi,cmap=cmap)
                    ax.set_xticks([]);ax.set_yticks([])
                    if i==0 and part==0:ax.set_title(names[j],fontsize=9)
                    if j==0:ax.set_ylabel(f'{parents[i]}\n'+('Re chi' if part==0 else 'Im chi'),fontsize=8)
                    if j:ax.text(.02,.02,f'MSE {np.mean(abs(image-values[0])**2):.4f}',transform=ax.transAxes,color='white',fontsize=7,bbox=dict(facecolor='black',alpha=.55,pad=1))
                if i==0:fig.colorbar(shown,ax=axes[part,:],shrink=.7,label='susceptibility chi')
        fig.suptitle('First four prespecified development parents | risk actions and unqualified model mean\nAll columns share material color scales; host pockets are not air',fontsize=11)
        fig.savefig(folder/'material_gallery.png',dpi=150);plt.close(fig)
        results=json.loads((risk/'scores_id_recovery/RESULTS.json').read_text())
        labels=['Physics','Frozen R2','Reoptimized R2','No geometry','Point phase8','Exact bin8','Exact bin16'];families=['anchor','r2','local_cnn','no_geometry_cnn','phase8','bin_phase8','bin_phase16'];means=[]
        for family in families:
            rows=[v['clipped']['full_mse']['mean'] for v in results['methods'].values() if v['family']==family];means.append(np.mean(rows))
        fig,ax=plt.subplots(figsize=(9.5,4.2),layout='constrained');ax.bar(labels,means,color=['#667580','#577e9b','#176d81','#96b5b8','#d8ad74','#b98b68','#c4a78c'])
        ax.set_ylabel('Full576 complex material MSE');ax.tick_params(axis='x',rotation=18);ax.set_title('128 new DEV parents; seed losses averaged, exploratory only')
        fig.savefig(folder/'risk_comparison.png',dpi=160);plt.close(fig)
        d=json.loads((density/'scores/RESULTS.json').read_text());labels=[];mse=[];energy=[]
        for name,v in d['methods'].items():
            labels.append(name.replace('_s1',''));mse.append(v['metrics']['raw_sample_mean_full_mse']['mean_all_parents']);energy.append(v['metrics']['material_energy_score']['mean_all_parents'])
        fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained');ax[0].bar(labels,mse,color='#648c9c');ax[1].bar(labels,energy,color='#9b876e')
        for a in ax:a.tick_params(axis='x',rotation=30)
        ax[0].set_ylabel('Model-sample mean material MSE');ax[1].set_ylabel('Joint material energy score');fig.suptitle('Full66 raw-y models | no posterior qualification from NLL')
        fig.savefig(folder/'density_material_scores.png',dpi=160);plt.close(fig)
        write(folder/'PROVENANCE.json',dict(parents=parents,selection='first four IDs before image reading; not best-of-K',risk_sources=str(risk),density_sources=str(density),
            GT='OFFLINE visualization after separate scoring gates',material='susceptibility chi, not epsilon',geometry_or_world_claim=False,new_RHS=0))
    except BaseException as e:status='FAILED';error=repr(e);raise
    finally:b.close(status,error)
if __name__=='__main__':main()
