"""Exploratory offline material-risk mechanisms, after frozen fresh scoring.

Phase strata describe the original host/inclusion fraction, never air/holes.
All confidence intervals resample parents, never individual pixels.
"""
from pathlib import Path
import json,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior
from p3r2.core import OUT,write,save
from infer_p3r2_locked import frozen_contract


def main():
    for folder in ('evaluation','trajectory_evaluation','offline_geometry'):
        if json.loads((OUT/'L'/folder/'COMPLETE.json').read_text()).get('status')!='COMPLETE':
            raise RuntimeError('All fresh scoring and isolated geometry must finish before error attribution')
    m,parents,_,_=frozen_contract(OUT,OUT/'METHODS_FROZEN_FOR_LOCKED.json')
    dest=OUT/'L/material_risk_mechanism'
    if dest.exists():raise FileExistsError('No overwrite or post-test estimator selection')
    dest.mkdir(parents=True)
    cfg=json.loads((OUT/'configs/TARGET.json').read_text())
    points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
    prior=MaterialPrior(cfg,points);roi=np.flatnonzero(prior.mask)
    names=['probit','fused','direct','residual','classical_postprocess','matched_point176']
    phase_names=['exact_host_plateau','transition_fraction','exact_inclusion_plateau']
    counts=[];sums=[];full=[];real=[];imag=[];confusions=[];groups=[]
    error_maps=[];gain_maps=[]
    for parent in parents:
        with np.load(OUT/'inputs/offline_only'/f'{parent}.npz') as z:truth=z['chi'];qt=z['q']
        with np.load(OUT/'L/contexts'/f'{parent}.npz') as z:q0=z['q0']
        with np.load(OUT/'L/estimates'/f'{parent}.npz') as z:images=np.stack([z[name] for name in names])
        gtfrac=prior.decode(qt)[1]['fraction'][roi];basefrac=prior.decode(q0)[1]['fraction'][roi]
        def phase(x):return np.where(x<=1e-12,0,np.where(x>=1-1e-12,2,1))
        labels=phase(gtfrac);base_labels=phase(basefrac)
        e=abs(images-truth)**2;co=np.bincount(labels,minlength=3)
        counts.append(co);sums.append(np.stack([e[:,roi][:,labels==k].sum(1) for k in range(3)],axis=1))
        full.append(e.mean(1));real.append(((images-truth).real**2).mean(1));imag.append(((images-truth).imag**2).mean(1))
        confusion=np.zeros((3,3),int)
        np.add.at(confusion,(labels,base_labels),1);confusions.append(confusion)
        groups.append('mostly_host' if co[0]/256>=.95 else 'mostly_inclusion' if co[2]/256>=.95 else 'mixed')
        error_maps.append(e[[0,1]]);gain_maps.append(e[0]-e[1])
    counts=np.asarray(counts);sums=np.asarray(sums);full=np.asarray(full);real=np.asarray(real);imag=np.asarray(imag)
    idx=np.random.default_rng(2026102292).integers(0,256,(2000,256))
    phase_results={}
    for k,name in enumerate(phase_names):
        denominator=counts[:,k].sum();bootden=counts[idx,k].sum(1)
        if denominator==0:phase_results[name]=dict(status='NO_TEST_CELLS');continue
        risks=sums[:,:,k].sum(0)/denominator
        delta=(sums[:,1,k]-sums[:,0,k]);valid=bootden>0
        boot=delta[idx].sum(1)[valid]/bootden[valid]
        phase_results[name]=dict(test_cells=int(denominator),parents_with_cells=int((counts[:,k]>0).sum()),
            mean_cell_MSE={n:float(risks[j]) for j,n in enumerate(names)},
            fused_minus_probit=float(delta.sum()/denominator),
            paired_parent_CI95=np.quantile(boot,[.025,.975]).tolist(),
            scope='Exploratory phase-stratified error, not independent pixels or posterior events')
    group_results={}
    for group in ('mostly_host','mostly_inclusion','mixed'):
        mask=np.asarray(groups)==group
        group_results[group]=dict(parents=int(mask.sum()),fullgrid_MSE={n:float(full[mask,j].mean()) if mask.any() else None for j,n in enumerate(names)})
    result=dict(status='COMPLETE',interpretation='Post-scoring exploratory mechanism audit; no model/hyperparameter reselection',
        GT_used_online=False,new_RHS=0,parents=parents,methods=names,
        component_MSE={n:dict(real=float(real[:,j].mean()),imag=float(imag[:,j].mean()),complex=float(full[:,j].mean())) for j,n in enumerate(names)},
        phase_strata=phase_results,parent_groups=group_results,
        fraction_definition='Original P1 C2 decoder fraction: host/inclusion effective material, not air or topology probability',
        fraction_thresholds=dict(plateau_tolerance=1e-12,parent_phase_fraction=.95),
        GT_phase_vs_point_phase_cell_counts=np.sum(confusions,axis=0),
        improved_parents=int((full[:,1]<full[:,0]).sum()),worsened_parents=int((full[:,1]>full[:,0]).sum()),
        bootstrap=dict(draws=2000,seed=2026102292,unit='independent parent'))
    write(dest/'RESULTS.json',result)
    save(dest/'RAW.npz',parent=np.array(parents),method=np.array(names),phase_counts=counts,phase_error_sums=sums,
         fullgrid_MSE=full,real_MSE=real,imag_MSE=imag,parent_group=np.array(groups),
         mean_squared_error_maps=np.mean(error_maps,axis=0),mean_material_risk_gain_map=np.mean(gain_maps,axis=0))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    maps=np.mean(error_maps,axis=0);vmax=float(np.max(maps))
    for ax,image,title in zip(axes[0,:2],maps,['Probit squared material error','Learned squared material error']):
        im=ax.imshow(image.reshape(24,24),origin='lower',vmin=0,vmax=vmax);ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.7)
    gain=np.mean(gain_maps,axis=0);lim=float(np.max(abs(gain)))
    im=axes[0,2].imshow(gain.reshape(24,24),origin='lower',cmap='RdBu_r',vmin=-lim,vmax=lim);axes[0,2].set_title('Point error minus learned error');fig.colorbar(im,ax=axes[0,2],shrink=.7)
    for j,name in enumerate(('probit','fused')):
        axes[1,0].bar(np.arange(3)+j*.35,[phase_results[k].get('mean_cell_MSE',{}).get(name,np.nan) for k in phase_names],width=.35,label=name)
    axes[1,0].set_xticks(np.arange(3)+.175,['host','transition','inclusion']);axes[1,0].set_title('Original decoder phase strata');axes[1,0].legend()
    axes[1,1].hist(full[:,0]-full[:,1],bins=30);axes[1,1].axvline(0,color='black');axes[1,1].set_title('Per-parent MSE gain; all256')
    im=axes[1,2].imshow(np.sum(confusions,axis=0));axes[1,2].set(xlabel='Point phase',ylabel='GT phase',title='Descriptive phase counts');fig.colorbar(im,ax=axes[1,2],shrink=.7)
    for ext in ('png','pdf'):fig.savefig(dest/('MATERIAL_RISK_MECHANISM.'+ext),dpi=160)
    plt.close(fig);write(dest/'COMPLETE.json',dict(status='COMPLETE',new_RHS=0,offline_only=True))


if __name__=='__main__':main()
