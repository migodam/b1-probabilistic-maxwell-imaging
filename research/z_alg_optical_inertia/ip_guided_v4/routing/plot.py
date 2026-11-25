#!/usr/bin/env python3
"""Create the display figures directly from the paired benchmark tables."""
from csv import DictReader
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def rows(name):
    with (HERE/name).open() as f:
        return list(DictReader(f))


def vec(data, name):
    return np.array([float(a[name]) if a.get(name) else np.nan for a in data])


def draw():
    native = rows('principal_28.csv')
    cie = rows('cie_diagnostics.csv')
    x = np.arange(1,len(native)+1)
    fig,axes = plt.subplots(2,2,figsize=(11,7),layout='constrained')
    ax = axes[0,0]
    for key, label, style, color in [('l1_lower_over_radius','Born cone','o','#17658c'),
                                     ('l2_cap05_distance_over_radius','Full wave, cap 0.50','^','#16804a'),
                                     ('l2_cap075_distance_over_radius','Full wave, cap 0.75','v','#a15a35')]:
        ax.plot(x,vec(native,key),style,color=color,markersize=4,label=label)
    ax.axhline(1,color='black',ls='--',lw=.9)
    ax.axvline(6.5,color='0.7',ls=':',lw=.9)
    ax.set(xlabel='Independent material (6 optfresh + 22 SRAP)',ylabel='Strict lower bound / radius',
           title='A  Model exclusions on 28 native objects')
    ax.legend(frameon=False,fontsize=8,loc='center right')
    ax = axes[0,1]
    m = vec(native,'zero_true_residual_over_radius')
    y = vec(native,'warm_true_residual_over_radius')
    ax.plot([1,11],[1,11],color='0.5',ls='--',lw=1)
    ax.scatter(m,y,c=np.where(y<m,'#17658c','#a15a35'),s=23,zorder=3)
    idx = next(i for i,r in enumerate(native) if r['object_id']=='s2d_018')
    ax.annotate('s2d_018',(m[idx],y[idx]),xytext=(-43,12),textcoords='offset points',fontsize=8)
    ax.set(xlim=(4,10.6),ylim=(4,10.6),xlabel='Zero start: true residual / radius',
           ylabel='Born warm start: true residual / radius',title='B  Common 16-call physical budget')
    ax.set_aspect('equal',anchor='C')
    for ax,object_id,letter in zip(axes[1],['optfresh_01','optfresh_05'],['C','D'],strict=True):
        for tag,color,label in [('standard','#17658c','Standard beta 8'),
                                ('equal_width_negative_swap','#a15a35','13-column negative-mode swap'),
                                ('beta_inverse_excluded_peak','#16804a','Beta from excluded cap'),
                                ('beta2_no_certificate','#78528c','Fixed beta 2')]:
            rr = [r for r in cie if r['object_id']==object_id and r['arm']==tag]
            if rr:
                ax.plot(vec(rr,'step'),vec(rr,'native_residual_over_radius'),'.-',lw=1.4,color=color,label=label)
        ax.axhline(1,color='black',ls='--',lw=.9)
        ax.set(xlabel='Sourcewise PR current steps',ylabel='True full-wave residual / radius',
               title=f'{letter}  Native Dual-CIE (32 tx, 32 rx): {object_id}',xticks=[4,8,16],ylim=(0,70))
    axes[1,1].legend(frameon=False,fontsize=7.8,loc='lower center')
    axes[1,1].text(8,41,'Two GMRES audit failures at step 8',ha='center',fontsize=8,color='0.35')
    fig.savefig(HERE/'routing_results.png',dpi=210)
    fig.savefig(HERE/'routing_results.pdf')


def draw_maps():
    from routing_run import native_suite, other_suite, read
    from score import reference_loader

    groups = [next((g,o) for g,o in native_suite() if o.object_id=='s2d_018'),
              next((g,o) for g,o in other_suite('mild_mid') if o.object_id=='guided_mildmid_04')]
    choices = [('Reference',None),('Born','born_only'),('Zero start','full_zero'),
               ('Born warm start','born_then_full')]
    fig,axes = plt.subplots(2,4,figsize=(10,5.6),layout='constrained')
    for line,((geo,obs),axs) in enumerate(zip(groups,axes,strict=True)):
        suite = 'native' if line==0 else 'mild_mid'
        if obs.id not in read(HERE/'online'/(suite+'__seal.json'))['ids']:
            raise ValueError('Reference plots require sealed predictions.')
        truth,_ = reference_loader(suite)(obs)
        path = HERE/'online'/(obs.id+'.npz')
        with np.load(path,allow_pickle=False) as pred:
            images = [truth]+[pred[name+'_material'] for _,name in choices[1:]]
        scores = read(HERE/'offline'/(obs.id+'.json'))['arms']
        max_scale = 4. if line==0 else 1.
        for ax,(label,name),values in zip(axs,choices,images,strict=True):
            latest = ax.imshow(values.reshape(geo.side,geo.side),origin='lower',vmin=0,vmax=max_scale,cmap='viridis')
            if name is None:
                ax.set_title(f'{label}\n{obs.object_id}',fontsize=9)
            else:
                ax.set_title(f'{label}\nMaterial relative L2 {scores[name]["full_relative_l2"]:.3f}',fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
        fig.colorbar(latest,ax=axs,shrink=.74,label='Contrast',location='right')
    fig.savefig(HERE/'native_maps.png',dpi=210)
    fig.savefig(HERE/'native_maps.pdf')


if __name__=='__main__':
    draw()
    draw_maps()
