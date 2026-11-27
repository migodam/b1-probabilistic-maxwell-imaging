"""Native material pairs consistent with each of the two prospective one-pulse records."""
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np


HERE=Path(__file__).resolve().parent
DATA=HERE.parent/'environment'/'fluence_modes'/'reattack'


def main():
    audit=json.loads((DATA/'dual_verification.json').read_text())['samples']
    fig=plt.figure(figsize=(8.9,5.6),layout='constrained')
    gs=fig.add_gridspec(3,3,width_ratios=(1,1,1.22),height_ratios=(1,1,.085))
    for row,(sample,resultfile) in enumerate([('unit_415','dual_unit_415'),('unit_862','continued_unit_862')]):
        d=audit[sample]
        if d['diameter_lower']!=2:
            raise ValueError('Two distinct uniform-core materials must have passed the strict audit')
        with np.load(DATA/f'{sample}_candidate.npz',allow_pickle=False) as f:
            ref=f['material_public_same_class']
        with np.load(DATA/f'{resultfile}.npz',allow_pickle=False) as f:
            alt,ix=f['material'],f['roi_indices']
        mask=np.zeros((96,96));mask.ravel()[ix]=1
        coord=(np.arange(96)+.5)*2/96-1
        for col,x in enumerate([ref,alt]):
            a=fig.add_subplot(gs[row,col])
            im=a.imshow(x.reshape(96,96),origin='lower',extent=(-1,1,-1,1),cmap='magma',
                        vmin=0,vmax=4.5,interpolation='nearest')
            a.contour(coord,coord,mask,levels=[.5],linewidths=.85,colors=['#62e8d7'])
            a.set(xticks=[-.75,0,.75],yticks=[-.75,0,.75])
            a.tick_params(labelsize=8)
            a.set_title(f'{"Public reference" if col==0 else "Refitted exterior"}  ·  core {x[ix[0]]:.0f}',fontsize=10)
            if col==0:a.set_ylabel(f'{"Record A" if row==0 else "Record B"}\ny (m)')
            if row==1:a.set_xlabel('x (m)')
        ax=fig.add_subplot(gs[row,2])
        for v,(key,label,color) in enumerate([('public_reference','Reference','#3b6983'),
                                               (resultfile,'Refit','#b25433')]):
            vals=d['materials'][key]['residual_upper_over_r']
            ax.barh(np.array([1.,0.])+(v-.5)*.26,[vals['400'],vals['460']],height=.24,
                    color=color,label=label)
        ax.axvline(1,color='#25262a',ls='--',lw=1)
        ax.set(xlim=(0,1.11),ylim=(-.4,1.4),yticks=[0,1],yticklabels=['460 MHz','400 MHz'],
               xticks=[0,.5,1.0])
        ax.tick_params(labelsize=9)
        if row==1:ax.set_xlabel('Residual upper bound / noise radius')
        if row==0:
            ax.set_title('Same recorded data',fontsize=10)
            ax.legend(ncol=2,loc='center',fontsize=8,frameon=False)
    cax=fig.add_subplot(gs[2,:2])
    cbar=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,1,2,3,4.5])
    cbar.set_label('Dielectric contrast',fontsize=9)
    for ext in ('pdf','png'):
        fig.savefig(HERE/f'native_wide_pairs.{ext}',dpi=195,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    main()
