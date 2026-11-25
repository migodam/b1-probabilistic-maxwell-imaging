"""Display the same observations and truth denominators for all matched methods."""
from pathlib import Path
import csv
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

p=Path(__file__).resolve().parents[1]
A=p/'experiments/A'
def records(name):
    with (A/name).open() as src:
        return list(csv.DictReader(src))
def key(row):
    return tuple(row[s] for s in ('pair_id','target_id','condition_id'))
rows=records('results/all_nodes.csv')
positive={key(v):v['status'] for v in records('psd_followup/positive_tail_nodes.csv')}
widths=[.5,.9,1.,1.5,2.,3.,5.]
background=[0.,1/32,1/16,1/8]
totals=defaultdict(lambda:np.zeros((4,7),int))
groups=[[defaultdict(int) for k in range(4)] for j in range(2)]
for v in rows:
    if v['C']!=v['D']:raise ValueError('unequal full and reduced tests')
    if v['contains_generating_material']=='True' or v['generating']=='True':continue
    fam=int(v['family']=='structure')
    i=background.index(float(v['background_cap']));j=widths.index(float(v['width_ratio']))
    totals[fam][i,j]+=1
    for col,st in enumerate((v['A'],positive[key(v)],v['C'],v['E'])):
        groups[fam][col][i,j]+=int(st=='EXCLUDED')
palette=LinearSegmentedColormap.from_list('exclusion',['#f8f6ea','#a0c9c3','#4d9395','#164359'])
fig,axes=plt.subplots(2,4,figsize=(10.25,4.15),layout='constrained',sharex=True,sharey=True)
for f in range(2):
    den=totals[f]
    for k,ax in enumerate(axes[f]):
        counts=np.array([[groups[f][k][i,j] for j in range(7)] for i in range(4)])
        val=np.divide(counts,den,out=np.full_like(counts,np.nan,dtype=float),where=den>0)
        im=ax.imshow(val,cmap=palette,vmin=0,vmax=1,aspect='auto')
        for (i,j),v in np.ndenumerate(val):
            if np.isfinite(v):
                ax.text(j,i,str(counts[i,j]),ha='center',va='center',size=8,color='white' if v>.72 else '#1b3545')
        ax.set_xticks(range(7),['.5','.9','1','1.5','2','3','5'],fontsize=8)
        ax.set_yticks(range(4),['0','1/32','1/16','1/8'],fontsize=8)
        ax.tick_params(length=0)
        if f==0:ax.set_title(['Directional perturbation','Positive source tail','Full = fixed traces','Additional shared rows'][k],size=10)
        if k==0:ax.set_ylabel(('Amplitude only' if f==0 else 'Internal structure')+'\nBackground cap',size=9)
        if f:ax.set_xlabel('Relative interval width',size=8)
bar=fig.colorbar(im,ax=axes,shrink=.87,pad=.015,aspect=26)
bar.set_label('Fraction excluded',size=9);bar.ax.tick_params(labelsize=8)
for ext in ['pdf','png']:
    fig.savefig(p/f'figures/phase_comparison.{ext}',dpi=190,bbox_inches='tight',pad_inches=.035)
