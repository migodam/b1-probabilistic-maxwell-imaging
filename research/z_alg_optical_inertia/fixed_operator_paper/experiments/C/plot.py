"""Render the prospective menu and paired-policy comparisons from receipts."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
BATCH=(1,2,4)
POLICIES=('random','uniform_gap','fisher_nuisance','conditional_ellipsoid')
LABELS=('Random','Uniform gap','Fisher with nuisance','Conditional ellipsoid')
BLUE='#246e88'
ORANGE='#a74825'


def loaded(path):
    return json.loads(path.read_text())


def menu(plan,contract):
    rows=[x for x in plan['source_rows'] if x['energy_contract']==contract]
    k=[sum(row['batches'][str(m)]['strict_separating_actions'] for row in rows) for m in BATCH]
    n=[sum(row['batches'][str(m)]['finite_menu_size'] for row in rows) for m in BATCH]
    return np.array(k), np.array(n)


def render(out):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'axes.labelsize':9})
    fig,axes=plt.subplots(1,3,figsize=(11.9,3.85),width_ratios=(1.4,1.,1.))
    primary=loaded(ROOT/'proposal.json')
    diagnostic=loaded(ROOT/'diagnostic/proposal.json')
    result=loaded(ROOT/'diagnostic/results.json')['by_policy_contract_batch']
    ax=axes[0]
    x=np.arange(len(BATCH))
    for i,(contract,display,color) in enumerate((('added','Added energy',BLUE),('fixed_total','Fixed total energy',ORANGE))):
        p,n=menu(primary,contract)
        v,_=menu(diagnostic,contract)
        delta=(-1 if i==0 else 1)*.08
        ax.plot(x+delta,100*v/n,'o-',color=color,lw=1.6,ms=5,label=display)
        for j,(a,b) in enumerate(zip(v,n,strict=True)):
            if j==0 and i==1:
                continue
            offsets=((12,9),(0,9),(5,-12)) if i==0 else ((0,0),(-12,-13),(0,-13))
            ax.annotate(f'{a}/{b}',(x[j]+delta,100*a/b),xytext=offsets[j],
                        textcoords='offset points',ha='center',fontsize=7.5,color=color)
        assert np.all(p==n)
    ax.axhline(100,color='#5a636f',lw=1,ls=(0,(3,2)),label='Original experiment')
    ax.set(title='A   Strict pre-acquisition edges',ylabel='Share of the finite action menu (%)',
           xlabel='New beams',xticks=x,xticklabels=BATCH,ylim=(15,108),yticks=(20,40,60,80,100))
    ax.grid(axis='y',color='#e6e9ec',lw=.7)
    ax.legend(loc='lower left',frameon=False,fontsize=7.7)

    for ax,metric,title in zip(axes[1:],('strict_predicted_edge','wrong_class_excluded'),
                               ('B   Planned strict edges','C   Wrong class excluded'),strict=True):
        conds=[(c,m) for c in ('added','fixed_total') for m in BATCH]
        grid=np.array([[result[f'{name}|{c}|{m}'][metric] for c,m in conds] for name in POLICIES])
        ax.imshow(grid,cmap='Blues',vmin=0,vmax=8,aspect='auto')
        for (i,j),k in np.ndenumerate(grid):
            ax.text(j,i,f'{k}/6',ha='center',va='center',color='white' if k>=5 else '#182431',fontsize=8.7)
        ax.axvline(2.5,color='white',lw=2)
        ax.set(title=title,xticks=range(6),xticklabels=BATCH+BATCH,yticks=range(4),
               yticklabels=LABELS if metric=='strict_predicted_edge' else [])
        ax.tick_params(length=0,pad=4)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.text(1,4.32,'Added energy',ha='center',va='top',fontsize=8.5,color=BLUE)
        ax.text(4,4.32,'Fixed total energy',ha='center',va='top',fontsize=8.5,color=ORANGE)
        ax.set_ylim(4.7,-.5)
    fig.subplots_adjust(left=.083,right=.994,bottom=.17,top=.85,wspace=.65)
    fig.text(.083,.938,'Prospective active illumination in six outgoing-TM targets',weight='bold',fontsize=13)
    fig.text(.083,.047,'Solid curves and policy panels show the paired fourfold-wider-class, 2% old-signal noise setting. '
             'Each policy cell counts physical targets.',fontsize=8,color='#344150')
    fig.savefig(out.with_suffix('.pdf'),bbox_inches='tight',pad_inches=.16)
    fig.savefig(out.with_suffix('.png'),dpi=240,bbox_inches='tight',pad_inches=.16)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'active_policy')
    render(parser.parse_args().output)
