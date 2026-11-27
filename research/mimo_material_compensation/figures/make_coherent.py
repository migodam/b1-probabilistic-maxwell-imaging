"""Comparison of complete-class projected margins and the actual low-energy single shots."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT=Path(__file__).resolve().parent
D=OUT.parent/'environment'
M=D/'fluence_modes'
def read(path):return json.loads(path.read_text())
data=read(D/'decision_limits.json')['records']
colors=['#116a71','#ad473b']
fig,(left,right)=plt.subplots(1,2,figsize=(11.2,4.5),layout='constrained',gridspec_kw={'width_ratios':[1.12,1]})
labels=['Eight sequential pulses','One pulse, four-phase feeds','One pulse, reference pair','One pulse, full-class direction']
for i,n in enumerate((32,96)):
    rec=read(M/f'classes_{n}.json')
    stages={v['key']:v['fullclass_lower_over_2r'] for v in rec['rows']}
    orig=next(b['new_frequency_full_classes_lower_over_2r'] for b in data if b['side']==n and b['contrast_increment']==2)
    control=read(D/'coherent_readout'/f'proposal_{n}.json')
    four=control['best']['four_by_four']['scalar_gap_over_two_r']
    a=[orig,four,stages['center_1'],stages['class_1']]
    y=np.arange(4)+(-.16 if i==0 else .16)
    left.barh(y,a,height=.28,color=colors[i],label=f'{n} × {n}',zorder=2)
    for val,yy in zip(a,y):
        left.text(val+.07,yy,f'{np.floor(100*val)/100:.2f}',va='center',fontsize=8,color=colors[i])
    p=read(M/f'n{n}'/'proposal.json')
    est=read(M/f'n{n}'/'decision.json')['readings']
    lo=p['low_noisy_upper']; hi=p['high_noisy_lower']; med=(lo+hi)/2; r=p['detector_radius']
    xx=i
    right.plot([xx-.29,xx+.29],[(lo-med)/r]*2,color=colors[i],linestyle='--',linewidth=1.2)
    right.plot([xx-.29,xx+.29],[(hi-med)/r]*2,color=colors[i],linestyle='--',linewidth=1.2)
    right.fill_between([xx-.29,xx+.29],[(lo-med)/r]*2,[(hi-med)/r]*2,color=colors[i],alpha=.15)
    for z in est.values():
        val=(z['score_lower']-med)/r
        right.scatter(xx,val,marker='o' if z['label']=='low' else '^',s=70,color=colors[i],zorder=3,edgecolor='white')
left.axvline(1,linestyle=':',color='#142b34',label='Separation threshold')
left.set_yticks(np.arange(4),labels);left.invert_yaxis();left.set_xlim(0,7.2)
left.set_xlabel('Whole-class distance lower bound / (2r)')
left.set_title('A   Readout at 460 MHz, unit energy',loc='left',fontsize=12,fontweight='bold');left.legend(loc='upper right',fontsize=8)
right.axhline(0,color='#60706e',lw=.8)
right.set_xlim(-.55,1.55);right.set_xticks([0,1],['32 × 32\nEnergy 9/256','96 × 96\nEnergy 9/64'])
right.set_ylabel('Recorded score minus threshold, in units of r')
right.set_title('B   One new reading per specimen',loc='left',fontsize=12,fontweight='bold')
right.legend(handles=[Line2D([],[],marker='o',linestyle='',color='#384841',label='Low specification'),Line2D([],[],marker='^',linestyle='',color='#384841',label='High specification'),Line2D([],[],ls='--',color='#384841',label='Limits including noise')],fontsize=8,loc='best')
for ax in (left,right):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=9);ax.grid(axis='x' if ax is left else 'y',lw=.4,alpha=.3,zorder=0)
for suff in ['pdf','png']:fig.savefig(OUT/f'one_pulse.{suff}',dpi=190)
