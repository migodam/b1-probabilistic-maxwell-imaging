"""Plot the largest tested sector threshold surviving four reference/blank choices."""
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

HERE=Path(__file__).resolve().parent
TARGETS=('nylon_large','nylon_small','nylon_pair','polyethylene_e')
LABELS=('Large nylon cylinder','Small nylon cylinder','Two nylon cylinders','Polyethylene E')


def main():
    with (HERE/'transfer_margins.csv').open() as source:
        rows=[r for r in csv.DictReader(source) if r['model']=='passive_sector' and r['positive']=='True']
    frequencies=np.arange(3,6.1,.5)
    cap=np.zeros((4,7)); remainder=np.zeros_like(cap)
    for r in rows:
        i=TARGETS.index(r['target']);j=list(frequencies).index(float(r['GHz']))
        value=float(r['cap'])
        if value>cap[i,j]:
            cap[i,j]=value;remainder[i,j]=100*float(r['relative_margin_lower'])
    print(cap.tolist())
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'font.size':10})
    fig,ax=plt.subplots(figsize=(8.2,3.25))
    colors=ListedColormap(['#ecefee','#d1e6e9','#9ec9cf','#5e9a9f','#2b6269'])
    ax.imshow(cap,cmap=colors,norm=BoundaryNorm([-.1,.1,.35,.60,.85,1.15],5),aspect='auto')
    for i in range(4):
        for j in range(7):
            value=cap[i,j]
            if value:
                label=f'{value:.2f}\n+{remainder[i,j]:.2f}%'
            else:label='—'
            ax.text(j,i,label,ha='center',va='center',color='white' if value>=.75 else '#172a2d',fontsize=9.5)
    ax.set_yticks(range(4),LABELS);ax.set_xticks(range(7),[f'{x:g}' for x in frequencies])
    ax.set_xticks(np.arange(-.5,7,1),minor=True);ax.set_yticks(np.arange(-.5,4,1),minor=True)
    ax.grid(which='minor',color='white',lw=1.8);ax.tick_params(which='both',length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set_xlabel('Frequency (GHz)')
    ax.set_title('Largest tested sector cap with positive calibration margin',loc='left',fontsize=12,pad=10)
    fig.text(.19,.023,'Second line: remaining discrepancy as a percentage of the primary field norm. Dash: no tested cap passed.',size=8.6)
    fig.subplots_adjust(left=.20,right=.985,top=.86,bottom=.20)
    for ext in ('pdf','png'):fig.savefig(HERE/('calibration_headroom.'+ext),dpi=215)
    plt.close(fig)

if __name__=='__main__':
    main()
