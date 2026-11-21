"""Single common-scale contact sheet per family subset for human reading after scoring."""
import argparse, csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

ROOT = Path(__file__).resolve().parent


def render(cohort, array, view, folder):
    with (ROOT/'results'/f'{cohort}_metrics.csv').open() as f:
        rows = list(csv.DictReader(f))
    by_case = {}
    for r in rows:
        if r['array']==array and r['view']==view and r['budget']=='30' and r['method'] in ('srap','dual'):
            by_case.setdefault(r['case'],{})[r['method']] = r
    data = ROOT / ('data_followup' if cohort=='followup' else 'data')
    size = 8
    cases = sorted(by_case)
    for pg in range(0,len(cases),size):
        fig = plt.figure(figsize=(9,7.8),layout='constrained')
        spec = GridSpec(1,2,figure=fig,wspace=.03)
        for col in range(2):
            group = fig.add_subfigure(spec[0,col])
            axes=group.subplots(4,3,gridspec_kw={'wspace':.02,'hspace':.16})
            for row in range(4):
                i=pg+col*4+row
                for ax in axes[row]:ax.set(xticks=[],yticks=[]);ax.spines[:].set_visible(False)
                if i>=len(cases): continue
                case=cases[i]
                with np.load(data/'truth'/(case+'.npz')) as z: truth=z['m']
                out=[truth]
                for m in ('dual','srap'):
                    with np.load(ROOT/'results'/'predictions'/f'{case}__{m}.npz') as z:out.append(z[f'm_30_{view}'])
                for ax,im,name in zip(axes[row],out,('Truth','Dual','SRAP')):
                    ax.imshow(im.reshape(32,32),origin='lower',interpolation='nearest',vmin=0,vmax=6,cmap='viridis')
                    if row==0:ax.set_title(name,fontsize=10)
                axes[row,0].set_ylabel(by_case[case]['dual']['object_id'],fontsize=9)
        fig.suptitle(f'{array}   |   30 s {view}   |   material contrast 0–6, shared scale',fontsize=12)
        fig.savefig(folder/f'{cohort}_{array}_{view}_{1+pg//size}.png',dpi=135,facecolor='white')
        plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('cohort');p.add_argument('array');p.add_argument('--view',default='light')
    p.add_argument('--output',type=Path,default=ROOT/'figures'/'contact');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    render(a.cohort,a.array,a.view,a.output)
