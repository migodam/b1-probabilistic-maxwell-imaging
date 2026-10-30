"""Actual endpoint quality/cost plots from post-seal scores. No solver or truth access."""
import json, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=pathlib.Path(__file__).resolve().parent
rows=json.loads((P/'RESULTS_SUMMARY.json').read_text())['rows']
methods=('native_trf','reciprocal_trf','data_lbfgs','globaldual')
names=('Native TRF','Reciprocal TRF','Data adjoint','GlobalDual')
colors=('#346c9c','#d4823f','#508c65','#b6427c')
P.joinpath('figures').mkdir(exist_ok=True)

for kind in ('contrast_quality','cost_vs_image_quality'):
    fig,axes=plt.subplots(2,3,figsize=(14,8.2),squeeze=False)
    for iy,backend in enumerate(('cpu','cuda')):
        for ix,amplitude in enumerate((.4,2.5,12.)):
            ax=axes[iy,ix]
            selected=[r for r in rows if r['backend']==backend and r['amplitude']==amplitude and r['finite_endpoint']]
            ax.set_title(f'{backend.upper()} | amplitude {amplitude:g} | {len(selected)} finite endpoints',fontsize=10)
            if not selected:
                ax.text(.5,.5,'NOT_RUN / no finite endpoints',ha='center',transform=ax.transAxes)
                continue
            for j,(method,name,color) in enumerate(zip(methods,names,colors)):
                for noise,marker in (('clean','o'),('noisy','^')):
                    group=[r for r in selected if r['method']==method and r['noise']==noise]
                    y=[max(r['NRMSE'],1e-12) for r in group]
                    if kind=='contrast_quality':
                        offset=(-.13 if noise=='clean' else .13)
                        x=j+offset+np.linspace(-.065,.065,len(group))
                    else:x=[r['wall_seconds'] for r in group]
                    ax.scatter(x,y,color=color,marker=marker,s=25,alpha=.75,label=name+' '+noise)
                if kind=='contrast_quality':
                    group=[r['NRMSE'] for r in selected if r['method']==method]
                    if group:ax.plot([j-.32,j+.32],[max(float(np.median(group)),1e-12)]*2,color=color,lw=2)
            ax.axhline(.15,color='#666666',ls='--',lw=.8,label='NRMSE gate .15')
            ax.set_yscale('log');ax.grid(alpha=.18)
            ax.set_ylabel('Online top1 NRMSE = RMSE / amplitude',fontsize=9)
            if kind=='contrast_quality':
                ax.set_xticks(range(4),names,rotation=17,ha='right',fontsize=8)
            else:ax.set_xlabel('Actual solver wall seconds (shared budget ceiling 45 s)',fontsize=9)
            if iy==0 and ix==0:ax.legend(fontsize=6,ncol=2,loc='best')
    title='Paired material recovery across contrast: dots clean, triangles noisy; bars median' if kind=='contrast_quality' else 'Actual computation versus online material image quality: endpoint comparisons'
    fig.suptitle(title,fontsize=13)
    fig.text(.01,.012,'CPU:12 shape clusters. CUDA:paired 3-shape subset. Variants are correlated; no hardware speedup or statistical significance claim. Display floor 1e-12 only.',fontsize=8)
    fig.tight_layout(rect=(0,.035,1,.95))
    for ext in ('png','pdf'):fig.savefig(P/'figures'/f'{kind}.{ext}',dpi=180)
    plt.close(fig)
