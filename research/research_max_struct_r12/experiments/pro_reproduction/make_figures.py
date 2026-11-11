from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).parent
OUT=ROOT.parent/'figures'
OUT.mkdir(exist_ok=True)
records=json.loads((ROOT/'toy_results.json').read_text())['cases']
names=['ring_high','split_shell_high','bridge_high','disk_low']
labels=['High contrast ring','High contrast split shell','High contrast bridge','Low contrast disk']
cols=['truth','direct','alternating_projection','douglas_rachford','douglas_rachford_polished']
titles=['Ground truth','Direct physics','Alternating projections','Reflection','Reflection + 200 physics']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(4,5,figsize=(14,10.8),constrained_layout=True)
for i,name in enumerate(names):
    a=np.load(ROOT/(name+'.npz'))
    for j,key in enumerate(cols):
        ax=axs[i,j]
        im=ax.imshow(a[key].reshape(10,10),origin='lower',vmin=0,vmax=12,cmap='viridis',extent=[-.5,.5,-.5,.5])
        ax.set_xticks([]);ax.set_yticks([])
        if i==0:ax.set_title(titles[j],fontsize=11,pad=14)
        if j==0:ax.set_ylabel(labels[i],fontsize=11,labelpad=13)
        if key!='truth':
            err=np.linalg.norm(a[key]-a['truth'])/np.linalg.norm(a['truth'])
            ax.set_xlabel('Relative material L2 = %.3f'%err,fontsize=10)
fig.colorbar(im,ax=axs,location='right',shrink=.75,label='Material contrast m')
fig.suptitle('Four fixed finite-model TM cases | 10 x 10 cells, 6 sources, 24 receivers, 1% noise\nExploratory mechanism test; no CIE comparison or equal-cost performance claim',fontsize=14)
fig.savefig(OUT/'toy_materials.png',dpi=180)
fig.savefig(OUT/'toy_materials.pdf')
plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(10.5,7.3),constrained_layout=True)
for ax,name,label,rec in zip(axs.flat,names,labels,records):
    a=np.load(ROOT/(name+'.npz'))
    for key,title,color in [('alternating_projection','Alternating projections','#8f709a'),('douglas_rachford','Reflection','#006e6c')]:
        h=a[key+'_history']
        ax.semilogy(h[:,0],np.sqrt(2*h[:,1]),label=title,color=color,lw=1.6)
    ax.axhline(rec['relative_noise_bound'],color='#c87520',ls='--',lw=1.2,label='99% noise radius')
    ax.axhline(rec['methods']['direct']['relative_data_l2'],color='0.5',ls=':',lw=1,label='Direct endpoint')
    ax.set_title(label);ax.set_xlabel('Projection iterations');ax.set_ylabel('True Maxwell relative data residual')
    ax.grid(alpha=.2);ax.set_ylim(5e-3,3)
handles,labs=axs[0,0].get_legend_handles_labels()
fig.legend(handles,labs,loc='outside upper center',ncol=2,frameon=False)
fig.savefig(OUT/'toy_trajectories.png',dpi=180)
plt.close(fig)

controls=json.loads((ROOT/'toy_controls.json').read_text())['cases']
fig,axs=plt.subplots(1,2,figsize=(9.6,3.7),constrained_layout=True)
for ax,rec,label in zip(axs,controls[:2],labels[:2]):
    pts=np.asarray(rec['endpoint_chord']['path'])
    ax.plot(pts[:,0],pts[:,1],'-o',color='#365e82',ms=4)
    ax.axhline(.5,color='.6',ls='--',lw=1,label='Air local minimum')
    ax.set_title(label);ax.set_xlabel('t on air-to-reflection-endpoint segment');ax.set_ylabel('Original data objective')
    ax.grid(alpha=.2)
fig.suptitle('Posthoc endpoint paths: the reflected solver did not use these segments',fontsize=12)
fig.savefig(OUT/'toy_endpoint_barriers.png',dpi=180)
plt.close(fig)
print(str(OUT.resolve()))
