"""Counts of retained centres at three paired noise radii."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
P=Path(__file__).resolve().parents[1]
conditions=[(s,b,n) for s in (32,96) for b in ('vacuum','weak') for n in (.01,.03,.1)]
rows=[r for f in sorted((P/'maps/results/stress').glob('stress_*.json')) for r in json.loads(f.read_text())['rows'].values()]
fig,axes=plt.subplots(1,2,figsize=(7.05,2.55),layout='constrained')
for s,ax in zip((32,96),axes):
    arr=np.zeros((2,3),int)
    for i,b in enumerate(('vacuum','weak')):
        for j,n in enumerate((.01,.03,.1)):
            selected=[r for r in rows if r['n_centres']==s*s and r['actual_scene']==b and np.isclose(r['noise'],n)]
            arr[i,j]=sum(r['n_centres']-r['excluded'] for r in selected)
    meshes=ax.imshow(arr,norm=LogNorm(vmin=12,vmax=380),cmap='YlGnBu',aspect='auto')
    for (i,j),value in np.ndenumerate(arr):
        ax.text(j,i,str(value),ha='center',va='center',color='white' if value>=115 else '#183347',fontsize=12)
    ax.set_xticks([0,1,2],['1%','3%','10%']);ax.set_xlabel('Receiver noise radius')
    ax.set_yticks([0,1],['Vacuum','[0, 1/32]'])
    ax.set_title('Four 32 by 32 targets' if s==32 else 'Two 96 by 96 targets',fontsize=10,loc='left')
    ax.tick_params(length=0)
fig.savefig(P/'figures/noise_phase.pdf',bbox_inches='tight')
fig.savefig(P/'figures/noise_phase.png',dpi=165,bbox_inches='tight')
