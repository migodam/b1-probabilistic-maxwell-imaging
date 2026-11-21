"""Side-by-side actual array coordinates in units of the free-space wavelength."""
from pathlib import Path
import json,itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
root=Path(__file__).resolve().parent;geom=json.loads((root/'acquisition.json').read_text());fig=plt.figure(figsize=(10,4.8),layout='constrained')
corners=np.array(list(itertools.product([-.3,.3],repeat=3)));edges=[(a,b) for i,a in enumerate(corners) for b in corners[i+1:] if np.count_nonzero(a!=b)==1]
for gi,g in enumerate(['plane_grid','hemi']):
    ax=fig.add_subplot(1,2,gi+1,projection='3d')
    for aa,bb in edges:
        p=np.stack([aa,bb]);ax.plot(*p.T,c='#54696b',lw=1,alpha=.75)
    for role,color,mark in [('transmit','#bd582c','^'),('receive','#186887','o')]:
        xd=geom[g][role];xyz=np.array(xd['positions']);pol=np.array(xd['polarizations'])
        ax.scatter(*xyz.T,s=34 if role=='transmit' else 20,marker=mark,c=color,depthshade=False)
        for d in range(2):ax.quiver(*xyz.T,*pol[:,d,:].T,length=.058,color=color,linewidth=.6,arrow_length_ratio=.2,alpha=.85)
    ax.set(xlim=(-.95,.95),ylim=(-.95,.95),zlim=(-.35,1),xlabel='x / λ',ylabel='y / λ',zlabel='z / λ')
    ax.set_box_aspect((1.9,1.9,1.35));ax.view_init(23,-55);ax.set_title(['Rectangular planar array','Hemispherical array'][gi],fontsize=12)
handles=[Line2D([],[],lw=0,marker=m,c=c,label=l) for m,c,l in [('^','#bd582c','8 transmit ports × 2 dipoles'),('o','#186887','16 receive ports × 2 components')]]
fig.legend(handles=handles,loc='outside lower center',ncol=2,frameon=False)
fig.suptitle('Acquisition geometry · reconstruction cube width 0.6 λ',fontsize=13)
(root/'figures').mkdir(exist_ok=True);fig.savefig(root/'figures/acquisition.png',dpi=175);plt.close(fig)
