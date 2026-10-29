"""Public-geometry-only larger-support conditioning; no targets or scoring."""
import sys,json
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model

def main():
 out={'scope':'GEOMETRY_ONLY_NO_TARGET_NO_IMAGING_RUN','presets':{}}
 for n,nx,ny,rad,ylim in [(24,6,4,.04,.24),(48,8,6,.025,.30)]:
  g={'centers':[[float(x),float(y)] for y in np.linspace(-ylim,ylim,ny) for x in np.linspace(-.35,.35,nx)],'sources':[[-.65,-.70],[.72,-.48],[-.67,.42]],'receivers':[[float(.95*np.cos(2*np.pi*j/(2*n)+.071)),float(.95*np.sin(2*np.pi*j/(2*n)+.071))] for j in range(2*n)],'radius':rad}
  model=Model(None,**g);s=model.S_s;d=np.linalg.norm(model.centers[:,None]-model.centers[None,:],axis=2)+np.eye(n)*1e10
  out['presets'][f'INFO{n}_EXTERNAL']={'geometry':g,'dimensions':{'cells':n,'sources':3,'receivers':2*n},'minimum_center_spacing':float(d.min()),'disjoint_disks':bool(d.min()>2*rad),'receiver_singular_values':s.tolist(),'receiver_singular_values_relative':(s/s[0]).tolist(),'receiver_condition':float(s[0]/s[-1]),'ranks':{str(t):int(sum(s>s[0]*t)) for t in [1e-8,1e-10,1e-12,1e-14]},'setup_counts':dict(model.counts),'setup_wall_seconds':dict(model.timings)}
 (P/'experiments/larger_geometry_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:{x:v[x] for x in ('ranks','receiver_condition','disjoint_disks')} for k,v in out['presets'].items()},indent=2))
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 from matplotlib.patches import Circle
 fig,axes=plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
 for ax,(name,r) in zip(axes[:2],out['presets'].items()):
  g=r['geometry'];centers=np.array(g['centers']);rec=np.array(g['receivers']);src=np.array(g['sources'])
  for x,y in centers:ax.add_patch(Circle((x,y),g['radius'],fc='#4477aa',alpha=.8))
  ax.scatter(*rec.T,s=13,c='#ee7733',label='receivers');ax.scatter(*src.T,s=45,c='#228833',marker='*',label='sources');ax.set(xlim=(-1.05,1.05),ylim=(-1.05,1.05),aspect='equal',title=name,xlabel='x / wavelength',ylabel='y / wavelength');ax.legend(loc='lower right',fontsize=8)
 for name,r in out['presets'].items():axes[2].semilogy(np.arange(1,len(r['receiver_singular_values'])+1),r['receiver_singular_values_relative'],'o-',ms=3,label=name)
 for t in [1e-8,1e-10,1e-12]:axes[2].axhline(t,color='.5',ls=':',lw=.6)
 axes[2].set(xlabel='receiver singular index',ylabel='sigma_i / sigma_max',title='Exterior receiver conditioning');axes[2].legend(fontsize=8);axes[2].grid(alpha=.2)
 fig.suptitle('Geometry only: no material targets, no larger imaging experiment')
 fig.savefig(P/'experiments/larger_geometry_audit.svg');plt.close(fig)
if __name__=='__main__':main()
