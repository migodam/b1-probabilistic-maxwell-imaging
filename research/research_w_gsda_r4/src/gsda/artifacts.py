"""Render an already selected estimate; this module never opens truth files."""
from pathlib import Path
import numpy as np,json,csv

def save_outputs(out,config=None):
 out=Path(out)
 with np.load(out/'endpoint.npz',allow_pickle=False) as f:e={k:f[k] for k in f.files}
 seal=json.loads((out/'ONLINE_SEALED.json').read_text())
 np.save(out/'target_material.npy',e['target_contrast']+1j*e['target_loss'])
 np.save(out/'wall_material.npy',e['wall_eps']+1j*e['wall_loss'])
 np.savez_compressed(out/'predicted_measurement.npz',data=e['data'])
 (out/'reconstruction.json').write_text(json.dumps(seal,indent=2)+'\n')
 trace=[json.loads(s) for s in (out/'TRACE.jsonl').read_text().splitlines()] if (out/'TRACE.jsonl').exists() else []
 rows=[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in trace]
 if rows:
  with (out/'optimization_trace.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));w.writeheader();w.writerows(rows)
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 c=config or {};tx=e.get('target_x_edges',np.asarray(c.get('target_x_edges',[-.7,.7])));ty=e.get('target_y_edges',np.asarray(c.get('target_y_edges',[1.4,2.8])))
 for name,keys,extent,ranges in [('target',('target_contrast','target_loss'),[tx[0],tx[-1],ty[0],ty[-1]],[(0,8),(0,.25)]),('wall',('wall_eps','wall_loss'),[-1.5,1.5,.5,1.25],[(1,6),(0,.25)])]:
  fig,axs=plt.subplots(1,2,figsize=(9,4),layout='constrained')
  for ax,key,lim in zip(axs,keys,ranges):
   im=ax.imshow(e[key],origin='lower',extent=extent,aspect='equal',vmin=lim[0],vmax=lim[1],interpolation='nearest');fig.colorbar(im,ax=ax);ax.set_title(key)
  fig.suptitle('Reconstructed '+name+' — actual material cells')
  fig.savefig(out/('reconstructed_'+name+'.png'),dpi=160);plt.close(fig)
 return dict(outputs='selected estimate only',new_physics=0,truth_reads=0)
