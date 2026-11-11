"""Offline scoring after all twelve remote predictions and prefixes are sealed."""
from pathlib import Path
import json,csv,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'RESULTS/basin_remote'
CASES=['r10_twod_0001','r10_twod_0004','r10_twod_0022']
METHODS=['direct','multistart','alternating_projection','douglas_rachford']
def metrics(m,t):
 rmse=float(np.linalg.norm(m-t)/np.sqrt(len(t)));amp=float(t.max()-.1)
 union=(m>1)|(t>1);iou=float(np.sum((m>1)&(t>1))/max(1,union.sum()))
 return dict(RMSE=rmse,NRMSE=rmse/amp,relative_material_L2=float(np.linalg.norm(m-t)/np.linalg.norm(t)),
  fixed_superlevel_IoU=iou,descriptive_material_success=bool(rmse/amp<=.15 and iou>=.7))
def main():
 started=time.process_time();seal=json.loads((BASE/'PREFIX_SEAL.json').read_text())
 assert len({r['run'] for r in seal['records']})==12
 assert all((BASE/f'{c}__{m}'/'RECEIPT.json').exists() for c in CASES for m in METHODS)
 rows=[];fig,axes=plt.subplots(3,6,figsize=(15,8));longfig,longaxes=plt.subplots(3,6,figsize=(15,8))
 old=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'
 for i,c in enumerate(CASES):
  with np.load(old/'sealed_truth/twod'/f'{c}.npz') as z:t=z['m'].copy()
  with np.load(ROOT/'data_public/basin'/f'{c}_initial.npz') as z:m0=z['m'].copy()
  for axset in (axes,longaxes):
   axset[i,0].imshow(t.reshape(24,24),vmin=.02,vmax=20,cmap='viridis');axset[i,0].set_title(c+' GT',fontsize=9)
   axset[i,1].imshow(m0.reshape(24,24),vmin=.02,vmax=20,cmap='viridis');axset[i,1].set_title('same bad start',fontsize=9)
  for j,method in enumerate(METHODS):
   folder=BASE/f'{c}__{method}';r=json.loads((folder/'RECEIPT.json').read_text())
   with np.load(folder/'prediction.npz') as z:m=z['m'].copy()
   row=dict(case=c,method=method,view='full_endpoint',**metrics(m,t),residual=r['relative_data_residual'],CPU=r['cpu_seconds'])
   rows.append(row)
   longaxes[i,j+2].imshow(m.reshape(24,24),vmin=.02,vmax=20,cmap='viridis');longaxes[i,j+2].set_title(f'{method}\nNRMSE {row["NRMSE"]:.3f}',fontsize=9)
   with np.load(folder/'paid_prefix.npz') as z:
    for cap in [30,60,120,240]:
     key=f'm_{cap}'
     if key not in z:continue
     mm=z[key];r2=dict(case=c,method=method,view=f'paid_CPU_{cap}',**metrics(mm,t))
     rows.append(r2)
     if cap==240:
      axes[i,j+2].imshow(mm.reshape(24,24),vmin=.02,vmax=20,cmap='viridis');axes[i,j+2].set_title(f'{method}\nNRMSE {r2["NRMSE"]:.3f}',fontsize=9)
 for axset in (axes,longaxes):
  for ax in axset.flat:ax.set_xticks([]);ax.set_yticks([])
 fig.suptitle('R11 actual bad endpoints: paid 240-body-CPU candidates (offline scoring)')
 longfig.suptitle('All actual full endpoints; some exceed nominal240 CPU, not equalcost')
 out=ROOT/'FIGURES';out.mkdir(exist_ok=True)
 fig.tight_layout();fig.savefig(out/'basin_remote_240.png',dpi=180);plt.close(fig)
 longfig.tight_layout();longfig.savefig(out/'basin_remote_full.png',dpi=180);plt.close(longfig)
 (BASE/'OFFLINE_SCORES.json').write_text(json.dumps(dict(scope='EXPOSED_DIAGNOSTIC after12sealed',rows=rows,scoringCPU=time.process_time()-started),indent=2))
 print(json.dumps([r for r in rows if r['view']=='paid_CPU_240'],indent=2))
if __name__=='__main__':main()
