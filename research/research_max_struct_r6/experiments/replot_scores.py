"""Plot-only replay of sealed R6 predictions/GT. Never scores or forwards physics."""
from pathlib import Path
from datetime import datetime, timezone
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
NAMES={'cie_single':'CIE','cie_dual':'Dual-CIE','global_dual':'Global dual','a_material':'A (material)'}
LABELS={'obs_b20d5a55677d4eca':'Historical curved bar','obs_42cb233e4f224bb3':'Confirmation curved bar 1','obs_63abafa750fc483c':'Confirmation curved bar 2'}
def read(p):return json.loads(Path(p).read_text())
def resolve(p):return Path(p) if Path(p).is_absolute() else REPO/p
def main():
 wall,cpu=time.perf_counter(),time.process_time();generated=[]
 mapping={r['case_id']:r for r in read(ROOT/'sealed_truth/scoring/historical_mapping_PRIVATE.json')+read(ROOT/'sealed_truth/confirmation/mapping_PRIVATE.json')['records']}
 for batch,seal_name in [('primary','PRIMARY'),('fft','FFT')]:
  directory=ROOT/'RESULTS'/('scoring_'+batch);scores=read(directory/'SCORES.json')['rows'];seal=read(ROOT/'audit'/(seal_name+'_PREDICTION_SEAL.json'))
  paths={r['task']['task_id']:resolve(r['task']['output_dir'])/'prediction.npz' for r in seal['tasks']}
  for stratum in ('DEV9_DIGITS','DEV6_DIFFICULT','NEW6_CONFIRMATION'):
   rows=[r for r in scores if r['stratum']==stratum]
   ids=list(dict.fromkeys(r['case_id'] for r in rows));arms=[a for a in ('cie_single','cie_dual','global_dual','a_material') if any(r['arm']==a for r in rows)]
   if not ids:continue
   # Retain complete original case inventory/order; fixed cell sizes keep titles inside columns.
   fig,axes=plt.subplots(len(ids),len(arms)+1,figsize=(3.1*(len(arms)+1)+.8,2.95*len(ids)+.8),squeeze=False)
   fig.subplots_adjust(left=.025,right=.93,bottom=.025,top=1-1.2/(2.95*len(ids)+.8),wspace=.17,hspace=.32)
   for i,ident in enumerate(ids):
    meta=mapping[ident]
    with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:m=z['m'].reshape(24,24)
    raw=next(r['label'] for r in rows if r['case_id']==ident)
    label=LABELS.get(ident,raw.replace('_',' ').title())
    axes[i,0].imshow(m,vmin=.02,vmax=20.,cmap='viridis');axes[i,0].set_title(label+' (GT)\nID …'+ident[-8:],fontsize=9.5,pad=6)
    for k,arm in enumerate(arms,1):
     row=next((r for r in rows if r['case_id']==ident and r['arm']==arm),None)
     if row and row['status']=='SEALED_PUBLIC_PREDICTION':
      with np.load(paths[row['task_id']],allow_pickle=False) as z:m=z['m'].reshape(24,24)
      axes[i,k].imshow(m,vmin=.02,vmax=20.,cmap='viridis')
      title=f"{NAMES[arm]}\nNRMSE {row['material_nrmse']:.3f} · IoU {row['otsu_iou']:.2f}\nResidual {row['relative_residual']:.3g}"
     else:title=NAMES[arm]+'\n'+str(row['status'] if row else 'MISSING')
     axes[i,k].set_title(title,fontsize=9.5,pad=6)
    for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
   fig.suptitle(f"R6 {stratum.replace('_',' ')} · {batch.upper()} · {len(ids)} objects",fontsize=12,y=.992)
   colorax=fig.add_axes([.947,.15,.012,.70]);bar=fig.colorbar(ScalarMappable(norm=Normalize(.02,20.),cmap='viridis'),cax=colorax)
   bar.set_label('Real material contrast m · shared range [0.02, 20]',fontsize=9);bar.ax.tick_params(labelsize=8);bar.set_ticks([.02,5,10,15,20])
   target=directory/(stratum+'_all_objects.png');fig.savefig(target,dpi=150,bbox_inches="tight",pad_inches=.15);plt.close(fig)
   generated.append(dict(path=str(target.relative_to(REPO)),batch=batch,stratum=stratum,objects=len(ids),methods=arms))
 receipt=dict(status='PLOT_ONLY_COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),generated=generated,normalization=[.02,20],scientific_arrays_or_scores_changed=False,extra_LU=0,extra_complex_RHS=0,source='experiments/replot_scores.py',wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
 receipt['visual_revision']='Reserve 1.2 inches above top axes; short labels and full shared color range'
 old=ROOT/'audit/REPLOT_RECEIPT.json'
 if old.exists():receipt['previous_plot_attempts']=read(old).get('previous_plot_attempts',[])+[dict(wall_seconds=read(old)['wall_seconds'],extra_LU=0,extra_complex_RHS=0,reason='Prior plot-only layout iteration retained; final uses reserved header and tight export bounding box')]
 with old.open('w') as f:json.dump(receipt,f,indent=2)
 print(json.dumps(receipt))
if __name__=='__main__':main()
