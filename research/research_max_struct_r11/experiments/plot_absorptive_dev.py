"""Saved-material absorptive DEV plate only; no operator/solver imports."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[k]='1'
from pathlib import Path
import json,csv,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).absolute().parents[1]
R10=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'
BASE=ROOT/'delegated/diagnostics/absorptive_dev'
OUT=ROOT/'FIGURES/absorptive_dev'
CASES=('r10_twod_0001','r10_twod_0004','r10_twod_0022','r11_bothblocks16_0000')
METHODS=('GT','CIE','Original G19','Real BothBlocks','Absorptive Green','Absorptive K0')
def main():
 wall=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True)
 scores=json.loads((BASE/'COMPARISON_WITH_EXISTING_REFERENCES.json').read_text())
 new=json.loads((BASE/'ALL8_PREDICTIONS_SEALED_BEFORE_GT.json').read_text())
 oldseal=json.loads((ROOT/'delegated/diagnostics/bothblocks_confirmation/ALL80_TERMINAL_PREDICTION_SEAL.json').read_text())
 oldpaths={(e['task']['case_id'],e['task']['method']):Path(e['prediction_path']) for e in oldseal['entries']}
 newpaths={(e['task']['case_id'],e['task']['material_readout']['label']):Path(e['prediction']) for e in new['entries']}
 rows=[];fig,axs=plt.subplots(4,6,figsize=(15.2,11),constrained_layout=True)
 for i,cid in enumerate(CASES):
  isold=cid.startswith('r10_');gt=(R10/'sealed_truth/twod' if isold else ROOT/'sealed_truth/bothblocks16')/(cid+'.npz')
  paths=[gt]
  if isold:
   paths += [R10/'RESULTS/twod'/f'{cid}__{m}__noise_ball_99/prediction.npz' for m in ('cie_single','matrix_green19')]
   paths += [ROOT/'RESULTS/exposed_dev'/f'{cid}__fullGreenBothBlocks__480CPU/R11_CAP_PREDICTION.npz']
   labels=('GT','cie_single','matrix_green19','real_GreenBothBlocks','absorptiveGreenBothBlocks','absorptiveScalarBothBlocks')
  else:
   paths += [oldpaths.get((cid,m)) for m in ('cie_single','fullG19','fullGreenBothBlocks')]
   labels=('GT','cie_single','fullG19','fullGreenBothBlocks','absorptiveGreenBothBlocks','absorptiveScalarBothBlocks')
  paths += [newpaths.get((cid,m)) for m in ('absorptiveGreenBothBlocks','absorptiveScalarBothBlocks')]
  for j,(name,label,p) in enumerate(zip(METHODS,labels,paths)):
   ax=axs[i,j];metric=next((r for r in scores if r['case_id']==cid and r['method']==label),None)
   exists=p is not None and p.is_file();row=dict(case_id=cid,column=name,method=label,raw_material_path=str(p) if p else None,exists=exists,amplitude_NRMSE=metric['amplitude_NRMSE'] if metric else None,material_success=metric['material_success'] if metric else None,CPU_comparison_scope='Eight new absorptive arms: same Mac; old references may cross hosts, CPU not horizontally comparable',exposure='EXPOSED_DEV_NOT_CONFIRMATION');rows.append(row)
   if exists:
    with np.load(p,allow_pickle=False) as z:m=z['m'].copy()
    im=ax.imshow(m.reshape(24,24),origin='lower',vmin=0,vmax=20,cmap='viridis',interpolation='nearest')
   else:
    ax.set_facecolor('#eeeeee');ax.text(.5,.5,'REFERENCE MISSING\nNO RERUN',ha='center',va='center',transform=ax.transAxes)
   title=name+('\nGT' if label=='GT' else '\nNRMSE '+f"{metric['amplitude_NRMSE']:.3f}" if metric else '\nUnavailable')
   ax.set_title(title,fontsize=10,color='#9f2020' if metric and not metric['material_success'] else 'black');ax.set_xticks([]);ax.set_yticks([])
   if j==0:ax.set_ylabel(cid.replace('r10_twod_','R10 ').replace('r11_bothblocks16_','16 DEV '),fontsize=10)
 fig.colorbar(im,ax=axs,shrink=.75,label='Real material contrast chi: shared 0–20 scale')
 fig.suptitle('Four exposed DEV objects | sealed material endpoints, no new solver\nOld references may cross hosts: CPU cannot be compared across columns; eight new absorptive runs share one Mac',fontsize=13)
 fig.supxlabel('Amplitude NRMSE labels; red = inherited material gate FAIL. Fixed absorptive phase, no tuning. Exposed DEV, not confirmation.',fontsize=10)
 files=[]
 for ext in ('png','pdf'):
  p=OUT/f'ABSORPTIVE_DEV_4x6.{ext}';fig.savefig(p,dpi=180);files.append(str(p.relative_to(ROOT)))
 plt.close(fig)
 (OUT/'MATERIAL_SOURCE_TABLE.json').write_text(json.dumps(rows,indent=2)+'\n')
 with (OUT/'MATERIAL_SOURCE_TABLE.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 rec=dict(status='COMPLETE',files=files,cases=4,columns=6,missing_references=sum(not r['exists'] for r in rows),display_range=[0,20],source_table='MATERIAL_SOURCE_TABLE.json',new_solver=0,new_FW=0,new_adjoint=0,GT48_reads=0,GT_scope='Only4 already exposed authorized DEV GT after all8 predictions sealed',CPU_scope='Absolute process_time near write, nested in external whole-child render measurement',CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall)
 (OUT/'PLOT_RECEIPT.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
if __name__=='__main__':main()
