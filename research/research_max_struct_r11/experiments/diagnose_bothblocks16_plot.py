"""Postscore full16 material plates and paid CPU/physics curves; no solver/FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time,csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation';FIG=ROOT/'FIGURES/bothblocks_confirmation'
METHODS=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks');LABELS=('CIE','Dual','Global','G19 material','BothBlocks');COLORS=plt.get_cmap('tab10').colors

def main():
 wall=time.perf_counter()
 if not (OUT/'SCORING_RECEIPT.json').exists() or not (OUT/'FRONTIER_SCORING_RECEIPT.json').exists():raise RuntimeError('Root-authorized sealed scoring/frontier required before GT plots')
 FIG.mkdir(parents=True,exist_ok=True);cfg=json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text());cases=cfg['cases'];scores=json.loads((OUT/'SINGLE_ARM_SCORES.json').read_text());ix={(x['case_id'],x['method']):x for x in scores};seal=json.loads((OUT/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text());sources={(x['task']['case_id'],x['task']['method']):x['prediction_path'] for x in seal['entries']};files=[]
 for regime in ('high','mixed'):
  selected=[x for x in cases if x['regime']==regime];fig,axes=plt.subplots(8,6,figsize=(17,23),layout='constrained')
  for ri,c in enumerate(selected):
   cid=c['case_id']
   with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:gt=z['m'].copy()
   mats=[gt]
   for method in METHODS:
    with np.load(sources[cid,method]) as z:mats.append(z['m'].copy())
   for j,m in enumerate(mats):
    ax=axes[ri,j];im=ax.imshow(m.reshape(24,24),origin='lower',cmap='viridis',vmin=0,vmax=20,interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
    if j==0:ax.set_title('GT',fontsize=11);ax.set_ylabel(cid[-4:]+' '+c['family'].replace('_',' ')+'\nhigh chi '+str(c['high_chi']),fontsize=10)
    else:
     x=ix[cid,METHODS[j-1]];ax.set_title(f"{LABELS[j-1]}\nN {x['amplitude_NRMSE']:.3f} | IoU {x['fixed_support_IoU']:.3f}",fontsize=10,color='black' if x['material_success'] else '#a51b1b',fontweight='bold' if j==5 else 'normal')
  fig.colorbar(im,ax=axes,shrink=.55,label='Real material contrast chi; common 0–20 scale');fig.suptitle('Independent BothBlocks16: '+regime+' objects | all5 arms frozen before generation',fontsize=15);fig.supxlabel('N = amplitude NRMSE. Red = material gate FAIL. All predictions selected by public Phi and sealed before GT scoring.',fontsize=10)
  for ext in ('png','pdf'):
   p=FIG/f'MATERIALS_{regime.upper()}.{ext}';fig.savefig(p,dpi=150);files.append(p.name)
  plt.close(fig)
 groups=defaultdict(list)
 with (OUT/'PAID_INCUMBENT_FRONTIER.csv').open() as f:
  for r in csv.DictReader(f):groups[r['case_id'],r['method']].append((float(r['CPU_seconds']),float(r['LU']),float(r['amplitude_NRMSE'])))
 terminals={(x['case_id'],x['method']):x for x in json.loads((OUT/'ACTUAL_TERMINAL_SPEND_POINTS.json').read_text())};common=json.loads((OUT/'COMMON_PHYSICS_CAP_SELECTION_BEFORE_GT.json').read_text())['rows'];caps={x['case_id']:x['common_physical_LU_cap'] for x in common};choices={(x['case_id'],x['method']):x for x in common};common_scores={(x['case_id'],x['method']):x for x in json.loads((OUT/'COMMON_PHYSICS_CAP_SCORES.json').read_text())}
 for axis in ('CPU','LU'):
  fig,axes=plt.subplots(4,4,figsize=(19,17),layout='constrained')
  for ci,c in enumerate(cases):
   ax=axes.flat[ci];cid=c['case_id']
   for mi,method in enumerate(METHODS):
    g=groups[cid,method];xidx=0 if axis=='CPU' else 1;g.sort(key=lambda p:p[xidx]);ax.plot([p[xidx] for p in g],[p[2] for p in g],color=COLORS[mi],lw=1.1,label=LABELS[mi])
    end=terminals.get((cid,method),{});ex=end.get('CPU_seconds') if axis=='CPU' else end.get('LU');ey=end.get('amplitude_NRMSE')
    if ex is not None and ey is not None:ax.scatter([ex],[ey],color=COLORS[mi],marker='x',s=38,zorder=4)
    if axis=='LU':
     point=choices.get((cid,method),{});score=common_scores.get((cid,method),{})
     if point.get('selected_LU') is not None and score.get('amplitude_NRMSE') is not None:ax.scatter([point['selected_LU']],[score['amplitude_NRMSE']],color=COLORS[mi],marker='*',s=70,zorder=5)
   ax.set_xscale('log');ax.axhline(.15,color='gray',lw=.6,ls='--');ax.grid(alpha=.2);ax.set_title(cid[-4:]+' '+c['family'].replace('_',' ')+' '+c['regime'],fontsize=10);ax.set_ylabel('Amplitude NRMSE');ax.set_xlabel('Paid whole-child CPU seconds' if axis=='CPU' else 'Paid physical Maxwell LU count')
   if axis=='LU':ax.axvline(caps[cid],color='gray',lw=.8,ls=':')
  handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=5,fontsize=11)
  note='Saved paid states; guide lines only. X = actual full child spend, eligible cap prediction. '+('Star = common-LU choice; dotted line = common cap. LU excludes auxiliary factors.' if axis=='LU' else 'CPU includes setup; 480 is allowed cap, not actual spend.')
  fig.suptitle('Independent16 paid '+axis+' frontier | all five arms on same local Mac\n'+note,fontsize=12)
  for ext in ('png','pdf'):
   p=FIG/f'PAID_{axis}_FRONTIERS.{ext}';fig.savefig(p,dpi=150);files.append(p.name)
  plt.close(fig)
 (FIG/'PLOT_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',objects=16,methods=5,files=files,new_forward=0,new_adjoint=0,new_inverse=0,GT_access='Only after root-authorized all80/portfolio/common-cap sealed scoring',near_exit_child_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall),indent=2)+'\n');print(json.dumps(dict(status='COMPLETE',files=files)))
if __name__=='__main__':main()
