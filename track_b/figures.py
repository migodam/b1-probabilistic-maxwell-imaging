import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run import OUT,SEEDS
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 s=json.loads((OUT/'SUMMARY.json').read_text());rows=s['records'];fig,ax=plt.subplots(figsize=(10,4.8));x=np.arange(len(rows));width=.35
 ax.bar(x-width/2,[r['mean_frozen'] for r in rows],width,label='Frozen head');ax.bar(x+width/2,[r['mean_matched'] for r in rows],width,label='Matched head')
 ax.axhline(.00746135,color='gray',ls=':',label='Historical strong CNN');ax.set_xticks(x,[r['arm'].replace('_','\n') for r in rows],fontsize=9);ax.set_ylabel('VAL256 mean seed full complex MSE');ax.set_ylim(.0065,.009);ax.legend(fontsize=9);ax.set_title('R6 finite-task training: all development/selection data, no fresh test');fig.tight_layout();fig.savefig(OUT/'RISK_COMPARISON.png',dpi=180);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(10,4))
 for ax,arm in zip(axs,['direct_gt','gt_action']):
  for seed in SEEDS:
   loss=[json.loads(v) for v in (OUT/'models'/f'{arm}_s{seed}'/'provider_losses.jsonl').read_text().splitlines()];u=[r['update'] for r in loss];ax.plot(u,[r['validation_material_component_MSE']*8/9 for r in loss],label=f'VAL seed{seed%100}');train=[r for r in loss if r['update']];ax.plot([r['update'] for r in train],[r['train_GT_action_MSE']*8/9 for r in train],ls='--',alpha=.5)
  ax.set_title(arm);ax.set_xlabel('Provider update');ax.set_ylabel('Full-grid equivalent complex MSE');ax.legend(fontsize=8)
 fig.suptitle('Solid=VAL; dashed=training minibatch; early selected providers preserved');fig.tight_layout();fig.savefig(OUT/'GT_ACTION_LEARNING_CURVES.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
