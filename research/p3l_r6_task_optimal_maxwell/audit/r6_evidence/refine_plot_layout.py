from pathlib import Path
import csv,json,hashlib,shutil,os
import numpy as np
ROOT=Path.cwd();F=ROOT/'track_d/figures';S=F/'sources';old=ROOT/'audit/r6_evidence/VISUAL_QA_INITIAL';old.mkdir(exist_ok=True)
for name in ['01_RISK_COST_PARETO','04_ANCHOR_PHASE_RISK','06_SEED_DISAGREEMENT']:
 for suffix in ['.png','.pdf']:shutil.copy2(F/(name+suffix),old/(name+suffix))
shutil.copy2(F/'FIGURE_MANIFEST.json',old/'FIGURE_MANIFEST.json');shutil.copy2(S/'plot_results.py',old/'plot_results.py')
os.environ['MPLCONFIGDIR']=str(F/'mpl_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42})
labels={'full_deterministic':'Full Det','selected_deterministic':'Low KKT-q','low_actual_deterministic':'Low actual','full_kkt_deterministic':'Full KKT-q','full_distill':'Full Distill','full_strong':'Full Strong'}
colors={'full_deterministic':'#334155','selected_deterministic':'#007f73','low_actual_deterministic':'#cf7b19','full_kkt_deterministic':'#6871bd','full_distill':'#c0517d','full_strong':'#8a949e'}
methods=list(labels)
def read(name):return list(csv.DictReader((S/name).open()))
def finish(fig,name):
 fig.savefig(F/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(F/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
rows=read('RISK_COST.csv');fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for ax,group in zip(axes,['DLOW','DHIGH']):
 pts=[r for r in rows if r['group']==group]
 for r in pts:
  m=r['method'];ax.scatter(float(r['mean_anchor_plus_learning_seconds']),float(r['mean_seed_loss_MSE']),s=52,color=colors[m],marker='*' if m=='selected_deterministic' else 'o',label=labels[m],zorder=3)
 front=sorted([r for r in pts if not any(float(t['mean_seed_loss_MSE'])<=float(r['mean_seed_loss_MSE']) and float(t['mean_anchor_plus_learning_seconds'])<=float(r['mean_anchor_plus_learning_seconds']) and (float(t['mean_seed_loss_MSE'])<float(r['mean_seed_loss_MSE']) or float(t['mean_anchor_plus_learning_seconds'])<float(r['mean_anchor_plus_learning_seconds'])) for t in pts)],key=lambda r:float(r['mean_anchor_plus_learning_seconds']))
 ax.plot([float(r['mean_anchor_plus_learning_seconds']) for r in front],[float(r['mean_seed_loss_MSE']) for r in front],ls='--',lw=.8,color='#8c979f');ax.set(title=f'{group}: noise x{.5 if group=="DLOW" else 2:g}',xlabel='Actual anchor + warm learning time (s / parent)',ylabel='Primary raw MSE: mean of three seed losses');ax.grid(alpha=.15)
fig.suptitle('Original execution-window risk / cost — frozen candidate, descriptive Pareto front');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=6,frameon=False);fig.tight_layout(rect=[0,.10,1,.93]);finish(fig,'01_RISK_COST_PARETO')
rows=read('ANCHOR_PHASE_RISK.csv');fig,axes=plt.subplots(1,2,figsize=(13,4.6))
for ax,group in zip(axes,['DLOW','DHIGH']):
 for j,m in enumerate(methods):
  rr=[next(r for r in rows if r['group']==group and r['method']==m and r['phase']==p) for p in ['host','transition','pure']];values=[float(r['cell_weighted_mse']) if r['cell_weighted_mse'] else np.nan for r in rr];ax.bar(np.arange(3)+(j-2.5)*.12,values,.12,label=labels[m],color=colors[m])
 ax.set_xticks(np.arange(3),[f'{r["phase"]}\n{int(r["cell_count"]):,} cells' for r in rr]);ax.set(title=group,ylabel='Raw MSE per phase cell');ax.grid(axis='y',alpha=.15)
fig.legend(*axes[1].get_legend_handles_labels(),loc='lower center',ncol=6,frameon=False);fig.suptitle('Risk by FULL-anchor phase — seed-mean losses, fixed ROI, no truth phase labels');fig.tight_layout(rect=[0,.10,1,.93]);finish(fig,'04_ANCHOR_PHASE_RISK')
rows=read('SEED_DISAGREEMENT.csv');fig,axes=plt.subplots(1,2,figsize=(13,4.8))
for ax,group in zip(axes,['DLOW','DHIGH']):
 ens=[];dis=[]
 for m in methods:
  rr=[r for r in rows if r['group']==group and r['method']==m];ens.append(np.mean([float(r['ensemble_MSE_from_identity']) for r in rr]));dis.append(np.mean([float(r['finite_three_seed_disagreement']) for r in rr]))
 x=np.arange(6);ax.bar(x,ens,color='#b8c6d5',label='Ensemble MSE');ax.bar(x,dis,bottom=ens,color='#526d8d',label='Fixed-seed disagreement');ax.set_xticks(x,[labels[m] for m in methods],rotation=30,ha='right');ax.set(title=group,ylabel='MSE: seed-mean = ensemble + disagreement')
fig.suptitle('Three frozen checkpoints: exact finite-ensemble decomposition\nNot population training variance; ensemble is a separate descriptive action.');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False);fig.tight_layout(rect=[0,.10,1,.93]);finish(fig,'06_SEED_DISAGREEMENT')
def identity(p):return dict(path=str(p.resolve()),bytes=p.stat().st_size,blake2b=hashlib.blake2b(p.read_bytes()).hexdigest())
shutil.copy2(ROOT/'track_d/plot_results.py',S/'plot_results.py');d=json.loads((F/'FIGURE_MANIFEST.json').read_text());d['source_script']=identity(S/'plot_results.py');d['presentation_revision']='Only Pareto/phase/disagreement legends and original-window label changed after view_image; plotted CSV values unchanged; initial renders archived under audit/r6_evidence/VISUAL_QA_INITIAL'
for r in d['figures']:
 r['png']=identity(Path(r['png']['path']));r['pdf']=identity(Path(r['pdf']['path']))
(F/'FIGURE_MANIFEST.json').write_text(json.dumps(d,indent=2))
