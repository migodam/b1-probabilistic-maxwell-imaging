"""Post-freeze/post-score supplementary figures and paired attribution data."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
 for split in ['holdout','budget_filled','data_adjoint']:
  if not (ROOT/f'audit/EXPERIMENT_PREDICTIONS_FREEZE_{split.upper()}.json').exists():raise RuntimeError('All63 predictions must freeze before comparison')
 score={s:json.loads((ROOT/f'results/{s}_scores.json').read_text()) for s in ['holdout','budget_filled','data_adjoint']}
 lookup={(s,r['observation_id'],r['method']):r for s,data in score.items() for r in data['rows']};ids=sorted(set(r['observation_id'] for r in score['holdout']['rows']),key=lambda ident:lookup[('holdout',ident,'fullj_native')]['digit'])
 methods=[('budget_filled','fullj_native','Native + remaining time'),('budget_filled','fullj_reciprocal','Reciprocal + remaining time'),('budget_filled','globaldual','Gram dual + Full-J'),('data_adjoint','data_adjoint','Data adjoint + Full-J')]
 material={};runs={}
 for split,method,label in methods:
  for ident in ids:
   r=json.loads((ROOT/f'experiments/frozen_predictions/{split}/{method}/{ident}.json').read_text());runs[(split,ident,method)]=r;material[(split,ident,method)]=np.array(r['m'],float)
 fig,axes=plt.subplots(9,5,figsize=(11.4,18.5),squeeze=False,layout='constrained');im=None
 for i,ident in enumerate(ids):
  digit=lookup[('holdout',ident,'fullj_native')]['digit'];truth=np.load(ROOT/f'truth_after_freeze/holdout/{ident}.npz')['m']
  for col,(m,label) in enumerate([(truth,'Ground truth')]+[(material[(s,ident,a)],label) for s,a,label in methods]):
   ax=axes[i,col];im=ax.imshow(m.reshape(24,24),vmin=.5,vmax=20,cmap='viridis',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
   if i==0:ax.set_title(label,fontsize=9)
   if col==0:ax.set_ylabel(f'Digit {digit}',fontsize=10)
   else:
    s,a,_=methods[col-1];r=lookup[(s,ident,a)];ax.set_xlabel(f'RMSE {r["RMSE"]:.2f} | IoU {r["IoU"]:.2f}',fontsize=8)
 fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.4,label='Real contrast m; common [0.5,20]');fig.suptitle('Supplement: same nine targets, 90-second total caps, ordinary attribution control',fontsize=12)
 fig.savefig(ROOT/'figures/holdout_with_ordinary_control_supplement.png',dpi=170);fig.savefig(ROOT/'figures/holdout_with_ordinary_control_supplement.pdf');plt.close(fig)
 # Replace the single-method scorer's generic four-column rendering with a
 # dedicated two-column control rendering; metrics and solver files unchanged.
 fig,axes=plt.subplots(9,2,figsize=(5,18.5),squeeze=False,layout='constrained')
 for i,ident in enumerate(ids):
  t=np.load(ROOT/f'truth_after_freeze/holdout/{ident}.npz')['m'];m=material[('data_adjoint',ident,'data_adjoint')];r=lookup[('data_adjoint',ident,'data_adjoint')]
  for col,x in enumerate([t,m]):
   axes[i,col].imshow(x.reshape(24,24),vmin=.5,vmax=20,cmap='viridis',interpolation='nearest');axes[i,col].set_xticks([]);axes[i,col].set_yticks([])
   if i==0:axes[i,col].set_title(['Ground truth','Data adjoint + Full-J'][col],fontsize=10)
  axes[i,0].set_ylabel(f'Digit {r["digit"]}');axes[i,1].set_xlabel(f'RMSE {r["RMSE"]:.2f} | IoU {r["IoU"]:.2f}',fontsize=8)
 fig.suptitle('Ordinary optimizer attribution control; same frozen targets',fontsize=10);fig.savefig(ROOT/'figures/data_adjoint_all_digits_main.png',dpi=170);fig.savefig(ROOT/'figures/data_adjoint_all_digits_main.pdf');plt.close(fig)
 paired=[];coverage=[]
 for ident in ids:
  gd=lookup[('budget_filled',ident,'globaldual')];da=lookup[('data_adjoint',ident,'data_adjoint')]
  paired.append(dict(observation_id=ident,digit=gd['digit'],Gram_dual_RMSE=gd['RMSE'],ordinary_data_adjoint_RMSE=da['RMSE'],Gram_minus_ordinary_RMSE=gd['RMSE']-da['RMSE'],Gram_dual_IoU=gd['IoU'],ordinary_data_adjoint_IoU=da['IoU'],Gram_minus_ordinary_IoU=gd['IoU']-da['IoU']))
  ordinary=runs[('data_adjoint',ident,'data_adjoint')];native=json.loads((ROOT/f'experiments/frozen_predictions/holdout/fullj_native/{ident}.json').read_text());st=next(s for s in ordinary['stages'] if s['stage']==24);ref=next(s for s in native['stages'] if s['stage']==24)
  init=np.array(st['initial_material']);refinit=np.array(ref['initial_material']);c=st.get('counts',{})
  coverage.append(dict(observation_id=ident,digit=gd['digit'],control_final576_initial_exact_match_native=bool(np.array_equal(init,refinit)),control_final576_initial_RMS_difference=float(np.sqrt(np.mean((init-refinit)**2))),final576_data_adjoint_gradients=c.get('data_adjoint_gradients',0),final576_data_adjoint_accepted_iterations=c.get('data_adjoint_optimizer_accepted_iterations',0),final576_source_adjoint_RHS=c.get('adjoint_RHS',0),final576_Gram_factorizations=c.get('Gram',0)))
 out=dict(scope='supplementary attribution control frozen after Y generation and before GT; original27 remain primary',paired=paired,coverage=coverage,aggregate_means={s:data['aggregate_means'] for s,data in score.items()},selector='common physical objective only; no GT selection of displayed endpoints')
 (ROOT/'results/ordinary_attribution_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(paired=paired,all_control_initials_match=all(r['control_final576_initial_exact_match_native'] for r in coverage))))
if __name__=='__main__':main()
