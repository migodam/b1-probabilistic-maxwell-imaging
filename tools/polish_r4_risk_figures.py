"""Display-only append to closed plots; use saved metrics, preserve original figs."""
from pathlib import Path
import sys,csv,json,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p3r4.evaluation import write_new_json,checkpoint_identity
O=ROOT/'research/p3l_r4_conditional_transport/figures/final_risk'
COLORS=['#6f7881','#353a40','#7a9dc4','#356a9b','#d47a24','#368777','#ad7a9b','#7597a8']
LABELS=['Probit\nanchor','Frozen R2\nCNN','Reoptimized CNN\n512 contexts','CNN\n1024 contexts','CNN + Spline\nfeatures','CNN + Diag\nfeatures','CNN + permuted\nSpline features','CNN + field\nmoment features']
def load(p):
 with p.open() as f:return list(csv.DictReader(f))
def save(fig,name):
 fig.savefig(O/(name+'.png'),dpi=180,bbox_inches='tight',facecolor='white');fig.savefig(O/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
def main():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 bars=load(O/'risk_bar_source_metrics.csv');phase=load(O/'phase_gain_summary.csv')
 fig,axes=plt.subplots(1,2,figsize=(17,5.7),sharey=True)
 ymax=max(float(r['pointwise_mean_CI_high']) for r in bars)*1000+.65
 for ax,group,title in zip(axes,['locked','confirmation'],['L1 · 512 parents','Independent confirmation · 256 parents']):
  rows=[r for r in bars if r['cohort']==group];m=np.array([float(r['mean_full_complex_mse']) for r in rows])*1000;lo=np.array([float(r['pointwise_mean_CI_low']) for r in rows])*1000;hi=np.array([float(r['pointwise_mean_CI_high']) for r in rows])*1000
  ax.bar(range(8),m,color=COLORS,yerr=np.stack((m-lo,hi-m)),capsize=3,error_kw=dict(lw=1));ax.set_xticks(range(8),LABELS,fontsize=8,rotation=25,ha='right');ax.set_title(title);ax.set_ylim(0,ymax);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
  for i,value in enumerate(m):ax.text(i,hi[i]+.12,f'{value:.3f}',ha='center',va='bottom',fontsize=8)
 axes[0].set_ylabel(r'Full 576-cell complex MSE ($\times10^{-3}$)');fig.suptitle('Frozen clipped risk actions · seed losses averaged within each parent',fontsize=14);fig.subplots_adjust(bottom=.29,top=.85,wspace=.08);fig.text(.5,.025,'Pointwise 98.333% parent-mean bootstrap intervals, 4000 draws. Values above whiskers are mean risks. No extra familywise claim.',ha='center',fontsize=10);save(fig,'two_cohort_risk_bars_display_v2')
 fig,axes=plt.subplots(1,2,figsize=(12,5.5),sharey=True)
 for ax,group,title in zip(axes,['locked','confirmation'],['L1 · 512','Confirmation · 256']):
  for x,ref in [(0,'frozen_r2'),(1,'strong_spectral')]:
   row=next(r for r in phase if r['cohort']==group and r['reference']==ref);a=float(row['mean_transition_gain_full576'])*1000;b=float(row['mean_other_gain_full576'])*1000;total=a+b
   ax.bar(x-.15,a,width=.3,color='#cc7728',label='Anchor-transition cells' if x==0 else None);ax.bar(x+.15,b,width=.3,color='#658d9e',label='All other cells' if x==0 else None);ax.plot([x-.34,x+.34],[total,total],color='#333333',lw=1.5)
   ax.text(x,max(a,b,total)+.045,f'total {total:.4f}',ha='center',fontsize=10);ax.text(x-.15,a*.52,f'{a:+.4f}',ha='center',va='center',fontsize=10,color='white',rotation=90);ax.text(x+.17,.035,f'other\n{b:+.6f}',ha='center',va='bottom',fontsize=9,color='#3d6576')
  ax.set_xticks([0,1],['Frozen R2 → Spline','Strong CNN → Spline']);ax.set_title(title);ax.set_ylim(-.04,1.20);ax.axhline(0,color='#444444',lw=.8);ax.grid(axis='y',alpha=.12);ax.set_axisbelow(True)
 axes[0].set_ylabel(r'Mean MSE gain contribution ($\times10^{-3}$; denominator576)');axes[0].legend(frameon=False,fontsize=9,loc='upper right');fig.suptitle('Anchor-phase gain decomposition · observed risk, not a mechanism proof',fontsize=13);fig.subplots_adjust(bottom=.24,top=.84,wspace=.10);fig.text(.5,.08,'Online-anchor transition mask −1 < s < 1, versus its full576-cell complement. Numeric labels expose very small complement contributions.',ha='center',fontsize=9);fig.text(.5,.035,'Per-cell seed loss means, not averaged images. Contributions add to total gain; no truth-mask selection or posterior claim.',ha='center',fontsize=9);save(fig,'anchor_transition_gain_decomposition_display_v2')
 shutil.copyfile(__file__,O/'display_source.py')
 write_new_json(O/'DISPLAY_REVISION_V2.json',dict(status='DISPLAY_ONLY_APPEND',changes=['mean labels above CI whiskers','explicit numeric labels for verysmall complement phase contributions'],data_sources=[checkpoint_identity(O/'risk_bar_source_metrics.csv'),checkpoint_identity(O/'phase_gain_summary.csv')],original_figures_preserved=True,statistics_changed=False,source_predictions_changed=False,model_inference=0,physics_RHS=0))
if __name__=='__main__':main()
