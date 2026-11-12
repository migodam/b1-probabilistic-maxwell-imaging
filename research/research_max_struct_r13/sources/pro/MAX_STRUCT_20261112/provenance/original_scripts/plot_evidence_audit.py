from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'work/evidence_analysis'
BASE=ROOT/'research_repo/research/research_max_struct_r12/supplement_paired'
x=pd.read_csv(OUT/'projection_split_outcomes.csv')
labels=[r.cohort+str(int(r.case_id.split('_')[-1])) for r in x.itertuples()]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140})
fig,axs=plt.subplots(1,3,figsize=(15,4.9),gridspec_kw={'width_ratios':[1.15,1.15,1.15]})
pos=np.arange(len(x))
axs[0].barh(pos,x.fourier_truncation_fraction*100,color='#416c9b',label='Initial Fourier truncation')
axs[0].barh(pos,x.fixed_DPCS_fraction*100,left=x.fourier_truncation_fraction*100,color='#d99d3f',label='Fixed DPCS mismatch')
axs[0].set_yticks(pos,labels);axs[0].invert_yaxis();axs[0].set_xlim(0,100)
axs[0].set_xlabel('Share of initial lost current energy (%)')
axs[0].set_title('(a) Most current loss is in truncated directions',fontsize=11)
axs[0].legend(loc='lower left',bbox_to_anchor=(0,-.28),fontsize=9,frameon=False)
axs[1].scatter(x.fourier_port_change_relative,pos,color='#416c9b',s=28,label='Initial Fourier truncation')
axs[1].scatter(x.fixed_DPCS_port_change_relative,pos,color='#d99d3f',s=28,label='Fixed DPCS mismatch')
for i,r in enumerate(x.itertuples()):
    axs[1].plot([r.fourier_port_change_relative,r.fixed_DPCS_port_change_relative],[i,i],color='#b5b5b5',lw=.8,zorder=0)
axs[1].set_xscale('log');axs[1].set_yticks(pos,labels);axs[1].invert_yaxis()
axs[1].set_xlabel(r'Induced port change / $\|Y\|_F$')
axs[1].set_title('(b) DPCS change dominates received-field shift',fontsize=11)
for group,marker,color in [('A','o','#33788a'),('B','s','#b75c37'),('C','^','#6e55a3')]:
    z=x[x.cohort==group]
    axs[2].scatter(z.relative_current_loss,z.dr_advantage,marker=marker,color=color,s=45,label=group)
for label,r in zip(labels,x.itertuples()):
    dx,dy=(4,4)
    if label=='C4':dx,dy=4,12
    if label=='C27':dx,dy=4,-14
    if label=='B3':dx,dy=-22,-14
    if label=='A1':dx,dy=-18,5
    axs[2].annotate(label,(r.relative_current_loss,r.dr_advantage),xytext=(dx,dy),textcoords='offset points',fontsize=8)
axs[2].axhline(0,color='#999999',lw=.8,ls='--')
axs[2].set_xlabel('Relative initial current loss')
axs[2].set_ylabel('G19 NRMSE − DR NRMSE')
axs[2].set_title('(c) Total loss does not order the recovery gain',fontsize=11)
axs[2].legend(title='Cohort',frameon=False,loc='upper left')
fig.suptitle('Retrospective R12 analysis: energy loss, data change and recovery are different quantities',fontsize=13,y=1.01)
fig.text(.5,-.055,'11 exposed starts; unchanged saved predictions. All runs later reached the full mf=12 Fourier stage. Positive panel (c) values favor DR.',ha='center',fontsize=9)
fig.tight_layout()
fig.savefig(OUT/'projection_mechanism.png',bbox_inches='tight')
fig.savefig(OUT/'projection_mechanism.pdf',bbox_inches='tight')

rows=[]
for r in x.itertuples():
    with np.load(BASE/'inputs'/f'{r.cohort}__{r.case_id}.npz') as z:
        Y=z['Y']; delta=float(z['sigma'])*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
    with np.load(BASE/'offline_truth'/r.cohort/f'{r.case_id}.npz') as z: truth=z['m']; amplitude=float(z['amplitude'])
    for method in ['g19','dr','ap','direct']:
        folder=BASE/'RESULTS'/f'{r.cohort}__{r.case_id}__{method}'
        with np.load(folder/'raw.npz') as z: m0=z['m']; y0=z['Y_prediction']
        with np.load(folder/'primary.npz') as z: m1=z['m']; y1=z['Y_prediction']
        dm=m1-m0
        rows.append(dict(cohort=r.cohort,case_id=r.case_id,method=method,
            material_step_nrmse=np.linalg.norm(dm)/(np.sqrt(len(dm))*amplitude),
            data_step_over_delta=np.linalg.norm(y1-y0)/delta,
            raw_nrmse=np.linalg.norm(m0-truth)/(np.sqrt(len(dm))*amplitude),
            primary_nrmse=np.linalg.norm(m1-truth)/(np.sqrt(len(dm))*amplitude),
            raw_data_over_delta=np.linalg.norm(y0-Y)/delta,primary_data_over_delta=np.linalg.norm(y1-Y)/delta,
            error_step_cosine=float(np.dot(m0-truth,dm)/(np.linalg.norm(m0-truth)*max(np.linalg.norm(dm),1e-300)))))
pd.DataFrame(rows).to_csv(OUT/'corrector_observation_vs_material.csv',index=False)
