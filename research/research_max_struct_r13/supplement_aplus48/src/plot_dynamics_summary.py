"""Render existing authorized CSV evidence only; no model, truth arrays or scores."""
import os,time,json,csv
from pathlib import Path
C=time.process_time();T=time.perf_counter()
os.environ['MPLBACKEND']='Agg'
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import numpy as np
import matplotlib.pyplot as plt
S=Path(__file__).resolve().parents[1]
O=S/'analysis/dynamics_summary';O.mkdir(exist_ok=True)
F=O/'FIGURES';F.mkdir(exist_ok=True)
def rows(path):return list(csv.DictReader(path.open()))
def number(x):
    try:return float(x)
    except (TypeError,ValueError):return float('nan')
colors=dict(cie='#857348',dual_cie='#2569a1',g19='#a64e9d',a_plus='#dc6b20',dr='#26916a')
labels=dict(cie='CIE',dual_cie='Dual-CIE',g19='G19',a_plus='A+ R13-MM',dr='DR')
def save(fig,name):
    for ext in ('png','pdf'):fig.savefig(F/(name+'.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
r=rows(S/'analysis/native_states/READOUT_TRANSITIONS.csv')
g=rows(S/'analysis/native_states/CURRENT_FORCE_TRANSPORT.csv')
fig,ax=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
for m in ('cie','dual_cie','g19','a_plus'):
    rr=[x for x in r if x['method']==m]
    ax[0].scatter([number(x['counterfactual_onehop_loss_change']) for x in rr],[number(x['actual_delta_Phi']) for x in rr],s=15,alpha=.6,c=colors[m],label=labels[m])
    gg=[x for x in g if x['method']==m]
    ax[1].scatter([number(x['gradient_transport_relative_to_prereadout']) for x in gg],[number(x['cosine_postreadout_prereadout_current_gradients']) for x in gg],s=15,alpha=.6,c=colors[m],label=labels[m])
ax[0].set(xscale='symlog',yscale='symlog',xlabel='Fixed-field one-hop loss change (counterfactual for G19)',ylabel='Actual original-Maxwell Phi change')
ax[0].set_xscale('symlog',linthresh=1e-7);ax[0].set_yscale('symlog',linthresh=1e-7)
ax[0].set_xticks([-1,-1e-2,-1e-4,-1e-6,0,1e-4,1e-2])
ax[0].axhline(0,c='.4',lw=.7);ax[0].axvline(0,c='.4',lw=.7)
ax[1].set_xscale('symlog',linthresh=.01)
ax[1].set(xlabel='Material-readout force change / pre-readout force norm',ylabel='Cosine of pre/post readout current gradients',ylim=(-1.05,1.05))
for a in ax:a.grid(alpha=.15);a.legend(fontsize=8,loc='best')
fig.suptitle('Retained native states of 8 post-hoc cases; points are not independent objects',fontsize=11)
save(fig,'readout_physics_and_current_force')
q=rows(S/'analysis/comparison48/ALL48_SCORES.csv')
lookup={(x['case_id'],x['method'],x['view']):x for x in q}
fig,ax=plt.subplots(1,2,figsize=(12,4.7),layout='constrained');records=[]
for m in colors:
    rr=[]
    for cid in sorted({x['case_id'] for x in q}):
        a=lookup[cid,m,'raw'];b=lookup[cid,m,'corrected']
        d=dict(case_id=cid,method=m,delta_Phi=number(b['Phi'])-number(a['Phi']),delta_NRMSE=number(b['NRMSE'])-number(a['NRMSE']),raw_NRMSE=number(a['NRMSE']),corrected_NRMSE=number(b['NRMSE']))
        rr.append(d);records.append(d)
    ax[0].scatter([d['raw_NRMSE'] for d in rr],[d['corrected_NRMSE'] for d in rr],s=22,alpha=.65,c=colors[m],label=labels[m])
    ax[1].scatter([d['delta_Phi'] for d in rr],[d['delta_NRMSE'] for d in rr],s=22,alpha=.65,c=colors[m],label=labels[m])
ax[0].plot([0,1],[0,1],c='.4',lw=.8);ax[0].set(xlabel='Raw material NRMSE',ylabel='Corrected material NRMSE')
ax[1].set_xscale('symlog',linthresh=1e-6);ax[1].set(xlabel='Corrected minus raw: original Phi',ylabel='Corrected minus raw: material NRMSE')
ax[1].set_xticks([-1,-1e-2,-1e-4,-1e-6,0])
ax[1].axhline(0,c='.4',lw=.7);ax[1].axvline(0,c='.4',lw=.7)
for a in ax:a.legend(fontsize=8);a.grid(alpha=.15)
fig.suptitle('Common corrector on all48: physical descent and material improvement are separate',fontsize=11)
save(fig,'corrector_material_and_objective')
(O/'CORRECTOR_CHANGES.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
(O/'RECEIPT.json').write_text(json.dumps(dict(status='EXISTING_AUTHORIZED_TABLE_PLOTS_COMPLETE',CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,new_forward=0,new_scores=0,new_GT_arrays=0,objects=48,native_objects=8,native_scope='posthoc retained snapshots, not independent rows or all iterations'),indent=2)+'\n')
