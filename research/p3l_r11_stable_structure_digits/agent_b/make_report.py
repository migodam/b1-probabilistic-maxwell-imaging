"""Zero-solve reporting for Sol-B, using existing kernel and inference records."""
import os
os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import csv,json
import numpy as np
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent
s=json.loads((D/'KERNEL_SUMMARY.json').read_text());u=json.loads((D/'R128_SUMMARY.json').read_text());rows=list(csv.DictReader((D/'KERNEL_NODE_RESULTS.csv').open()));r128=list(csv.DictReader((D/'R128_NODE_RESULTS.csv').open()))
methods=list(s['methods']);labels=['Anchor 32','Anchor 64','Anchor 96','Green 32','Green 64','Green 96','Weighted old 64','Weighted rebuild 64']
colors=['#789bd0','#4063ad','#163d7d','#7dcac5','#389e99','#176b67','#b57042','#884e9e']
fig,ax=plt.subplots(1,3,figsize=(16.2,5.2),layout='constrained')
y=np.arange(len(methods)+1)
qual=[s['methods'][m]['qualified'] for m in methods]+[u['qualified']]
actual=[s['methods'][m]['actual_error_below_point01'] for m in methods]+[u['actual_error_below_point01']]
lab=labels+['Anchor 128*'];col=colors+['#4063ad']
ax[0].barh(y,actual,color='#d7dce3',label='Actual error <= 0.01 (evaluation only)')
ax[0].barh(y,qual,color=col,label='Residual + potential gate accepted')
for j,(q,a) in enumerate(zip(qual,actual)):ax[0].text(max(q,a)+2,j,f'{q} / {a}',va='center',fontsize=8)
ax[0].set(yticks=y,yticklabels=lab,xlim=(0,230),xlabel='Nodes / 192: dark = gate; pale = actual error <= 0.01',title='Few fixed-space proposals qualify');ax[0].invert_yaxis();ax[0].xaxis.label.set_size(9)
for m,l,c in zip(methods,labels,colors):
 er=np.sort([float(r['actual_approx_white_error']) for r in rows if r['method']==m and r['basis_available']=='True'])
 ax[1].plot(er,np.arange(1,len(er)+1)/len(er),label=l,color=c,ls='-' if 'ANCHOR' in m or 'REBUILD' in m else '--',lw=1.4)
e=np.sort([float(r['actual_approx_white_error'])for r in r128]);ax[1].plot(e,np.arange(1,193)/192,label='Anchor 128*',color='#172f70',lw=2)
ax[1].axvline(.01,color='black',ls=':',lw=1);ax[1].set(xscale='log',xlim=(1e-14,30),ylim=(0,1.03),xlabel='Actual whitened response error',ylabel='Fraction of computed nodes',title='Total-field coordinates help, but are insufficient');ax[1].legend(fontsize=7,loc='upper left');ax[1].grid(alpha=.2)
basis=np.array([s['methods'][m]['fixed_basis_setup_seconds'] for m in methods]);shape=np.array([s['methods'][m]['current_shape_setup_seconds']for m in methods]);pure=np.array([s['methods'][m]['pure_online_seconds']for m in methods]);fallback=np.array([s['methods'][m]['online_including_fallback_seconds']-s['methods'][m]['pure_online_seconds']for m in methods]);bottom=np.zeros(8)
for a,c,l in [(basis,'#c8ceda','Fixed basis'),(shape,'#8fa9c4','Current shape + bound setup'),(pure,'#4d789e','Reduced solve + gate'),(fallback,'#cf816e','Independent exact fallback')]:
 ax[2].barh(np.arange(8),a,left=bottom,color=c,label=l);bottom+=a
ax[2].axvline(s['exact_A']['cold_task_seconds'],color='black',ls='--',lw=1.2,label='Exact A reference')
ax[2].set(yticks=np.arange(8),yticklabels=labels,xlabel='Measured cold task seconds',title='Core matched window: full paid cost');ax[2].invert_yaxis();ax[2].legend(fontsize=7,loc='lower right');ax[2].grid(axis='x',alpha=.2)
fig.suptitle('R11 Sol-B: 24 exposed R10 path shapes x 8 stored amplitude nodes',fontsize=14)
fig.text(.5,-.035,'* Optional rank 128 is a separate supplemental timing window. Eight rebuild nodes have rank 61 and were not approximated. No GT; 0 new 576-state solves.',ha='center',fontsize=9)
fig.savefig(D/'FIXED_TOTAL_FIELD_QUALIFICATION.png',dpi=190,bbox_inches='tight');fig.savefig(D/'FIXED_TOTAL_FIELD_QUALIFICATION.pdf',bbox_inches='tight')
plt.close(fig)
# Standalone short evidence table for parent assembly.
with (D/'METHOD_COMPARISON.csv').open('w',newline='')as f:
 keys=['method','nodes','qualified','actual_error_below_point01','fallback_nodes','max_actual_approx_white_error','median_actual_approx_white_error','cold_task_seconds','timing_window']
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
 for m,v in s['methods'].items():w.writerow(dict(method=m,nodes=192,qualified=v['qualified'],actual_error_below_point01=v['actual_error_below_point01'],fallback_nodes=v['fallback_nodes'],max_actual_approx_white_error=v['max_actual_approx_white_error'],median_actual_approx_white_error=v['median_actual_approx_white_error'],cold_task_seconds=v['cold_task_including_fallback_seconds'],timing_window='matched_core'))
 w.writerow(dict(method='TOTAL_ANCHOR_R128',nodes=192,qualified=u['qualified'],actual_error_below_point01=u['actual_error_below_point01'],fallback_nodes=u['fallback_nodes'],max_actual_approx_white_error=u['max_actual_approx_white_error'],median_actual_approx_white_error=u['median_actual_approx_white_error'],cold_task_seconds=u['comparable_component_timing']['cold_task_seconds'],timing_window='separate_supplement'))
print(D/'FIXED_TOTAL_FIELD_QUALIFICATION.png')
