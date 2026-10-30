"""Plot existing score/diagnostic JSON only. No physics or private-data IO."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
def read(name):return json.loads((P/name).read_text())
def save(fig,name):
    fig.tight_layout();fig.savefig(P/(name+'.png'),dpi=155);fig.savefig(P/(name+'.pdf'));plt.close(fig)

control=read('SAMEJ_VS_REFRESHED_CONTROL.json');names=[r['parent'][-3:] for r in control]
fig,axs=plt.subplots(1,2,figsize=(10.5,4.6))
axs[0].bar(names,[r['sameJ_nested_span_violation'] for r in control],color='#667a8b')
axs[0].set_ylabel('Distance outside original weak8 span');axs[0].set_title('Same J, nested SVD: 8 → 4 → 2')
axs[0].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
axs[1].bar(names,[r['secant']['jacobian_relative_change'] for r in control],color='#b5633c')
axs[1].set_ylabel('||J(new state) − J(anchor)||F / ||J(anchor)||F');axs[1].set_title('Actually paid candidate J changes')
for ax in axs:ax.set_xlabel('Exposed R9 parent ID suffix');ax.grid(axis='y',alpha=.2)
fig.suptitle('Algebraic reuse versus a new local solve; same observed data y throughout',fontsize=12)
save(fig,'SAMEJ_VS_REFRESHED_CONTROL')

layers=read('p1_diagnostics/LAYER_COSTS.json');methods=['H0_RML','R9_STYLE_FIXED8','H1_FULL16','H2_FULL16','HA_FULL16','H1_RECEIVER16','H1_PAIR16','H1_FULL32']
fig,ax=plt.subplots(figsize=(10,5.4))
for method in methods:
    rr=[r for r in layers if r['method']==method]
    xx=sorted(set(r['layer'] for r in rr))
    ax.plot(xx,[np.mean([r['RHS_layer'] for r in rr if r['layer']==i]) for i in xx],'.-',label=method)
ax.set_xticks([1,2,3],['1 (.02 RMS)','2 (.05 RMS)','3 (.10 RMS)'])
ax.set_ylabel('Mean paid RHS per branch per layer');ax.set_xlabel('Predictor layer')
ax.set_title('Actual trace costs; shared COMMON and final checks reported separately')
ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2)
save(fig,'P1_LAYER_COSTS')

candidates=read('p1_score/CANDIDATE_SCORES.json');fig,axs=plt.subplots(1,2,figsize=(13,5.4))
for ax,key,title in zip(axs,['parameter_span_outside_initial_weak8','material_tangent_span_outside_initial_weak8'],
                         ['Parameter-coordinate span','Original576 material-tangent span']):
    data=[[r[key]['relative'] for r in candidates if r['method']==m and r['qualified'] and r['anchor_material_RMS']>.002 and key in r] for m in methods]
    ax.boxplot(data,tick_labels=methods,showfliers=False);ax.set_xticklabels(methods,rotation=45,ha='right',fontsize=8)
    ax.set_ylabel('Relative component outside original R9 weak8');ax.set_ylim(0,1.04);ax.grid(axis='y',alpha=.2);ax.set_title(title)
fig.suptitle('Qualified nonanchor candidates leave weak8; ε=.075 coverage remains 0/4 for every arm',fontsize=12)
save(fig,'P1_OUTSIDE_WEAK8')
