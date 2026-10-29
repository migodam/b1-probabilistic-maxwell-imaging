"""Scientific figure from paid same-state dev Jacobian records only."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
ROOT=Path(__file__).resolve().parents[3]
OWN=Path(__file__).resolve().parent
OUT=OWN/'CONSEQUENCE_FIGURE'
def main():
 OUT.mkdir(exist_ok=True)
 data=json.loads((OWN/'ENDPOINT_CONSEQUENCE_RUN/SUMMARY.json').read_text())
 rows=[r for r in data['rows'] if r.get('chart_pair_status')=='SAME_STATE_PAIR']
 values=[]
 for r in rows:
  category='refresh_trigger_audit' if 'CURRENT_AUDIT_' in r['source'] else ('fixed_dev' if r['method']=='FIXED_CURRENT_CHART' else 'adaptive_dev')
  values.append(dict(category=category,method=r['method'],source=r['source'],selected_endpoint=r['source'].endswith('/endpoint_exact_jacobian.npz'),whole_tangent_relative=r['tangent_errors']['whole']['relative_frobenius'],wall_tangent_relative=r['tangent_errors']['wall']['relative_frobenius'],target_tangent_relative=r['tangent_errors']['target']['relative_frobenius'],targetSchur_relative=r['dataSchur']['relative_frobenius_error'],Schur_leading_relative_remainder=r['first_order_target_effective']['relative_linearization_remainder'],target_effective_tangent_relative=r['first_order_target_effective']['relative_to_profiled_response'],exact_projected_stationarity=r['exact_projected_stationarity'],chart_projected_stationarity=r['chart_projected_stationarity_same_full_diagonal'],full_dataSchur_min_eigenvalue=r['dataSchur']['full_eigenvalues'][0]))
 with (OUT/'same_state_consequences.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':10,'pdf.fonttype':42,'ps.fonttype':42})
 fig=plt.figure(figsize=(10,6.6));gs=fig.add_gridspec(2,2,height_ratios=[1,.26],hspace=.30,wspace=.28)
 axes=[fig.add_subplot(gs[0,i]) for i in range(2)]
 groups=[('fixed_dev','Fixed chart dev records','#0072B2','s'),('adaptive_dev','Adaptive chart dev records','#D55E00','o'),('refresh_trigger_audit','Adaptive refresh-trigger audits','#009E73','^')]
 for cat,label,color,marker in groups:
  rs=[r for r in values if r['category']==cat]
  x=[r['whole_tangent_relative'] for r in rs]
  for ax,key in zip(axes,['targetSchur_relative','Schur_leading_relative_remainder']):
   ax.scatter(x,[r[key] for r in rs],s=44,c=color,marker=marker,alpha=.82,edgecolors='white',linewidths=.5,label=f'{label} (n={len(rs)})',zorder=3)
 for ax in axes:
  ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(.012,.85);ax.xaxis.set_major_locator(FixedLocator([.02,.05,.1,.2,.5]));ax.xaxis.set_major_formatter(FixedFormatter(['0.02','0.05','0.1','0.2','0.5']));ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(True,which='major',alpha=.25);ax.grid(True,which='minor',alpha=.08);ax.set_xlabel(r'Scaled tangent error $\|\widehat J-J\|_F/\|J\|_F$')
 axes[0].set_ylabel(r'Target Schur error $\|\widehat S-S\|_F/\|S\|_F$');axes[0].set_title('(a) Actual finite target-information change',loc='left');axes[0].set_ylim(.02,1.5)
 axes[1].set_ylabel(r'Leading-term remainder $\|\Delta S-dS\|_F/\|\Delta S\|_F$');axes[1].set_title('(b) First-order prediction remainder',loc='left');axes[1].set_ylim(.02,2.5);axes[1].axhline(1,color='.45',linestyle='--',linewidth=.8,zorder=1);axes[1].text(.015,1.08,'remainder = actual change',fontsize=8,color='.4')
 endpoints=[r for r in values if r['selected_endpoint']]
 for r in endpoints:
  for i,(ax,key) in enumerate(zip(axes,['targetSchur_relative','Schur_leading_relative_remainder'])):
   xy=(r['whole_tangent_relative'],r[key]);ax.scatter(*xy,s=110,facecolors='none',edgecolors='black',linewidths=1.0,zorder=5)
   if i==0:
    name='Fixed selected' if r['category']=='fixed_dev' else 'Adaptive selected';offset=(-82,-28) if r['category']=='fixed_dev' else (10,15)
    ax.annotate(name,xy,xytext=offset,textcoords='offset points',fontsize=8,arrowprops=dict(arrowstyle='-',lw=.6,color='.3'))
 tableax=fig.add_subplot(gs[1,:]);tableax.axis('off')
 table=tableax.table(cellText=[[('Fixed' if r['category']=='fixed_dev' else 'Adaptive'),f"{r['exact_projected_stationarity']:.6f}",f"{r['chart_projected_stationarity']:.6f}",('PASS' if r['exact_projected_stationarity']<.005 else 'FAIL')] for r in endpoints],colLabels=['Selected dev endpoint','Exact projected stationarity','Chart value, same normalization','Exact < 0.005'],cellLoc='center',loc='center',colWidths=[.22,.26,.34,.18]);table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.35)
 for (rr,cc),cell in table.get_celld().items():
  cell.set_linewidth(.4);cell.set_edgecolor('.75');cell.set_facecolor('#eaf0f4' if rr==0 else 'white')
 handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.98),ncol=3,frameon=False,fontsize=8)
 fig.suptitle('A local chart check is not a stationarity or target-accuracy certificate',fontsize=12,y=1.025)
 fig.text(.07,.035,'33 saved same-state pairs from one development case; endpoint/final records may overlap.\nReal–imag stacked J, noise scale 0.01 and public parameter scales; no new PDE or truth reads.',fontsize=8,color='.3')
 fig.subplots_adjust(left=.085,right=.98,top=.89,bottom=.13)
 fig.savefig(OUT/'same_state_consequences.png',dpi=240,bbox_inches='tight');fig.savefig(OUT/'same_state_consequences.pdf',bbox_inches='tight');plt.close(fig)
 (OUT/'CAPTION.md').write_text('Same-state development chart consequences. Left: whole tangent relative error versus actual target Schur relative error. Right: remainder of the first-order Schur prediction, normalized by actual finite Schur change. Refresh-trigger audit pairs compare the pre-refresh chart to the already-paid full Jacobian at the same state. Selected endpoint circles and stationarity table are descriptive; target accuracy is not evaluated. 33 saved records are not independent scenes and include overlapping endpoint/final entries. All matrices use sigma=.01, public parameter scaling, fixed wall prior, zero target prior, and no numerical damping in the data Schur. No PDE or truth reads.\n')
 print(json.dumps(dict(records=len(values),outputs=[str(p) for p in OUT.iterdir()])))
if __name__=='__main__':main()
