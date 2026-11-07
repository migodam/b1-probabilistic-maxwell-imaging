"""Plot only saved27 endpoint-segment Phi samples; no physics/GT module import."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).absolute().parents[1];DATA=ROOT/'RESULTS/exposed_endpoint_paths';FIG=ROOT/'FIGURES/exposed_development'
LABELS={'cie_single':'CIE','matrix_green19':'G19'};FAMILIES={'r10_twod_0001':'annulus','r10_twod_0004':'slotted shell','r10_twod_0022':'bridge'}
def main():
 wall=time.perf_counter();rows=json.loads((DATA/'ROWS.json').read_text());summary=json.loads((DATA/'SUMMARY.json').read_text());cfg=json.loads((DATA/'FROZEN_CONFIG.json').read_text());FIG.mkdir(parents=True,exist_ok=True)
 if len(rows)!=27 or len(summary)!=3:raise ValueError('Use the exact saved3x9 table, no extra evaluations')
 fig,axes=plt.subplots(1,3,figsize=(17.5,6.3),layout='constrained')
 for ax,s in zip(axes,summary):
  cid=s['case_id'];r=sorted((x for x in rows if x['case_id']==cid),key=lambda x:x['t'])
  if [x['t'] for x in r]!=cfg['t']:raise ValueError('Frozen sample locations differ')
  t=[x['t'] for x in r];v=[x['original_Phi'] for x in r];first=r[0];peak=max(r,key=lambda x:x['original_Phi'])
  ax.plot(t,v,color='#245ea8',marker='o',ms=5,ls=':',lw=1.1,label='9 paid Phi samples; dotted guide only')
  ax.axhline(first['original_Phi'],color='#a33b20',ls='--',lw=1,label='Initial Phi reference')
  ax.set_xscale('symlog',linthresh=.001,linscale=1.);ax.set_xlim(0,1);ax.set_xticks([0,.001,.01,.1,.5,1],['0','.001','.01','.1','.5','1']);ax.set_yscale('log');ax.grid(alpha=.2)
  ax.set_xlabel('Segment t (symlog; 0–.001 linear)',fontsize=10);ax.set_ylabel('Original physical Phi including R',fontsize=10)
  ax.set_title(cid[-4:]+' '+FAMILIES[cid]+'\n'+LABELS[s['from_method']]+' → '+LABELS[s['to_method']],fontsize=12)
  note=(f"At t=0: dPhi/dt = +{first['directional_derivative']:.6g}\n"
        f"Projected-gradient inf = {first['projected_gradient_inf']:.3g}\n"
        f"Active lower / upper = {first['lower_active_cells']} / {first['upper_active_cells']} cells\n"
        f"Phi(0) = {v[0]:.6g}; Phi(1) = {v[-1]:.6g}")
  ax.text(.025,.97,note,transform=ax.transAxes,va='top',fontsize=9,bbox=dict(boxstyle='round,pad=.35',facecolor='white',edgecolor='#cccccc',alpha=.95))
  ax.annotate(f"Sampled max {peak['original_Phi']:.6g}\nat t={peak['t']:g}",xy=(peak['t'],peak['original_Phi']),xytext=(.56,.76),textcoords='axes fraction',fontsize=9,arrowprops=dict(arrowstyle='->',lw=.8,color='#555555'),bbox=dict(boxstyle='round,pad=.25',fc='white',ec='#cccccc',alpha=.95))
  ax.set_ylim(min(v)*.65,max(v)*4.)
 handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=2,fontsize=10)
 fig.suptitle('One saved legal endpoint segment per exposed case: m(t)=(1−t)m_highPhi+t m_lowPhi\nOnly9 measured samples per case; guides are not a measured continuous curve, all-path barrier or local-minimum certificate',fontsize=13)
 files=[]
 for ext in ('png','pdf'):
  p=FIG/f'EXPOSED_ENDPOINT_PATH_PHI.{ext}';fig.savefig(p,dpi=160);files.append(p.name)
 plt.close(fig)
 receipt=dict(status='COMPLETE_SAVED_ROWS_ONLY',sources=[str((DATA/x).relative_to(ROOT)) for x in ('ROWS.json','SUMMARY.json','FROZEN_CONFIG.json')],objects=3,points_per_object=9,plotted_points=27,files=files,new_forward=0,new_adjoint=0,new_GT_access=0,new_optimizer=0,interpolated_physical_samples=0,visual_guides='Dotted segments connect existing samples only, no smoothing or numerical interpolation',source_matrix_dimension=576,legacy_counter_warning='full_factorization_dimension is an accumulated dimension counter, not matrix order; plots do not create or relabel its values',near_exit_child_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall)
 (FIG/'EXPOSED_ENDPOINT_PATH_PLOT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
