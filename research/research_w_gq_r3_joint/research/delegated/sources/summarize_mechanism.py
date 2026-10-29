"""Offline figures and explicit differences from frozen paid public 2x2 audit."""
from pathlib import Path
import json,csv,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];OWN=Path(__file__).resolve().parent;RUN=OWN/'CHART_MECHANISM_RUN';OUT=RUN/'FIGURES'

def main():
 cpu=time.process_time();wall=time.perf_counter();OUT.mkdir(exist_ok=True);summary=json.loads((RUN/'SUMMARY.json').read_text());keys=['uniform_absent','changed_absent','uniform_strong','changed_strong'];loaded={k:np.load(RUN/(k+'.npz')) for k in keys};rows=[];arrays={}
 for k in keys:
  r=summary['states'][k];d=loaded[k];row=dict(state=k,whole_J_relative=r['jacobian_errors']['whole']['scaled_white_relative'],wall_J_relative=r['jacobian_errors']['wall']['scaled_white_relative'],target_J_relative=r['jacobian_errors']['target']['scaled_white_relative'],whole_J_absolute=r['jacobian_errors']['whole']['scaled_white_absolute'],wall_J_absolute=r['jacobian_errors']['wall']['scaled_white_absolute'],target_J_absolute=r['jacobian_errors']['target']['scaled_white_absolute'],dataSchur_relative=r['data_schur']['relative_frobenius_error'],stepSchur_relative=r['step_schur']['relative_frobenius_error'],damped_target_step_relative=r['step']['relative_target_error'],damped_target_step_absolute=r['step']['target_error_norm']);rows.append(row)
  arrays[k+'_J_full_scaled_white']=d['jacobian_scaled_white_full'];arrays[k+'_J_fixed_scaled_white']=d['jacobian_scaled_white_fixed'];arrays[k+'_defect_scaled_white']=d['jacobian_scaled_white_fixed']-d['jacobian_scaled_white_full'];arrays[k+'_data_full']=d['data_full'];arrays[k+'_dataSchur_full']=d['data_schur_full'];arrays[k+'_dataSchur_fixed']=d['data_schur_fixed']
 with (OUT/'FOUR_STATE_TABLE.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 pairs={'wall_at_absent':('changed_absent','uniform_absent'),'wall_at_strong':('changed_strong','uniform_strong'),'target_at_uniform':('uniform_strong','uniform_absent'),'target_at_changed':('changed_strong','changed_absent')}
 norms={}
 for objectname in ['J_full_scaled_white','J_fixed_scaled_white','defect_scaled_white','data_full','dataSchur_full','dataSchur_fixed']:
  for effect,(a,b) in pairs.items():
   v=arrays[a+'_'+objectname]-arrays[b+'_'+objectname];arrays[effect+'_'+objectname]=v;norms[effect+'_'+objectname]=float(np.linalg.norm(v))
  interaction=arrays['changed_strong_'+objectname]-arrays['changed_absent_'+objectname]-arrays['uniform_strong_'+objectname]+arrays['uniform_absent_'+objectname];arrays['interaction_'+objectname]=interaction;norms['interaction_'+objectname]=float(np.linalg.norm(interaction))
 np.savez_compressed(OUT/'RAW_MATCHED_DIFFERENCES.npz',**arrays)
 absolute={k:summary['states'][k]['jacobian_errors']['whole']['scaled_white_absolute'] for k in keys};abs_effect={e:absolute[a]-absolute[b] for e,(a,b) in pairs.items()};abs_effect['interaction']=absolute['changed_strong']-absolute['changed_absent']-absolute['uniform_strong']+absolute['uniform_absent']
 effects=summary['relative_error_factorial_contrasts'];(OUT/'MATCHED_EFFECTS.json').write_text(json.dumps(dict(relative_error_contrasts=effects,absolute_whole_error_contrasts=abs_effect,matrix_difference_norms=norms,caveat='Difference of error norms is not norm of matrix differences. Relative contrasts change denominators; matrix/raw absolute differences retained. All contrasts are for these declared public states only.'),indent=2)+'\n')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
 fig,axes=plt.subplots(2,2,figsize=(11.2,7.8),gridspec_kw={'hspace':.50,'wspace':.52})
 for ax,metric,title in zip(axes.ravel()[:3],['whole_J_relative','target_J_relative','dataSchur_relative'],['(a) Whole tangent relative error','(b) Target tangent relative error','(c) Target Schur relative error']):
  lookup={r['state']:r[metric] for r in rows};a=np.array([[lookup['uniform_absent'],lookup['uniform_strong']],[lookup['changed_absent'],lookup['changed_strong']]])
  im=ax.imshow(a,cmap='Blues',vmin=0,vmax=np.max(a)*1.10,aspect='auto');ax.set_xticks([0,1],['Absent target','Strong target']);ax.set_yticks([0,1],['Uniform wall','Changed wall']);ax.set_title(title,loc='left',fontsize=10)
  for i in range(2):
   for j in range(2):ax.text(j,i,'reference\n≈ 0' if (i,j)==(0,0) else f'{a[i,j]:.4f}',ha='center',va='center',fontsize=12,color='white' if a[i,j]>.55*np.max(a) else 'black')
  fig.colorbar(im,ax=ax,fraction=.045,pad=.025)
 ax=axes[1,1];vals=[effects[k] for k in ['wall_effect_absent','wall_effect_strong','target_effect_uniform','target_effect_changed','interaction']];labels=['Wall\nX absent','Wall\nX strong','Target\nW uniform','Target\nW changed','Interaction'];ax.bar(np.arange(5),vals,color=['#56B4E9','#56B4E9','#E69F00','#E69F00','#999999'],width=.65);ax.axhline(0,color='.3',lw=.7);ax.set_xticks(np.arange(5),labels,fontsize=8);ax.set_ylim(-.8,1.75);ax.set_ylabel('Contrast of whole tangent relative errors');ax.set_title('(d) Matched effects are not additive',loc='left',fontsize=10);ax.grid(axis='y',alpha=.2)
 for i,v in enumerate(vals):ax.text(i,v+(.04 if v>=0 else -.06),f'{v:+.3f}',ha='center',va='bottom' if v>=0 else 'top',fontsize=9)
 fig.suptitle('Public 2×2 control: strong-target change alone can break the fixed chart',fontsize=13,y=1.01)
 fig.text(.08,.02,'One frozen rank-16 chart from uniform wall + absent target; same geometry and three frequencies throughout.\nReal-white, public-scaled tangent errors; Schur uses fixed wall prior and zero target prior. No algorithm damping in panels (a–c).\nFinite error contrasts describe this public example; they are not additive independent scattering amplitudes.',fontsize=8,color='.3')
 fig.subplots_adjust(bottom=.17,top=.94);fig.savefig(OUT/'public_2x2_mechanism.png',dpi=240,bbox_inches='tight');fig.savefig(OUT/'public_2x2_mechanism.pdf',bbox_inches='tight');plt.close(fig)
 receipt=dict(physical_cost=dict(full_fe=0,rhs=0),process_cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall);(OUT/'OFFLINE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
