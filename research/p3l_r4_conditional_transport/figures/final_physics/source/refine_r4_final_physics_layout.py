#!/usr/bin/env python3
"""Layout-only visual QA pass while final_physics_figures_v1 Book stays open."""
import os,time,resource,json,csv
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
t0=time.perf_counter();cpu0=resource.getrusage(resource.RUSAGE_SELF)
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
D=Path(__file__).resolve().parents[1]/'figures/final_physics'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':'#26374a','text.color':'#26374a','axes.edgecolor':'#b1bac5','pdf.fonttype':42,'savefig.facecolor':'white'})
def save(fig,name):
 fig.savefig(D/f'{name}.png',dpi=220);fig.savefig(D/f'{name}.pdf');plt.close(fig)
with (D/'chain_costs.csv').open() as f:chain=list(csv.DictReader(f))
C=['#526579','#2b867e','#538fb4','#b68b42'];labels=['Full','Richardson','Global GMRES','GMRES fixed accuracy']
fig,ax=plt.subplots(1,2,figsize=(12.5,5.1));fig.subplots_adjust(left=.185,right=.98,bottom=.24,top=.82,wspace=.32);fig.suptitle('Same-cache fixed chain: measured cost on eight new locked parents',y=.96,fontsize=14,fontweight='bold')
yy=np.arange(4);vals=[float(r['whole_process_seconds']) for r in chain]
ax[0].barh(yy,vals,color=C,height=.58);ax[0].scatter([float(r['controlled_seconds']) for r in chain],yy,marker='|',s=160,color='white',linewidth=2,label='Book controlled wall')
for i,v in enumerate(vals):ax[0].text(v+1,i,f'{v:.2f} s',va='center',fontsize=10)
ax[0].set(yticks=yy,yticklabels=labels,xlabel='Whole arm process wall (seconds)',xlim=(0,85));ax[0].invert_yaxis();ax[0].set_title('A  Includes setup, I/O and postprocessing',loc='left',fontsize=10);ax[0].legend(loc='lower right',fontsize=8,frameon=False)
left=np.zeros(4)
for key,label,col in [('full_forward_RHS','Full forward','#526579'),('R0_setup_RHS','R0 setup','#2b867e'),('reduced_LS_RHS','Reduced complex least squares','#b68b42')]:
 val=np.array([int(r[key]) for r in chain]);ax[1].barh(yy,val,left=left,height=.58,label=label,color=col);left+=val
for i,v in enumerate(left):ax[1].text(v+1000,i,f'{int(v):,}',va='center')
ax[1].set(yticks=yy,yticklabels=[],xlabel='All charged complex right-hand sides',xlim=(0,82000));ax[1].invert_yaxis();ax[1].xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'));ax[1].set_title('B  Reduced-system RHS are included',loc='left',fontsize=10);ax[1].legend(loc='lower right',fontsize=8,frameon=False)
fig.text(.04,.105,'8 parents × 1,024 pCN transitions; beta = .002. Same proposals/uniforms: zero saved path, decision or proof misses.',fontsize=9)
fig.text(.04,.055,'Cached anchor/receiver rights shared by all arms; R0 setup is paid in each approximate arm. One quiet serial run; no ESS or posterior claim.',fontsize=8.5);save(fig,'fig01_locked_chain_cost')
with (D/'risk_fit_means.csv').open() as f:means=list(csv.DictReader(f))
anchor=means[-1];fig,ax=plt.subplots(1,2,figsize=(11.6,5.2));fig.subplots_adjust(left=.085,right=.97,bottom=.23,top=.82,wspace=.28);fig.suptitle('Material risk and nonlinear data fit are different quantities',y=.96,fontsize=15,fontweight='bold')
for j,clip in enumerate(('raw','clipped')):
 a=ax[j];a.axvline(282.873792,color='#a97842',linestyle='--',linewidth=1,label='Original noise-ball threshold')
 for fam,col in [('spline','#2b867e'),('r2','#538fb4')]:
  key=fam+'_'+clip;rr=next(r for r in means if r['key']==key and r['kind']=='raw_risk');re=next(r for r in means if r['key']==key and r['kind']=='representative');rx,ry=float(rr['mean_Phi']),float(rr['mean_full_mse']);ex,ey=float(re['mean_Phi']),float(re['mean_full_mse'])
  a.annotate('',xy=(ex,ey),xytext=(rx,ry),arrowprops=dict(arrowstyle='->',color=col,lw=1.6));a.scatter(rx,ry,s=76,marker='s',color=col,zorder=3);a.scatter(ex,ey,s=76,marker='o',facecolor='white',edgecolor=col,lw=2,zorder=3)
  a.annotate(f"{fam.title()} risk ({rr['noiseball_count']}/64 in ball)",(rx,ry),xytext=(-7,-17 if fam=='spline' else 9),textcoords='offset points',fontsize=8,color=col,ha='right')
  a.annotate(f'{fam.title()} representative',(ex,ey),xytext=(7,-16 if fam=='spline' else 9),textcoords='offset points',fontsize=8,color=col)
 a.scatter(float(anchor['mean_Phi']),float(anchor['mean_full_mse']),s=70,color='#526579',marker='D',zorder=4);a.annotate('Original anchor',(float(anchor['mean_Phi']),float(anchor['mean_full_mse'])),xytext=(7,4),textcoords='offset points',fontsize=8)
 a.set_title(('A  Raw estimators','B  Component-clipped estimators')[j],loc='left',fontsize=11);a.set_xlabel('Mean nonlinear Full Phi');a.set_ylabel('Mean full 576-cell material MSE');a.grid(alpha=.13);a.set_ylim(.0056,.0114);a.margins(x=.30)
ax[0].legend(loc='lower left',fontsize=8,frameon=False)
fig.text(.085,.115,'Identical first 64 L1 parents; fixed Spline seed 2026102501 and frozen R2. Squares: material risk actions; circles: original bounded representatives.',fontsize=8.2)
fig.text(.085,.075,'Arrows link fixed outputs, not an optimization trajectory. Raw/clipped Full evaluation does not prove decoder realizability; representatives are not posterior means.',fontsize=8.2)
fig.text(.085,.035,'All 320 original SLSQP cases stopped at budget. Callback-only selection limitation is preserved; independent Full-J controls are shown separately.',fontsize=8.2);save(fig,'fig03_risk_fit_tradeoff')
cpu=resource.getrusage(resource.RUSAGE_SELF);profile=dict(status='LAYOUT_ONLY_COMPLETE',parent_Book='IV/jobs/final_physics_figures_v1',wall_seconds=time.perf_counter()-t0,CPU_seconds=cpu.ru_utime+cpu.ru_stime-cpu0.ru_utime-cpu0.ru_stime,RSS_peak_bytes=cpu.ru_maxrss,physics_RHS=0,source_data_unchanged=True,scope='Same open Book controlled wall covers layout process and visual QA; supplemental process CPU disclosed separately, no second wall debit')
(D/'LAYOUT_QA_PROFILE.json').write_text(json.dumps(profile,indent=2)+'\n');print(json.dumps(profile))
