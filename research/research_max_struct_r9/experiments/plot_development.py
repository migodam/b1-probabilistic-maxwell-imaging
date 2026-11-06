"""Plot only sealed/scored four-case development outputs; zero new physics.

Truth is read only through score_development.truth_records. No confirmation
truth, metric recomputation, GT iterate choice or solver import is permitted.
"""
from pathlib import Path
import os
os.environ['MPLBACKEND']='Agg'
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import json,time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import score_development as scoring

ROOT=Path(__file__).resolve().parents[1]
SCORE=ROOT/'RESULTS/development_scoring'
FIG=ROOT/'FIGURES'
CASES=('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3','r8_test_0000','r8_test_0027')
NAMES=dict(zip(CASES,['b20d','40ac','0000','0027']))
METHODS=[
 ('cie_dual','fixed_001','Dual fixed','#0072B2','--'),
 ('global_dual','fixed_001','Global fixed','#009E73','--'),
 ('a_plus','fixed_001','A+ fixed','#D55E00','--'),
 ('a_plus_sourcewise','fixed_001','A+ source fixed','#CC79A7','--'),
 ('cie_dual','noise_ball_99','Dual noise','#0072B2','-'),
 ('global_dual','noise_ball_99','Global noise','#009E73','-'),
 ('a_plus','noise_ball_99','A+ noise','#D55E00','-'),
 ('a_plus_sourcewise','noise_ball_99','A+ source noise','#CC79A7','-'),
 ('profile_homotopy',None,'Profile continuation','#6A3D9A','-'),
 ('profile_fixed',None,'Profile fixed','#6A3D9A',':'),
 ('data_homotopy',None,'Data continuation','#A6761D','-.'),
 ('physical_lbfgs','fixed_001','Direct physical','#111111','-')]


def load(path):return json.loads(Path(path).read_text())


def select(rows,case,method,rule):
    found=[r for r in rows if r['case_id']==case and r['method']==method
           and (r.get('dpcs_rule')==rule or (rule=='fixed_001' and r.get('dpcs_rule') is None))
           and r.get('prediction_path')]
    fresh=[r for r in found if r['origin'].startswith('NEW_R9_')]
    chosen=fresh if fresh else [r for r in found if r['origin']=='REUSED_R8_UNIFORM_ENDPOINT']
    if len(chosen)!=1:
        raise ValueError(f'Expected exactly one protocol-selected endpoint: {case}/{method}/{rule}; got{len(chosen)}')
    return chosen[0]


def main():
    wall,cpu=time.perf_counter(),time.process_time()
    scored=load(SCORE/'SCORES.json');budget=load(SCORE/'BUDGET_SCORES.json')
    seal=load(ROOT/'audit/DEVELOPMENT_PUBLIC_SEAL.json')
    if scored['status']!='SCORED_FROZEN_FOUR_CASE_DEVELOPMENT' or seal['status']!='ALL_REQUESTED_R9_DEVELOPMENT_PUBLIC_OUTPUTS_SEALED':
        raise ValueError('Complete four-case development seal and scoring required')
    if tuple(seal['cases'])!=CASES or set(r['case_id'] for r in scored['rows'])!=set(CASES):
        raise ValueError('Plot scope differs from the four authorized DEV cases')
    if scored['new_tasks']!=43:raise ValueError('Expected43 newly sealed/scored development tasks')
    if tuple(scoring.CASES)!=CASES:raise ValueError('Truth helper case scope changed')
    truths=scoring.truth_records(Path(seal['r8_repo']))
    if set(truths)!=set(CASES):raise ValueError('Unexpected truth scope')
    rows=scored['rows'];selected={(case,method,rule):select(rows,case,method,rule)
        for case in CASES for method,rule,*_ in METHODS}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':10,
        'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,'pdf.fonttype':42})
    FIG.mkdir(exist_ok=True,parents=True)
    matrix=[('GT',None,None),('Dual fixed','cie_dual','fixed_001'),
        ('Global fixed','global_dual','fixed_001'),('A+ fixed','a_plus','fixed_001'),
        ('Dual noise','cie_dual','noise_ball_99'),('A+ noise','a_plus','noise_ball_99'),
        ('Profile continuation','profile_homotopy',None),('Direct physical','physical_lbfgs','fixed_001')]
    fig=plt.figure(figsize=(23,14))
    gs=fig.add_gridspec(4,9,width_ratios=[1]*8+[.065],left=.035,right=.975,bottom=.065,top=.86,
                       hspace=.28,wspace=.065)
    matrix_records=[]
    for i,case in enumerate(CASES):
        amp=truths[case]['amplitude'];arrays=[truths[case]['m']]
        for _,method,rule in matrix[1:]:
            with np.load(selected[case,method,rule]['prediction_path'],allow_pickle=False) as z:
                arrays.append(z['m'].ravel().copy())
        # One amplitude-normalized scale per case, shared by every column.
        # Its upper bound includes all displayed predictions, retaining spikes.
        vmax=max(1.,max(float(np.max(a))/amp for a in arrays))
        for j,((label,method,rule),a) in enumerate(zip(matrix,arrays)):
            ax=fig.add_subplot(gs[i,j]);im=ax.imshow((a/amp).reshape(24,24),origin='lower',
                vmin=0,vmax=vmax,cmap='viridis',interpolation='nearest')
            ax.set_xticks([]);ax.set_yticks([])
            if method is None:
                title='NRMSE 0 (reference)\nPhi not evaluated'
            else:
                r=selected[case,method,rule];old=r['origin']=='REUSED_R8_UNIFORM_ENDPOINT'
                title=f"NRMSE {r['material_nrmse']:.3f}\nPhi {r['original_Phi']:.2e}"+('  [old Mac]' if old else '')
                matrix_records.append(dict(case_id=case,column=label,task_id=r['task_id'],origin=r['origin'],prediction_path=r['prediction_path'],scores_reused=True))
            ax.set_title(title,pad=5)
            if i==0:ax.text(.5,1.31,label,ha='center',va='bottom',transform=ax.transAxes,fontsize=12,fontweight='bold')
            if j==0:ax.set_ylabel(f'{NAMES[case]}\namplitude {amp:g}',fontsize=12,fontweight='bold',labelpad=10)
        cax=fig.add_subplot(gs[i,8]);cb=fig.colorbar(im,cax=cax);cb.set_label('m / amplitude',fontsize=9)
    fig.suptitle('PUBLIC_DEVELOPMENT | Four fixed exposed cases; sealed endpoint selection',fontsize=19,y=.982)
    fig.text(.5,.025,'Fresh local fixed controls take priority. Case0000 fixed/direct reuse old Mac results. Colors share one amplitude-normalized range within each row; spikes are not clipped.\nNRMSE and Phi are existing endpoint scores; GT Phi is not evaluated here. No confirmation truth or GT-selected iterates.',ha='center',fontsize=10)
    fig.savefig(FIG/'development_materials.png',dpi=180);fig.savefig(FIG/'development_materials.pdf');plt.close(fig)

    grouped={}
    for r in budget['rows']:grouped.setdefault(r['task_id'],[]).append(r)
    fig,axs=plt.subplots(4,2,figsize=(16,17))
    fig.subplots_adjust(left=.075,right=.98,bottom=.185,top=.93,hspace=.40,wspace=.20)
    curve_records=[]
    for i,case in enumerate(CASES):
        upper_labels=[];old_labels=[];cpu_max=0.
        for method,rule,label,color,style in METHODS:
            r=selected[case,method,rule];rr=sorted(grouped.get(r['task_id'],[]),key=lambda x:x['index'])
            if not rr:raise ValueError('Missing paid budget scores for '+r['task_id'])
            if any(x['CPU_seconds'] is None for x in rr):raise ValueError('Unmeasured prefix cannot be placed on a CPU curve')
            x=np.array([a['CPU_seconds'] for a in rr]);err=np.array([a['material_nrmse'] for a in rr]);phi=np.array([a['original_Phi'] for a in rr])
            if np.any(np.diff(x)<-1e-9):raise ValueError('Nonmonotone paid cost')
            upper=any('UPPER_SHIFT' in a['CPU_semantics'] for a in rr);old=r['origin']=='REUSED_R8_UNIFORM_ENDPOINT'
            if upper:upper_labels.append(label)
            if old:old_labels.append(label)
            cpu_max=max(cpu_max,float(x.max()))
            width=2.1 if method=='profile_homotopy' else 1.45
            for ax,y in zip(axs[i],[err,phi]):
                ax.step(x,y,where='post',color=color,linestyle=style,linewidth=width,alpha=.90)
                # Mark paid endpoint samples, never a minimum selected by GT.
                if old or upper:ax.plot(x[-1],y[-1],marker='^' if upper else 's',markersize=5,
                    markerfacecolor='none',markeredgecolor=color,linestyle='None')
            curve_records.append(dict(case_id=case,label=label,task_id=r['task_id'],origin=r['origin'],samples=len(rr),
                CPU_semantics=sorted(set(a['CPU_semantics'] for a in rr)),old_reuse=old,startup_prefix_upper_bound=upper))
        for ax in axs[i]:
            ax.set_xlim(0,cpu_max*1.03);ax.grid(alpha=.18);ax.set_xlabel('Recorded paid CPU seconds')
        axs[i,0].set_ylim(bottom=0);axs[i,0].set_ylabel('Amplitude-normalized material NRMSE')
        axs[i,1].set_yscale('log');axs[i,1].set_ylabel('Original Phi (log scale)')
        subtitle='old Mac fixed/direct reuse' if old_labels else 'fresh local controls'
        axs[i,0].set_title(f'{NAMES[case]} | {subtitle}',fontsize=12,loc='left')
        scope='upper-shift: '+', '.join(upper_labels) if upper_labels else 'recorded startup-inclusive CPU'
        axs[i,1].set_title(f'{NAMES[case]} | {scope}',fontsize=10.5,loc='left')
    handles=[Line2D([],[],color=c,linestyle=st,label=label,linewidth=1.8) for _,_,label,c,st in METHODS]
    handles += [Line2D([],[],color='black',marker='s',markerfacecolor='none',linestyle='None',label='Old Mac reuse (0000)'),
                Line2D([],[],color='black',marker='^',markerfacecolor='none',linestyle='None',label='Startup prefix upper bound')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.065),ncol=4,frameon=False,fontsize=10)
    fig.suptitle('PUBLIC_DEVELOPMENT | Paid original-Phi incumbents: material error and objective',fontsize=17,y=.982)
    fig.text(.5,.955,'All12 fixed method/rule controls retained per case, including profile-fixed, data-continuation and direct negatives.',ha='center',fontsize=11)
    fig.text(.5,.022,'Steps show the last saved PAID original-Phi incumbent; material NRMSE is not minimized or made monotone for display.\nLandscape upper-shift CPU includes the unseparated output tail conservatively. Historical Mac reuse is not a fresh isolated timing control.\nAuxiliary dense Cholesky work is included in CPU; physical LU alone is not total compute. No interpolation, GT iterate selection or new physics.',ha='center',fontsize=10)
    fig.savefig(FIG/'development_cost_curves.png',dpi=180);fig.savefig(FIG/'development_cost_curves.pdf');plt.close(fig)
    path=ROOT/'audit/development_plot_receipt.json';past=load(path) if path.exists() else []
    if isinstance(past,dict):past=[past]
    receipt=dict(status='RENDERED_AWAITING_VISUAL_REVIEW',scope='PUBLIC_DEVELOPMENT_FOUR_SPECIFIED_CASES_ONLY',
        cases=list(CASES),truth_paths_opened=[truths[c]['truth_path'] for c in CASES],truth_loader='score_development.truth_records',
        confirmation_truth_access=False,new_GT_scoring_calls=0,new_physical_calls=0,new_LU=0,new_complex_RHS=0,
        endpoint_selection='Protocol-only: fresh local same method/rule preferred; otherwise exact reused R8 uniform endpoint. No GT ranking.',
        endpoint_columns=matrix_records,curves=curve_records,
        files=['FIGURES/development_materials.png','FIGURES/development_materials.pdf','FIGURES/development_cost_curves.png','FIGURES/development_cost_curves.pdf'],
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    past.append(receipt);path.write_text(json.dumps(past,indent=2)+'\n')
    print(json.dumps(dict(status=receipt['status'],files=receipt['files'],CPU_seconds=receipt['CPU_seconds'])))


if __name__=='__main__':main()
