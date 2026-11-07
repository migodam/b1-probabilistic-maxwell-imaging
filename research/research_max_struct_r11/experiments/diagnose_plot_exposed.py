"""Eight-column exposed material plate and saved-stage/corrector trajectories; zero FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,time,json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
R10=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10';OUT=ROOT/'FIGURES/exposed_development';DIAG=ROOT/'delegated/diagnostics'
CASES=('r10_twod_0001','r10_twod_0004','r10_twod_0022');LABELS=('CIE','Dual','Global','R10 G19','freeze A','contracted T','BothBlocks');COLORS=plt.get_cmap('tab10').colors

def sources(cid):
 return [(LABELS[i],R10/'RESULTS/twod'/f'{cid}__{method}__noise_ball_99','prediction') for i,method in enumerate(('cie_single','cie_dual','global_dual','matrix_green19'))]+[(LABELS[i+4],ROOT/'RESULTS/exposed_dev'/f'{cid}__{method}__480CPU','R11_CAP_PREDICTION') for i,method in enumerate(('lagged','contractedLagged','fullGreenBothBlocks'))]
def load(path):
 with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def main():
 wall=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True);meta={x['case_id']:x for x in json.loads((R10/'configs/twod_generator.json').read_text())['cases']};old_scores=json.loads((R10/'delegated/twod_score/scored/SCORES_FINAL.json').read_text());index={(r['case_id'],r['method']):r for r in old_scores};rows=[];truth={};images=[]
 fig,axes=plt.subplots(3,8,figsize=(22,8.6),layout='constrained')
 for ci,cid in enumerate(CASES):
  gt=load(R10/'sealed_truth/twod'/f'{cid}.npz')['m'];truth[cid]=gt;amp=float(gt.max()-.1);mats=[gt]
  info=[('GT',None,None)]
  for mi,(name,d,predname) in enumerate(sources(cid)):
   if mi<4:
    method=('cie_single','cie_dual','global_dual','matrix_green19')[mi]
    if index[cid,method]['task_id']!=d.name:raise RuntimeError('Historical SCORES_FINAL source identity differs')
   arr=load(d/(predname+'.npz'));metric=material_metrics(arr['m'],gt,amp,meta[cid]['high_chi']);mats.append(arr['m']);info.append((name,metric,d));rows.append(dict(case_id=cid,method=name,prediction=str(d/(predname+'.npz')),**metric))
  for col,(m,(name,metric,d)) in enumerate(zip(mats,info)):
   ax=axes[ci,col];im=ax.imshow(m.reshape(24,24),origin='lower',vmin=0,vmax=20,cmap='viridis',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
   title=name if metric is None else f"{name}\nNRMSE {metric['amplitude_NRMSE']:.3f} | IoU {metric['fixed_support_IoU']:.3f}"
   ax.set_title(title,fontsize=10,color='black' if metric is None or metric['material_success'] else '#9f2020',fontweight='bold' if name=='BothBlocks' else 'normal')
   if col==0:ax.set_ylabel(cid[-4:]+' '+meta[cid]['family'].replace('_',' ')+'\nchi high '+str(meta[cid]['high_chi']),fontsize=10)
 fig.colorbar(im,ax=axes,shrink=.75,label='Real material contrast chi; shared scale 0–20')
 fig.suptitle('Exposed DEV endpoints: historical baselines and declared 480-CPU paid incumbents',fontsize=15)
 fig.supxlabel('Amplitude NRMSE = RMSE / (max GT − 0.1). Red labels fail inherited material gate. Shell failure remains. No new inverse/FW; no new16/48 truth read.',fontsize=10)
 fig.savefig(OUT/'MATERIAL_COMPARISON_3x8.png',dpi=160);fig.savefig(OUT/'MATERIAL_COMPARISON_3x8.pdf');plt.close(fig)
 # Keep auxiliary and physical-corrector iteration axes separate: their counts have distinct meaning.
 fig,axes=plt.subplots(3,4,figsize=(20,12),layout='constrained');traj=[]
 for ci,cid in enumerate(CASES):
  gt=truth[cid];amp=float(gt.max()-.1)
  for mi,(name,d,predname) in enumerate(sources(cid)):
   rec=json.loads((d/'RECEIPT.json').read_text());color=COLORS[mi];aux=[];corr=[]
   for stage in range(4):
    n=f'stage_{stage:02d}_physical';p=d/(n+'.npz')
    if not p.exists():continue
    m=load(p)['m'];r=json.loads((d/(n+'.json')).read_text());metric=material_metrics(m,gt,amp,meta[cid]['high_chi']);aux.append(((500,1000,1500,3000)[stage],r['objective'],metric['amplitude_NRMSE']))
    traj.append(dict(case_id=cid,method=name,phase='auxiliary',iteration=(500,1000,1500,3000)[stage],stage=stage,mf=(3,6,10,12)[stage],original_Phi=r['objective'],amplitude_NRMSE=metric['amplitude_NRMSE'],fixed_IoU=metric['fixed_support_IoU']))
   points=[(0,'selected_before_corrector')]+[(int(p.stem.rsplit('_',1)[1]),p.stem) for p in sorted(d.glob('corrector_accepted_*.npz'))]
   endit=(rec.get('corrector') or {}).get('iterations')
   if endit is not None:points.append((endit,predname))
   for it,n in points:
    p=d/(n+'.npz')
    if not p.exists():continue
    m=load(p)['m'];r=json.loads((d/(n+'.json')).read_text());metric=material_metrics(m,gt,amp,meta[cid]['high_chi']);corr.append((it,r['objective'],metric['amplitude_NRMSE']))
    traj.append(dict(case_id=cid,method=name,phase='corrector',iteration=it,state=n,original_Phi=r['objective'],amplitude_NRMSE=metric['amplitude_NRMSE'],fixed_IoU=metric['fixed_support_IoU'],selected_aux_stage=rec.get('auxiliary_selected_stage'),corrector_stop=(rec.get('corrector') or {}).get('stop')))
   aux.sort();corr.sort()
   for j,value in [(0,1),(1,2)]:
    if aux:axes[ci,j].plot([x[0] for x in aux],[x[value] for x in aux],':o',color=color,label=name,ms=4)
   for j,value in [(2,1),(3,2)]:
    if corr:axes[ci,j].plot([x[0] for x in corr],[x[value] for x in corr],':o',color=color,label=name,ms=4)
   selected=rec.get('auxiliary_selected_stage')
   if selected is not None and selected>=0:
    match=next((x for x in aux if x[0]==(500,1000,1500,3000)[selected]),None)
    if match:
     for j,v in ((0,1),(1,2)):axes[ci,j].scatter([match[0]],[match[v]],color=color,marker='*',s=90,zorder=5)
  for j in range(4):
   ax=axes[ci,j];ax.grid(alpha=.2)
   if j%2==0:ax.set_yscale('log');ax.set_ylabel('Original physical Phi')
   else:ax.set_ylabel('Amplitude NRMSE');ax.axhline(.15,color='gray',lw=.7,ls='--')
   if j<2:ax.set_xticks([500,1000,1500,3000],['mf3\n500','mf6\n1000','mf10\n1500','mf12\n3000']);ax.set_xlabel('Auxiliary cumulative iterations; star = Phi-selected start')
   else:ax.set_xscale('symlog',linthresh=100);ax.set_xlim(left=0);ax.set_xlabel('Physical corrector accepted iteration / saved stop endpoint')
   if ci==0:ax.set_title(('Auxiliary Phi','Auxiliary material error','Same inherited corrector: Phi','Same inherited corrector: material error')[j],fontsize=12)
   if j==0:ax.text(.01,.02,cid[-4:]+' '+meta[cid]['family'].replace('_',' '),transform=ax.transAxes,fontsize=10)
 handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=7,fontsize=11)
 fig.suptitle('Saved native-stage and corrector trajectories | complete-trajectory intervention, not same-state comparison',fontsize=15)
 fig.savefig(OUT/'STAGE_CORRECTOR_TRAJECTORIES.png',dpi=150);fig.savefig(OUT/'STAGE_CORRECTOR_TRAJECTORIES.pdf');plt.close(fig)
 (OUT/'MATERIAL_METRICS.json').write_text(json.dumps(rows,indent=2)+'\n');(OUT/'TRAJECTORY_SAVED_POINTS.json').write_text(json.dumps(traj,indent=2)+'\n')
 with (OUT/'TRAJECTORY_SAVED_POINTS.csv').open('w') as f:
  keys=list(dict.fromkeys(k for r in traj for k in r));w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(traj)
 rec=dict(status='COMPLETE',cases=3,columns=8,figure_files=['MATERIAL_COMPARISON_3x8.png','MATERIAL_COMPARISON_3x8.pdf','STAGE_CORRECTOR_TRAJECTORIES.png','STAGE_CORRECTOR_TRAJECTORIES.pdf'],new_inverse=0,new_physical_FW=0,new_adjoint=0,new16_GT_reads=0,new48_GT_reads=0,truth_scope='Previously exposed R10 0001/0004/0022 only',historical_score_index=str(R10/'delegated/twod_score/scored/SCORES_FINAL.json'),CPU_seconds_near_exit=time.process_time(),wall_seconds=time.perf_counter()-wall,cost_scope='Absolute plot child process_time including imports; exit tail not externally sampled',trajectory_scope='Native auxiliary stage and saved accepted-corrector iteration axes kept separate; discrete point guide lines; no cross-platform CPU axis')
 (OUT/'PLOT_RECEIPT.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
if __name__=='__main__':main()
