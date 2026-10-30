"""Final figures from frozen endpoint arrays and actual paid material trajectories."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
METHODS=['trf','data_lbfgs','globaldual','uniform_bridge','uniform_best_only','uniform_multistart']
COLORS=['#0072b2','#009e73','#d55e00','#882255','#44aabb','#332288']

def main():
    scores=json.loads((ROOT/'results/fresh_v3/SCORES.json').read_text())['rows']
    cases=sorted(set(r['case_id'] for r in scores),key=lambda c:next(r['digit'] for r in scores if r['case_id']==c))
    totals=json.loads((ROOT/'COST_TOTALS.json').read_text())['fresh_method_totals']
    figures=ROOT/'figures'
    for stem in ('reconstruction_and_cost','maxwell_branch_trajectories'):
        for ext in ('png','pdf','svg'):
            old=figures/(stem+'.'+ext)
            dest=figures/(stem+'_v1.'+ext)
            if old.exists() and not dest.exists():shutil.copy2(old,dest)
    fig=plt.figure(figsize=(19,12),layout='constrained')
    grid=fig.add_gridspec(4,7,height_ratios=[1,1,1,1.0])
    for i,case in enumerate(cases):
        truth=np.load(ROOT/'truth_after_freeze/fresh'/(case+'.npz'))['m']
        digit=next(r['digit'] for r in scores if r['case_id']==case)
        for j,key in enumerate(['truth']+METHODS):
            ax=fig.add_subplot(grid[i,j])
            if key=='truth':material=truth
            else:
                endpoint=json.loads((ROOT/'published_results/fresh_v3'/(case+'_'+key+'.json')).read_text())
                material=np.asarray(endpoint['top1']['m'])
            image=ax.imshow(material.reshape(24,24),vmin=.5,vmax=20,cmap='viridis',interpolation='nearest')
            ax.set_xticks([]);ax.set_yticks([])
            if i==0:ax.set_title('GT' if key=='truth' else key,fontsize=10)
            if j==0:ax.set_ylabel('Fresh digit '+str(digit))
            else:
                row=next(r for r in scores if r['case_id']==case and r['method']==key)
                ax.set_xlabel('RMSE %.2f | IoU %.2f\nr %.2g | %.1fs | K%d'%(row['RMSE'],row['IoU'],row['relative_data_residual'],row['wall_seconds'],row['K']),fontsize=8)
    for j,(field,title,span) in enumerate([('wall_seconds','Solver wall, three targets (s)',slice(0,2)),('LU','Full 576 LU, three targets',slice(2,4)),('complex_LU_RHS','Physical LU/kernel complex RHS',slice(4,7))]):
        ax=fig.add_subplot(grid[3,span]);ax.bar(np.arange(6),[totals[m][field] for m in METHODS],color=COLORS)
        ax.set_xticks(np.arange(6),METHODS,rotation=35,ha='right',fontsize=8);ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.2)
    fig.suptitle('R4 frozen V3: all fresh reconstructions and paid costs\nCommon material color range [0.5,20]; all six methods fail the joint recovery gate; no GT-selected image',fontsize=13)
    for ext in ('png','pdf'):fig.savefig(figures/('reconstruction_and_cost.'+ext),dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(17,5.6),layout='constrained')
    for ax,case in zip(axes,cases):
        arrays={m:np.load(ROOT/'results/fresh_v3'/(case+'_'+m+'.trajectory.npz'))['m'] for m in METHODS}
        samples=np.vstack([a[np.linspace(0,len(a)-1,min(60,len(a))).astype(int)] for a in arrays.values()])
        reference=np.ones(576)
        _,sv,Vh=np.linalg.svd(samples-reference,full_matrices=False);basis=Vh[:2].T
        out=ROOT/'diagnostics'/('fixed_plot_projection_'+case+'.npz')
        np.savez_compressed(out,basis=basis,reference=reference,singular_values=sv,definition='PCA of frozen paid material trials only; GT excluded from basis construction')
        for method,color in zip(METHODS,COLORS):
            z=(arrays[method]-reference)@basis
            ax.plot(z[:,0],z[:,1],color=color,lw=.8,alpha=.7,label=method)
            endpoint=json.loads((ROOT/'published_results/fresh_v3'/(case+'_'+method+'.json')).read_text())
            selected=(np.asarray(endpoint['top1']['m'])-reference)@basis
            ax.scatter(*selected,color=color,s=30,marker='o')
        truth=np.load(ROOT/'truth_after_freeze/fresh'/(case+'.npz'))['m'];z=(truth-reference)@basis
        ax.scatter(*z,color='black',marker='*',s=130,label='GT, posthoc')
        digit=next(r['digit'] for r in scores if r['case_id']==case)
        ax.set_title('Fresh digit %d | PCA energy %.1f%%'%(digit,100*np.sum(sv[:2]**2)/np.sum(sv**2)),fontsize=11)
        ax.set_xlabel('Fixed material PC1');ax.set_ylabel('Fixed material PC2');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Actual paid material-trial trajectories in a fixed, non-invertible projection\nLines include rejected optimizer trials; circles = online top1; this is not a basin boundary, solution fiber or posterior',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(figures/('maxwell_branch_trajectories.'+ext),dpi=180)
    plt.close(fig)
    print('Rendered final reconstruction/cost and fixed-projection figures; zero new Maxwell solves.')

if __name__=='__main__':main()
