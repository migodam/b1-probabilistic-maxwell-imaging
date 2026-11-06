"""Offline fixed-scale figures after all nine seals and material scoring."""
import os
os.environ['MPLBACKEND']='Agg'
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d';FIG=BASE/'FIGURES'
NAMES={'sourceCIE':'Source CIE','sourceDual':'Source Dual','GlobalDual':'Global Dual','AplusGlobal':'R9 vector A+','LSIEsource':'Source LSIE'}
ORDER=['sourceCIE','sourceDual','GlobalDual','AplusGlobal','LSIEsource']

def main():
    wall=time.perf_counter();scores=json.loads((OUT/'NINE_OFFLINE_SCORES.json').read_text())
    if scores['status']!='OFFLINE_SCORED_AFTER_ALL_NINE_SEALED':raise RuntimeError('All9 seals and offline scores required')
    if len(json.loads((OUT/'SEALED_ALL_NINE.json').read_text())['terminals'])!=9:raise RuntimeError('Missing seals')
    protocol=json.loads((BASE/'configs/original_cm_nine.json').read_text())
    truth=np.load(OUT/'sealed_truth/truth.npz')['inverse_material'];n=30;z=(np.arange(n)-(n-1)/2)*.1
    fixed_z=[int(np.argmin(abs(z-v))) for v in [-.7,.7]];FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11})
    colors={a:plt.rcParams['axes.prop_cycle'].by_key()['color'][i] for i,a in enumerate(ORDER)}
    artifacts=[]
    for cohort in protocol['cohorts']:
        name=cohort['name'];arms=cohort['arms'];nc=len(arms)+1
        fig,axes=plt.subplots(3,nc,figsize=(3*nc,8.6),sharex=True,sharey=True,layout='constrained')
        materials=[truth]+[np.load(OUT/'cohorts'/name/'inversions'/a/'stage2.npz')['material'] for a in arms]
        for col,(label,m) in enumerate(zip(['Truth']+[NAMES[a] for a in arms],materials)):
            cube=m.reshape(n,n,n);views=[cube[:,:,i].T for i in fixed_z]+[cube.max(axis=2).T]
            for row,v in enumerate(views):
                im=axes[row,col].imshow(v,origin='lower',extent=[-1.5,1.5,-1.5,1.5],vmin=0,vmax=4,cmap='viridis',interpolation='nearest')
                axes[row,col].set_xticks([-1,0,1]);axes[row,col].set_yticks([-1,0,1])
                if row==2:axes[row,col].set_xlabel('x / λ')
            if col:
                score=next(r for r in scores['primary_rows'] if r['cohort']==name and r['arm']==arms[col-1])
                label+='\nErr %.3f; χ-NRMSE %.3f\nΦ %.4f'%(score['paper_Err'],score['contrast_NRMSE'],score['Phi'])
            axes[0,col].set_title(label)
        for row,label in enumerate(['z / λ = %.2f'%z[fixed_z[0]],'z / λ = %.2f'%z[fixed_z[1]],'Maximum along z']):axes[row,0].set_ylabel(label+'\ny / λ')
        fig.colorbar(im,ax=axes,shrink=.85,label='Real contrast χ = εᵣ−1 (common0–4 box)')
        caption='60³ data /30³ inverse; model discrepancy0.535' if name=='cross_grid' else '30³ data /30³ inverse; matched finite control'
        fig.suptitle('PUBLIC DEVELOPMENT · CM+RR3D DDA · original3λ four cubes\n'+caption+' · final-stage endpoints, no GT iterate selection',fontsize=14)
        for ext in ['png','pdf']:
            dest=FIG/('original_CM_'+name+'_materials.'+ext);fig.savefig(dest,dpi=190);artifacts.append(str(dest.relative_to(BASE)) if dest.is_relative_to(BASE) else str(dest))
        plt.close(fig)
    fig,axes=plt.subplots(4,2,figsize=(13,14),layout='constrained')
    for column,cohort in enumerate(protocol['cohorts']):
        name=cohort['name']
        for arm in cohort['arms']:
            dest=OUT/'cohorts'/name/'inversions'/arm
            eps=json.loads((dest/'ENDPOINTS.json').read_text());paid=[r for r in scores['all_paid_stage_rows'] if r['cohort']==name and r['arm']==arm]
            xx=[e['cumulative_wall_seconds']/60 for e in eps]
            yy=[[next(r for r in paid if r['endpoint']==e['name'])[k] for e in eps] for k in ['paper_Err','contrast_NRMSE','physical_data_RMS']]
            for row,v in enumerate(yy):axes[row,column].plot(xx,v,'-o',color=colors[arm],label=NAMES[arm],linewidth=1.5,markersize=4)
            hist=[json.loads(x) for x in (dest/'history.jsonl').read_text().splitlines()]
            good=[r for r in hist if r['coefficient_delta_RMS'] is not None and r['coefficient_delta_RMS']>0]
            axes[3,column].semilogy([r['cumulative_wall_seconds']/60 for r in good],[r['coefficient_delta_RMS'] for r in good],color=colors[arm],linewidth=1.1,label=NAMES[arm])
        axes[0,column].set_title('Cross-grid CM60→30' if name=='cross_grid' else 'Matched finite CM30→30')
        gt=scores['witnesses'][name]['Phi_at_same_known_feasible_material']
        axes[2,column].axhline(np.sqrt(gt),color='black',linestyle=':',label='Known GT feasible response RMS')
        axes[3,column].axhline(.001,color='black',linestyle=':',label='Eq16 stop0.001')
        for row,label in enumerate(['Paper Err (ε relative voxel RMS)','Contrast NRMSE','Physical data RMS','Coefficient change RMS']):
            axes[row,column].set_ylabel(label);axes[row,column].set_xlabel('Actual wall minutes (setup amortized)');axes[row,column].grid(alpha=.2)
        axes[0,column].legend(fontsize=10);axes[2,column].legend(fontsize=8)
    fig.suptitle('PUBLIC DEVELOPMENT · stored stage endpoints and paid coefficient traces\nActual shared queue wall cost; host CPU excludes GPU arithmetic; no standalone speed claim',fontsize=13)
    for ext in ['png','pdf']:
        dest=FIG/('original_CM_nine_convergence.'+ext);fig.savefig(dest,dpi=180);artifacts.append(str(dest.relative_to(BASE)) if dest.is_relative_to(BASE) else str(dest))
    plt.close(fig)
    (OUT/'NINE_PLOT_COST.json').write_text(json.dumps(dict(CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,
        platform='Mac',physics_evaluations=0,scope='Offline plots from all9sealed predictions/scores; common0..4scale; fixedz±.7 and max projection',
        additive=True,included_in=None,artifacts=artifacts),indent=2)+'\n')

if __name__=='__main__':main()
