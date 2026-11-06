"""Offline figures from saved predictions and receipts; never runs physics or optimization."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[_k]='1'
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]
OLD=HERE.parent/'research_z_alg_r4_diversity'
OUT=HERE/'figures'
ARMS=['continue_trf','haar_unfiltered','green_pullback','soft_jacobian']
LABELS={'continue_trf':'Continue TRF','haar_unfiltered':'Haar','green_pullback':'Green','soft_jacobian':'Soft J'}
COLORS=['#686868','#d68a24','#278e72','#6259ad']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
def csvrows(path):
    return list(csv.DictReader(path.open())) if path.exists() else []
def records(folder):
    out=[]
    for p in sorted((HERE/'results'/folder).glob('*.json')):
        if not p.with_suffix('.npz').exists():continue
        r=json.loads(p.read_text())
        if 'result' in r:r=r['result']
        if 'task' in r and 'top1' in r:out.append((p,r))
    return out

def save(fig,name):
    OUT.mkdir(exist_ok=True);fig.savefig(OUT/name,bbox_inches='tight',facecolor='white');plt.close(fig)

def geometry():
    folder=HERE/'results/e1'; dims=csvrows(folder/'dimensions.csv');cov=csvrows(folder/'coverage.csv');cases=csvrows(folder/'cases.csv')
    if not (dims and cov and cases):return
    fig,axs=plt.subplots(2,2,figsize=(13,9))
    for label,color in [('Born','#aaa'),('public_mref','#3988b0'),('anchor','#d88036'),('offline_GT','#32956f')]:
        for eps,style in [(1e-5,':'),(1e-4,'-'),(1e-3,'--')]:
            rs=[r for r in dims if r['linearization']==label and np.isclose(float(r['RMS']),.03) and np.isclose(float(r['epsilon']),eps)]
            x=sorted(set(int(r['sources']) for r in rs));y=[np.median([float(r['d_weak']) for r in rs if int(r['sources'])==n]) for n in x]
            axs[0,0].plot(x,y,style,color=color,label=f'{label}, ε={eps:g}')
    axs[0,0].set(xlabel='Illuminations (fixed receiver array)',ylabel='Tolerance-weak dimensions / 576',title='A  |  Weak dimension depends on noise tolerance (RMS=0.03)',ylim=(0,590));axs[0,0].legend(fontsize=7,ncol=2)
    for label,color in [('Born','#aaa'),('public_mref','#3988b0'),('anchor','#d88036'),('offline_GT','#32956f')]:
        rs=[r for r in dims if r['linearization']==label and float(r['RMS'])==.03 and float(r['epsilon'])==1e-4]
        x=sorted(set(int(r['sources']) for r in rs));y=[np.median([float(r['numerical_null_dimension']) for r in rs if int(r['sources'])==n]) for n in x]
        axs[0,1].plot(x,y,'o-',label=label,color=color)
    subset=csvrows(folder/'current_subset_dimensions.csv')
    if subset:
        x=sorted(set(int(r['sources'])for r in subset));y=[np.median([float(r['numerical_kernel_dimension'])for r in subset if int(r['sources'])==n])for n in x]
        axs[0,1].plot(x,y,'s--',color='black',label='Anchor current-leakage material kernel')
    axs[0,1].set(xlabel='Illuminations',ylabel='Numerical null/kernel dimension / 576',title='B  |  Distinct maps, separate relative rank threshold 1e−12',ylim=(-10,590));axs[0,1].legend(fontsize=8)
    chosen={r['task_id'] for r in cases if r['method']=='native_trf' and np.isclose(float(r['amplitude']),.4)}
    cumulative=HERE/'results/coverage_cumulative.npz'
    coverage_dimensions=csvrows(HERE/'results/coverage_dimension.csv')
    if not cumulative.exists():
        raise FileNotFoundError('Full cumulative coverage is required; sparse-rank interpolation is not used.')
    with np.load(cumulative,allow_pickle=False)as z:
        curves={str(key):fraction for key,fraction in zip(z['key'],z['fraction'])}
    table=[]
    for basis,label,color in [('anchor_J','Anchor J','#d88036'),('anchor_Green','Anchor Green','#32956f'),('reference_J','Fixed-reference J','#3988b0'),('reference_Green','Fixed-reference Green','#8058a3')]:
        ys=[curves[tid+'|'+basis]for tid in sorted(chosen)if tid+'|'+basis in curves]
        if not ys:continue
        axs[1,0].plot(np.arange(1,577),np.median(ys,axis=0),label=label,color=color,linewidth=1.8)
        rs=[r for r in coverage_dimensions if r['task_id']in chosen and r['basis']==basis]
        if rs:table.append([label,f"{np.median([float(r['r90'])for r in rs]):g}",f"{np.median([float(r['r95'])for r in rs]):g}"])
    axs[1,0].plot([0,576],[0,1],':',color='#999',label='Rank/576 reference')
    axs[1,0].set(xlabel='Retained material directions (all ranks 1–576)',ylabel='Fraction of squared GT error covered',title=f'C  |  Full cumulative coverage, native {len(chosen)} anchors (A=0.4)',ylim=(0,1.02))
    axs[1,0].legend(fontsize=7,loc='lower right')
    if table:
        tab=axs[1,0].table(cellText=table,colLabels=['Native cohort medians','r90','r95'],cellLoc='center',colWidths=[.64,.18,.18],bbox=[.035,.60,.58,.33])
        tab.auto_set_font_size(False);tab.set_fontsize(7)
        for cell in tab.get_celld().values():cell.set_edgecolor('#dddddd');cell.set_linewidth(.5)
    for field,label,color in [('material_approx_rank','Green ≤0.01 leakage','#32956f'),('numerical_material_kernel_dimension','Material numerical kernel','#ba5252')]:
        vals=[int(r[field]) for r in cases];axs[1,1].plot(range(len(vals)),vals,'.',label=f'{label}: {min(vals)}–{max(vals)}',color=color)
    axs[1,1].axhline(509,color='#888',linestyle=':',label='Fixed current dark dimension=509')
    axs[1,1].set(xlabel='Frozen anchor index (both historical methods)',ylabel='Material dimension',title=f'D  |  Current-space darkness versus material pullback ({len(cases)} anchors)',ylim=(-10,590));axs[1,1].legend(fontsize=8)
    fig.suptitle('GN–Green R2: saved geometry and offline error labels',fontsize=14);fig.tight_layout();save(fig,'search_space_and_gt_coverage.png')

def finite_and_correction():
    rec=records('e3'); proposals=[(r['task']['arm'],p) for _,r in rec for p in r.get('proposals',[])]
    if not proposals:return
    fig,axs=plt.subplots(2,2,figsize=(13,9))
    for arm,c in zip(ARMS,COLORS):
        ps=[p for a,p in proposals if a==arm]
        if not ps:continue
        x=np.array([p['raw_RMS'] for p in ps]);y=np.array([p['filtered_RMS'] for p in ps]);axs[0,0].scatter(x,y,s=13,alpha=.35,color=c,label=f'{LABELS[arm]} ({len(ps)} rays)')
    axs[0,0].set(yscale='symlog',xlabel='Raw proposal RMS',ylabel='After filter + scalar box shrink: RMS',title='A  |  Collapsed box-feasible rays are retained');axs[0,0].set_yscale('symlog',linthresh=1e-12);axs[0,0].legend(fontsize=8)
    labels=[];rates=[];med=[]
    for arm in ARMS[1:]:
        ps=[p for a,p in proposals if a==arm]
        if ps:labels.append(LABELS[arm]);rates.append(np.mean([p['filtered_RMS']<1e-12 for p in ps]));med.append(np.median([p['ray_scale'] for p in ps]))
    axs[0,1].bar(labels,rates,color=COLORS[1:]);axs[0,1].set(ylabel='Fraction with feasible RMS < 1e−12',ylim=(0,1.05),title='B  |  Box collapse rate, all completed proposal rows')
    for i,(r,m) in enumerate(zip(rates,med)):axs[0,1].text(i,r+.025,f'median ray scale={m:.1g}',ha='center',fontsize=8)
    retention=csvrows(HERE/'results/TRANSVERSE_RETENTION.csv')
    if not retention:raise FileNotFoundError('Saved transverse-retention rows are required for panel C.')
    xs=np.arange(3);ticks=[]
    for i,arm in enumerate(ARMS[1:]):
        counts=[]
        for component,offset,color in [('weak',-.18,'#376fa0'),('strong',.18,'#c47c35')]:
            vals=[]
            for row in retention:
                if row['arm']!=arm or float(row['seed_RMS'])<=1e-12:continue
                value=row[component+'_retention']
                if not value or float(row[component+'_before'])<=1e-12:continue
                if np.isfinite(float(value)):vals.append(float(value))
            counts.append(len(vals))
            if vals:
                q1,median,q3=np.quantile(vals,[.25,.5,.75])
                axs[1,0].bar(i+offset,median,width=.32,color=color,label=component.capitalize()+' J component'if i==0 else None)
                axs[1,0].errorbar(i+offset,median,yerr=np.array([[median-q1],[q3-median]]),fmt='none',color='#333',capsize=4,linewidth=1)
                axs[1,0].text(i+offset,median+.035,f'{median:.3f}',ha='center',fontsize=8)
        ticks.append(f'{LABELS[arm]}\nn weak/strong={counts[0]}/{counts[1]}')
    axs[1,0].axhline(1,color='#999',linestyle=':')
    axs[1,0].set_xticks(xs,ticks)
    axs[1,0].set(ylabel='After/seed anchor-displacement component norm',title='C  |  Actual TRF retention: median and interquartile range',ylim=(0,1.45))
    axs[1,0].legend(fontsize=8,loc='upper right')
    axs[1,0].text(.02,.95,'Exclude seed RMS ≤ 1e−12 and component denominator ≤ 1e−12',transform=axs[1,0].transAxes,fontsize=7,va='top')
    finite=csvrows(HERE/'results/finite/finite_steps.csv')
    if finite:
        for basis in sorted(set(r['basis'] for r in finite)):
            rs=[r for r in finite if r['basis']==basis];x=sorted(set(float(r['requested_RMS']) for r in rs))
            y=[np.median([float(r['relative_data_drift'])/float(r['epsilon']) for r in rs if float(r['requested_RMS'])==radius]) for radius in x]
            axs[1,1].plot(x,y,'o-',label=basis)
        axs[1,1].axhline(1,color='#999',linestyle=':',label='Drift budget')
        axs[1,1].set(xscale='log',yscale='log',xlabel='Requested projected GT-error RMS',ylabel='Nonlinear relative data drift / ε',title='D  |  Offline GT-directed finite-step ceiling (median)');axs[1,1].legend(fontsize=6)
    else:
        axs[1,1].axis('off');axs[1,1].text(.05,.6,'Finite-step diagnostic not yet available.\nPanels A–C use saved E3 proposal receipts only.',transform=axs[1,1].transAxes)
    fig.suptitle('GN–Green R2: finite proposals, box shrink, and nonlinear correction',fontsize=14);fig.tight_layout();save(fig,'finite_steps_and_correction.png')

def reconstructions():
    rec=records('e3');clean=[(p,r) for p,r in rec if r['task']['noise']=='clean' and float(r['task']['amplitude'])==.4]
    if len(clean)<48:return
    shapes=sorted(set(r['task']['shape'] for _,r in clean)); lookup={(r['task']['shape'],r['task']['arm']):(p,r)for p,r in clean}
    fig,axs=plt.subplots(len(shapes),7,figsize=(14,len(shapes)*1.9),squeeze=False)
    titles=['GT','Original native']+[LABELS[a]+' top1' for a in ARMS]+['Green bestK\nOFFLINE oracle']
    for i,shape in enumerate(shapes):
        p,r=lookup[(shape,'green_pullback')];truth=np.load(OLD/'truth_after_freeze'/(r['task']['observation_id']+'.npz'))['m'];z=np.load(p.with_suffix('.npz'));anchor=z['anchor']
        candidates=np.concatenate([z['anchor'][None],z['proposal_before'],z['proposal_after'],z['top1'][None]])
        best=candidates[np.argmin(np.mean((candidates-truth)**2,axis=1))]
        arrays=[truth,anchor]+[np.load(lookup[(shape,a)][0].with_suffix('.npz'))['top1']for a in ARMS]+[best]
        for k,m in enumerate(arrays):
            image=(m-.02)/.4;im=axs[i,k].imshow(image.reshape(24,24),cmap='viridis',vmin=0,vmax=1,interpolation='nearest');axs[i,k].set_xticks([]);axs[i,k].set_yticks([])
            if i==0:axs[i,k].set_title(titles[k],fontsize=9)
            if k==0:axs[i,k].set_ylabel(shape,fontsize=9)
            else:axs[i,k].text(.02,.02,f'RMSE {np.sqrt(np.mean((m-truth)**2)):.4f}',transform=axs[i,k].transAxes,color='white',fontsize=7,bbox={'facecolor':'black','alpha':.45,'pad':1})
    fig.suptitle('All 12 low-amplitude clean shapes | (m − 0.02)/0.4 fixed display [0,1]\nTop1 uses public nonlinear objective; bestK uses disclosed GT offline only',fontsize=12,y=.995)
    fig.subplots_adjust(left=.05,right=.91,bottom=.01,top=.96,wspace=.08,hspace=.16);fig.colorbar(im,ax=axs.ravel().tolist(),fraction=.012,pad=.012,label='(m − background) / amplitude');save(fig,'reconstruction_and_cost.png')

def quality_cost():
    rec=records('e3')
    if not rec:return
    fig,axs=plt.subplots(1,3,figsize=(15,4.6));all_ratios=[]
    for arm,c in zip(ARMS,COLORS):
        vals=[]
        for p,r in rec:
            if r['task']['arm']!=arm:continue
            gt=np.load(OLD/'truth_after_freeze'/(r['task']['observation_id']+'.npz'))['m'];z=np.load(p.with_suffix('.npz'));base=np.sqrt(np.mean((z['anchor']-gt)**2));top=np.sqrt(np.mean((z['top1']-gt)**2))
            candidates=[z['anchor'],z['top1']]
            for key in ('proposal_before','proposal_after'):
                if key in z:candidates.extend(z[key])
            best=min(np.sqrt(np.mean((v-gt)**2))for v in candidates);vals.append((r['wall_seconds'],r['cpu_seconds'],top/base,best/base,r['costcounts'].get('LU',0)))
        if not vals:continue
        v=np.array(vals);all_ratios.extend(v[:,2:4].tolist());axs[0].scatter(v[:,0],v[:,2],s=28,color=c,alpha=.65,label=LABELS[arm]);axs[1].scatter(v[:,1],v[:,2],s=28,color=c,alpha=.65);axs[2].scatter(v[:,2],v[:,3],s=28,color=c,alpha=.65)
    for ax in axs:ax.axhline(1,color='#aaa',linestyle=':')
    axs[0].set(xlabel='Per-task additional wall seconds',ylabel='Top1 RMSE / frozen-anchor RMSE',title='Public-objective top1 versus wall cost');axs[0].legend(fontsize=8)
    axs[1].set(xlabel='Per-task CPU seconds',ylabel='Top1 RMSE / frozen-anchor RMSE',title='CPU cost (operator setup included)')
    rr=np.array(all_ratios);xmin,xmax=float(rr[:,0].min()),float(rr[:,0].max());ymin=float(rr[:,1].min());pad=max((xmax-xmin)*.12,.003)
    axs[2].plot([min(xmin,ymin)-pad,xmax+pad],[min(xmin,ymin)-pad,xmax+pad],':',color='#999');axs[2].set_xlim(xmin-pad,xmax+pad);axs[2].set_ylim(ymin-.01,1.01);axs[2].set(xlabel='Top1 RMSE / anchor RMSE',ylabel='Offline bestK RMSE / anchor RMSE',title='Oracle candidate diversity ≠ online success')
    fig.suptitle(f'GN–Green R2: {len(rec)} saved E3 tasks, clean + noisy (offline scoring)',fontsize=13);fig.tight_layout();save(fig,'quality_vs_cost.png')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--only',choices=['geometry','finite','reconstruction','cost','all'],default='all');args=parser.parse_args()
    funcs={'geometry':geometry,'finite':finite_and_correction,'reconstruction':reconstructions,'cost':quality_cost}
    for name,fn in funcs.items():
        if args.only in ('all',name):fn()
if __name__=='__main__':main()
