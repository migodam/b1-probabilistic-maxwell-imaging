"""Publication panels from already-sealed TEST scoring; no candidate selection."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ARMS=['v_image','v_physics','j_full','j_null','j_null_current']
LABELS=['Image NN','Voxel + physics','Full current','Null current','Null current loss']

def main():
    p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    s=json.loads(a.scores.read_text());study=a.scores.resolve().parents[1];a.output.mkdir(parents=True,exist_ok=True)
    seal=json.loads(Path(s['seal_path']).read_text());data=Path(seal['data_root'])
    records=json.loads(Path(seal['manifest']['path']).read_text())['records']
    rows={(r['run'],r['obs_id']):r for r in s['rows']}
    runs=['continue_trf','weakreg45_anchor']+[f'confirm1024_{x}_seed7'for x in ARMS]
    labels=['Continue TRF','Weak-reg TRF']+LABELS
    chosen={}
    for r in records:
        if r.get('level','low')!='low' or r['noise']!=0:continue
        key=f"digit{r['digit']}"if r['family']=='digits'else r['family']
        if key not in chosen:chosen[key]=r['parent_id']
    keys=[f'digit{x}'for x in range(1,10)]+['primitive_composites','islands','smooth']
    selection=[]
    for noise in (0,.001):
        selected=[(key,next(r for r in records if r['parent_id']==chosen[key]and r['noise']==noise and r.get('level','low')=='low'))for key in keys if key in chosen]
        fig,axes=plt.subplots(len(selected),2+len(runs),figsize=(18,1.9*len(selected)),squeeze=False,layout='constrained')
        for i,(key,r)in enumerate(selected):
            obs=r['obs_id']
            with np.load(data/r['truth_path'])as z:gt=z['m'].reshape(24,24)
            with np.load(study/'anchors'/f'{obs}.npz')as z:gn=z['m0'].reshape(24,24)
            panels=[gt,gn];captions=['','']
            for run in runs:
                metric=rows.get((run,obs),{}).get('top1')
                path=study/'evaluation'/run/f'{obs}.npz'
                arr=None
                if path.exists():
                    with np.load(path)as z:
                        if 'top1'in z:arr=z['top1'].reshape(24,24)
                panels.append(arr)
                captions.append('NOT RUN / FAIL'if not metric else f"N {metric['NRMSE']:.3f} | IoU {metric['IoU_Otsu']:.2f}")
            for j,(arr,caption)in enumerate(zip(panels,captions)):
                ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                if arr is not None:
                    im=ax.imshow(arr,vmin=.02,vmax=.55,cmap='viridis',interpolation='nearest')
                    sat=int(np.count_nonzero((arr<.02)|(arr>.55)))
                    if sat:ax.text(.02,.02,f'sat {sat}',transform=ax.transAxes,fontsize=6,color='white',backgroundcolor='black')
                else:ax.text(.5,.5,'FAIL',ha='center',transform=ax.transAxes)
                if caption:ax.set_xlabel(caption,fontsize=7)
                if i==0:ax.set_title((['GT','GN anchor']+labels)[j],fontsize=9)
                if j==0:ax.set_ylabel(key,fontsize=9)
            selection.append(dict(obs_id=obs,parent_id=r['parent_id'],key=key,noise=noise))
        fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.45,label='Material contrast; common [.02,.55], saturation marked',extend='both')
        fig.suptitle(f"Independent TEST | {'clean'if noise==0 else'noise 0.001'} | first manifest parent per digit/family | online top1",fontsize=14)
        fig.savefig(a.output/('test_clean_reconstructions.png'if noise==0 else'test_noisy_reconstructions.png'),dpi=155);plt.close(fig)
    (a.output/'REPORT_PANEL_SELECTION.json').write_text(json.dumps(dict(rule='first manifest low-clean parent per digit/family; paired noisy same parent; no metric selection',cases=selection),indent=2))
    fig,axes=plt.subplots(2,3,figsize=(17,8),layout='constrained')
    for col,(level,noise,title)in enumerate([('low',0,'Low contrast / clean'),('low',.001,'Low contrast / noisy'),('stress',None,'A=2.5 extrapolation')]):
        vals=[];success=[]
        for run in runs:
            rr=[r for r in s['rows']if r['run']==run and r.get('level','low')==level and(noise is None or r['noise']==noise)]
            vals.append([r['top1']['NRMSE']for r in rr if r.get('top1')]);success.append(sum(r.get('joint_success',False)for r in rr)/max(1,len(rr)))
        axes[0,col].boxplot(vals,tick_labels=labels,showfliers=True)
        axes[0,col].axhline(.15,ls='--',color='gray');axes[0,col].set(title=title,ylabel='Selected material NRMSE / amplitude',yscale='log')
        axes[1,col].bar(labels,success);axes[1,col].set(ylabel='Joint success rate',ylim=(0,1))
        for ax in axes[:,col]:ax.tick_params(axis='x',rotation=45,labelsize=8)
    fig.savefig(a.output/'quality_and_joint_success.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.5),layout='constrained')
    for run,label in zip(runs,labels):
        rr=[r for r in s['rows']if r['run']==run and r.get('level','low')=='low'and r.get('top1')]
        scope='45s independent anchor'if run=='weakreg45_anchor'else'post shared GN'
        axes[0].scatter([r['wall_seconds']for r in rr],[r['top1']['NRMSE']for r in rr],s=13,alpha=.5,label=label+' ('+scope+')')
        axes[1].scatter([r['costcounts'].get('LU',0)for r in rr],[r['top1']['NRMSE']for r in rr],s=13,alpha=.5)
        axes[2].scatter([r['top1']['relative_residual']for r in rr],[r['top1']['NRMSE']for r in rr],s=13,alpha=.5)
    for ax in axes:ax.set_yscale('log');ax.set_ylabel('Material NRMSE / amplitude')
    axes[0].set_xlabel('Measured charged wall seconds (different scopes labeled)');axes[1].set_xlabel('Charged LU');axes[2].set(xlabel='Full VIE relative data residual',xscale='log')
    axes[0].legend(fontsize=6);fig.savefig(a.output/'cost_quality_physical.png',dpi=175);plt.close(fig)

if __name__=='__main__':main()
