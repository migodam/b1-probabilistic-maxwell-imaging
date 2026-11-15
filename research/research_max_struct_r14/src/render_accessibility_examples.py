"""Post-seal explanatory material figures; full24 statistics remain unchanged."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,csv,json,time
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve()
    auth=json.loads((r/'audit/PARENT_INDEPENDENT24_AUTHORIZATION.json').read_text())
    if not auth['authorized_scoring']:raise PermissionError('Post-seal parent approval required')
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    c=time.process_time();w=time.perf_counter();out=r/'FIGURES/independent24_examples';out.mkdir(exist_ok=False)
    rows=json.loads((r/'analysis/independent24_v1/ALL_VIEWS.json').read_text());scores={(x['case_id'],x['solver']):x for x in rows if x['view']=='primary'}
    arrays={};truths={}
    for i in range(24):
        cid=f'r14_new_{i:04d}'
        with np.load(r/'sealed_truth/independent24'/(cid+'.npz'),allow_pickle=False) as z:truths[cid]=z['material_coarse'].copy()
        for m in ['cie','dual_cie','pair01','pair11','g19','dr']:
            with np.load(r/'RESULTS/independent24_v1'/(cid+'__'+m)/'primary.npz',allow_pickle=False) as z:arrays[(cid,m)]=z['m'].copy()
    def gallery(ids,methods,name,title):
        fig=plt.figure(figsize=(3*(len(methods)+1)+.7,3.25*len(ids)),layout='constrained')
        grid=fig.add_gridspec(len(ids),len(methods)+2,width_ratios=[1]*(len(methods)+1)+[.055])
        labels={'cie':'CIE','dual_cie':'Dual-CIE','pair01':'A+ (01)','pair11':'Nonlocal consistent (11)','g19':'G19','dr':'DR'}
        for row,i in enumerate(ids):
            cid=f'r14_new_{i:04d}';items=[truths[cid]]+[arrays[(cid,m)] for m in methods];hi=max(float(x.max()) for x in items)
            for col,m in enumerate(items):
                ax=fig.add_subplot(grid[row,col]);im=ax.imshow(m.reshape(24,24),vmin=.02,vmax=hi,cmap='viridis',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title('Offline GT' if col==0 else labels[methods[col-1]],fontsize=11)
                if col==0:ax.set_ylabel(cid,fontsize=10)
                else:
                    s=scores[(cid,methods[col-1])];ax.set_xlabel(f"NRMSE {s['NRMSE']:.3f} | IoU {s['IoU']:.3f}\nresidual / delta {s['residual_over_delta']:.3f}",fontsize=9)
            cb=fig.colorbar(im,cax=fig.add_subplot(grid[row,-1]));cb.ax.tick_params(labelsize=9)
        fig.suptitle(title+'\nPost-seal explanatory examples; all24 cases and failures retained separately. Row-shared chi scale.',fontsize=13)
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=160)
        plt.close(fig)
    gallery([12,15,21,9,20],['cie','dual_cie','pair01','pair11','g19','dr'],'ACCESSIBILITY_COUNTEREXAMPLES','Different metric placements and reflection access different materials')
    gallery([2,19],['pair01','dr'],'THRESHOLD_SENSITIVE_CASES','Continuous recovery versus a fixed support threshold')
    differences=[]
    for cid,t in truths.items():
        x=arrays[(cid,'pair01')];y=arrays[(cid,'dr')];amp=float(t.max()-.1)
        differences.append(dict(case_id=cid,pair_material_RMSE=float(np.sqrt(np.mean((x-y)**2))),
            pair_amplitude_normalized_RMSE=float(np.sqrt(np.mean((x-y)**2)))/amp,
            pair_material_Linf=float(np.max(np.abs(x-y))),support_disagreement_cells=int(np.count_nonzero((x>1)!=(y>1))),
            Aplus_NRMSE=scores[(cid,'pair01')]['NRMSE'],DR_NRMSE=scores[(cid,'dr')]['NRMSE'],
            Aplus_IoU=scores[(cid,'pair01')]['IoU'],DR_IoU=scores[(cid,'dr')]['IoU']))
    with (out/'A_PLUS_DR_MATERIAL_DIFFERENCES.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(differences[0]));wr.writeheader();wr.writerows(differences)
    rec=dict(status='POST_SEAL_EXPLANATORY_FIGURES_COMPLETE',selection_changed=False,all_statistics_denominator=24,
        case_choice_scope='Posthoc mechanism and threshold illustrations,not a newtest cohort',new_forward=0,new_optimizer=0,
        material_views='Original minPhi primary; local explicit m-only transport value unchanged; fullfields in raw release',
        CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w)
    (out/'RECEIPT.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))


if __name__=='__main__':main()
