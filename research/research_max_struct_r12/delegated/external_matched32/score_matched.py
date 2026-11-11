"""Post-seal scorer for explicitly matched32 physics diagnostic."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).absolute().parent;sys.path.insert(0,str(HERE.parent/'external'));from run_external import dump
cpu=time.process_time();cid='austria_matched32';methods=['direct','Dual','DR'];paths=[HERE/'RESULTS'/cid/k/'derived.npz'for k in methods]
if not all(p.exists()and p.with_name('RECEIPT.json').exists()for p in paths):raise RuntimeError('All3 matched predictions required before scorerGT')
dump(HERE/'PREDICTION_SEAL.json',dict(status='ALL3_MATCHED_PREDICTIONS_SEALED_BEFORE_GT_SCORE',paths=[str(p.relative_to(HERE))for p in paths],qualification='OFFLINE_PHYSICS_DIAGNOSTIC generated after original exposed scores; not freshconfirmation',selection='registered original physical Phi only'))
truth=np.load(HERE/'sealed_truth/GT.npz')['m'];generation=json.loads((HERE/'GENERATION_RECEIPT.json').read_text());rows=[]
for k,p in zip(methods,paths):
    z=np.load(p);r=json.loads(p.with_name('RECEIPT.json').read_text());rows.append(dict(method=k,final_NRMSE=float(np.linalg.norm(z['m']-truth)/np.linalg.norm(truth)),raw_NRMSE=float(np.linalg.norm(z['raw']-truth)/np.linalg.norm(truth)),train_residual=r['train_relative_complex_residual'],final_Phi=r['original_Phi'],raw_Phi=r['raw_Phi'],CPU_seconds=r['CPU_seconds'],auxiliary=r['auxiliary'],corrector=r['corrector'],peak_chi=float(z['m'].max()),counts=r['physical_counts']))
dump(HERE/'SCORES.json',dict(status='MATCHED_OFFLINE_PHYSICS_DIAGNOSTIC_SCORED',rows=rows,coarse_GT_true_residual=generation['coarse_GT_observed_relative_residual'],new_physical_solves=0,original150to32_unchanged=True))
fig,axs=plt.subplots(2,4,figsize=(13,7),constrained_layout=True)
for i,key in enumerate(['raw','m']):
    for j,k in enumerate(['truth']+methods):
        m=truth if k=='truth'else np.load(HERE/'RESULTS'/cid/k/'derived.npz')[key];im=axs[i,j].imshow(m.reshape(32,32),origin='lower',vmin=0,vmax=20,cmap='magma');axs[i,j].set_xticks([]);axs[i,j].set_yticks([]);axs[i,j].set_title(k+' / '+('selected auxiliary'if i==0 else'common corrector'))
        if k!='truth':r=next(x for x in rows if x['method']==k);axs[i,j].set_xlabel(f"NRMSE {r['raw_NRMSE'if i==0 else'final_NRMSE']:.3f}")
cb=fig.colorbar(im,ax=axs,label='chi full0–20 box; targetchi3',shrink=.7);fig.suptitle('Matched32 Austria offline physics diagnostic; fresh10%fixednormnoise, no radius/scale tuning')
fig.savefig(HERE/'MATCHED32_MATERIALS_FULL_RANGE.png',dpi=160)
for ax in axs.flat:
    for im in ax.images:im.set_clim(0,4)
cb.set_label('chi display0–4; inversion box0–20, targetchi3')
fig.savefig(HERE/'MATCHED32_MATERIALS_DETAIL.png',dpi=160);plt.close(fig)
queue=json.loads((HERE/'QUEUE.json').read_text());dump(HERE/'COST_LEDGER.json',dict(status='COMPLETE',external_child_CPU_seconds=sum(x['external_child_CPU_seconds']for x in queue['rows']),generation_CPU_seconds=generation['CPU_seconds'],generation_new_LU=0,generation_reused_clean_forward='external_physics/RESULTS.json sole1LU20RHS accounted upstreamonce',rows=queue['rows'],physical_counts=[dict(method=x['method'],counts=x['counts'])for x in rows],scoring_CPU_near_exit=time.process_time()-cpu))
lines=['# Matched32 Austria: bounded offline physics diagnostic','', 'Original Austria data were actually generated on150² and inverted on32². Its coarseGT residual.2252565 exceeds the unchanged.12 ball. This separately authorized matched32 control removes that discretization mismatch using the same32²GT and operators with fresh10%fixednorm noise, seed2026111232. The clean prediction was reused from the sole authorized offlineGT forward; no redundant forward generation was performed. Public noisy observations and generatorGT are separate. The true matched target residual is'+f"{generation['coarse_GT_observed_relative_residual']:.6f}"+' and lies inside.12.','', 'This is OFFLINE_PHYSICS_DIAGNOSTIC after exposed original results, not original-paper reproduction or new confirmation. Same realbox0..20,chi=epsilon−1,targetchi3,initial0,prior,scale6,rawball.12 and180aux/60corrector CPU caps. Full native Dual uses inherited300perstage allowances,mf5/8/10,beta8/2/.5,currentbeta8,DPCS15. Direct/DR use the same previous2000physical/5000projection allowances; commoncorrector max1000. No contrast/radius/DPCS/scale retuning. All3 predictions were sealed before matchedscorerGT.','', '| Method | Raw NRMSE | Corrected NRMSE | Residual | Phi | CPU s |','|---|---:|---:|---:|---:|---:|']
for r in rows:lines.append(f"|{r['method']}|{r['raw_NRMSE']:.6f}|{r['final_NRMSE']:.6f}|{r['train_residual']:.6f}|{r['final_Phi']:.7g}|{r['CPU_seconds']:.2f}|")
lines+=['', 'This finite diagnostic separates model mismatch from the remaining search/objective ambiguity. Feasible true material does not imply a nonconvex method will reach it, and lower registeredPhi alone does not establish material recovery. All original150→32 outcomes remain unchanged, including bad material. All raw/corrected matched images and exact stops/costs are preserved.']
(HERE/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status='COMPLETE',rows=rows),indent=2))
