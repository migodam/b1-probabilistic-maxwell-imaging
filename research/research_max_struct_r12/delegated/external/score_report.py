"""Seal all external predictions, then evaluate heldout predictive fit/ Austria GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).absolute().parent;sys.path.insert(0,str(HERE));from run_external import CASES,METHODS,R8,load_file,dump
cpu=time.process_time();wall=time.perf_counter();paths=[HERE/'RESULTS'/c/k/'derived.npz' for c in CASES for k in METHODS]
if not all(p.exists() and p.with_name('RECEIPT.json').exists() for p in paths):raise RuntimeError('All12 frozen predictions must exist before heldout/GT evaluation')
dump(HERE/'PREDICTION_SEAL.json',dict(status='ALL12_PREDICTIONS_SEALED_BEFORE_HELDOUT_GT_SCORE',paths=[str(p.relative_to(HERE)) for p in paths],selection='Training original Phi only',heldout_selection=False))
rows=[]
for cid in CASES:
    austria=cid.startswith('austria')
    if austria:
        truthpath=R8/'data/paper_benchmarks_private'/(cid+'.npz');truth=np.load(truthpath)['m']
    else:a,Y,I=load_file(R8/'data_local/fresnel'/(cid+'.exp'),2.)
    for method in METHODS:
        out=HERE/'RESULTS'/cid/method;r=json.loads((out/'RECEIPT.json').read_text());z=np.load(out/'derived.npz');row=dict(case_id=cid,method=method,CPU_seconds=r['CPU_seconds'],train_residual=r['train_relative_complex_residual'],original_Phi=r['original_Phi'],raw_Phi=r['raw_Phi'],bounds=[float(z['m'].min()),float(z['m'].max())],corrector=r['corrector'])
        if austria:
            row.update(NRMSE=float(np.linalg.norm(z['m']-truth)/np.linalg.norm(truth)),raw_NRMSE=float(np.linalg.norm(z['raw']-truth)/np.linalg.norm(truth)),truth_source=str(truthpath),chi_target=3.)
        else:
            mask=z['withheld_mask'];sourcehold=np.zeros(mask.shape,bool);sourcehold[:,z['held_sources']]=True;sourcehold&=np.isfinite(Y);receiverhold=mask&~sourcehold;pred=z['Y_prediction']
            rel=lambda q:float(np.linalg.norm((pred-Y)[q])/np.linalg.norm(Y[q]))
            row.update(heldout_all_residual=rel(mask),heldout_source_residual=rel(sourcehold),heldout_receiver_residual=rel(receiverhold),heldout_channels=int(mask.sum()),calibration_empty_validation_median=float(np.median([x['validation_relative'] for x in r['calibration']])),pixel_GT='UNAVAILABLE')
        rows.append(row)
dump(HERE/'SCORES.json',dict(status='POSTSEAL_OFFLINE_SCORED',rows=rows,new_physical_solves=0,truth_scope='Austria offline scoring only; Fresnel has no pixelGT',heldout_selection=False))
fig,axs=plt.subplots(6,4,figsize=(14,19),constrained_layout=True)
for i,cid in enumerate(CASES):
    side=32 if cid.startswith('austria') else 24
    for j,method in enumerate(METHODS):
        z=np.load(HERE/'RESULTS'/cid/method/'derived.npz');r=next(x for x in rows if x['case_id']==cid and x['method']==method)
        for phase,key in [(0,'raw'),(1,'m')]:
            ax=axs[2*i+phase,j];im=ax.imshow(z[key].reshape(side,side),origin='lower',vmin=0,vmax=4,cmap='magma');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f"{cid} / {method}\n{'auxiliary' if phase==0 else 'common corrector'}",fontsize=9)
            if phase==1:ax.set_xlabel(f"train {r['train_residual']:.3f}; "+(f"NRMSE {r['NRMSE']:.3f}" if 'NRMSE'in r else f"heldout {r['heldout_all_residual']:.3f}"),fontsize=9)
cb=fig.colorbar(im,ax=axs,label='real material contrast chi (display0–4; inversion box0–20)',shrink=.5)
fig.suptitle('External anchors: retain all auxiliary and corrected material images; Fresnel has no pixel GT',fontsize=13)
fig.savefig(HERE/'EXTERNAL_MATERIALS.png',dpi=150)
for ax in axs.flat:
    for img in ax.images:img.set_clim(0,20)
cb.set_label('real material contrast chi; full inversion range0–20')
fig.savefig(HERE/'EXTERNAL_MATERIALS_FULL_RANGE.png',dpi=150);plt.close(fig)
fig,axs=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
for i,cid in enumerate(CASES[:2]):
    a,Y,I=load_file(R8/'data_local/fresnel'/(cid+'.exp'),2.)
    for j,method in enumerate(METHODS):
        z=np.load(HERE/'RESULTS'/cid/method/'derived.npz');mask=z['withheld_mask'];scale=np.linalg.norm(Y[mask])/np.sqrt(mask.sum());p=z['Y_prediction'][mask]/scale;q=Y[mask]/scale;ax=axs[i,j]
        ax.scatter(q.real,p.real,s=3,alpha=.5,label='real');ax.scatter(q.imag,p.imag,s=3,alpha=.5,label='imag');lims=np.array([-3,3]);ax.plot(lims,lims,'k--',lw=.7);ax.set_title(cid+' / '+method,fontsize=9);ax.set_xlabel('heldout measured field / RMS');ax.set_ylabel('predicted / measured RMS');ax.legend(fontsize=7)
fig.savefig(HERE/'FRESNEL_HELDOUT_PREDICTIONS.png',dpi=150);plt.close(fig)
truth=np.load(R8/'data/paper_benchmarks_private/austria_eps4_400MHz.npz')['m'];fig,axs=plt.subplots(1,5,figsize=(15,3.5),constrained_layout=True)
for j,method in enumerate(['truth']+METHODS):
    a=truth if method=='truth' else np.load(HERE/'RESULTS/austria_eps4_400MHz'/method/'derived.npz')['m']
    im=axs[j].imshow(a.reshape(32,32),origin='lower',vmin=0,vmax=4,cmap='magma');axs[j].set_title(method);axs[j].set_xticks([]);axs[j].set_yticks([])
fig.colorbar(im,ax=axs,label='chi=epsilon_r−1; targetchi3',shrink=.75);fig.suptitle('Existing Austria epsilon4/400MHz anchor;32² inverse,150² observation generation')
fig.savefig(HERE/'AUSTRIA_MATERIALS.png',dpi=160);plt.close(fig)
costs=[]
for r in rows:
    x=json.loads((HERE/'RESULTS'/r['case_id']/r['method']/'RECEIPT.json').read_text());costs.append(dict(case_id=r['case_id'],method=r['method'],CPU_seconds=x['CPU_seconds'],physical_counts=x['physical_counts'],full_prediction_counts=x['full_prediction_counts'],auxiliary=x['auxiliary']))
dump(HERE/'COST_LEDGER.json',dict(status='COMPLETE_WITH_LOWER_BOUND_TIMING',rows=costs,solver_CPU_lower_bound=sum(r['CPU_seconds']for r in rows),postseal_report_CPU_near_exit=time.process_time()-cpu,scope='Solver timing includes model/calibration/projector setup, auxiliary steps/audits, full audit-state output, common corrector and full36source Fresnel prediction. Imports/final receipts excluded in per-arm clocks. QUEUE external child CPU overlaps these and is not additive. Existing Austria observation generation is reused, not new free generation. No original measured raw files copied.'))
ledger=json.loads((HERE/'COST_LEDGER.json').read_text());queue=json.loads((HERE/'QUEUE.json').read_text());ledger['external_child_CPU_seconds_sum']=sum(r['external_child_CPU_seconds'] for r in queue['rows']);ledger['external_child_CPU_scope']='all12 solver children including imports and final output writes; overlaps solver_CPU_lower_bound, do not add them';ledger['reused_Austria_generation_receipt']=json.loads((R8/'RESULTS/paper_benchmarks/austria_eps4_400MHz/DATA_RECEIPT.json').read_text());dump(HERE/'COST_LEDGER.json',ledger)
lines=['# External measured Fresnel and Austria anchors','', 'Actual inherited R8 measured source files and32² Austria observations were reused. All12 predictions were sealed before heldout/Austria truth scoring. No original measured raw files were copied into this deliverable.','', 'Fresnel uses original2001 single/twin dielectric measurements at2GHz,24²/0.16m square,27 training sources with39 actual receiver channels each,711 heldout channels, source multipole illumination fitted only to empty-target measurements. Exact source-specific receive/adjoint masks are retained. Empty-field calibration validation medians are about.02515/.02408; unknown absolute units, receiver-horn response and noise covariance remain unresolved. The fixed .15||Ytrain|| ball is explicitly MODEL_DISCREPANCY_DIAGNOSTIC, not a known-noise feasibility certificate. No pixelGT or material superiority claim is available.','', 'Austria reuses actual R8 austria_eps4_400MHz PUBLIC_DATA.npz:150² forward generation,32² inversion,20 plane-wave sources,40 receiver channels,10% fixed-norm complex noise. Material is real chi=epsilon_r−1=3 on the target, not epsilon4 used as chi. The predetermined .12||Y|| ball is a discrepancy diagnostic, not iid CN99. The unchanged R8 equal-area disk-pulse acquisition/Maxwell backend is used; no new backend or full-author-reproduction claim.','', 'All arms use the same real[0,20] box, inherited smoothed-TV/squared-difference prior, and common scaled L-BFGS-B corrector max1000/60CPU. Fresnel starts.1, Austria starts0. Direct uses up to2000/180CPU physical auxiliary iterations. Dual reuses actual R8 mature sourcewise CIE-Dual with bounded300 total PR iterations; Fresnel beta1/.2/.2/.2,mf3/6/10/12, Austria beta8/2/.5,mf5/8/10,currentbeta8. AP/DR use groupedv2,scale6,max5000/audit50/180CPU. These are declared bounded-work comparisons, not equal CPU.','', '| Case | Arm | Train residual | Heldout residual / Austria NRMSE | Phi | Paid CPU s |','|---|---|---:|---:|---:|---:|']
for r in rows:lines.append(f"|{r['case_id']}|{r['method']}|{r['train_residual']:.5f}|{r.get('heldout_all_residual',r.get('NRMSE')):.5f}|{r['original_Phi']:.6g}|{r['CPU_seconds']:.2f}|")
lines += ['', 'All bad auxiliary and corrected images are retained in EXTERNAL_MATERIALS.png and per-arm derived arrays. FRESNEL_HELDOUT_PREDICTIONS.png shows prediction against heldout measurements without fitted gain. Receipt corrector stops/costs distinguish budget caps from numerical convergence. No ball radius, calibration, prior or material bounds were selected by heldout/GT. Physical affine reachability of a data ball alone is not Maxwell-material feasibility.','', 'Source files/API: sibling r8_worktree/research/research_max_struct_r8/data_local/fresnel/{dielTM_dec8f,twodielTM_8f}.exp; sibling RESULTS/paper_benchmarks/austria_eps4_400MHz/PUBLIC_DATA.npz and DATA_RECEIPT.json; R8 src/fresnel_adapter.py, fresnel_solver.py, paper_physics.py, paper_cie_fft.py; current R12 src/field_current.py groupedv2. Context diagnostics read from R8 FRESNEL_REPORT.md and INCIDENT_FIT_AUDIT/REFERENCE_POSE records; no claim to have read all previous reports.']
(HERE/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status='COMPLETE',rows=len(rows),scoring_CPU_seconds=time.process_time()-cpu)))
