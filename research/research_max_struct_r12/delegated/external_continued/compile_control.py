"""Postscore native baseline budget control and inherited calibration evidence."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,time,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).absolute().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO));sys.path.insert(0,str(HERE.parent/'external'))
from run_external import CASES,R8,load_file,dump
from research.research_z_alg_r4_diversity.src.physics import Regularizer
paths=[HERE/'RESULTS'/c/'Dual/derived.npz' for c in CASES]
if not all(p.exists() and p.with_name('RECEIPT.json').exists() for p in paths):raise RuntimeError('All3fullallowance predictions required before control scoring')
dump(HERE/'PREDICTION_SEAL.json',dict(status='ALL3_CONTROL_PREDICTIONS_SEALED_BEFORE_CONTROL_SCORE',paths=[str(p.relative_to(HERE))for p in paths],qualification='prior v1scores already known; current runs are explicitly postscore budget diagnostics, selection does not use truth or heldout'))
rows=[]
for cid,p in zip(CASES,paths):
    r=json.loads(p.with_name('RECEIPT.json').read_text());z=np.load(p);row=dict(case_id=cid,method='Dual_full_allowance',Phi=r['original_Phi'],raw_Phi=r['raw_Phi'],train_residual=r['train_relative_complex_residual'],CPU_seconds=r['CPU_seconds'],auxiliary=r['auxiliary'],corrector=r['corrector'])
    if cid.startswith('austria'):
        truth=np.load(R8/'data/paper_benchmarks_private'/(cid+'.npz'))['m'];row.update(NRMSE=float(np.linalg.norm(z['m']-truth)/np.linalg.norm(truth)),raw_NRMSE=float(np.linalg.norm(z['raw']-truth)/np.linalg.norm(truth)))
    else:
        a,Y,I=load_file(R8/'data_local/fresnel'/(cid+'.exp'),2.);mask=z['withheld_mask'];row['heldout_all_residual']=float(np.linalg.norm((z['Y_prediction']-Y)[mask])/np.linalg.norm(Y[mask]))
    rows.append(row)
dump(HERE/'SCORES.json',dict(status='POSTSCORE_BUDGET_DIAGNOSTIC_SCORED',rows=rows,new_physical_solves=0))
diag=[];cid='austria_eps4_400MHz';reg=Regularizer((32,32));truth=np.load(R8/'data/paper_benchmarks_private'/(cid+'.npz'))['m']
for cohort,root,methods in [('v1',HERE.parent/'external',['direct','Dual','AP','DR']),('full_allowance',HERE,['Dual'])]:
    for method in methods:
        out=root/'RESULTS'/cid/method;r=json.loads((out/'RECEIPT.json').read_text());z=np.load(out/'derived.npz')
        for phase,key,phi in [('auxiliary','raw',r['raw_Phi']),('corrected','m',r['original_Phi'])]:
            m=z[key];prior=float(reg.evaluate(m)[0]);diag.append(dict(cohort=cohort,method=method,phase=phase,Phi=phi,prior=prior,data_term=phi-prior,prior_fraction=prior/phi,NRMSE=float(np.linalg.norm(m-truth)/np.linalg.norm(truth)),peak_chi=float(m.max()),near_upper_cells=int(np.sum(m>20-1e-8)),scope='Pure saved-array regularizer/data decomposition; no new forward solves; NRMSE offline only'))
dump(HERE/'AUSTRIA_OBJECTIVE_MATERIAL_DIAGNOSTIC.json',dict(status='OFFLINE_DIAGNOSTIC',rows=diag,new_physical_solves=0))
old=REPO/'research/research_max_struct_r8/RESULTS/fresnel';cal=json.loads((old/'INCIDENT_FIT_AUDIT.json').read_text());ref=json.loads((old/'REFERENCE_VALIDATION.json').read_text());pose=json.loads((old/'REFERENCE_POSE_DIAGNOSTIC.json').read_text());transfer=json.loads((old/'ILLUMINATION_STABILITY_AND_REFERENCE_TRANSFER.json').read_text())
evidence=dict(status='EXACT_INHERITED_R8_EVIDENCE_READOUT_NO_RERUN',noise_covariance='not supplied/unknown, empty calibration residual is not a noise bound',source_files=[str(old/x)for x in ['INCIDENT_FIT_AUDIT.json','REFERENCE_VALIDATION.json','REFERENCE_POSE_DIAGNOSTIC.json','ILLUMINATION_STABILITY_AND_REFERENCE_TRANSFER.json']],cases=[])
for c in cal:
    cid=c['case'];rr=next(x for x in ref['cases']if x['case']==cid);pp=next(x for x in pose['cases']if x['case']==cid);sel=next(x for x in pp['choices']if x['degrees']==pp['selected_degrees'])
    evidence['cases'].append(dict(case_id=cid,empty_target_validation_median=float(np.median([x['validation_relative']for x in c['calibration']])),old_cylindrical_validation_median=float(np.median([x['cylindrical_validation_relative']for x in c['calibration']])),independent_interface_reference_residual=rr['interface_relative_residual'],discretization=rr['discretization'],original_reference_pose='Initial axis interpretation was a retained failure (>100% measured residual); discretization comparison is only VIE vs independent solution at that same nominal geometry',selected_nominal_pose=pp['selected_degrees'],nominal_train_residual_no_gain=sel['train_residual'],diagnostic_gain_heldout_residual=sel.get('selected_pose_heldout_residual'),qualification='Nominal known geometry is not pixelGT; pose/gain diagnostics were not used by inversions'))
evidence['known_reference_transfer']=transfer['known_reference_transfer'];dump(HERE/'FRESNEL_CALIBRATION_FORWARD_EVIDENCE.json',evidence)
fig,axs=plt.subplots(3,4,figsize=(13,10),constrained_layout=True)
for i,cid in enumerate(CASES):
    side=32 if cid.startswith('austria')else 24
    for j,(root,key,title)in enumerate([(HERE.parent/'external','raw','v1aux'),(HERE.parent/'external','m','v1corrected'),(HERE,'raw','fullallowanceaux'),(HERE,'m','fullallowancecorrected')]):
        z=np.load(root/'RESULTS'/cid/'Dual/derived.npz');im=axs[i,j].imshow(z[key].reshape(side,side),origin='lower',vmin=0,vmax=20,cmap='magma');axs[i,j].set_title(cid+'\n'+title,fontsize=9);axs[i,j].set_xticks([]);axs[i,j].set_yticks([])
fig.colorbar(im,ax=axs,label='chi, full box0–20',shrink=.7);fig.savefig(HERE/'FULL_DUAL_CONTROL.png',dpi=160);plt.close(fig)
queue=json.loads((HERE/'QUEUE.json').read_text());dump(HERE/'COST_LEDGER.json',dict(status='COMPLETE',external_child_CPU_seconds=sum(r['external_child_CPU_seconds']for r in queue['rows']),rows=queue['rows'],physical_counts=[dict(case_id=c,counts=json.loads((HERE/'RESULTS'/c/'Dual/RECEIPT.json').read_text())['physical_counts'])for c in CASES],new_forward_solves_in_report=0,prior_v1_CPU_separate=True))
lines=['# Native Dual full-allowance budget control','', 'This is a post-v1-score budget diagnostic requested by root, not a fresh preregistered confirmation. Allv1 endpoints remain immutable. Three new Dual arms use full inherited allowances: Fresnel500/500/500/1500 from actual R8 fresnel_invert.py full stages; Austria300perstage from paper_benchmark.py default --iterations300 and actual case mf5/8/10,beta8/2/.5,currentbeta8. Both use the same180proposal+60corrector CPU caps, realbox0..20, original starts/prior, unchanged DPCS and original physical-Phi selection. Austria literalDPCS15 remains the inherited benchmark contract; R9 noise-adaptive selection is a separate previously studied policy and was not silently substituted.','', '| Case | Train residual | Heldout residual / Austria NRMSE | Raw Austria NRMSE | CPU s |','|---|---:|---:|---:|---:|']
for r in rows:lines.append(f"|{r['case_id']}|{r['train_residual']:.6f}|{r.get('heldout_all_residual',r.get('NRMSE')):.6f}|{r.get('raw_NRMSE','NA')}|{r['CPU_seconds']:.2f}|")
lines += ['', 'AUSTRIA_OBJECTIVE_MATERIAL_DIAGNOSTIC.json decomposes each saved auxiliary/corrected Phi into its exact inherited prior and data term. A physical corrector can lower the registered Phi and increase material error even from a recognizable Dual auxiliary image. This implicates data/prior/selection adequacy in addition to search behavior; it does not prove a universal optimal replacement prior. There is no GT-based early stopping or alternate endpoint selection.','', 'FRESNEL_CALIBRATION_FORWARD_EVIDENCE.json gives exact inherited source paths and numbers: empty-target validation medians.0251475723/.0240806883; equal-area VIE vs independent circular-interface reference errors at24² are.0243024323/.0189421309, dropping at36² to.0129375430/.0120232929. Those forward differences quantify discretization at nominal geometry, not measured-material accuracy. The initial schematic-axis interpretation failed (>100% measurement error) and is preserved; a later training-only pose diagnostic gave nominal measured train errors.1889438838/.3060280655. Optional known-reference gain transfer gave heldout.1607348202/.2961524435 and was not used by inversions. Source/receive horns and absolute field units remain unseparated; covariance unknown. Therefore similar common Fresnel heldout fit is predictive evidence, not physical material recovery.','', 'All bad images and original v1 scores remain preserved. No delta, box, DPCS or calibration retuning occurred. Full-range material comparisons are in FULL_DUAL_CONTROL.png. New3 predictions were sealed before control scores, with explicit acknowledgement that earlier v1 scores were already known.']
(HERE/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status='COMPLETE',rows=rows),default=str))
