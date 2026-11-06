"""Post-six-seal offline historical DEV scoring. Never calls an inverse solver."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import csv,json,sys,time,importlib.util
import numpy as np
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;STUDY=O.parents[1];REPO=STUDY.parents[1]
R9=REPO.parent/'r9_worktree';R8=REPO.parent/'r8_worktree'
CASES=('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3');BASES=('twofold','green_equal_rank','green_strong')
def js(p):return json.loads(p.read_text())
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def stamp(p):
    s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def dump(name,v):(O/name).write_text(json.dumps(v,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def csvwrite(name,rows):
    if not rows:return
    with (O/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def main():
    begin=time.perf_counter();cpu=time.process_time();runs=[];evidence=[]
    # Complete all-six public gate BEFORE this process opens any GT file.
    for case in CASES:
        for basis in BASES:
            d=STUDY/'RESULTS/denominator_dev'/f'{case}__{basis}__galerkin'
            r=js(d/'RECEIPT.json')
            if r['status']!='SEALED_PUBLIC_PREDICTION' or not (d/'prediction.npz').exists():raise RuntimeError('Not all six sealed: '+str(d))
            runs.append((case,'R10_'+basis,d,r));evidence.extend(stamp(d/n) for n in ('RECEIPT.json','prediction.npz','prediction.json','PUBLIC_TASK.json','FROZEN_CONFIG.json','DENOMINATOR_SETUP.json'))
    dump('ALL_SIX_PUBLIC_SEAL.json',dict(status='ALL_SIX_PUBLIC_OUTPUTS_SEALED',truth_opened_by_this_process=False,candidates=[r[1] for r in runs],evidence=evidence))
    metricpath=R8/'research/research_max_struct_r7/experiments/score_primary.py'
    spec=importlib.util.spec_from_file_location('original_r8_material_metrics',metricpath);metrics=importlib.util.module_from_spec(spec);spec.loader.exec_module(metrics)
    truth={};public={};sources=[stamp(metricpath)]
    case_meta={x['case_id']:x for x in js(R8/'research/research_max_struct_r6/data_public/cases.json')['records']}
    for case in CASES:
        tp=R9/'research/research_max_struct_r9/sealed_truth/development'/f'{case}.npz';truth[case]=load(tp);sources.append(stamp(tp))
        pp=R9/'research/research_max_struct_r9/data_public/development'/f'{case}.npz'
        pub=load(R9/'research/research_max_struct_r9/data_public/development/operators.npz');pub.update(load(pp));public[case]=pub;sources.append(stamp(pp))
        for arm in ('cie_dual','global_dual','a_plus'):
            d=R9/'research/research_max_struct_r9/RESULTS/development'/f'{case}__{arm}__fixed_001';r=js(d/'RECEIPT.json')
            if r['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError(d)
            runs.append((case,'R9_'+arm,d,r))
        d=R8/'research/research_max_struct_r6/RESULTS/primary'/f'{case}__cie_single';r=js(d/'RECEIPT.json')
        if r['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError(d)
        runs.append((case,'R6_CIE_single_reference',d,r))
    rows=[];costs=[];pictures={};newcurrents=[]
    for case,method,d,r in runs:
        amplitude=float(case_meta[case]['amplitude']);gt=truth[case];pub=public[case]
        counts=r.get('physical_model_cost',{});ops=r.get('auxiliary_execution_receipt',{}).get('operation_counts',{})
        platform=r.get('execution_platform',{})
        cost=dict(case_id=case,method=method,source=str(d),recorded_device=r.get('device'),recorded_execution_platform=repr(platform),
            executable=r.get('executable'),CPU_seconds=r.get('cpu_seconds'),wall_seconds=r.get('wall_seconds'),
            physical_LU=counts.get('LU'),source_RHS=counts.get('source_RHS'),adjoint_RHS=counts.get('adjoint_RHS'),complex_RHS=counts.get('actual_complex_RHS'),
            auxiliary_CPU_seconds=r.get('auxiliary_cpu_seconds'),corrector_CPU_seconds=r.get('corrector',{}).get('cpu_seconds'),
            denominator_small_LU=ops.get('denominator_small_LU'),denominator_small_triangular_complex_RHS=ops.get('denominator_small_triangular_complex_RHS'),
            CPU_comparability='NOT_A_MATCHED_HOST_TIMING_COMPARISON; retain original receipt hardware/implementation',
            implementation='R6 mature CIE reference; different historical implementation' if method.startswith('R6') else 'R9 inherited FFT/current/corrector' if method.startswith('R9') else 'R10 material metric only; inherited FFT/current/corrector')
        costs.append(cost);sources.append(stamp(d/'RECEIPT.json'))
        paths=[d/'selected_before_corrector.npz',d/'prediction.npz']
        if method.startswith('R10'):paths=sorted(d.glob('stage_*_physical.npz'))+paths
        for p in paths:
            if not p.exists():continue
            pred=load(p);m=pred['m'];sources.append(stamp(p));met=metrics.material_metrics(m,gt['m'],gt['mask'],amplitude)
            state_meta=js(p.with_suffix('.json')) if p.with_suffix('.json').exists() else {}
            yp=pred.get('Y_prediction',pub['S']@pred['j'] if 'j' in pred else None)
            row=dict(case_id=case,method=method,state=p.stem,amplitude=amplitude,truth_status='EXPOSED_OFFLINE_DIAGNOSTIC',
                source=str(p),original_Phi=state_meta.get('objective'),relative_residual=float(np.linalg.norm(yp-pub['Y'])/np.linalg.norm(pub['Y'])) if yp is not None else None,
                material_error_energy=float(np.linalg.norm(m-gt['m'])**2),**met)
            rows.append(row)
            if p.stem=='prediction':pictures[case,method]=m
            if method.startswith('R10') and p.stem in ('prediction','selected_before_corrector'):
                newcurrents.append((case,method,p.stem,m,pred.get('j')))
    csvwrite('SCORES_ALL_SAVED_STATES.csv',rows);csvwrite('SCORES_FINAL.csv',[r for r in rows if r['state']=='prediction']);csvwrite('COST_LU_RHS.csv',costs)
    # New candidates only; existing 449-state diagnostics are not touched or recomputed.
    sys.path.insert(0,str(REPO));from research.research_max_struct_r10.src.energy_diagnostics import Geometry,physical_current,energy
    G,S=public[CASES[0]]['G'],public[CASES[0]]['S'];geom=Geometry(G,S);single=[];two=[];rem=[];cm=[];LU=0;RHS=0;jtruth={}
    for case,method,state,m,j in newcurrents:
        pub=public[case]
        if not np.array_equal(G,pub['G']) or not np.array_equal(S,pub['S']):raise RuntimeError('Geometry mismatch')
        if case not in jtruth:jtruth[case]=physical_current(G,pub['E0'],truth[case]['m']);LU+=1;RHS+=pub['E0'].shape[1]
        if j is None:j=physical_current(G,pub['E0'],m);LU+=1;RHS+=pub['E0'].shape[1]
        dm=m-truth[case]['m'];dj=j-jtruth[case]
        base=dict(case_id=case,method=method,state=state)
        cm.append(dict(**base,material_error_energy=energy(dm),true_current_error_energy=energy(dj),receiver_true_error_energy=energy(S@dj),internal_error_energy=energy(G@dj)))
        for quantity,x,family in [('material',dm,'complex'),('material_real',dm,'real'),('true_current',dj,'complex')]:
            a,b,c=geom.decompose(x,family)
            for rr,target in [(a,single),(b,two),(c,rem)]:target.extend(dict(**base,quantity=quantity,**v) for v in rr)
    csvwrite('NEW_CANDIDATE_CURRENT_ENERGY.csv',cm);csvwrite('NEW_CANDIDATE_GREEN_ENERGY.csv',single);csvwrite('NEW_CANDIDATE_TWOFOLD_ENERGY.csv',two);csvwrite('NEW_CANDIDATE_CROSS_REMAINDER.csv',rem)
    methods=['R6_CIE_single_reference','R9_cie_dual','R9_global_dual','R9_a_plus']+['R10_'+b for b in BASES]
    fig,axs=plt.subplots(2,len(methods)+1,figsize=(19,5),constrained_layout=True)
    for i,case in enumerate(CASES):
        for k,method in enumerate(['GT']+methods):
            ax=axs[i,k];m=truth[case]['m'] if method=='GT' else pictures[case,method]
            im=ax.imshow(m.reshape(24,24),vmin=.02,vmax=20,cmap='viridis',origin='lower');ax.set_xticks([]);ax.set_yticks([])
            if method=='GT':label=case[4:8]+' GT'
            else:
                x=next(x for x in rows if x['case_id']==case and x['method']==method and x['state']=='prediction')
                label=method.replace('_reference','').replace('green_equal_rank','G28').replace('green_strong','G19')+'\nNRMSE %.4f | IoU %.3f'%(x['material_nrmse'],x['otsu_iou'])
            ax.set_title(label,fontsize=8)
    fig.colorbar(im,ax=axs,shrink=.78,label='physical material m; shared scale [.02,20]')
    fig.suptitle('Disclosed historical DEV | post-seal offline scoring | no unseen-data claim',fontsize=12)
    fig.savefig(O/'MATERIAL_COMPARISON.png',dpi=180);fig.savefig(O/'MATERIAL_COMPARISON.pdf');plt.close(fig)
    dump('SOURCE_PROVENANCE.json',sources)
    dump('RECEIPT.json',dict(status='COMPLETE',new_candidates=6,historical_references=8,all_six_seal=str(O/'ALL_SIX_PUBLIC_SEAL.json'),
        truth_access='EXPOSED_OFFLINE_DIAGNOSTIC_AFTER_ALL_SIX_SEALED',GT_selection=False,new_inverse_runs=0,scoring_CPU_seconds=time.process_time()-cpu,
        scoring_wall_seconds=time.perf_counter()-begin,new_offline_physical_LU=LU,new_offline_complex_RHS=RHS,
        metric_source=str(metricpath),normalization='original material_metrics: RMSE/declared amplitude from R6 public cases.json',
        primary_CIE_reference='R6 mature CIE-single old implementation; no exact R9 historical same-code fixed001 CIE-single available',
        scope='Historical exposed DEV only; material metric candidate. Current algorithm, affine FFT family, acquisition and corrector inherited. Twofold vs G28 near-equivalence cannot establish receiver-specific information gain.',
        timing='Receipt CPU values retained with platform/device; no Mac/Windows/unknown-host speedup claim; nested costs overlap'))
    print(json.dumps(dict(status='COMPLETE',new_candidates=6,rows=len(rows),output=str(O))),flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
