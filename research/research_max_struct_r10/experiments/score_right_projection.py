"""Root-authorized right3 scoring. Exposed posthoc structural diagnostic only."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,csv,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics,dump,csvwrite,stamp
CASES=('r10_twod_0001','r10_twod_0004','r10_twod_0022')
METHODS=('cie_single','a_plus','matrix_green19','matrix_green19_right')
LABELS=dict(cie_single='Frozen CIE single',a_plus='Frozen A+ global',matrix_green19='Frozen PGP19',matrix_green19_right='Posthoc GP19 (right)')
def js(p):return json.loads(p.read_text())
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--authorization-note',required=True)
    ap.add_argument('--queue',type=Path,default=ROOT/'delegated/twod/fetched_right_queue/RECEIPT.json')
    ap.add_argument('--output',type=Path,default=ROOT/'delegated/right_projection_score');args=ap.parse_args()
    start=time.perf_counter();cpu=time.process_time();q=js(args.queue)
    if q['status']!='COMPLETED' or len(q['jobs'])!=3 or q.get('running') or q.get('pending',0):raise RuntimeError('Right3 queue not complete')
    jobs={r['case_id']:r for r in q['jobs']};new={};evidence=[stamp(args.queue)]
    for case in CASES:
        d=ROOT/'RESULTS/right_projection_diagnostic'/f'{case}__matrix_green19_right__noise_ball_99';r=js(d/'RECEIPT.json');setup=js(d/'DENOMINATOR_SETUP.json')
        if r['status']!='SEALED_PUBLIC_PREDICTION' or not (d/'prediction.npz').exists() or jobs[case]['status']!='FINISHED':raise RuntimeError('Right3 missing public seal')
        if (setup['projection'],setup['actual_denominator_rank'],setup['matrix_inner_maxiter'])!=('right',19,8):raise RuntimeError('Wrong frozen right/rank/inner contract')
        task=js(d/'PUBLIC_TASK.json')
        if task['material_readout']['label']!='matrix_green19_right' or task['case_id']!=case:raise RuntimeError('Wrong candidate identity')
        new[case]=(d,r,setup);evidence.extend(stamp(d/n) for n in ('RECEIPT.json','prediction.npz','PUBLIC_TASK.json','DENOMINATOR_SETUP.json','NUMERICAL_SCHEDULE.json'))
    if args.output.exists():raise FileExistsError('Write-once right3 score directory exists')
    out=args.output;out.mkdir(parents=True)
    dump(out/'ALL3_PUBLIC_SEAL.json',dict(status='ALL3_RIGHT_OUTPUTS_SEALED',scope='EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC',root_authorization=args.authorization_note,evidence=evidence))
    # Only previously authorized/unlocked original24 GT cache. Conditional4 is never referenced.
    reference=ROOT/'delegated/twod_score/scored';refs=js(reference/'SCORES_FINAL.json');ix={(r['case_id'],r['method']):r for r in refs}
    with (reference/'SCORES_STAGES.csv').open() as f:stages=list(csv.DictReader(f))
    stageix={(r['case_id'],r['method'],r['state']):r for r in stages};rows=[];costs=[];pics={};truth={}
    for case in CASES:
        cached=load(ROOT/'delegated/twod_energy/gt_current_cache'/(case+'.npz'));truth[case]=cached['m_gt'];jt=cached['j_gt'];meta=ix[case,'cie_single'];amp=meta['amplitude'];high=meta['high_chi']
        with np.load(ROOT/'data_public/twod'/(case+'.npz'),allow_pickle=False) as z:Y=z['Y'].copy();sigma=float(z['sigma'])
        for method in METHODS:
            right=method=='matrix_green19_right'
            if right:d,r,setup=new[case];job=jobs[case]
            else:d=ROOT/'RESULTS/twod'/ix[case,method]['task_id'];r=js(d/'RECEIPT.json');setup={};job={}
            c=r.get('physical_model_cost') or {};ops=(r.get('auxiliary_execution_receipt') or {}).get('operation_counts') or {}
            ext=job.get('child_process_CPU_seconds') if right else ix[case,method]['external_full_child_CPU_seconds']
            costs.append(dict(case_id=case,method=method,cost_role='NEW_POSTHOC_INVERSE_ALREADY_RUN' if right else'REUSED_FROZEN_REFERENCE_NOT_NEW_SPENDING',
                external_full_child_CPU_seconds=ext,external_CPU_scope=job.get('CPU_scope') if right else ix[case,method]['external_CPU_scope'],
                nested_R10_CPU_seconds=r.get('r10_total_cpu_seconds'),nested_R8_CPU_seconds=r.get('cpu_seconds'),nested_auxiliary_CPU_seconds=r.get('auxiliary_cpu_seconds'),
                nested_corrector_CPU_seconds=(r.get('corrector') or {}).get('cpu_seconds'),physical_LU=c.get('LU'),physical_LU_dimension=576,
                physical_source_RHS=c.get('source_RHS'),physical_adjoint_RHS=c.get('adjoint_RHS'),physical_complex_RHS=c.get('actual_complex_RHS'),
                denominator_small_LU=ops.get('denominator_small_LU'),denominator_small_LU_dimension_sum=ops.get('denominator_small_LU_dimension_sum'),
                denominator_small_triangular_complex_RHS=ops.get('denominator_small_triangular_complex_RHS'),corrector_stop=(r.get('corrector') or {}).get('stop'),
                corrector_success=(r.get('corrector') or {}).get('success'),corrector_iterations=(r.get('corrector') or {}).get('iterations'),source=str(d/'RECEIPT.json')))
            for p in sorted(d.glob('stage_*_physical.npz'))+[d/'selected_before_corrector.npz',d/'prediction.npz']:
                if not p.exists():continue
                z=load(p);m=z['m'];sm=js(p.with_suffix('.json'));state=p.stem
                if right:metrics=material_metrics(m,truth[case],amp,high);metric_source='NEW_RIGHT3_POSTSEAL_SCORE'
                else:
                    src=ix[case,method] if state=='prediction' else stageix[case,method,state]
                    metrics={k:(float(src[k]) if src[k] not in ('',None) and k not in ('material_success','material_legal') else src[k]) for k in ('amplitude_NRMSE','fixed_support_IoU','legacy_Otsu_IoU','Otsu_full_support_IoU','material_success','material_legal')}
                    for k in ('material_success','material_legal'):
                        if isinstance(metrics[k],str):metrics[k]=metrics[k]=='True'
                    metric_source='REUSED_ORIGINAL168_SCORE_UNCHANGED'
                yp=z['Y_prediction'];rr=float(np.linalg.norm(yp-Y)/np.linalg.norm(Y));we=float(np.linalg.norm(yp-Y)/sigma)
                pg=float(np.max(abs(m-np.clip(m-z['gradient'],.02,20)))) if 'gradient' in z else None
                row=dict(case_id=case,method=method,state=state,scope='EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC',metric_source=metric_source,
                    auxiliary_selected_stage=r.get('auxiliary_selected_stage'),original_Phi=sm.get('objective'),relative_physical_residual=rr,
                    whitened_data_residual=we,CN99_radius=meta['CN99_radius'],CN99_raw_noise_pass=bool(we<=meta['CN99_radius']),
                    saved_projected_gradient_inf=pg,corrector_stop=(r.get('corrector') or {}).get('stop'),corrector_success=(r.get('corrector') or {}).get('success'),
                    corrector_iterations=(r.get('corrector') or {}).get('iterations'),external_full_child_CPU_seconds=ext,**metrics)
                if 'j' in z:row['physical_current_error_squared_energy']=float(np.linalg.norm(z['j']-jt)**2)
                rows.append(row)
                if state=='prediction':pics[case,method]=m
    csvwrite(out/'ALL_STAGE_START_FINAL_SCORES.csv',rows);csvwrite(out/'FINAL_COMPARISON.csv',[r for r in rows if r['state']=='prediction']);csvwrite(out/'COST_LU_RHS.csv',costs)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(3,5,figsize=(12,7.8),constrained_layout=True)
    for i,case in enumerate(CASES):
        for k,method in enumerate(('GT',)+METHODS):
            ax=axs[i,k];m=truth[case] if method=='GT' else pics[case,method];im=ax.imshow(m.reshape(24,24),origin='lower',cmap='viridis',vmin=.02,vmax=20)
            ax.set_xticks([]);ax.set_yticks([])
            title=case.removeprefix('r10_twod_')+' GT' if method=='GT' else LABELS[method]
            if method!='GT':x=next(r for r in rows if r['case_id']==case and r['method']==method and r['state']=='prediction');title+='\nN %.4f | IoU %.3f'%(x['amplitude_NRMSE'],x['fixed_support_IoU'])
            ax.set_title(title,fontsize=8)
    fig.colorbar(im,ax=axs,shrink=.8,label='real contrast chi; shared [.02,20]');fig.suptitle('EXPOSED POSTHOC STRUCTURAL DIAGNOSTIC | right GP19 vs frozen PGP19 | not new confirmation',fontsize=11)
    fig.savefig(out/'MATERIAL_COMPARISON.png',dpi=180);fig.savefig(out/'MATERIAL_COMPARISON.pdf');plt.close(fig)
    own=[c for c in costs if c['cost_role']=='NEW_POSTHOC_INVERSE_ALREADY_RUN']
    dump(out/'RECEIPT.json',dict(status='COMPLETE',scope='EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC',root_authorization=args.authorization_note,
        right_predictions_scored=3,original168_rescored=False,original168_results_modified=False,new_inverse_runs=0,new_Maxwell_forward_calls=0,new_LU=0,
        GT_source='Original24 previously authorized GT-current cache only; conditional4 truth not read',
        new_right_inverse_full_child_CPU_seconds=sum(c['external_full_child_CPU_seconds'] for c in own),new_right_inverse_physical_LU=sum(c['physical_LU'] for c in own),
        new_right_inverse_physical_complex_RHS=sum(c['physical_complex_RHS'] for c in own),cost_semantics='These three inversions were already run by root/producer; scorer adds none. Reference inverse costs reused, not charged as new; internal timers nested.',
        scoring_CPU_seconds=time.process_time()-cpu,scoring_wall_seconds=time.perf_counter()-start,source_provenance=evidence))
    print(json.dumps([r for r in rows if r['method']=='matrix_green19_right' and r['state']=='prediction'],indent=2))
if __name__=='__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):main()
