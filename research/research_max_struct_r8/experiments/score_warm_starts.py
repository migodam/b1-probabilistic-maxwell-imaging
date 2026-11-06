"""Offline independent warm-start scoring. No Maxwell calls or online GT use.

Seal first, after every frozen task is terminal. Then score from the seal.
The scorer only reuses R7's pure metric definitions, never its solver.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,csv,time,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.score_primary import material_metrics,noise_metrics,fixed_prior

def load(p):return json.loads(Path(p).read_text())
def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def resolve(p):return Path(p) if Path(p).is_absolute() else REPO/p
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def csvwrite(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def seal():
    manifest=ROOT/'configs/warm_all_manifest.json';tasks=load(manifest)['tasks'];entries=[]
    if len(tasks)!=750 or len({t['task_id'] for t in tasks})!=750:
        raise ValueError('Frozen user-amended experiment requires 750 distinct task identities')
    evidence={manifest,ROOT/'configs/warm_primary_manifest.json',ROOT/'configs/solver.json',ROOT/'configs/warm_portfolio_protocol.json',ROOT/'data_public/warm_test/test_manifest.json',ROOT/'audit/generation_test.json',ROOT/'audit/HARDWARE_ASSIGNMENT.json',ROOT/'audit/PRIMARY_AMENDMENT_02.json'}
    if (ROOT/'audit/HARDWARE_REALLOCATION.json').exists():evidence.add(ROOT/'audit/HARDWARE_REALLOCATION.json')
    for name in ('primary_physics','primary_born_network'):
        evidence.add(ROOT/f'data_public/initializers/{name}/RECEIPT.json');evidence.add(ROOT/f'data_public/initializers/{name}/MANIFEST.json')
    for name in ('physics','born'):
        evidence.add(ROOT/f'RESULTS/training/{name}/RECEIPT.json');evidence.add(ROOT/f'RESULTS/training/{name}/weights.pt')
    for name in ('XINAN_FEATURE_INFERENCE.json','MAC_FEATURE_INFERENCE.json'):
        evidence.add(ROOT/'audit'/name)
    for t in tasks:
        d=resolve(t['output_dir']);rec=load(d/'RECEIPT.json')
        if (rec['task_id'],rec['case_id'],rec['arm'])!=(t['task_id'],t['case_id'],t['arm']):
            raise ValueError('Receipt identity mismatch '+t['task_id'])
        if rec['status'] not in ('SEALED_PUBLIC_PREDICTION','FAILED'):raise RuntimeError('Nonterminal frozen task '+t['task_id'])
        names=['RECEIPT.json','PUBLIC_TASK.json']
        if rec['status']=='SEALED_PUBLIC_PREDICTION':names+=['prediction.npz','prediction.json','BUDGET_STATES.npz','initial.npz','initial.json']
        evidence.update(resolve(p) for p in t['public_files'])
        if t.get('initialization_file'):evidence.add(resolve(t['initialization_file']))
        entries.append(dict(task=t,status=rec['status'],files=[stamp(d/n) for n in names]))
    # All two-branch choices are frozen using only paid public objectives.
    protocol=load(ROOT/'configs/warm_portfolio_protocol.json');taskmap={(t['case_id'],t['arm'],t.get('initialization','uniform')):t for t in tasks};portfolios=[]
    for ident in sorted({t['case_id'] for t in tasks}):
        pairs=[]
        for arm in protocol['methods']:
            for init in ('nn','uniform'):pairs.append((arm+' / '+init+' / retained_physics',[(ident,'physical_lbfgs',init),(ident,arm,init)]))
        pairs.append(('analytic_DPCS_physics + uniform_Dual',[(ident,'physical_lbfgs','dpcs'),(ident,'cie_dual','uniform')]))
        for name,keys in pairs:
            candidates=[taskmap[k] for k in keys];records=[load(resolve(t['output_dir'])/'RECEIPT.json') for t in candidates]
            good=[(r['physical_objective'],t['task_id']) for r,t in zip(records,candidates) if r['status']=='SEALED_PUBLIC_PREDICTION']
            portfolios.append(dict(case_id=ident,method=name,candidate_ids=[t['task_id'] for t in candidates],selected_id=min(good)[1] if good else None,selection='Minimum final original Phi, then task_id; before private reads',truth_access='NONE'))
    pp=ROOT/'audit/WARM_PORTFOLIO_PUBLIC_SELECTION.json';dump(pp,dict(records=portfolios,new_physical_calls=0));evidence.add(pp)
    out=ROOT/'audit/WARM_ALL_PUBLIC_SEAL.json'
    if out.exists():raise FileExistsError(out)
    dump(out,dict(status='ALL_FROZEN_WARM_TASKS_SEALED',created_utc=datetime.now(timezone.utc).isoformat(),objects=len({t['case_id'] for t in tasks}),tasks=len(tasks),manifest=stamp(manifest),evidence=[stamp(p) for p in sorted(evidence)],entries=entries,truth_access='NONE',hash_checks='NONE'))
    print(json.dumps(dict(status='SEALED',tasks=len(tasks),objects=len({t['case_id'] for t in tasks}))))

def public_observation(t):
    d={}
    for f in t['public_files']:
        with np.load(resolve(f),allow_pickle=False) as z:
            for k in ('Y','sigma','noise_radius'):
                if k in z.files:d[k]=z[k].copy()
    return d

def initialization_costs():
    gen=load(ROOT/'audit/generation_test.json');feature={r['case_id']:f for r,f in zip(gen['per_object'],gen['feature_cost'])}
    setup={tuple(v['key']):v for v in gen['operator_feature_setup']};inf={}
    for name in ('primary_physics','primary_born_network'):
        rec=load(ROOT/f'data_public/initializers/{name}/RECEIPT.json');inf[name]=(rec,{r['case_id']:r for r in rec['records']})
    return feature,setup,inf

def score():
    begin,cpu=time.perf_counter(),time.process_time();seal=load(ROOT/'audit/WARM_ALL_PUBLIC_SEAL.json')
    if seal['status']!='ALL_FROZEN_WARM_TASKS_SEALED':raise ValueError('A complete public seal is required')
    # The later affine-current mechanism was frozen while all test truth was
    # still private. Finish and seal that extension before the first GT read.
    if (ROOT/'configs/warm_anchor_confirmation_manifest.json').exists():
        anchor=load(ROOT/'audit/WARM_ANCHOR_PUBLIC_SEAL.json')
        if anchor['status']!='ALL_ANCHOR_TASKS_SEALED' or anchor['tasks']!=126:
            raise ValueError('The frozen126 affine-anchor predictions must also be sealed before test scoring')
        expected={t['task_id'] for t in load(ROOT/'configs/warm_anchor_confirmation_manifest.json')['tasks']}
        if len(anchor['entries'])!=126 or {e['task']['task_id'] for e in anchor['entries']}!=expected:
            raise ValueError('Affine-anchor seal does not cover the frozen task identities')
        for s in anchor['evidence']+[f for e in anchor['entries'] for f in e['files']]:
            if stamp(resolve(s['path']))!=s:raise RuntimeError('Changed sealed anchor file '+s['path'])
    for s in seal['evidence']+[f for e in seal['entries'] for f in e['files']]:
        if stamp(resolve(s['path']))!=s:raise RuntimeError('Changed sealed file '+s['path'])
    out=ROOT/'RESULTS/warm_scoring';out.mkdir(exist_ok=False)
    # Private reads begin only after the complete public seal above.
    scoremeta={r['case_id']:r for r in load(ROOT/'sealed_truth/warm_test/score_metadata.json')['records']}
    cases={r['case_id']:r for r in load(ROOT/'data_public/warm_test/test_manifest.json')['records']}
    primary={t['task_id'] for t in load(ROOT/'configs/warm_primary_manifest.json')['tasks']}
    feature,setup,inference=initialization_costs();truths={};publics={};rows=[];budgets=[];witnesses=[]
    hardware_overhead={host:load(ROOT/f'audit/{name}_FEATURE_INFERENCE.json') for host,name in [('Mac','MAC'),('XINAN','XINAN')]}
    hardware_rows={host:{r['case_id']:r for r in rec['records']} for host,rec in hardware_overhead.items()}
    for entry in seal['entries']:
        t=entry['task'];ident=t['case_id'];case=cases[ident];cond=case['condition'];d=resolve(t['output_dir']);rec=load(d/'RECEIPT.json');corr=rec.get('corrector') or {};cost=rec.get('physical_model_cost') or {}
        if ident not in truths:
            with np.load(ROOT/f'sealed_truth/warm_test/{ident}.npz') as z:truths[ident]=(z['m'].ravel().copy(),z['mask'].ravel().astype(bool),z['Y_clean'].copy())
            publics[ident]=public_observation(t)
        truth,support,clean=truths[ident];public=publics[ident];amp=scoremeta[ident]['amplitude'];init=t.get('initialization','uniform');isnn=init.startswith('nn');usesfeature=isnn or init in ('born','dpcs')
        fcpu=feature[ident]['cpu_seconds'] if usesfeature else 0.;fwall=feature[ident]['wall_seconds'] if usesfeature else 0.
        fsetup=setup[(24,cond['sources'],cond['aperture_fraction'])] if usesfeature else {};icpu=iwall=loadcpu=loadwall=0.
        if isnn:
            ir,records=inference['primary_born_network' if init=='nn_born_features' else 'primary_physics'];ir0=records[ident];icpu=ir0['inference_cpu_seconds'];iwall=ir0['inference_wall_seconds'];loadcpu=ir['model_load_cpu_seconds'];loadwall=ir['model_load_wall_seconds']
        actual_prefix=fcpu+icpu;label=t['arm']+' / '+init
        platform=rec.get('execution_platform') or dict(host='Mac',system='Darwin');host='XINAN' if platform.get('system')=='Windows' else 'Mac'
        hr=hardware_rows[host][ident];hrec=hardware_overhead[host];mode='born' if init=='nn_born_features' else 'physics'
        local_feature_cpu=hr['feature']['cpu_seconds'] if usesfeature else 0.;local_feature_wall=hr['feature']['wall_seconds'] if usesfeature else 0.
        local_inference_cpu=hr['models'][mode]['cpu_seconds'] if isnn else 0.;local_inference_wall=hr['models'][mode]['wall_seconds'] if isnn else 0.
        prefix=local_feature_cpu+local_inference_cpu
        cold_extra=(hrec['operator_setup'][hr['operator_key']]['cpu_seconds'] if usesfeature else 0.)+(hrec['model_load'][mode]['cpu_seconds'] if isnn else 0.)+(hrec['imports_cpu_seconds'] if usesfeature else 0.)
        row=dict(task_id=t['task_id'],case_id=ident,cohort=case['cohort'],primary=t['task_id'] in primary,method=label,arm=t['arm'],initialization=init,status=rec['status'],hardware=host,amplitude=amp,sources=cond['sources'],aperture_fraction=cond['aperture_fraction'],relative_complex_RMS_noise=cond['relative_complex_RMS_noise'],generation_grid=cond['generation_grid'],inverse_grid=24,material_image_success=False,material_nrmse=None,otsu_iou=None,original_Phi=rec.get('physical_objective'),relative_residual=rec.get('physical_relative_residual'),solver_cpu_seconds=rec.get('cpu_seconds'),solver_wall_seconds=rec.get('wall_seconds'),feature_cpu_seconds=fcpu,inference_cpu_seconds=icpu,warm_deployment_cpu_seconds=rec.get('cpu_seconds',0)+prefix,cold_deployment_cpu_seconds=rec.get('cpu_seconds',0)+prefix+loadcpu+fsetup.get('cpu_seconds',0),warm_deployment_wall_seconds=rec.get('wall_seconds',0)+fwall+iwall,cold_model_load_cpu_seconds=loadcpu,feature_operator_setup_cpu_seconds=fsetup.get('cpu_seconds',0),LU=cost.get('LU'),complex_RHS=cost.get('actual_complex_RHS'),source_RHS=cost.get('source_RHS'),adjoint_RHS=cost.get('adjoint_RHS'),auxiliary_cpu_seconds=rec.get('auxiliary_cpu_seconds'),auxiliary_iterations=rec.get('auxiliary_iterations'),selected_auxiliary_stage=rec.get('auxiliary_selected_stage'),corrector_stop=corr.get('stop'),corrector_success=corr.get('success'),corrector_iterations=corr.get('iterations'),corrector_evaluations=corr.get('evaluations'),truth_access='OFFLINE_POST_COMPLETE_PUBLIC_SEAL')
        row.update(actual_workflow_cpu_seconds=rec.get('cpu_seconds',0)+actual_prefix,actual_preprocessing_hardware='Mac',deployment_preprocessing_hardware=host,same_hardware_feature_cpu_seconds=local_feature_cpu,same_hardware_inference_cpu_seconds=local_inference_cpu,deployment_prefix_cpu_seconds=prefix,cold_deployment_cpu_seconds=rec.get('cpu_seconds',0)+prefix+cold_extra,warm_deployment_wall_seconds=rec.get('wall_seconds',0)+local_feature_wall+local_inference_wall,hardware_overhead_scope='Independent batched timing of frozen functions; initializers unchanged. Separate-process cold cost includes measured module imports. Training separate.',random_initializer_CPU_not_isolated=(init=='random'))
        if rec['status']=='SEALED_PUBLIC_PREDICTION':
            with np.load(d/'prediction.npz') as z:m=z['m'].ravel();yp=z['Y_prediction'];grad=z['gradient'].ravel() if 'gradient' in z.files else None
            row.update(material_metrics(m,truth,support,amp));row.update(noise_metrics(yp,public,cond['generation_grid']==24))
            row['same_grid_material_noise_success']=bool(row['material_image_success'] and row['noise_CN99_pass']) if cond['generation_grid']==24 else None
            row['independent_saved_Phi']=.5*row['relative_residual']**2+fixed_prior(m);row['Phi_receipt_difference']=row['independent_saved_Phi']-rec['physical_objective']
            row['projected_gradient_inf']=float(np.linalg.norm(m-np.clip(m-grad,.02,20),np.inf)) if grad is not None else None
            with np.load(d/'initial.npz') as z:mi=z['m'].ravel()
            im=material_metrics(mi,truth,support,amp);row['initial_material_nrmse']=im['material_nrmse'];row['initial_otsu_iou']=im['otsu_iou']
            row['initial_original_Phi']=(rec.get('auxiliary_candidate_scores') or [{}])[0].get('objective')
            with np.load(d/'BUDGET_STATES.npz') as z:cs=z['costs'].copy();ms=z['m'].copy()
            for k,(c,mat) in enumerate(zip(cs,ms)):
                met=material_metrics(mat,truth,support,amp)
                budgets.append(dict(task_id=t['task_id'],case_id=ident,cohort=case['cohort'],method=label,hardware=host,index=k,CPU_seconds=float(c[0])+prefix,LU=float(c[1]),complex_RHS=float(c[2]),original_Phi=float(c[3]),material_nrmse=met['material_nrmse'],otsu_iou=met['otsu_iou'],material_image_success=met['material_image_success']))
            if cond['generation_grid']==24:
                gtphi=.5*float(np.linalg.norm(clean-public['Y'])/np.linalg.norm(public['Y']))**2+fixed_prior(truth)
                row['offline_feasible_GT_Phi']=gtphi;row['prediction_minus_feasible_GT_Phi']=rec['physical_objective']-gtphi
                witnesses.append(dict(task_id=t['task_id'],known_better_feasible_objective=bool(rec['physical_objective']>gtphi+1e-12),GT_Phi=gtphi,prediction_Phi=rec['physical_objective'],stop=corr.get('stop'),projected_gradient_inf=row['projected_gradient_inf'],scope='OFFLINE same finite model feasible witness; not uniqueness or basin barrier proof'))
        rows.append(row)
    dump(out/'SCORES.json',dict(rows=rows,objects=42,tasks=len(rows),material_success='R7 inherited amplitude NRMSE<=.15 and Otsu IoU>=.7 and legal box; not an experimental stopping gate',support_semantics='Truth threshold .1+.5*declared amplitude; prediction-only Otsu256. Smooth/nested topology metrics diagnostic.',cost_semantics='Warm deployment includes measured feature application and inference; cold adds feature-operator SVD plus model load. Training/generation separate. Compare absolute CPU within hardware.',new_physical_calls=0))
    csvwrite(out/'SCORES.csv',rows);dump(out/'BUDGET_SCORES.json',dict(rows=budgets,selection='Paid original-Phi incumbent only; no GT iterate choice; sparse checkpoints give lower-resolution budget observation, not interpolation',new_physical_calls=0));dump(out/'FEASIBLE_WITNESSES.json',witnesses)
    dump(out/'RECEIPT.json',dict(status='SCORED_ALL_FROZEN_PREDICTIONS',objects=42,tasks=len(rows),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin,truth_access='OFFLINE_POST_COMPLETE_PUBLIC_SEAL',LU=0,complex_RHS=0,training_not_changed=True))
    print(json.dumps(dict(status='SCORED',tasks=len(rows),output=str(out))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['seal','score']);a=p.parse_args()
    seal() if a.command=='seal' else score()
