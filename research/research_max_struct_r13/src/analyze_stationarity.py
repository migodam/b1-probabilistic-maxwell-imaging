"""Post-joint-authorization original-Phi stationarity diagnostics, source only.

No NumPy/model/material read occurs before joint authorization, ALL required
authorized seals and frozen public-selection metadata have been validated.
This is a charged same-discretization M4 diagnostic, not independent physical
validation, optimization, a new candidate, or a local-minimum certificate.
"""
import os,time
_CPU_START=time.process_time();_WALL_START=time.perf_counter()
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,csv,json,re,sys
from pathlib import Path
from collections import defaultdict

JOINT_CAMPAIGNS={
    'confirmation_core48','confirmation_mechanism16','physics_v1',
    'physics_recovery_v1','confirmation_online_controls16',
    'confirmation_rescue48','confirmation_local_reference16','confirmation_local_trf16'}
CORE_METHODS=('cie','dual_cie','g19','dr','event_integer','direct')
LOCAL_METHODS=('cie','dual_cie','g19','dr','direct','physical_globaldual')
CORE_CASES=tuple('r13_fresh_%04d'%i for i in range(48))
SUBSET_CASES=CORE_CASES[::3]


def read(path):return json.loads(Path(path).read_text())


def relative_path(root,value):
    p=Path(str(value).replace('\\','/'))
    return p if p.is_absolute() else root/p


def authorize(root,path):
    # This function and seal_results are stdlib-only. Do not import a numerical
    # helper here: the full joint guard must finish before NumPy or models.
    try:from .seal_results import validate_seal
    except ImportError:from seal_results import validate_seal
    auth=read(path)
    if auth.get('authorized_joint_offline_scoring') is not True:
        raise PermissionError('Parent joint offline authorization is required')
    if auth.get('root') and Path(auth['root']).resolve()!=root:
        raise ValueError('Joint authorization root differs from requested root')
    required=auth['required_campaigns']
    if not JOINT_CAMPAIGNS.issubset(set(required)):
        raise PermissionError('Joint authorization omits a required registered campaign')
    selected={'confirmation_core48','confirmation_local_reference16','confirmation_local_trf16'}
    if not selected.issubset(set(auth['score_campaigns'])):
        raise PermissionError('Fixed stationarity campaigns are outside the joint score scope')
    seals={}
    for campaign in required:
        location=auth.get('seals',{}).get(campaign)
        seal_path=relative_path(root,location) if location else root/'audit'/(campaign+'_seal.json')
        seal=read(seal_path);validate_seal(seal,require_authorization=True)
        if Path(seal['root']).resolve()!=root or seal['campaign']!=campaign:
            raise ValueError('Required seal root/campaign mismatch: '+campaign)
        seals[campaign]=dict(path=str(seal_path.resolve()),seal=seal)
    for item in auth['public_selection_files']:
        p=relative_path(root,item['relative_path']);actual=p.stat()
        if actual.st_size!=item['bytes'] or actual.st_mtime_ns!=item['mtime_ns']:
            raise RuntimeError('Public selection index changed after joint freeze: '+str(p))
    return auth,seals


def fixed_requests(seals):
    specifications=(
        ('confirmation_core48',CORE_CASES,CORE_METHODS),
        ('confirmation_local_reference16',SUBSET_CASES,LOCAL_METHODS),
        ('confirmation_local_trf16',SUBSET_CASES,('physical_trf',)))
    requests=[]
    for campaign,cases,methods in specifications:
        tasks=seals[campaign]['seal']['manifest_snapshot']['tasks']
        for case in cases:
            for method in methods:
                matches=[t for t in tasks if t.get('case_id')==case and t.get('method')==method]
                requests.append(dict(campaign=campaign,case_id=case,method=method,
                    requested_task_id=case+'__'+method,
                    task=matches[0] if len(matches)==1 else None,
                    scope_issue=None if len(matches)==1 else 'MISSING_REGISTERED_TASK' if not matches else 'AMBIGUOUS_REGISTERED_TASK'))
    if len(requests)!=400:raise RuntimeError('Internal fixed-scope count differs from400')
    return requests


def reject_nn(task):
    text=' '.join(str(task.get(k,'')) for k in ('method','algorithm','initialization','initializer','input_file')).lower()
    return bool(re.search(r'\bnn\b|nn_|_nn|neural|network',text))


def delta(after,before):
    return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before)}


def norm_pair_metrics(np,gd,gr):
    nd,nr=float(np.linalg.norm(gd)),float(np.linalg.norm(gr))
    return dict(data_gradient_L2=nd,prior_gradient_L2=nr,
        data_gradient_inf=float(np.linalg.norm(gd,np.inf)),
        prior_gradient_inf=float(np.linalg.norm(gr,np.inf)),
        data_prior_gradient_cosine=float(np.dot(gd,gr)/(nd*nr)) if nd>0 and nr>0 else None,
        cosine_unknown_reason=None if nd>0 and nr>0 else 'A zero gradient has no defined cosine')


def evaluate_material(np,model,m,Y,lo,hi):
    value,g=model.data_value_gradient(m,Y)  # beta=None: original physical Phi
    state=model.state(m)                  # retained state, no second forward
    prior,gr,_=model.reg.evaluate(m)
    residual=state['Y']-Y
    rr=float(np.linalg.norm(residual));yn=float(np.linalg.norm(Y))
    data=.5*rr*rr/(yn*yn)
    gd=g-gr
    mapping=m-np.clip(m-g,lo,hi)  # exact declared unit step, unscaled gradient
    metrics=dict(Phi=float(value),D=data,R=float(prior),
        Phi_minus_sum_D_R=float(value-data-prior),
        original_data_residual_norm=rr,original_data_relative_residual=rr/yn,
        original_gradient_inf=float(np.linalg.norm(g,np.inf)),
        original_gradient_L2=float(np.linalg.norm(g)),
        unit_step_box_gradient_mapping_inf=float(np.linalg.norm(mapping,np.inf)),
        unit_step_box_gradient_mapping_L2=float(np.linalg.norm(mapping)),
        unit_step=1.,gradient_scale='original unscaled Phi gradient; no optimizer1e6 scale',
        lower_active_exact=int(np.count_nonzero(m==lo)),
        upper_active_exact=int(np.count_nonzero(m==hi)),
        below_lower=int(np.count_nonzero(m<lo)),above_upper=int(np.count_nonzero(m>hi)),
        box_violation_max=float(max(np.max(lo-m,initial=0.),np.max(m-hi,initial=0.))),
        material_in_box_exact=bool(np.all(m>=lo) and np.all(m<=hi)),
        material_min=float(m.min()),material_max=float(m.max()),
        physical_state_residual=float(state['state_relative_residual']),
        **norm_pair_metrics(np,gd,gr))
    if not np.isfinite(g).all() or not all(np.isfinite(metrics[k]) for k in ('Phi','unit_step_box_gradient_mapping_inf')):
        raise FloatingPointError('Nonfinite original objective/gradient diagnostic')
    return metrics


def same_operator(np,a,b):
    return all(a[k].shape==b[k].shape and np.array_equal(a[k],b[k]) for k in ('G','E0','S','Y','grid_shape'))


def load_public(np,path):
    with np.load(path,allow_pickle=False) as z:
        if 'components' in z and int(z['components'])!=1:
            raise ValueError('Fixed scope is scalar TM only')
        result={k:z[k].copy() for k in ('G','S','Y','grid_shape')}
        result['E0']=z['E0' if 'E0' in z else 'E'].copy()
        for k in ('lower','upper'):
            if k in z:result[k]=float(z[k])
    return result


def table(path,rows):
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for row in rows:
            writer.writerow({k:'null' if row.get(k) is None else json.dumps(row[k]) if isinstance(row[k],(dict,list)) else row[k] for k in columns})


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--authorization',type=Path,required=True)
    p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();root=args.root.resolve()
    auth,seals=authorize(root,args.authorization.resolve())
    requests=fixed_requests(seals)
    out=args.out.resolve()
    if any(out.is_relative_to(root/name) for name in ('RESULTS','queues','sealed_truth')):
        raise ValueError('Analysis output must be outside immutable/private data directories')
    for item in seals.values():
        if any(out.is_relative_to(Path(d['path'])) for d in item['seal']['directories']):
            raise ValueError('Analysis output overlaps a frozen prediction directory')
    out.mkdir(parents=True,exist_ok=False)
    # ALL joint/seal/policy guards above precede NumPy, model imports, public
    # numerical arrays and every primary material read.
    import numpy as np
    legacy=args.legacy_repo.resolve();sys.path.insert(0,str(legacy))
    import research.research_fbg_r3_experiments.src.physics as physical_module
    if not Path(physical_module.__file__).resolve().is_relative_to(legacy):
        raise RuntimeError('Imported Maxwell model differs from the requested legacy repository')
    MaxwellModel,Regularizer=physical_module.MaxwellModel,physical_module.Regularizer
    cases=defaultdict(list)
    for request in requests:cases[request['case_id']].append(request)
    rows=[];evaluations=[];setups=[];unique_counter=0
    for case,case_requests in cases.items():
        # Bound array retention to one public case at a time. Exact material
        # duplicates reuse one measured evaluation, including across campaigns.
        public_cache={};operator_groups=[]
        for request in case_requests:
            row={k:v for k,v in request.items() if k!='task'}
            row.update(task_id=None,status='UNAVAILABLE',unique_evaluation_id=None,
                metrics=None,saved_primary_Phi=None,Phi_minus_saved_primary=None,
                cost_reference=None,truth_access='NONE')
            task=request['task'];rows.append(row)
            if task is None:continue
            row['task_id']=task['task_id']
            if reject_nn(task):
                row['status']='OUT_OF_FIXED_SCOPE_NN';continue
            campaign=request['campaign'];manifest=seals[campaign]['seal']['manifest_snapshot']
            folder=root/'RESULTS'/campaign/task['task_id']
            public_path=relative_path(root,task['input_file'])
            row['public_input_path']=str(public_path.resolve())
            row['primary_material_path']=str((folder/'primary.npz').resolve())
            row['source_seal']=seals[campaign]['path']
            try:
                primary_meta=read(folder/'primary.json')
                with np.load(folder/'primary.npz',allow_pickle=False) as z:m=z['m'].copy()
                if np.iscomplexobj(m) or m.ndim!=1 or not np.isfinite(m).all():
                    raise ValueError('Frozen primary material must be a finite real vector')
                row['saved_primary_Phi']=primary_meta.get('Phi')
                row['saved_candidate_id']=primary_meta.get('candidate')
                row['saved_primary_eligible']=primary_meta.get('eligible')
                row['saved_primary_legal']=primary_meta.get('legal')
                if public_path not in public_cache:public_cache[public_path]=load_public(np,public_path)
                public=public_cache[public_path]
                protocol=manifest['protocol'];box=protocol.get('material_box',[.02,20.])
                lo,hi=float(public.get('lower',box[0])),float(public.get('upper',box[1]))
                prior=protocol.get('prior')
                if not isinstance(prior,dict) or not all(k in prior for k in ('smooth','tv','epsilon')):
                    raise ValueError('Exact frozen scalar prior parameters are required')
                prior_key=tuple(float(prior[k]) for k in ('smooth','tv','epsilon'))
                row['material_box']=[lo,hi];row['prior_parameters']=dict(zip(('smooth','tv','epsilon'),prior_key))
                group=next((g for g in operator_groups if g['prior_key']==prior_key and g['box']==(lo,hi) and same_operator(np,g['public'],public)),None)
                if group is None:
                    setup_cpu,setup_wall=time.process_time(),time.perf_counter()
                    model=MaxwellModel(public)
                    model.reg=Regularizer(model.shape,smooth=prior_key[0],tv=prior_key[1],epsilon=prior_key[2])
                    group=dict(public=public,model=model,prior_key=prior_key,box=(lo,hi),materials=[])
                    operator_groups.append(group)
                    setups.append(dict(case_id=case,operator_group=len(operator_groups),
                        model_source=str(Path(physical_module.__file__).resolve()),
                        CPU_seconds=time.process_time()-setup_cpu,wall_seconds=time.perf_counter()-setup_wall,
                        model_setup_cost=dict(model.cost),scope='public model/identical regularizer setup; not added to per-evaluation deltas'))
                old=next((e for e in group['materials'] if np.array_equal(e['m'],m)),None)
                reference=dict(campaign=campaign,task_id=task['task_id'],method=request['method'])
                if old is None:
                    unique_counter+=1;eid=case+'/material_%04d'%unique_counter
                    model=group['model'];before=dict(model.cost);cpu,wall=time.process_time(),time.perf_counter()
                    entry=dict(unique_evaluation_id=eid,case_id=case,status='FAILED',metrics=None,error=None,references=[reference])
                    try:
                        if m.shape!=(model.n,):raise ValueError('Frozen material shape differs from original public model')
                        entry['metrics']=evaluate_material(np,model,m,public['Y'],lo,hi)
                        entry['status']='DIAGNOSTIC_EVALUATED'
                    except Exception as ex:entry['error']=repr(ex)
                    finally:
                        count=delta(model.cost,before)
                        entry['cost']=dict(physical_model_delta=count,LU=count.get('LU'),
                            forward_RHS=count.get('source_RHS'),adjoint_RHS=count.get('adjoint_RHS'),
                            actual_complex_RHS=count.get('actual_complex_RHS'),
                            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                            additional_prior_gradient_evaluations=1 if entry['status']=='DIAGNOSTIC_EVALUATED' else None,
                            cost_scope='charged once for this exact operator/case/material; no optimization/Hessian or candidate generation')
                    evaluations.append(entry);old=dict(m=m.copy(),entry=entry);group['materials'].append(old)
                else:old['entry']['references'].append(reference)
                entry=old['entry'];row['unique_evaluation_id']=entry['unique_evaluation_id']
                row['cost_reference']=entry['unique_evaluation_id'];row['status']=entry['status'];row['metrics']=entry['metrics']
                row['error']=entry['error']
                if row['metrics'] is not None and row['saved_primary_Phi'] is not None:
                    row['Phi_minus_saved_primary']=row['metrics']['Phi']-row['saved_primary_Phi']
                row['difference_scope']='saved-Phi discrepancy is a diagnostic number, never a pass threshold or a new candidate'
            except Exception as ex:row['status']='UNAVAILABLE_OR_INVALID_FROZEN_INPUT';row['error']=repr(ex)
    # NumPy scalar conversion is applied only after authorization.
    def convert(v):
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(type(v).__name__)
    metrics_rows=[]
    for row in rows:
        flat={k:v for k,v in row.items() if k!='metrics'}
        if row['metrics']:flat.update(row['metrics'])
        metrics_rows.append(flat)
    receipt=dict(status='POST_SEAL_SAME_DISCRETIZATION_STATIONARITY_DIAGNOSTIC',
        root=str(root),truth_access='NONE',joint_authorization=str(args.authorization.resolve()),
        validated_required_seals={name:item['path'] for name,item in seals.items()},
        requested_rows=400,retained_rows=len(rows),evaluated_rows=sum(r['status']=='DIAGNOSTIC_EVALUATED' for r in rows),
        unique_material_evaluations=len(evaluations),failed_unique_evaluations=sum(e['status']!='DIAGNOSTIC_EVALUATED' for e in evaluations),
        physical_LU_known_sum=sum(e['cost']['LU'] for e in evaluations if e['cost']['LU'] is not None),
        physical_complex_RHS_known_sum=sum(e['cost']['actual_complex_RHS'] for e in evaluations if e['cost']['actual_complex_RHS'] is not None),
        physical_forward_RHS_known_sum=sum(e['cost']['forward_RHS'] for e in evaluations if e['cost']['forward_RHS'] is not None),
        physical_adjoint_RHS_known_sum=sum(e['cost']['adjoint_RHS'] for e in evaluations if e['cost']['adjoint_RHS'] is not None),
        unknown_unique_evaluation_costs=sum(e['cost']['LU'] is None for e in evaluations),
        CPU_seconds=time.process_time()-_CPU_START,wall_seconds=time.perf_counter()-_WALL_START,
        timer_scope='script entry, all guards/imports/public/primary reads/model setup/unique diagnostics; final output serialization/exit excluded',
        host='parent-declared Apple_M4 analysis',platform=sys.platform,
        physics_scope='same finite public physics, not independent validation',
        model_source=str(Path(physical_module.__file__).resolve()),new_optimization=0,new_Hessian=0,new_candidates=0,
        scope='All fixed400 requested rows retained, including missing/invalid/failing inputs. No GT/case selection by quality; exact duplicates share references and measured cost.',
        interpretation='Unit-step box gradient mappings and gradients are descriptive at frozen primary materials; no stationarity gate, local-minimum/recovery certificate, or added inference information.')
    for name,value in (('STATIONARITY.json',dict(receipt=receipt,rows=rows)),
                       ('UNIQUE_FORWARD_COSTS.json',evaluations),('MODEL_SETUP_COSTS.json',setups),('RECEIPT.json',receipt)):
        with (out/name).open('x') as f:json.dump(value,f,indent=2,default=convert,allow_nan=False);f.write('\n')
    table(out/'STATIONARITY.csv',metrics_rows);table(out/'UNIQUE_FORWARD_COSTS.csv',evaluations)
    print(json.dumps({k:receipt[k] for k in ('status','requested_rows','retained_rows','evaluated_rows','unique_material_evaluations','truth_access')}),flush=True)


if __name__=='__main__':main()
