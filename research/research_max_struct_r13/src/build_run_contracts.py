"""Metadata-only immutable-receipt wrappers; no truth, physics or candidates.

Only terminal queues yield execution fields. Every expected task is retained;
pending/unknown execution fields stay null. Public G/E/S/m_init headers and
small whitelisted metadata are read, never Y or material-array values.
"""
import os,time
_CPU_START=time.process_time();_WALL_START=time.perf_counter()
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
from zipfile import ZipFile
from datetime import datetime,timezone
import argparse,csv,json
import numpy as np

TERMINAL={'COMPLETED','COMPLETED_WITH_ERRORS','COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY'}
HASH_REASON='userforbidsSHAchecks'
C2={'confirmation_core48','confirmation_mechanism16','physics_v1'}


def read(path):return json.loads(Path(path).read_text())


def encode(v):
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,Path):return str(v)
    raise TypeError(type(v).__name__)


def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,default=encode,allow_nan=False);f.write('\n')


def safe_json(path):
    try:return read(path),None
    except FileNotFoundError:return None,'metadata file absent'
    except Exception as ex:return None,repr(ex)


def npz_headers(path,keys):
    answer={}
    with ZipFile(path) as archive:
        available={x[:-4] for x in archive.namelist() if x.endswith('.npy')}
        for key in keys:
            if key not in available:continue
            with archive.open(key+'.npy') as f:
                version=np.lib.format.read_magic(f)
                if version==(1,0):shape,fortran,dtype=np.lib.format.read_array_header_1_0(f)
                elif version==(2,0):shape,fortran,dtype=np.lib.format.read_array_header_2_0(f)
                else:
                    answer[key]=dict(shape=None,reason='unsupported NPY header version '+str(version));continue
            answer[key]=dict(shape=list(shape),dtype=str(dtype),fortran_order=fortran)
    return answer,sorted(available)


def public_metadata(path):
    headers,keys=npz_headers(path,['G','E0','E','S','S_ports','source_receiver_indices','m_init','centers','coords'])
    whitelist=('grid_shape','components','lower','upper','k','frequency','frequency_Hz','frequency_GHz',
        'frequency_units','wavenumber_units','coordinate_units','length_units','width','h','spacing',
        'radius','cell_radius','inverse_side','forward_side','kind','material_convention','time_convention',
        'noise_sigma_known','sigma_role','noise_model','material_objective')
    scalar={}
    with np.load(path,allow_pickle=False) as z:
        for key in whitelist:
            if key in keys:
                value=z[key]
                if value.size<=64 and value.dtype.kind!='O':scalar[key]=value.item() if value.shape==() else value.tolist()
    e=headers.get('E0',headers.get('E',{})).get('shape');g=headers.get('G',{}).get('shape')
    s=headers.get('S_ports',headers.get('S',{})).get('shape');indices=headers.get('source_receiver_indices',{}).get('shape')
    m=headers.get('m_init',{}).get('shape');components=scalar.get('components')
    if components is None and g and m and len(m)==1 and m[0] and g[0] in (m[0],3*m[0]):components=g[0]//m[0]
    nr=s[1] if s and len(s)==3 else indices[1] if indices and len(indices)==2 else s[0] if s and len(s)==2 else None
    frequency_Hz=scalar.get('frequency_Hz')
    if frequency_Hz is None and scalar.get('frequency_GHz') is not None:frequency_Hz=scalar['frequency_GHz']*1e9
    if frequency_Hz is None and scalar.get('frequency') is not None and scalar.get('frequency_units')=='Hz':frequency_Hz=scalar['frequency']
    return dict(path=str(path),public_keys=keys,array_headers=headers,metadata=scalar,
        material_cells=m[0] if m and len(m)==1 else None,components=components,
        mesh_shape=scalar.get('grid_shape'),sources=e[1] if e and len(e)==2 else None,receivers_per_source=nr,
        frequency_Hz=frequency_Hz,frequency_unknown_reason=None if frequency_Hz is not None else 'No explicitly Hz/GHz-qualified public frequency; not inferred from k/normalized units',
        access_scope='headers and whitelist metadata only; Y/G/E/S/m_init numeric payloads not read')


def provenance(root,campaign,manifest,terminal):
    protocol=manifest.get('protocol',{});snapshot_path=root/'queues'/campaign/'SOURCE_SNAPSHOT.json'
    snapshot,_=safe_json(snapshot_path)
    declared=protocol.get('source_commit') or manifest.get('source_freeze')
    shipping=None;algorithm=declared;evidence=[];qualification='Unknown executed source; current HEAD never substituted'
    if campaign in C2:
        p=root/'audit/CONFIRMATION_FINAL_BUNDLE.json';bundle,_=safe_json(p)
        if bundle:
            shipping=bundle.get('shipping_commit');algorithm=bundle.get('algorithm_freeze',declared)
            evidence.append(str(p));qualification='C2 final shipping record distinguished from algorithm freeze; no byte/hash validation'
    elif campaign=='confirmation_online_controls16':
        shipping='4bd';evidence.append('parent session declaration: C3 controls shipping source4bd')
        qualification='C3 session-declared commit prefix; full commit not invented; separate from inherited algorithm freeze3a91'
    if snapshot:
        algorithm=snapshot.get('source_freeze',algorithm);evidence.append(str(snapshot_path))
        qualification='terminal local source snapshot record; source byte/mtime inventory is not cryptographic verification'
    if not evidence and declared:
        evidence.append('manifest protocol.source_commit / source_freeze declaration')
        qualification='Manifest algorithm/config declaration only; actual executed source unknown'
    executed=shipping or ((snapshot or{}).get('source_freeze'))
    return dict(executed_source_commit=executed if terminal else None,
        shipping_commit=shipping if terminal else None,algorithm_freeze_commit=algorithm,
        declared_config_source_commit=declared,source_evidence=evidence or None,
        qualification=qualification if terminal else 'Execution pending/queue nonterminal; declarations retained but no executed source asserted',
        code_hash=None,config_hash=None,code_hash_reason=HASH_REASON,config_hash_reason=HASH_REASON)


def costs(receipt,job):
    p=receipt.get('physical_cost') or{};lu_key='physical_LU' if 'physical_LU' in p else 'LU' if 'LU' in p else None
    if 'actual_complex_RHS' in p:rhs=p['actual_complex_RHS'];rhs_keys=['actual_complex_RHS']
    else:
        rhs_keys=[k for k in ('source_RHS','adjoint_RHS','tangent_RHS','receiver_transpose_RHS') if k in p]
        rhs=sum(p[k] for k in rhs_keys) if rhs_keys and all(isinstance(p[k],(int,float)) for k in rhs_keys) else None
    external=(job.get('external') or{}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))
    return dict(internal_process_CPU_seconds=receipt.get('whole_process_CPU_seconds'),external_child_CPU_seconds=external,
        wall_seconds=receipt.get('wall_seconds'),physical_LU=p.get(lu_key) if lu_key else None,physical_complex_RHS=rhs,
        physical_cost_source_keys=dict(LU=lu_key,RHS=rhs_keys),physical_cost=p or None,
        projector_cost=receipt.get('projector_cost'),auxiliary_operation_counts=(receipt.get('auxiliary_receipt') or{}).get('operation_counts'),
        cost_scope='Internal/external CPU are alternatives; auxiliary physical deltas are nested, never added twice; unknown cost remains null')


def baseline(algorithm,protocol,aux):
    family='material-only physical optimizer' if algorithm in ('direct','physical_trf','physical_globaldual') else \
        'LSIE complete-current/material joint L-BFGS-B adaptation' if algorithm=='nfft_equations' else \
        'CIE auxiliary current/material alternating' if algorithm in ('cie','dual_cie','global_dual','a_plus','g19','g19_full_apcs','g19_full_current') else \
        'projection/event two-field iterator'
    qualifier='NOT_ORIGINAL_NFFT_SOM' if algorithm=='nfft_equations' else 'global-PR Dual-CIE, not physical GlobalDual' if algorithm=='global_dual' else \
        'physical total-field source-Gram GlobalDual beta0, not feedback beta.2' if algorithm=='physical_globaldual' else \
        'real-material/common-prior project adaptation; no invented author reproduction'
    return dict(name=algorithm,family=family,qualification=qualifier,actual_current_spaces=aux.get('current_spaces'),
        actual_arm=aux.get('arm'),current_solver=aux.get('current_solver'),
        completed_PR_steps=(aux.get('operation_counts') or{}).get('PR_CG_steps'),
        declared_schedule_policy=protocol.get('schedule_policy'),
        missing_stage_settings_reason='Unretained per-stage beta/iteration fields are not reconstructed from current source; see BASELINE_SOURCE_AUDIT.md for source-qualified native settings',
        effective_stage_coefficient_tolerances=aux.get('effective_stage_coefficient_tolerances'),
        stopping_policy=aux.get('stopping_policy'),current_control=aux.get('current_control'),
        sourcewise_cg=aux.get('sourcewise_cg'),initialize_physical_current=aux.get('initialize_physical_current'),
        actual_beta=aux.get('beta'),termination_scope='Native stage setup/reached rank is not proof of full schedule completion',
        physical_audit_CPU_interval=protocol.get('physical_audit_CPU_interval'),audit_every_fallback=protocol.get('audit_every'),
        periodic_audit_policy=aux.get('periodic_physical_audit_policy'))


def make_row(root,manifest_path,manifest,task,queue,terminal,public_cache,validation_rows):
    campaign=manifest['campaign'];protocol=manifest.get('protocol',{});tid=task['task_id'];algorithm=task.get('algorithm',task['method'])
    folder=root/'RESULTS'/campaign/tid;public_path=(root/task['input_file']).resolve()
    row=dict(run_id=campaign+'/'+tid,task_id=tid,case_id=task['case_id'],campaign=campaign,method=task['method'],algorithm=algorithm,
        manifest_path=str(manifest_path),data_id=task['case_id'],public_filepath=str(public_path),
        execution_available=False,execution_status=None,queue_status=queue.get('status') if queue else None,
        provenance=provenance(root,campaign,manifest,terminal),public_metadata=None,material_contract=None,
        box=None,prior=None,initialization=None,state_representation=None,baseline=None,termination=None,
        selection=None,physical_residual=None,auxiliary_defect=None,constraint=None,cost=None,candidate_ids=None,
        online_truth_access=None,validation=None,artifact_coverage=None,unknown_reasons={},
        declared_contract=dict(CPU_cap=task.get('CPU_cap',protocol.get('CPU_cap')),proposal_fraction=protocol.get('proposal_fraction'),
            material_box=protocol.get('material_box'),prior=protocol.get('prior'),initialization=task.get('initialization'),
            selection=protocol.get('selection'),cohort=task.get('cohort'),physics_kind=task.get('physics_kind')))
    row.update(source_commit=row['provenance']['executed_source_commit'],
        shipping_commit=row['provenance']['shipping_commit'],algorithm_freeze_commit=row['provenance']['algorithm_freeze_commit'],
        code_hash=None,config_hash=None,code_hash_reason=HASH_REASON,config_hash_reason=HASH_REASON)
    if not terminal:
        row['unknown_reasons']['execution']='Queue missing/nonterminal; receipt, result metadata and public payload not opened'
        return row
    receipt,error=safe_json(folder/'RECEIPT.json')
    if receipt is None:
        row['unknown_reasons']['execution']='Terminal expected task has no readable frozen receipt: '+str(error)
        row['source_commit']=None;row['provenance']['executed_source_commit']=None
        row['unknown_reasons']['source_commit']='Shipped/config context may be known, but no per-run execution receipt exists'
        return row
    row['execution_available']=True;row['execution_status']=receipt.get('status');row['receipt_path']=str(folder/'RECEIPT.json')
    if public_path not in public_cache:
        try:public_cache[public_path]=(public_metadata(public_path),None)
        except Exception as ex:public_cache[public_path]=(None,repr(ex))
    public,public_error=public_cache[public_path];row['public_metadata']=public
    if public_error:row['unknown_reasons']['public_metadata']=public_error
    pm=public.get('metadata',{}) if public else{};aux=receipt.get('auxiliary_receipt') or{}
    components=public.get('components') if public else receipt.get('components')
    row['material_contract']=dict(material_class='real shared isotropic susceptibility' if components==3 else 'real scalar susceptibility' if components==1 else None,
        convention=pm.get('material_convention','chi=epsilon_r-1; project-declared real susceptibility'),
        time_convention=pm.get('time_convention','exp(-i omega t)'),
        convention_source='public metadata when present; otherwise project THEORY_DERIVATIONS.md opening declaration',
        current_convention='J=chi Q; vector physical electric-current conversion is a separate backend convention',
        components=components)
    row['box']=[pm['lower'],pm['upper']] if 'lower' in pm and 'upper' in pm else protocol.get('material_box_by_family',{}).get(task.get('physics_kind'),protocol.get('material_box'))
    row['prior']=dict(actual_receipt=receipt.get('prior'),declared_protocol=protocol.get('prior_by_family',{}).get(task.get('physics_kind'),protocol.get('prior')))
    row['initialization']=dict(declaration=task.get('initialization'),public_array_header=(public or{}).get('array_headers',{}).get('m_init'),
        numeric_values_read=False,numeric_summary=None,numeric_unknown_reason='Initializer values are not read by metadata-only wrapper',
        auxiliary_current_warm=aux.get('initialize_physical_current'))
    row['state_representation']=dict(material_real_DOF=receipt.get('material_DOF'),
        independent_algorithm_state_DOF=receipt.get('independent_algorithm_state_DOF'),
        actual_per_stage_current_spaces=aux.get('current_spaces'),stored_field_real_scalars=receipt.get('stored_field_real_scalars'),
        legacy_state_real_DOF=receipt.get('state_real_DOF'),
        qualification=receipt.get('state_DOF_qualification') or 'Legacy state_real_DOF, if present, is stored two-field ambient scalars; never an independent CIE/direct variable count')
    row['baseline']=baseline(algorithm,protocol,aux);row['termination']=receipt.get('termination')
    row['selection']=dict(primary=protocol.get('selection'),secondary=protocol.get('secondary'),
        incumbent_scope='original paid physical Phi; auxiliary values do not replace primary selector',prefixes=protocol.get('prefixes'))
    views={}
    for name in ('initial','raw','corrected','primary','secondary'):
        meta,_=safe_json(folder/(name+'.json'))
        views[name]=None if meta is None else {k:meta.get(k) for k in ('candidate','CPU_seconds','phase','relative_residual','residual_over_delta','noise_feasible','legal','eligible','deadline','auxiliary')}
    row['physical_residual']=dict(saved_physical_snapshot_metadata=views,
        PDE_constitutive_propagation_defects=None,scope='Existing saved data-residual metadata only; no Y read, new field solve or residual recomputation')
    projection,_=safe_json(folder/'PROJECTION_HISTORY.json');events,_=safe_json(folder/'EVENTS.json')
    row['auxiliary_defect']=dict(last_saved_projection_metadata=projection[-1] if isinstance(projection,list) and projection else None,
        saved_event_metadata=events if events else None,
        complete_CIE_G19_PR_state_history=None,coverage_reason='Not persisted by main runner; see AUXILIARY_STATE_COVERAGE.md; not reconstructed')
    primary=views['primary'] or{}
    row['constraint']=dict(material_box=row['box'],delta=receipt.get('delta'),sigma_complex_RMS=receipt.get('sigma_complex_RMS'),
        qualification=receipt.get('noise_qualification'),noise_sigma_known=pm.get('noise_sigma_known'),
        primary_legal=primary.get('legal'),primary_discrepancy_feasible=primary.get('noise_feasible'),
        interpretation='Discrepancy feasibility is not material truth or cross-grid/model/calibration certification')
    jobs={x.get('task_id'):x for x in (queue or{}).get('jobs',[]) if isinstance(x,dict)}
    row['cost']=costs(receipt,jobs.get(tid,{}));row['online_truth_access']=receipt.get('truth_access')
    row['candidate_ids']=dict(primary=receipt.get('primary_candidate'),secondary=receipt.get('secondary_candidate'),raw=receipt.get('raw_candidate'),
        corrected=(views['corrected'] or{}).get('candidate'),initial=(views['initial'] or{}).get('candidate'))
    row['artifact_coverage']=dict(physical_snapshots={name:(folder/(name+'.npz')).is_file() for name in views},
        material_candidate_pool=(folder/'CANDIDATES.npz').is_file(),final_auxiliary_state=(folder/'FINAL_AUXILIARY_STATE.npz').is_file(),
        exact_optimizer_continuation_state=None,continuation_unknown_reason='Full PR/L-BFGS/TRF internal state is not generally persisted')
    audit=validation_rows.get(tid)
    if audit:
        row['validation']=dict(source=audit['source'],issues=audit['row'].get('issues'),selector_checks=audit['row'].get('selector_checks'),
            snapshot_defects={k:v.get('defects') for k,v in (audit['row'].get('snapshots') or{}).items() if isinstance(v,dict)})
        row['physical_residual']['PDE_constitutive_propagation_defects']=row['validation']['snapshot_defects']
    else:row['unknown_reasons']['validation']='No existing public acceptance audit row; wrapper does not create one'
    for key in ('box','termination','online_truth_access'):
        if row[key] is None:row['unknown_reasons'][key]='Absent in frozen receipt/declared metadata'
    return row


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve();out=a.output.resolve()
    for protected in ('RESULTS','queues','sealed_truth'):
        if out.is_relative_to(root/protected):raise ValueError('Output cannot modify immutable/private '+protected+' directories')
    if any(part.lower() in ('seal','seals','sealed') for part in out.parts):raise ValueError('Choose a new non-seal audit directory')
    out.mkdir(parents=True,exist_ok=False)
    rows=[];campaigns=[];public_cache={}
    for manifest_path in a.manifest:
        manifest_path=manifest_path.resolve();manifest=read(manifest_path);campaign=manifest['campaign']
        queue,queue_error=safe_json(root/'queues'/campaign/'RECEIPT.json')
        terminal=bool(queue and queue.get('status') in TERMINAL and not queue.get('pending') and not queue.get('running'))
        audit_path=root/'audit'/(campaign+'_PUBLIC_ACCEPTANCE')/'AUDIT.json';validation_rows={}
        if terminal:
            audit,_=safe_json(audit_path)
            if audit:
                validation_rows={r['task_id']:dict(source=str(audit_path),row=r) for r in audit.get('tasks',[]) if 'task_id' in r}
        campaigns.append(dict(campaign=campaign,expected_tasks_from_list=len(manifest['tasks']),
            declared_expected_tasks=manifest.get('expected_tasks'),queue_terminal=terminal,queue_error=queue_error,
            queue_status=(queue or{}).get('status'),source_manifest=str(manifest_path)))
        for task in manifest['tasks']:
            try:row=make_row(root,manifest_path,manifest,task,queue,terminal,public_cache,validation_rows)
            except Exception as ex:
                row=dict(run_id=campaign+'/'+task['task_id'],task_id=task['task_id'],case_id=task.get('case_id'),campaign=campaign,
                    execution_available=False,execution_status=None,unknown_reasons={'contract_build':repr(ex)},
                    code_hash=None,config_hash=None,code_hash_reason=HASH_REASON,config_hash_reason=HASH_REASON)
            rows.append(row)
    report=dict(schema='R13_public_run_contract_wrappers_v1',created_at_utc=datetime.now(timezone.utc).isoformat(),
        root=str(root),truth_access='NONE',expected_run_rows=len(rows),known_terminal_execution_rows=sum(r.get('execution_available',False) for r in rows),
        campaigns=campaigns,runs=rows,code_hash=None,config_hash=None,code_hash_reason=HASH_REASON,config_hash_reason=HASH_REASON,
        metadata_build_cost=dict(CPU_seconds=time.process_time()-_CPU_START,wall_seconds=time.perf_counter()-_WALL_START,
            new_physics=0,new_factorization=0,new_candidate_generation=0,scope='metadata/header reads only; final serialization/exit excluded'),
        scope='Frozen receipts wrapped, never edited. Pending/missing expected rows retained with null execution fields. No GT/Y/material numeric payload or new validation physics.')
    write(out/'RUN_CONTRACTS.json',report)
    columns=list(dict.fromkeys(k for r in rows for k in r))
    with (out/'RUN_CONTRACTS.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for r in rows:
            writer.writerow({k:'null' if r.get(k) is None else json.dumps(r[k],default=encode) if isinstance(r[k],(dict,list)) else r[k] for k in columns})
    print(json.dumps({k:report[k] for k in ('schema','expected_run_rows','known_terminal_execution_rows','truth_access')}),flush=True)


if __name__=='__main__':main()
