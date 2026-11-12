"""Public-only postqueue prediction/selector acceptance; never loads truth.

Only saved states, public arrays, JSON metadata and matrix products are used.
No Maxwell model import, new factorization, optimization, or numerical defect
pass threshold. Missing tasks stay in the full expected-task denominator.
"""
import os,time
_CPU_START=time.process_time();_WALL_START=time.perf_counter()
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,json,re
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from zipfile import ZipFile

TERMINAL={'COMPLETED','COMPLETED_WITH_ERRORS','COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY'}
VIEWS=('initial','raw','corrected','primary','secondary')


def read(path):return json.loads(Path(path).read_text())


def encode(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,Path):return str(value)
    raise TypeError(type(value).__name__)


def write(path,value):path.write_text(json.dumps(value,indent=2,default=encode,allow_nan=False)+'\n')


def issue(rows,kind,**details):rows.append(dict(kind=kind,**details))


def maxdiff(a,b):
    if a.shape!=b.shape:return dict(same_shape=False,shape_a=list(a.shape),shape_b=list(b.shape),exact_equal=False,max_abs_difference=None)
    finite=bool(np.isfinite(a).all() and np.isfinite(b).all())
    return dict(same_shape=True,exact_equal=bool(np.array_equal(a,b)),finite=finite,
        max_abs_difference=float(np.max(abs(a-b),initial=0.)) if finite else None)


def norm_record(value,denominator,units,denominator_name):
    absolute=float(np.linalg.norm(value));den=float(denominator)
    return dict(Frobenius_norm=absolute if np.isfinite(absolute) else None,
        relative_norm=absolute/den if np.isfinite(absolute) and np.isfinite(den) and den>0 else None,
        denominator_norm=den if np.isfinite(den) else None,denominator=denominator_name,units=units,
        interpretation='numerical defect only; no pass/error/recovery threshold')


def initial_family(task):
    text=' '.join(str(task.get(k,'')) for k in ('initialization','initializer','initial_state','input_file')).lower()
    return 'NN' if re.search(r'\bnn\b|nn_|_nn|neural|network',text) else 'public_non_NN'


def npz_preflight(path,keys,max_bytes):
    with ZipFile(path) as archive:
        total=sum(archive.getinfo(k+'.npy').file_size for k in keys)
    if total>max_bytes:raise MemoryError('Requested NPZ arrays exceed declared memory cap before loading')


def load_public(path,max_bytes):
    # Deliberately load only physical operators/initial material/bounds. No
    # observation, truth-spec path, hidden material, training label or GT key.
    with np.load(path,allow_pickle=False) as z:
        names=set(z.files)
        keys=[k for k in ('G','m_init','components','lower','upper','source_receiver_indices','grid_shape') if k in names]
        ekey='E0' if 'E0' in names else 'E';skey='S_ports' if 'S_ports' in names else 'S'
        npz_preflight(path,keys+[ekey,skey],max_bytes)
        p={k:z[k].copy() for k in keys}
        p['E']=z[ekey].copy();p['S']=z[skey].copy()
    for k in ('G','E','S','m_init'):
        if k not in p:raise ValueError('Missing public '+k)
    return p


def same_start(a,b):
    return all(k in a and k in b and np.array_equal(a[k],b[k]) for k in ('G','E','S','m_init')) and \
        all(np.array_equal(a.get(k,np.array(1 if k=='components' else [])),b.get(k,np.array(1 if k=='components' else []))) for k in ('components','source_receiver_indices'))


def canonical_cost(cost):
    lu_key='physical_LU' if 'physical_LU' in cost else 'LU' if 'LU' in cost else None
    if 'actual_complex_RHS' in cost:rhs=cost['actual_complex_RHS'];keys=['actual_complex_RHS']
    else:
        keys=[k for k in ('source_RHS','adjoint_RHS','tangent_RHS','receiver_transpose_RHS') if k in cost]
        rhs=sum(cost[k] for k in keys) if keys and all(isinstance(cost[k],(int,float)) for k in keys) else None
    return dict(physical_LU=cost.get(lu_key) if lu_key else None,physical_complex_RHS=rhs,
        source_keys=dict(LU=lu_key,RHS=keys),unknown_reason=None if lu_key and rhs is not None else 'counter absent; never imputed as zero')


def state_defects(state,p,actions):
    m,J,Q,yp=(state[k] for k in ('m','j','Q','Y_prediction'))
    components=int(p.get('components',1));d=np.repeat(m.ravel(),components)
    if J.shape!=Q.shape or J.shape!=p['E'].shape or d.shape!=(J.shape[0],):
        raise ValueError('Saved fields/material/components do not match public E grid')
    GJ=p['G']@J;actions['Green_forward_calls']+=1;actions['Green_forward_complex_RHS']+=J.shape[1]
    S=p['S']
    if S.ndim==3:
        if S.shape[0]!=J.shape[1]:raise ValueError('Source-specific port count mismatch')
        receive=np.column_stack([S[s]@J[:,s] for s in range(J.shape[1])])
        actions['receiver_forward_calls']+=J.shape[1];actions['receiver_forward_complex_RHS']+=J.shape[1]
    elif 'source_receiver_indices' in p:
        indices=p['source_receiver_indices']
        if len(indices)!=J.shape[1]:raise ValueError('Source-specific receiver-index count mismatch')
        receive=np.column_stack([S[np.asarray(indices[s],int)]@J[:,s] for s in range(J.shape[1])])
        actions['receiver_forward_calls']+=J.shape[1];actions['receiver_forward_complex_RHS']+=J.shape[1]
    else:
        receive=S@J;actions['receiver_forward_calls']+=1;actions['receiver_forward_complex_RHS']+=J.shape[1]
    if yp.shape!=receive.shape:raise ValueError('Saved prediction shape differs from actual public receiver map')
    return dict(constitutive=norm_record(J-d[:,None]*Q,np.linalg.norm(J),'current units','||saved J||F'),
        propagation=norm_record(Q-p['E']-GJ,max(np.linalg.norm(Q),np.linalg.norm(p['E']),np.linalg.norm(GJ)),
            'field units','max(||saved Q||F,||E||F,||GJ||F)'),
        reception=norm_record(yp-receive,max(np.linalg.norm(yp),np.linalg.norm(receive)),
            'receiver field units','max(||saved Y_prediction||F,||public receive(J)||F)'),
        components=components,material_cells=m.size,source_columns=J.shape[1],
        receiver_map='per-source explicit matrices' if S.ndim==3 else 'common S with per-source original receiver indices' if 'source_receiver_indices' in p else 'common S',
        numeric_finiteness=bool(all(np.isfinite(state[k]).all() for k in state)))


def audit_task(root,task,protocol,p,job,actions,max_bytes):
    started,wall=time.process_time(),time.perf_counter();issues=[]
    folder=root/'RESULTS'/task['campaign']/task['task_id'];receipt={};hist=[];mats=None
    row=dict(task_id=task['task_id'],case_id=task['case_id'],method=task['method'],algorithm=task.get('algorithm',task['method']),
        initialization_family=initial_family(task),issues=issues,snapshots={},selector_checks={},truth_access='NONE')
    for filename in ('RECEIPT.json','HISTORY.json','CANDIDATES.npz','PREFIXES.json'):
        if not(folder/filename).is_file():issue(issues,'MISSING_ARTIFACT',file=filename)
    try:
        receipt=read(folder/'RECEIPT.json');hist=read(folder/'HISTORY.json')
        npz_preflight(folder/'CANDIDATES.npz',['m'],max_bytes)
        with np.load(folder/'CANDIDATES.npz',allow_pickle=False) as z:mats=z['m'].copy()
        if mats.nbytes>max_bytes:raise MemoryError('Candidate material array exceeds declared audit memory cap')
        if len(hist)!=len(mats):issue(issues,'HISTORY_CANDIDATE_LENGTH_MISMATCH',history=len(hist),materials=len(mats))
        if receipt.get('truth_access')!='NONE':issue(issues,'ONLINE_TRUTH_DECLARATION_NOT_NONE',value=receipt.get('truth_access'))
    except Exception as ex:issue(issues,'CANDIDATE_LOAD_FAILED',error=repr(ex));mats=None
    cap=float(task.get('CPU_cap',protocol['CPU_cap']));bounds=protocol.get('material_box',[.02,20.])
    lower=float(p.get('lower',bounds[0])) if p is not None else float(bounds[0])
    upper=float(p.get('upper',bounds[1])) if p is not None else float(bounds[1])
    valid=[];all_indices={};last_time=None
    for position,(h,m) in enumerate(zip(hist,mats if mats is not None else [])):
        index=h.get('candidate')
        if isinstance(index,bool) or not isinstance(index,int):
            issue(issues,'INVALID_CANDIDATE_INDEX_TYPE',position=position);continue
        all_indices[index]=position
        if index!=position:issue(issues,'NONCONTIGUOUS_CANDIDATE_INDEX',position=position,stored=index)
        valid_index=not isinstance(index,bool) and isinstance(index,int) and index==position
        if isinstance(index,bool) or not isinstance(index,int):issue(issues,'INVALID_CANDIDATE_INDEX_TYPE',position=position)
        if np.iscomplexobj(m) or not np.isfinite(m).all():issue(issues,'INVALID_MATERIAL_ARRAY',candidate=index);continue
        legal=bool(np.min(m)>=lower-1e-12 and np.max(m)<=upper+1e-12)
        if h.get('legal') is not legal:issue(issues,'LEGAL_FLAG_MISMATCH',candidate=index,computed=legal,stored=h.get('legal'))
        t,deadline=h.get('CPU_seconds'),h.get('deadline')
        if t is None or deadline is None:issue(issues,'MISSING_PAID_TIME_OR_DEADLINE',candidate=index);continue
        if last_time is not None and t<last_time:issue(issues,'NONMONOTONE_CANDIDATE_CPU',candidate=index)
        last_time=t;eligible=bool(t<=deadline and t<=cap)
        if h.get('eligible') is not eligible:issue(issues,'ELIGIBILITY_FLAG_MISMATCH',candidate=index,computed=eligible,stored=h.get('eligible'))
        if not all(isinstance(h.get(k),(int,float)) and np.isfinite(h[k]) for k in ('Phi','R')):
            issue(issues,'MISSING_OR_NONFINITE_SELECTOR_VALUE',candidate=index);continue
        if valid_index and legal and eligible:valid.append(h)
    def minimum(pool,key):return min(pool,key=key)['candidate'] if pool else None
    expected_primary=minimum(valid,lambda h:h['Phi'])
    expected_secondary=minimum([h for h in valid if h.get('noise_feasible') is True],lambda h:(h['R'],h['Phi']))
    direct=task.get('algorithm',task['method'])=='direct'
    expected_raw=minimum([h for h in valid if h.get('phase')=='initial' or not direct and h.get('phase')=='proposal'],lambda h:h['Phi'])
    expected_corrected=expected_primary if direct else minimum([h for h in valid if h.get('phase')=='corrector'],lambda h:h['Phi'])
    if expected_corrected is None:expected_corrected=expected_primary
    for name,expected in (('primary',expected_primary),('secondary',expected_secondary)):
        stored=receipt.get(name+'_candidate');row['selector_checks'][name]=dict(stored=stored,expected=expected,exact_match=stored==expected)
        if stored!=expected:issue(issues,'SELECTOR_INDEX_MISMATCH',view=name,stored=stored,expected=expected)
    row['selector_checks']['raw']=dict(stored=receipt.get('raw_candidate'),expected=expected_raw)
    if receipt.get('raw_candidate')!=expected_raw and receipt.get('status')=='COMPLETE':
        issue(issues,'RAW_SELECTOR_INDEX_MISMATCH',stored=receipt.get('raw_candidate'),expected=expected_raw)
    prefix_rows=[]
    try:stored_prefix=read(folder/'PREFIXES.json')
    except Exception as ex:stored_prefix=[];issue(issues,'PREFIX_LOAD_FAILED',error=repr(ex))
    mapping={float(x['CPU_cap']):x.get('candidate') for x in stored_prefix}
    valid_by_index={h['candidate']:h for h in valid}
    for point in protocol.get('prefixes',[]):
        expected=minimum([h for h in valid if h['CPU_seconds']<=point],lambda h:h['Phi']);stored=mapping.get(float(point))
        check=dict(CPU_prefix=point,stored=stored,expected=expected,exact_match=stored==expected,stored_candidate_CPU=valid_by_index.get(stored,{}).get('CPU_seconds'))
        if float(point) not in mapping:issue(issues,'MISSING_PREFIX',CPU_prefix=point)
        if stored!=expected:issue(issues,'PREFIX_SELECTOR_MISMATCH',**check)
        if stored is not None and (stored not in valid_by_index or valid_by_index[stored]['CPU_seconds']>point):
            issue(issues,'PREFIX_FUTURE_OR_INELIGIBLE_CANDIDATE',**check)
        prefix_rows.append(check)
    row['selector_checks']['prefixes']=prefix_rows
    initial=None
    for view in VIEWS:
        path=folder/(view+'.npz');meta_path=folder/(view+'.json')
        if not path.is_file():
            row['snapshots'][view]=dict(available=False,reason='missing snapshot; optional secondary may have no feasible candidate')
            if view=='primary' and expected_primary is not None:issue(issues,'MISSING_SELECTED_PRIMARY_SNAPSHOT')
            continue
        try:
            npz_preflight(path,['m','j','Q','Y_prediction'],max_bytes)
            with np.load(path,allow_pickle=False) as z:state={k:z[k].copy() for k in ('m','j','Q','Y_prediction')}
            meta=read(meta_path) if meta_path.is_file() else {};index=meta.get('candidate')
            entry=dict(available=True,candidate=index,metadata_present=bool(meta))
            if isinstance(index,bool) or not isinstance(index,int):
                issue(issues,'INVALID_SNAPSHOT_CANDIDATE_INDEX',view=view,value=index);index=None
            if index not in all_indices:issue(issues,'SNAPSHOT_INDEX_NOT_IN_HISTORY',view=view,candidate=index)
            elif mats is not None:
                comparison=maxdiff(state['m'],mats[all_indices[index]])
                entry['material_matches_candidate']=comparison
                if not comparison['exact_equal']:issue(issues,'SNAPSHOT_MATERIAL_INDEX_MISMATCH',view=view,candidate=index)
            if view in ('primary','secondary') and index!=receipt.get(view+'_candidate'):
                issue(issues,'SNAPSHOT_RECEIPT_INDEX_MISMATCH',view=view,snapshot=index,receipt=receipt.get(view+'_candidate'))
            if view=='raw' and index!=receipt.get('raw_candidate'):
                issue(issues,'RAW_SNAPSHOT_RECEIPT_INDEX_MISMATCH',snapshot=index,receipt=receipt.get('raw_candidate'))
            if view=='corrected':
                row['selector_checks']['corrected']=dict(stored=index,expected=expected_corrected)
                if index!=expected_corrected:issue(issues,'CORRECTED_SELECTOR_INDEX_MISMATCH',snapshot=index,expected=expected_corrected)
            if not all(np.isfinite(value).all() for value in state.values()):
                issue(issues,'NONFINITE_SNAPSHOT_FIELDS',view=view)
            if view=='initial' and p is not None:
                entry['public_m_init_comparison']=maxdiff(state['m'],p['m_init'])
                if not entry['public_m_init_comparison']['exact_equal']:issue(issues,'INITIAL_MATERIAL_DIFFERS_FROM_PUBLIC_M_INIT')
                initial=state
            entry['defects']=state_defects(state,p,actions) if p is not None else None
            entry['unknown_defect_reason']=None if p is not None else 'public operators unavailable'
            row['snapshots'][view]=entry
        except Exception as ex:issue(issues,'SNAPSHOT_AUDIT_FAILED',view=view,error=repr(ex));row['snapshots'][view]=dict(available=False,reason=repr(ex))
    physical=receipt.get('physical_cost') or{};external=(job.get('external') or{}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))
    auxiliary=receipt.get('auxiliary_receipt') or{}
    row.update(run_status=receipt.get('status'),candidate_rows=len(hist),candidate_material_rows=len(mats) if mats is not None else None,
        valid_candidate_rows=len(valid),termination=receipt.get('termination'),
        cost_coverage=dict(**canonical_cost(physical),internal_process_CPU_seconds=receipt.get('whole_process_CPU_seconds'),
            external_child_CPU_seconds=external,external_unknown_reason=None if external is not None else 'external exit CPU absent; not imputed',
            physical_cost=physical,projector_cost=receipt.get('projector_cost'),
            auxiliary_operation_counts=auxiliary.get('operation_counts'),
            nested_scope='auxiliary physical deltas already included in final physical_cost; never added twice'),
        stage_coverage=dict(current_spaces=auxiliary.get('current_spaces'),PR_steps=(auxiliary.get('operation_counts') or{}).get('PR_CG_steps'),
            status=auxiliary.get('status'),qualification='space setup/reached stage is not evidence that its iterations or physical stage audit completed'),
        DOF_metadata={k:receipt.get(k) for k in ('material_DOF','stored_field_real_scalars','independent_algorithm_state_DOF','state_real_DOF','state_DOF_qualification')},
        audit_CPU_seconds=time.process_time()-started,audit_wall_seconds=time.perf_counter()-wall)
    return row,initial


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--max-array-MiB',type=float,default=512.)
    args=ap.parse_args();root=args.root.resolve();manifest=read(args.manifest);campaign=manifest['campaign']
    queue_path=root/'queues'/campaign/'RECEIPT.json';queue=read(queue_path)
    if queue.get('status') not in TERMINAL or queue.get('pending') or queue.get('running'):
        raise RuntimeError('Public audit requires an explicit terminal queue with no running/pending jobs')
    output=args.output.resolve()
    if output.is_relative_to(root/'RESULTS'/campaign) or output.is_relative_to(root/'queues'/campaign):
        raise ValueError('Audit output must be outside immutable prediction/queue directories')
    output.mkdir(parents=True,exist_ok=False)
    max_bytes=int(args.max_array_MiB*1024**2)
    if max_bytes<=0:raise ValueError('Positive array memory cap required')
    groups=defaultdict(list)
    for task in manifest['tasks']:groups[task['case_id']].append(task)
    jobs={x['task_id']:x for x in queue.get('jobs',[])};rows=[];comparisons=[];actions=Counter();peak_case_public_bytes=0
    public_groups=0
    for case,tasks in groups.items():
        # Cache/deduplicate only this case's public operators, including separate
        # NN initializers; discard the entire case cache before the next case.
        cache={};starts=[];pooled=[];public_bytes=0
        for task in tasks:
            path=(root/task['input_file'].replace('\\','/')).resolve();public=None;public_error=None
            if path in cache:public,public_error=cache[path]
            else:
                try:
                    public=load_public(path,max_bytes)
                    for key,array in list(public.items()):
                        shared=next((old for old in pooled if old.shape==array.shape and old.dtype==array.dtype and np.array_equal(old,array)),None)
                        if shared is not None:public[key]=shared
                        else:
                            if public_bytes+array.nbytes>max_bytes:raise MemoryError('Per-case public array memory cap exceeded')
                            pooled.append(array);public_bytes+=array.nbytes
                    cache[path]=(public,None)
                except Exception as ex:public=None;public_error=repr(ex);cache[path]=(None,public_error)
            task_cpu,task_wall=time.process_time(),time.perf_counter()
            try:
                if task.get('campaign')!=campaign:
                    raise ValueError('Task campaign differs from the terminal queue/manifest campaign; no other result directory opened')
                row,initial=audit_task(root,task,manifest['protocol'],public,jobs.get(task['task_id'],{}),actions,max_bytes)
            except Exception as ex:
                # One malformed/missing task never removes the remaining
                # expected denominator or resets already spent action counts.
                row=dict(task_id=task['task_id'],case_id=task['case_id'],method=task['method'],
                    issues=[dict(kind='TASK_AUDIT_EXCEPTION',error=repr(ex))],snapshots=None,selector_checks=None,
                    cost_coverage=None,stage_coverage=None,run_status=None,truth_access='NONE',
                    audit_CPU_seconds=time.process_time()-task_cpu,audit_wall_seconds=time.perf_counter()-task_wall)
                initial=None
            if public_error:issue(row['issues'],'PUBLIC_OPERATOR_LOAD_FAILED',error=public_error)
            if public is not None and initial is not None:
                family=initial_family(task)
                match=next((g for g in starts if g['family']==family and same_start(g['public'],public)),None)
                if match is None:
                    public_groups+=1;match=dict(group=public_groups,family=family,public=public,initial=initial,task_id=task['task_id']);starts.append(match)
                comparison=dict(case_id=case,group=match['group'],initialization_family=family,
                    reference_task=match['task_id'],task_id=task['task_id'],group_rule='exact public G/E/receiver map/components/m_init; NN family kept separate',
                    fields={k:maxdiff(initial[k],match['initial'][k]) for k in ('m','j','Q','Y_prediction')},
                    interpretation='actual differences recorded; floating field defects have no arbitrary pass threshold')
                comparisons.append(comparison);row['initial_comparison_group']=match['group']
            else:row['initial_comparison_group']=None
            rows.append(row)
        peak_case_public_bytes=max(peak_case_public_bytes,public_bytes)
    report=dict(status='PUBLIC_PREDICTION_ACCEPTANCE_AUDIT_COMPLETED',campaign=campaign,truth_access='NONE',
        expected_tasks=len(manifest['tasks']),audited_task_rows=len(rows),tasks_with_issues=sum(bool(r['issues']) for r in rows),
        queue_status=queue['status'],tasks=rows,initial_state_comparisons=comparisons,
        validation_cost=dict(new_Maxwell_factorizations=0,new_optimization=0,matvec_actions=dict(actions),
            CPU_seconds=time.process_time()-_CPU_START,wall_seconds=time.perf_counter()-_WALL_START,
            timer_scope='script-entry timer before NumPy import, input/metadata/state reads, comparisons and matrix actions; final output serialization/exit excluded',
            peak_cached_case_public_bytes=peak_case_public_bytes,declared_per_case_public_and_per_candidate_array_cap_bytes=max_bytes,
            memory_scope='per-case deduplicated public arrays plus one task candidate array and saved snapshots; not measured process peak RSS'),
        qualification='No GT, fresh quality score, Maxwell solve, inverse iteration or numerical error threshold. Exact index/wiring issues and numerical defects are separate.')
    write(output/'AUDIT.json',report)
    lines=['# Public prediction acceptance audit','',f'Campaign `{campaign}`; {len(rows)}/{len(manifest["tasks"])} expected task rows retained; {report["tasks_with_issues"]} rows have reported wiring/availability issues.',
        '', 'No truth, new Maxwell factorization, optimization or numerical-defect pass gate. Numerical field residuals and exact selector/index issues remain separate. See AUDIT.json for per-view units, relative norms, initial-state max differences, selector/prefix checks, missing counters and stage termination metadata.',
        '', 'Validation matrix actions/CPU are a separate cost, never added to solver physical LU/RHS or used to retroactively create candidates. Shared starts require exact public operators and material initializer; NN groups are separate. Missing tasks/costs are retained as missing/null, never as zero.', '', 'Reported issues:']
    for row in rows:
        for item in row['issues']:lines.append('- `'+row['task_id']+'`: '+json.dumps(item,default=encode))
    if not any(r['issues'] for r in rows):lines.append('- No exact wiring/index/availability issue detected in these inputs; numerical defects are recorded without a threshold.')
    (output/'README.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:report[k] for k in ('status','campaign','expected_tasks','audited_task_rows','tasks_with_issues','truth_access')}),flush=True)


if __name__=='__main__':main()
