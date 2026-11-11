"""Copy allowlisted public inputs only; no physics, optimization, trigger or truth reads."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
REPO=STUDY.parents[1]
R10=REPO.parent/'r10_worktree/research/research_max_struct_r10'
ALLOW=('G','E0','S','Y','sigma','noise_radius','noise_model','grid_shape','material_bounds','coords','centers','sources','receivers','k','radius')
PRIOR=dict(name='inherited_original_Maxwell_Phi',smooth=1e-4,tv=1e-4,epsilon=.01,source='research/research_z_alg_r4_diversity/src/physics.py:Regularizer')
READS=[]
def public(path):
    path=Path(path)
    if 'sealed_truth' in path.parts:raise ValueError('Private input forbidden')
    READS.append(str(path))
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files if k in ALLOW}
def material(path,key):
    path=Path(path)
    if 'sealed_truth' in path.parts:raise ValueError('Private input forbidden')
    READS.append(str(path))
    with np.load(path,allow_pickle=False) as z:return z[key].copy()
def read(path):
    READS.append(str(path));return json.loads(Path(path).read_text())
def dump(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def costs(receipt):
    return {k:receipt[k] for k in ('cpu_seconds','whole_child_CPU_seconds','wall_seconds','physical_model_cost','auxiliary_physical_model_cost','setup_cpu_seconds','auxiliary_cpu_seconds','physical_operation_cpu','physical_operation_wall','cost_scope') if k in receipt}
def save(cohort,cid,pub,m,provenance):
    assert m.shape==(576,) and np.isfinite(m).all()
    assert pub['G'].shape==(576,576) and pub['E0'].shape[0]==576
    assert pub['Y'].shape==(pub['S'].shape[0],pub['E0'].shape[1])
    arrays=dict(pub,m_init=m)
    if 'radius' in pub:arrays['cell_radius']=pub['radius'].copy()
    # Store only allowlisted public physical inputs; material bytes unchanged.
    dest=HERE/'inputs'/f'{cohort}__{cid}.npz'
    np.savez_compressed(dest,**arrays)
    with np.load(dest,allow_pickle=False) as z:
        assert set(z.files)==set(arrays)
        for key,val in arrays.items():assert np.array_equal(z[key],val)
    return dict(cohort=cohort,case_id=cid,input_file=str(dest.relative_to(HERE)),
                source_count=int(pub['E0'].shape[1]),receiver_count=int(pub['S'].shape[0]),
                material_cells=576,material_unchanged=True,keys=list(arrays),
                public_sigma=float(pub['sigma']),public_noise_radius=float(pub['noise_radius']),
                raw_CN99_radius=float(pub['sigma']*pub['noise_radius']),
                prior=PRIOR,trigger_evaluated=False,**provenance)
def main():
    rows=[]
    for cid,arm,label in [('r10_twod_0001','cie_single','ring'),('r10_twod_0004','matrix_green19','shell'),('r10_twod_0022','matrix_green19','bridge')]:
        base=STUDY/'data_public/basin';op=base/'operators.npz';obs=base/f'{cid}.npz';init=base/f'{cid}_initial.npz'
        pub=public(op);pub.update(public(obs));m=material(init,'m')
        parent=R10/'RESULTS/twod'/f'{cid}__{arm}__noise_ball_99'
        assert np.array_equal(m,material(parent/'prediction.npz','m'))
        receipt=read(parent/'RECEIPT.json')
        rows.append(save('A',cid,pub,m,dict(label=label,role='three fixed R11-exposed bad endpoints inherited from actual R10 predictions',public_sources=[str(op),str(obs)],initialization_source=str(init),parent_task_id=receipt['task_id'],parent_method=arm,parent_prediction=str(parent/'prediction.npz'),parent_receipt=str(parent/'RECEIPT.json'),parent_original_actual_cost=costs(receipt),parent_material_exact_equality=True,trigger_rule='ALL_THREE_FIXED_EXPOSED_ENDPOINTS',trigger_required=True)))
    c=STUDY/'delegated/confirmation';seal=read(c/'SEALED_CONFIRM24_PUBLIC_COHORT.json')
    parents=[r for r in seal['rows'] if r['method']=='fullG19']
    assert [r['case_id'] for r in parents]==[f'r12_confirm24_{i:04}' for i in range(24)]
    for row in parents:
        cid=row['case_id'];base=c/'selected_public_outputs/confirm24'/row['task_id'];init=base/'primary.npz'
        op=c/'data_public/confirm24/operators.npz';obs=op.parent/f'{cid}.npz'
        pub=public(op);pub.update(public(obs));m=material(init,'m');receipt=read(base/'RECEIPT.json')
        cap=read(base/'R11_CAP_PREDICTION.json')
        rows.append(save('B',cid,pub,m,dict(role='all24 original primary fullG19 endpoints; no post-trigger filtering',public_sources=[str(op),str(obs)],initialization_source=str(init),parent_task_id=row['task_id'],parent_method='fullG19',parent_prediction_original=row['prediction'],parent_receipt_original=row['receipt'],parent_receipt_local=str(base/'RECEIPT.json'),compact_export_source=str(c/'export_selected_public.py'),parent_primary_Phi=row['Phi'],parent_paid_state=cap,parent_original_actual_cost=dict(**costs(receipt),completed_external_child_CPU_seconds=row['completed_external_child_CPU_seconds']),trigger_rule='||Y_forward(m_init)-Y||_F / sigma > noise_radius',trigger_required='DEFERRED_TO_RUNNER_PUBLIC_FORWARD',noise_convention='iid circular complex RMS sigma; noise_radius=sqrt(chi2.ppf(.99,2*Y.size)/2); use stored original radius; strict >',parent_secondary_not_used=True)))
    n=STUDY/'delegated/neural';prep=read(n/'PREPROCESS_RECEIPT.json');protocol=read(n/'FROZEN_PROTOCOL.json')
    for cid in ['r8_test_0001','r8_test_0004','r8_test_0022','r8_test_0027']:
        obs=n/'public'/f'{cid}.npz';init=obs.with_name(f'{cid}_nn.npz');pub=public(obs);m=material(init,'m_init');p=next(r for r in prep['records'] if r['case_id']==cid)
        rows.append(save('C',cid,pub,m,dict(role='four fixed exposed frozen R8 NN initial materials',public_sources=[str(obs)],initialization_source=str(init),parent_task_id=cid+'__frozen_nn_initializer',parent_method='R8 frozen physics NN',parent_original_actual_cost={'preprocessing_record':p,'shared_preprocess_summary':{k:v for k,v in prep.items() if k!='records'}},original_frozen_task=next(t for t in protocol['tasks'] if t['case_id']==cid) if 'tasks' in protocol else None,trigger_rule='ALL_FOUR_FROZEN_NN_INITIALIZATIONS',trigger_required=True)))
    for row in rows:
        row['file']=row['input_file'].removeprefix('inputs/')
        row['sigma']=row['public_sigma'];row['delta']=row['raw_CN99_radius']
        row['parentcost']=row['parent_original_actual_cost']
        if row['cohort']=='B':
            with np.load(HERE/row['input_file'],allow_pickle=False) as z:ynorm=float(np.linalg.norm(z['Y']))
            residual=row['parent_paid_state']['physical_relative_residual']*ynorm
            row['original_primary_raw_residual_from_metadata']=residual
            row['original_primary_whitened_residual_from_metadata']=residual/row['sigma']
            row['trigger_from_original_primary_metadata']=bool(residual>row['delta'])
            row['trigger_evaluated']=True
            row['trigger_evaluation_kind']='saved original primary metadata only; no new forward'
            row['trigger_metadata_source']='original R11_CAP_PREDICTION.json physical_relative_residual times exact public normY; no forward; original primary m retained'
    dump(HERE/'inputs/manifest.json',dict(cases=rows,counts={'A':3,'B':24,'C':4,'total':31},path_base='inputs directory',truth_access='NONE',new_forward=0))
    result=dict(status='PUBLIC_INPUTS_PREPARED_NO_NEW_PHYSICS',counts={'A':3,'B':24,'C':4,'total':31},truth_access='NONE',new_forward=0,new_optimization=0,trigger_evaluation='B_CLASSIFIED_FROM_SAVED_PRIMARY_RESIDUAL_METADATA_NO_FORWARD',array_policy='Exact material and public arrays; allowlist only; no private arrays; no material clipping/casting',rows=rows,source_reads=READS,uncertainties=['A are R11-exposed bad endpoints with exact producing R10 prediction match; original nominal R11 wording must not replace actual provenance.','B public compact primary export carries exact original primary material; original full remote arrays remain separate.','Original parent CPU views overlap; completed external child CPU is inclusive parent deployment cost for B.','C preprocessing shared cost is a cohort cost; do not add four copies of the shared summary.'])
    dump(HERE/'audit/INPUTS_PROVENANCE.json',result)
    text=['# Public input preparation','', '31 fixed materials: A=3, B=24, C=4. No physics, optimization, GT or sealed-truth reads. B triggers classified from saved original primary residual metadata; no new forward.','', 'NPZ mandatory: G, E0, S, Y, m_init. Public sigma/noise_radius, grid/material bounds and available geometry are preserved. Array equality checked after serialization.','', 'B retains all24 independent confirmation original fullG19 primary materials. Strict ||Y(m_init)-Y||/sigma > stored noise_radius classification is provided from original primary residual metadata; preserve untriggered rows as NOT_RUN. Original prior smooth=tv=1e-4, epsilon=.01.','', 'A inputs exactly match producing R10 cie_single ring0001 and matrix_green19 shell0004/bridge0022 predictions; R11 exposed these endpoints.','', 'Detailed source paths, task IDs, original parent cost views, public noise and initialization provenance: INPUTS_PROVENANCE.json. No SHA checks.']
    (HERE/'audit/INPUTS_README.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(dict(status=result['status'],counts=result['counts'],GT='NONE',new_forward=0)))
if __name__=='__main__':main()
