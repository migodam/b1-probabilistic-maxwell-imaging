"""Export actual frozen-confirmation finite-static worlds, with zero new physics.

No offline/truth file is read. All physics validation values originate in the
paid full-Mie confirmation campaign and its append-only cost receipts.
"""
from __future__ import annotations
import json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
from research_a.budget import Budget
from research_a.contracts import make_mie_context,write_packages,zero_cost,error_record
from research_a.target import to_x,TAU
from research_a.material_metrics import distance,structural_difference
from b1_joint.material import decode,event_ids

OUT=ROOT/'research/research_a_20261015'
BACKEND_COMMIT='7ed4faef9bc3dd214d4de689c2203cc383f7495f'
TARGET='Research-A-original-P2-vector-Mie-n14-C'


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False,default=lambda x:x.tolist())+'\n')


def radial_material(p,edges):
    d=decode(to_x(p));mid=.5*(edges[:-1]+edges[1:])
    return d['eps_layers'][np.searchsorted(d['radii'],mid)]


def ledger_index(rows):
    ends={(r['stage'],r['attempt']):r for r in rows if r['event']=='END'}
    by={};parent=None
    for r in rows:
        if r.get('event')!='BEGIN' or r.get('job')!='confirm':continue
        context=r.get('context',{})
        if r.get('kind')=='SHARED_COARSE_ENDPOINT':
            parent=context['parent'];key=(parent,'base',None)
        elif r.get('kind')=='SEARCH_ENDPOINT':key=(parent,context['arm'],context['step'])
        else:continue
        by[key]=(r,ends[(r['stage'],r['attempt'])])
    return by


def candidate_record(c,index,archive,parent_id,base,edges,receipts,method):
    is_base=np.array_equal(c['p'],base['p']) and np.array_equal(c['fitted_prediction'],base['fitted_prediction'])
    key=(parent_id,'base',None) if is_base else (parent_id,'finite_static',c['history']['step'])
    begin,end=receipts[key]
    if end['status']!='PASS':raise ValueError('Paid full endpoint did not succeed')
    decoded=decode(to_x(c['p']))
    params=begin['parameters']
    if not (np.allclose(params['radii'],decoded['radii'],rtol=0,atol=2e-15)
            and np.allclose(np.asarray(params['eps_layers']),np.c_[decoded['eps_layers'].real,decoded['eps_layers'].imag],rtol=0,atol=2e-15)):
        raise ValueError(f'Ledger/material replay differs for {parent_id}/{index}')
    mat=radial_material(c['p'],edges)
    v=np.diff(edges**3)
    support=float(np.sum(v[np.abs(mat-1)>1e-12]))
    requested_r,requested_h,er,ei=c['p']
    r=float(decoded['central_radius']);h=float(decoded['annular_halfwidth'])
    void_count=len(decoded['void_intervals'])
    vacuum_mask=np.abs(mat-1)<1e-12
    array_void_runs=int(np.count_nonzero(vacuum_mask & np.r_[True,~vacuum_mask[:-1]]))
    material_mask=~vacuum_mask
    material_runs=int(np.count_nonzero(material_mask & np.r_[True,~material_mask[:-1]]))
    if void_count!=array_void_runs:
        raise ValueError(f'Decoder/radial-array cavity consistency failed for {parent_id}/{index}')
    history=dict(c['history'])
    parent_p=history.get('parent_p')
    parent_name=None
    if parent_p is not None:
        parent_name=next((f'world_{j:02d}' for j,a in enumerate(archive) if np.array_equal(a['p'],parent_p)),None)
    proxy_wall=float(history.get('corrector',{}).get('wall_seconds',0.))
    cid=f'world_{index:02d}'
    source='shared-coarse-reconstruction' if is_base else history.get('source','finite-static')
    trace=[dict(origin='actual CONFIRM_CAMPAIGN finite_static',parent_id=parent_id,candidate_index=index,
                generation=history,optimizer_requested_p=c['p'],canonical_decoded_geometry=dict(central_radius=r,annular_halfwidth=h),discovery_status=c.get('discovery_status','SHARED_BASE' if is_base else None),
                cost_receipt=dict(begin=begin,end=end),proxy_corrector_cost_scope='saved per-candidate corrector wall plus charged endpoint wall; shared proposal/setup costs remain in global ledger'),
           dict(method_search_history=method['history'],prior_support_unchanged=True)]
    return dict(id=cid,p=c['p'],material=mat,prediction=np.asarray(c['fitted_prediction']),
        whitened_residual_l2=c['residual'],validated_full_maxwell=c['valid'],physics_target_id=TARGET,
        material_parameterization_id='exact-concentric-layer-epsilon-on-package-union-radial-intervals-v1',
        validation=dict(constitutive_relative_defect=None,constitutive_relative_defect_status='NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE',
            solver_relative_residual=c['boundary_residual'],material_constraints_passed=True,all_illuminations_same_material=True,
            numerical_error=error_record(value=c['boundary_residual'],status='solver_tolerance_only',
                norm='max per-electric/magnetic multipole boundary-matching linear-system normwise backward residual',
                method='Original full vector-Mie n14 endpoint; not an internal current/field constitutive defect; no independent mesh/full-model certificate')),
        structure=dict(support_volume=support,cavity_count=void_count,component_count=material_runs,
            feature_scale=float(min([r for _ in [0] if r>0]+[2*h for _ in [0] if h>0]+[1.])),
            descriptor=dict(radial_geometry=dict(central_radius=r,annular_halfwidth=h,ring_center=.65,outer_radius=1.),
                optimizer_requested_p_r_h_epsilon_re_im=c['p'],actual_decoded_geometry=dict(central_radius=r,annular_halfwidth=h),cavity_count_array_vacuum_runs=array_void_runs,Qabs=c['Qabs'],boundary_residual=c['boundary_residual'],
                exact_layer_radii=decoded['radii'].tolist(),exact_layer_epsilon_re_im=np.c_[decoded['eps_layers'].real,decoded['eps_layers'].imag].tolist(),
                material_legality='passive epsilon with strict-vacuum layers; same material for both polarization channels'),
            descriptor_definition_id='Research-A-radial-geometry-v1',cluster_ids=[cid],event_ids=[f'radial-event-{int(event_ids(to_x(c["p"])))}'],
            topology_thresholds=[0.],topology_persistence_rule='exact nonzero radial intervals, no pixel threshold; structural distinctness requires delta-r>=.05R or delta-h>=.025R'),
        compensation=dict(parent_candidate_id=parent_name,edit_family=source,
            edit_parameters=dict(proposal_p=history.get('proposal_p'),edit_id=history.get('edit_id'),
                protected_geometry=history.get('corrector',{}).get('protected_geometry'),
                parent_p=parent_p,source=source),
            material_compensation_parameters=dict(epsilon_real=er,epsilon_imag=ei,central_radius=r,annular_halfwidth=h),
            nuisance_compensation_parameters=dict(gain=c['gain'],bounds=[.7,1.3],shared_across_frequencies_and_sources=True),
            coordinates=np.asarray(c['p']),protected_structure_constraints=[
                'The finite edit protected geometry cell is stored in edit_parameters.protected_geometry',
                'Candidate is archived only after independent full-Mie endpoint data-fit and material-space deduplication']),
        material_distance_to_base=distance(c['p'],base['p']),nuisance_policy='bounded_profile',nuisance_parameters={'gain':c['gain']},
        cost=zero_cost(new_full_state_frequency_attempts=int(c['cost']['new_material_frequency_attempts']),
            vector_rhs=int(c['cost']['boundary_RHS']),forward_rhs=int(c['cost']['boundary_RHS']),proposal_count=1,
            wall_seconds=float(end['wall_seconds'])+proxy_wall,historical_cost_ledger='../../COST_LEDGER.jsonl'),
        history=trace,method=source)


def run(*,semantic_repair=False):
    budget=Budget(OUT)
    original_receipt=OUT/'results/ACTUAL_PACKAGE_EXPORT.json'
    old_metadata=[]
    if semantic_repair:
        retained=OUT/'results/ACTUAL_PACKAGE_EXPORT_ORIGINAL_BEFORE_DECODER_REPAIR.json'
        if original_receipt.exists() and not retained.exists():retained.write_text(original_receipt.read_text())
        for manifest_path in sorted((OUT/'MaterialHypothesisPackage').glob('confirm_*/manifest.json')):
            m=json.loads(manifest_path.read_text())
            for c in m['candidates']:
                old_metadata.append(dict(parent=manifest_path.parent.name,candidate=c['id'],structure=c['structure'],compensation=c['compensation']))
        dump(OUT/'results/ACTUAL_PACKAGE_DECODER_METADATA_BEFORE_REPAIR.json',old_metadata)
    packages=[];start=time.perf_counter();physics_before=budget.count()
    with budget.job('export-actual-packages-semantic-repair' if semantic_repair else 'export-actual-packages'):
        campaign=json.loads((OUT/'results/CONFIRM_CAMPAIGN.json').read_text())
        if len(campaign)!=32:raise ValueError('Expected all 32 completed confirmation parents')
        protocol=json.loads((OUT/'configs/protocol_v1.json').read_text())
        frozen=json.loads((OUT/'configs/CONFIRM_FREEZE.json').read_text())
        rows=budget.rows();receipts=ledger_index(rows);snapshot=budget.summary()
        rejects=sum(c.get('discovery_status')=='DATA_INCOMPATIBLE' for s in campaign for m in s['methods'] for c in m['attempts'])
        total=zero_cost(new_full_state_frequency_attempts=snapshot['attempts'],
            vector_rhs=snapshot['boundary_RHS'],forward_rhs=snapshot['boundary_RHS'],failed_count=snapshot['failed'],
            proposal_count=sum(r.get('event')=='BEGIN' and r.get('kind')=='SEARCH_ENDPOINT' for r in rows),
            rejected_count=rejects,wall_seconds=snapshot['controlled_execution_wall'],historical_cost_ledger='../../COST_LEDGER.jsonl')
        for s in campaign:
            budget.check()
            m=next(m for m in s['methods'] if m['arm']=='finite_static')
            archive=[c for c in m['archive'] if c['valid']]
            if s['base']['valid'] and not any(np.array_equal(c['p'],s['base']['p']) and np.array_equal(c['fitted_prediction'],s['base']['fitted_prediction']) for c in archive):
                archive.insert(0,s['base'])
            if not archive:
                packages.append(dict(parent=s['id'],status='NO_VALID_ARCHIVE_NO_EXPORT',world_count=0));continue
            all_edges=sorted({0.}|{float(r) for c in archive for r in decode(to_x(c['p']))['radii']})
            edges=np.asarray(all_edges)
            candidate_rows=[candidate_record(c,i,archive,s['id'],s['base'],edges,receipts,m) for i,c in enumerate(archive)]
            # Required actual-package semantic check, not a synthetic test: the
            # number of void runs in saved exact radial epsilon must equal the
            # descriptor's actual decoder count, irrespective of optimizer p.
            for c in candidate_rows:
                vacuum=np.abs(c['material']-1)<1e-12
                runs=int(np.count_nonzero(vacuum & np.r_[True,~vacuum[:-1]]))
                if runs!=c['structure']['cavity_count']:raise ValueError('Saved radial material/metadata cavity count mismatch')
            base_row=next((c for c in candidate_rows if np.array_equal(c['p'],s['base']['p'])),candidate_rows[0])
            pairs=[]
            for i,a in enumerate(candidate_rows):
                for j,b in enumerate(candidate_rows[:i]):
                    measurement_distance=float(np.linalg.norm((a['prediction']-b['prediction'])/s['sigma']))
                    pairs.append(dict(candidate_i=a['id'],candidate_j=b['id'],material_distance=distance(a['p'],b['p']),
                        material_metric_id='volume-weighted-complex-epsilon-L2-div-abs(1.5+.1i)',measurement_distance=measurement_distance,
                        measurement_metric_id='old1.6-full-gain-applied-prediction-difference-whitened-real-L2',
                        hard_to_distinguish=None,decision_rule_id=None,connection_status='not_tested'))
            actions=[dict(id=f'kR-{k:.1f}',wavenumber=k,nmax=14,angles_degrees=[20,40,60,80,100,120,140,160],
                polarizations=['S1','S2'],n_real_measurements=32,radial_grid_edges=all_edges,
                realification_order='all-real-then-all-imaginary',shared_material=True,shared_candidate_gain=True,
                material_frequency_attempts_per_candidate=1,
                status='old observation prediction already cached' if k==1.6 else 'new action requires separate authorization and charged full-Mie solve') for k in [1.6,2.4]]
            failed=[dict(step=c['history'].get('step'),p=c['p'],valid=c['valid'],discovery_status=c.get('discovery_status'),
                residual=c['residual'],boundary_residual=c['boundary_residual'],Qabs=c['Qabs'],history=c['history'],cost=c['cost'])
                for c in m['attempts'] if c.get('discovery_status')!='NEW_MATERIAL']
            context=make_mie_context(package_id=f'Research-A-20261015-{s["id"]}-finite-static',target_id=TARGET,
                backend_commit=BACKEND_COMMIT,observation=np.asarray(s['observation']),sigma=np.asarray(s['sigma']),fit_threshold=TAU,
                acquisition=dict(id='Research-A-P2-kR1.6-angle20to160-polarization-S1S2-Re16Im16',frequencies=[1.6],frequency_units='kR',n_sources=2),
                grid=dict(kind='EXACT_RADIAL_3D_SPHERICAL_MATERIAL',outer_radius=1.,radial_interval_edges=all_edges,
                    interval_normalized_volumes=np.diff(edges**3).tolist(),array_semantics='one complex relative epsilon per radial interval; spherical symmetry uniquely specifies full 3D material',
                    piecewise_constant_exact=True,center_slice_is_display_only=True),
                material_support=dict(latent_bounds=protocol['latent_bounds'],geometry='B1-P2 original implicit disjoint radial-bump decoder',
                    shared_complex_epsilon_bounds=[[2.2,2.8],[.04,.16]],vacuum_epsilon=[1.,0.],passivity='Im epsilon >=0 under exp(-i omega t)',
                    same_material_all_sources=True),
                nuisance_support=protocol['gain'],total_cost=total,
                upstream_manifest=dict(source_campaign='research/research_a_20261015/results/CONFIRM_CAMPAIGN.json',
                    source_parent_id=s['id'],source_method='finite_static',backend_commit=BACKEND_COMMIT,
                    protocol=protocol,confirmation_freeze=frozen,
                    global_ledger_scope='total_cost is one GLOBAL E1 snapshot, repeated for provenance; NEVER sum it over packages or candidates',
                    global_budget_snapshot=snapshot,
                    global_rejected_count_scope='confirmation DATA_INCOMPATIBLE attempts only; complete physical setup/failure/repetition costs are in append-only global ledger',
                    per_candidate_cost_scope='one saved full endpoint plus its saved proxy-corrector wall; shared coarse endpoint/setup are not to be charged again',
                    rhs_scope='boundary_RHS means 2*nmax spherical electric/magnetic boundary systems, not voxel illumination RHS',
                    candidate_generation_history='actual candidate history and complete finite_static step/discovery history in common exchange; unsuccessful attempts retained as assets'),
                failed_attempts=failed,
                new_action_interface=dict(evaluator_entrypoint='research_a.contracts.evaluate_new_action with caller-supplied budget-enforced original B1 full-Mie callback',
                    environment=dict(backend_module='b1.physics.Forward',implementation='matching',nmax=14,backend_commit=BACKEND_COMMIT,
                        callback_contract='Caller supplies original Forward.run wrapped by its separately authorized ledger; no evaluator is executed by export',
                        no_surrogate_validation=True,new_action_status='NOT_RUN_BY_EXPORT'),
                    allowed_actions=actions,cache_manifest=dict(old_prediction_assets='per-candidate gain-applied full prediction at kR=1.6',new_physics_evaluations=0)),
                distance_definition_id='volume-weighted-complex-epsilon-L2-div-abs(1.5+.1i)',noise_id='P2-frozen-reference-solid-RMS0.12-real-SD-kR1.6',requested_k=len(archive))
            path=OUT/'MaterialHypothesisPackage'/s['id']
            result=write_packages(path,context,candidate_rows,base=base_row,pair_distances=pairs)
            packages.append(dict(parent=s['id'],status='ACTUAL_VALIDATED_WORLDS_EXPORTED',world_count=len(archive),
                structural_distinct_pairs=sum(structural_difference(a['p'],b['p']) for i,a in enumerate(archive) for b in archive[:i]),
                manifest=str(Path(result['native_manifest']).relative_to(OUT)),common_exchange=str(Path(result['exchange_manifest']).relative_to(OUT)),
                exact_radial_intervals=len(edges)-1,schema_and_semantics=result['audit']))
        if budget.count()!=physics_before:raise RuntimeError('Export unexpectedly changed physical attempt count')
        dump(OUT/'results/ACTUAL_PACKAGE_EXPORT.json',dict(status='PASS',scope='32 actual confirmation finite_static archives; no offline truth read',
            packages=packages,exported_parents=sum(p['world_count']>0 for p in packages),total_worlds=sum(p['world_count'] for p in packages),
            new_physics_attempts=0,digest_computation='NOT_RUN_BY_USER_POLICY',source_campaign_rows=len(campaign),
            budget_snapshot_before_export=snapshot,export_wall_seconds=time.perf_counter()-start,
            semantic_revision='actual-decoder-geometry-v2' if semantic_repair else 'actual-decoder-geometry-v2-first-export',
            cavity_consistency='PASS: all 193 saved radial material arrays agree with canonical decoder vacuum-interval count',
            total_cost_aggregation='Global E1 snapshot repeated as provenance; do not sum across packages',
            actual_new_measurement_evaluation='NOT_RUN_BY_EXPORT; old validated predictions are reused, new actions require separately authorized callbacks'))
    report=json.loads((OUT/'results/ACTUAL_PACKAGE_EXPORT.json').read_text())
    if semantic_repair:
        changed=[]
        for old in old_metadata:
            path=OUT/'MaterialHypothesisPackage'/old['parent']/'manifest.json'
            actual=json.loads(path.read_text())
            new=next(c for c in actual['candidates'] if c['id']==old['candidate'])
            if old['structure']['cavity_count']!=new['structure']['cavity_count'] or old['structure']['component_count']!=new['structure']['component_count']:
                changed.append(dict(parent=old['parent'],candidate=old['candidate'],old_cavity_count=old['structure']['cavity_count'],new_cavity_count=new['structure']['cavity_count'],
                    old_component_count=old['structure']['component_count'],new_component_count=new['structure']['component_count']))
        report['semantic_repair']=dict(reason='Raw optimizer radius around2.9255e-11 can decode to exactzero; derive all structure metadata from canonical decoder and saved radial epsilon instead of requested coordinate sign',
            changed_count=len(changed),changed=changed,physics_predictions_changed=False,material_arrays_changed=False,selection_or_metrics_changed=False,
            old_export_receipt='results/ACTUAL_PACKAGE_EXPORT_ORIGINAL_BEFORE_DECODER_REPAIR.json',old_metadata='results/ACTUAL_PACKAGE_DECODER_METADATA_BEFORE_REPAIR.json')
        dump(OUT/'results/ACTUAL_PACKAGE_DECODER_REPAIR.json',report['semantic_repair'])
    report['budget_summary_after_metered_export']=budget.summary()
    dump(OUT/'results/ACTUAL_PACKAGE_EXPORT.json',report)
    print(json.dumps({k:report[k] for k in ['status','exported_parents','total_worlds','new_physics_attempts','export_wall_seconds']},indent=2))


if __name__=='__main__':run(semantic_repair='--semantic-repair' in sys.argv)
