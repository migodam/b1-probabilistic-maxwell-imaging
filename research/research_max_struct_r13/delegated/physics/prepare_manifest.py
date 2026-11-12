"""Declare the finite physics campaign; no computation, result or truth reads."""
from pathlib import Path
import json
from datetime import datetime,timezone
HERE=Path(__file__).resolve().parent;STUDY=HERE.parents[1]


def main():
    ordinary=['cie','dual_cie','global_dual','dr','ap','event_integer','direct']
    classical=['cie','dual_cie','global_dual','g19','g19_full_apcs','g19_full_current','dr','event_integer','direct']
    vector=['slotted_high_frame','unequal_coupled_steps','mixed_crossed_bracket','slotted_high_frame_crossgrid8']
    measured=['dielTM_dec8f','twodielTM_8f'];classic=['austria_150to32','austria_matched32']
    tasks=[];truth_specs={};heldout_specs={}
    for index,cid in enumerate(vector+measured+classic):
        if cid in vector:
            family='vector3d';input_file=f'delegated/physics/new3d/runner_public/{cid}.npz';cap=180;methods=ordinary
            spec=dict(path=f'delegated/physics/new3d/sealed_truth/{cid}.npz',key='chi',
                support_rule='fixed_chi_gt_1',support_threshold=1.,amplitude=8.,host=0.,
                amplitude_rule='predeclared_characteristic_chi8',object_id=cid.split('_crossgrid')[0],
                material_truth_available=True,truth_scope='generator coarse volume-average chi, offline after seal only',
                crossgrid_pair_not_independent=cid.endswith('_crossgrid8'))
        elif cid in measured:
            family='fresnel';input_file=f'delegated/physics/fresnel_public/{cid}.npz';cap=180;methods=ordinary
            spec=dict(path=None,material_truth_available=False,support_rule='NO_PIXEL_GT',object_id=cid)
            heldout_specs[cid]=dict(path=f'delegated/physics/fresnel_public/{cid}_HELDOUT_OFFLINE_ONLY.npz',
                measured_values_key='Y',mask_key='mask',full_incident_key='E0_all_sources',
                full_receiver_key='S_all_receivers',heldout_channels=711,
                use='prediction evaluation only after all physics_v1 predictions and choices sealed',
                pixel_GT=False,online_selection=False,
                physical_evaluation_requires_new_forward=True,physical_costs_must_be_recorded=True)
        else:
            family='classic';input_file=f'delegated/physics/classic_public/{cid}.npz';cap=240;methods=classical
            spec=dict(path=f'delegated/physics/classic_sealed_truth/{cid}_OFFLINE_ONLY.npz',key='m',
                support_rule='fixed_chi_gt_1',support_threshold=1.,amplitude=3.,host=0.,
                amplitude_rule='known physical chi=epsilon_r-1=3',object_id='austria_eps4_400MHz',
                material_truth_available=True,truth_scope='existing coarse material copied opaquely; offline after seal only',
                matched_and_crossgrid_same_object_not_independent=True)
        truth_specs[cid]=spec
        for method in methods:
            task=dict(task_id=f'{cid}__{method}',case_id=cid,method=method,campaign='physics_v1',
                input_file=input_file,CPU_cap=cap,seed=13013,family=family,truth_spec=spec)
            if cid in heldout_specs:task['heldout_spec']=heldout_specs[cid]
            tasks.append(task)
    assert len(tasks)==60 and len({t['task_id'] for t in tasks})==60
    missing=[t['input_file'] for t in tasks if not(STUDY/t['input_file']).is_file()]
    if missing:raise FileNotFoundError(sorted(set(missing)))
    # Path availability is checked without opening any truth/heldout archive.
    for spec in list(truth_specs.values())+list(heldout_specs.values()):
        if spec.get('path') and not(STUDY/spec['path']).is_file():raise FileNotFoundError(spec['path'])
    manifest=dict(campaign='physics_v1',status='READY',
        created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_status='FIXED_BEFORE_PHYSICS_V1_INVERSION',
        source_freeze='parent freezes C2 with shared2CPU-second proposal audit clock before ANYphysics_v1 runs; existing development runtimes unchanged',
        execution_performed_by_preparer=False,
        protocol=dict(CPU_cap=180,task_overrides='classic240; all other families180',
            proposal_fraction=.6,steps=3000,audit_every=50,AP_preparation=128,scale=6,omega=1,
            physical_audit_CPU_interval=2.,
            physical_audit_scope='proposal periodic candidates use common absolute process CPU audit_due clock, replacing iteration cadence for all physics families; original native initial/stage/final audits retained; corrector objective evaluations audited as paid',
            corrector='common original physical L-BFGS-B, same finite allowances and caller whole-process eligibility',
            selection='minimum original paid physical Phi; secondary minimum prior among candidates in the predeclared ball',
            primary='all eligible paid physical candidates including original initial material',
            material_box=[0,20],material_box_by_family=dict(vector3d=[0,8],fresnel=[0,20],classic=[0,20]),
            prior_by_family=dict(vector3d=dict(kind='R8quadratic',edge_mean_square=2e-4,material_mean_square=1e-6),
                fresnel=dict(smooth=1e-4,tv=1e-4,epsilon=.01),classic=dict(smooth=1e-4,tv=1e-4,epsilon=.01)),
            prior='per unchanged qualified model: R8vector edge/l2; scalar smooth1e-4/tv1e-4/epsilon.01',
            random_seed=13013,prefixes=[15,30,60,90,120,180,240],
            scope='three independent new vector geometries plus crossgrid pair; two actual measured inputs; two same-object classic diagnoses',
            outcome='retain all losses/failures and finite budget stops; no pass threshold, winner quota or adaptive case selection',
            ball_qualification='matched vector known iidCN99; crossgrid vector includes unresolved model mismatch; Fresnel/classic fixed discrepancy balls are not CN99 certificates',
            truth_access='NONE online; no material truth or heldout observations read until complete expected-task seal reviewed and authorized'),
        queue=dict(implementation='existing finite_queue.py',one_campaign=True,
            maximum_concurrent_processes=4,CPU_cores_only=True,extra_GPU=False,
            retries=False,task_count=60,total_declared_CPU_cap_seconds=sum(t['CPU_cap'] for t in tasks),
            queue_order='fixed family,case,method order; finite ready queue; no outcome-dependent adaptation'),
        family_counts=dict(vector3d=dict(cases=4,tasks=28,independent_geometries=3),
            fresnel=dict(cases=2,tasks=14,pixel_GT=False),
            classic=dict(cases=2,tasks=18,independent_physical_object=1)),
        adapter_contracts=dict(vector3d='VectorPublicModel; native3D VectorPR3D CIE/Dual/global; source-normalized physical objective',
            fresnel='SourceSpecificPublicModel; projector_S list; R8 source-specific calibrated measured CIE/Dual/global',
            classic='classic_auxiliary_config(public,method,model,Y,out); physics_family equal_area_disk; original stage/rank/box/prior',
            classic_G19='physical warm R8loop with EqualAreaGreen/PaperOneHop; explicit global stop convention; full masks16 on32²',
            classic_mature='cold paperloop; original perstage coefficient stops; DPCS15,mf5/8/10,materialbeta8/2/.5,Dualcurrent8',
            event_integer='128AP preparation then finite integer normal event and ordinaryDR; common steps/audits/model'),
        qualification_report='delegated/physics/QUALIFICATION_DIFFERENCES.md',
        truth_specs=truth_specs,heldout_specs=heldout_specs,tasks=tasks)
    path=STUDY/'configs/physics_v1.json'
    with path.open('x') as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
    ready=dict(status='READY',manifest='configs/physics_v1.json',task_count=len(tasks),
        unique_task_count=len({t['task_id'] for t in tasks}),public_inputs_available=True,
        truth_and_heldout_path_metadata_only=True,truth_arrays_read=False,
        new_physical_solves=0,execution_performed=False,
        parent_integration_required='classic scalar branch calls classic_auxiliary_config; final shared sources frozen before one4corequeue',
        tasks=tasks)
    (HERE/'PHYSICS_V1_READY.json').write_text(json.dumps(ready,indent=2)+'\n')
    print(json.dumps(dict(status='READY',manifest=str(path),tasks=len(tasks),declared_CPU_cap_seconds=sum(t['CPU_cap'] for t in tasks))))
if __name__=='__main__':main()
