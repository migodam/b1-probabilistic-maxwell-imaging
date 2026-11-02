"""Read existing artifacts and standardize receipts; performs no physics.

Never add nested solver/setup/per-case timers to their enclosing process timer.
Raw numerical outputs are preserved. Re-running this file updates annotations
from current saved raw metadata but does not rerun any experiment.
"""
from pathlib import Path
from copy import deepcopy
import json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'audit/receipts'
BASE='7317601b43b914525c21f09e5a60020552ea677b'
FIELDS=['run_id','case_id','branch','git_sha','code_hash','config_hash','seed','data_id',
 'material_class','physical_convention','mesh','sources','receivers','frequencies','noise_model',
 'whitening','prior','initialization','state_representation','state_dof','material_dof','baseline_id',
 'stage','termination_reason','online_selection_rule','physical_residual','auxiliary_residual',
 'constitutive_defect','constraint_violation','cost','candidate_ids','truth_access','validation','status']

def read(p):return json.loads((ROOT/p).read_text())

def write(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False))

def normalize(r,reason):
    r=deepcopy(r)
    for k in FIELDS:r.setdefault(k,None)
    r['git_sha']=BASE;r['code_hash']=None;r['config_hash']=None
    r['git_provenance_note']='git_sha identifies the frozen input/base tree. Newly written study source is retained by referenced paths; this is not a claim that the new code was already committed at execution.'
    explanations=r.setdefault('null_reasons',{})
    explanations.update(code_hash='User forbids hash checks; retained source paths and frozen Git base identify provenance',
                        config_hash='User forbids hash checks; configs/initial_plan.json and raw run configs retained')
    def scan(x,p=''):
        if isinstance(x,dict):
            for k,v in x.items():
                if k=='null_reasons':continue
                q=f'{p}.{k}' if p else k
                if v is None:explanations.setdefault(q,reason)
                else:scan(v,q)
        elif isinstance(x,list):
            for i,v in enumerate(x):
                if v is None:explanations.setdefault(f'{p}[{i}]',reason)
                else:scan(v,f'{p}[{i}]')
    scan(r)
    return r

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    f=read('results/families/SUMMARY.json');b=read('results/backend_qualification/QUALIFICATION.json')
    d=read('research/delegated/audit/dynamics_cost.json')
    index=[]
    def save(r,reason,group='required'):
        r=normalize(r,reason);fn=r['run_id']+'.json';write(OUT/fn,r)
        index.append(dict(run_id=r['run_id'],path='audit/receipts/'+fn,status=r['status'],group=group))
    frows=f['rows']+[f['interpolation'],f['finite_noise']]
    for row in frows:
        r=deepcopy(row)
        r.update(source_references=['src/tensor_families.py','experiments/review_f_metadata.py','results/families/SUMMARY.json'],
          cost_accounting_annotation='Per-case wall includes its nested setup and compute, but excludes final NPZ I/O for four primary cases. All per-case timings and raw solver receipts are nested within SUMMARY process timer. Do not add them to SUMMARY.',
          blind_status='OFFLINE_DIAGNOSTIC_DISCLOSED_ANCHOR',imaging_claim=False)
        save(r,'Not measured or not applicable to this disclosed finite-DDA branch diagnostic; see raw scientific_scope/cost scope')
    for row in b['rows']:
        h=row['h'];sr=row['solve_receipt'];N=row['N']
        r=dict(run_id=f'BACKEND_SPHERE_h{h:.2f}',case_id=f'sphere_r.12_eps3_h{h:.2f}',branch='BACKEND',
          data_id=f'analytic_Mie_sphere_r.12_eps3_h{h:.2f}',material_class='Fixed isotropic real epsilon=3; mu=1',
          physical_convention='exp(-i omega t); lambda0=epsilon_b=mu_b=Z0=1',
          mesh=dict(N=N,h=h,cell_volume=h**3,discretization='CM/RR DDA staircase sphere',PML='not applicable: outgoing free-space Green'),
          sources=dict(directions=[[0,0,1],[0,0,1]],polarizations=[[1,0,0],[0,1,0]]),
          receivers=dict(type='far electric vector amplitudes; analytic Mie; separate large-radius near check',source='experiments/qualify_backend.py'),
          frequencies=[1.],noise_model=dict(type='clean analytic reference'),whitening='none; relative physical error',prior='none',
          initialization='fixed known homogeneous sphere; no inverse fit',state_representation='E,H,j,J,p,E_exc at occupied cell centers',
          state_dof=12*N,material_dof=0,material_dof_reason='Known fixed benchmark material; no inverse variables optimized',
          baseline_id='independent analytic Mie sphere',stage='validation',termination_reason='PREDECLARED_TWO_GRID_SPHERE_CONTROL_COMPLETE',
          online_selection_rule='none; both qualification states retained',physical_residual=sr['linear_relative_residual'],
          constitutive_defect=sr['constitutive_relative_residual'],constraint_violation=0,
          cost=dict(LU=sr['LU'],forward_RHS=sr['RHS'],additional_triangular_RHS=0,
            setup=sr['setup'],solve_wall_s=sr['wall_s'],solve_cpu_s=sr['cpu_s'],whole_state_validation_wall_s=None,
            whole_state_validation_reason='Per-state far/near/Mie/I/O overhead not separately timed; enclosing QUALIFICATION process total retains it',
            FFT=0,elliptic_solves=0,Green_actions=None,small_SVD=None,GPU_s=0,currency=None),
          truth_access='Known analytic benchmark, no hidden target',validation=dict(qualification=b['qualification'],metrics={k:v for k,v in row.items() if k!='solve_receipt'},Hcurl='NOT_RUN'),
          status='COMPLETED',scientific_scope='Bounded DDA qualification, no continuum or Hcurl certificate',
          source_references=['experiments/qualify_backend.py','src/vector_backend.py','results/backend_qualification/QUALIFICATION.json'],
          cost_accounting_annotation='Setup and solver are disjoint suboperations, both nested in QUALIFICATION process total; do not add any nested time to that outer total')
        save(r,'Not separately instrumented or inapplicable to this fixed clean benchmark; no inverse candidate selection')
    for name in ('D1','D3'):
        raw=read(f'research/delegated/audit/{name}.json');rr=raw['receipt']
        r=dict(run_id=rr['run_id'],case_id=rr['run_id'],branch=name,seed=1951 if name=='D1' else 829,
          data_id='frozen_dev_calibration_4x4_operators_plus_disclosed_chosen_root',material_class='16 shared independent real scalar contrasts in [.02,20]',
          physical_convention='2D scalar-TM disk-pulse VIE; k=2pi; complex128',mesh=dict(shape=[4,4],N=16),
          sources=dict(incident_operator='G/E0/S/grid_shape from frozen dev_calibration.npz',columns=6),receivers=dict(receiver_operator='S',count=16),
          frequencies=[1.],noise_model=dict(type='clean chosen development root; no random noise observation'),
          whitening='fixed native source weights; D3 includes full unitary-transported source metrics',prior='D1 primary none; targeted common-prior control separate; D3 none',
          initialization='disclosed nonconstant physical root for D1; disclosed uniform m=.1 and public backprojection J for D3',
          state_representation='fixed full auxiliary current C^(16x6), plus shared real material',state_dof=192,material_dof=16,
          baseline_id='frozen R3 alternating structure; exact-current control and full-metric covariance diagnostic',stage='OFFLINE_DIAGNOSTIC',
          termination_reason='BOUNDED_DIAGNOSTIC_COMPLETE',online_selection_rule='none; all predefined controls retained',
          cost=dict(wall_s=rr['wall_seconds'],cpu_s=rr['cpu_seconds'],LU=rr['physical_LU'],forward_RHS=rr['physical_complex_RHS'],
            operation_counts=rr['operation_counts'],setup_included=False,FFT=0,elliptic_solves=0,GPU_s=0,small_SVD=None,Green_actions=None,currency=None),
          truth_access=rr['truth_access'],validation=dict(source=f'research/delegated/audit/{name}.json',scope=raw.get('interpretation',raw.get('stopping'))),
          status='COMPLETED',imaging_claim=False,raw_receipt=rr,
          source_references=['src/dynamics_audit.py',f'research/delegated/audit/{name}.json','audit/D_INPUT_PROVENANCE.json'])
        save(r,'Iteration-map diagnostic: no imaging residual/candidate score; omitted product/factorization subcounts not instrumented')
    setup=d['shared_setup']
    save(dict(run_id='D_SHARED_SETUP',case_id='D_disclosed_root',branch='D_SETUP',stage='setup',status='COMPLETED',
      data_id='disclosed_root_plus_frozen_operator_load',material_class='16 shared real scalars',
      termination_reason='SHARED_ROOT_PREPARED',truth_access='Chosen OFFLINE_DIAGNOSTIC root',
      cost=dict(wall_s=setup['wall_seconds'],cpu_s=setup['cpu_seconds'],LU=setup['physical_LU'],forward_RHS=setup['physical_complex_RHS']),
      source_references=['research/delegated/audit/dynamics_cost.json','audit/D_INPUT_PROVENANCE.json'],
      cost_accounting_annotation='Included once in dynamics_cost outer timer; do not add to it'),
      'Shared setup only; no standalone inverse, candidate selection or gate')
    stops=[('D2_FULL_HISTORY','D2','NOT_RUN','Missing per-iterate full currents, gradients, PR memory, line steps and materials; only scalar histories exist',
            ['BASELINE_AUDIT.md','research/delegated/audit/D2_historical_evidence.json']),
           ('M_PRIMARY','M','NOT_RUN','No qualified compatible Hcurl spaces/mass/weak residual/boundary contract; finite DDA sphere does not supply them',
            ['research/delegated/implementation/M_QUALIFICATION.md']),
           ('T_CERTIFICATE','T','INCONCLUSIVE','No usable certified uniform exception count k<p or complete delta budget; no finite certificate matrix evaluated',
            ['research/delegated/theory/THEORY_QUALIFICATION.md']),
           ('ROM_PRIMARY','ROM','NOT_RUN','No qualified selfadjoint reciprocal colocated temporal observations supplied; frozen plan stops extension',
            ['configs/initial_plan.json','DATA_CONTRACT.md'])]
    for run,branch,status,reason,refs in stops:
        save(dict(run_id=run,case_id=run,branch=branch,status=status,stage='qualification',termination_reason=reason,
          truth_access='NO_NEW_TRUTH_OPENED',source_references=refs,online_selection_rule='No inverse execution',
          cost=dict(LU=0,forward_RHS=0,additional_triangular_RHS=0,known_zero_reason='No new numerical solver launched for this unqualified extension',
                    wall_s=None,cpu_s=None,setup=None,online=None,validation=None,scoring=None,currency=None)),reason)
    aux=[('BACKEND_PLOT','figures/BACKEND_PLOT_RECEIPT.json'),('F_PLOT','figures/F_PLOT_RECEIPT.json'),('F_METADATA_REVIEW','audit/F_METADATA_REVIEW.json')]
    for name,path in aux:
        raw=read(path)
        save(dict(run_id=name,case_id=name,branch='ARTIFACT_REVIEW',status='COMPLETED',stage=raw.get('stage','scoring'),
          termination_reason='SAVED_OUTPUT_PROCESSING_COMPLETE',truth_access='SAVED_DISCLOSED_OUTPUTS_ONLY',cost=raw,
          source_references=[path],scientific_scope='No new forward physical states or imaging predictions'),
          'Saved-output rendering/metadata stage; physical/inverse task fields inapplicable',group='auxiliary')
    components=[dict(name='backend_qualification',source='results/backend_qualification/QUALIFICATION.json',wall_s=b['wall_s'],cpu_s=b['cpu_s']),
                dict(name='F_all_tasks_and_IO',source='results/families/SUMMARY.json',wall_s=f['wall_s'],cpu_s=f['cpu_s']),
                dict(name='D_all_tasks_and_setup',source='research/delegated/audit/dynamics_cost.json',wall_s=d['wall_seconds'],cpu_s=d['cpu_seconds'])]
    scientific=dict(wall_s=sum(x['wall_s']for x in components),cpu_s=sum(x['cpu_s']for x in components))
    for name,path in aux:
        raw=read(path);components.append(dict(name=name,source=path,wall_s=raw['wall_s'],cpu_s=raw['cpu_s']))
        for i,prior in enumerate(raw.get('previous_plot_runs',[])):
            components.append(dict(name=name+f'_prior_render_{i+1}',source=path+f'#previous_plot_runs/{i}',wall_s=prior['wall_s'],cpu_s=prior['cpu_s']))
    totals=dict(status='MEASURED_PROCESS_SUMS_WITH_EXPLICIT_GAPS',components=components,
      scientific_numerical_process_sum=scientific,
      all_measured_process_sum=dict(wall_s=sum(x['wall_s']for x in components),cpu_s=sum(x['cpu_s']for x in components)),
      physical_LU=b['LU']+sum(x['cost']['LU']for x in frows)+d['physical_LU'],
      forward_complex_RHS=b['RHS']+sum(x['cost']['forward_RHS']for x in frows)+d['physical_complex_RHS'],
      additional_retained_LU_triangular_RHS=sum(x['cost']['additional_triangular_RHS']for x in frows),
      current_block_Cholesky=d['operation_counts']['current_block_Cholesky'],root_JVP_Cholesky=d['operation_counts']['root_JVP_Cholesky'],
      current_and_root_JVP_complex_RHS=d['operation_counts']['current_block_complex_RHS']+d['operation_counts']['root_JVP_complex_RHS'],
      small_SVD=None,small_SVD_reason='Local-rank counters and metadata review are recorded, but spectral norms, legality checks and LAPACK subcalls are not exhaustively instrumented; cannot claim aggregate zero or exact total',
      Green_actions=None,Green_actions_reason='Dense Green/receiver products included in timers, not comprehensively counted',
      FFT=0,elliptic_solves=0,GPU_s=0,
      known_nested_timers_not_added=dict(F_per_case_wall_s=sum(x['cost']['wall_s']for x in frows),
        F_per_case_cpu_s=sum(x['cost']['cpu_s']for x in frows),
        F_setup_wall_s=sum(x['cost']['setup']['wall_s']for x in frows),
        backend_solve_wall_s=sum(x['solve_receipt']['wall_s']for x in b['rows']),
        backend_setup_wall_s=sum(x['solve_receipt']['setup']['wall_s']for x in b['rows']),
        D_setup_wall_s=setup['wall_seconds']),
      stage_partition=dict(setup='Nested measured setup suboperations retained; not added to parent process timers',
        online=None,online_reason='No online blind inversion executed',validation=None,
        validation_reason='Physics/validation/root work included in outer numerical timers; exact phase split not instrumented',
        scoring='F metadata review and plot process timers listed separately; other report/scoring overhead not independently measured'),
      unmeasured=dict(startup_failures=d['failed_attempts'],tool_research_and_report_wall_s=None,
        reason='Interpreter imports/startups, shell/tool orchestration, literature retrieval, read-only audit, source preparation and report composition not comprehensively timed',
        currency=None,currency_reason='No billed-service pricing meter or amortized local-hardware price evaluated'),
      counting_rule='Each enclosing process timer added once. Per-case/solver/setup/current timers and standardized receipt copies are references, not additional runs. This is process-time sum, not task elapsed time.')
    write(ROOT/'audit/FINAL_COST_TOTALS.json',totals)
    write(OUT/'INDEX.json',dict(schema_fields=FIELDS,receipts=index,
      null_policy='Every null field has a dotted-path explanation in null_reasons; no hashes computed',
      aggregate_cost_source='audit/FINAL_COST_TOTALS.json',raw_results_preserved=True))
    write(ROOT/'audit/D_INPUT_PROVENANCE.json',dict(frozen_git_base=BASE,
      verified_in_git_tree=[BASE+':research/research_fbg_r3_experiments/src/physics.py',BASE+':research/research_fbg_r3_experiments/data/dev_calibration.npz'],
      in_study_dependency_snapshot='research/delegated/audit/frozen_diversity_physics.py',
      dependency_original_path='/Volumes/migodam\'s-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_diversity/src/physics.py',
      in_study_adapter_snapshot='research/delegated/audit/frozen_r3_physics.py',
      adapter_change='Only import route changed to included frozen diversity module alias; numerical code unchanged',
      in_study_public_operators='data_public/dynamics_operators_4x4.npz',
      copied_operator_fields=['G','E0','S','grid_shape'],
      unexported_fields='No Y or material/truth arrays opened or copied',
      execution='Standalone study: python src/dynamics_audit.py. No external historical study modules/data required; dependencies NumPy/SciPy/threadpoolctl',
      hashes=None,hashes_reason='Forbidden by user; Git tree source existence inspected without recomputing hashes'))
    print(json.dumps(dict(receipts=len(index),totals=totals['all_measured_process_sum'],LU=totals['physical_LU'],RHS=totals['forward_complex_RHS'])))

if __name__=='__main__':main()
