import sys,json,time
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
from research_w.geometry import antennas,probes
tx,rx=antennas()
config=dict(study='Research W P0',frozen_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
 scope=['P0-A','P0-B'],not_authorized=['P1 computation','NN training','3D'],
 backend='A2 outgoing Sommerfeld Green VIE; exact direct multiple scattering, square-cell integration',
 domain=[[-2,2],[0,3]],tx=tx.tolist(),rx=rx.tolist(),
 frequencies=dict(wide=[.8,.9,1.,1.1,1.2],narrow=[.9975,.99875,1.,1.00125,1.0025]),
 source_normalization='Unit line-source Green, constant multiplier 1 at every frequency; known calibrated source spectrum',
 noise=dict(each_real_quadrature_std=.01,covariance='0.0001 I; one fixed absolute field unit',realification='source-major [Re(receiver),Im(receiver)] per frequency',noise_seed=2026102203),
 materials=dict(real_epsilon=[1,9],conductivity_over_omega0=[0,.15],dispersion='epsilon=epsilon_infinity+i conductivity/f; same two real spatial fields at all frequencies',maximum_occupied_area=3.0,air_phase='exact epsilon=1',minimum_feature=.10),
 truth_seeds=[2026102201,2026102202],
 query=[[-1.2,1.2],[1.6,2.8]],probes=probes(),construction_probe_indices=[8],validation_probe_indices=[i for i in range(24) if i!=8],
 inverse_representation=dict(grid=[32,48],coordinates='full-domain continuous passive cells plus thresholded multiscale occupancy proposals; no known wall mask or fixed cavity count'),
 forward_grid_ladder=[32,64,96,128,192],forward_grid_choice='Before search, finest affordable adjacent grids meeting <=0.1 whitened full-data mesh differences; no claim of continuum certification',
 quadrature_order=4,coefficient_subsamples=4,
 benchmark=dict(nx=32,frequency=1.,world='homogeneous dielectric disk qualification only; no wall search outcome'),
 budgets=dict(complex_rhs=8000,material_frequency_factorizations=2000,charged_cpu_seconds=7200,allocation=dict(search=.4,endpoint=.3,validation=.2,repair=.1)),
 search=dict(methods=['FULL_J','RANDOM_MATCHED','MULTISTART_FWI','TASK_AWARE'],seed=2026102204,shared_initial_derivatives=True,per_context_rhs=800,per_method_rhs=160,contexts=4,corrector='same regularized full-J trust-region corrector, passive and occupancy constraints',maximum_corrector_trials=4),
 nuisance=dict(primary='calibrated gamma=1',stress_gain_amplitude=[.99,1.01],stress_shared_phase_rad=[-.01,.01],per_frequency_free_gain=False),
 task=dict(calibration_max_mean_distance=1.,target_min_response_distance=3.,joint_max_mean_distance=1.,joint_min_centroid_difference=.1,centroid='positive total dielectric contrast within fixed T; mass within +/-20% of reference target'),
 numerics=dict(mesh_whitened_tolerance=.1,backward_residual_tolerance=1e-10,margin='mesh differences are consistency checks, not rigorous PDE bounds; no GO at a threshold inside measured numerical uncertainty',max_active_dense_cells=4500,peak_memory_target_gb=1.5),
 perturbations=dict(epsilon_infinity_relative=[-.002,.002],target_position_lambda=[-.005,.005]),
 stops=dict(analytic_failure='BLOCKED_PHYSICS',qualification_failure='UNQUALIFIED_SOLVER; do not run search',no_witness='NOT_FOUND_WITHIN_BUDGET, not uniqueness',P1='NOT_RUN_NOT_AUTHORIZED'))
dest=root/'w_config.json'
if dest.exists():raise RuntimeError('Frozen config already exists; inspect instead of overwriting')
dest.write_text(json.dumps(config,indent=2,ensure_ascii=False));print(dest)
