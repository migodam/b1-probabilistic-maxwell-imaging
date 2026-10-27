"""Inspect pinned local assets and write proposed contracts; performs no PDE solve."""
from pathlib import Path
import importlib.metadata as md
import json, os, platform, subprocess, sys, time

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
W='958a478739d9553269316ded8125fda9df225c0c'
SLAM='52351aee6ddd41b09851fe5e52d1a8861891ef13'

def write(name,obj):
    (ROOT/name).parent.mkdir(parents=True,exist_ok=True)
    (ROOT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def main():
    t0=time.process_time()
    auth=ROOT/'AUTHORIZATION.json'
    if auth.exists() and json.loads(auth.read_text()).get('new_physics_authorized'):
        raise RuntimeError('Preflight may not overwrite an already granted/frozen physical campaign')
    vendor=ROOT/'vendor'/'wslam';vendor.mkdir(parents=True,exist_ok=True)
    names=('__init__','fem_backend','forward','physics','world','geometry','accounting','factors','information')
    files=[]
    for name in names:
        path=f'research/research_w_slam/src/wslam/{name}.py'
        data=subprocess.run(['git','show',f'{SLAM}:{path}'],cwd=REPO,capture_output=True,check=True).stdout
        (vendor/f'{name}.py').write_bytes(data)
        files.append(dict(commit=SLAM,path=path,bytes=len(data),operation='git show; no hash computation'))
    cp=dict(schema='BGGA_STATIC_V1',status='FROZEN_DESIGN_PENDING_NEW_PHYSICS_AUTHORIZATION',
        source_commits=dict(W=W,WSLAM=SLAM),static_wall=True,static_target=True,
        acquisition_index='frequency/source/receiver, never motion time',
        physics=dict(polarization='2D scalar TM Ez',time='exp(-i omega t)',mu_r=1.,f0_hz=1e9,
            relative_frequencies=[1.0],dispersion='epsinf+i ell/f',dtype='complex128/float64',
            threads=1,noise_std_each_real_imag=.01,complex_noise_variance=.0002,
            packing='complex frequency/source/receiver C order flattened, then [all Re, all Im]; identical y/J ordering',
            whitening='divide each real coordinate by .01; no sqrt(channel count)',
            source='unit i H0^(1)(2pi f r)/4; no fitted gain',
            tx=[[-1.5,0],[1.5,0]],rx=dict(count=16,x_min=-2,x_max=2,y=0),
            public_wall_roi=[-1.5,1.5,.5,1.25],target_material_integration_roi=[-1.2,1.2,1.4,2.8],
            pml=dict(physical_box=[-2.75,2.75,-.8,3.5],thickness=1.,strength=4.,power=3),
            field='scalar H1 P3 FEM, not 3D H(curl)',max_dofs=60000,
            max_factor_bytes=1500000000),
        fixed_space=dict(h=.17,wall_macro_edge=.125,wall_subdivision=2,target_step=.125,
            circle_polygon_n=384,cut_cell_integration='P3 polynomial on fixed triangles; physical circle polygon intersections',
            no_target_in_mesh_key=True,independent='pinned interface-fitted P3, h=.17, n=384, wall refinement=3'),
        statistical_scope='new replay uses same explicit finite counting prior; no continuous posterior claim',
        stage_order=['preflight/code/zero-PDE verification','G0 new physical qualification','G1','G2'],
        proposed_not_authorized=dict(G0=dict(rhs=4000,full_fe=500,cpu_seconds=7200),
            G1=dict(rhs=20000,full_fe=5000,cpu_seconds=28800,replay_cpu_seconds=3600),
            G2=dict(rhs=80000,full_fe=12000,cpu_seconds=86400)),
        resource_policy='Old W/WSLAM caps are campaign-specific and cannot be transferred. Authorization must be explicit.',
        exclusions=['locked test access','NN','3D','new A17 modules','active sensing campaign','hash checks'])
    write('CONTRACT.json',cp)
    write('AUTHORIZATION.json',dict(schema='BGGA_AUTH_V1',new_physics_authorized=False,
        authorized_caps={},reason='User explicitly labels budgets as suggestions; no BGGA-specific resource grant located.',
        allowed=['read-only audit','implementation','zero-PDE tests','cost projection'],
        old_budgets_not_transferable=True))
    write('PRIOR.json',dict(prior_id='FINITE_PRIOR_ENUM4_V1',measure='finite discrete counting prior, exactly 2950 atoms',
        support='14 legal 2x2 occupancies; nonempty 4-connected, no diagonal-only checkerboard; includes empty',
        nonempty_wall_groups=117,empty_wall_groups=1,pose_atoms_per_group=25,
        material=dict(epsinf=[1.5,3.75,6],ell=[0,.125,.25],weights=[1/6,4/6,1/6],shared_across_occupied_cells=True),
        occupancy_p=.65,pose='pinned 5x5 nodes and counting masses; not continuous quadrature certification',
        target=dict(radius=.12,contrast=3,ell=.02,intrinsic_known=True),gamma_known=True,
        no_holes_in_2x2=True,continuous_A0='different prior; never label finite table as its quadrature certificate'))
    write('TASKS.json',dict(T_obs=dict(kind='fixed full nonlinear total response',sources=[0,1],receivers=list(range(16)),frequency=1.),
        T_future=dict(kind='same target disk at 25 prior pose nodes, plus G0 predeclared radius/contrast tests',
            norm='max_j whitened response norm, not concatenation'),
        posterior_task=dict(action='posterior mean position',loss='squared distance in lambda0^2',
            events=['pose_x>=0','pose_y>=2.1'],mode_partition='joint (wall support ID, x sign, y>=2.1); support alone separately'),
        history_stress=dict(identity='oracle conditional discrete witness, separate from new-prior draws',
            joint_distance=.7987161019689071,cal_distance=3.5605119112114565,
            target_x_displacement=.1,probe_max=4.81819338103861,probe_argmax=6,
            prohibited='not a calibration omega(1) bound or posterior mass')))
    write('GATES.json',dict(status='new engineering thresholds from user suggestion, frozen for proposed run; no result yet',
        G0=dict(relative_same_discrete_identity=1e-9,linear_backward_residual=1e-10,
            independent_whitened_max=.1,small_margin_rule='independent numerical error must be below decision margin; .1 is not a posterior certificate'),
        G1=dict(joint_TV=.02,coordinate_CDF_sup=.02,event_probability=.02,important_joint_mode_min_mass=.05,
            important_joint_mode_error=.02,position_Bayes_regret='<=.05*R_ref+(.005 lambda0)^2',
            exact_bound_violation='FAIL; isolate bug/assumptions',empirical_bound='EMPIRICAL_INTERVAL only'),
        G2=dict(mean_end_to_end_cpu_reduction=.20,paired_parent_interval='must support positive saving',
            physical_reduction='FE or RHS lower without shifting cost to other physical work or memory'),
        no_automatic_authorization=True))
    write('DATA_SPLIT.json',dict(development=dict(parents=12,noise_per_parent=2,
        source='all 24 W-SLAM cases previously exposed; replay/development only'),
        stress='12 G0 parents constructed before physics, 3 intended examples/class, observations may refute classification',
        locked=dict(parents_min=30,noise_per_parent=2,status='NOT_CREATED_NOT_OPENED',
            conditional_60='requires predeclared variance-only rule and separate resource grant'),
        truth_access='independent scorer only; online API has prior, data, own cache and physics query'))
    env={p:md.version(p) for p in ['numpy','scipy','matplotlib','scikit-fem','triangle']}
    write('source_review/PREFLIGHT.json',dict(source_files=files,python=sys.version,executable=sys.executable,
        platform=platform.platform(),architecture=platform.machine(),hardware='Apple M4, 24 GB; CPU only',
        dependencies=env,head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        process_cpu_seconds=time.process_time()-t0,new_fe=0,new_rhs=0,
        source_pdf=dict(path='/Volumes/migodam\'s-external-brain/archive/user-home/downloads/Research_W_BGGA_Chinese_Report.pdf',pages=32),
        no_sha256_performed=True))
    print('Pinned sources and proposed static contracts prepared; new physics authorization remains false.')

if __name__=='__main__':main()
