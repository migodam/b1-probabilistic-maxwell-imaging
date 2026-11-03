"""Fixed two-cell vector-DDA residual-spectrum example using the R1 backend."""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import json
import sys
import time
import traceback
from datetime import datetime,timezone

WALL,CPU=time.perf_counter(),time.process_time()
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
cfg=json.loads((ROOT/'configs/vector_pair_frozen.json').read_text())
out=ROOT/'RESULTS/vector_pair'
out.mkdir(parents=True,exist_ok=False)
def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False))
def cj(z):return dict(real=float(z.real),imag=float(z.imag))
def cm(z):return dict(real=z.real.tolist(),imag=z.imag.tolist())
dump(out/'FROZEN_CONFIG.json',cfg)
receipt=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
    configuration=cfg,executable=sys.executable,truth_access='NONE',scoring='NOT_RUN',
    inversion='NOT_RUN',derivative_validation='NOT_RUN',new_formula_tests='NOT_RUN',
    seed=None,failed_attempts=[],imaging_status='NOT_AN_IMAGING_EXPERIMENT')
dump(out/'RECEIPT.json',receipt)
model=state=None
try:
    import numpy as np
    from threadpoolctl import threadpool_limits,threadpool_info
    from research.research_max_struct_r1.src.vector_backend import TensorDDA
    with threadpool_limits(limits=1):
        receipt['threadpools']=threadpool_info()
        points=np.asarray(cfg['centers'],float)
        k=float(cfg['k']);h=float(cfg['h']);volume=float(cfg['volume'])
        start,ct=time.perf_counter(),time.process_time()
        model=TensorDDA(points,volume,k=k)
        receipt['operator_setup_wall_seconds']=time.perf_counter()-start
        receipt['operator_setup_cpu_seconds']=time.process_time()-ct
        epsilon=np.broadcast_to(cfg['epsilon']*np.eye(3),(2,3,3)).copy()
        directions=np.asarray(cfg['directions'],float)
        polarizations=np.asarray(cfg['polarizations'],float)
        # One actual backend call: shared factorization, both legal waves as RHS.
        state=model.solve(epsilon,directions,polarizations,receipt_path=out/'SOLVE_RECEIPT.json')
        start,ct=time.perf_counter(),time.process_time()
        receivers=np.asarray(cfg['receivers'],float)
        near=state.near(receivers)
        receipt['external_near_wall_seconds']=time.perf_counter()-start
        receipt['external_near_cpu_seconds']=time.process_time()-ct
        np.savez_compressed(out/'FIELDS.npz',centers=points,epsilon=epsilon,
            directions=directions,polarizations=polarizations,incident_E=state.incident,
            E=state.E,j=state.j,J=state.J,p=state.p,
            receivers=receivers,E_near=near['E'],E_near_scattered=near['E_scattered'])
        # G acts on density j. The backend's Goff property alone is the
        # unscaled dipole kernel; it must be multiplied by cell volume here.
        goff=model.G-model.self_density*np.eye(model.n)
        block=model.G[:3,3:6]
        self_block=model.G[:3,:3]
        m=float(cfg['epsilon']-1)
        rows=[];operators={}
        for beta in cfg['beta']:
            a=m/(1+(float(beta)-model.self_density)*m)
            operators[f'T_beta_{beta}']=np.eye(model.n)+a*goff
            for component,label in enumerate(('x_longitudinal','y_transverse','z_transverse')):
                e=np.eye(3)[component]
                for spatial,sign in (('coherent',1),('antisymmetric',-1)):
                    mode=np.r_[e,sign*e]/np.sqrt(2)
                    eigenvalue=sign*block[component,component]
                    action=a*eigenvalue
                    rows.append(dict(beta=float(beta),polarization=label,spatial_mode=spatial,
                        eigenvalue=cj(eigenvalue),a=cj(a),a_lambda=cj(action),
                        interference_term_2Re_a_lambda=float(2*action.real),
                        quadratic_term_abs_a_lambda_squared=float(abs(action)**2),
                        residual_weight_abs_1_plus_a_lambda_squared=float(abs(1+action)**2),
                        mode_vector=mode.tolist(),
                        physically_excited_by_fixed_waves=bool(spatial=='coherent' and component in (0,2))))
        static_block=volume/(4*np.pi*h**3)*np.diag([2.,-1.,-1.])
        report=dict(status='COMPLETED_FIXED_VECTOR_PAIR_SPECTRUM',m_susceptibility=m,
            convention='exp(-i omega t), lambda0=Z0=1; j=(epsilon-I)E is polarization density; J=-ik*j',
            G_density_offdiag_block_01=cm(block),G_density_self_block_00=cm(self_block),
            backend_self_density=cj(model.self_density),
            raw_dipole_offdiag_block_01=cm(model.Goff[:3,3:6]),
            volume_scaling='T uses density Goff=volume*backend.Goff, not the unscaled backend.Goff property',
            static_sign_explanation_only=dict(block=static_block.tolist(),longitudinal_sign='positive',
                transverse_sign='negative',used_in_actual_fields_or_weights=False),
            spatial_polarization_modes=rows,
            degrees_of_freedom=dict(electric_complex_dof=6,known_isotropic_scalar_material_cells=2,
                material_unknowns_fitted=0,mu_fixed=1),
            electric_sizes=dict(h_over_lambda=h,kh=k*h,inside_material_kh=np.sqrt(cfg['epsilon'])*k*h,
                occupied_extent_over_lambda=[2*h,h,h],total_volume_over_lambda_cubed=2*volume),
            scope='Coarse fixed two-cell CM+RR point-dipole model; not an analytic sphere, continuous PDE qualification, inverse recovery or 3D imaging PASS.',
            historical_backend_qualification='R1 results/backend_qualification/QUALIFICATION.json: bounded DDA qualification only; continuum Hcurl NOT_RUN; epsilon3 sphere radius.12 at h.04/.03 had far-vector errors .104679/.036554. These historical sphere errors do not certify the present two-cell geometry.',
            comparison_to_2D='Vector dyadic CM+RR DDA has longitudinal/transverse polarization structure; it is a different backend from 2D scalar-TM disk-pulse VIE. No common-backend or imaging comparison.')
        dump(out/'SPECTRUM.json',report)
        np.savez_compressed(out/'OPERATORS_AND_MODES.npz',G=model.G,
            Goff_density=goff,Goff_raw=model.Goff,Gself_block=self_block,
            offdiag_block=block,**operators)
        lines=['# 固定两体素 vector Green 残差谱例','','使用R1原TensorDDA：lambda=1、h=.04、epsilon=3、mu=1；两相邻立方体素中心沿x。电场自由度6，已知各向同性材料标量2；没有反演。两个+y传播平面波分别x/z极化，一次6×6 LU共享2 RHS。外部电场仅取沿三轴±.5的六个点。','',
            f"读取的self={model.self_density.real:.12g}{model.self_density.imag:+.12g}i。密度Green offdiag纵向x={block[0,0].real:.12g}{block[0,0].imag:+.12g}i；横向y/z={block[1,1].real:.12g}{block[1,1].imag:+.12g}i。实际retarded Green的纵/横实部分别正/负；静电diag(2,-1,-1)仅作符号解释，未用于场解或权重。",'',
            'm=epsilon−1=2；a=m/[1+(beta−gself)m]；T=I+a Goff，其中Goff是乘过体积的密度算子。以下权重均为实际复Green的 |1+a lambda|²。','',
            '| beta | 极化 | coherent | antisymmetric |','|---:|---|---:|---:|']
        for beta in cfg['beta']:
            for label in ('x_longitudinal','y_transverse','z_transverse'):
                selected=[row for row in rows if row['beta']==beta and row['polarization']==label]
                lines.append(f"| {beta} | {label} | {selected[0]['residual_weight_abs_1_plus_a_lambda_squared']:.12g} | {selected[1]['residual_weight_abs_1_plus_a_lambda_squared']:.12g} |")
        lines.extend(['','两固定平面波只激励coherent x与coherent z；其余mode是完整残差空间的代数方向。纵向coherent与横向coherent的权重方向不同，不能据2D scalar的同相解释直接统一描述该vector算子。','',
            'kh='+format(k*h,'.8g')+'，材料内部kh='+format(np.sqrt(cfg['epsilon'])*k*h,'.8g')+'。这是粗离散两cell模型，不是解析球或连续PDE。继承R1 bounded-DDA适用界；历史Mie球结果不认证这两个体素。没有新增恒等式/导数测试，没有GT评分，不宣布3D imaging PASS。科学解释与新颖性由父线程审查。'])
        (out/'README.md').write_text('\n'.join(lines)+'\n')
        receipt.update(status='COMPLETED_FIXED_VECTOR_PAIR_MECHANISM',solve_receipt=state.receipt,
            LU=state.receipt['LU'],actual_complex_RHS=state.receipt['RHS'],solve_calls=1,
            external_electric_operator_builds=1,external_electric_matrix_apply_RHS=2,
            external_magnetic_operator_builds=1,external_magnetic_matrix_apply_RHS=2,
            near_API_scope='Existing state.near computes E and H together; only E requested/output. No additional linear solve.',
            dof=6,known_material_scalar_cells=2,h=h,kh=k*h)
except Exception as exc:
    receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
    receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
finally:
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds=time.perf_counter()-WALL,cpu_seconds=time.process_time()-CPU,
        backend_setup_receipt=model.setup_receipt if model is not None else None,
        cost_scope='Process imports,existing backend operator assembly,one shared6x6 LU/two RHS,existing near E/H operator applications,six-mode formula table at exactly two beta values,all artifacts; final receipt write excluded.')
    dump(out/'RECEIPT.json',receipt)
print(json.dumps(dict(status=receipt['status'],LU=receipt.get('LU'),RHS=receipt.get('actual_complex_RHS'),wall_seconds=receipt['wall_seconds'],cpu_seconds=receipt['cpu_seconds'])),flush=True)
if receipt['status']=='FAILED':raise SystemExit(1)
