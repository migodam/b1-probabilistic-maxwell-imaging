"""Read actual local measured/benchmark packages; no inverse or truth reads."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.physics_inputs import prepare_fresnel,prepare_austria
from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,UnreachableDataBall
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump


def main():
    start,cpu=time.perf_counter(),time.process_time();rows=[]
    for case in ('dielTM_dec8f','twodielTM_8f'):
        d=prepare_fresnel(case);m=d['model'];Y=d['Y'];t=time.process_time()
        p=DenseFieldCurrentProjectors(d['G'],d['E'],d['S'],Y,d['delta'],
            components=1,scale=6,lower=0,upper=20,grouped_ports=True)
        initial=np.full(m.n,.1);st=m.state(initial);f,g=m.data_value_gradient(initial,Y)
        X,Q=p.affine_ball(st['j']/6,st['Q'])
        error=np.sqrt(sum(np.linalg.norm(6*S@X[:,j]-Y[:,j])**2 for j,S in enumerate(d['S'])))
        try:
            tiny=DenseFieldCurrentProjectors(d['G'],d['E'],d['S'],Y,0,
                components=1,scale=6,lower=0,upper=20,grouped_ports=True)
            zero='AFFINE_REACHABLE_WITHIN_NUMERICAL_RANK; not physical-material feasible'
        except UnreachableDataBall as exc:zero=str(exc)
        cal=d['metadata']['calibration']
        rows.append(dict(case=case,source_path=d['metadata']['source_path'],
            G_shape=list(d['G'].shape),E_shape=list(d['E'].shape),
            source_port_shapes=[list(S.shape) for S in d['S']],
            training_channels=Y.size,withheld_channels=d['metadata']['withheld_channels'],
            calibration_validation_median=float(np.median([r['validation_relative'] for r in cal])),
            delta=d['delta'],delta_relative=.15,affine_rank_floor=p.minimum_reachable_residual,
            affine_projected_residual=float(error),zero_ball_diagnostic=zero,
            material_physical_initial_relative_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),
            original_initial_Phi=f,gradient_norm=float(np.linalg.norm(g)),
            physical_costs=dict(m.cost),projector_setup_CPU=p.setup_cpu_seconds,
            verification_CPU=time.process_time()-t,
            qualification='actual measured source-dependent ports and empty-field source fit; unknown receive-horn response and noise covariance; discrepancy ball only'))
    try:
        a=prepare_austria()
        austria=dict(status='PUBLIC_LOADED_NO_TRUTH',source_path=a['metadata']['source_path'],
            Y_shape=list(a['Y'].shape),G_shape=list(a['G'].shape),metadata=a['metadata'])
    except FileNotFoundError as exc:austria=dict(status='PUBLIC_PACKAGE_ABSENT',error=str(exc))
    dump(HERE/'EXTERNAL_INPUT_RECEIPT.json',dict(status='COMPLETE',Fresnel=rows,
        Austria=austria,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
        original_raw_copied=False,inversion_run=False,material_truth_access=False))
    print('Actual Fresnel adapters checked; receipt saved',flush=True)
if __name__=='__main__':main()
