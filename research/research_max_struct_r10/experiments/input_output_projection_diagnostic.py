"""Public-operator/state diagnosis of the PGP versus GP structural correction.

No truth and no full-wave solve. Reuse the already saved strong right basis;
quantify the actual omitted output, not an identity error-floor check.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import json,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve,svdvals

ROOT=Path(__file__).resolve().parents[1]

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'delegated/input_output_projection';out.mkdir(parents=True,exist_ok=False)
    with np.load(ROOT/'data_public/twod/operators.npz') as z:G=z['G']
    first=ROOT/'RESULTS/twod/r10_twod_0001__matrix_green19__noise_ball_99'
    with np.load(first/'DENOMINATOR_BASIS.npz') as z:U=z['U']
    B=G@U;compressed=U@(U.conj().T@B);lost=B-compressed
    n=lambda x:float(np.linalg.norm(x))
    result=dict(scope='EXPOSED_POSTHOC_PUBLIC_STATE_DIAGNOSTIC',truth_access='NONE',
        rank=U.shape[1],basis_source=str((first/'DENOMINATOR_BASIS.npz').relative_to(ROOT)),
        strong_input_full_output_Frobenius=n(B),discarded_left_output_Frobenius=n(lost),
        discarded_left_output_relative=n(lost)/n(B),
        strong_output_fraction_retained_by_left_projection=1-(n(lost)/n(B))**2,
        rows=[])
    for suffix in ('0001','0004','0022'):
        folder=ROOT/'RESULTS/twod'/('r10_twod_'+suffix+'__matrix_green19__noise_ball_99')
        with np.load(folder/'stage_03_joint.npz') as z:
            m=z['m'];J=z['J'];Q=z['Q']
        beta=.2;a=1+beta*m;z=m/a;e=J-m[:,None]*Q;residuals={};spectra={}
        for name,basis in (('PGP',compressed),('GP',B)):
            F=z[:,None]*basis;T=np.eye(U.shape[1])-U.conj().T@F
            factor=lu_factor(T,check_finite=False);b=e/a[:,None]
            residuals[name]=beta*(b+F@lu_solve(factor,U.conj().T@b,check_finite=False))
            s=svdvals(T,check_finite=False);spectra[name]=dict(min_singular=float(s[-1]),condition_2=float(s[0]/s[-1]))
        result['rows'].append(dict(case_id='r10_twod_'+suffix,state='stage_03_joint',beta=beta,
            PGP_residual_norm=n(residuals['PGP']),GP_residual_norm=n(residuals['GP']),
            relative_residual_change=n(residuals['GP']-residuals['PGP'])/n(residuals['PGP']),
            small_systems=spectra))
    result.update(physical_LU=0,physical_RHS=0,auxiliary_small_LU=6,auxiliary_small_dimension=19,
        auxiliary_small_RHS=72,small_SVD=6,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        limitation='Smaller Green approximation error or residual difference is not a material recovery guarantee. PGP dissipation proof does not apply automatically to GP.')
    (out/'RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
