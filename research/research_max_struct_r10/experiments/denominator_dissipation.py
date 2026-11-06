"""New auxiliary-denominator qualification; public compressed Green only.

Do not interpret floating near-zero eigenvalues as a rigorous positive margin.
No material truth, inverse solve, parameter selection or new forward.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import json,time
import numpy as np
from scipy.linalg import eigvalsh,norm
ROOT=Path(__file__).resolve().parents[1]


def main():
    wall,cpu=time.perf_counter(),time.process_time()
    output=ROOT/'delegated/denominator_dissipation';output.mkdir(exist_ok=False)
    with np.load(ROOT/'data_public/twod/operators.npz',allow_pickle=False) as z:G=z['G']
    rows=[];arrays={}
    for basis in ('green_strong','green_equal_rank','twofold'):
        path=ROOT/'RESULTS/denominator_dev'/f'obs_b20d5a55677d4eca__{basis}__galerkin'/'DENOMINATOR_BASIS.npz'
        with np.load(path,allow_pickle=False) as z:U=z['U']
        H=U.conj().T@G@U
        ih=(H-H.conj().T)/(2j);rh=(H+H.conj().T)/2
        ei=eigvalsh(ih);er=eigvalsh(rh)
        arrays.update({basis+'_H':H,basis+'_imag_spectrum':ei,basis+'_real_spectrum':er})
        rows.append(dict(basis=basis,rank=U.shape[1],ImH_min=float(ei[0]),ImH_max=float(ei[-1]),
            ReH_min=float(er[0]),ReH_max=float(er[-1]),H_norm=float(norm(H,2)),
            real_coercivity_margin_float_only={str(beta):float(beta+1/20-er[-1]) for beta in (1.,.2)},
            certificate='NONE; floating radiative minima at rounding scale; no interval bound',
            source=str(path.relative_to(ROOT))))
    np.savez_compressed(output/'COMPRESSED_SPECTRA.npz',**arrays)
    result=dict(status='COMPLETED_PUBLIC_OPERATOR_DIAGNOSTIC_NOT_CERTIFICATE',rows=rows,
        truth_access='NONE',physical_LU=0,physical_RHS=0,new_inverse_runs=0,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        initial_exploration_CPU_seconds=.0072470000000000034,
        initial_exploration_wall_seconds=None,
        limitations='Positive/negative eigenvalues at about1e-16 do not establish strict dissipation or its failure. No formula/rank/budget changed.')
    (output/'RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
