"""Verified completion for a uniformly coercive material form without an observation."""
from pathlib import Path
from time import process_time,perf_counter
import json,sys

import numpy as np
from scipy.sparse.linalg import LinearOperator,cg
from flint import acb,acb_mat,arb,ctx
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[2]
V2=STUDY/'ip_submission_v2'
sys.path.insert(0,str(V2))
from energy_port_control import eigen_interval
sys.path.insert(0,str(STUDY/'submission/single'))
from local_moment import ExactGeometry,NumericGeometry,blocks,hi,lo,ball,star
from general_fourier_verify import scale_rows,conjugate_values
from coherent_dual import expand_block
from core.certificate import frobenius_squared


def main():
    ctx.prec,ctx.threads=112,1
    wall,cpu=perf_counter(),process_time()
    with np.load(STUDY/'synthetic/public/fine_geometry.npz') as arr:
        S,B,symbol=(arr[k] for k in ('PW360_10x20__S','PW360_10x20__B','spectrum'))
    idx=blocks(4);n=len(B);num=NumericGeometry(S,B,symbol,idx);exact=ExactGeometry(S,B,symbol,idx)
    transfer=json.loads((HERE.parent/'receipts/radiation_transfer.json').read_text()); gamma=arb(transfer['radiation_coefficient']); eta=arb(transfer['isotropic_remainder_upper']);
    profiles=json.loads((STUDY/'submission/single/coherent_general_profiles.json').read_text());data={}
    # A single valid pre-existing geometry profile suffices for all materials and all a.
    z=profiles[2];mu=arb(z['mu']);w=expand_block(z['sector']);t=expand_block(z['optical']);
    strong_path=HERE.parent/'receipts/energy_bound.json'
    shape_key='optfresh_01__pw360_10x20__n02';delta=arb(1)/1024
    if strong_path.exists():
        row=json.loads(strong_path.read_text())['cases'][shape_key]
        if row['gap_replay'] and row['profile_index']==2:
            delta=arb(row['unshifted_scaled_K_gap_lower'])
    floor=(delta-eta/mu).lower();shift=arb(float(floor)) * 15/16
    if not floor>shift>0: raise ArithmeticError('The shifted whole-box form has no certified floor.')
    wz=[arb(v) for v in w];za=[acb(arb(.75)*ww,arb(tt)+gamma/mu)/2 for ww,tt in zip(wz,t,strict=True)]
    zn=np.asarray([complex(v) for v in za]);rhs=zn[:,None]*B
    def mv(v):return (w-float(shift))*v-zn*num.green(v)-num.green(zn.conj()*v,num.symbol.conj())
    op=LinearOperator((n,n),matvec=mv,dtype=complex);cols=[];its=[]
    for j in range(B.shape[1]):
        calls=[];val,info=cg(op,rhs[:,j],rtol=2e-11,atol=0.,maxiter=950,callback=lambda _:calls.append(1));its.append(len(calls));cols.append(val)
        if info: raise ArithmeticError(f'CG proposal stop {info}.')
    X=ball(np.column_stack(cols));L=scale_rows(exact.b,za)
    kx=(scale_rows(X,[v-shift for v in wz])-scale_rows(exact.green(X),za)-exact.green(scale_rows(X,conjugate_values(za)),adjoint=True))
    defect=frobenius_squared(kx-L);lin=star(L)*X+star(X)*L-star(X)*kx
    _,bound=eigen_interval(lin);correction=(defect/(floor-shift)).upper();energy=((bound+correction)/shift).upper()
    data={'scope':'all admissible real materials 0<=m<=.75 and unit coherent a, no Y and no noise radius, 96x96 finite TM','profile_index':2,'optical_addition_normalized':hi(gamma/mu),'common_K_floor':lo(delta),'whole_box_form_floor':lo(floor),'shift':float(shift),'stationary_cg_calls':its,'fft_verified_rhs':exact.green_columns,'defect_Fro_upper':hi(defect.sqrt()),'source_correction_upper':hi(correction),'energy_upper':hi(energy),'wall_seconds':perf_counter()-wall,'cpu_seconds':process_time()-cpu,'precision_bits':ctx.prec};(HERE.parent/'receipts/unconditional_completion.json').write_text(json.dumps(data,indent=2)+'\n');print(data)

if __name__=='__main__':
    with threadpool_limits(1):main()
