"""A data-independent current-energy bound using a verified radiation-to-port inequality."""
from pathlib import Path
import json
import sys
from time import process_time, perf_counter

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import svd
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
V2 = HERE.parents[1]
STUDY = V2.parent
sys.path.insert(0, str(V2))
from energy_port_control import eigen_interval
sys.path.insert(0, str(STUDY / 'submission/single'))
from local_moment import SIDE, N, hi, lo
from coherent_dual import expand_block
from core.certificate import ball, frobenius_squared, star


def inverse_symbol(symbol):
    rows = [acb.dft([acb(float(z.real),float(z.imag)) for z in row],inverse=True) for row in symbol]
    return [acb.dft(column,inverse=True) for column in zip(*rows,strict=True)]  # kernel[x][y]


def main():
    start, cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = 112, 1
    with np.load(STUDY/'synthetic/public/fine_geometry.npz') as data:
        symbol, s, b = (data[k] for k in ('spectrum','PW360_10x20__S','PW360_10x20__B'))
    k = 2*np.pi*4e8/3e8*2/SIDE
    q = 64
    angle = np.arange(q)*2*np.pi/q
    kx, ky = k*np.cos(angle), k*np.sin(angle)
    c = arb(.007606209054719168)
    xs = range(1-SIDE,SIDE)
    phase = lambda theta: [[acb(0,arb(float(a))*i).exp() for a in theta] for i in xs]
    fx,fy = phase(kx),phase(ky)
    equiv = acb_mat(fx)*acb_mat(fy).transpose()* (c/q)
    ker = inverse_symbol(symbol)
    pad=len(ker)
    worst = arb(0)
    for i,x in enumerate(xs):
        for j,y in enumerate(xs):
            rij=(ker[x%pad][y%pad]-ker[-x%pad][-y%pad].conjugate())/acb(0,2)
            err=abs(rij-equiv[i,j]).upper()
            if err>worst: worst=err
    rho=(N*worst).upper()
    xx,yy=np.meshgrid(np.arange(SIDE),np.arange(SIDE))
    U=np.sqrt(float(c)/q)*np.exp(1j*(xx.ravel()[:,None]*kx + yy.ravel()[:,None]*ky))
    up,vals,vh=svd(U,full_matrices=False,check_finite=False)
    take=vals>1.e-4
    coeff=vh[take].conj().T@((up[:,take].conj().T@s.conj().T)/vals[take,None])
    cf=ball(coeff)
    _, op_sq=eigen_interval(star(cf)*cf)
    first=SIDE-1
    root=(c/q).sqrt()
    err_sq=arb(0)
    for y in range(SIDE):
        row=acb_mat([[fx[x+first][h]*fy[y+first][h]*root for h in range(q)] for x in range(SIDE)])
        delta=star(ball(s[:,y*SIDE:(y+1)*SIDE]))-row*cf
        err_sq+=frobenius_squared(delta)
    # The directed endpoints themselves are fixed real upper bounds.
    gamma=arb(1)/128
    if not (arb(9)/8)*op_sq < gamma:
        raise ArithmeticError('The proposed rational radiation coefficient is too small.')
    eta=(9*err_sq + gamma*rho).upper()
    profiles=json.loads((STUDY/'submission/single/coherent_general_profiles.json').read_text())
    caps={}
    for idx in (0,2):
        z=profiles[idx];mu=arb(z['mu']); gap=arb(1)/1024
        floor=(gap-eta/mu).lower()
        if not floor>0:raise ArithmeticError('The radiation inequality does not absorb receiver power.')
        # The exact coefficients are taken from the archived, normalized K.
        w=expand_block(z['sector']);t=expand_block(z['optical']);two_l=[]
        for real,opt,row in zip(w,t,b,strict=True):
            sgn=acb(arb(.75)*arb(float(real)),arb(float(opt))+gamma/mu)
            two_l.extend(sgn*acb(float(v.real),float(v.imag)) for v in row)
        rhs=acb_mat(N,b.shape[1],two_l);_,a=eigen_interval(star(rhs)*rhs)
        cap=(a/floor**2).upper()
        caps[str(idx)]={'whole_box_energy_upper':hi(cap),'normalized_state_floor_lower':lo(floor),'linear_operator_norm_squared_upper':hi(a),'data_access':'no Y, no r, no truth'}
    out={'scope':'All currents for every unit coherent wave and any real pixelwise contrast in [0,.75], independently of its readings, in the exact stored native discrete model.','native_cells':N,'angular_columns':q,'numeric_fit_columns':int(take.sum()),'gram_kernel_max_entry_error_upper':hi(worst),'radiation_gram_spectral_error_upper':hi(rho),'coefficient_norm_squared_upper':hi(op_sq),'receiver_residual_fro_squared_upper':hi(err_sq),'radiation_coefficient':1/128,'isotropic_remainder_upper':hi(eta),'inherited_scaled_floor':1/1024,'profiles':caps,'wall_seconds':perf_counter()-start,'cpu_seconds':process_time()-cpu,'precision':ctx.prec}
    path=HERE.parent/'receipts/radiation_transfer.json';path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':
    with threadpool_limits(1):main()
