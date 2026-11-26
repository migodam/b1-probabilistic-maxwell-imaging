"""One-sided 96² replay of a saved physical candidate without an Arb dense LU.

Arb encloses the exact lag entries of the released binary64 FFT symbol.
Dense binary64 calculations are proposals with conservative normwise errors.
Streamed inverse columns keep the largest two matrices in ordinary binary64.
The unmodified native active block is used throughout the final check.
"""
from __future__ import annotations

import argparse

from shared import Model, OUT, elapsed, now, write_json
from certify import gamma, norm, up

from flint import acb, ctx
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits


def kernel_enclosure(symbol):
    # Separable inverse DFT of exact dyadic samples. Row transforms run over
    # the x coordinate and column transforms over y.
    cols = symbol.shape[1]
    rows = [acb.dft([acb(complex(v)) for v in row], inverse=True) for row in symbol]
    lag_columns = [acb.dft([row[x] for row in rows], inverse=True) for x in range(cols)]
    lag = np.array([[complex(v) for v in col] for col in lag_columns]).T
    max_error = max((abs(v-acb(complex(lag[i,j]))).upper()
                     for j, col in enumerate(lag_columns) for i,v in enumerate(col)))
    return lag, up(np.nextafter(float(max_error),np.inf))


def active_matrix(model, x, active):
    started = now()
    lo, lag_err = kernel_enclosure(model.green.spectrum)
    x_act = x[active]
    n = len(active)
    gy, gx = np.divmod(active, model.side)
    G_l = model.green.left[active]
    G_r = model.green.right[:,active]
    A = np.empty((n,n), dtype=complex, order="F")
    green_norm_sq = 0.
    for start in range(0,n,80):
        stop = min(start+80,n)
        sel = slice(start,stop)
        dy = (gy[sel,None]-gy[None,:]) % lo.shape[0]
        dx = (gx[sel,None]-gx[None,:]) % lo.shape[1]
        approx = lo[dy,dx] + G_l[sel] @ G_r
        green_norm_sq = up(green_norm_sq + up(norm(approx)**2))
        A[sel,:] = -x_act[sel,None]*approx
        jj = np.arange(start,stop)
        A[jj,jj] += 1.
    g_hat_norm = up(green_norm_sq**.5)
    lo_fro = up(n * max(abs(lo).ravel()))
    prod_fro = up(norm(G_l)*norm(G_r))
    lr_error = up(gamma(16*G_l.shape[1]+32)*prod_fro)
    matrix_error_g = up(n*lag_err + lr_error + gamma(8)*(lo_fro+prod_fro+lr_error))
    x_max = float(np.max(x_act))
    matrix_error_a = up(x_max * matrix_error_g +
                        gamma(12)*((n**.5) + x_max*g_hat_norm))
    return A, matrix_error_a, dict(lag_entry_error=lag_err, green_frobenius_error=matrix_error_g,
                                   state_matrix_frobenius_error=matrix_error_a, **elapsed(started))


def run(case, roi_num, source="zero", stationary=False):
    t = now()
    model = Model.read(case)
    p = OUT/case
    if stationary:
        p /= "stationary"
    with np.load(p/"baseline.npz", allow_pickle=False) as file:
        roi = file["rois"][roi_num-1].ravel()
    name = f"roi_{roi_num:02d}"
    if source == "zero":
        with np.load(p/"b1_predictions.npz", allow_pickle=False) as file:
            x = file["x_"+name].ravel()
            current = file["j_"+name]
        if np.any(x[roi] != 0.):
            raise ValueError("Illegal zero-region candidate")
    else:
        with np.load(p/"b1_controls.npz", allow_pickle=False) as file:
            x = file["physical_chi"].ravel()
            current = file["physical_j"]
    if np.any(x < 0.) or np.any(x > 4.5) or np.any(~np.isfinite(x)):
        raise ValueError("Candidate violates the shared material box")
    active = np.flatnonzero(x)
    n = len(active)
    a, a_err, setup = active_matrix(model,x,active)
    a_n = norm(a)
    print(case, name, "assembled", n, setup, flush=True)
    stage = now()
    factor = lu_factor(a.copy(order="F"),overwrite_a=True,check_finite=False)
    factor_time = elapsed(stage)
    print(case, name, "factor", factor_time, flush=True)
    inv_sq, d_sq = 0., 0.
    stage = now()
    step = 112
    for count, start in enumerate(range(0,n,step)):
        stop = min(start+step,n)
        width = stop-start
        row = np.arange(start,stop)
        col = np.arange(width)
        rhs = np.zeros((n,width),complex,order='F')
        rhs[row,col] = 1.
        rj = lu_solve(factor,rhs,overwrite_b=True,check_finite=False)
        inv_sq = up(inv_sq + up(norm(rj)**2))
        d = a @ rj
        d[row,col] -= 1.
        d_sq = up(d_sq + up(norm(d)**2))
        if count % 12 == 0:
            print(case, name, "inverse columns", stop, n, elapsed(stage), flush=True)
    stream_time = elapsed(stage)
    nr = up(inv_sq**.5)
    nd = up(d_sq**.5)
    dot = gamma(16*n+32)
    defect = up(nd + a_err*nr + dot*(a_n*nr+n**.5))
    res = dict(case=case, region=name, material_active=n, roi_cells=int(roi.sum()),
               source=source, baseline_family="stationary" if stationary else "dual",
               roi_peak=float(x[roi].max()), roi_mean=float(x[roi].mean()),
               precision_arb_bits=ctx.prec, setup=setup, factor_time=factor_time,
               inverse_stream_time=stream_time, inverse_defect_upper=defect,
               contract="Exact released binary64 nominal B,S,Y,r and FFT symbol plus exact low-rank operator")
    if defect >= 1.:
        res.update(verified=False,reason="inverse defect exceeds one")
    else:
        inv_bound = up(nr/np.nextafter(1-defect,0.))
        ja = current[active]
        ba = model.B[active]
        sa = model.S[:,active]
        f = x[active,None] * ba
        f_err = up(gamma(8) * float(np.max(x))*norm(ba))
        jn, fn, sn = norm(ja), norm(f), norm(sa)
        state_upper = up(norm(a @ ja - f) + a_err*jn + f_err + dot*(a_n*jn+fn))
        data_obs = float(np.linalg.norm(sa @ ja - model.Y))
        data_upper = up(up(data_obs*(1+gamma(8*model.Y.size+32))) +
                        dot*(sn*jn+norm(model.Y)) + sn*inv_bound*state_upper)
        res.update(verified=data_upper<model.radius, radius=model.radius, residual_float=data_obs,
                   residual_upper=data_upper, residual_upper_over_radius=up(data_upper/model.radius),
                   state_defect_upper=state_upper, exact_inverse_norm_upper=inv_bound)
    res.update(**elapsed(t))
    suffix = "" if source == "zero" else "_matched"
    write_json(p/f"native_one_sided_{name}{suffix}.json",res)
    print(case,name,res,flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("case")
    cli.add_argument("roi",type=int)
    cli.add_argument("--source",choices=("zero","matched"),default="zero")
    cli.add_argument("--stationary",action="store_true")
    a = cli.parse_args()
    ctx.prec = 112
    with threadpool_limits(limits=1,user_api="blas"):
        run(a.case,a.roi,a.source,a.stationary)
