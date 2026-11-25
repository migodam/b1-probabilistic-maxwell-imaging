"""Public-data witness proposals and checked physical replay on the full grid."""
import argparse
from collections import Counter
import json
import os

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.linalg import solve
from scipy.optimize import lsq_linear
from scipy.sparse.linalg import LinearOperator, gmres
from threadpoolctl import threadpool_limits

from certified_map import ARITH, add_bound, fnorm, op_upper
from common import (FREEZE, HERE, Geometry, clock, elapsed, geom_path, pattern, place,
                    save_json)
from sector_map import SectorGeometry, matrix_epsilon


def forecast(g, x, cells, profile, Y, r, background_cap=0., extra_slack=False):
    measured, current, info = g.forward(x)
    if info:
        return None
    residual = fnorm(Y-measured)
    if residual > 0.985*r and residual < 2.75*r:
        local = x[cells]
        block = g.restrict(cells)
        nom = np.eye(len(cells))-local[:,None]*block
        sensor = solve(nom.T, g.S[:, cells].T, check_finite=False).T
        field = (g.B + g.apply(current))[cells]
        differential = np.einsum("ri,is->rsi", sensor, field).reshape(Y.size, -1)
        data = Y-measured
        design = np.r_[differential.real, differential.imag]
        design = np.vstack((design,0.0005*np.eye(len(cells))))
        target = np.r_[data.real.ravel(), data.imag.ravel(), np.zeros(len(cells))]
        lo,hi = FREEZE["material_classes"]["primary_core_range"]
        lower, upper = lo*profile-local, hi*profile-local
        try:
            correction = lsq_linear(design,target,bounds=(lower,upper),
                                    method="bvls",max_iter=32).x
        except ValueError:
            return x,current,residual
        for step in (1., .5):
            trial = x.copy()
            trial[cells] = np.minimum(hi*profile,np.maximum(lo*profile,local+step*correction))
            model, j, info = g.forward(trial)
            err = fnorm(Y-model)
            if info == 0 and err < residual:
                x, current, residual = trial, j, err
            if residual < 0.985*r:
                break
    if background_cap and (0.95 if extra_slack else 1.)*r <= residual < 1.8*r:
        right = Y-g.S@current
        field = g.B + g.apply(current)
        # The shift-invariant Green operator is symmetric to the working
        # precision of the proposal. The final check uses exact coefficients.
        adj=LinearOperator((g.n,g.n),matvec=lambda v:v-g.apply(x*v),dtype=complex)
        jac_sensor=np.empty_like(g.S)
        good=True
        for i,row in enumerate(g.S):
            jac_sensor[i],code=gmres(adj,row,rtol=2e-8,restart=60,maxiter=25)
            good &= code==0
        if good:
            jac=np.einsum("ri,is->rsi",jac_sensor,field).reshape(Y.size,g.n)
            A=np.r_[jac.real,jac.imag]
            lo,hi=FREEZE["material_classes"]["primary_core_range"]
            lower=np.zeros(g.n)
            upper=np.full(g.n,background_cap)
            lower[cells],upper[cells]=lo*profile,hi*profile
            width=upper-lower
            A*=width
            delta=np.zeros(g.n)
            lb,ub=(lower-x)/width,(upper-x)/width
            target=np.r_[right.real.ravel(),right.imag.ravel()]
            effect=target.copy()
            for _ in range(32):
                gradient=A.T@effect
                gradient[(delta <= lb)&(gradient < 0)] = 0
                gradient[(delta >= ub)&(gradient > 0)] = 0
                Av=A@gradient
                denom=Av@Av
                if denom <= 1e-30:
                    break
                step=max(0.,effect@Av/denom)
                candidate=np.minimum(ub,np.maximum(lb,delta+step*gradient))
                delta=candidate
                effect=target-A@delta
            for q in (1.,.5):
                trial=np.minimum(upper,np.maximum(lower,x+q*width*delta))
                model,j,info=g.forward(trial)
                err=fnorm(Y-model)
                if info==0 and err<residual:
                    x,current,residual=trial,j,err
                if residual<0.985*r:
                    break
    return x,current,residual


def certify_witness(sector, x, cells, current, Y, r, b_uniform=0.):
    g = sector.g
    start = clock()
    xref=x.copy()
    xref[cells]=b_uniform
    reference_max=float(np.max(xref))
    margin=1-reference_max*sector.full_g
    if margin <= 0:
        return None
    active_delta = x[cells]-b_uniform
    n = g.n
    # Validate the resolvent of the known uniform reference. Their
    # complements remain exact incidence grids, not fitted exterior pixels.
    if reference_max:
        xy = g.xy
        lag = (xy[:,None,:]-xy[cells][None,:,:])%sector.pad
        Z0 = g.lag[lag[:,:,1],lag[:,:,0]]
        operator = json.loads((HERE/"results"/"operator_bounds.json").read_text())[str(g.side)]
        q_error = np.sqrt(n*len(cells))*operator["kernel_entry_error"]
        Z = np.zeros_like(Z0)
        z_error_sq = 0.
        z_norm_sq = 0.
        for block in range(0,len(cells),32):
            part = slice(block,block+32)
            source = Z0[:,part]
            proposal = source.copy()
            for _ in range(13):
                proposal=source+g.apply(xref[:,None]*proposal)
            ref_current=xref[:,None]*proposal
            qhat=sector.apply(ref_current,g.symbol)
            epsfft=(sector.rel_fft_bound * fnorm(ref_current)*sector.full_g +
                    matrix_epsilon(sector.pad)*reference_max*fnorm(proposal)*sector.full_g)
            residual=fnorm(source-proposal+qhat)
            out=(residual+epsfft+matrix_epsilon(sector.pad)*
                 (fnorm(source)+fnorm(proposal)+fnorm(qhat)))
            z_error_sq += out*out
            z_norm_sq += fnorm(proposal)**2
            Z[:,part] = proposal
        z_err = (np.sqrt(z_error_sq) + q_error) / margin
        zv = Z[cells]
        znorm = np.sqrt(z_norm_sq) + z_err
    else:
        zv = g.restrict(cells)
        operator = json.loads((HERE/"results"/"operator_bounds.json").read_text())[str(g.side)]
        z_err = len(cells)*operator["kernel_entry_error"]
        znorm = sector.full_g
    a = np.eye(len(cells)) - zv*active_delta
    ae = z_err*np.max(np.abs(active_delta)) + (
        matrix_epsilon(len(cells)) * (1+fnorm(zv)*np.max(np.abs(active_delta)))
    )
    inv = solve(a,np.eye(len(cells)),check_finite=False)
    invnorm = op_upper(inv)
    rho = (fnorm(np.eye(len(cells))-inv@a) + invnorm*ae +
           matrix_epsilon(len(cells))*(np.sqrt(len(cells))+fnorm(inv)*fnorm(a)))
    if rho >= 1:
        return None
    # Exact Woodbury factor of I-GX. Convert the raw current residual
    # to a field residual via e_hat=B+G j_hat, then back to current.
    core = invnorm/(1-rho)
    norm_inverse = (1+znorm*np.max(np.abs(active_delta))*core)/margin
    wave = sector.apply(current,g.symbol)
    e = g.B+wave
    error_gj = sector.rel_fft_bound*sector.full_g*fnorm(current)
    linear_residual = current - x[:,None]*e
    product_arithmetic = matrix_epsilon(n)*(1+fnorm(current)+np.max(x)*(fnorm(g.B)+fnorm(wave)))
    current_equation_upper = fnorm(linear_residual)+np.max(x)*error_gj+product_arithmetic
    current_bound = (1+np.max(x)*norm_inverse*sector.full_g)*current_equation_upper
    observed = g.S@current
    predicted_residual = (fnorm(Y-observed)+norm_sensor(sector)*current_bound+
                          matrix_epsilon(n)*(1+fnorm(Y)+fnorm(g.S)*fnorm(current)))
    if not predicted_residual < r:
        return None
    return dict(checked_residual_upper=float(predicted_residual),noise_radius=r,
                residual_ratio_upper=float(predicted_residual/r),
                core_inverse_defect_upper=float(rho),full_field_inverse_upper=float(norm_inverse),
                physical_state_residual_upper=float(current_equation_upper),
                reference_value_on_active=b_uniform,
                real_exterior_min=float(np.min(np.delete(x,cells))),
                real_exterior_max=float(np.max(np.delete(x,cells))),
                **elapsed(start))


def norm_sensor(sector):
    # Frobenius bounds the full spectral norm and suffices in this audit.
    return add_bound(sector.normS)


def run(side, shape):
    started=clock()
    g = Geometry(geom_path(side))
    sector = SectorGeometry(g,FREEZE["model"]["background_upper"])
    amplitudes=(4.0,3.875,4.125,3.9375,4.0625,3.8125,4.1875)
    offsets,p = pattern(side,shape)
    acq=json.loads((HERE/"public"/"acquisition.json").read_text())
    roster={}
    for oid,info in acq["observations"].items():
        if info["side"]!=side or info["shape"]!=shape:
            continue
        data=np.load(HERE/"public"/f"{oid}.npz",allow_pickle=False)
        Y,r=data["Y"],float(data["r"])
        roster[oid]={}
        for b in ((.03125,) if "__weak__" in oid else (0.,.03125)):
            one=clock()
            manifest=HERE/"results"/f"{oid}__b{b:g}__primary.npz"
            if not manifest.exists():
                print("waiting for",manifest.name,flush=True)
                continue
            possibilities=np.load(manifest,allow_pickle=False)
            reference=possibilities["reference_error"]
            local_step=0
            proposals=0
            best_margin=np.inf
            audited=[]
            candidate_centres=np.argsort(reference)[:FREEZE["material_classes"][
                "witness_proposals_per_map"]]
            accepted=None
            for node in candidate_centres:
                xy=g.xy[node]
                cells,profile=place(side,xy,offsets,p)
                for bg in ((0.,) if not b else (b/2,0.)):
                    for scale in amplitudes:
                        proposals+=1
                        x=np.full(g.n,bg)
                        x[cells]=scale*profile
                        tentative=forecast(g,x,cells,profile,Y,r,b if bg else 0.)
                        if tentative is None:
                            continue
                        xx,j,pres=tentative
                        best_margin=min(best_margin,pres)
                        if pres>=r:
                            continue
                        receipt=certify_witness(sector,xx,cells,j,Y,r,bg)
                        if receipt is not None:
                            witness=HERE/"results"/f"{oid}__b{b:g}__witness.npz"
                            np.savez_compressed(witness,material=xx,current=j,centre=xy,
                                                geometric_profile=profile,active_cells=cells)
                            receipt.update(centre=[int(c) for c in xy],candidate_scale=scale,
                                           proposal_residual=float(pres),file=witness.name)
                            accepted=receipt
                            break
                    if accepted is not None:
                        break
                if accepted is not None:
                    break
            roster[oid][str(b)]=dict(witness=accepted,proposals=proposals,
                   minimum_proposal_residual=None if not np.isfinite(best_margin) else float(best_margin),
                   **elapsed(one))
            print("witness",side,shape,oid,b,
                  None if accepted is None else accepted["residual_ratio_upper"],
                  "proposals",proposals,elapsed(one),flush=True)
    save_json(HERE/"results"/f"witnesses_{side}_{shape}.json",dict(witnesses=roster,elapsed=elapsed(started)))


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("side",type=int)
    a.add_argument("shape")
    args=a.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        run(args.side,args.shape)
