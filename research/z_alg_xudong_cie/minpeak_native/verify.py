"""Independently certify box candidates with original native Green coefficients."""

import argparse
import json
import sys

import numpy as np
from scipy.linalg import lu_solve
from threadpoolctl import threadpool_limits

from search import HERE, STUDY, N, full_ring_inputs, elapsed, clocks

sys.path.insert(0, str(STUDY/'synthetic'))
from nonvacuity import ExactFourierGeometry, StableActiveModel, acb_mat, arb, ball, ctx, frobenius_squared, hi, lo
sys.path.insert(0, str(STUDY/'submission/single'))
from local_moment import ExactGeometry


def norm(z):
    return frobenius_squared(z).sqrt()


def embed(a, ids):
    temp=np.zeros((N,a.shape[1]),dtype=a.dtype)
    temp[ids]=a
    return temp


def active_rows(a, ids):
    return acb_mat(len(ids),a.ncols(),[a[int(i),j] for i in ids for j in range(a.ncols())])


def verified_ports(geometry,fft,material):
    """Exact stored DFT; residual correction with a verified approximate inverse."""
    t0=clocks()
    model=StableActiveModel(geometry,material)
    inv_bound,details=model.inverse_bound()
    j=lu_solve(model.lu,model.m[:,None]*model.B,check_finite=False)
    w=lu_solve(model.lu,model.S.T,trans=1,check_finite=False).T
    x,wq,s,b=map(ball,(j,w,model.S,model.B))
    ms=list(map(arb,map(float,model.m)))
    g=active_rows(fft.green(ball(embed(j,model.active))),model.active)
    state_defect=acb_mat(model.n,b.ncols(),[ms[i]*(b[i,k]+g[i,k])-x[i,k]
                       for i in range(model.n) for k in range(b.ncols())])
    # Exact G^T applied with a conjugated input and Hermitian-adjoint FFT.
    wm=acb_mat(N,wq.nrows())
    for k,i in enumerate(model.active):
        for v in range(wq.nrows()):
            wm[int(i),v]=wq[v,k].conjugate()*ms[k]
    wt=active_rows(fft.green(wm,adjoint=True),model.active).conjugate().transpose()
    receiver_defect=s-wq+wt
    correction = (norm(receiver_defect).upper()*inv_bound*norm(state_defect).upper()).upper()
    ports = s*x + wq*state_defect
    return ports,correction,dict(active_cells=model.n, strict_inverse=details,
        primal_defect_upper=hi(norm(state_defect)),receiver_adjoint_defect_upper=hi(norm(receiver_defect)),
        corrected_port_frobenius_remainder=hi(correction),active_dense_lu=1,full_wave_source_columns=b.ncols(),
        receiver_adjoint_rows=s.nrows(),arb_green_columns=b.ncols()+s.nrows(),**elapsed(t0))


def run(args):
    start=clocks()
    s,p,b,observations=full_ring_inputs()
    ctx.prec,ctx.threads=128,1
    geo=ExactFourierGeometry(STUDY/'synthetic/public/fine_geometry.npz')
    geo.S,geo.B=p,b
    fft=ExactGeometry(p,b,s,[])
    shared=elapsed(start)
    details_by_case={}
    for obs in observations:
        if obs.noise_pct not in (args.noise or [1.,2.]) or (args.object and obs.object_id not in args.object):
            continue
        rows=[]
        for folder in args.inputs:
            path=folder/(obs.id+'.npz')
            if not path.exists():
                continue
            pay=json.loads(path.with_suffix('.json').read_text())
            if not pay['eligible_for_independent_validation'] and not args.include_outside:
                rows.append(dict(source=folder.name,status='NUMERICAL_RESIDUAL_ABOVE_RADIUS',
                                 numerical_residual=pay['selected_numeric_residual']))
                continue
            t0=clocks()
            with np.load(path,allow_pickle=False) as raw:
                material=raw['m_selected']
                if not np.array_equal(obs.y,raw['Y']):
                    raise ArithmeticError(f'Observed centers differ at {obs.id}.')
            if material.shape!=(N,) or not np.isrealobj(material) or not np.all(np.isfinite(material)):
                raise ArithmeticError('Malformed or non-real native material.')
            if material.min()<0 or material.max()>4.:
                raise ArithmeticError('Saved material violates the certified nonnegative ≤4 search envelope.')
            try:
                port,error,details=verified_ports(geo,fft,material)
                distance=norm(port-ball(obs.y))
                lower=max(0.,lo(distance.lower()-error))
                upper=hi(distance.upper()+error)
                status='CERTIFIED_FEASIBLE' if upper<=obs.padded_radius else (
                    'CERTIFIED_OUTSIDE' if lower>obs.padded_radius else 'UNRESOLVED_RADIUS')
                rows.append(dict(source=folder.name,status=status,peak=float(material.max()),
                                 distance_lower=lower,distance_upper=upper,radius=obs.padded_radius,
                                 details=details,**elapsed(t0)))
            except (ArithmeticError,ValueError) as err:
                rows.append(dict(source=folder.name,status='UNRESOLVED_NATIVE_INVERSE',error=str(err),**elapsed(t0)))
            print(obs.id,rows[-1]['source'],rows[-1]['status'],rows[-1].get('distance_upper'),flush=True)
        if rows:
            details_by_case[obs.id]=rows
    args.output.write_text(json.dumps({'interval_bits':ctx.prec,'geometry_prep':shared,
                                      'solver':'Active-state LU + inverse Neumann proof + independent Arb DFT residuals',
                                      'cases':details_by_case,**elapsed(start)},indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',nargs='+',type=type(HERE),default=[HERE/'runs',HERE/'gn'])
    parser.add_argument('--output',type=type(HERE),default=HERE/'verified.json')
    parser.add_argument('--object',action='append')
    parser.add_argument('--noise',type=float,action='append')
    parser.add_argument('--include-outside',action='store_true')
    with threadpool_limits(1):
        run(parser.parse_args())
