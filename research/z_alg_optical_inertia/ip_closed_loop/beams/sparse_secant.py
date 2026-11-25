"""Try closing the native two/three-port gap using exact real-box secant support."""
import argparse
import json

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, elapsed, geometry, herm, normalize, save
from interval_certify import ProofGrid, ball, hi, lo, midpoint, norm_upper, tr
from native_paired import green_kernel_roundoff, pad_all
from native_residual import ResidualBoundBox


def bound(box,a,q):
    e=pad_all(box.E*a,(box.e_error*norm_upper(a)).upper())
    d=pad_all(tr(q)*box.D,(box.d_error*norm_upper(q)).upper())
    nominal=(tr(q)*box.Y*a)[0,0].real
    dw=[d[0,j]*acb(v) for j,v in enumerate(box.width)]
    first=sum((abs((t*e[k,0]).real).upper() for k,t in enumerate(dw)),arb(0))
    dn=sum((abs(u).upper()**2 for u in dw),arb(0)).sqrt().upper()
    second=dn*box.eta*norm_upper(e)/(1-box.eta)
    offset=norm_upper(q)*norm_upper(a)*box.y_error
    rad=first+second+offset
    return (nominal-rad).lower(),(nominal+rad).upper(),hi(first),hi(second)


def receive_direction(difference, references, a):
    """Find a receive direction reducing projected real-box nuisance."""
    initial=np.r_[difference.real,difference.imag]/np.linalg.norm(difference)
    center=np.r_[difference.real,difference.imag]
    penalty=[]
    curvature=[]
    for _,d,e,width,eta in references:
        columns=d*width[None,:]
        delta_e=e@a
        p=columns*delta_e[None,:]
        penalty.append(np.vstack((p.real,p.imag)))
        r=np.vstack((columns.real,columns.imag))
        i=np.vstack((columns.imag,-columns.real))
        curvature.append((r@r.T+i@i.T,np.linalg.norm(delta_e)*eta/(1-eta)))

    def objective(x):
        cost=-float(x@center)
        gradient=-center.copy()
        for p,(matrix,scale) in zip(penalty,curvature,strict=True):
            z=x@p
            w=np.sqrt(z**2+2.e-21)
            cost+=float(np.sum(w))
            gradient+=p@(z/w)
            v=matrix@x
            root=np.sqrt(float(x@v)+2.e-21)
            cost+=root*scale
            gradient+=v*(scale/root)
        return cost,gradient

    constraint=dict(type='ineq',fun=lambda x:1-float(x@x),jac=lambda x:-2*x)
    fitted=minimize(objective,initial,jac=True,method='SLSQP',constraints=[constraint],
                    options={'maxiter':80,'ftol':2.e-13})
    direct=fitted.x/np.linalg.norm(fitted.x)
    x=min((initial,direct),key=lambda x:objective(x)[0])
    p=len(difference)
    return x[:p]+1j*x[p:]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=(96,))
    parser.add_argument('--optimize-receive',action='store_true')
    args=parser.parse_args()
    src=HERE.parent/'pair'/f'power_{args.side}'
    p=json.loads((src/'proposal.json').read_text())
    count=json.loads((HERE/f'sparse_physical_power_{args.side}.json').read_text())['rows']
    bank=np.load(src/'public_classes.npz')
    gram=bank['power']
    start=clocks()
    ctx.prec=128
    ctx.threads=3
    geo=geometry(args.side)
    grid=ProofGrid(geo)
    eps=green_kernel_roundoff(grid)
    b=ball(geo.B)
    exact_power=tr(b)*b
    exact_power=acb(1/exact_power[0,0].real)*exact_power
    results={}
    stored={}
    for i,meta in p['rows'].items():
        tick=clocks()
        key=i+'__'
        idx=bank[key+'cells']
        numeric,exact=grid.restricted(idx)
        forms=[ResidualBoundBox(exact,numeric,c,bank[key+'widths'],eps*len(idx)) for c in bank[key+'centers']]
        D=midpoint(forms[1].Y)-midpoint(forms[0].Y)
        allowed=[r['indices'] for r in count[i]['caps']['2']['rows'] if not r['common_reading_implied_for_any_coefficients']]
        rows=[]
        act={}
        for inds in allowed:
            M=gram[np.ix_(inds,inds)]
            H=herm(D[:,inds].conj().T@D[:,inds])
            d,v=eigh(H,M,check_finite=False)
            weight=np.linalg.cholesky(M).conj().T
            z=weight@v[:,-1]
            theta=np.arctan2(abs(z[1]),abs(z[0]))
            phi=np.angle(z[1])-np.angle(z[0])
            def action(x):
                y=np.zeros(len(gram),complex)
                y[inds]=np.linalg.solve(weight,np.r_[np.cos(x[0]),np.exp(1j*x[1])*np.sin(x[0])])
                return normalize(y,gram,meta['fluence'])
            # Cache the original object response matrices once: optimization
            # uses only those old-data public class models.
            original=[(midpoint(x.Y),midpoint(x.D),midpoint(x.E),x.width_np,float(x.eta)) for x in forms]
            def approximate(a,q=None):
                ydelta=D@a
                norm=float(np.linalg.norm(ydelta))
                q=ydelta/norm if q is None else q
                r=[]
                for _,d,e,width,eta in original:
                    dw=(q.conj()@d)*width
                    ex=e@a
                    r.append(np.sum(abs((dw*ex).real))+np.linalg.norm(dw)*np.linalg.norm(ex)*eta/(1-eta))
                return np.vdot(q,ydelta).real-sum(r),q
            def fast_score(x):
                return -approximate(action(x))[0]
            fitted=minimize(fast_score,np.r_[theta,phi],method='Powell',bounds=((0.,np.pi/2),(phi-np.pi,phi+np.pi)),
                            options={'maxiter':35,'xtol':5e-5,'ftol':1e-7})
            x=min((fitted.x,np.r_[theta,phi]),key=fast_score)
            raw=action(x)
            guess,q=approximate(raw)
            if args.optimize_receive:
                q=receive_direction(D@raw,original,raw)
                raw=action(minimize(lambda x:-approximate(action(x),q)[0],x,method='Powell',
                                    bounds=((0.,np.pi/2),(phi-np.pi,phi+np.pi)),
                                    options={'maxiter':35,'xtol':5e-5,'ftol':1e-7}).x)
                q=receive_direction(D@raw,original,raw)
                guess,q=approximate(raw,q)
            exact_a=ball(raw[:,None])
            exact_a=acb((arb(meta['fluence'])/(tr(exact_a)*exact_power*exact_a)[0,0].real).sqrt())*exact_a
            q=ball(q[:,None])
            ends=[bound(m,exact_a,q) for m in forms]
            noise=arb(meta['new_radius'])*norm_upper(q)
            gain=ends[1][0]-ends[0][1]-2*noise
            rows.append(dict(indices=inds,strict=bool(gain>0),whole_box_gap_over_noise_lower=lo(gain/arb(meta['new_radius'])),
                             projected_noise=hi(noise),estimated_clean_distance_over_noise=guess/meta['new_radius'],
                             support_first_order=[x[2] for x in ends],support_higher_order=[x[3] for x in ends]))
            act[tuple(inds)]=(raw,midpoint(q))
        best=max(rows,key=lambda r:r['whole_box_gap_over_noise_lower'])
        stored[key+'best_two_actions'],stored[key+'best_two_receive_direction']=act[tuple(best['indices'])]
        print(i,best,flush=True)
        results[i]=dict(best=best,candidates=rows,**elapsed(tick))
    suffix='_receiver' if args.optimize_receive else ''
    name=HERE/f'sparse_secant_power_{args.side}{suffix}'
    atomic_write(name.with_suffix('.json'),dict(rows=results,physical_inputs='Old-compatible declared class material boxes only, no acquired new data read',**elapsed(start)))
    save(name.with_suffix('.npz'),**stored)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
