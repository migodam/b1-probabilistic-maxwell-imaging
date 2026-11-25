"""One-sided necessary dark-space constraints for the whole two-group PSD cone."""
import argparse
from fractions import Fraction
import json
import sys

import numpy as np
from flint import acb,acb_mat,arb,ctx
from scipy.linalg import eigh,svd
from scipy.optimize import minimize_scalar

from search import Geometry, HERE, box, clock, elapsed, geom_path, read_spec, save_json, threadpool_limits
from operator_bounds import ball,inverse_dft_table
sys.path.insert(0,str(HERE.parent))
from interval_backend import adj,ball as matrix_ball,row_scaled


def trial_region(g,cells,padding):
    lo=g.xy[cells].min(axis=0)-padding
    hi=g.xy[cells].max(axis=0)+padding
    return np.flatnonzero(np.all((g.xy>=lo)&(g.xy<=hi),axis=1))


def exact_kernel(g):
    table=[[ball(v) for v in row] for row in g.symbol]
    return [tuple(row) for row in inverse_dft_table(table)]


def exact_restricted(g,kernel,restriction):
    xy=g.xy[restriction]
    pad=g.symbol.shape[0]
    delta=(xy[:,None,:]-xy[None,:,:])%pad
    return acb_mat(len(xy),len(xy),[kernel[int(y)][int(x)] for x,y in delta.reshape(-1,2)])


def herm(a):
    return (a+a.conj().T)/2


def upper(a):
    return float(np.nextafter(float(a.upper()),np.inf))


def exact_energies(v,g,s,tempered_out,tempered_in,l,u):
    # S*S projection is computed symbolically through a verified Gram solve.
    raw=matrix_ball(v[:,None])
    h=adj(s)
    w=raw-h*(s*h).solve(s*raw)
    field=g*w
    term=lambda z,k: z[k,0].real**2+z[k,0].imag**2
    results=[]
    for group in (tempered_out,tempered_in):
        val=arb(0)
        for k in group:
            wk=w[k,0];ek=field[k,0]
            val += term(w,k) - (l[k]+u[k])*(wk.conjugate()*ek).real + l[k]*u[k]*term(field,k)
        results.append(upper(val))
    residual_upper = upper(sum((abs(x).upper()**2 for x in (s*w).entries()),arb(0)).upper().sqrt())
    proven_nonzero = any(abs(x).lower()>0 for x in w.entries())
    return results, residual_upper,proven_nonzero


def proof_bounds(entries):
    low,high=Fraction(0),Fraction(1)
    direct=None
    for i,entry in enumerate(entries):
        if not entry['vector_certified_nonzero']:
            continue
        a,b=map(Fraction,entry['strict_energy_upper_outside_inside'])
        if a<0 and b<0:
            direct=i
            break
        if a<0 and b>=0:
            low=max(low,-a/(b-a))
        if b<0 and a>=0:
            high=min(high,a/(a-b))
    return {'lower_tau_exact_ratio':[low.numerator,low.denominator],
            'upper_tau_exact_ratio':[high.numerator,high.denominator],
            'necessary_tau_interval':[float(low),float(high)],'direct_both_negative_trial':direct,
            'entire_nonzero_two_group_psd_cone_obstructed':direct is not None or low>high}


def one(g,kernel,node,full=False):
    started=clock()
    cells,profile,l_arg,u_arg=box(node)
    p=2 if g.side==32 else 4
    j=np.arange(g.n) if full else trial_region(g,cells,p)
    inside=np.isin(j,cells)
    lu,uu=l_arg*profile,u_arg*profile
    l,u=np.zeros(len(j)),np.full(len(j),1/32)
    # Both index arrays are sorted in raster order.
    l[inside],u[inside]=lu,uu
    g_num=g.restrict(j)
    s=g.S[:,j]
    svals=svd(s,full_matrices=True,check_finite=False)[2]
    z=svals[s.shape[0]:].conj().T
    e=g_num@z
    jm,jo=np.flatnonzero(inside),np.flatnonzero(~inside)
    def group(rows):
        zz,ee=z[rows],e[rows]
        mid=zz.conj().T@(((l[rows]+u[rows])/2)[:,None]*ee)
        return herm(zz.conj().T@zz-mid-mid.conj().T+ee.conj().T@((l[rows]*u[rows])[:,None]*ee))
    out,in_ = group(jo),group(jm)
    evaluations=[]
    def weakest(tau, vector=False):
        h=(1-tau)*out+tau*in_
        if vector:
            val,v=eigh(h,subset_by_index=[0,0],check_finite=False)
            return float(val[0]),z@v[:,0]
        val=eigh(h,subset_by_index=[0,0],eigvals_only=True,check_finite=False)
        evaluations.append((float(tau),float(val[0])))
        return -float(val[0])
    ans=minimize_scalar(weakest,bounds=(0.,1.),method='bounded',options={'xatol':1.e-7,'maxiter':40})
    tau=float(ans.x)
    weak=-float(ans.fun)
    ee=exact_restricted(g,kernel,j)
    se=matrix_ball(s)
    l_exact=[arb(float(v)) for v in l]
    u_exact=[arb(float(v)) for v in u]
    attempts=[]
    used=[]
    # At an interior concave maximum, two neighboring weak directions can
    # bracket an empty feasible interval; endpoints cover near-zero groups.
    coords=(tau,max(0.,tau-1.e-5),min(1.,tau+1.e-5),max(0.,tau-2.e-4),min(1.,tau+2.e-4),0.,1.,
            *[float(2.**w/(1+2.**w)) for w in (-20,-16,-12,-8,-4,0)])
    for t in dict.fromkeys(coords):
        t_clock=clock()
        value,v=weakest(t,True)
        try:
            energies,err,nonzero=exact_energies(v,ee,se,jo,jm,l_exact,u_exact)
        except (ValueError,ZeroDivisionError) as exc:
            attempts.append(dict(tau=t,numeric_dark_value=value,unresolved=str(exc),cost=elapsed(t_clock)))
            continue
        used.append(v)
        attempts.append(dict(tau=t,numeric_dark_value=value,strict_energy_upper_outside_inside=energies,
                             receiver_null_by_exact_verified_projection=True,
                             receiver_zero_enclosure_radius_diagnostic=err,vector_certified_nonzero=nonzero,
                             cost=elapsed(t_clock)))
        if proof_bounds([r for r in attempts if 'vector_certified_nonzero' in r])['entire_nonzero_two_group_psd_cone_obstructed']:
            break
    valid=[r for r in attempts if 'vector_certified_nonzero' in r]
    result=proof_bounds(valid)
    label='full_' if full else ''
    file=HERE/f'dark_frame_{label}{g.side}_{node["index"]}_{node["family"]}.npz'
    np.savez_compressed(file,ids=j,trial_vectors=np.asarray(used),support_cells=cells)
    return {**node,**result,'candidate_cells':len(cells),'trial_current_cells':len(j),
            'numerical_largest_minimum_dark_value':weak,'numerical_tau':tau,'optimizer_success':bool(ans.success),
            'optimizer_trials':evaluations,'one_sided_trials':attempts,'vectors':file.name,'cost':elapsed(started)}


def run(full32=False):
    spec=read_spec()
    results=[]
    scales=(32,) if full32 else (32,96)
    for side in scales:
        g=Geometry(geom_path(side)); kernel=exact_kernel(g)
        for node in (x for x in spec['nodes'] if x['side']==side):
            rec=one(g,kernel,node,full=full32); results.append(rec)
            dest='dark_cone_full_32.json' if full32 else 'dark_cone.json'
            save_json(HERE/dest,results)
            print(side,node['index'],rec['numerical_largest_minimum_dark_value'],rec['necessary_tau_interval'],
                  rec['entire_nonzero_two_group_psd_cone_obstructed'],rec['cost'],flush=True)


def replay_full_32():
    """Recompute the signed energies from the stored trials and exact original operator."""
    started=clock()
    rows=json.loads((HERE/'dark_cone_full_32.json').read_text())
    g=Geometry(geom_path(32))
    ids=np.arange(g.n)
    G=exact_restricted(g,exact_kernel(g),ids)
    S=matrix_ball(g.S)
    results=[]
    for row in rows:
        k=row['direct_both_negative_trial']
        if k is None:
            continue
        cells,profile,l,u=box(row)
        lower=np.zeros(g.n); upper=np.full(g.n,1/32)
        lower[cells],upper[cells]=l*profile,u*profile
        outside=np.setdiff1d(ids,cells)
        with np.load(HERE/row['vectors'],allow_pickle=False) as record:
            if not (np.array_equal(ids,record['ids']) and np.array_equal(cells,record['support_cells'])):
                raise ValueError('stored trial geometry mismatch')
            trial=record['trial_vectors'][k]
        bounds,residual,nonzero=exact_energies(trial,G,S,outside,cells,
                                              [arb(float(x)) for x in lower],[arb(float(x)) for x in upper])
        if not (nonzero and max(bounds)<0):
            raise ArithmeticError('whole-cone obstruction replay failed')
        results.append(dict(index=row['index'],observation=row['oid'],energy_uppers=bounds,
                            exact_receiver_null=True,receiver_residual_enclosure=residual))
    save_json(HERE/'dark_replay_full_32.json',{'rows':results,'precision_bits':ctx.prec,'cost':elapsed(started)})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    switches=parser.add_mutually_exclusive_group()
    switches.add_argument('--full32',action='store_true')
    switches.add_argument('--replay-full32',action='store_true')
    args=parser.parse_args()
    ctx.prec,ctx.threads=192,1
    with threadpool_limits(limits=1):
        if args.replay_full32:
            replay_full_32()
        else:
            run(args.full32)
