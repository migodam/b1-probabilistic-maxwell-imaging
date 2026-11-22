"""Apply six previously proved native profiles to the publicly frozen SRAP22 views."""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, cg, gmres
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
SRAP = STUDY.parent / 'srap_cie_array_study_20261121'
CURVE = STUDY / 'submission/curve'
SINGLE = STUDY / 'submission/single'
FINE = STUDY / 'synthetic/public/fine_geometry.npz'
ARRAY = 'PW360_10x20'
SCOPE = tuple(i for i in range(1,25) if i not in (7,10))
sys.path.insert(0,str(CURVE))
from fixed_blocks import SavedProfile, vertical
from native import Geometry, Native, Weights, ball, frobenius_squared, gershgorin
from native import high, low, mid, mult_rows, star
from search import h, save_json, times


def read(path):
    return json.loads(path.read_text())


def load_public(bundle=None):
    source = SRAP / 'data/public'
    if bundle is not None:
        with np.load(bundle,allow_pickle=False) as raw:
            if set(raw.files)!={'ids','Y','original_delta'}:
                raise ValueError('The compact public pack may contain only IDs, observations and sample radii.')
            ids,yy,delta_array=raw['ids'],raw['Y'],raw['original_delta']
        expected=[f's2d_{i:03d}__pw360_10x20' for i in SCOPE]
        if ids.tolist()!=expected or yy.shape!=(22,20,10) or delta_array.shape!=(22,):
            raise ValueError('Compact public pack does not match the registered study IDs and shape.')
    result = []
    for j,i in enumerate(SCOPE):
        object_id = f's2d_{i:03d}'
        name = object_id+'__pw360_10x20'
        if bundle is None:
            with np.load(source/f'{name}.npz',allow_pickle=False) as f:
                y = f['Y']
                delta = float(f['delta'])
                p = json.loads(f['case_json'].item())
            if (p['fine'],p['frequency'],p['array'],p['noise_relative']) != (96,400e6,ARRAY,.1):
                raise ValueError(f'Frozen public observation contract disagrees for {name}.')
        else:
            y,delta=yy[j],float(delta_array[j])
        y_norm = frobenius_squared(ball(y)).sqrt()
        radius = high(arb(delta)+y_norm/10**6)
        result.append((y,dict(id=name,object_id=object_id,stored_delta=delta,padded_radius=radius,
                              observation_norm_upper=high(y_norm))))
    return result


def load_geometry():
    with np.load(FINE,allow_pickle=False) as f:
        return f['spectrum'],f[ARRAY+'__S'],f[ARRAY+'__B']


def choices():
    for file,index,noise,tag in (
        ('global.5.json',4,2,'global_real_050'),
        ('global.75.json',1,2,'global_real_075'),
        ('global.75.json',5,2,'global_optical_075'),
        ('local.75.json',4,1,'quadrant_optical_075'),
    ):
        record=read(CURVE/file)
        source=next(r for r in record['rows'] if r['object_id']==f'optfresh_{index:02d}' and r['noise_pct']==noise)
        key=source['profile_key']
        proof=next(r for r in record['profiles'] if all(r[k]==v for k,v in key.items()))
        if (not source['strictly_excludes'] or proof['status']!='native_gap_and_inertia_verified'
                or not proof['positive_gap_lower']>0):
            raise ValueError('The original frozen native strict profile is required.')
        yield dict(label=tag,case=source['id'],block_side=2,cap=key['cap'],
                   inherited_proof=f'submission/curve/{file}',nu=proof['exact_nu'],
                   positive_gap_lower=proof['positive_gap_lower'],parameters=key,
                   inherited_profile_seconds=proof['profile_wall_seconds'],
                   inherited_profile_cpu_seconds=proof['profile_cpu_seconds'])
    fixed=read(SINGLE/'coherent_general_profiles.json')
    for obj,proof in read(SINGLE/'coherent_general_verified.json')['cases'].items():
        index=proof['profile_index']
        yield dict(label=f'fixed_4x4_{index}_075',case=obj,block_side=4,cap=.75,nu=0,
                   positive_gap_lower=proof['proved_K_scaled_gap'],
                   inherited_proof='submission/single/coherent_general_verified.json',
                   parameters=fixed[index],profile_index=index,
                   inherited_profile_seconds=proof.get('wall_seconds'),
                   inherited_profile_cpu_seconds=proof.get('cpu_seconds'))


def saved_operator(g,trial):
    p=trial['parameters']
    if 'profile_index' in trial:
        return SavedProfile(g,p,trial['cap'])
    return Native(g,Weights(p['cap'],p['scale'],tuple(p['weights']),tuple(p['optical'])))


def exact_common_response(p,nu):
    start,cpu=perf_counter(),process_time()
    g=p.g
    rec=len(g.sensor)
    mu=p.w.scale
    rhs=np.column_stack((g.sensor_star/mu,p.z[:,None]*g.incident))
    pre=LinearOperator((g.n,g.n),matvec=p.precondition,dtype=np.complex128)
    columns,counts=[],[]
    for k,r in enumerate(rhs.T):
        track=[]
        if nu==0:
            x,info=cg(p.operator,r,M=pre,rtol=3e-13,atol=0.,maxiter=700,
                      callback=lambda _:track.append(1))
        else:
            x,info=gmres(p.operator,r,M=pre,rtol=3e-13,atol=0.,restart=95,maxiter=12,
                         callback=track.append,callback_type='pr_norm')
        if info:
            raise ArithmeticError(f'Common RHS solver error {info} at index {k}.')
        columns.append(x)
        counts.append(len(track))
    x=ball(np.column_stack(columns))
    ax=p.exact_action(x)
    _,z=p.exact_scaled_coefficients()
    s=star(g.sensor_ball)/arb(mu)
    zb=mult_rows(g.incident_ball,z)
    l=acb_mat(g.n,rec+g.s,[v for i in range(g.n)
                            for v in [s[i,j] for j in range(rec)]+[zb[i,j] for j in range(g.s)]])
    r=ax-l
    m=star(x)*ax-star(x)*l-star(l)*x
    defect=star(r)*r
    return m,defect,dict(algorithm='CG' if nu==0 else 'GMRES',shared_rhs=rec+g.s,
                         iterations=counts,defect_norm_upper=high(frobenius_squared(r).sqrt()),
                         **times(start,cpu))


def bound_source(y,p,nu,gap,ngap,m,defect,meta):
    ns=y.shape[1]
    b=ball(y)
    eye=acb_mat(ns,ns)
    for i in range(ns):
        eye[i,i]=1
    c=vertical(b,eye)
    h0=star(b)*b+arb(p.w.scale)*(star(c)*m*c)
    projected_defect=star(c)*defect*c
    e=arb(p.w.scale)*sum((projected_defect[i,i].real for i in range(ns)),arb(0))
    below=(e/arb(gap)).upper()
    above=(e/arb(ngap)).upper() if nu else arb(0)
    vals,v=eigh(h(mid(h0)),check_finite=False,driver='evd')
    vals,v=vals[::-1],v[:,::-1]
    v=ball(v)
    gram=star(v)*v
    congruence=star(v)*h0*v
    lower=np.zeros(ns)
    for k in range(1,ns+1):
        hh=acb_mat(k,k,[congruence[i,j] for i in range(k) for j in range(k)])
        gg=acb_mat(k,k,[gram[i,j] for i in range(k) for j in range(k)])
        denom=gershgorin(gg,upper=True)
        floor=gershgorin(hh,upper=False)-below*denom
        if floor<=0:
            break
        lower[k-1]=max(0.,low(floor/denom))
    gmin=gershgorin(gram,upper=False)
    if gmin<=0:
        raise ArithmeticError('No certified invertible source basis.')
    for i,vv in enumerate(vals):
        congruence[i,i]-=arb(float(vv))
    bump=frobenius_squared(congruence).sqrt().upper()+above*gershgorin(gram,upper=True)
    upper=[]
    for v0 in vals:
        eig_upper=(arb(float(v0))+bump).upper()
        upper.append(high(eig_upper/gmin) if eig_upper>0 else 0.)
    jointlo=max(0.,low(sum((arb(float(q)) for q in lower[nu:]),arb(0)).sqrt()))
    jointhi=high(sum((arb(float(q)) for q in upper[nu:]),arb(0)).sqrt()) if any(upper[nu:]) else 0.
    coherentlo=max(0.,low(arb(float(lower[0])).sqrt())) if nu==0 else 0.
    coherenthi=high(arb(float(upper[0])).sqrt()) if nu==0 and upper[0]>0. else 0.
    radius=arb(meta['padded_radius'])
    return dict(**meta,exact_nu=nu,native_retained_modes_lower=int(np.count_nonzero(lower[nu:])),
                joint_distance_lower=jointlo,joint_distance_upper=jointhi,
                joint_margin_lower=low(arb(jointlo)-radius),
                joint_margin_upper=high(arb(jointhi)-radius),
                joint_status='excluded' if arb(jointlo)>radius else
                             ('verified_below_radius' if arb(jointhi)<=radius else 'unresolved'),
                coherent_distance_lower=coherentlo,coherent_distance_upper=coherenthi,
                coherent_status='inertia_exhausts_single_source' if nu else
                                ('excluded' if arb(coherentlo)>radius else
                                 ('verified_below_radius' if arb(coherenthi)<=radius else 'unresolved')),
                source_eigenvalue_positive_parts_lower=lower.tolist(),
                source_eigenvalue_positive_parts_upper=upper,
                stationary_error_lower_correction=high(below),
                stationary_error_upper_correction=high(above))


def main(bundle=None):
    ctx.prec,ctx.threads=112,2
    start,cpu=perf_counter(),process_time()
    y=load_public(bundle)
    args=load_geometry()
    geometries={}
    receipts=[]
    out=dict(scope='SRAP22 original full-ring 10% observations; no fit or access to any fine truth',
             radius='stored delta + 1e-6 ||Y||_F, outward rounded',
             native_dimension=96**2,joint_channel_shape=[20,10],precision_bits=ctx.prec,
             objects=[row for _,row in y],profiles=receipts,setup=times(start,cpu),sealed=False)
    for trial in choices():
        t,c=perf_counter(),process_time()
        side=trial['block_side']
        if side not in geometries:
            geometries[side]=Geometry(*args,block_side=side)
        g=geometries[side]
        p=saved_operator(g,trial)
        nu=trial['nu']
        prior=g.cost.copy()
        ngap=None
        negative={}
        if nu:
            tt,tc=perf_counter(),process_time()
            p.modes(nu)
            negative_upper,ngap=p.ensure_negative(nu,return_gap=True)
            negative=dict(negative_gap_lower=ngap,negative_congruence_upper=negative_upper,
                          **times(tt,tc))
        m,defect,work=exact_common_response(p,nu)
        stones=[bound_source(v,p,nu,trial['positive_gap_lower'],ngap,m,defect,row)
                for v,row in y]
        receipts.append(dict(**trial,negative_proof=negative,shared_stationary_work=work,
                             operations=dict(g.cost-prior),rows=stones,**times(t,c)))
        out.update(times(start,cpu))
        save_json(HERE/'online.json',out)
        print(trial['label'],sum(r['joint_status']=='excluded' for r in stones),
              sum(r['coherent_status']=='excluded' for r in stones),
              'strongest ratios',max(r['joint_distance_upper']/r['padded_radius'] for r in stones),flush=True)
    out.update(sealed=True,**times(start,cpu))
    save_json(HERE/'online.json',out)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--public-bundle', type=Path,
                     help='Use the portable pack of 22 observations and original radii.')
    args=cli.parse_args()
    with threadpool_limits(1):
        main(args.public_bundle)
