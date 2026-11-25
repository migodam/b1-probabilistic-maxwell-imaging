"""Independent strict original-operator verification of paired quadrics."""
import argparse
import json

import numpy as np
from flint import acb, acb_mat, arb, ctx
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, elapsed, geometry, save
from interval_certify import (ExactBox, ProofGrid, ball, column, conditional_q, endpoints,
                              eye, hi, joined, lo, midpoint, norm_upper,
                              positive_lower, tr, trust_direction)


def dimblock(q,rows,cols):
    rows,cols=list(rows),list(cols)
    return acb_mat(len(rows),len(cols),[q[i,j] for i in rows for j in cols])


def exact_form(box, form):
    if form['nu']==0:
        base,l,c=box.optical_base(form['mu'])
        return c-tr(l)*base.solve(l),dict(state_negative=0,normal_base_lower=lo(arb(form['mu'])*(1-box.eta_reverse**2)))
    return box.optical_q(form,base=box.optical_base(form['mu']))


def projected_separation(q, action, gamma, sigma, total):
    m=action.nrows()
    p=q[0].nrows()-m
    recv=range(m,m+p)
    parts=[]
    for f,weight in zip(q,gamma,strict=True):
        w=acb(float(weight))
        cross=w*dimblock(f,recv,range(m))*action
        r=eye(p)+w*dimblock(f,recv,recv)
        c=w*(tr(action)*dimblock(f,range(m),range(m))*action)[0,0]
        parts.append((r,cross,c))
    one,two=parts
    r=acb_mat(2*p,2*p,joined(one[0],-eye(p)).entries()+joined(-eye(p),two[0]).entries())
    b=acb_mat(2*p,1,one[1].entries()+two[1].entries())
    const=one[2]+two[2]
    strict_floor=positive_lower(r)
    if strict_floor is None:
        return dict(receive_pair_verified_positive=False,strict_separation=False)
    obj=(const-(tr(b)*r.solve(b))[0,0]).real
    threshold=arb(4)*arb(sigma)**2
    qualified=bool(obj>threshold)
    nn=one[1].nrows()+two[1].nrows()
    bordered=acb_mat(nn+1,nn+1,joined(r,b).entries()+tr(b).entries()+[const-acb(threshold)])
    border_floor=positive_lower(bordered) if qualified else None
    return dict(receive_pair_verified_positive=True,
                receiving_joint_gap_lower=lo(strict_floor),
                distance_squared_enclosure=endpoints(obj),
                lower_distance_over_radius=lo(obj.lower().sqrt()/arb(sigma)) if obj>0 else 0.,
                sufficient_fluence_for_fixed_direction_upper=hi(total*threshold/obj) if obj>0 else None,
                certificate_strict_margin_squared_lower=lo(obj-threshold),
                bordered_direct_positive_floor_lower=None if border_floor is None else lo(border_floor),
                strict_separation=qualified)


def matrix_values(q):
    mid=midpoint(q)
    radius=np.array([hi(abs(z-acb(float(v.real),float(v.imag))).upper()) for z,v in zip(q.entries(),mid.flat,strict=True)])
    return mid,radius.reshape(mid.shape)


def checked_row(models,meta,arrays,main_arrays,index,old,new,gram,numeric,bank,truth):
    tick=clocks()
    xold=old[index+'__Y']
    actions=arrays[index+'__actions']
    V=ball(actions)
    power=tr(V)*gram*V
    target=arb(float(meta['fluence']))
    v=[acb((target/power[j,j].real).sqrt())*column(V,j) for j in range(actions.shape[1])]
    d0=np.eye(actions.shape[0],dtype=complex)[:,[0]]
    zold=acb_mat(actions.shape[0]+len(xold),1,ball(d0).entries()+ball(xold).entries())
    sigma=meta['new_radius']
    ro=meta['old_radius']
    policies=meta['policies']
    outcomes={}
    arr={}
    witnesses={}
    midpoint_tests={}
    receiver_center=[getattr(model,'center_y',model.Y) for model in models]
    noiseless_olds=[hi(norm_upper(y*ball(d0)-ball(xold))/arb(ro)) for y in receiver_center]
    original=np.array(main_arrays[index+'__actions'])
    original_names=('uniform_angle','random_angle','fisher_angle','direct_coherent','source_eigenbeam','source_shaped')
    standard=ball(original)
    pow_orig=tr(standard)*gram*standard
    standard_arms={name:acb((target/pow_orig[j,j].real).sqrt())*column(standard,j)
                   for j,name in enumerate(original_names)}
    all_arms=dict(zip(policies,v,strict=True))
    all_arms.update({'original_'+n:x for n,x in standard_arms.items()})
    for name,action in all_arms.items():
        centers=[y*action for y in receiver_center]
        p=ball(midpoint(acb(arb(1)/2)*(centers[0]+centers[1])))
        reaches=[hi(norm_upper(y-p)/arb(sigma)) for y in centers]
        midpoint_tests[name]=dict(old_center_residual_over_radius_upper=noiseless_olds,
                                  new_center_residual_over_radius_upper=reaches,
                                  both_classes_have_same_legal_old_and_new_reading=bool(max(noiseless_olds+reaches)<1.0))
        arr[name+'__common_new_witness']=midpoint(p)
    if new is not None:
        vectors=np.column_stack([d0,actions])
        exact=joined(ball(d0),*v)
        observations=np.column_stack([xold,new[index+'__Y']])
        radii=[ro]+[sigma]*len(v)
        fg=bank[index+'__fg']
        w=bank[index+'__widths']
        residual_rows=[]
        for j,box in enumerate(models):
            center=box.center_np
            alt=w*np.where(fg,-.46,.56)
            phys=np.stack((center,center+alt,center-alt))
            residual=[]
            for point in phys:
                val,_,_,_=ExactBox.physical_candidate(box,numeric,point,vectors,observations,radii,exact_vectors=exact)
                residual.append(val)
            arr[f'{j}__candidate_materials']=phys
            arr[f'{j}__candidate_residual_over_radius_upper']=np.array(residual)
            residual_rows.append(np.array(residual))
        for j,name in enumerate(policies):
            witnesses[name]=[int(np.count_nonzero((rr[:,0]<=1)&(rr[:,j+1]<=1))) for rr in residual_rows]

    for label, candidate in meta['cases'].items():
        if not candidate['found']:
            continue
        if any(w<0 for w in candidate['dual_weights']):
            raise ArithmeticError('Distance multipliers must be nonnegative.')
        t=clocks()
        raw=[]
        scripts=[]
        for j,(model,spec) in enumerate(zip(models,candidate['forms'],strict=True)):
            raw_q,nu=exact_form(model,spec) if isinstance(model,ExactBox) else model.optical_form(spec)
            if raw_q is None:
                scripts.append(nu)
                continue
            q,pivot,k=conditional_q(raw_q,zold,ro,actions.shape[0])
            saturated=k==nu['state_negative']==spec['nu']
            margin=abs(pivot).lower()
            scripts.append(dict(**nu,old_pivot_interval=endpoints(pivot),old_positive=k,
                                old_saturation=saturated,old_pivot_absolute_margin_lower=lo(margin)))
            if saturated:
                raw.append(q)
                arr[f'{label}__{j}__q_midpoint'],arr[f'{label}__{j}__q_entry_radii']=matrix_values(q)
        results={}
        threshold_tests={}
        if len(raw)==2:
            for j,name in enumerate(policies):
                results[name]=projected_separation(raw,v[j],candidate['dual_weights'],sigma,target)
            for name,a in standard_arms.items():
                results['original_'+name]=projected_separation(raw,a,candidate['dual_weights'],sigma,target)
            if new is not None:
                if not np.array_equal(new[index+'__V'],actions):
                    raise ArithmeticError('Acquired actions differ from prospective frozen proposal.')
                Y=ball(new[index+'__Y'])
                for j,name in enumerate(policies):
                    cls=[]
                    for q in raw:
                        trust=trust_direction(q,actions.shape[0],v[j],column(Y,j),sigma)
                        cls.append(dict(excludes=bool(trust is not None and trust[0]>0),
                                        strict_lower_bound=None if trust is None else trust[0],
                                        multiplier=None if trust is None else trust[1]))
                    for k,decision in enumerate(cls):
                        n_wit=witnesses[name][k]
                        if decision['excludes'] and n_wit:
                            raise ArithmeticError('Physical witness contradicts strict quadratic exclusion.')
                        decision['physical_witness_count']=n_wit
                        decision['state']='EXCLUDED' if decision['excludes'] else ('WITNESSED' if n_wit else 'UNRESOLVED')
                    threshold_tests[name]=cls
        outcomes[label]=dict(forms=scripts,certified=results,observed_wrong_class_checks=threshold_tests,
                             **elapsed(t))
    audited={}
    if truth is not None:
        true_label=int(truth[index+'__label'])
        actual=truth[index+'__chi']
        bounds,_,_,_=ExactBox.physical_candidate(models[true_label],numeric,actual,vectors,observations,radii,exact_vectors=exact)
        if max(bounds)>=1:
            raise ArithmeticError('Actual simulated material lies outside the registered data ball.')
        audited=dict(offline_generator_label=true_label,offline_generator_all_actions_residual_upper=bounds)
    return dict(outcomes=outcomes,exact_common_fluence=float(meta['fluence']),
                verification_after_frozen_data=audited,per_action_common_reading_witness=midpoint_tests,
                **elapsed(tick)),arr


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('study')
    parser.add_argument('indices',type=int,nargs='*')
    parser.add_argument('--forecast-only',action='store_true')
    parser.add_argument('--suffix',default='')
    args=parser.parse_args()
    ctx.prec=128
    ctx.threads=3
    side=int(args.study.split('_')[-1])
    src=HERE.parent/'pair'/args.study
    old=np.load(src/'old_readings.npz')
    prop=json.loads((src/'proposal.json').read_text())
    arrays=np.load(src/'proposal_arrays.npz')
    main_arrays=np.load(HERE/prop['physical_source']/'proposal_arrays.npz')
    bank=np.load(src/'public_classes.npz')
    new=np.load(src/'future_readings.npz') if (not args.forecast_only and (src/'future_readings.npz').exists()) else None
    truth=np.load(src/'future_offline.npz') if new is not None else None
    stem=HERE/f'paired_replay_{args.study}{args.suffix}'
    tick=clocks()
    geo=geometry(side)
    pgrid=ProofGrid(geo)
    p=tr(ball(geo.B))*ball(geo.B)
    p=acb(1/p[0,0].real)*p
    green_err_per_entry=None
    if side==96:
        from native_paired import FastOpticalBox, green_kernel_roundoff
        green_err_per_entry=green_kernel_roundoff(pgrid)
    rows={}
    data={}
    for i,meta in prop['rows'].items():
        if args.indices and int(i) not in args.indices:
            continue
        t=clocks()
        cells=bank[i+'__cells']
        ctr=bank[i+'__centers']
        widths=bank[i+'__widths']
        numeric,exact=pgrid.restricted(cells)
        if side==96:
            models=[FastOpticalBox(exact,numeric,center,widths,green_err_per_entry*len(cells)) for center in ctr]
        else:
            models=[ExactBox(exact,numeric,center,widths) for center in ctr]
        setup=elapsed(t)
        row,stored=checked_row(models,meta,arrays,main_arrays,i,old,new,p,numeric,bank,truth)
        rows[i]=dict(row,box_and_original_operator_setup=setup,full_row_cost=elapsed(t))
        data.update({i+'__'+k:v for k,v in stored.items()})
        atomic_write(stem.with_suffix('.json'),
                     dict(side=side,stage='VERIFIED_ORIGINAL_OPERATOR',has_fresh_readings=new is not None,
                          precision_bits=ctx.prec,
                          method='verified original-operator state/congruence, exact common fluence, old saturation, paired distance dual, independent trust-region new observation',
                          rows=rows,**elapsed(tick)))
        print(i,[(label,[v.get('state_negative') for v in r['forms']],
                  {name:round(z.get('lower_distance_over_radius',0),5) for name,z in r['certified'].items() if 'original_' not in name})
                 for label,r in row['outcomes'].items()],flush=True)
    save(stem.with_suffix('.npz'),**data)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
