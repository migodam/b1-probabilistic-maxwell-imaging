"""Select and seal future actions from public data, with outward support proofs."""
from dataclasses import dataclass
from itertools import combinations
from time import time

from flint import acb, acb_mat, arb
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from model import (ANGLES, BATCHES, CONTRACTS, HERE, OLD, POLICIES, REGISTRY, WEIGHTS,
                   adj, ball, clocks, diag, elapsed, energy, eye, forward_jacobian,
                   geometry, hash_arrays, herm, high, hjoin, inertia, init_arithmetic,
                   low, mid, minor, noise_radius, norm, norm2, radius_squared, read_json, require_p0,
                   selection_matrix, spectrum, star, vjoin, write_json)


@dataclass
class Border:
    weight: float
    negative: int
    q: object
    receipt: dict


@dataclass
class Envelope:
    weight: float
    center: object
    gram: object
    receive: object
    inverse: object
    receipt: dict


def form_border(g,lo,hi,weight):
    started=clocks()
    intervals=list(zip([arb(float(x)) for x in lo],[arb(float(x)) for x in hi],strict=True))
    m=diag([(x+y)/2 for x,y in intervals])
    h2=diag([((y-x)/2)**2 for x,y in intervals])
    w=arb(float(weight))
    U=eye(len(lo))-m*g.ex_g
    B,S,G=g.ex_b,g.ex_s,g.ex_g
    K=star(S)*S+w*(star(U)*U-star(G)*h2*G)
    witness,_,_=inertia(K)
    if not witness['full_count_verified']:
        return None,dict(weight=weight,stage='K_INERTIA_UNRESOLVED',**witness,**elapsed(started))
    NB=w*(star(U)*m+star(G)*h2)*B
    RB=hjoin(NB,star(S))
    n,p=B.ncols(),S.nrows()
    D=vjoin(hjoin(w*star(B)*(m*m-h2)*B,acb_mat(n,p)),
            hjoin(acb_mat(p,n),eye(p)))
    R=K.solve(RB)
    defect=norm(K*R-RB)
    Q=D-star(RB)*R
    row=dict(weight=weight,stage='K_FULL_COUNT_AND_SCHUR_ENCLOSED',**witness,
             interval_solve_residual_upper=high(defect),**elapsed(started))
    return Border(weight,witness['negative'],Q,row),row


def source_pair(ids,amplitude,y):
    return vjoin(selection_matrix(ids,amplitude),ball(y))


def joint_tail(border,z):
    H=star(z)*border.q*z
    lo,hi,err=spectrum(H)
    num=max(0,H.nrows()-border.negative)
    lower_sq=sum((max(arb(0),e.lower()) for e in lo[:num]),arb(0)).lower()
    upper_sq=sum((max(arb(0),e.upper()) for e in hi[:num]),arb(0)).upper()
    return lower_sq,upper_sq,high(err)


def envelope(border,g,z,nr):
    started=clocks()
    d,p=len(ANGLES),g.s.shape[0]
    Q=border.q
    Po=star(z)*Q*z-radius_squared(len(OLD),nr)*eye(len(OLD))
    po,_,_=inertia(Po)
    status=dict(weight=border.weight,negative_state=border.negative,old_Inertia=po,
                old_positive=po['positive'],base='positive' if border.negative==0 else 'indefinite',
                saturated=False,receiver_positive=False)
    if not po['full_count_verified']:
        return None,dict(status,stage='OLD_INERTIA_UNRESOLVED',**elapsed(started))
    if po['positive']!=border.negative:
        relation='INDEX_DEFICIT' if po['positive']<border.negative else 'OLD_ALREADY_EXCLUDED'
        return None,dict(status,stage=relation,missing_old_index=max(0,border.negative-po['positive']),
                         **elapsed(started))
    status['saturated']=True
    Wq=Q*z
    qc=Q-Wq*Po.solve(star(Wq))
    Ry=minor(qc,range(d,d+p),range(d,d+p))
    Ry_rng,rylo,_=inertia(Ry)
    if Ry_rng['positive']!=p:
        return None,dict(status,stage='CONDITIONAL_RECEIVER_NOT_PROVED_POSITIVE',
                         receiver_inertia=Ry_rng,**elapsed(started))
    status['receiver_positive']=True
    off=minor(qc,range(d,d+p),range(d))
    ca=-Ry.solve(off)
    th=-star(off)*ca-minor(qc,range(d),range(d))
    thlo,_,terr=spectrum(th)
    # An outward isotropic repair makes the carried source Gram positive.
    # If roundoff already establishes positivity, no shift is applied.
    shift=max(arb(0),(-thlo[0]+arb("1e-24")).upper())
    th=th+shift*eye(d)
    RI=Ry.solve(eye(p))
    traceh=sum(float(th[i,i].real.mid()) for i in range(d))
    tracer=sum(float(RI[i,i].real.mid()) for i in range(p))
    info=dict(status,stage='BOUNDED_CONDITIONAL_OPERATOR',receiver_inertia=Ry_rng,
              receiver_gap_lower=low(rylo[0]),source_psd_shift=high(shift),
              source_spectral_error_upper=high(terr),
              average_width_proxy=float(np.sqrt(max(0.,traceh*tracer/d))),**elapsed(started))
    return Envelope(border.weight,ca,th,Ry,RI,info),info


def best_option(g, z, lo, hi, nr):
    borders=[]
    ok=[]
    rates=[]
    held=[]
    for w in WEIGHTS:
        b,rec=form_border(g,lo,hi,w)
        if b is None:
            rates.append(rec)
            continue
        borders.append(b)
        v,u,er=joint_tail(b,z)
        held.append((v,u,b))
        obj,rate=envelope(b,g,z,nr)
        rate['K_proof']=rec
        rate['old_joint_distance_lower']=max(0.,low(v.sqrt()))
        rate['old_joint_distance_upper']=high(u.sqrt())
        rates.append(rate)
        if obj:
            ok.append(obj)
    winner=min(ok,key=lambda e:e.receipt['average_width_proxy'],default=None)
    if held:
        v,_,best=max(held,key=lambda e:float(e[0]))
        old=dict(excluded=(v>radius_squared(len(OLD),nr)),best_distance_lower=max(0.,low(v.sqrt())),
                 weight=best.weight)
    else:
        old=dict(excluded=False,best_distance_lower=0.,weight=None)
    return borders,winner,rates,old


def sqrt_num(a, delta):
    vals,vec=eigh(herm(a),check_finite=False)
    return (vec*np.sqrt(np.maximum(vals,delta)))@adj(vec)


def prospective_numeric(option,ids,per_beam):
    c=mid(option.center)[:,list(ids)]*np.sqrt(per_beam)
    t=mid(option.gram)[np.ix_(ids,ids)]*per_beam
    return c,t,mid(option.inverse)


def real_trace(a):
    return sum((a[i,i].real for i in range(a.nrows())),arb(0))


def trace_support(H,C):
    """A precise enough, independently checked T-factor bound on the nuclear support."""
    hn,cn=mid(H),mid(C)
    floor=max(1.e-20,np.linalg.norm(hn,2)*1.e-12)
    h=sqrt_num(hn,floor)
    ih=np.linalg.solve(h,np.eye(h.shape[0]))
    test=herm(ih@sqrt_num(h@cn@h,max(1.e-30,np.linalg.norm(h@cn@h,2)*1.e-12))@ih)
    tb=ball(test)
    tlo,_,te=spectrum(tb)
    if not tlo[0]>0:
        # Identity keeps the inequality sound; a wide bound need not separate.
        tb=eye(tb.nrows())
        tlo=[arb(1)]
    upper=(real_trace(H*tb)+real_trace(C*tb.solve(eye(tb.nrows()))))/2
    if upper.upper()<0:
        raise ArithmeticError('Computed negative upper support.')
    return upper.upper(),dict(test_positive_gap_lower=low(tlo[0]),test_eigensystem_error_upper=high(te))


def report_support(option,w,ids,amplitude,nr):
    wb=ball(w)
    cb=minor(option.center,range(w.shape[0]),ids)*amplitude
    c=minor(option.gram,ids,ids)*amplitude**2
    h=star(wb)*option.receive.solve(wb)
    diameter,check=trace_support(h,c)
    center=real_trace(star(wb)*cb)
    noise=radius_squared(len(ids),nr).sqrt()*norm(wb)
    pad=diameter+noise
    return dict(weight=option.weight,center_interval=[low(center),high(center)],
        clean_support_upper=high(diameter),new_noise_support_upper=high(noise),
        noisy_support_interval=[low(center-pad),high(center+pad)],
        shared_source_gram_real=mid(c).real.tolist(),shared_source_gram_imag=mid(c).imag.tolist(),**check)


def split_query(options,ids,amplitude,nr):
    masses=float(amplitude**2)
    a,b=[prospective_numeric(o,ids,masses) for o in options]
    diff=b[0]-a[0]
    scale=np.linalg.norm(diff)
    if scale<=1.e-16:
        return 0.,None
    w=diff/scale
    penalties=0.
    for _,c,inv in (a,b):
        c12=sqrt_num(c,0.)
        P=herm(c12@adj(w)@inv@w@c12)
        penalties+=np.sqrt(np.maximum(0.,eigh(P,eigvals_only=True,check_finite=False))).sum()
    return float(scale-penalties-2*float(radius_squared(len(ids),nr).sqrt())),w


def certify_pair(options,w,ids,amplitude,nr):
    checked=[report_support(e,w,ids,amplitude,nr) for e in options]
    a,b=[c['noisy_support_interval'] for c in checked]
    gap=low(arb(b[0])-arb(a[1]))
    return dict(noisy_separation_gap_lower=gap,strict_edge=(gap>0),supports=dict(zip(('box_a','box_b'),checked)),
                worst_survivors_bound=(1 if gap>0 else 2),
                support_query_index=None)


def angular_gap(ids):
    v=np.sort([ANGLES[i] for i in (*OLD,*ids)])
    g=np.diff(np.r_[v,v[0]+360])
    return float(-g@g)


def fisher_terms(g,centers,old_energy):
    d=centers[1]-centers[0]; d=d/np.linalg.norm(d)
    database=[]
    p=g.s.shape[0]
    for m in centers:
        _,J=forward_jacobian(g,m,tuple(range(len(ANGLES))))
        J=J.reshape(2,p,len(ANGLES),len(m))
        physical=np.einsum('rpan,rpam->anm',J,J)
        old=float(old_energy)*physical[list(OLD)].sum(axis=0)
        ridge=max(1.e-14,1.e-5*np.trace(old)/len(m))
        database.append((physical,old+ridge*np.eye(len(m)),ridge))
    return d,database


def fisher_score(ids, per_beam, d, database):
    # 1/(d*F^-1 d) is the information for d*x after profiling all
    # orthogonal material directions, under a disclosed weak ridge.
    return float(min(1/(d@np.linalg.solve(old+per_beam*phys[list(ids)].sum(axis=0),d))
                     for phys,old,_ in database))


def design_scene(scene,info,pub,old,contract,arrays,nr,audit=None):
    start=clocks()
    g=geometry(info['layout'])
    assert hash_arrays(g.g,g.s,g.b)==info['model_hash']
    old_power,_=energy(contract,1)
    z=source_pair(OLD,old_power.sqrt(),old)
    independent_candidates = {}
    for i in range(len(ANGLES)):
        if i in OLD:
            continue
        proposed=minor(g.ex_b,range(g.g.shape[0]),(*OLD,i))
        gap=low(spectrum(star(proposed)*proposed)[0][0])
        if gap>0:
            independent_candidates[i]=gap
    if len(independent_candidates)<max(BATCHES):
        raise ArithmeticError('Too few physical menu beams are certified independent of old beams.')
    lo,hi,centers=[pub[scene+'__'+name] for name in ('lo','hi','center')]
    if hash_arrays(lo,hi)!=info['box_pair_hash']:
        raise ArithmeticError('The preregistered material boxes changed.')
    options={}
    old_classes={}
    for label,l,h in zip(info['classes'],lo,hi,strict=True):
        _,e,rates,old_check=best_option(g,z,l,h,nr)
        options[label]=e
        old_classes[label]=dict(profiles=rates,old_joint=old_check,
                                chosen_envelope=e.receipt if e else None)
    # An overlapping broad box is diagnostic, never leads to a candidate removal.
    _,_,wide,_=best_option(g,z,np.zeros(g.g.shape[0]),np.full(g.g.shape[0],3.),nr)
    survivors=[c for c in info['classes'] if not old_classes[c]['old_joint']['excluded']]
    d,priors=fisher_terms(g,centers,old_power)
    allowed=tuple(independent_candidates)
    available=all(options[c] is not None for c in info['classes'])
    envelopes=list(options.values())
    out={}
    for size in BATCHES:
        subset_started=clocks()
        _,power=energy(contract,size)
        amp=power.sqrt()
        all_combos=list(combinations(allowed,size))
        if len(survivors)<2:
            out[str(size)]=dict(status='STOPPED_OLD_SINGLETON',actions={},actual_old=len(OLD),actual_new=0)
            continue
        vals=dict(uniform_gap=[],fisher_nuisance=[],numeric_envelope=[])
        query={}
        for j,ids in enumerate(all_combos):
            vals['uniform_gap'].append(angular_gap(ids))
            vals['fisher_nuisance'].append(fisher_score(ids,float(power),d,priors))
            if available:
                gap,w=split_query(envelopes,ids,amp,nr)
                vals['numeric_envelope'].append(gap)
                if w is not None:query[j]=w
            else:
                vals['numeric_envelope'].append(None)
        rng=np.random.default_rng(np.random.SeedSequence([sum(scene.encode()),len(contract),size,2097]))
        uniform=int(np.argmax(vals['uniform_gap']))
        fisher=int(np.argmax(vals['fisher_nuisance']))
        selected=dict(random=int(rng.integers(len(all_combos))),uniform_gap=uniform,fisher_nuisance=fisher)
        candidate_proofs={}
        positive_numeric=sorted(query,key=lambda j:(-vals['numeric_envelope'][j],j))
        # Full interval support for every numerically positive proposed edge,
        # plus the controlled policy actions so post-readout uses the same test.
        attempted=sorted(set(selected.values())|{j for j in positive_numeric if vals['numeric_envelope'][j]>0}
                         |set(positive_numeric[:1]))
        if available:
            for j in attempted:
                if j in query:
                    candidate_proofs[j]=certify_pair(envelopes,query[j],all_combos[j],amp,nr)
        edges=[j for j,o in candidate_proofs.items() if o['strict_edge']]
        if audit is not None:
            audited=[]
            for j,ids in enumerate(all_combos):
                proof=candidate_proofs.get(j)
                if proof:
                    key=f'{scene}__{contract}__{size}__{j}'
                    arrays[key]=query[j]
                    proof=dict(proof,support_query_index=key,support_query_hash=hash_arrays(query[j]))
                status=('NOT_INTERVAL_CHECKED' if proof is None else
                        'STRICT_EDGE' if proof['strict_edge'] else 'NO_STRICT_EDGE_IN_THIS_CHECK')
                source=minor(g.ex_b,range(g.g.shape[0]),(*OLD,*ids))
                source_gap=low(spectrum(star(source)*source)[0][0])
                if source_gap<=0:
                    raise ArithmeticError('No strict source Gram rank for an archived action.')
                audited.append(dict(source_indices=ids,status=status,verified_support=proof,
                                    all_old_and_batch_source_gram_floor=source_gap))
            audit[str(size)]=audited
        if edges:
            selected['conditional_ellipsoid']=max(edges,key=lambda j:(candidate_proofs[j]['noisy_separation_gap_lower'],-j))
        else:
            selected['conditional_ellipsoid']=uniform
        actions={}
        for policy,j in selected.items():
            ids=all_combos[j]
            pr=candidate_proofs.get(j)
            V=minor(g.ex_b,range(g.g.shape[0]),(*OLD,*ids))
            chosen_batch_floor=low(spectrum(star(V)*V)[0][0])
            if pr:
                key=f'{scene}__{contract}__{size}__{j}'
                arrays[key]=query[j]
                pr=dict(pr,support_query_index=key,support_query_hash=hash_arrays(query[j]))
            actions[policy]=dict(source_indices=ids,source_angles_degrees=[ANGLES[k] for k in ids],
                score_j=j,all_old_and_batch_source_gram_floor=chosen_batch_floor,source_energy_each_exact=str(power),source_energy_each=float(power),
                total_new_energy=float(power*size),old_total_energy=float(old_power*len(OLD)),
                radius_old=high(radius_squared(len(OLD),nr).sqrt()),radius_new=high(radius_squared(size,nr).sqrt()),
                radius_joint=high(radius_squared(len(OLD)+size,nr).sqrt()),decision_rule=(
                    'strict_class_separation' if policy=='conditional_ellipsoid' and j in edges else
                    'uniform_gap_fallback' if policy=='conditional_ellipsoid' else policy),
                verified_support=pr,guarantee='CONDITIONAL_TWO_BOX_EDGE' if pr and pr['strict_edge'] else 'NO_CERTIFIED_SEPARATION')
        out[str(size)]=dict(status='SEALED_BEFORE_ACQUISITION',actions=actions,
                finite_menu_size=len(all_combos),attempted_interval_actions=len(candidate_proofs),
                strict_separating_actions=len(edges),strict_separation_fraction_lower=len(edges)/len(all_combos),
                all_action_menu=[list(j) for j in all_combos],numerical_scores=vals,
                fisher_ridges=[float(k[2]) for k in priors],**elapsed(subset_started))
    return dict(scene=scene,layout=info['layout'],operator_model_hash=info['model_hash'],
        public_candidate_cover=info['box_pair_hash'],old_data_hash=hash_arrays(old),energy_contract=contract,
        matrix_noise_upper='r_old^2 I',physical_noise_radius_per_beam=nr,candidate_vs_old_source_gram_floors=independent_candidates,
        coverage_before_update=survivors,old_classes=old_classes,indefinite_stress_profiles=wide,
        batches=out,**elapsed(start))


def main():
    start=clocks()
    gate=require_p0()
    reg=read_json(REGISTRY)
    if (HERE/'proposal.json').exists():
        print('Prospective action receipt already frozen. Clear only in an isolated replay.')
        return
    if (HERE/'acquired').exists():
        raise RuntimeError('Future acquisition already exists before proposal freeze.')
    all_rows=[]
    arrays={}
    with np.load(HERE/'public_model.npz',allow_pickle=False) as pub:
        for scene,info in reg['scenes'].items():
            for contract in CONTRACTS:
                old_file=HERE/'old'/f'{scene}__{contract}.npz'
                with np.load(old_file,allow_pickle=False) as f:
                    if set(f.files)!={'Y'}:
                        raise ArithmeticError('Old design input contains non-public arrays.')
                    row=design_scene(scene,info,pub,f['Y'],contract,arrays,noise_radius(reg,scene))
                all_rows.append(row)
                print(scene,contract,[(size,x['strict_separating_actions']) for size,x in row['batches'].items()
                                       if 'strict_separating_actions' in x],flush=True)
    np.savez_compressed(HERE/'proposal_arrays.npz',**arrays)
    write_json(HERE/'proposal.json',dict(status='SEALED_READ_ONLY_BEFORE_ANY_NEW_DATA',
       protocol='two-class exhaustive fixed basis angle menu, separate contracts, equal policy spend',
       preregistration='registration.json',source_arrays='proposal_arrays.npz',gate=gate,
       wall_time_unix=time(),precision_bits=128,source_rows=all_rows,**elapsed(start)))
    print('Complete pre-acquisition proposal sealed.',flush=True)


if __name__=='__main__':
    init_arithmetic()
    with threadpool_limits(1):
        main()
