"""Factor a public class-separating port observable into physical, fluence-limited shots."""
import argparse
from pathlib import Path
import json
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
from frequency_classes import (HERE, LAYOUT, AuditedState, VerifiedWave, check_direction,
                              explicit_difference, down, elapsed, fro, material, now, rd, released_operator, up)
from frequency_classes_open import empty_network_coefficients, uniform_wave_error
from matched_amplitude import read_pair
from flint import arb

TARGET=HERE/'fluence_modes'
DEN=4096


def safe_vectors(columns):
    ticks=np.stack([np.trunc(columns.real*DEN),np.trunc(columns.imag*DEN)]).astype(np.int64)
    if np.any(np.sum(ticks*ticks,axis=(0,1)) > DEN**2):
        raise ArithmeticError('Unnormalized power wave')
    return ticks


def physical(observable,count):
    left,s,vh=svd(observable,full_matrices=False)
    left=safe_vectors(left[:,:count])
    right=safe_vectors(vh.conj().T[:,:count])
    alpha=np.floor(DEN*np.sqrt(s[:count]/s[:count].sum())).astype(np.int64)
    if int(alpha@alpha)>DEN**2:
        raise ArithmeticError('Source and receiver allocations exceed the total fluence')
    # Exact integer products: entries have modulus at most 2^48.
    qr=(left[0]*alpha**2)@right[0].T + (left[1]*alpha**2)@right[1].T
    qi=(left[1]*alpha**2)@right[0].T - (left[0]*alpha**2)@right[1].T
    if max(np.max(abs(qr)),np.max(abs(qi)))>DEN**4:
        raise ArithmeticError('Observable bound violated')
    q=np.r_[qr.ravel(),qi.ravel()]/DEN**4
    return q,dict(rx_ticks=left,tx_ticks=right,amp_ticks=alpha),s


def run(side):
    clock=now()
    if side==96:
        mats,_,_,r,_=read_pair()
        data_file=HERE/'matched_amplitude.npz'
        rec=json.loads((HERE/'matched_amplitude.json').read_text())
        item=next(x for x in rec['outputs']['460']['budgets'] if x['beta']==.01)
    else:
        _,_,r,_,_,_,plus,minus=material(32)
        mats=[plus,minus]
        data_file=HERE/'frequency_classes_32.npz'
        rec=json.loads((HERE/'frequency_classes_32.json').read_text())
        item=next(x for x in rec['frequency']['460']['budget'] if x['beta']==.01)
    with np.load(data_file,allow_pickle=False) as stored: lam=stored[item['lambda_key']]
    dual=lam[:64].reshape(8,8)+1j*lam[64:].reshape(8,8)
    _,op=released_operator(side,460)
    wave=VerifiedWave(op)
    states=[AuditedState(wave,op['S'],op['E0'],x) for x in mats]
    jacs=[p.jacobian(LAYOUT) for p in states]
    coeff,empty=empty_network_coefficients(op)
    transfers=[uniform_wave_error(p,p.k_total,.01,coeff)['response_error_fro_upper'] for p in states]
    d,ed=explicit_difference(states[0].p,states[1].p)
    ed=up(ed+sum(p.response_error for p in states))
    # Bounds on the true complete total wave at every allowed material, not only the centers.
    fullnorm=[]
    for st,transfer in zip(states,transfers):
        q=up(.01*st.k_total)
        if q>=1:raise ArithmeticError('No uniform inverse')
        edit=up((fro(st.u)+st.u_error)*.01*(fro(st.e)+st.e_error)/(1-q))
        fullnorm.append(up(fro(st.p)+st.response_error+edit+fro(empty)+coeff['empty_total_Fro_error_upper']+transfer))
    center=d[:64].reshape(8,8)+1j*d[64:].reshape(8,8)
    rows=[]
    arrays={}
    root=8**.5
    sqrt_gap=up(float(abs(arb(8).sqrt()-arb(root)).upper()))
    for algo,reference in [('center',center),('class',dual)]:
        for m in [1,2,3,4,8]:
            q,factors,s=physical(reference,m)
            _,qn,projection,terms=check_direction(d,ed,q,states,jacs,.01)
            if not all(x['nonlinear'] is not None for x in terms):
                raise ArithmeticError('Remainder undefined')
            qdiff=up(fro(q-qn)+rd(3)*(fro(q)+fro(qn)))
            scores=[]
            for j,(p,term,transfer,fn) in enumerate(zip(states,terms,transfers,fullnorm)):
                values=p.p+empty
                realval=np.r_[values.real.ravel(),values.imag.ravel()]
                a=float(qn@realval)
                center_error=up(p.response_error+coeff['empty_total_Fro_error_upper']+
                                rd(132)*(fro(p.p)+fro(empty))+transfer)
                penalty=up(term['safe_support']+term['nonlinear']+term['derivative_term_error']+qdiff*fn+center_error)
                # Include the effect of using sqrt(8) at the physical transmitter.
                error=up(root*penalty+sqrt_gap*fn+rd(5)*(abs(a)*root+penalty*root))
                scores.append(down(a*root-error) if j==0 else -down(-a*root-error))
            sep=max(0.,down((scores[0]-scores[1]-rd(4)*(abs(scores[0])+abs(scores[1])))/(2*r)))
            key=f'{algo}_{m}'
            arrays.update({key+'_'+k:v for k,v in factors.items()})
            rows.append(dict(key=key,shots=m,complex_readings_after_analog_combining=m,
                             source_design=algo,fullclass_lower_over_2r=sep,
                             sufficient_energy_strictly_above=up(1/(sep*sep)) if sep>0 else None,
                             unit_high_clean_score_lower=scores[0],unit_low_clean_score_upper=scores[1],
                             feasible_unit_energy=sep>1,nominal_l2_spectrum=s.tolist(),
                             q_rounding_guard=qdiff,analog_phase_grid_denominator=DEN,
                             tx_fluence_upper=float(int(factors['amp_ticks']@factors['amp_ticks'])/DEN**2),
                             two_physical_endpoint_terms=terms))
    TARGET.mkdir(exist_ok=True)
    np.savez_compressed(TARGET/f'actions_{side}.npz',**arrays)
    result=dict(side=side,radius=r,freq_MHz=460,beta=.01,core_gap=2.,hardware='8 transmit and 8 receive ports. One separable coherent sample each pulse, ideal calibrated passive combining at receiver plus digital weighting. One joint absolute Frobenius radius per action. Analytical whole-class check using only public specifications.',
                guard_exactness='Complex coefficients are products/sums of dyadic 12-bit tick weights. integer entries <=2^48; exact in binary64. Error of floating normalization and open-circuit transfer audited per direction.',
                rows=rows,**elapsed(clock))
    (TARGET/f'classes_{side}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(side,[(v['key'],round(v['fullclass_lower_over_2r'],5)) for v in rows],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('sides',type=int,choices=(32,96),nargs='+')
    opts=parser.parse_args()
    with threadpool_limits(1):
        for count in opts.sides:run(count)
