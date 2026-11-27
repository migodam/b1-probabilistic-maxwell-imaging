"""Enclose the minimum exterior edit needed by each opposite core for the actual acquired data."""
from fractions import Fraction
import json

import numpy as np
from flint import arb
from threadpoolctl import threadpool_limits

from run import HERE, SIDE, controller, read_pair
from coherent_modes import AuditedState, DEN, LAYOUT, VerifiedWave, down, elapsed, fro, now, rd, released_operator, up
from frequency_classes import directional
from frequency_classes_open import empty_network_coefficients, uniform_wave_error
from verify_dual import upper


def scalar_cover(state,jac,port,empty,coeff,width,sign,roi):
    """One noisy halfspace endpoint for the exact frozen receiver and fixed core."""
    j,jerror,u,e=jac
    lo=-state.x.copy()
    hi=np.nextafter(4.5-state.x,np.inf)
    hi[state.x==4.5]=0.
    lo[roi]=hi[roi]=0.
    bound,qn=directional(np.zeros_like(port),port,sign*j,jerror,u,e,state.u_error,
                          state.e_error,state.k_total,width,'max',lo,hi,0.)
    if bound['nonlinear'] is None:
        return None
    try:
        transfer=uniform_wave_error(state,state.k_total,width,coeff)['response_error_fro_upper']
    except ArithmeticError:
        return None
    k=up(width*state.k_total)
    if k>=1:
        return None
    edit=up((fro(state.u)+state.u_error)*width*(fro(state.e)+state.e_error)/(1-k))
    fn=up(fro(state.p)+state.response_error+edit+fro(empty)+coeff['empty_total_Fro_error_upper']+transfer)
    values=state.p+empty
    v=np.r_[values.real.ravel(),values.imag.ravel()]
    score=float(qn@v)
    guard=up(state.response_error+coeff['empty_total_Fro_error_upper']+rd(132)*(fro(state.p)+fro(empty))+transfer)
    direrr=up(fro(port-qn)+rd(3)*(fro(port)+fro(qn)))
    penalty=up(bound['safe_support']+bound['nonlinear']+bound['derivative_term_error']+direrr*fn+guard)
    root=8**.5
    cerr=up(float(abs(arb(8).sqrt()-arb(root)).upper()))
    error=up(root*penalty+cerr*fn+rd(5)*(abs(score)*root+penalty*root))
    endpoint=down(score*root-error) if sign==1 else -down(-score*root-error)
    return Fraction(endpoint)


def main():
    clock=now()
    proposal,dec,data,_,_,_,(tx,rx,weights)=controller()
    amp=Fraction(proposal['amplitude'])
    refs,_,roi,r,_=read_pair()
    # Integer common denominator exactly reproduces both analog weights and radiated amplitude split.
    assert tx.shape==rx.shape==(2,8,1)
    qr=(rx[0]*weights**2)@tx[0].T+(rx[1]*weights**2)@tx[1].T
    qi=(rx[1]*weights**2)@tx[0].T-(rx[0]*weights**2)@tx[1].T
    q=np.r_[qr.ravel(),qi.ravel()]/DEN**4
    _,op=released_operator(SIDE,460)
    wave=VerifiedWave(op)
    coeff,empty=empty_network_coefficients(op)
    replay=json.loads((HERE/'dual_verification.json').read_text())
    rows={}
    for name,measure in data.items():
        index=1 if dec[name]['label']=='high' else 0
        sign=1 if index==0 else -1
        x=refs[index]
        y=Fraction(int(weights[0]),DEN)*Fraction(float(measure[1][0].real))
        st=AuditedState(wave,op['S'],op['E0'],x)
        jac=st.jacobian(LAYOUT)
        history=[]
        safe,other=0,4096
        for _ in range(12):
            k=(safe+other)//2
            b=k/16384
            endpoint=scalar_cover(st,jac,q,empty,coeff,b,sign,roi)
            margin=None if endpoint is None else sign*(amp*endpoint-y)-Fraction(r)
            passed=margin is not None and margin>0
            history.append(dict(width=b,passed=passed,strict_scalar_gap_after_noise=None if margin is None else -upper(-margin)))
            if passed:safe=k
            else:other=k
        assert safe>0
        att=next(k for k,v in replay['samples'][name]['materials'].items()
                 if k!='public_reference' and v['both_balls'])
        with np.load(HERE/(att+'.npz'),allow_pickle=False) as a:
            chosen=a['material']
        assert np.array_equal(chosen[roi],x[roi])
        delta=[abs(Fraction(f)-Fraction(g)) for f,g in zip(chosen[~roi],x[~roi])]
        area=Fraction(4,SIDE*SIDE)
        l1=area*sum(delta,Fraction())
        l2sq=area*sum((v*v for v in delta),Fraction())
        rows[name]=dict(opposite_core=float(x[roi][0]),excluded_exterior_Linf_through=safe/16384,
                        verified_witness_Linf_upper=upper(max(delta)),area_L1_witness_upper_m2=upper(l1),
                        area_L2_squared_witness_upper_m2=upper(l2sq),witness=att,
                        bracket_interpretation='The lower endpoint is excluded by this actual scalar even with the old data omitted. The upper endpoint has a common two-ball witness.',history=history)
    (HERE/'background_budget.json').write_text(json.dumps(dict(side=96,materials_box=[0,4.5],
        core_pixels=int(sum(roi)),observed_added='stored one scalar at 460 MHz',rows=rows,**elapsed(clock)),indent=2)+'\n')
    for name,r in rows.items():
        print(name,r['excluded_exterior_Linf_through'],r['verified_witness_Linf_upper'],flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
