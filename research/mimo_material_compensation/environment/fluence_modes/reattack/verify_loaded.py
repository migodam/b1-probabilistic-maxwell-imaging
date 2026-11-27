"""Independent, loaded-24-port proof of both actual data balls, without the inverse-search code."""
from pathlib import Path
import json
import sys

import numpy as np
from flint import acb, arb, ctx
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
ENV=HERE.parents[1]
sys.path.insert(0,str(ENV))
from frequency_classes import AuditedState, VerifiedWave, elapsed, fro, now, rd, released_operator, up
from frequency_classes_open import empty_network_coefficients, uniform_wave_error


def ball_complex(z):
    return acb(float(z.real),float(z.imag))


def check():
    start=now()
    project=ENV.parents[1]/'z_alg_optical_inertia'/'blind_uq_mimo_theory'
    frozen=HERE.parent/'n96'
    proposal=json.loads((frozen/'proposal.json').read_text())
    radius=proposal['detector_radius']
    with np.load(project/'ports/public/cases/n96_c01.npz',allow_pickle=False) as case:
        assert float(case['radius'])==radius
    with np.load(ENV/'wave_96_2.0_1.npz',allow_pickle=False) as src:
        assert float(src['r'])==radius
        footprint=src['ix']
        roi=np.sort(footprint[src['group']==0])
        references={}
        for key in ('plus','minus'):
            x=np.zeros(96**2)
            x[footprint]=src[key]
            references[int(x[roi[0]])]=x
    assert len(roi)==164 and set(references)=={0,2}
    with np.load(frozen/'controls.npz',allow_pickle=False) as c:
        old=c['old_observed']
        tr,rr,m=[c[k] for k in ('tx_ticks','rx_ticks','amp_ticks')]
    den=4096
    assert old.shape==(8,8) and tr.shape==rr.shape==(2,8,1) and m.shape==(1,)
    for ticks in (tr,rr,m):
        assert sum(int(v)**2 for v in ticks.ravel())<=den**2
    vectors=[]
    for ticks in (tr,rr):
        vectors.append([acb(int(a),int(b))/den for a,b in zip(ticks[0,:,0],ticks[1,:,0])])
    tx,rx=vectors
    specimens={}
    with np.load(frozen/'new_readings.npz',allow_pickle=False) as readings:
        for tag in readings.files:
            observed=readings[tag]
            assert observed.shape==(1,)
            with np.load(HERE/f'{tag}_candidate.npz',allow_pickle=False) as c:
                assert np.array_equal(c['roi_indices'],roi)
                assert np.array_equal(c['observed_old_400'],old) and np.array_equal(c['observed_added_scalar_460'],observed)
                ref=c['material_public_same_class']
                assert np.array_equal(ref,references[int(ref[roi[0]])])
            name=('continued_' if tag=='unit_862' else 'dual_')+tag
            with np.load(HERE/(name+'.npz'),allow_pickle=False) as c:
                assert np.array_equal(c['roi_indices'],roi)
                assert np.array_equal(c['observed_400'],old) and np.array_equal(c['observed_460'],observed)
                other=c['material']
            for label,x in [('reference',ref),('refit',other)]:
                assert x.shape==(9216,) and np.all((x>=0)&(x<=4.5))
                assert np.all(x[roi]==x[roi[0]])
                specimens[(tag,label)]=(x,observed[0])
    result={};operators={}
    for band in (400,460):
        op_id,op=released_operator(96,band)
        wave=VerifiedWave(op)
        coeff,empty=empty_network_coefficients(op)
        operators[str(band)]=dict(operator=op_id,open_network_error_terms=coeff)
        for (tag,label),(x,observed) in specimens.items():
            tick=now()
            state=AuditedState(wave,op['S'],op['E0'],x)
            transfer=uniform_wave_error(state,state.k_total,0.,coeff)
            err=up(state.response_error+transfer['response_error_fro_upper']+
                   coeff['empty_total_Fro_error_upper']+rd(4)*(fro(state.p)+fro(empty)))
            t=state.p+empty
            ctx.prec=160
            if band==400:
                z=[ball_complex(a)-ball_complex(b) for a,b in zip(t.ravel(),old.ravel())]
                discrepancy=sum((abs(v)**2 for v in z),arb(0)).sqrt()+err
            else:
                scale=arb(8).sqrt()*arb(proposal['amplitude'])*int(m[0])/den
                scalar=sum((rx[i].conjugate()*ball_complex(t[i,j])*tx[j]
                            for i in range(8) for j in range(8)),acb(0))
                discrepancy=abs(scale*scalar-ball_complex(observed))+abs(scale)*err
            ratio=up(float((discrepancy/arb(radius)).upper()))
            branch=result.setdefault(tag,{}).setdefault(label,dict(core=int(x[roi[0]]),bands={}))
            branch['bands'][str(band)]=dict(residual_upper_over_r=ratio,current_response_error=state.response_error,
                                           model_transfer_error=transfer['response_error_fro_upper'],**elapsed(tick))
            if ratio>=1:raise ArithmeticError(f'Loaded VIE could not establish both balls at {band} MHz for {tag}/{label}: {ratio}')
    for branches in result.values():
        if {v['core'] for v in branches.values()}!={0,2}:
            raise ValueError('The two validated material tasks must differ by two')
    report=dict(model='24-port loaded native VIE, analytic empty-port wave and bounded coefficient transfer to the open probe',
                common_assumptions='Same frozen native kernel, discrete power-wave port law and Frobenius noise balls. The two verifiers share the kernel interval bounds.',
                dependencies='No search-side controller, proof, read_pair, 25-port inverse or stored predicted wave. Inputs are raw material files, frozen weights and actual acquired waves.',
                projection_precision_bits=160,core_indices_checked_in_reference_and_refit_files=True,
                radius=radius,operators=operators,samples=result,**elapsed(start))
    (HERE/'loaded_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    for tag,parts in result.items():
        print(tag,[(k,{f:w['residual_upper_over_r'] for f,w in v['bands'].items()}) for k,v in parts.items()],flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        check()
