"""Read the already acquired two balls and independently enclose each proposed physical material."""
from fractions import Fraction
import json

import numpy as np
from threadpoolctl import threadpool_limits

from run import HERE, SIDE, controller, elapsed, proof, read_pair, timer


def upper(x):
    value=float(x)
    return float(np.nextafter(value,np.inf)) if Fraction(value)<x else value


def verify():
    clock=timer()
    proposal,labels,observations,tx,rx,amplitude,_=controller()
    _,_,roi_mask,r,_=read_pair()
    roi=np.flatnonzero(roi_mask)
    assert len(roi)==164 and proposal['detector_radius']==r
    samples={}
    for name,data in observations.items():
        candidates=sorted(HERE.glob(f'*_{name}.npz'))
        candidates=[p for p in candidates if p.stem.startswith(('dual_','continued_'))]
        if not candidates:
            raise FileNotFoundError(name)
        with np.load(HERE/f'{name}_candidate.npz',allow_pickle=False) as p:
            materials={'public_reference':p['material_public_same_class']}
        for file in candidates:
            with np.load(file,allow_pickle=False) as p:
                assert np.array_equal(roi,p['roi_indices'])
                assert np.array_equal(p['observed_400'],data[0])
                assert np.array_equal(p['observed_460'],data[1])
                materials[file.stem]=p['material']
        tested={}
        for key,x in materials.items():
            assert x.size==SIDE*SIDE and np.all((x>=0)&(x<=4.5))
            core=sum((Fraction(p) for p in x[roi]),Fraction(0))/len(roi)
            assert core in (0,2) and np.all(x[roi]==int(core))
            res,_,cost=proof(x,data,r,tx,rx,amplitude,proposal['amplitude'])
            bounds={b:v['residual_upper_over_r'] for b,v in res.items()}
            # Treat the 8-column old transfer error conservatively as spectral.
            # Its Frobenius norm is at most sqrt(8) times as large; 3 > sqrt(8).
            b=res['400']
            bounds['400']=upper(Fraction(bounds['400'])+2*Fraction(b['response_error_fro_upper'])/Fraction(r))
            tested[key]=dict(core_mean=int(core),residual_upper_over_r=bounds,
                             both_balls=max(bounds.values())<=1,material_nonzero=int(np.count_nonzero(x)),
                             native_audit=res,replay_cost=cost)
        feasible=[a['core_mean'] for a in tested.values() if a['both_balls']]
        samples[name]=dict(sealed_label=labels[name]['label'],materials=tested,
                           diameter_lower=max(feasible)-min(feasible) if feasible else None)
    report=dict(grid=SIDE,area_m2=len(roi)*(2/SIDE)**2,observations='one common 8x8 at 400 MHz and one recorded complex scalar per specimen at 460 MHz',
                material_box=[0,4.5],absolute_radius=r,added_energy=proposal['available_energy_upper'],
                old_fro_conversion='three times the original response error instead of one; 3 exceeds sqrt(8)',
                samples=samples,**elapsed(clock))
    (HERE/'dual_verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for key,s in samples.items():
        print(key,s['diameter_lower'],[(k,v['residual_upper_over_r']) for k,v in s['materials'].items()],flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        verify()
