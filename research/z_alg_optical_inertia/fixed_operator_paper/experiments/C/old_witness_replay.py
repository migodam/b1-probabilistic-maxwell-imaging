"""Physically verify each publicly registered box center against the old reading."""
import numpy as np
from threadpoolctl import threadpool_limits

from model import (CONTRACTS, HERE, OLD, REGISTRY, ball, clocks, elapsed, energy, exact_forward,
                   geometry, hash_arrays, high, init_arithmetic, low, noise_radius, norm,
                   radius_squared, read_json, require_p0, write_json)


def main():
    start=clocks()
    gate=require_p0()
    reg=read_json(REGISTRY)
    results=[]
    with np.load(HERE/'public_model.npz',allow_pickle=False) as pub:
        for scene,meta in reg['scenes'].items():
            geo=geometry(meta['layout'])
            lo,hi,centers=[pub[scene+'__'+q] for q in ('lo','hi','center')]
            for contract in CONTRACTS:
                old_power,_=energy(contract,1)
                with np.load(HERE/'old'/f'{scene}__{contract}.npz',allow_pickle=False) as raw:
                    y=raw['Y']
                radius=radius_squared(len(OLD),noise_radius(reg,scene)).sqrt()
                for label,l,u,c in zip(meta['classes'],lo,hi,centers,strict=True):
                    if not np.all((l<=c)&(c<=u)):
                        raise ArithmeticError('Preregistered nominal is outside own box.')
                    f=exact_forward(geo,c,OLD,old_power.sqrt())
                    rn=norm(ball(y)-f)
                    results.append(dict(scene=scene,layout=meta['layout'],contract=contract,candidate=label,
                        class_center_hash=hash_arrays(c),old_data_hash=hash_arrays(y),
                        status='WITNESS_FOUND' if rn<radius else 'UNRESOLVED',
                        forward_residual_upper=high(rn),old_radius_upper=high(radius),
                        radius_margin_lower=low(radius-rn)))
    summary=dict(status='Post-freeze physical proof of old ambiguity, not used to select actions',
        gate=gate,rows=results,verified_witnesses=sum(r['status']=='WITNESS_FOUND' for r in results),
        **elapsed(start))
    write_json(HERE/'old_box_witnesses.json',summary)
    print(summary['verified_witnesses'],'of',len(results),'public box centers fit the old ball.',flush=True)


if __name__=='__main__':
    init_arithmetic()
    with threadpool_limits(1):
        main()
