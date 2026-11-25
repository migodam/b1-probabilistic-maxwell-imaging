"""Preregister a paired diagnostic using only old publicly measured signal scales."""
from time import time

from flint import arb
import numpy as np

from model import COMPONENT, OLD, ball, clocks, elapsed, hash_arrays, high, init_arithmetic, norm, read_json, write_json


def main():
    init_arithmetic()
    start=clocks()
    base=COMPONENT
    target=base/'diagnostic'
    target.mkdir(exist_ok=True)
    if (target/'registration.json').exists():
        raise RuntimeError('Diagnostic already preregistered.')
    registration=read_json(base/'registration.json')
    target_pub={}
    radii={}
    scales={}
    with np.load(base/'public_model.npz',allow_pickle=False) as pub:
        target_pub.update({k:pub[k] for k in pub.files})
        for name,scene in registration['scenes'].items():
            with np.load(base/'old'/f'{name}__added.npz',allow_pickle=False) as old:
                y=old['Y']
            # The same absolute error applies to both energy contracts. The
            # lower old amplitude in the fixed-total contract gets no noise discount.
            old_rms=norm(ball(y))/arb(len(OLD)).sqrt()
            r=high(arb("0.02")*old_rms)
            radii[name]=repr(r)
            scales[name]=dict(old_added_public_hash=hash_arrays(y),public_old_rms_norm_per_beam=high(old_rms),
                              noise_radius_decimal=radii[name])
            center=pub[name+'__center']
            for side,sgn in (('lo',-1),('hi',1)):
                target_pub[name+'__'+side]=center+sgn*.0032
            lo,hi=[target_pub[name+'__'+b] for b in ('lo','hi')]
            if not np.any(hi[0]<lo[1]):
                raise ArithmeticError('Fourfold enlargement merged whole classes.')
            scene['box_pair_hash']=hash_arrays(lo,hi)
            scene['material_halfwidth']=.0032
    np.savez_compressed(target/'public_model.npz',**target_pub)
    registration.update(phase='predeclared_paired_followup_no_diagnostic_read_yet',
         protocol='One diagnostic at nominal 2% old-ref receive noise and fourfold halfwidth. Original same six physically sealed media.',
         exposure='Paired follow-up on primary cohort. Chosen by fixed rule from old added-energy public signal; no new response, material truth or separation-optimized noise tuning.',
         noise_radius_by_scene=radii, public_scale_receipt=scales,
         recorded_unix_wall=time(), diagnostic_supplement_cost=elapsed(start))
    registration['noise_contract']['bound_per_normalized_receiver_column']='see noise_radius_by_scene; same absolute radius across both physical energy contracts'
    registration['noise_contract']['realization']='independent new diagnostic per scene and angle, shared across diagnostic policies at each physical setting'
    write_json(target/'registration.json',registration)
    print('Diagnostic noise and widened classes frozen.',radii,flush=True)


if __name__=='__main__':
    main()
