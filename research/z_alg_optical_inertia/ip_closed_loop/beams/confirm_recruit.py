"""Independent material confirmation cohort with old-channel-only recruitment."""
import argparse

import numpy as np
from threadpoolctl import threadpool_limits

from beam_common import (HERE, OLD, atomic_write, clocks, elapsed, energy_gram,
                         forward, geometry, pair_jacobian, profiles, require_absent,
                         save, shapes)

CONFIRM_SEED=202611260830
FIXED_FLUENCE=.014
CHANGE=.9
OLD_REL_CUTOFF=.76
WIDTH_BG=.012
WIDTH_FG=.0025
GROUP={32:tuple(range(300,312)),96:tuple(range(400,402))}
KINDS=('solid','ring','double','austria')


def recruit(geo, idx):
    start=clocks()
    kind=KINDS[idx%4]
    cells,fg,bg,base,xy=shapes(geo.side,idx,kind)
    g,s,b=geo.restricted(cells)
    # There is exactly one old column here. Neither future dictionary columns
    # nor future simulated readings enter the profile, screen or scale.
    old=b[:,OLD,None]
    y,e,d,*_=forward(g,s,old,base,auxiliary=True)
    funcs=profiles(xy,fg)
    orth,_,_=np.linalg.svd(funcs,full_matrices=False)
    j=pair_jacobian(e,d,orth)[:,0,:]
    real=np.vstack((j.real,j.imag))
    _,_,vh=np.linalg.svd(real,full_matrices=True)
    radius=.04*np.linalg.norm(y)
    rng=np.random.default_rng(np.random.SeedSequence([CONFIRM_SEED,idx]))
    rec=dict(index=idx,kind=kind,old_radius=radius,candidate_new_noise_radius=.65*radius,
             fixed_fluence=FIXED_FLUENCE,selected=False,trials=[])
    arrays=None
    for candidate in range(6):
        coefficients=rng.normal(size=5)
        contrast=orth@(vh[-5:,:].T@coefficients)
        contrast*=CHANGE/max(abs(contrast))
        other=base+contrast
        predicted=forward(g,s,old,other)
        distance=np.linalg.norm(predicted-y)/radius
        legal=bool(np.min(other)>=0)
        rec['trials'].append(dict(candidate=candidate,old_center_distance_over_radius=distance,legal=legal))
        if arrays is None and legal and distance<=OLD_REL_CUTOFF:
            rec.update(selected=True,candidate=candidate,old_relative=distance,
                       material_difference_rmse=float(np.linalg.norm(contrast)/np.sqrt(len(contrast))))
            widths=np.where(bg,WIDTH_BG,WIDTH_FG)
            arrays=dict(cells=cells,fg=fg,bg=bg,centers=np.stack((base,other)),widths=widths)
    return dict(rec,**elapsed(start)),arrays


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=sorted(GROUP))
    args=parser.parse_args()
    folder=HERE/f'confirm_{args.side}'
    folder.mkdir(exist_ok=True)
    path=folder/'registration.json'
    require_absent(path,folder/'public_classes.npz')
    start=clocks()
    geo=geometry(args.side)
    arrays=dict(power=energy_gram(geo),rx_indices=np.arange(16),tx_indices=np.arange(16))
    records=[]
    for idx in GROUP[args.side]:
        rec,arr=recruit(geo,idx)
        records.append(rec)
        if arr is not None:
            arrays.update({f'{idx}__{k}':v for k,v in arr.items()})
        print(idx,rec['kind'],rec['selected'],rec.get('old_relative'),flush=True)
    save(folder/'public_classes.npz',**arrays)
    atomic_write(path,dict(side=args.side,stage='REGISTERED_BEFORE_DATA',
                           origin='New disjoint seeds and heterogeneous deterministic material families. Candidate perturbations use old Jacobian only, independently of every future direction and reading. First old-only passing item out of six is kept, six recorded regardless of screen. Class labels in each recruited pair remain sealed until readout.',
                           material_change_maximum=CHANGE,old_only_cutoff=OLD_REL_CUTOFF,
                           foreground_contrast_uncertainty_half_width=WIDTH_FG,
                           background_contrast_half_width=WIDTH_BG,
                           background_support='Two weak unknown lobes with independent per-cell values; exact vacuum outside allowed supports.',
                           receivers='half ring (16 out of 32, 0 through 168.75 deg)',
                           illuminators='half dictionary (16 out of 32, 0 through 168.75 deg)',
                           frequency_hz=400000000.,side_m=2.,old_indices=[OLD],
                           energy='One fixed new shot with 0.014 full-square fluence, in beam-0 energy units; chosen once after engineering pilot, no per-object design calibration.',
                           noise='Old radius 0.04 times class-0 nominal old response norm; new radius 0.65 old radius. Fixed before prospective action design and data acquisition.',
                           cases=records,**elapsed(start)))


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
