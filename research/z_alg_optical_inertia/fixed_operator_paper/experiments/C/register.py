"""Freeze public competing boxes and the acquisition protocol before any truth."""
from time import time

import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits

from model import (ANGLES, BATCHES, COLUMN_RADIUS, CONTRACTS, HERE, OLD, POLICIES,
                   PORTS, REGISTRY, SEEDS, SIDE, WEIGHTS, WIDTH_WAVELENGTHS,
                   WAVELENGTH, clocks, elapsed, forward_jacobian, geometry,
                   hash_arrays, init_arithmetic, numeric_forward, write_json)


def template(index, rng):
    xy=(np.indices((SIDE,SIDE))+.5)/SIDE-.5
    centers=((-0.15,.13),(.16,-.18),(-.17,-.15))
    z=xy-np.asarray(centers[index])[:,None,None]
    g=.14+.31*rng.uniform(size=(SIDE,SIDE))+1.26*np.exp(-np.sum(z*z,axis=0)/.102)
    return g.ravel()


def public_pair(geo, reference):
    start=clocks()
    y0,jac0=forward_jacobian(geo,reference,OLD)
    _,singular,vh=np.linalg.svd(jac0,full_matrices=True)
    rank=int(np.sum(singular>singular.max()*1.e-10))
    null=vh[rank:].T
    options=np.sum(null*null,axis=1)
    anchor=int(np.argmax(options))
    if options[anchor] < .1:
        raise ArithmeticError('Public ambiguity design failed to find a visible anchor.') 
    perturb=null@(null[anchor]/options[anchor])
    delta=.23
    new=reference+delta*perturb
    free=np.arange(len(reference))!=anchor
    lo=np.maximum(.025,reference-.65)
    hi=np.minimum(2.95,reference+.65)
    target_norm=np.linalg.norm(y0)
    calls=0

    def evaluate(x):
        nonlocal calls
        calls+=1
        value=new.copy()
        value[free]=x
        f,j=forward_jacobian(geo,value,OLD)
        error=(f-y0).ravel()
        return np.r_[error.real,error.imag]/target_norm,j[:,free]/target_norm

    result=least_squares(lambda x:evaluate(x)[0],np.clip(new[free],lo[free]+1.e-7,hi[free]-1.e-7),
                         jac=lambda x:evaluate(x)[1],bounds=(lo[free],hi[free]),max_nfev=150,
                         gtol=2.e-14,ftol=2.e-14,xtol=2.e-14)
    new[free]=result.x
    discrepancy=float(np.linalg.norm(numeric_forward(geo,new,OLD)-y0))
    return np.stack((reference,new)),dict(anchor_cell=anchor,public_anchor_change=delta,
       public_old_sensitivity_rank=rank,public_old_nominal_residual=discrepancy,
       optimizer_calls=calls,stopped=str(result.message),**elapsed(start))


def main():
    start=clocks()
    if REGISTRY.exists():
        print('Public preregistration already frozen.')
        return
    public={}
    rows={}
    for layout,seed in SEEDS.items():
        geo=geometry(layout)
        for index in range(2):
            scene=f'{layout}_{index+1:02d}'
            rng=np.random.default_rng(np.random.SeedSequence([seed,index,9]))
            cen,boundary=public_pair(geo,template(index+(1 if layout=='occluded' else 0),rng))
            width=.0008
            lo,hi=cen-width,cen+width
            if not np.any(hi[0]<lo[1]):
                raise ArithmeticError('Public classes overlap.')
            rows[scene]=dict(layout=layout, classes=['box_a','box_b'],nominal_fibre_design=boundary,
                material_hash_scheme='BLAKE2b of shape, dtype and contiguous arrays',
                coverage='exact two-box union task prior, independent per-cell intervals',
                model_hash=hash_arrays(geo.g,geo.s,geo.b),
                box_pair_hash=hash_arrays(lo,hi))
            public[scene+'__lo']=lo
            public[scene+'__hi']=hi
            public[scene+'__center']=cen
        public[layout+'__G']=geo.g
        public[layout+'__S']=geo.s
        public[layout+'__Bdict']=geo.b
    np.savez_compressed(HERE/'public_model.npz',**public)
    write_json(REGISTRY,dict(phase='frozen_public_design_no_true_material_no_observations',
       experiment='prospective active TM design',model='INTERVAL_DISCRETIZATION, outgoing 16 cell pulse disk',
       frequency_hz=400e6,wavelength_m=WAVELENGTH,width_wavelengths=WIDTH_WAVELENGTHS,side=SIDE,
       time_convention='exp(-i omega t)', port_coordinates='angle degrees and radius rational in wavelengths',
       arrays=PORTS,old_source_indices=OLD,candidate_angles_degrees=ANGLES,allowed_batch_sizes=BATCHES,
       same_paired_candidates=True,policies=POLICIES,energy_contracts=CONTRACTS,
       energy_units='Integral over one unit time of area-average incident |B|^2 on the test square per plane wave.',
       energy_rule=dict(added='2 old beams energy 1 each plus 1 per future beam',
                        fixed_total='independently measured two old beams energy 1/2 each, new energy totals 1'),
       noise_contract=dict(bound_per_normalized_receiver_column=COLUMN_RADIUS,
         realization='independently paired by physical beam, fixed absolute readout error budget; observed perturbation is 0.65 times per-column budget',
         joint_radius='sqrt(sum disjoint column noise squares); Omega_old = r_old^2 I_source'),
       weights=WEIGHTS, broad_indefinite_diagnostic_box=[0,3], scenes=rows,
       no_adaptive_budget_reuse='One independently planned follow-up batch per policy and size; stop before acquiring if only one class survived the old data. This is not a cumulative sequence of three batches.',
       recorded_unix_wall=time(),**elapsed(start)))
    print('Froze public boxes; no material truth or readings created.')
    for key,row in rows.items():
        print(key,row['nominal_fibre_design']['public_old_nominal_residual'])


if __name__=='__main__':
    init_arithmetic()
    with threadpool_limits(1):
        main()
