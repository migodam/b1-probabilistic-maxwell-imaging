"""Isolated simulated instrument process. Future reads require the sealed proposal."""
import argparse
from time import time

from flint import arb
import numpy as np
from threadpoolctl import threadpool_limits

from model import (COMPONENT, CONTRACTS, HERE, OLD, REGISTRY, ball, clocks,
                   elapsed, energy, exact_forward, geometry, hash_arrays, hash_file, high,
                   init_arithmetic, mid, noise_radius, norm, radius_squared, read_json,
                   require_p0, write_json)

PRIVATE = COMPONENT/'.sealed'


def material_truth(scene,public):
    f=PRIVATE/f'{scene}.npz'
    if f.exists():
        with np.load(f,allow_pickle=False) as data:
            return data['material'],int(data['owner'])
    rng=np.random.default_rng(np.random.SeedSequence([862410, *scene.encode()]))
    owner=int(rng.integers(2))
    lo,hi=[public[scene+'__'+key][owner] for key in ('lo','hi')]
    q=rng.uniform(.26,.74,size=len(lo))
    m=lo+(hi-lo)*q
    f.parent.mkdir(mode=0o700,exist_ok=True)
    np.savez_compressed(f,material=m,owner=np.array(owner))
    f.chmod(0o600)
    return m,owner


def measure(scene,geo,x,angles,amplitude,nr):
    start=clocks()
    clean=exact_forward(geo,x,angles,amplitude)
    values=mid(clean)
    gaussian=[]
    for a in angles:
        random=np.random.default_rng(np.random.SeedSequence([831774 if HERE!=COMPONENT else 831773,a,*scene.encode()]))
        noise=random.standard_normal(geo.s.shape[0])+1j*random.standard_normal(geo.s.shape[0])
        gaussian.append(noise/np.linalg.norm(noise))
    error=float(arb(nr))*0.65*np.column_stack(gaussian)
    got=values+error
    r=radius_squared(len(angles),nr).sqrt()
    verified=norm(ball(got)-clean)
    if not verified<r:
        raise ArithmeticError('Simulated acquisition violated the fixed port noise contract.')
    return got,dict(exact_forward_distance_upper=high(verified),declared_radius_upper=high(r),
                    operator_scope='OUTWARD_ANALYTIC_TM_DISCRETIZATION', **elapsed(start))


def old(gate,reg):
    start=clocks()
    dst=HERE/'old'
    dst.mkdir(exist_ok=True)
    observations={}
    with np.load(HERE/'public_model.npz',allow_pickle=False) as pub:
        for scene,meta in reg['scenes'].items():
            x,_=material_truth(scene,pub)
            geo=geometry(meta['layout'])
            for contract in CONTRACTS:
                old_power,_=energy(contract,1)
                y,proof=measure(scene,geo,x,OLD,old_power.sqrt(),noise_radius(reg,scene))
                name=f'{scene}__{contract}'
                path=dst/(name+'.npz')
                if not path.exists():
                    np.savez_compressed(path,Y=y)
                else:
                    with np.load(path,allow_pickle=False) as f:
                        if not np.array_equal(f['Y'],y):
                            raise ArithmeticError('Existing old data cannot be changed.')
                observations[name]=dict(scene=scene,energy_contract=contract,
                       old_data_hash=hash_arrays(y),old_beams=OLD,old_energy=float(old_power*len(OLD)),**proof)
    write_json(dst/'old_acquisition.json',dict(scope='Old reads only. No truth or future arrays in old design files.',
                gate=gate,observations=observations,**elapsed(start)))
    print(f'Released {len(observations)} old observations only.',flush=True)


def future(gate,reg):
    start=clocks()
    sealed=HERE/'proposal.json'
    if not sealed.is_file():
        raise RuntimeError('No future read before proposal freeze.')
    record=read_json(sealed)
    if record['status']!='SEALED_READ_ONLY_BEFORE_ANY_NEW_DATA':
        raise RuntimeError('Proposal seal state invalid.')
    dst=HERE/'acquired'
    dst.mkdir(exist_ok=True)
    old_cert=read_json(HERE/'old/old_acquisition.json')['observations']
    rows={}
    files={}
    with np.load(HERE/'public_model.npz',allow_pickle=False) as pub:
        for row in record['source_rows']:
            scene,contract=row['scene'],row['energy_contract']
            name=f'{scene}__{contract}'
            if old_cert[name]['old_data_hash']!=row['old_data_hash']:
                raise ArithmeticError('Frozen old array hash discrepancy.')
            x,expected=material_truth(scene,pub)
            if not np.all((pub[scene+'__lo'][expected]<=x)&(x<=pub[scene+'__hi'][expected])):
                raise ArithmeticError('Witness is outside its exact predeclared box.')
            geo=geometry(reg['scenes'][scene]['layout'])
            for size,proposal in row['batches'].items():
                new_s=int(size)
                _,power=energy(contract,new_s)
                for policy,action in proposal['actions'].items():
                    angles=tuple(action['source_indices'])
                    key=name+'__m'+size+'__'+'-'.join(map(str,angles))
                    action_id=name+'__m'+size+'__'+policy
                    if key not in files:
                        y,proof=measure(scene,geo,x,angles,power.sqrt(),noise_radius(reg,scene))
                        np.savez_compressed(dst/(key+'.npz'),Y=y)
                        files[key]=dict(filename=key+'.npz',data_hash=hash_arrays(y),scene=scene,
                              angles_indices=angles,**proof)
                    rows[action_id]=dict(acquisition_key=key,witness_class=reg['scenes'][scene]['classes'][expected],
                       material_hash=hash_arrays(x),contract=contract,scene=scene,
                       new_total_energy=float(power*new_s),old_data_hash=row['old_data_hash'])
    write_json(dst/'acquisition.json',dict(gate=gate,frozen_proposal_stat_mtime=sealed.stat().st_mtime, frozen_proposal_hash=hash_file(sealed),
       proposal_timestamp=record['wall_time_unix'],acquisition_end_unix=time(),
       reading_files=files,policy_rows=rows,unique_physical_evaluations=len(files),
       counted_policy_paths=len(rows),**elapsed(start)))
    print(f'{len(files)} unique prospective new observations, {len(rows)} policy-priced paths.',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('old','future'))
    args=p.parse_args()
    gate=require_p0()
    registration=read_json(REGISTRY)
    if args.stage=='old':
        old(gate,registration)
    else:
        future(gate,registration)


if __name__=='__main__':
    init_arithmetic()
    with threadpool_limits(1):
        main()
