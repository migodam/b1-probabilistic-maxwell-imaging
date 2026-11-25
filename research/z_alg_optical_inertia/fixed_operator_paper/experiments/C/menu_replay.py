"""Recreate all archived pre-acquisition menu checks from public inputs only."""
import numpy as np
from threadpoolctl import threadpool_limits

from design import design_scene
from model import (HERE, REGISTRY, clocks, elapsed, hash_arrays, hash_file, init_arithmetic,
                   noise_radius, read_json, require_p0, write_json)


def main():
    start=clocks()
    gate=require_p0()
    proposal=read_json(HERE/'proposal.json')
    registration=read_json(REGISTRY)
    results=[]
    query_arrays={}
    with np.load(HERE/'public_model.npz',allow_pickle=False) as pub:
        for frozen in proposal['source_rows']:
            scene,contract=frozen['scene'],frozen['energy_contract']
            with np.load(HERE/'old'/f'{scene}__{contract}.npz',allow_pickle=False) as reader:
                y=reader['Y']
            if hash_arrays(y)!=frozen['old_data_hash']:
                raise ArithmeticError('Archived old input changed.')
            audits={}
            replay=design_scene(scene,registration['scenes'][scene],pub,y,contract,query_arrays,
                                noise_radius(registration,scene),audit=audits)
            for size,batch in frozen['batches'].items():
                copy=replay['batches'][size]
                for key in ('all_action_menu','attempted_interval_actions','strict_separating_actions'):
                    if batch[key]!=copy[key]:
                        raise ArithmeticError(f'Archived menu mismatch: {scene}, {contract}, {size}, {key}')
                for policy,action in batch['actions'].items():
                    repeated=copy['actions'][policy]
                    if list(action['source_indices'])!=list(repeated['source_indices']):
                        raise ArithmeticError('Archived policy action changed.')
                    a,b=action['verified_support'],repeated['verified_support']
                    if (a is None)!=(b is None) or a and a['support_query_hash']!=b['support_query_hash']:
                        raise ArithmeticError('Archived query differs.')
            results.append(dict(scene=scene,contract=contract,by_batch=audits))
    np.savez_compressed(HERE/'menu_query_replay.npz',**query_arrays)
    write_json(HERE/'full_menu_replay.json',dict(
        status='RECREATED_FROM_OLD_DATA_AND_CHECKED_AGAINST_PRE_ACQUISITION_ARCHIVE',
        observed_future_data_read=False, interval_bits=128, gate=gate,
        proposal_hash=hash_file(HERE/'proposal.json'),query_arrays='menu_query_replay.npz',
        rows=results,**elapsed(start)))
    print(f'Reproduced every menu summary and selected action on {len(results)} old conditions.')


if __name__=='__main__':
    init_arithmetic()
    with threadpool_limits(1):
        main()
