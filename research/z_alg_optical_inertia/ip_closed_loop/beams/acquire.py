"""Instrument simulator: record genuine original-operator fields after freeze."""
import argparse
from datetime import datetime, timezone

import numpy as np
from threadpoolctl import threadpool_limits

from beam_common import (HERE, OLD, NOISE_DRAW, atomic_write, clocks, digest, elapsed,
                         forward, geometry, require_absent, save)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def rng_for(index,purpose):
    return np.random.default_rng(np.random.SeedSequence([NOISE_DRAW,index,purpose]))


def noise(index,purpose,radius,ports):
    rng=rng_for(index,purpose)
    val=rng.normal(size=ports)+1j*rng.normal(size=ports)
    return radius*val/np.linalg.norm(val)


def material(data, index):
    actual=rng_for(index,8)
    center=data[f'{index}__centers']
    widths=data[f'{index}__widths']
    offset=actual.uniform(-.65,.65,len(widths))
    label=int(actual.integers(2))
    return label,center[label] + widths*offset


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['old','future'])
    parser.add_argument('side',type=int,choices=[32,96])
    parser.add_argument('--tag',default='main')
    parser.add_argument('--folder',type=type(HERE),help='Use an independently staged instrument bundle.')
    args=parser.parse_args()
    folder=args.folder or HERE/f'{args.tag}_{args.side}'
    old=args.stage=='old'
    stage_path=folder/(args.stage+'_acquisition.json')
    out=folder/(args.stage+'_readings.npz')
    audit=folder/(args.stage+'_offline.npz')
    require_absent(stage_path,out,audit)
    if old:
        require_absent(folder/'proposal.json')
    if not old and not (folder/'proposal.json').is_file():
        raise FileNotFoundError('A frozen proposal is required before new acquisitions.')
    import json
    reg=json.loads((folder/'registration.json').read_text())
    bank=np.load(folder/'public_classes.npz')
    t=clocks()
    geom=geometry(args.side)
    records=[]
    readings={}
    offline={}
    if not old:
        prop=json.loads((folder/'proposal.json').read_text())
        bank_design=np.load(folder/'proposal_arrays.npz')
        if prop['old_file_digest']!=digest(np.frombuffer((folder/'old_readings.npz').read_bytes(),dtype=np.uint8)):
            raise ValueError('Old data changed after freeze.')
        if prop['design_binary_digest']!=digest(np.frombuffer((folder/'proposal_arrays.npz').read_bytes(),dtype=np.uint8)):
            raise ValueError('Design arrays changed after freeze.')
    for record in reg['cases']:
        if not record['selected']:
            continue
        tick=clocks()
        idx=record['index']
        i=f'{idx}__'
        cells=bank[i+'cells']
        label,chi=material(bank,idx)
        g,s,b=geom.restricted(cells)
        if old:
            v=np.eye(b.shape[1],dtype=complex)[:,[OLD]]
            names=['initial']
            e=noise(idx,13,.28*record['old_radius'],s.shape[0])
        else:
            if str(idx) not in prop['rows']:
                continue
            names=prop['rows'][str(idx)]['policies']
            v=bank_design[i+'actions']
            e=noise(idx,15,.62*record['candidate_new_noise_radius'],s.shape[0])
        clean, _, _, current, born, _=forward(g,s,b@v,chi,auxiliary=True)
        y=clean+e[:,None]
        readings[i+'Y']=y
        if not old:
            readings[i+'V']=v
        offline[i+'label']=np.array(label)
        offline[i+'chi']=chi
        offline[i+'clean']=clean
        offline[i+'currents']=current
        offline[i+'error']=e
        offline[i+'born']=born
        ident=f'ipcl_beam_{idx:03d}__s{args.side}__tx180_rx180_16x16__{"old1" if old else "added1"}'
        records.append(dict(observation_id=ident,index=idx,policies=names,
                            full_wave_state_dimension=len(chi), physical_rhs=v.shape[1],
                            float_y_digest=digest(y),**elapsed(tick)))
    save(out,**readings)
    save(audit,**offline)
    atomic_write(stage_path,dict(stage=args.stage,completed_utc=stamp(),records=records,
                                data_visibility='offline file withheld from design',
                                source='independent full-wave solve at material within sealed boxes',
                                paired_noise='same 16 complex receiver errors on every policy, fixed absolute radius',
                                **elapsed(t)))
    print(args.stage,args.side,len(records),elapsed(t),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
