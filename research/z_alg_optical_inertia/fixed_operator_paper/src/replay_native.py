"""Replay the two fixed 96 x 96 measurements, using the published proof primitives."""
import json
import sys
import time
from pathlib import Path
import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT.parent
BOXES = STUDY / 'ip_guided_v4/boxes'
RING = STUDY / 'ip_final_v5/strong_control'
sys.path[:0] = [str(BOXES), str(RING)]
from run import CONTROL_ID, CONTROL, load, Observation
from verify import ProofGrid, prove
from run_strong_ring import verify_generator_in_ball

def main():
    ctx.prec = 128
    ctx.threads = 4
    res = BOXES/'results'
    disk = load(res, [CONTROL_ID])[CONTROL_ID]
    ring = np.load(RING/'strong_ring_arrays.npz')
    meta = json.loads((RING/'strong_ring_verified.json').read_text())
    saved = json.loads((res/'followup_disk.json').read_text())['cases'][CONTROL_ID]['strict_checks']['rank_tail_disk']
    inputs = np.load(res/'followup_disk_witnesses.npz')
    key = CONTROL_ID+'__rank_tail_disk__'
    ids, lower, upper = [inputs[key+name] for name in ('active','lower','upper')]
    material = np.load(CONTROL/(CONTROL_ID+'.npz'))['material'].ravel()
    assert np.all(material[np.setdiff1d(np.arange(material.size),ids)]==0)
    solid = ring['c_solid']
    w = meta['chosen_solid_width']
    obs2 = Observation(meta['observation_id'], disk.geometry,ring['Y'],float(ring['radius'][0]),'stored')
    cases = [(disk, ids, lower, upper, saved['mu'],material[ids]),
             (obs2,ring['ids'],np.maximum(0,solid-w),solid+w,meta['solid_peak4_arb_certificate']['mu'],ring['chi_active'])]
    records = []
    for obs, region, l, u, mu, physical in cases:
        start=time.perf_counter()
        grid=ProofGrid(obs.geometry)
        cert,_=prove(obs,grid,region,l,u,mu,try_spectral=True,try_dark=False)
        witness=verify_generator_in_ball(grid,obs,region,physical)
        assert cert['strictly_excludes'] and witness['fits_noise_ball']
        records.append(dict(id=obs.id,precision=128,allowable_cells=len(region),actual_nonzero_cells=int(np.count_nonzero(physical)),peak=float(max(physical)),cap=cert,witness=witness,seconds=time.perf_counter()-start))
        (ROOT/'results/native_replay.json').write_text(json.dumps(records,indent=2))
        print(obs.id,cert['distance_over_radius_lower'],witness['residual_over_radius_upper'],flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=3):main()
