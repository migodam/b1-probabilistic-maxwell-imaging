"""After the batch is fixed, score generating fields and verify separate legal controls."""
import json

from benchmark import (HERE, STUDY, ProofGrid, SMALL_ID, clock, elapsed, endpoints,
                       freeze_config, load_model, prepare_candidates, write)
from interval_backend import ball, diagonal, norm, row_scaled
from flint import arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits


def run():
    frozen = freeze_config()
    records = []
    public = {}
    for side in (32,96):
        key = str(side)
        started = clock()
        spec = frozen['scales'][key]
        geom, ids, (_,s,b), y, r, base, sites, name = load_model(side,spec['candidate_base_upper'])
        public[f'{key}_Y'] = y
        public[f'{key}_radius'] = np.array([r])
        public[f'{key}_state_ids'] = ids
        if side==32:
            source = STUDY/'ip_guided_v4/boxes/results'/f'{SMALL_ID}__offline.npz'
            with np.load(source,allow_pickle=False) as z:
                truth, clean = z['material'],z['Y_clean']
            support = np.flatnonzero(truth)
            (g,s,b),(eg,es,eb) = ProofGrid(geom).restricted(support)
            x = [arb(float(v)) for v in truth[support]]
            j = (diagonal([1]*len(support))-row_scaled(eg,x)).solve(row_scaled(eb,x))
            proof = norm(es*j-ball(y)).upper()/arb(r)
            ratio = float(np.nextafter(float(proof.upper()),np.inf))
            passed = bool(proof < 1)
            linear = s@(truth[support,None]*b)
            analysis = 'new 128-bit pulse-scalar-TM replay on all 76 active source cells with exact exterior vacuum of this known generating material'
        else:
            source = STUDY/'ip_final_v5/strong_control/strong_ring_arrays.npz'
            with np.load(source,allow_pickle=False) as z:
                truth,clean = z['chi_active'],z['Y_clean']
            record = json.loads((source.parent/'strong_ring_verified.json').read_text())
            ratio = record['generator_check']['residual_over_radius_upper']
            passed = record['generator_check']['fits_noise_ball']
            linear = s@(truth[:,None]*b)
            support = np.flatnonzero(truth)
            analysis = 'archival 128-bit original-TM forward proof reused on native stored Fourier operator'
        born = np.linalg.norm(linear-clean)/np.linalg.norm(clean)
        weights = np.full(len(ids),spec['weight'])
        zero = np.zeros(len(ids))
        candidates = prepare_candidates(side,ids,sites,base,weights,1000,spec['seed'])
        hits = 0
        for candidate in candidates:
            lo,up = endpoints(candidate,zero,base)
            hits += bool(np.all((truth>=lo) & (truth<=up)))
        records.append(dict(observation=name,source=str(source.relative_to(STUDY)),generating_peak=float(max(truth)),
                            nonzero_material_cells=len(support),truth_in_base=bool(np.all(truth<=base)),
                            query_boxes_containing_truth=hits,
                            truth_residual_over_radius_upper=ratio,known_physical_control=passed,
                            born_relative_readout_error_float=float(born),forward_evidence=analysis,cost=elapsed(started)))
    np.savez_compressed(HERE/'public_readings.npz',**public)
    write(HERE/'physics_controls.json',records)


if __name__ == '__main__':
    ctx.prec,ctx.threads = 128,1
    with threadpool_limits(limits=1):
        run()
