"""Offline SRAP22 exact true-material audit, gated on the sealed online result."""

import argparse
import gc
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from scipy.linalg import lu_factor
from threadpoolctl import threadpool_limits

from online import ARRAY, FINE, HERE, SCOPE, SRAP, STUDY, ball, frobenius_squared
from online import high, load_geometry, load_public, low, read, save_json, times


class ActiveSRAP:
    """Reuse the exact retained material solve with general nonnegative contrast."""

    def __init__(self, geometry, material):
        self.geo = geometry
        self.active = np.flatnonzero(material)
        self.m = material[self.active]
        if not (np.all(np.isfinite(self.m)) and np.min(self.m)>0):
            raise ValueError('The offline SRAP material must be finite and nonnegative.')
        self.n = len(self.active)
        self.I = geometry.active_block_indices(self.active)
        self.S = geometry.S[:,self.active]
        self.B = geometry.B[self.active]
        self.A = -self.m[:,None]*geometry.centers[self.I]
        self.A.flat[::self.n+1] += 1.
        self.lu = lu_factor(self.A,check_finite=False)


def main(args):
    receipt=read(HERE/'online.json')
    if not receipt['sealed'] or len(receipt['profiles'])!=6 or len(receipt['objects'])!=22:
        raise ValueError('All six public-only profile replays for all 22 objects must be sealed first.')
    ctx.prec,ctx.threads=128,1
    start,cpu=perf_counter(),process_time()
    sys.path.insert(0,str(STUDY/'synthetic'))
    from nonvacuity import ExactFourierGeometry, StableActiveModel
    sys.path.insert(0,str(SRAP))
    import data2d

    class ExactActive(ActiveSRAP,StableActiveModel):
        pass

    source=data2d.SyntheticOperators(96,data2d.CONFIRM_ARRAYS)
    spectrum,s,b=load_geometry()
    matches=dict(spectrum=np.array_equal(spectrum,source.green.spectrum),
                 incident=np.array_equal(b,source.fine_e[:,source.tx_ids[ARRAY]]),
                 sensor=np.array_equal(s,source.fine_s[source.rx_ids[ARRAY]]))
    if not all(matches.values()):
        raise ArithmeticError('Saved 96² operator disagrees with original SRAP generator.')
    geometry=ExactFourierGeometry(FINE)
    geometry.S,geometry.B=s,b
    labels=np.random.default_rng(data2d.FAMILY_SEED).permutation(np.repeat(np.arange(6),4))
    yy={x['object_id']:(y,x) for y,x in load_public()}
    online={q['object_id']:q for q in receipt['objects']}
    rows=[]
    out=dict(scope='Independent original fine SRAP material, checked after public-only online certificate seal',
             generator='research/srap_cie_array_study_20261121/data2d.py:random_material',
             cohort='srap_main24',geometry_replay_binary64_equal=matches,precision_bits=ctx.prec,
             original_fine_side=96,algorithm='Arb inverse FFT and active-row residual with exact-port adjoint correction',
             setup=times(start,cpu),rows=rows)
    for index in SCOPE:
        oid=f's2d_{index:03d}'
        if args.object and oid not in args.object:
            continue
        t0,c0=perf_counter(),process_time()
        family=labels[index-1]
        y,meta=yy[oid]
        if meta!=online[oid]:
            raise ArithmeticError('Reopened public observation metadata differs from the sealed ledger.')
        material=data2d.random_material(source.fine_centers,index,family)
        with np.load(SRAP/'data/truth'/f"{meta['id']}.npz",allow_pickle=False) as ar:
            score_m=ar['m']
            name=ar['family'].item()
        coarse=data2d.area_average(material,96,32)
        discrepancy=float(np.max(np.abs(score_m-coarse)))
        if name!=data2d.FAMILY_NAMES[family] or discrepancy>5.e-14:
            raise ArithmeticError('Offline reconstruction does not match the archived scoring material.')
        model=ExactActive(geometry,material)
        inv,inv_note=model.inverse_bound()
        corrected,error,note=model.enclose_field(inv,args.chunk)
        distance=frobenius_squared(ball(y)-corrected).sqrt()
        truth_lower=max(0.,low(distance.lower()-error))
        truth_upper=high(distance.upper()+error)
        radius=meta['padded_radius']
        inside=arb(truth_upper)<arb(radius)
        rows.append(dict(**meta,family=name,active_cells=model.n,original_fine_peak=float(np.max(material)),
                         stored_coarse_peak=float(np.max(score_m)),
                         scoring_material_difference_max=discrepancy,
                         true_residual_lower=truth_lower,true_residual_upper=truth_upper,
                         inside_padded_ball=bool(inside),
                         admissibility_margin_lower=low(arb(radius)-arb(truth_upper)),
                         **times(t0,c0),inverse=inv_note,forward=note))
        out.update(times(start,cpu),completed_rows=len(rows))
        save_json(args.out,out)
        print(oid,'peak',float(material.max()),'cells',model.n,'distance',truth_upper,
              'radius',radius,'admitted',inside,flush=True)
        del model,corrected
        gc.collect()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--object',action='append')
    parser.add_argument('--chunk',default=128,type=int)
    parser.add_argument('--out',default=HERE/'offline.json',type=Path)
    with threadpool_limits(1):
        main(parser.parse_args())
