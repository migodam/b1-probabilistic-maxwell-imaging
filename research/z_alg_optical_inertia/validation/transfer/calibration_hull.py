"""Bound the displacement of four deterministic calibrations of each EIGOR reading."""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter, process_time
from zipfile import ZipFile

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'manitoba'))
sys.path.insert(0, str(HERE.parents[1]/'core'))
from physical import FILE_MAP, FREQUENCIES, RAW, PREPARED, TARGETS, USED, load_file, pec_scattered
from certificate import ball, frobenius_squared, safe_float


def interval_norm(matrix):
    value = frobenius_squared(matrix).sqrt()
    return [max(0.0, safe_float(value, upper=False)), safe_float(value, upper=True)]


def main(raw, prepared):
    start_wall, start_cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = 128, 1
    with ZipFile(raw) as archive:
        data = {key: load_file(archive, file)[0] for key, file in FILE_MAP.items()}
    comparisons = []
    index = np.full(USED.shape, -1, dtype=int)
    index[USED] = np.arange(USED.sum())
    reverse = index.T[USED]
    with np.load(prepared, allow_pickle=False) as primary:
        for i, frequency in enumerate(FREQUENCIES):
            k = 2*np.pi*frequency*1.e9/299792458.
            standards = [('large', 'cal_metal_3p5in', .04445), ('small', 'held_metal_2in', .0254)]
            calibrations = {}
            for name, standard, radius in standards:
                theory = pec_scattered(k, radius)[USED]
                for blank_name, blank_file in [('before','empty_first'),('after','empty_last')]:
                    empty = data[blank_file][i][USED]
                    gain = theory/(data[standard][i][USED]-empty)
                    calibrations[name+'_'+blank_name] = (gain, empty)
            for target in TARGETS:
                reading = data[target][i][USED]
                fields = {name: gain*(reading-empty) for name,(gain,empty) in calibrations.items()}
                reference = fields['large_before']
                if not np.array_equal(reference, primary[f'f{i}_{target}'][USED]):
                    raise ValueError('Primary processing differs from the input of the certified solver')
                ref_ball = ball(reference[:,None])
                gaps = {key: interval_norm(ball(values[:,None])-ref_ball) for key,values in fields.items()}
                opposite = ball(reference[reverse,None])
                anti_norm = interval_norm((ref_ball-opposite)/2)
                # This is exactly the binary centre passed to the separate replay.
                sym_primary = ball((.5*(reference+reference[reverse]))[:,None])
                sym_gaps = {}
                for name,values in fields.items():
                    sym_value = (ball(values[:,None])+ball(values[reverse,None]))/2
                    sym_gaps[name] = interval_norm(sym_value-sym_primary)
                comparisons.append(dict(target=target, frequency_GHz=float(frequency),
                                        primary_norm=interval_norm(ref_ball),
                                        displacements=gaps, hull_radius_upper=max(x[1] for x in gaps.values()),
                                        reciprocity_error=anti_norm, symmetric_displacements=sym_gaps,
                                        symmetric_hull_radius_upper=max(x[1] for x in sym_gaps.values())))
    result=dict(scope='Convex hull of four saved-model preprocessing choices on the same 288 complex channels.',
                underlying_physical_transfer_bound='unknown', precision_bits=ctx.prec,
                timestamps='before/after are separate empty-chamber acquisitions, not independent dielectric objects',
                elapsed_wall_seconds=perf_counter()-start_wall, elapsed_cpu_seconds=process_time()-start_cpu,
                rows=comparisons)
    (HERE/'four_calibrations.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,default=RAW)
    p.add_argument('--prepared',type=Path,default=PREPARED)
    args=p.parse_args()
    with threadpool_limits(1):
        main(args.raw,args.prepared)
