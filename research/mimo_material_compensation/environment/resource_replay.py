"""Measure peak process RSS of one native, independently enclosed two-material replay."""
import argparse
import json
import platform
import resource
from time import perf_counter, process_time

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits
from passive_replay import Circuit, HERE, VerifiedWave, fro, stored_response


def main(side):
    start, start_cpu = perf_counter(), process_time()
    file = "phase_32_two_shells" if side == 32 else "wave_96_2.0_1"
    with np.load(HERE / f"{file}.npz") as raw:
        ids = raw["ix"]
        x = [raw[key] if key in raw else raw[f"{key}_2.0"] for key in ("plus","minus")]
    w = Circuit(side, 400, ids)
    n = side**2
    auditor = VerifiedWave(dict(spectrum=w.symbol, grid_shape=[side,side], S=w.c, E0=w.b,
                                low_rank_left=np.zeros((n,1)),low_rank_right=np.zeros((1,n))))
    errors, defects = [], []
    with np.load(HERE / f"verified_responses_{side}_2.0_{file}_400.npz") as frozen:
        for material, key, error_key in zip(x, ("s1", "s0"), ("e1","e0")):
            sc, err, evidence = stored_response(w, auditor, material)
            observed = fro(sc - frozen[key])
            if observed > err+float(frozen[error_key]):
                raise ArithmeticError("Disjoint forward enclosures")
            errors.append(float(err))
            defects.append(float(evidence["inverse_defect"]))
    usage = resource.getrusage(resource.RUSAGE_SELF)
    report = dict(side=side, area_cells=side**2, current_support_size=len(ids), band_MHz=400,
                  measure="separate process; interpreter, input, 128-bit Green audit, two material responses, residual inverse checks included",
                  platform=platform.platform(), arb_bits=ctx.prec,
                  cpu_s=process_time()-start_cpu, wall_s=perf_counter()-start,
                  peak_RSS_MiB=usage.ru_maxrss/(2**20 if platform.system()=="Darwin" else 1024),
                  material_lu=2, material_RHS=2*(len(ids)+25), terminal_factorizations=2,
                  action="explicit current support and circuit; verified convolution kernel construction",
                  inverse_defects=defects, endpoint_response_errors=errors)
    (HERE/f"resource_replay_{side}.json").write_text(json.dumps(report,indent=2)+"\n")


if __name__ == "__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("side",type=int,choices=(32,96))
    with threadpool_limits(1):
        main(arg.parse_args().side)
