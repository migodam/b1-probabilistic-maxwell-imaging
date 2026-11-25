"""Re-run the original point-ranked direct full-wave control on all 12 CIE inputs."""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from paired import DEVELOPMENT, HERE, POOL, RX, TX, get_inputs
from pool_common import FullWave, clock, elapsed, write_arrays, write_json

spec = importlib.util.spec_from_file_location("original_guided_pool_replay", POOL/"replay.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def geometry():
    with np.load(DEVELOPMENT/"public/geometry.npz",allow_pickle=False) as a:
        return a["green"],a["sensor"],a["incident"]


def run(out):
    begin = clock()
    ctx.prec,ctx.threads=112,1
    if (out/"online.json").exists():
        raise FileExistsError(out/"online.json")
    out.mkdir(parents=True,exist_ok=True)
    g,s,b = geometry()
    table, arrays = {}, {}
    for cohort, seed, name, observations, _ in get_inputs("all"):
        with np.load(observations,allow_pickle=False) as a:
            y = a["train"] if cohort=="development" else a["Y"]
            radius = float(a["radius"])
        usage, result = prior.run_one("top_early",g,s[RX],b[:,TX],y,radius,seed+1007)
        table[name] = dict(cohort=cohort, radius=radius, **usage)
        arrays[f"{name}__material"], arrays[f"{name}__train_prediction"] = result
        print(name, "direct point-to-fit", f"{usage['wall_seconds']:.4f}s",
              f"{usage['verified_residual_upper']/radius:.4f}r", flush=True)
    write_arrays(out/"predictions.npz",**arrays)
    write_json(out/"online.json",dict(cases=table,total=elapsed(begin),
               scope="Original all-64 point ranking, within-box best ranked sample initializes direct 42-call full-wave corrector. "
                     "Stop at first independently verified feasible endpoint. Re-run on same exposed CIE objects/boxes."))


def score(out):
    table=json.loads((out/"online.json").read_text())["cases"]
    g,s,b=geometry()
    rx_held=np.ix_(np.setdiff1d(np.arange(s.shape[0]),RX),TX)
    tx_held=np.ix_(RX,np.setdiff1d(np.arange(b.shape[1]),TX))
    results,receipts={},{}
    with np.load(out/"predictions.npz",allow_pickle=False) as predictions:
        for cohort,_,name,_,private in get_inputs("all"):
            with np.load(private,allow_pickle=False) as reference:
                truth = reference["material"] if cohort=="development" else reference["truth"]
                ports = reference["noiseless_ports"] if cohort=="development" else reference["ports"]
            estimated = predictions[f"{name}__material"]
            refiner = FullWave(g,s,b,ports)
            field, _, _ = refiner.evaluate(estimated)
            receipts[f"{name}__ports"] = field
            error = float(np.linalg.norm(estimated-truth)/np.linalg.norm(truth))
            results[name] = dict(full_relative_l2=error, inside_relative_l2=error,
                 background_rms=0., heldout_receivers=float(np.linalg.norm((field-ports)[rx_held])/np.linalg.norm(ports[rx_held])),
                 heldout_sources=float(np.linalg.norm((field-ports)[tx_held])/np.linalg.norm(ports[tx_held])),
                 feasible=bool(table[name]["verified_residual_upper"]<=table[name]["radius"]))
    write_arrays(out/"offline_ports.npz",**receipts)
    write_json(out/"offline.json",results)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase",choices=("online","score"))
    parser.add_argument("--out",type=Path,default=HERE/"direct_results")
    args=parser.parse_args()
    with threadpool_limits(1):
        (run if args.phase=="online" else score)(args.out)
