"""A real batch of 24 observations per physical array. Charge preparation to all 24."""
from concurrent.futures import ProcessPoolExecutor
import json

from run import Benchmark, HERE, collection, clocks, spent, threadpool_limits
from screen import Checks

def batch(aperture):
    with threadpool_limits(1):
        cfg=json.loads((HERE/"FREEZE.json").read_text())
        names=[r["case"] for r in collection.planned_cases() if r["array"]==aperture]
        t=clocks()
        _,p=collection.get_native_case(names[0])
        cert=Checks(p["spectrum"],p["S"],p["E"])
        cert.uniform.prepare()
        cert.tile_projection.prepare()
        cost=spent(t)
        out=HERE/"predictions"
        collection.write_json(HERE/f"preparation_{aperture}.json",dict(
            aperture=aperture,n_observations=len(names),n_independent_materials=12,
            uniform=cert.uniform.prepare_cost,tile=cert.tile_projection.setup_cost,**cost))
        for name in names:
            if (out/f"{name}__pruned_reuse.json").exists():continue
            e=Benchmark(name,"pruned_reuse",cfg,shared=cert,allocated_wall=cost["wall"]/len(names),
                        allocated_cpu=cost["cpu"]/len(names))
            print(*e.run(out),flush=True)
    return aperture

if __name__=="__main__":
    with ProcessPoolExecutor(2) as ex:
        for key in ex.map(batch,collection.APERTURES):
            print("complete",key,flush=True)
