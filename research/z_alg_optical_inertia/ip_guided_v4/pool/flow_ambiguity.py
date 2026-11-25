"""Offline, certify separated A13 materials explaining the same archived reading."""

from itertools import combinations
import json

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from box_certificates import ball, lower, norm
from flow_audit import OUTPUT
from flow_verify import SignedNative
from pool_common import FullWave, clock, elapsed, write_json


def difference(a, b, scale):
    aa, bb = ball(a[:, None]), ball(b[:, None])
    return lower((norm(aa-bb).lower()/scale.upper()).lower())


def main():
    ctx.prec, ctx.threads = 96, 1
    begin = clock()
    sources = json.loads((OUTPUT/"online.json").read_text())["cases"]
    checks = json.loads((OUTPUT/"strict_point_witnesses.json").read_text())["cases"]
    with np.load(OUTPUT/"native_geometry.npz", allow_pickle=False) as geom:
        g, s, b = geom["green"], geom["sensor"][geom["used_rx"]], geom["incident"]
    verifier = SignedNative(g,s,b)
    results = {}
    with np.load(OUTPUT/"offline_geometry_scores.npz", allow_pickle=False) as offline, \
            np.load(OUTPUT/"continuation_arrays.npz", allow_pickle=False) as discovered:
        for case, r in sources.items():
            truth = offline[f"{case}__truth"]
            scale = norm(ball(truth[:,None]))
            radius = r["radius"]
            with np.load(OUTPUT/f"{case}.npz") as evidence:
                observed, lo, hi = evidence["train"], evidence["lows"], evidence["highs"]
            model = FullWave(g,s,b,observed)
            _, j, _ = model.evaluate(truth)
            archived = verifier.point(truth,j,observed)
            if archived["strict_residual_upper"] > radius:
                raise ArithmeticError("Archived material is not within the declared observation ball.")
            physical = checks[case]["boxes"]
            feasible = sorted(int(i) for i,v in physical.items() if v["strictly_feasible"])
            exact_original = [i for i,h in enumerate(hi) if np.all((lo[i]<=truth)&(truth<=h))]
            observed_materials = {i: discovered[f"{case}__box{i:02}__material"] for i in feasible}
            paired = []
            for i,a in observed_materials.items():
                paired.append(dict(verified_box=i, truth_material_separation_relative_lower=difference(a,truth,scale),
                                   residual_radius_ratio=physical[str(i)]["strict_residual_upper"]/radius))
            alternatives = []
            for i,j in combinations(feasible, 2):
                disjoint = np.flatnonzero((hi[i]<lo[j])|(hi[j]<lo[i]))
                alternatives.append(dict(pair=[i,j], disjoint_boxes=bool(disjoint.size),
                    disjoint_cells=int(disjoint.size),
                    material_separation_relative_lower=difference(observed_materials[i],observed_materials[j],scale)))
            best = min(paired, key=lambda v:v["residual_radius_ratio"])
            results[case] = dict(radius=radius, verified_original=archived, original_inside_boxes=exact_original,
                                 best_box=best, alternatives_to_truth=paired, pairs_between_boxes=alternatives)
            print(case,"truth/r",archived["strict_residual_upper"]/radius, "best relative separation",
                  best["truth_material_separation_relative_lower"], "disjoint feasible pairs",sum(v["disjoint_boxes"] for v in alternatives),
                  flush=True)
    write_json(OUTPUT/"ambiguity_offline.json", dict(cases=results,setup=verifier.setup,total=elapsed(begin),
        scope="Strict finite observation-ball ambiguity with distinct, explicitly verified physical materials. "
             "The observed training channels, grid, operator, and radius are shared; truth is loaded only offline."))


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
