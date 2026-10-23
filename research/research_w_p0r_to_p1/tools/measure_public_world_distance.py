"""Decode the frozen public joint endpoints on one common cell-centre grid.

No Book, FEM, GT, optimization or new response. These are sampled material
distances, not continuous geometry certificates or posterior probabilities.
"""
import json
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    cpu, wall = time.process_time(), time.perf_counter()
    sys.path[:0] = [str(ROOT/'src'), str(ROOT.parent/'research_w_continuation/src'),
                   str(ROOT.parent/'research_w_continuation/vendor')]
    import numpy as np
    from wr.probability_family import Family
    from wr.endpoint_geometry import linear_material_values
    from p0.geometry import wall_values
    import types
    family = Family.load(ROOT)
    nx, ny = 480, 320
    xs = -2 + (np.arange(nx)+.5)*4/nx
    ys = (np.arange(ny)+.5)*3/ny
    xx, yy = np.meshgrid(xs, ys)
    grid = np.stack([xx, yy], axis=-1)
    material = np.stack([family.material(b, 0., 3., grid, f=1.) for b in (0, 1)])
    decoder = types.SimpleNamespace(wall_values=wall_values)
    walls = []
    for q in family.walls:
        er, sigma = linear_material_values(2026102202, q, grid.reshape(-1, 2), decoder)
        walls.append((er+1j*sigma).reshape(xx.shape))
    walls = np.stack(walls)
    union = np.any(abs(walls-1)>1e-12, axis=0)
    def metrics(fields, mask=None):
        delta = fields[1]-fields[0]
        if mask is not None: delta = delta[mask]
        return dict(complex_RMS=float(np.sqrt(np.mean(abs(delta)**2))),
                    complex_MAE=float(np.mean(abs(delta))),
                    real_RMS=float(np.sqrt(np.mean(delta.real**2))),
                    conductivity_RMS=float(np.sqrt(np.mean(delta.imag**2))))
    out = ROOT/'tables/public_world_distance'
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(out/'MATERIALS.npz', x=xs, y=ys, canonical_joint=material,
                        calibration_walls=walls, wall_union_mask=union)
    result = dict(identity='PUBLIC_FROZEN_JOINT_ENDPOINT_DECODER_ONLY',
        frequency=1., grid=dict(shape=[ny,nx], domain=[-2,2,0,3],
                               sampling='equal-area cell centres; no continuous quadrature bound'),
        canonical_coordinates=[family.coordinates(b,0.,3.) for b in (0,1)],
        joint_full_domain=metrics(material), wall_full_domain=metrics(walls),
        wall_union_ROI=metrics(walls,union), wall_union_sampled_area=float(union.mean()*12),
        target_translation=.1, unchanged_target_mass=True,
        original_minimum_feature='UNRESOLVED', geometry_support_unchanged_for_wall=True,
        new_Maxwell_RHS=0, new_LU=0, GT_read=False,
        CPU_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        costs='pure decoder/postprocessing; separate from paid physics ledger')
    (out/'RESULTS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,ensure_ascii=False))


if __name__ == '__main__':
    main()
