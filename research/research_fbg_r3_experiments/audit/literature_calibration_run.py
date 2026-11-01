"""Reproduce TIME-ONLY baseline calibration without opening material arrays."""
import json
from pathlib import Path
import numpy as np
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import run_cie


def main():
    root = Path(__file__).resolve().parents[1]
    with np.load(root/'data/dev_calibration.npz', allow_pickle=False) as f:
        public = {k: f[k] for k in ('G', 'E0', 'S', 'grid_shape')}
        Y = f['Y']
    model = MaxwellModel(public)
    spectrum = np.linalg.svd(model.S, compute_uv=False)
    base = dict(dpcs_rtol=1e-3, physical_sigma=np.linalg.norm(Y),
                regularizer=model.regularizer, miniter=2,
                coefficient_tol=1e-3, material_box=(.02, 20.))
    m0 = np.full(model.n, .1)
    timing = []
    for _ in range(3):
        out = run_cie(model, Y, m0, {**base, 'stages': [dict(mf=1, beta=1., maxiter=50)],
                                    'coefficient_tol': 0., 'miniter': 50})
        timing.append(dict(wall=out['cost']['wall_seconds'], cost=out['cost'],
                           iterations=sum('iteration' in h for h in out['history'])))
    definitions = [
        ('cie_single', 'cie', [dict(mf=1, beta=1., maxiter=300), dict(mf=2, beta=.2, maxiter=300)]),
        ('cie_dual', 'cie', [dict(mf=1, beta=1., maxiter=300), dict(mf=2, beta=.2, beta_current=1., maxiter=300)]),
        ('tsom', 'tsom', [dict(mf=1, maxiter=300), dict(mf=2, maxiter=300)]),
        ('csi', 'csi', [dict(mf=1, maxiter=600)]),
        ('cie_single_strict', 'cie', [dict(mf=1, beta=1., maxiter=500), dict(mf=2, beta=.2, maxiter=500)]),
        ('cie_dual_strict', 'cie', [dict(mf=1, beta=1., maxiter=500), dict(mf=2, beta=.2, beta_current=1., maxiter=500)]),
    ]
    rows = []
    for name, method, stages in definitions:
        stopping = dict(coefficient_tol=1e-8, miniter=50) if name.endswith('_strict') else {}
        out = run_cie(model, Y, m0, {**base, **stopping, 'stages': stages, 'method': method})
        rows.append(dict(name=name, stages=stages, status=out['status'], cost=out['cost'],
                         coefficient_tol=stopping.get('coefficient_tol', base['coefficient_tol']),
                         miniter=stopping.get('miniter', base['miniter']),
                         stage_ends=[h for h in out['history'] if h.get('event') == 'stage_end'],
                         candidates=[{k:v for k,v in c.items() if k != 'm'} for c in out['candidates']],
                         iterations=sum('iteration' in h for h in out['history'])))
    result = dict(dataset='data/dev_calibration.npz',
                  arrays_opened=['G', 'E0', 'S', 'grid_shape', 'Y'],
                  material_arrays_opened=False, purpose='TIME ONLY, no parameter optimization',
                  source_rank=int(np.count_nonzero(spectrum > 1e-3*spectrum[0])),
                  dpcs_rtol=1e-3, receiver_singular_values=spectrum.tolist(),
                  timing_50=timing, T0=float(np.median([t['wall'] for t in timing])),
                  full_schedule=rows)
    (root/'audit/literature_calibration_timing.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
