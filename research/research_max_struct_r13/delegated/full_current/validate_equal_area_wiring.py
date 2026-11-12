"""Actual Austria public observations; at most three new-wiring updates.

No truth load and no recovery experiment. Distinct test stopping thresholds
explicitly exercise paper-stage stopping versus original G19 global stopping.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[k] = '1'
from pathlib import Path
import sys
import json
import numpy as np
from threadpoolctl import threadpool_limits
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from research.research_max_struct_r13.src.physics_inputs import existing_r8_root
from research.research_max_struct_r8.src.paper_physics import acquisition, make_model
from research.research_max_struct_r13.src.full_current import run_full_current


def main():
    out = Path(__file__).resolve().parent/'equal_area_wiring_v2'
    out.mkdir(exist_ok=True)
    source = existing_r8_root()/'RESULTS/paper_benchmarks/austria_eps4_400MHz/PUBLIC_DATA.npz'
    with np.load(source, allow_pickle=False) as z:
        Y = z['Y'].copy()
        public = acquisition(int(z['inverse_side']), float(z['frequency']))
    reports = []
    for method, mode in [('cie', 'restricted'), ('g19', 'restricted'), ('g19', 'full_apcs'), ('g19', 'full_current')]:
        model = make_model(public)
        m0 = np.full(model.n, .1)
        sub = out/(method+'_'+mode)
        cfg = dict(arm='cie_single' if method == 'cie' else 'a_plus', method='cie',
            material_box=[0., 20.], sourcewise_cg=method == 'cie', physics_family='equal_area_disk',
            dpcs_rank=15, coefficient_tol=0., miniter=2,
            stages=[dict(mf=5, beta_material=8., beta_current=8., maxiter=3, coefficient_tol=1e10)],
            fft_geometry=dict(coords=public['centers'], radius=public['radius'], k=public['k']),
            regularizer=lambda m:model.reg.evaluate(m)[:2], aplus_regularizer=model.reg,
            physical_sigma=float(np.linalg.norm(Y)), initialize_physical_current=method == 'g19',
            receipt_path=str(sub/'EXECUTION.json'))
        result = run_full_current(model, Y, m0, cfg, mode=mode, g19=method == 'g19', audit_every=1)
        receipt = result['execution_receipt']
        assert receipt['operation_counts']['PR_CG_steps'] == (2 if method == 'cie' else 3)
        assert receipt['initialize_physical_current'] == (method == 'g19')
        assert ('paper_cie_fft' in receipt['inherited_loop_entrypoint']) == (method == 'cie')
        if method == 'cie':
            assert result['initial_current_projection'] is None
        else:
            assert result['initial_current_projection'] is not None
            if mode == 'full_current':
                assert result['initial_current_projection']['relative_current_loss'] < 1e-10
        reports.append(dict(method=method, mode=mode, receipt=receipt))
    (out/'VALIDATION.json').write_text(json.dumps(dict(status='NEW_EQUAL_AREA_WIRING_ONLY',
        source=str(source), truth_access='NONE', radius=public['radius'],
        qualification='Two/three deliberately shortened updates; no reconstruction or recovery claim', runs=reports), indent=2)+'\n')
    print(json.dumps(dict(status='NEW_EQUAL_AREA_WIRING_ONLY',
        runs=[dict(method=r['method'], mode=r['mode'], loop=r['receipt']['inherited_loop_entrypoint'],
            PR_steps=r['receipt']['operation_counts']['PR_CG_steps'],
            warm=r['receipt']['initialize_physical_current']) for r in reports]), indent=2))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
