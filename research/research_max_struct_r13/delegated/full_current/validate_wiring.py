"""One known public-state wiring exercise; no material truth or scoring.

Two updates per family do not constitute an inversion experiment. This only
checks new basis/audit plumbing, paid projection loss and callback interruption.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import sys
import json
import numpy as np
from threadpoolctl import threadpool_limits

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r8.src.benchmark_solver import arm_config
from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank
from research.research_max_struct_r13.src.full_current import run_full_current, run_audited_cie


def main():
    out = Path(__file__).resolve().parent/'validation'
    out.mkdir(exist_ok=True)
    source = REPO/'research/research_max_struct_r12/supplement_paired/inputs/A__r10_twod_0004.npz'
    with np.load(source, allow_pickle=False) as z:
        public = {k:z[k].copy() for k in z.files}
    m0 = public.pop('m_init')
    Y = public['Y']
    decision = noise_ball_rank(public['S'], Y, float(public['sigma']))
    reports = []
    for mode, arm in [('restricted', 'a_plus'), ('full_apcs', 'a_plus'), ('full_current', 'a_plus'), ('restricted', 'cie_dual')]:
        model = MaxwellModel(public)
        cfg = arm_config(arm, dict(host=.1, material_box=[.02, 20.], dpcs_rank=decision['rank']), Y, model.E)
        # Explicitly short NON-SCIENTIFIC wiring exercise; never presented as
        # a completed/frozen-budget reconstruction or recovery comparison.
        cfg['stages'] = [dict(cfg['stages'][0], maxiter=2)]
        cfg.update(fft_geometry=dict(coords=public['coords'], radius=float(public['radius']), k=float(public['k'])),
            regularizer=lambda m:model.reg.evaluate(m)[:2], aplus_regularizer=model.reg,
            physical_sigma=float(np.linalg.norm(Y)), receipt_path=str(out/(arm+'_'+mode)/'EXECUTION.json'))
        audits = []
        def callback(m, state, value, metadata):
            expected = model.physical_score(m, Y)  # cache hit: no extra solve
            assert abs(expected-value) <= 1e-14
            audits.append(metadata)
        fn = run_full_current if arm == 'a_plus' else run_audited_cie
        result = fn(model, Y, m0, cfg, mode=mode, audit_every=1, audit_callback=callback)
        assert len(audits) == 4  # initial, each update, inherited stage-end
        assert result['execution_receipt']['operation_counts']['PR_CG_steps'] == 2
        assert len(result['physical_audit_candidates']) == 4
        assert np.isfinite(result['m']).all()
        p = result['initial_current_projection']
        if mode == 'full_current':
            assert p['relative_current_loss'] < 1e-10
        elif p is not None:
            assert abs(p['relative_current_loss']-p['DPCS_projection_loss']) < 1e-10 if mode == 'full_apcs' else True
        space = result['current_spaces'][0]
        expected_rank = 576 if mode == 'full_current' else 576-decision['rank'] if mode == 'full_apcs' else 36
        assert space['independent_complex_spatial_rank'] == expected_rank
        if arm == 'cie_dual':
            assert p is None  # preserve original CIE's absent warm-start flag
        reports.append(dict(arm=arm, mode=mode, status=result['status'], spaces=result['current_spaces'],
            projection=p, audit_count=len(audits), receipt=result['execution_receipt'],
            selection=result['selection']))
    # Caller interruption after an already-paid periodic audit: no new loop or
    # PR state, and the finally receipt still carries the retained audit prefix.
    class Stop(BaseException):
        pass
    saved = []
    cfg['cost_callback'] = saved.append
    cfg['receipt_path'] = str(out/'interrupted'/'EXECUTION.json')
    def interrupt(m, state, value, meta):
        if meta['phase'] == 'periodic':
            raise Stop()
    try:
        run_audited_cie(model, Y, m0, cfg, audit_every=1, audit_callback=interrupt)
    except Stop:
        pass
    else:
        raise AssertionError('Caller interruption not propagated')
    assert saved[0]['status'] == 'INTERRUPTED_BY_CALLER'
    assert saved[0]['eligible_physical_audits'] == 2
    summary = dict(status='WIRING_VALIDATED_NOT_AN_EXPERIMENT', source=str(source.relative_to(REPO)),
        truth_access='NONE', dpcs_rule=decision, runs=reports, interruption_receipt=saved[0],
        qualification='One known 576-cell public state; two updates per family; no recovery/performance/novelty inference')
    (out/'WIRING_VALIDATION.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=summary['status'], dpcs_rank=decision['rank'],
        families=[dict(arm=r['arm'], mode=r['mode'], rank=r['spaces'][0]['independent_complex_spatial_rank'],
                       projection_loss=r['projection']['relative_current_loss'] if r['projection'] else None,
                       audits=r['audit_count']) for r in reports]), indent=2))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
