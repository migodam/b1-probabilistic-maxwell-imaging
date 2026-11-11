"""Short public AP-pair → actual DR forecasts, with observed error balls.

The ball radius comes from the ACTUAL ray/trajectory deviation. Its branch
certificate is therefore a posterior diagnostic, not a future guarantee.
No event proposal, optimizer repair, truth scoring, or inherited source edit.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from pathlib import Path
import json
import time
import numpy as np
from scipy.stats import chi2
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r13.runtime.field_current import DenseFieldCurrentProjectors
from research.research_max_struct_r13.src.material_events import finite_threshold, angular_branch_margin, norm_pair


def _cost_delta(after, before):
    return {k:after.get(k, 0)-before.get(k, 0) for k in set(after)|set(before)}


def forecast(P, x0, q0, steps=30):
    wall, cpu = time.perf_counter(), time.process_time()
    before = dict(P.costs)
    cx, cq, m0 = P.material(x0, q0)
    dx, dq = x0-cx, q0-cq
    ax, aq = P.affine_ball(cx, cq)
    gap = norm_pair(dx, dq)
    pair_error = norm_pair(ax-x0, aq-q0)
    th = finite_threshold(cx, cq, dx, dq, m0, P.scale, P.lower, P.upper, P.components)
    threshold = th['first_t']
    predicted_step = max(0, int(np.floor(threshold))) if np.isfinite(threshold) else None
    mid = (np.arctan(P.lower/P.scale)+np.arctan(P.upper/P.scale))/2
    original_high = np.arctan(m0/P.scale) > mid
    # Descriptive, predeclared size: one quarter of the legal material span.
    # This is not an acceptance gate or the finite-threshold event itself.
    large_jump = .25*(P.upper-P.lower)
    x, q = x0.copy(), q0.copy()
    materials, ray_materials, margins, radii, rows = [], [], [], [], []
    first_flip = first_large = first_ray_flip = None
    for k in range(steps+1):
        ux, uq, material = P.material(x, q)
        rx, rq = cx+(k+1)*dx, cq+(k+1)*dq
        _, _, ray_material = P.material(rx, rq)
        error_x, error_q = x-rx, q-rq
        radius = np.sqrt(np.sum(abs(error_x.reshape(P.cells, -1))**2, axis=1)+
                         np.sum(abs(error_q.reshape(P.cells, -1))**2, axis=1))
        cert = angular_branch_margin(rx, rq, P.lower, P.upper, P.scale, P.components,
                                     perturbation_radius=radius)
        actual_high = np.arctan(material/P.scale) > mid
        ray_high = np.arctan(ray_material/P.scale) > mid
        robust = cert['robust_high']|cert['robust_low']
        disagreements = (cert['robust_high'] & ~actual_high)|(cert['robust_low'] & actual_high)
        flips = actual_high != original_high
        ray_flips = ray_high != original_high
        jumps = abs(material-m0) >= large_jump
        if k and first_flip is None and np.any(flips):
            first_flip = k
        if k and first_ray_flip is None and np.any(ray_flips):
            first_ray_flip = k
        if k and first_large is None and np.any(jumps):
            first_large = k
        rows.append(dict(DR_step=k, ideal_ray_t=k+1,
            ray_deviation=norm_pair(error_x, error_q),
            ray_deviation_over_initial_normal_gap=norm_pair(error_x, error_q)/max(gap, 1e-300),
            actual_branch_flips=int(np.count_nonzero(flips)), ray_branch_flips=int(np.count_nonzero(ray_flips)),
            large_material_jump_cells=int(np.count_nonzero(jumps)),
            max_material_change=float(np.max(abs(material-m0))),
            robust_branch_cells=int(np.count_nonzero(robust)), posterior_certificate_disagreements=int(np.count_nonzero(disagreements)),
            first_threshold_cell_actual_material=float(material[th['first_cell']]),
            first_threshold_cell_ray_material=float(ray_material[th['first_cell']]),
            first_threshold_cell_margin=float(cert['margin'][th['first_cell']]),
            first_threshold_cell_uncertainty=float(cert['uncertainty'][th['first_cell']]),
            first_threshold_cell_actual_deviation_radius=float(radius[th['first_cell']]),
            first_threshold_cell_robust=bool(robust[th['first_cell']])))
        materials.append(material)
        ray_materials.append(ray_material)
        margins.append(cert['margin'])
        radii.append(radius)
        if k < steps:
            arx, arq = P.affine_ball(2*ux-x, 2*uq-q)
            x, q = x+arx-ux, q+arq-uq
    finite = th['thresholds'][np.isfinite(th['thresholds'])]
    report = dict(steps=steps, first_ray_threshold_t=float(threshold) if np.isfinite(threshold) else None,
        threshold_first_cell=int(th['first_cell']), predicted_first_strict_posttie_DR_step=predicted_step,
        prediction_scope='floor(t*) future steps only on exact stationary normal ray; diagnostic hypothesis for approximate AP pair',
        actual_first_angular_half_branch_flip=first_flip, ray_first_angular_half_branch_flip=first_ray_flip,
        actual_first_large_material_jump=first_large,
        large_jump_definition='abs(m_k-m_0) >= one quarter of legal material interval; descriptive only',
        normal_gap=gap, pair_defect=pair_error, relative_pair_defect=pair_error/max(gap, 1e-300),
        amplitude_normal_defect=float(np.linalg.norm(th['amplitude_normal_defect'])),
        interior_stationarity_defect=float(np.linalg.norm(th['interior_stationarity_defect'])),
        finite_threshold_cells=len(finite), min_finite_threshold=float(np.min(finite)) if len(finite) else None,
        rows=rows, projector_cost=_cost_delta(P.costs, before),
        local_margin_evaluations=(steps+1)*P.cells,
        certificate_scope='a posterior actual-deviation state balls; no future trajectory or physical-material guarantee',
        CPU_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall)
    arrays = dict(base_A_X=x0, base_A_Q=q0, cX=cx, cQ=cq, dX=dx, dQ=dq,
        final_DR_X=x, final_DR_Q=q, material_trace=np.array(materials), ray_material_trace=np.array(ray_materials),
        ray_margin_trace=np.array(margins), actual_deviation_cell_radius_trace=np.array(radii),
        material_thresholds=th['thresholds'], threshold_a=th['a'], threshold_b=th['b'], threshold_gamma=th['gamma'])
    return report, arrays


def run_public(source, out, *, checkpoints=(50, 128, 512), steps=30):
    wall, cpu = time.perf_counter(), time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    with np.load(source, allow_pickle=False) as z:
        if any(any(word in k.lower() for word in ('truth', 'm_true', 'ground', 'label')) for k in z.files):
            raise ValueError('Truth-like field in public diagnostic input')
        public = {k:z[k].copy() for k in z.files}
    m0 = public.pop('m_init')
    Y = public['Y']
    sigma = float(public['sigma'])
    delta = float(public.get('delta', sigma*np.sqrt(chi2.ppf(.99, 2*Y.size)/2)))
    model = MaxwellModel(public)
    initial = model.state(m0)
    P = DenseFieldCurrentProjectors(model.G, model.E, model.S, Y, delta,
        scale=6., lower=.02, upper=20., grouped_ports=True)
    x, q = initial['j']/P.scale, initial['Q'].copy()
    records = []
    for k in range(1, max(checkpoints)+1):
        cx, cq, _ = P.material(x, q)
        x, q = P.affine_ball(cx, cq)
        if k in checkpoints:
            # forecast works on copies; AP continuation remains untouched.
            report, arrays = forecast(P, x, q, steps)
            report['AP_checkpoint'] = k
            np.savez_compressed(out/('AP_'+str(k)+'_STATE_AND_TRACE.npz'), **arrays)
            records.append(report)
    # Two physical LUs per case: original initializer and final material of the
    # longest-checkpoint forecast. Other traces are saved without extra scoring.
    end_material = arrays['material_trace'][-1]
    end = model.state(end_material)
    def metrics(m, state):
        residual = state['Y']-Y
        return dict(Phi=.5*float(np.vdot(residual, residual).real/np.vdot(Y, Y).real)+float(model.reg.evaluate(m)[0]),
            relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)),
            residual_over_delta=float(np.linalg.norm(residual)/delta))
    result = dict(status='PUBLIC_EARLY_PAIR_FORECAST_DIAGNOSTIC', truth_access='NONE', source=str(source),
        checkpoints=list(checkpoints), actual_DR_steps_per_checkpoint=steps,
        current_state_initialization='same physical initializer; AP continuation does not consume forecast states',
        scope='3 already exposed public cases; early-pair hypothesis limits, not blind reconstruction or frozen identity retest',
        initial_physics=metrics(m0, initial), final_longest_forecast_physics=metrics(end_material, end),
        forecasts=records, physical_model_cost=dict(model.cost), projector_cost=dict(P.costs),
        CPU_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall,
        timer_scope='input/model/projector setup, 512AP plus three30DR forecasts, trace serialization and two physical audits; final report serialization excluded')
    (out/'FORECAST.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, action='append', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if len(a.input) > 3:
        raise ValueError('Bounded forecast diagnostic accepts at most three public cases')
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        for source in a.input:
            r = run_public(source.resolve(), a.output/source.stem)
            print(json.dumps(dict(case=source.stem, CPU_seconds=r['CPU_seconds'],
                forecasts=[dict(AP=f['AP_checkpoint'], t=f['first_ray_threshold_t'],
                    relative_pair_defect=f['relative_pair_defect'], first_actual_branch=f['actual_first_angular_half_branch_flip'],
                    last_ray_error=f['rows'][-1]['ray_deviation_over_initial_normal_gap'],
                    posterior_disagreements=sum(v['posterior_certificate_disagreements'] for v in f['rows'])) for f in r['forecasts']])), flush=True)


if __name__ == '__main__':
    main()
