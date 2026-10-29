"""Predeclared public 2x2 chart mechanism audit. No private inputs or PDE on import.

--freeze writes a public protocol only. Parent schedules --run through
run_suite.py / guarded_online.py after granting the serial PDE queue.
"""
from pathlib import Path
import argparse
import json
import os
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
OWN = ROOT / 'research/delegated/sources'
FREEZE = OWN / 'CHART_MECHANISM_FROZEN.json'


def protocol():
    return dict(
        schema='R3_PUBLIC_CHART_MECHANISM_2X2_V1',
        scope='Matched public states; mechanism attribution, not blind recovery accuracy',
        config=dict(version='R3_CONTINUOUS_FULL_MACRO_3X2_V1', wall_x=[-1.5, -.5, .5, 1.5],
                    wall_y=[.5, .875, 1.25], frequencies=[.85, 1., 1.15], sigma=.01,
                    h=.17, targetstep=.10, wall_refine=1, target_degree=8,
                    phase='R3_CHART_MECHANISM_PUBLIC', pixel_nx=241, pixel_ny=181),
        uniform_wall_eps=[3.] * 6, uniform_wall_loss=[.06] * 6,
        changed_wall_eps=[2.1, 3.8, 2.7, 4.2, 2.4, 3.6],
        changed_wall_loss=[.03, .09, .04, .08, .05, .10],
        target_geometry=dict(x=.13, y=2.17, rx=.137, ry=.093, angle=.27),
        absent_target=dict(contrast=0., loss=0.), strong_target=dict(contrast=6.8, loss=.075),
        states=['uniform_absent', 'changed_absent', 'uniform_strong', 'changed_strong'],
        reference='uniform_absent', chart_rank=16, chart_updates=False,
        prior=dict(wall_eps_mean=3., wall_eps_sd=1.5, wall_loss_mean=.06,
                   wall_loss_sd=.15, target_precision=0.),
        step_damping=1e-3, step_damping_scope='1e-3 I on ALL scaled parameter coordinates; not posterior precision',
        metrics=dict(jacobian='complex Frobenius and real-white parameter-scaled Frobenius, whole/wall/target',
                     data_schur='B^T Q(A,LambdaW) B; LambdaX=0; no algorithm damping',
                     step_schur='Schur of J^T J + Lambda + 1e-3 I, all scaled coordinates',
                     step='unconstrained diagnostic GN step; no feasible or accepted reconstruction claim',
                     normalization='relative norms divide by full counterpart; zero denominator represented null'),
        observation=str(ROOT / 'inputs/online/dev_r0/observation.npz'),
        observation_scope='measured development y only; no truth, scoring, generator, or metadata reads',
        budget=dict(max_full_fe=18, max_rhs=600, expected_full_fe=12, expected_rhs=432,
                    rss_bytes=4294967296, physical_concurrency=1),
        no_monte_carlo=True, no_rank_sweep=True, no_pde_until_parent_queue_grant=True,
        output=str(OWN / 'CHART_MECHANISM_RUN'))


def write_json(path, value):
    with Path(path).open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def freeze():
    OWN.mkdir(parents=True, exist_ok=True)
    p = protocol()
    if FREEZE.exists():
        if json.loads(FREEZE.read_text()) != p:
            raise RuntimeError('Existing protocol differs; preserve and seek parent decision')
        print('Existing identical protocol retained:', FREEZE)
        return
    with FREEZE.open('x') as f:
        json.dump(p, f, indent=2, allow_nan=False)
        f.write('\n')
    print('FROZEN; no physical calculation launched:', FREEZE)


def run():
    if not FREEZE.exists():
        raise RuntimeError('Freeze public protocol before any physical action')
    frozen = json.loads(FREEZE.read_text())
    if frozen != protocol():
        raise RuntimeError('Frozen protocol differs from this implementation')
    r2 = ROOT.parent / 'research_w_gq_r2'
    for p in [ROOT / 'src', r2 / 'src', r2 / 'reuse/previous_gq/src',
              r2 / 'reuse/fixed_vendor', ROOT.parents[2] / 'research-w-continuation/.deps']:
        sys.path.append(str(p))
    import numpy as np
    from r3.accounting import BudgetBook, BudgetExceeded, fem_lock
    from r3.physics import JointModel, NAMES_TARGET
    from r3.charts import ForcingChart
    from r3.inverse import prior, whiten

    class LocalCapBook(BudgetBook):
        def __init__(self):
            super().__init__(root=ROOT)
            self.local_rhs = self.local_fe = 0

        def reserve(self, phase, label, rhs=0, factorizations=0, **meta):
            limit = frozen['budget']
            if self.local_rhs + rhs > limit['max_rhs'] or self.local_fe + factorizations > limit['max_full_fe']:
                raise BudgetExceeded('Public 2x2 mechanism audit local physical budget exhausted')
            super().reserve(phase, label, rhs=rhs, factorizations=factorizations,
                            diagnostic='CHART_MECHANISM_2X2', **meta)
            self.local_rhs += int(rhs)
            self.local_fe += int(factorizations)

    def rel(a, b):
        denom = float(np.linalg.norm(b))
        return float(np.linalg.norm(a) / denom) if denom else None

    def spectral_summary(a):
        vals = np.linalg.eigvalsh((a + a.T) / 2)
        return dict(eigenvalues=vals.tolist(), min_eigenvalue=float(vals[0]),
                    max_eigenvalue=float(vals[-1]), frobenius_norm=float(np.linalg.norm(a)))

    def schur(J, lamw, damping=0.):
        A, B = J[:, :12], J[:, 12:]
        hw = A.T @ A + lamw + damping * np.eye(12)
        cross = A.T @ B
        hx = B.T @ B + damping * np.eye(7)
        return hx - cross.T @ np.linalg.solve(hw, cross), hw, cross

    out = Path(frozen['output'])
    out.mkdir(parents=True, exist_ok=False)
    write_json(out / 'FROZEN_PROTOCOL.json', frozen)
    book = LocalCapBook()
    before = book.snapshot(pid=os.getpid())
    started, cpu = time.perf_counter(), time.process_time()
    status, error = 'RUNNING', None
    summaries = {}
    try:
        book.check(new_job=True)
        model = JointModel(frozen['config'], book)
        scale = np.asarray(model.spec['scale'], float)
        y = np.load(frozen['observation'])['y']
        expected_shape = (len(frozen['config']['frequencies']), 4, 16)
        if y.shape != expected_shape:
            raise ValueError('Measured y shape differs from frozen public acquisition')
        states = {}
        for key in frozen['states']:
            changed, strong = key.startswith('changed'), key.endswith('strong')
            mat = frozen['strong_target' if strong else 'absent_target']
            geom = frozen['target_geometry']
            eps = frozen['changed_wall_eps' if changed else 'uniform_wall_eps']
            loss = frozen['changed_wall_loss' if changed else 'uniform_wall_loss']
            states[key] = np.asarray(eps + loss + [geom['x'], geom['y'], geom['rx'], geom['ry'],
                                                  mat['contrast'], mat['loss'], geom['angle']], float)
        # Reference construction has only the uniform wall + zero target material.
        # Geometry is held identical in all four states for matched attribution.
        refz = states[frozen['reference']]
        wall, target = model.decode(refz)
        with fem_lock(book.root):
            for ev in model.evaluators:
                st = ev.evaluate(wall, target, label='public 2x2 reference state')
                ch = ForcingChart(ev, st, rank=frozen['chart_rank'], max_rank=frozen['chart_rank'])
                model.charts.append(ch)
                np.savez_compressed(out / f'REFERENCE_f{ev.p.freq:g}.npz', z=refz,
                                    parameter_scales=scale, chart_U=ch.U,
                                    **ch.reference_audit)
        np.savez_compressed(out / 'PUBLIC_OBSERVATION.npz', y=y, parameter_scales=scale)
        for key in frozen['states']:
            z = states[key]
            mark = book.snapshot(pid=os.getpid())
            full = model.forward(z, jacobian=True, label='public 2x2 ' + key)
            if key == frozen['reference']:
                fixed_j = np.vstack([ch.reference_audit['approximate_jacobian'] for ch in model.charts])
                fixed_data = full.data.copy()
                fixed_reference_reused = True
            else:
                fixed = model.forward(z, jacobian='FIXED_CURRENT_CHART', label='public fixed ' + key)
                fixed_j, fixed_data = fixed.jacobian, fixed.data
                fixed_reference_reused = False
            # This guard is semantic: chart comparisons must share one full response.
            # It does not create an extra forward solve or a machine-precision test.
            if not np.array_equal(full.data, fixed_data):
                raise RuntimeError('Full/fixed chart response arrays do not refer to the same state')
            J = np.vstack((full.jacobian.real, full.jacobian.imag)) * scale[None, :] / model.config['sigma']
            Jh = np.vstack((fixed_j.real, fixed_j.imag)) * scale[None, :] / model.config['sigma']
            E = Jh - J
            rp, jp = prior(model, z)
            jp = jp * scale[None, :]
            precision = jp.T @ jp
            lamw = precision[:12, :12]
            Sd, _, _ = schur(J, lamw)
            Sdh, _, _ = schur(Jh, lamw)
            damp = frozen['step_damping']
            Ss, hw, cross = schur(J, lamw, damp)
            Ssh, _, _ = schur(Jh, lamw, damp)
            r = whiten(full.data - y, model.config['sigma'])
            H = J.T @ J + precision + damp * np.eye(len(z))
            Hh = Jh.T @ Jh + precision + damp * np.eye(len(z))
            g = J.T @ r + jp.T @ rp
            gh = Jh.T @ r + jp.T @ rp
            delta = -np.linalg.solve(H, g)
            deltah = -np.linalg.solve(Hh, gh)
            normal_residual = H @ deltah + g
            profiled_residual = normal_residual[12:] - cross.T @ np.linalg.solve(hw, normal_residual[:12])
            f29 = np.linalg.solve(Ss, profiled_residual)
            # F27 with fixed wall prior, no target prior or numerical damping.
            A, B = J[:, :12], J[:, 12:]
            EA, EB = E[:, :12], E[:, 12:]
            invsqrt_lamw = np.diag(1 / np.sqrt(np.diag(lamw)))
            abar, eabar = A @ invsqrt_lamw, EA @ invsqrt_lamw
            qbound = (2 * np.linalg.norm(abar, 2) + np.linalg.norm(eabar, 2)) * np.linalg.norm(eabar, 2)
            f27_bound = 2 * np.linalg.norm(B, 2) * np.linalg.norm(EB, 2) + np.linalg.norm(EB, 2)**2 + np.linalg.norm(B, 2)**2 * qbound
            errors = {}
            for block, sl in [('whole', slice(None)), ('wall', slice(0, 12)), ('target', slice(12, None))]:
                errors[block] = dict(scaled_white_relative=rel(E[:, sl], J[:, sl]),
                                     scaled_white_absolute=float(np.linalg.norm(E[:, sl])),
                                     complex_raw_relative=rel(fixed_j[:, sl] - full.jacobian[:, sl], full.jacobian[:, sl]))
            after = book.snapshot(pid=os.getpid())
            summary = dict(state=key, physical_cost_since_reference_setup={k: after[k] - mark[k] for k in mark},
                           fixed_reference_arrays_reused=fixed_reference_reused, jacobian_errors=errors,
                           data_schur=dict(full=spectral_summary(Sd), fixed=spectral_summary(Sdh),
                                           relative_frobenius_error=rel(Sdh - Sd, Sd),
                                           measured_operator_error=float(np.linalg.norm(Sdh - Sd, 2)),
                                           F27_operator_bound=float(f27_bound)),
                           step_schur=dict(full=spectral_summary(Ss), fixed=spectral_summary(Ssh),
                                           relative_frobenius_error=rel(Ssh - Ss, Ss), damping_all_coordinates=damp,
                                           interpretation='algorithmic update matrix, not posterior precision'),
                           step=dict(unconstrained=True, relative_whole_error=rel(deltah - delta, delta),
                                     relative_target_error=rel(deltah[12:] - delta[12:], delta[12:]),
                                     target_error_norm=float(np.linalg.norm(deltah[12:] - delta[12:])),
                                     F29_predicted_target_error_norm=float(np.linalg.norm(f29)),
                                     normal_residual_norm=float(np.linalg.norm(normal_residual))))
            np.savez_compressed(out / (key + '.npz'), z=z, data_full=full.data, data_fixed=fixed_data,
                                jacobian_complex_full=full.jacobian, jacobian_complex_fixed=fixed_j,
                                jacobian_scaled_white_full=J, jacobian_scaled_white_fixed=Jh,
                                parameter_scales=scale, prior_precision_scaled=precision,
                                data_schur_full=Sd, data_schur_fixed=Sdh,
                                step_schur_full=Ss, step_schur_fixed=Ssh,
                                step_hessian_full=H, step_hessian_fixed=Hh,
                                step_gradient_full=g, step_gradient_fixed=gh,
                                step_full=delta, step_fixed=deltah, residual_white=r,
                                normal_residual=normal_residual, F29_target_error=f29,
                                F27_operator_bound=np.asarray(f27_bound))
            summaries[key] = summary
            write_json(out / (key + '.json'), summary)
            write_json(out / 'SUMMARY_PARTIAL.json', summaries)
            print(key, 'saved; cost', after, flush=True)
        # Matched deltas are descriptive factorial contrasts, not a physical
        # decomposition into independently additive wall and target scattering.
        state_errors = {k: summaries[k]['jacobian_errors']['whole']['scaled_white_relative'] for k in summaries}
        contrasts = dict(wall_effect_absent=state_errors['changed_absent'] - state_errors['uniform_absent'],
                         wall_effect_strong=state_errors['changed_strong'] - state_errors['uniform_strong'],
                         target_effect_uniform=state_errors['uniform_strong'] - state_errors['uniform_absent'],
                         target_effect_changed=state_errors['changed_strong'] - state_errors['changed_absent'])
        contrasts['interaction'] = contrasts['wall_effect_strong'] - contrasts['wall_effect_absent']
        write_json(out / 'SUMMARY.json', dict(scope=frozen['scope'], states=summaries,
                                             relative_error_factorial_contrasts=contrasts,
                                             caveat='Contrasts of relative errors; not additive scattering amplitudes or blind accuracy'))
        status = 'COMPLETE'
    except BaseException as exc:
        status, error = 'FAILED_PARTIAL_RETAINED', str(exc)
        book.failure('R3_CHART_MECHANISM_PUBLIC', 'public 2x2 mechanism audit', exc)
        raise
    finally:
        after = book.snapshot(pid=os.getpid())
        receipt = dict(status=status, error=error, pid=os.getpid(),
                       physical_cost={k: after[k] - before[k] for k in before},
                       process_cpu_seconds=time.process_time() - cpu,
                       wall_seconds=time.perf_counter() - started,
                       completed_states=list(summaries), protocol=str(FREEZE),
                       rss_guard='BudgetBook checks + required run_suite.py polled 4GiB wrapper',
                       queue_guard='Parent serial scheduling + shared .fem.lock for all FE actions')
        write_json(out / 'RECEIPT.json', receipt)
        book.event('R3_CHART_MECHANISM_PUBLIC', 'mechanism_audit_finished',
                   **{k: v for k, v in receipt.items() if k != 'pid'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--freeze', action='store_true')
    group.add_argument('--run', action='store_true', help='Parent queue grant required; use guarded wrapper')
    args = parser.parse_args()
    freeze() if args.freeze else run()


if __name__ == '__main__':
    main()
