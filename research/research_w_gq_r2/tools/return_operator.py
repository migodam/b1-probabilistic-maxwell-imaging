"""Frozen E3 representation qualification; the default invocation makes no PDE.

Execute only after the parent releases the frozen experiment, using the R2
run_job wrapper and --execute-frozen.  This is independent of cold Bayesian
deployment.  It supplies no probability, operator-norm, angular-spectrum,
continuous-Maxwell, or out-of-family generalization certificate.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse import save_npz

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

from r2.accounting import BudgetBook, BudgetExceeded, plain
from r2.return_model import (
    PublicPorts, build_return_approximation, build_target_self, linear_defect,
    measure_field, relative_error, solve_return, solve_whole_same_basis,
)


PHASE = 'E3'
SEED = 28102931
COLUMNS = 32
RANKS = (4, 8, 16, 32)
WALL_IDS = (0, 2)
PLAN = dict(
    experiment='E3 return-only representation qualification',
    phase=PHASE, wall_ids=list(WALL_IDS), target_selection='2 train + first 2 heldout',
    frequency=1., sources=2, receivers=16,
    range_columns=COLUMNS, oversampling=0, seed=SEED, ranks=list(RANKS),
    public_ports='all target-ROI FE free DOFs, independent of current target',
    basis_metric='Euclidean FE coefficients; deterministic complex QR in Gaussian column order',
    target_self='exact free-target sparse FEM factor L0-VX',
    ordinary_diagnostic='whole Cw projected into the same H-fitted basis; not optimal whole-C ROM',
    mass_metric='fixed target-ROI P3 L2 Gram matrix for field errors only',
    full_C0_cache=False, operator_norm_certificate=False,
    periodic_angular_claim=False, cold_Bayesian_run=False,
    estimated_rhs=dict(range_and_core=256, wall_primal=4, receiver_adjoint=32,
                       exact_target_self=272, direct_full_reference=16,
                       full_feedback_validation=128, total=708),
    estimated_full_fe=dict(free=1, wall=2, target_self=8, direct_full_reference=8, total=19),
    target_basis_reuse='separate Q per wall; each exact target template is rebuilt and charged',
    output_default='results/e3_return_operator',
)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=plain, allow_nan=False)+'\n')


def cost_delta(before, after):
    # run_job writes its once-only CPU charge after this child exits.  A zero
    # in-job CPU-ledger delta must not be presented as zero actual CPU work.
    return {key: after[key]-before[key] for key in ('rhs', 'full_fe')}


def _roi_mass(ev, ports, box):
    """Fixed FE field L2 Gram matrix; no whitening/factorization/PDE solve."""
    from skfem import BilinearForm, asm

    @BilinearForm(dtype=np.float64)
    def mass(u, v, w):
        return w.roi*u*v

    points = ev.p.material_points
    xl, xr, yb, yt = box
    inside = ((points[:, 0] >= xl) & (points[:, 0] <= xr)
              & (points[:, 1] >= yb) & (points[:, 1] <= yt))
    matrix = asm(mass, ev.p.material_basis, roi=inside.reshape(ev.p._material_shape))
    matrix = matrix[ev.p._free_dofs][:, ev.p._free_dofs]
    return matrix[ports.indices][:, ports.indices].tocsc()


def _preflight(book):
    expected_rhs = PLAN['estimated_rhs']['total']
    expected_fe = PLAN['estimated_full_fe']['total']
    for scope, cap in ((PHASE, book.auth['caps'][PHASE]), (None, book.cap)):
        used = book.snapshot(scope)
        if used['rhs']+expected_rhs > cap['rhs'] or used['full_fe']+expected_fe > cap['full_fe']:
            raise BudgetExceeded('insufficient remaining budget for frozen E3 return study: '+str((scope, used)))


def _receiver_map(ev, wall_factor, ports, label):
    """Bw=R Gw P^T via exactly 16 shared receiver-adjoint RHS."""
    adjoint = wall_factor.solve(ev.R.conj().T.toarray(), trans='H',
                                label=label+' shared receiver adjoints',
                                purpose='Bw=R Gw P^T; no per-port Green cache')
    result = adjoint[ports.indices].conj().T.copy()
    del adjoint
    return result


def _response_metrics(data, reference, wall_data, from_field, sigma):
    error = np.asarray(data)-reference
    white = np.r_[error.real.ravel(), error.imag.ravel()]/sigma
    return dict(total_response_relative=relative_error(data, reference),
                target_response_relative=relative_error(data-wall_data, reference-wall_data),
                white_response_norm=float(np.linalg.norm(white)),
                white_response_RMS=float(np.sqrt(np.mean(white**2))),
                max_complex_response_error=float(np.max(np.abs(error))),
                receiver_map_vs_paid_reconstruction_relative=relative_error(data, from_field),
                receiver_map_vs_paid_reconstruction_absolute=float(np.linalg.norm(data-from_field)))


def _validate(ev, wall_factor, ports, mass_ports, wall_matrix, wall_rhs, wall_u,
              wall_data, Bw, assembly, template, reference, solution, label, sigma):
    """Paid true Gw action is used only for validation, never for fitting."""
    started_cpu, started_wall = time.process_time(), time.perf_counter()
    delta = wall_factor.solve(ports.inject(solution.current),
                               label=label+' true Gw current validation',
                               purpose='validation only; full true feedback, no refit')
    model_u = wall_u+delta
    data = wall_data+(Bw@solution.current).T
    from_field = ev.direct+np.asarray(ev.R@model_u).T
    full_rhs = wall_rhs+assembly.bX
    current_from_field = template.q+template.Vports@ports.restrict(delta)
    true_current = template.q+template.Vports@ports.restrict(reference.u-wall_u)
    full_defect = linear_defect(wall_matrix-assembly.Vtarget, model_u, full_rhs)
    metrics = dict(**_response_metrics(data, reference.data, wall_data, from_field, sigma),
                   **measure_field(ports, mass_ports, model_u, reference.u, wall_u),
                   current_relative=relative_error(solution.current, true_current),
                   true_feedback_current_relative=relative_error(solution.current, current_from_field),
                   true_feedback_current_absolute=float(np.linalg.norm(solution.current-current_from_field)),
                   full_FE_defect=full_defect,
                   validation_cpu_seconds=time.process_time()-started_cpu,
                   validation_wall_seconds=time.perf_counter()-started_wall,
                   validation_rhs=int(delta.shape[1]),
                   validation_role='true wall Green current action; no basis update or rank selection')
    arrays = dict(current=solution.current, u=model_u, delta=delta, data=data,
                  response_from_paid_field=from_field, current_from_full_feedback=current_from_field,
                  closure=solution.core, reduced_coefficients=solution.coefficients,
                  reduced_rhs=solution.reduced_rhs)
    return arrays, metrics


def _write_rows(path, rows):
    fields = ('wall_id', 'target_id', 'target_role', 'method', 'rank', 'status',
              'white_response_norm', 'white_response_RMS', 'target_response_relative',
              'ROI_target_induced_L2_relative', 'current_relative',
              'true_feedback_current_relative', 'full_FE_rhs_relative',
              'closure_condition_2', 'closure_inverse_norm_2',
              'algebra_cpu_seconds', 'validation_cpu_seconds', 'error')
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            metrics = row.get('metrics', {})
            closure = row.get('closure', {})
            flat = {**metrics, **closure, **{key: row.get(key) for key in
                    ('wall_id', 'target_id', 'target_role', 'method', 'rank', 'status',
                     'algebra_cpu_seconds', 'error')}}
            flat['full_FE_rhs_relative'] = metrics.get('full_FE_defect', {}).get('rhs_relative')
            writer.writerow({key: flat.get(key) for key in fields})


def execute(output):
    # Import the physical adapter only on the explicitly executing path.
    from r2.physics import Eval, S0_TX
    from r2.worlds import make_target, public_wall

    if not (ROOT/'CONTRACT.json').exists():
        raise RuntimeError('parent must freeze CONTRACT.json before E3 physics')
    contract = json.loads((ROOT/'CONTRACT.json').read_text())
    sigma = float(contract['physics']['sigma_real_imag'])
    definitions = ([(f'train_{i}', d) for i, d in enumerate(contract['train_targets'])]
                   + [(f'heldout_{i}', d) for i, d in enumerate(contract['heldout_targets'][:2])])
    if len(definitions) != 4:
        raise ValueError('frozen E3 requires exactly 2 train plus first 2 heldout targets')
    walls = [(key, public_wall(key)) for key in WALL_IDS]
    if any(not np.array_equal(w.x, walls[0][1].x)
           or not np.array_equal(w.y, walls[0][1].y) for _, w in walls):
        raise ValueError('all E3 walls must share one public common mesh')
    output = Path(output).resolve()
    if ROOT/'results' not in output.parents:
        raise ValueError('E3 artifacts must remain under this new campaign results directory')
    book = BudgetBook(ROOT)
    _preflight(book)
    output.mkdir(parents=True, exist_ok=False)
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    before = book.snapshot(PHASE)
    frozen_run = dict(**PLAN, walls=[dict(id=key, **wall.metadata()) for key, wall in walls],
                      targets=[dict(id=i, role=role, **d) for i, (role, d) in enumerate(definitions)],
                      contract_frozen_unix=contract['frozen_unix'],
                      E3_contract=contract['E3'], full_FE_reference='fresh direct solve for every target/wall',
                      PDE_status='STARTING', expected_cost_at_start=PLAN['estimated_rhs'])
    write_json(output/'RUN_DEFINITION.json', frozen_run)
    book.event(PHASE, 'e3_return_operator_start', output=str(output), plan=PLAN)
    rows, wall_summaries = [], []
    try:
        ev = Eval(book, PHASE, f=1., tx=S0_TX, wall_x=walls[0][1].x,
                  wall_y=walls[0][1].y, wall_refine=1, h=.17, targetstep=.125,
                  max_factor_bytes=1_000_000_000)
        L0, _ = ev.free_operator()
        ports = PublicPorts(ev.ass.public_roi_free_indices, L0.shape[0])
        mass_ports = _roi_mass(ev, ports, contract['physics']['target_ROI'])
        save_npz(output/'L0.npz', L0)
        save_npz(output/'ROI_L2_MASS.npz', mass_ports)
        np.savez_compressed(output/'PUBLIC_SPACE.npz', port_free_indices=ports.indices,
                            free_global_dofs=ev.p._free_dofs, mesh_vertices=ev.p.mesh.p,
                            mesh_triangles=ev.p.mesh.t, tx=ev.p.tx, rx=ev.p.rx,
                            direct=ev.direct)
        write_json(output/'SPACE_METADATA.json', dict(
            free_dofs=ports.nfree, global_dofs=ev.p.ndofs, public_target_ports=ports.size,
            public_target_box=contract['physics']['target_ROI'],
            mesh_version=ev.p.mesh.r2_public_key[0],
            mass='P3 FE L2 Gram on fixed public target ROI; field diagnostics only',
            QR_metric='Euclidean FE coefficients; no mass whitening',
            source_order='source-major receiver arrays', full_C0_cached=False,
            max_live_large_factors=3, per_factor_storage_guard_bytes=1_000_000_000))
        free_factor, L0, free_bytes = ev.factor(L0, 'E3 free background sparse LU')
        for wall_id, wall in walls:
            directory = output/f'wall_{wall_id}'
            directory.mkdir()
            wall_before = book.snapshot(PHASE)
            wall_cpu = time.process_time()
            Lw, bw = ev.background(wall)
            wall_factor, Lw, wall_bytes = ev.factor(Lw, f'E3 wall {wall_id} background sparse LU')
            uw = wall_factor.solve(bw, label=f'E3 wall {wall_id} background primal')
            d_wall = ev.direct+np.asarray(ev.R@uw).T
            Bw = _receiver_map(ev, wall_factor, ports, f'E3 wall {wall_id}')
            approximation = build_return_approximation(
                free_factor, wall_factor, ports, columns=COLUMNS, seed=SEED,
                label=f'E3 wall {wall_id}')
            save_npz(directory/'Lw.npz', Lw)
            np.savez_compressed(directory/'BACKGROUND.npz', bw=bw, uw=uw, d_wall=d_wall, Bw=Bw)
            np.savez_compressed(directory/'RAW_RETURN_ACTIONS.npz', **approximation.raw)
            write_json(directory/'RETURN_METADATA.json', approximation.metadata)
            rank_diagnostics = []
            for rank in RANKS:
                Q, Hcore, Ccore = approximation.rank(rank)
                rank_diagnostics.append(dict(
                    rank=rank,
                    in_sample_H_Omega_relative=relative_error(
                        Q@Hcore@(Q.conj().T@approximation.raw['Omega']),
                        approximation.raw['H_Omega']),
                    measured_H_Q_relative=relative_error(
                        Q@Hcore, approximation.raw['H_Q'][:, :rank]),
                    measured_Cw_Q_relative=relative_error(
                        Q@Ccore, approximation.raw['Gw_Q_ports'][:, :rank]),
                    diagnostic_scope='paid setup directions only; no independent validation or operator-norm certificate'))
            write_json(directory/'RANK_ACTION_DIAGNOSTICS.json', rank_diagnostics)
            for target_id, (role, definition) in enumerate(definitions):
                target_directory = directory/f'target_{target_id}_{role}'
                target_directory.mkdir()
                label = f'E3 wall {wall_id} target {target_id}'
                target = make_target(definition)
                assembly = ev.assemble_target(target)
                save_npz(target_directory/'Vtarget.npz', assembly.Vtarget)
                target_cpu = time.process_time()
                target_before = book.snapshot(PHASE)
                target_factor, _, target_bytes = ev.factor(
                    L0-assembly.Vtarget, label+' exact free-target sparse LU')
                template = build_target_self(ports, assembly, target_factor, approximation.Q,
                                               uw, label=label)
                target_cost = cost_delta(target_before, book.snapshot(PHASE))
                target_cpu = time.process_time()-target_cpu
                np.savez_compressed(target_directory/'EXACT_SELF_TEMPLATE.npz', **template.raw)
                write_json(target_directory/'TARGET_METADATA.json', dict(
                    target=target.metadata(), role=role, assembly=assembly.metadata,
                    exact_self=template.metadata, cost=target_cost,
                    kernel_cpu_seconds=target_cpu, factor_bytes=target_bytes))
                # The template arrays suffice for both methods; release its
                # full factor before constructing the independent reference.
                del target_factor
                reference = ev.evaluate(wall, target, label=label+' direct full reference')
                j_reference = template.q+template.Vports@ports.restrict(reference.u-uw)
                np.savez_compressed(target_directory/'DIRECT_REFERENCE.npz',
                                    u=reference.u, data=reference.data, rhs=reference.rhs,
                                    j=j_reference, bX=assembly.bX)
                write_json(target_directory/'REFERENCE_METADATA.json', dict(
                    residual=reference.residual, timings=reference.timings,
                    reference='fresh direct (Lw-VX)u=bw+bX sparse factor for this target/wall'))
                for rank in RANKS:
                    for method, solver in (
                            ('return_only_exact_target_self', solve_return),
                            ('whole_C_same_return_basis_diagnostic', solve_whole_same_basis)):
                        row = dict(wall_id=wall_id, target_id=target_id, target_role=role,
                                   method=method, rank=rank)
                        algebra_cpu = time.process_time()
                        try:
                            solution = solver(approximation, template, rank)
                            row['algebra_cpu_seconds'] = time.process_time()-algebra_cpu
                            arrays, metrics = _validate(
                                ev, wall_factor, ports, mass_ports, Lw, bw, uw,
                                d_wall, Bw, assembly, template, reference, solution,
                                label+f' {method} rank {rank}', sigma)
                            row.update(status='COMPUTED_UNCERTIFIED', closure=solution.metadata,
                                       metrics=metrics)
                            np.savez_compressed(target_directory/f'{method}_r{rank:02d}.npz', **arrays)
                            del arrays, solution
                        except np.linalg.LinAlgError as error:
                            row.update(status='FAILED_NUMERICAL_CLOSURE', error=str(error),
                                       algebra_cpu_seconds=time.process_time()-algebra_cpu)
                            book.failure(PHASE, label+f' {method} rank {rank}', error)
                        rows.append(row)
                        write_json(target_directory/f'{method}_r{rank:02d}.json', row)
                        write_json(output/'PER_RANK_ERRORS.json', rows)
                        _write_rows(output/'PER_RANK_ERRORS.csv', rows)
                ev.clear_cache()
                del reference, template, assembly
                gc.collect()
            wall_summaries.append(dict(wall_id=wall_id, background_factor_bytes=wall_bytes,
                                       kernel_cpu_seconds=time.process_time()-wall_cpu,
                                       cost=cost_delta(wall_before, book.snapshot(PHASE))))
            del wall_factor, approximation, Lw, bw, uw, d_wall, Bw
            gc.collect()
        after = book.snapshot(PHASE)
        status = ('COMPUTED_UNCERTIFIED' if all(r['status'] == 'COMPUTED_UNCERTIFIED' for r in rows)
                  else 'COMPUTED_WITH_NUMERICAL_FAILURES')
        summary = dict(status=status, scope=PLAN['experiment'], plan=PLAN,
                       free_dofs=ports.nfree, public_target_ports=ports.size,
                       free_factor_bytes=free_bytes, walls=wall_summaries,
                       phase_cost=cost_delta(before, after),
                       kernel_cpu_seconds=time.process_time()-start_cpu,
                       wall_seconds=time.perf_counter()-start_wall,
                       CPU_charge='outer R2 run_job receipt only; local CPU fields are informational',
                       rows=len(rows), failed_rows=sum(r['status'] != 'COMPUTED_UNCERTIFIED' for r in rows),
                       certificates=dict(operator_norm=False, probability=False,
                                         continuous_PDE=False, L2_operator_spectrum=False),
                       cold_Bayesian='NOT_RUN_SEPARATE_GATE',
                       interpretation='finite common-space representation diagnostic; final scientific review belongs to parent',
                       whole_C_limitation='same H-fitted QR basis, not optimal whole-C or all ordinary ROMs',
                       direct_full_comparison='direct references use 16 RHS/8 FE; setup and validation costs remain charged',
                       shared_setup=dict(rhs=292, full_fe=3),
                       return_self_setup=dict(rhs=272, full_fe=8),
                       direct_reference=dict(rhs=16, full_fe=8),
                       validation=dict(rhs=sum(r.get('metrics', {}).get('validation_rhs', 0) for r in rows)))
        write_json(output/'SUMMARY.json', summary)
        book.event(PHASE, 'e3_return_operator_completion', output=str(output), status=status,
                   cost=summary['phase_cost'], kernel_cpu_seconds=summary['kernel_cpu_seconds'])
        print(json.dumps(dict(status=status, output=str(output), cost=summary['phase_cost'], rows=len(rows))))
    except Exception as error:
        book.failure(PHASE, 'e3_return_operator', error, output=str(output))
        write_json(output/'FAILED.json', dict(status='FAILED', error=str(error),
                   phase_cost=cost_delta(before, book.snapshot(PHASE)),
                   kernel_cpu_seconds=time.process_time()-start_cpu,
                   wall_seconds=time.perf_counter()-start_wall, completed_rows=len(rows)))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute-frozen', action='store_true',
                        help='run the parent-released frozen study; use R2 run_job for CPU/RSS accounting')
    parser.add_argument('--output', type=Path, default=ROOT/'results/e3_return_operator')
    args = parser.parse_args()
    if not args.execute_frozen:
        print(json.dumps(dict(**PLAN, PDE_actions=0, status='PLAN_ONLY'), indent=2))
        return
    execute(args.output)


if __name__ == '__main__':
    main()
