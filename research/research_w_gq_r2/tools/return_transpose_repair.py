"""One parent-authorized complex-symmetric projection repair, without rank growth.

The first E3 study remains intact.  Reuse its paid Q, exact target templates,
receiver maps and direct references.  Pay for H conj(Q), then replace only
Q core Q* by Q [Q* H conj(Q)] Q^T.  No Bayesian cold run or certificate is
created.  Default invocation is PLAN_ONLY; execute with the R2 job wrapper.
"""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np
from scipy.sparse import load_npz

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

from r2.accounting import BudgetBook, BudgetExceeded
from r2.return_model import (
    ClosureSolution, PublicPorts, TargetSelfTemplate, linear_defect,
    ratio, relative_error,
)
from return_operator import _validate, _write_rows, cost_delta, write_json


PHASE = 'E3'
RANKS = (4, 8, 16, 32)
PLAN = dict(
    study='single complex-symmetric transpose repair of frozen E3 representation',
    changed_factor='right projection Q^T instead of Q*; new core Q* H conj(Q)',
    unchanged='same paid per-wall Q32, ranks, targets, walls, exact free-target templates, Bw and full references',
    return_model='H_r=Q_r (Q_r* H conj(Q_r)) Q_r^T',
    whole_diagnostic='Cw_r=Q_r (Q_r* Cw conj(Q_r)) Q_r^T; same H-fitted basis',
    expected_rhs=dict(new_Gw_G0_conjugate_basis=128, true_feedback_validation=128, total=256),
    expected_full_fe=3, phase=PHASE, ranks=list(RANKS),
    operator_norm_certificate=False, cold_Bayesian=False, new_rank_search=False,
    first_attempt_preserved=True,
)


def transpose_closure(Q, core, left_basis, base, rank, method):
    """j=base+left_basis core (I-Q^T left_basis core)^-1 Q^T base."""
    q, c = Q[:, :rank], core[:rank, :rank]
    left = left_basis[:, :rank]@c
    closure = np.eye(rank, dtype=complex)-q.T@left
    rhs = q.T@base
    coefficients = np.linalg.solve(closure, rhs)
    current = base+left@coefficients
    if not np.isfinite(current).all():
        raise np.linalg.LinAlgError('nonfinite repaired closure current')
    singular = np.linalg.svd(closure, compute_uv=False)
    metadata = dict(method=method, rank=rank,
                    reduced_closure='I-Q_r^T left_basis core; transpose, not adjoint',
                    closure_singular_values=singular.tolist(),
                    closure_condition_2=ratio(singular[0], singular[-1]),
                    closure_inverse_norm_2=ratio(1., singular[-1]),
                    closure_norms='numerical small-coordinate matrix values; not certified physical resolvent bounds',
                    reduced_backward=linear_defect(closure, coefficients, rhs),
                    current_relation_relative=relative_error(current, base+left@(q.T@current)),
                    changed_factor='complex-symmetric right projection',
                    full_feedback_validation='separate paid true Gw current action')
    return ClosureSolution(current, closure, coefficients, rhs, metadata)


def _preflight(book):
    for scope, cap in ((PHASE, book.auth['caps'][PHASE]), (None, book.cap)):
        used = book.snapshot(scope)
        if (used['rhs']+PLAN['expected_rhs']['total'] > cap['rhs']
                or used['full_fe']+PLAN['expected_full_fe'] > cap['full_fe']):
            raise BudgetExceeded('remaining resource cap cannot cover one frozen transpose repair: '+str((scope, used)))


def _archive(path):
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key].copy() for key in data.files}


def execute(origin, output):
    from r2.physics import Eval, S0_TX
    from r2.worlds import public_wall

    origin, output = Path(origin).resolve(), Path(output).resolve()
    if ROOT/'results' not in origin.parents or ROOT/'results' not in output.parents:
        raise ValueError('origin and output must be inside the new R2 campaign results')
    original_summary = json.loads((origin/'SUMMARY.json').read_text())
    original_run = json.loads((origin/'RUN_DEFINITION.json').read_text())
    if original_summary['status'] != 'COMPUTED_UNCERTIFIED':
        raise ValueError('repair requires the fully retained computed first attempt')
    if (original_run['range_columns'] != 32 or tuple(original_run['ranks']) != RANKS
            or len(original_run['targets']) != 4 or original_run['wall_ids'] != [0, 2]):
        raise ValueError('frozen first-attempt dimensions changed')
    contract = json.loads((ROOT/'CONTRACT.json').read_text())
    if original_run['contract_frozen_unix'] != contract['frozen_unix']:
        raise ValueError('contract changed after the first study')
    book = BudgetBook(ROOT)
    _preflight(book)
    output.mkdir(parents=True, exist_ok=False)
    before = book.snapshot(PHASE)
    cpu_start, wall_start = time.process_time(), time.perf_counter()
    write_json(output/'RUN_DEFINITION.json', dict(**PLAN, origin=str(origin),
               original_status=original_summary['status'],
               contract_frozen_unix=contract['frozen_unix'],
               authorization='parent explicitly authorized one mechanism-driven transpose repair <=300 additional RHS'))
    book.event(PHASE, 'e3_transpose_repair_start', origin=str(origin), output=str(output), plan=PLAN)
    rows, walls = [], []
    try:
        public = public_wall(0)
        ev = Eval(book, PHASE, f=1., tx=S0_TX, wall_x=public.x, wall_y=public.y,
                  wall_refine=1, h=.17, targetstep=.125, max_factor_bytes=1_000_000_000)
        space = _archive(origin/'PUBLIC_SPACE.npz')
        if (not np.array_equal(space['free_global_dofs'], ev.p._free_dofs)
                or not np.array_equal(space['port_free_indices'], ev.ass.public_roi_free_indices)
                or not np.array_equal(space['mesh_vertices'], ev.p.mesh.p)
                or not np.array_equal(space['mesh_triangles'], ev.p.mesh.t)
                or not np.array_equal(space['tx'], ev.p.tx)
                or not np.array_equal(space['rx'], ev.p.rx)
                or not np.array_equal(space['direct'], ev.direct)):
            raise ValueError('replayed evaluator does not equal the saved public FE/source/receiver space')
        L0 = load_npz(origin/'L0.npz').tocsc()
        ports = PublicPorts(space['port_free_indices'], L0.shape[0])
        mass_ports = load_npz(origin/'ROI_L2_MASS.npz').tocsc()
        fresh_L0, _ = ev.free_operator()
        if (L0-fresh_L0).nnz:
            raise ValueError('saved and fresh free operators are not identical')
        free_factor, _, free_bytes = ev.factor(L0, 'E3 transpose repair free background sparse LU')
        for wall_id in original_run['wall_ids']:
            source = origin/f'wall_{wall_id}'
            directory = output/f'wall_{wall_id}'
            directory.mkdir()
            wall_before = book.snapshot(PHASE)
            basis = _archive(source/'RAW_RETURN_ACTIONS.npz')
            Q = basis['Q']
            if Q.shape != (ports.size, 32):
                raise ValueError('saved per-wall range basis is not the frozen 32-column basis')
            Lw = load_npz(source/'Lw.npz').tocsc()
            background = _archive(source/'BACKGROUND.npz')
            wall_factor, _, wall_bytes = ev.factor(Lw, f'E3 transpose repair wall {wall_id} sparse LU')
            load = ports.inject(Q.conj())
            wall_conjugate = ports.restrict(wall_factor.solve(
                load, label=f'E3 transpose repair wall {wall_id} Gw conjQ',
                purpose='one-factor complex-symmetric return repair'))
            free_conjugate = ports.restrict(free_factor.solve(
                load, label=f'E3 transpose repair wall {wall_id} G0 conjQ',
                purpose='one-factor complex-symmetric return repair'))
            del load
            Hconjugate = wall_conjugate-free_conjugate
            Hcore, Ccore = Q.conj().T@Hconjugate, Q.conj().T@wall_conjugate
            np.savez_compressed(directory/'RAW_TRANSPOSE_ACTIONS.npz',
                                Q=Q, Gw_conjugate_Q_ports=wall_conjugate,
                                G0_conjugate_Q_ports=free_conjugate,
                                H_conjugate_Q=Hconjugate, Hcore=Hcore, Ccore=Ccore,
                                singular_values_Hcore=np.linalg.svd(Hcore, compute_uv=False),
                                singular_values_Ccore=np.linalg.svd(Ccore, compute_uv=False))
            diagnostics = []
            for rank in RANKS:
                q, h, c = Q[:, :rank], Hcore[:rank, :rank], Ccore[:rank, :rank]
                diagnostics.append(dict(
                    rank=rank,
                    in_sample_H_Omega_relative=relative_error(q@h@(q.T@basis['Omega']), basis['H_Omega']),
                    measured_H_conjugate_Q_relative=relative_error(q@h, Hconjugate[:, :rank]),
                    measured_Cw_conjugate_Q_relative=relative_error(q@c, wall_conjugate[:, :rank]),
                    Hcore_transpose_symmetry_relative=relative_error(h, h.T),
                    Ccore_transpose_symmetry_relative=relative_error(c, c.T),
                    diagnostic_scope='paid setup actions only; no operator-norm certification'))
            write_json(directory/'RANK_ACTION_DIAGNOSTICS.json', diagnostics)
            for target_definition in original_run['targets']:
                target_id, role = target_definition['id'], target_definition['role']
                target_source = source/f'target_{target_id}_{role}'
                target_output = directory/f'target_{target_id}_{role}'
                target_output.mkdir()
                raw_template = _archive(target_source/'EXACT_SELF_TEMPLATE.npz')
                raw_reference = _archive(target_source/'DIRECT_REFERENCE.npz')
                V = load_npz(target_source/'Vtarget.npz').tocsc()
                v, bx = ports.target_block(V, raw_reference['bX'])
                if (not np.array_equal(bx, raw_template['bX_ports'])
                        or not np.array_equal(background['bw']+raw_reference['bX'], raw_reference['rhs'])):
                    raise ValueError('saved analytic forcing/reference RHS changed')
                template = TargetSelfTemplate(
                    raw_template['q'], raw_template['j0'], raw_template['KQ'], raw_template['VQ'],
                    v, bx, raw_template, dict(reused_exact_free_target=True, no_new_template_RHS=True))
                reference = SimpleNamespace(u=raw_reference['u'], data=raw_reference['data'])
                assembly = SimpleNamespace(Vtarget=V, bX=raw_reference['bX'])
                for rank in RANKS:
                    for method, core, left_basis, base in (
                            ('return_transpose_exact_target_self', Hcore, template.KQ, template.j0),
                            ('whole_C_transpose_same_return_basis_diagnostic', Ccore, template.VQ, template.q)):
                        row = dict(wall_id=wall_id, target_id=target_id, target_role=role,
                                   method=method, rank=rank)
                        algebra_cpu = time.process_time()
                        label = f'E3 transpose repair wall {wall_id} target {target_id} {method} rank {rank}'
                        try:
                            solution = transpose_closure(Q, core, left_basis, base, rank, method)
                            row['algebra_cpu_seconds'] = time.process_time()-algebra_cpu
                            arrays, metrics = _validate(
                                ev, wall_factor, ports, mass_ports, Lw, background['bw'],
                                background['uw'], background['d_wall'], background['Bw'],
                                assembly, template, reference, solution, label,
                                float(contract['physics']['sigma_real_imag']))
                            row.update(status='COMPUTED_UNCERTIFIED', closure=solution.metadata, metrics=metrics)
                            np.savez_compressed(target_output/f'{method}_r{rank:02d}.npz', **arrays)
                            del arrays, solution
                        except np.linalg.LinAlgError as error:
                            row.update(status='FAILED_NUMERICAL_CLOSURE', error=str(error),
                                       algebra_cpu_seconds=time.process_time()-algebra_cpu)
                            book.failure(PHASE, label, error)
                        rows.append(row)
                        write_json(target_output/f'{method}_r{rank:02d}.json', row)
                        write_json(output/'PER_RANK_ERRORS.json', rows)
                        _write_rows(output/'PER_RANK_ERRORS.csv', rows)
                del template, reference, assembly, raw_template, raw_reference, V
                gc.collect()
            walls.append(dict(wall_id=wall_id, background_factor_bytes=wall_bytes,
                              cost=cost_delta(wall_before, book.snapshot(PHASE))))
            del wall_factor, basis, background, Lw
            gc.collect()
        cost = cost_delta(before, book.snapshot(PHASE))
        status = ('COMPUTED_UNCERTIFIED' if all(r['status'] == 'COMPUTED_UNCERTIFIED' for r in rows)
                  else 'COMPUTED_WITH_NUMERICAL_FAILURES')
        summary = dict(status=status, plan=PLAN, origin=str(origin), original_attempt_preserved=True,
                       phase_cost=cost, cumulative_E3_resources=book.snapshot(PHASE),
                       free_dofs=ports.nfree, public_target_ports=ports.size,
                       free_factor_bytes=free_bytes, walls=walls, rows=len(rows),
                       failed_rows=sum(r['status'] != 'COMPUTED_UNCERTIFIED' for r in rows),
                       kernel_cpu_seconds=time.process_time()-cpu_start,
                       wall_seconds=time.perf_counter()-wall_start,
                       CPU_charge='outer R2 run_job receipt only; local CPU is informational',
                       certificates=dict(operator_norm=False, probability=False, continuous_PDE=False),
                       cold_Bayesian='NOT_RUN_SEPARATE_GATE',
                       no_general_ROM_superiority='whole-C comparator uses exactly the same H-fitted basis',
                       scientific_review='parent owns final interpretation of this one-factor mechanism repair')
        write_json(output/'SUMMARY.json', summary)
        book.event(PHASE, 'e3_transpose_repair_completion', output=str(output), status=status, cost=cost)
        print(json.dumps(dict(status=status, output=str(output), cost=cost, rows=len(rows))))
    except Exception as error:
        book.failure(PHASE, 'e3_transpose_repair', error, output=str(output))
        write_json(output/'FAILED.json', dict(status='FAILED', error=str(error),
                   phase_cost=cost_delta(before, book.snapshot(PHASE)), completed_rows=len(rows),
                   kernel_cpu_seconds=time.process_time()-cpu_start,
                   wall_seconds=time.perf_counter()-wall_start))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--origin', type=Path, default=ROOT/'results/e3_return_operator')
    parser.add_argument('--output', type=Path, default=ROOT/'results/e3_return_transpose_repair')
    parser.add_argument('--execute-frozen', action='store_true')
    args = parser.parse_args()
    if not args.execute_frozen:
        print(json.dumps(dict(**PLAN, status='PLAN_ONLY', PDE_actions=0), indent=2))
        return
    execute(args.origin, args.output)


if __name__ == '__main__':
    main()
