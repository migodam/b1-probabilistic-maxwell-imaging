"""Audited legacy controls. Preview by default; physical work needs --execute --stage.

Original inputs remain immutable. Parent owns campaign physics and accounting.
"""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import runpy
import shutil
import sys
import time
for _var in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_var]='1'
import numpy as np
from accounting import Book, P, write

K = P / 'inputs/Research_Z_ALG_Research_Kit'


def _contract():
    c = json.loads((P / 'EXPERIMENT_CONTRACT.json').read_text())
    assert c['E4_legacy']['single_new_order'] == 32
    assert c['resources']['conservative_actual_plus_closed_form_source_RHS_equivalents_cap'] == 30000
    assert c['resources']['per_run_hard_cap'] == 10000
    return c


@contextlib.contextmanager
def _linear_accounting(book):
    solve, inv = np.linalg.solve, np.linalg.inv
    def charged_solve(a, b):
        a, b = np.asarray(a), np.asarray(b)
        complex_system = np.iscomplexobj(a)
        columns = 1 if b.ndim == 1 else b.shape[-1]
        with book.op('legacy_complex_solve' if complex_system else 'legacy_real_algebra_solve',
                     rhs=int(columns) if complex_system else 0,
                     LU=1 if complex_system else 0, algebraic_real_factorizations=0 if complex_system else 1, real_algebra_RHS=int(columns) if not complex_system else 0) as action:
            answer = solve(a, b)
            book.save(action, A=a, RHS=b, solution=answer)
            return answer
    def charged_inv(a):
        a = np.asarray(a)
        with book.op('legacy_complex_inverse' if np.iscomplexobj(a) else 'legacy_real_algebra_inverse',
                     rhs=int(a.shape[-1]) if np.iscomplexobj(a) else 0,
                     LU=1 if np.iscomplexobj(a) else 0, algebraic_real_factorizations=0 if np.iscomplexobj(a) else 1) as action:
            answer = inv(a)
            book.save(action, A=a, inverse=answer)
            return answer
    np.linalg.solve, np.linalg.inv = charged_solve, charged_inv
    try:
        yield
    finally:
        np.linalg.solve, np.linalg.inv = solve, inv


def _run_copy(path, argv=None, init_globals=None):
    old_argv = sys.argv
    stdout = io.StringIO()
    try:
        sys.argv = [str(path)] + (argv or [])
        with contextlib.redirect_stdout(stdout):
            namespace = runpy.run_path(str(path), run_name='__main__', init_globals=init_globals)
    finally:
        sys.argv = old_argv
    path.with_suffix('.stdout.txt').write_text(stdout.getvalue())
    return namespace


def e0(book):
    out = P / 'audit/E0_replay'
    out.mkdir(parents=True, exist_ok=True)
    paths = ['zalg_algebra_e0.py', 'zalgo_probability/analytic_probability_checks.py',
             'zalgo_maxwell/two_cell_witness.py']
    ns = None
    import scipy.special as special
    import scipy.integrate as integrate
    callbacks={'hankel1_calls':0,'hankel1_elements':0,'quad_calls':0,'quad_callback_calls':0}
    original_hankel,original_quad=special.hankel1,integrate.quad
    def counted_hankel(*args,**kwargs):
        callbacks['hankel1_calls']+=1
        callbacks['hankel1_elements']+=int(np.broadcast_arrays(*args[:2])[0].size)
        return original_hankel(*args,**kwargs)
    def counted_quad(fn,*args,**kwargs):
        callbacks['quad_calls']+=1
        def callback(*values):
            callbacks['quad_callback_calls']+=1
            return fn(*values)
        return original_quad(callback,*args,**kwargs)
    for relative in paths:
        dst = out / relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(K / relative, dst)
        with book.op('legacy_E0_script' if 'witness' not in relative else 'legacy_two_cell_script',
                     note='CPU/memory measured; Green adaptive quadrature not Maxwell RHS'):
            if 'witness' in relative:
                special.hankel1,integrate.quad=counted_hankel,counted_quad
                try:
                    with _linear_accounting(book):
                        ns = _run_copy(dst)
                finally:
                    special.hankel1,integrate.quad=original_hankel,original_quad
                book.append(dict(event='GREEN_KERNEL_CALLBACK_COUNTS',**callbacks,
                                 note='All supplied two-cell Hankel calls, including geometry and adaptive Green quadrature; zero extra states'))
                write(out/'GREEN_CALLBACK_COUNTS.json',callbacks)
            else:
                _run_copy(dst)
    from e0_extra_algebra import compute
    with book.op('E0_extra_algebra_zero_Maxwell'):
        extra=compute()
        write(out/'EXTRA_ALGEBRA_QA.json',extra)
    # Independent arithmetic checks reuse saved states; zero new forward evaluations.
    x, G, S, J, Jswap = [ns[key] for key in ['x', 'G', 'S', 'J', 'J_swap']]
    state_residual = np.linalg.norm((np.eye(2)-np.diag(x)@G)@J-np.diag(x)@ns['E0'])
    swap_difference = np.linalg.norm(S@J-S@Jswap)/np.linalg.norm(S@J)
    minors = ns['Jcandidate'][:,0]*ns['Ecandidate'][:,1]-ns['Jcandidate'][:,1]*ns['Ecandidate'][:,0]
    qa = dict(state_residual_fro=float(state_residual), swap_relative_difference=float(swap_difference),
              dark_shared_material_minors_abs=np.abs(minors).tolist(),
              dark_shared_material_residual=float(np.linalg.norm(ns['R'])),
              lagrange_identity_error=float(np.max(np.abs(ns['lagrange_lhs']-ns['lagrange_rhs']))),
              independent_checks='Arithmetic on freshly saved original states; no additional physical evaluation',
              extra_algebra_all_passed=extra['all_passed'],
              gate=bool(extra['all_passed'] and state_residual < 1e-10 and swap_difference < 1e-10 and np.linalg.norm(ns['R']) > 1e-3))
    np.savez_compressed(out/'two_cell_full_states.npz', x=x, G=G, S=S, E0=ns['E0'], J=J, J_swap=Jswap,
                        y=ns['y'], y_swap=ns['y_swap'], E0_dark=ns['E0two'], Jdark=ns['Jdark'],
                        Jcandidate=ns['Jcandidate'], Ecandidate=ns['Ecandidate'], xhat=ns['xhat'],
                        ydark=ns['ydark'], yhat=ns['yhat'], Jhat=ns['Jhat'], R=ns['R'])
    write(out/'INDEPENDENT_QA.json',qa)
    write(out/'REPLAY_STATUS.json',dict(stage='E0', actual_complex_RHS=20,
          all_calls_charged='8 two-source forward, one two-column sensitivity, one 2x2 inverse',
          complete_protocol_E0=False, completion_qualification='Parent must review all replay and extra algebra gates; script-level execution is not protocol acceptance', algebra_controls=extra, callback_counts=callbacks, scope='Supplied E0 scripts cover only a subset; two-cell is positive control',
          source_files=paths, source_immutable=True, qa=qa))
    return qa


def _replace_once(source, old, new):
    if source.count(old) != 1:
        raise RuntimeError('Auditable source patch anchor mismatch: '+old[:100])
    return source.replace(old, new, 1)


def e4(book):
    out = P / 'audit/E4_legacy_replay'
    out.mkdir(parents=True,exist_ok=True)
    source_path=K/'zalgo_probability/e4_two_cell_posterior.py'
    original=source_path.read_text()
    patched=original.replace('hashlib.sha256(', 'hashlib.blake2b(').replace('source_witness_sha256','source_witness_blake2b').replace('source_script_sha256','source_script_blake2b')
    if original.count('hashlib.sha256(') != 2:
        raise RuntimeError('Unexpected original identity code')
    patched=_replace_once(patched, '    def likelihood(pred):',
        '    response_direct = _formula_wrapper(response_direct, "direct_current_analytic_response")\n'
        '    response_quotient = _formula_wrapper(response_quotient, "quotient_rational_response")\n\n'
        '    def likelihood(pred):')
    old='    results = {str(n): {"original_box": integrate_box(n),\n                        "ordered_quotient": integrate_ordered_quotient(n)} for n in levels}'
    new='    results = {"32": {"original_box": integrate_box(32), "ordered_quotient": integrate_ordered_quotient(32)},\n               "64": _frozen_64}\n'
    patched=_replace_once(patched,old,new)
    adapted=out/'e4_two_cell_posterior_adapted.py'
    adapted.write_text(patched)
    import difflib
    (out/'SOURCE_PATCH.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile=str(source_path),tofile=str(adapted))))
    frozen=json.loads((K/'zalgo_probability/e4_two_cell_posterior.json').read_text())
    batch=0
    def formula_wrapper(fn, kind):
        def wrapped(x1,x2):
            nonlocal batch
            aa,bb=np.broadcast_arrays(x1,x2)
            with book.op(kind,formula=2*int(aa.size),states=int(aa.size),sources_per_state=2) as action:
                pred=fn(x1,x2)
                book.save(action,x1=aa,x2=bb,predictions=pred)
                batch+=1
                return pred
        return wrapped
    target=out/'e4_two_cell_posterior.json'
    with _linear_accounting(book):
        _run_copy(adapted,['--n','32','--witness',str(K/'zalgo_maxwell/two_cell_witness.json'),
                           '--output',str(target)],
                  {'_formula_wrapper':formula_wrapper,'_frozen_64':frozen['results']['64']})
    result=json.loads(target.read_text())
    goals=_contract()['E4_legacy']
    diffs=[*result['n_to_2n_differences'].values(),result['box_vs_quotient_at_2n']]
    # The legacy 'at_2n' field is cached64 versus cached64. Add NEW32 cross-check explicitly.
    newbox=result['results']['32']['original_box'];newq=result['results']['32']['ordered_quotient']
    cross32={key:abs(newbox[key]-newq[key])/(abs(newq[key]) if key=='scaled_evidence' else 1)
             for key in ['scaled_evidence','mean_x1','mean_x2','event_x1_gt_x2','psi_x1_squared']}
    limits={'scaled_evidence':goals['evidence_relative_goal'],'mean_x1':goals['mean_absolute_goal'],
            'mean_x2':goals['mean_absolute_goal'],'event_x1_gt_x2':goals['event_goal'],'psi_x1_squared':goals['second_moment_goal']}
    gate=all(d[key]<=limits[key] for d in [*diffs,cross32] for key in limits)
    metadata=dict(new_levels=[32],reused_levels=[64],reference_source=str(K/'zalgo_probability/e4_two_cell_posterior.json'),
                  certification='New32 against supplied cached64 regression; not fresh64 and not interval integration certificate',
                  cross_new32_direct_vs_quotient=cross32,new32_vs_cached64=result['n_to_2n_differences'],
                  fresh_actual_complex_RHS=10,fresh_closed_form_source_RHS_equivalent=6164,
                  expected_fresh_combined=6174,response_batches=batch,gate=gate,
                  frozen_reference_runtime='Historical cached64 time is not fresh runtime',source_immutable=True)
    result['replay_provenance']=metadata
    result['numerical_reference_gate_passed']=gate
    write(target,result)
    write(out/'REPLAY_STATUS.json',metadata)
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['E0','E4'])
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    _contract()
    if not args.execute:
        print(json.dumps(dict(mode='PREVIEW_NO_EXECUTION',E0_actual_complex_RHS=20,
              E4_actual_complex_RHS=10,E4_formula_equivalent=6164,E4_combined=6174,
              command='python src/legacy_regression.py --execute --stage E0|E4')))
        return
    if not args.stage:
        parser.error('--execute requires --stage')
    for var in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        os.environ[var]='1'
    book=Book('legacy_'+args.stage)
    start_actual,start_formula=book.actual,book.formula
    try:
        result=e0(book) if args.stage=='E0' else e4(book)
        measured=(book.actual-start_actual,book.formula-start_formula)
        expected=(20,0) if args.stage=='E0' else (10,6164)
        if measured != expected:
            raise RuntimeError(f'Unexpected billed counts: {measured}, expected {expected}')
    except BaseException as exc:
        book.end('FAILED',error=repr(exc));raise
    else:
        book.end('COMPLETE',result=result)
        print(json.dumps(result,indent=2))

if __name__=='__main__':main()
