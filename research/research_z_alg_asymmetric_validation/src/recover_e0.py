"""Recover E0 artifacts from already paid exact cached solve outputs; no new RHS.

Default preview. --execute replays immutable source with fail-closed cached
solve/inverse replacement, never falling back to a linear solver.
"""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
for _v in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_v]='1'
import numpy as np
from accounting import P, Book, write
from e0_extra_algebra import compute


def identical(a,b):
    a,b=np.asarray(a),np.asarray(b)
    return a.dtype==b.dtype and a.shape==b.shape and a.tobytes(order='C')==b.tobytes(order='C')


def recovery(book):
    out=P/'audit/E0_replay'
    old=json.loads((out/'zalgo_maxwell/two_cell_witness.json').read_text())
    rows=[json.loads(line) for line in (P/'COST_LEDGER.jsonl').read_text().splitlines()]
    paid={r['action']:r for r in rows if r['run']=='legacy_E0' and r['event']=='BEGIN'}
    completed={r['action'] for r in rows if r['run']=='legacy_E0' and r['event']=='END' and r['status']=='COMPLETE'}
    cache=[]
    for row in rows:
        if row['run']!='legacy_E0' or row['event']!='SAVE':continue
        action=row['action']
        if action not in completed:raise RuntimeError('Cannot reuse incomplete paid action')
        with np.load(P/row['path']) as z:
            entry={key:z[key].copy() for key in z.files}
        entry.update(action=action,kind=paid[action]['kind'],path=row['path'])
        cache.append(entry)
    if len(cache)!=11:raise RuntimeError(f'Expected 11 saved algebra/complex actions, found {len(cache)}')
    reused_actual=sum(paid[e['action']]['actual_complex_RHS'] for e in cache)
    if reused_actual!=20:raise RuntimeError('Expected exactly 20 previously paid RHS')
    position=0
    matches=[]
    def take(kind,a,b=None):
        nonlocal position
        if position>=len(cache):raise RuntimeError('Cached lookup exhausted; new solve forbidden')
        entry=cache[position]
        if entry['kind']!=kind or not identical(a,entry['A']):
            raise RuntimeError(f'Cached A/kind bytes differ at lookup {position}; no solver fallback')
        if b is not None and not identical(b,entry['RHS']):
            raise RuntimeError(f'Cached RHS bytes differ at lookup {position}; no solver fallback')
        answer=entry['inverse' if b is None else 'solution'].copy()
        matches.append(dict(original_paid_action=entry['action'],cache_path=entry['path'],
                            A_bytes_exact=True,RHS_bytes_exact=True if b is not None else 'identity implicit',kind=kind))
        book.append(dict(event='CACHED_PAID_RESULT_REUSE',original_run='legacy_E0',**matches[-1],new_RHS=0))
        position+=1
        return answer
    def cached_solve(a,b):
        return take('legacy_complex_solve' if np.iscomplexobj(a) else 'legacy_real_algebra_solve',a,b)
    def cached_inv(a):
        return take('legacy_complex_inverse' if np.iscomplexobj(a) else 'legacy_real_algebra_inverse',a)
    original_solve,original_inv=np.linalg.solve,np.linalg.inv
    source=P/'inputs/Research_Z_ALG_Research_Kit/zalgo_maxwell/two_cell_witness.py'
    isolation=out/'serialization_recovery'
    isolation.mkdir(parents=True,exist_ok=True)
    namespace={'__name__':'__main__','__file__':str(isolation/'two_cell_witness.py')}
    captured=io.StringIO()
    with book.op('cached_namespace_reconstruction_zero_new_RHS',reuses_paid_complex_RHS=20,
                 note='Immutable source replay; solve/inv return exact cached outputs only; Green coefficient algebra CPU paid'):
        try:
            np.linalg.solve,np.linalg.inv=cached_solve,cached_inv
            with contextlib.redirect_stdout(captured):
                exec(compile(source.read_text(),str(source),'exec'),namespace)
        finally:
            np.linalg.solve,np.linalg.inv=original_solve,original_inv
    if position!=len(cache):raise RuntimeError('Not all cached entries consumed')
    (isolation/'reconstruction.stdout.txt').write_text(captured.getvalue())
    if namespace['out']!=old:raise RuntimeError('Reconstructed original output differs from previously paid witness')
    with book.op('E0_extra_algebra_serialization_recovery_zero_RHS'):
        extra=compute();write(out/'EXTRA_ALGEBRA_QA.json',extra)
    n=namespace
    with book.op('E0_recovered_state_arithmetic_QA_zero_RHS'):
        state=np.linalg.norm((np.eye(2)-np.diag(n['x'])@n['G'])@n['J']-np.diag(n['x'])@n['E0'])
        swap=np.linalg.norm(n['S']@n['J']-n['S']@n['J_swap'])/np.linalg.norm(n['S']@n['J'])
        minor=n['Jcandidate'][:,0]*n['Ecandidate'][:,1]-n['Jcandidate'][:,1]*n['Ecandidate'][:,0]
        qa=dict(state_residual_fro=float(state),swap_relative_difference=float(swap),
                dark_shared_material_minors_abs=np.abs(minor).tolist(),
                dark_shared_material_residual=float(np.linalg.norm(n['R'])),
                lagrange_identity_error=float(np.max(abs(n['lagrange_lhs']-n['lagrange_rhs']))),
                extra_algebra_all_passed=bool(extra['all_passed']),original_output_exact_match=True,
                independent_checks='Arithmetic on recovered paid states; no independent new physical solve',
                gate=bool(extra['all_passed'] and state<1e-10 and swap<1e-10 and np.linalg.norm(n['R'])>1e-3))
        np.savez_compressed(out/'two_cell_full_states.npz',
            **{key:n[key] for key in ['x','G','S','E0','J','J_swap','y','y_swap','Jdark','Jcandidate','Ecandidate','xhat','ydark','yhat','Jhat','R']},
            E0_dark=n['E0two'])
        write(out/'INDEPENDENT_QA.json',qa)
    status=dict(status='RECOVERED' if qa['gate'] else 'RECOVERED_QA_FAILED',
                recovery_run=book.run,original_run='legacy_E0',new_actual_RHS=0,new_formula_equivalents=0,
                reused_paid_actual_RHS=reused_actual,cache_entries_reused=position,
                exact_cache_matches=matches,original_output_exact_match=True,qa=qa,
                original_failed_log_preserved=True,original_stdout_and_witness_preserved=True,
                protocol_acceptance='Parent scientific review required; recovery is not fresh physical validation')
    write(out/'RECOVERY_STATUS.json',status)
    return status


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');args=p.parse_args()
    if not args.execute:
        print('PREVIEW: zero-new-RHS cached E0 recovery; add --execute to recover artifacts');return
    book=Book('E0_SERIALIZATION_RECOVERY')
    start=(book.actual,book.formula)
    try:
        result=recovery(book)
        if (book.actual,book.formula)!=start:raise RuntimeError('Recovery billed unexpected new RHS')
    except BaseException as exc:
        book.end('FAILED',error=repr(exc));raise
    else:
        book.end('COMPLETE',qa_gate=result['qa']['gate'],new_RHS=0)
        print(json.dumps({'status':result['status'],'new_RHS':0,'qa':result['qa']},indent=2))

if __name__=='__main__':main()
