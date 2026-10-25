#!/usr/bin/env python3
"""Build paid explicit finite-prior bank; importing runs no forward action.

Run only through the serialized stage wrapper:
  python tools/run_job.py E2E3 finite-bank python tools/finite_bank.py --output PATH
Resume uses --resume; uncertain paid rows are never automatically retried.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import tempfile
import time
import numpy as np
from wslam.accounting import BudgetBook, BudgetExceeded, json_default
from wslam.finite_bayes import build_grid, PRIOR_ID
from wslam.forward import Evaluator
from wslam.physics import acquisition
from wslam.world import World

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'E2E3'
TRUTH_SEED = 2026102504
NOISE_SEED = 2026102505


def write_json(path,data):
    path = Path(path)
    fd,tmp = tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as stream:
            json.dump(data,stream,indent=2,allow_nan=False,default=json_default)
            stream.write('\n');stream.flush();os.fsync(stream.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def save_npz(path,**arrays):
    path = Path(path)
    fd,tmp = tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as stream:
            np.savez_compressed(stream,**arrays)
            stream.flush();os.fsync(stream.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def config(grid,tx,rx,ledger):
    supplement_path=ROOT/'FINITE_PILOT_SUPPLEMENT.json'
    supplement=json.loads(supplement_path.read_text()) if supplement_path.exists() else None
    return dict(finite_pilot_supplement=supplement,schema='FINITE_BANK_CHECKPOINT_V1',prior_id=PRIOR_ID,
        state_count=grid['state_count'],stage=STAGE,ledger=str(ledger),
        frequency=1.,tx=tx.tolist(),rx=rx.tolist(),h=.17,polygon_n=256,
        field_grid_refinement=2,order=3,target_intrinsic=dict(radius=.12,contrast=3.,loss=.02),
        known_gamma=[1.,0.],noise_sigma_real_imag=.01,parents=12,noise_repeats=2,
        truth_seed=TRUTH_SEED,noise_seed=NOISE_SEED,
        explicit_discrete_prior=True,continuous_A0_quadrature_certificate=False,
        planned_rhs=11800,planned_factorizations=2950)


def load_rows(directory,total):
    """Read the exact durable prefix; no extra solver action or reservation."""
    predictions=np.empty((total,4,16),complex)
    residuals=np.empty(total,float)
    prefix=0
    while (directory/f'{prefix:04d}.npz').exists() and prefix<total:
        with np.load(directory/f'{prefix:04d}.npz',allow_pickle=False) as row:
            if int(row['row_index'])!=prefix or row['data'].shape!=(4,16):
                raise RuntimeError('durable bank-row schema/index mismatch')
            if not np.isfinite(row['data']).all() or not np.isfinite(row['residual']):
                raise RuntimeError('durable bank row contains nonfinite values')
            predictions[prefix]=row['data'];residuals[prefix]=float(row['residual'])
        prefix+=1
    if any(int(path.stem)>=prefix for path in directory.glob('[0-9][0-9][0-9][0-9].npz')):
        raise RuntimeError('bank rows are not an exact completed prefix')
    return predictions,residuals,prefix


def checkpoint(online,data,residuals,prefix,grid,tx,rx,book):
    save_npz(online/'FINITE_BANK.npz',data=data[:prefix],residuals=residuals[:prefix],
        completed=np.array(prefix),total=np.array(grid['state_count']),tx=tx,rx=rx,
        frequency=np.array(1.),prior_id=np.array(PRIOR_ID))
    receipt=dict(completed=prefix,total=grid['state_count'],costs=book.totals(STAGE),
        prediction_only=True,truth_solve_count=0,
        statement='Explicit finite prior; no continuous A0 integration certificate or scientific gate')
    write_json(online/'COST_RECEIPT.json',receipt)
    print(json.dumps(receipt,default=json_default),flush=True)


def make_observations(online,offline,grid,data,cfg):
    """Draw prior worlds only after bank completion; reuse their paid means."""
    truth_rng=np.random.default_rng(TRUTH_SEED)
    indices=truth_rng.choice(grid['state_count'],size=12,p=grid['prior_weight'])
    gamma=complex(*cfg['known_gamma'])
    means=gamma*data[indices]
    noise_rng=np.random.default_rng(NOISE_SEED)
    noise=.01*(noise_rng.normal(size=(12,2,4,16))+1j*noise_rng.normal(size=(12,2,4,16)))
    observed=means[:,None,:,:]+noise
    # Truth indices and realized noise belong only to offline scoring.
    save_npz(offline/'FINITE_PRIOR_TRUTHS_AND_NOISE.npz',truth_index=indices,
        independent_noise=noise,clean_means=means,parent_id=np.arange(12),
        truth_seed=np.array(TRUTH_SEED),noise_seed=np.array(NOISE_SEED),
        known_gamma=np.array(gamma),prior_id=np.array(PRIOR_ID))
    base_ids=np.array([[f'E2E3-parent-{p:02d}-repeat-{r}-base' for r in range(2)] for p in range(12)])
    outer_ids=np.array([[f'E2E3-parent-{p:02d}-repeat-{r}-outer' for r in range(2)] for p in range(12)])
    sequential_ids=np.array([[[f'E2E3-parent-{p:02d}-repeat-{r}-t{t}' for t in (1,2,3)] for r in range(2)] for p in range(12)])
    remaining_rx=np.array([r for r in range(16) if r!=8],dtype=int)
    save_npz(online/'OBSERVATIONS.npz',data=observed,base=observed[:,:,:2,:],outer=observed[:,:,2:,:],
        sequential_factor_id=sequential_ids,sequence_t1_tx=np.array([0]),sequence_t1_rx=np.array([8]),
        sequence_t2_tx=np.array([1]),sequence_t2_rx=np.array([8]),sequence_t3_tx=np.array([0,1]),sequence_t3_rx=remaining_rx,
        parent_id=np.arange(12),repeat_id=np.arange(2),base_factor_id=base_ids,outer_factor_id=outer_ids,
        known_gamma=np.array(gamma),known_target_radius=np.array(.12),known_target_contrast=np.array(3.),
        known_target_loss=np.array(.02),sigma_real_imag=np.array(.01),prior_id=np.array(PRIOR_ID))
    write_json(online/'FACTOR_IDS.json',dict(base=base_ids.tolist(),outer=outer_ids.tolist(),
        sequential=sequential_ids.tolist(),sequential_channels=dict(t1=dict(tx=[0],rx=[8]),t2=dict(tx=[1],rx=[8]),t3=dict(tx=[0,1],rx=remaining_rx.tolist())),
        warning='Base full-block IDs and sequential t1/t2/t3 are alternative representations of the same base data; never insert both.',
        observation_dimensions_complex=dict(base=32,outer=32),
        rule='One likelihood per parent/repeat/block ID; refreshing a likelihood replaces its version, never inserts a second copy.',
        repeats='Independent noise replicates; analyze separately unless the declared task explicitly combines them.',
        outer_description='Additional source block available for declared sequential update; not a truth label.'))


def run(args):
    out=Path(args.output).resolve()
    if out.exists() and not args.resume:
        raise RuntimeError('output exists; --resume required to reuse paid-prefix artifacts')
    if args.resume and not out.exists():
        raise RuntimeError('--resume requires an existing output')
    out.mkdir(parents=True,exist_ok=True)
    online,offline=out/'online',out/'offline'
    rows=online/'paid_rows';fields=online/'representative_fields'
    for path in (online,offline,rows,fields):path.mkdir(parents=True,exist_ok=True)
    with (out/'BANK.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        book=BudgetBook(args.ledger)
        grid=build_grid();tx,rx=acquisition(outer=True)
        cfg=config(grid,tx,rx,book.path)
        if grid['state_count']!=2950 or tx.shape!=(4,2) or rx.shape!=(16,2):
            raise RuntimeError('frozen state/acquisition count changed')
        manifest_path=out/'CHECKPOINT.json'
        if manifest_path.exists():
            manifest=json.loads(manifest_path.read_text())
            if manifest.get('config')!=cfg:
                raise RuntimeError('resume config differs; no paid calls attempted')
        else:
            if args.resume:
                raise RuntimeError('missing original checkpoint; refuse inferred/free retry')
            manifest=dict(config=cfg,status='INITIALIZED',completed=0,pending_row=None,created_unix=time.time(),initial_costs=book.totals(STAGE))
            write_json(manifest_path,manifest)
            arrays={k:v for k,v in grid.items() if isinstance(v,np.ndarray)}
            arrays.update(prior_id=np.array(PRIOR_ID),explicit_discrete_prior=np.array(True),continuous_A0_quadrature_certificate=np.array(False))
            save_npz(online/'FINITE_GRID.npz',**arrays)
        data,residuals,prefix=load_rows(rows,2950)
        if manifest.get('pending_row') is not None and int(manifest['pending_row'])>=prefix:
            raise RuntimeError('uncertain or failed paid row without durable prediction; no automatic retry or refund')
        if prefix<int(manifest.get('completed',0)):
            raise RuntimeError('durable prefix shorter than recorded completed count; no automatic retry')
        manifest.update(completed=prefix,pending_row=None)
        if prefix==2950:
            # Offline algebra can be reconstructed deterministically without PDE work.
            checkpoint(online,data,residuals,prefix,grid,tx,rx,book)
            make_observations(online,offline,grid,data,cfg)
            manifest.update(status='COMPLETE',completed=prefix,final_costs=book.totals(STAGE))
            manifest.update(scientific_gate='NOT_RUN',truth_solve_count=0,representative_field_files=[str(p.relative_to(out)) for p in sorted(fields.glob('ROW_*.npz'))])
            write_json(manifest_path,manifest)
            write_json(out/'COMPLETE.json',manifest)
            return
        needed=2950-prefix
        available=book.snapshot(STAGE)[STAGE]['remaining']
        if available['rhs']<needed*4 or available['factorizations']<needed:
            raise BudgetExceeded('remaining shared E2E3 cap cannot finish the exact bank; no calls made')
        solver=Evaluator(book,STAGE,h=.17,n=256,order=3,grid_refinement=2,tx=tx,rx=rx)
        representatives={0,2950//2,2949}
        manifest.update(status='RUNNING')
        write_json(manifest_path,manifest)
        try:
            for i in range(prefix,2950):
                before=book.totals(STAGE)
                manifest.update(pending_row=i,completed=i,pre_action_costs=before)
                write_json(manifest_path,manifest)
                world=World(grid['support'][i].reshape(2,2),float(grid['epsinf'][i]),float(grid['ell'][i]),grid['pose'][i],radius=.12,contrast=3.,target_loss=.02)
                ev=solver.evaluate(world,f=1.,jacobian=False)
                prediction=np.asarray(ev.data,complex)
                residual=float(ev.state.residual)
                if prediction.shape!=(4,16) or not np.isfinite(prediction).all() or not np.isfinite(residual):
                    raise RuntimeError('invalid paid prediction/residual')
                after=book.totals(STAGE)
                receipt=dict(row_index=i,before=before,after=after,
                    rhs_charged=after['rhs']-before['rhs'],factorizations_charged=after['factorizations']-before['factorizations'],
                    all_physics_reserved_by='shared BudgetBook through Evaluator/FEM',metadata=ev.metadata)
                # Persist each completed paid row before the aggregate 25-row checkpoint.
                save_npz(rows/f'{i:04d}.npz',row_index=np.array(i),data=prediction,residual=np.array(residual),receipt_json=np.array(json.dumps(receipt,default=json_default)))
                data[i]=prediction;residuals[i]=residual;prefix=i+1
                if i in representatives:
                    save_npz(fields/f'ROW_{i:04d}.npz',row_index=np.array(i),data=prediction,
                        material_points=ev.state.material_points,epsilon=ev.state.epsilon,field=ev.state.field,
                        mesh_p=ev.backend.mesh.p,mesh_t=ev.backend.mesh.t)
                manifest.update(completed=prefix,pending_row=None,last_costs=after)
                write_json(manifest_path,manifest)
                if prefix%25==0 or prefix==2950:
                    checkpoint(online,data,residuals,prefix,grid,tx,rx,book)
                del ev
        except BaseException as exc:
            # Keep all reserved fees; do not retry an uncertain or failed row.
            manifest.update(status='PARTIAL_FAILURE',completed=prefix,error=repr(exc),failure_costs=book.totals(STAGE))
            write_json(manifest_path,manifest)
            checkpoint(online,data,residuals,prefix,grid,tx,rx,book)
            raise
        finally:
            solver.clear()
        make_observations(online,offline,grid,data,cfg)
        manifest.update(status='COMPLETE',completed=2950,pending_row=None,final_costs=book.totals(STAGE),
            truth_solve_count=0,scientific_gate='NOT_RUN',
            representative_field_rows=sorted(representatives),
            representative_field_files=[str(p.relative_to(out)) for p in sorted(fields.glob('ROW_*.npz'))])
        write_json(manifest_path,manifest)
        write_json(out/'COMPLETE.json',manifest)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default=str(ROOT/'results/e2e3/finite_bank'))
    parser.add_argument('--ledger',default=None,help='Same shared project ledger as the stage wrapper')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args(argv)
    run(args)


if __name__=='__main__':
    main()
