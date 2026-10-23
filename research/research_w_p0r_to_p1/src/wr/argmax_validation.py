"""Independent grid check of the EXISTING fixed-24-task argmax; ORACLE only.

No reoptimization, probability/prior editing, or physics at import. The original
endpoint validation files are read-only; a new argmax_validation folder refuses
overwrite. New work is at most two 5-frequency F calls = 40 complex RHS.
"""
from __future__ import annotations

import argparse
import json
import time
import uuid
from pathlib import Path

import numpy as np
from p0.geometry import probes

from .accounting import ROOT,ScopedBook,BudgetStop,write
from .native_world import LinearWorld
from .physical_witness import ModulatedWorld,mesh_for
from .physics import Forward

IDENTITY='ORACLE_FIXED_24_TASK_ARGMAX_EMPIRICAL_GRID_CHECK'


def _rows(book):
    if not book.path.exists():return []
    return [json.loads(line) for line in book.path.read_text().splitlines() if line.strip()]


def load_completed_validation(manifest):
    """Read/verify complete endpoint artifacts without any Book or solve."""
    path=Path(manifest)
    if not path.is_absolute():path=ROOT/path
    cfg=json.loads(path.read_text())
    identifier=cfg['id']
    if Path(identifier).name!=identifier or identifier in ('.','..'):
        raise ValueError('invalid validation ID')
    directory=ROOT/'results/B2/validation'/identifier
    if json.loads((directory/'MANIFEST.json').read_text())!=cfg:
        raise RuntimeError('manifest differs from completed validation manifest')
    gate=json.loads((directory/'RESULTS.json').read_text())
    tasks=json.loads((directory/'ALL_TASKS.json').read_text())
    if (not isinstance(tasks,list) or len(tasks)!=24
            or {r['probe'] for r in tasks}!=set(range(24))
            or not all((directory/f'probe{i}_PAIR.npz').exists() for i in range(24))
            or not isinstance(gate.get('numerical_pass'),bool)
            or 'd_T_sup_24' not in gate or 'task_argmax' not in gate):
        raise RuntimeError('complete full endpoint validation required before argmax check')
    public=probes()
    by_index={r['probe']:r for r in tasks}
    for index,row in by_index.items():
        expected=public[index]
        if (not np.array_equal(row['center'],expected[0]) or row['radius']!=expected[1]
                or row['contrast']!=expected[2] or not np.isfinite(row['distance']) or row['distance']<0):
            raise RuntimeError('fixed task class or saved task distance is inconsistent')
    # Frozen probe order resolves exact ties, matching endpoint_validation.
    index=max(range(24),key=lambda i:by_index[i]['distance'])
    if gate['task_argmax']!=index or not np.isclose(gate['d_T_sup_24'],by_index[index]['distance'],rtol=1e-12,atol=1e-12):
        raise RuntimeError('saved argmax/gate disagrees with complete fixed-task table')
    with np.load(directory/f'probe{index}_PAIR.npz',allow_pickle=False) as archive:
        fine0,fine1=archive['data0'].copy(),archive['data1'].copy()
    distance=float(np.linalg.norm(fine1-fine0)/.01)
    if not np.isclose(distance,by_index[index]['distance'],rtol=1e-12,atol=1e-12):
        raise RuntimeError('saved fine pair does not reproduce declared argmax distance')
    with np.load(ROOT/cfg['x_file'],allow_pickle=False) as archive:q=archive['x'].copy()
    reused=cfg['kind']=='cavity' and index==8
    label='joint8' if reused else f'probe{index}'
    for world,expected in ((0,np.zeros_like(q)),(1,q)):
        recorded=json.loads((directory/f'fields/{label}_world{world}_fine/endpoint_0000/PARAMETERS.json').read_text())
        if not np.array_equal(np.asarray(recorded['x']),expected):
            raise RuntimeError('current latent differs from original validated endpoint')
    return cfg,directory,index,public[index],q,fine0,fine1,reused


def main(manifest):
    cfg,source,index,probe,q,fine0,fine1,reused=load_completed_validation(manifest)
    frequencies=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())['frequencies'][cfg['band']]
    shape=(len(frequencies),4,12)
    if len(frequencies)!=5 or fine0.shape!=shape or fine1.shape!=shape:
        raise RuntimeError('argmax check requires original five-frequency 4x12 acquisition')
    if not np.all(np.isfinite(fine0)) or not np.all(np.isfinite(fine1)):
        raise RuntimeError('saved fine endpoint contains nonfinite data')
    out=source/'argmax_validation'
    out.mkdir(parents=True,exist_ok=False)
    book=ScopedBook('B2','validation')
    phase=f'B2_task_argmax_grid_{cfg["id"]}_{uuid.uuid4().hex}'
    cpu0,wall0=time.process_time(),time.perf_counter()
    failure=None
    coarse0=coarse1=None
    result=None

    def receipt():
        rows=[r for r in _rows(book) if r.get('stage')=='B2'
              and r.get('scope')=='validation' and r.get('phase')==phase]
        reserves=[r for r in rows if r.get('event')=='reserve']
        return dict(scope='validation',new_rhs=sum(r.get('rhs',0) for r in reserves),
                    new_factorizations=sum(r.get('factorizations',0) for r in reserves),
                    reserve_ids=[r['id'] for r in reserves],
                    cpu_seconds=time.process_time()-cpu0,wall_seconds=time.perf_counter()-wall0,
                    cached_source_work_was_previously_paid=True,
                    cpu_billing='external wrapper charges process CPU once')
    write(out/'MANIFEST.json',cfg)
    write(out/'PLAN.json',dict(identity=IDENTITY,fixed_task_argmax=index,probe=probe,
          source_validation=str(source),reused_existing_joint8_grid=reused,
          max_new_rhs=0 if reused else 40,probability=None,
          prior_files_read_only=True,original_min_feature_status='UNRESOLVED'))
    try:
        if reused:
            with np.load(source/'joint8_PAIR.npz',allow_pickle=False) as archive:
                coarse0,coarse1=archive['coarse0'].copy(),archive['coarse1'].copy()
                if not np.array_equal(fine0,archive['fine0']) or not np.array_equal(fine1,archive['fine1']):
                    raise RuntimeError('reused joint8 pair differs from fixed probe8 fine pair')
            grid=[r for r in json.loads((source/'GRID_COMPARISON.json').read_text()) if r['mode']=='joint8']
            if len(grid)!=1:raise RuntimeError('saved joint8 grid comparison missing or ambiguous')
            recorded_grid=grid[0]
        else:
            used=sum(r.get('rhs',0) for r in _rows(book) if r.get('event')=='reserve'
                     and r.get('stage')=='B2' and r.get('scope')=='validation')
            if used+40>book.scope_caps['rhs']:
                raise BudgetStop('whole argmax grid check needs 40 remaining validation RHS')
            cls=LinearWorld if cfg.get('chart','linear')=='linear' else ModulatedWorld
            cavities=tuple((tuple(c[0]),c[1]) for c in cfg.get('cavities',[]))
            data=[]
            for world in (0,1):
                model=cls(cfg['seed'],probe=probe,cavities=cavities if world else ())
                mesh=mesh_for(cfg['seed'],probe,h=.145,tol=5e-5,res=2,n=256,
                              cavities=cavities if world else ())
                forward=Forward(model,mesh,frequencies,book,
                                str((out/f'fields/world{world}_coarse').relative_to(ROOT/'results')),
                                phase=phase,identity=IDENTITY,h=.145,order=3,pml_strength=6,
                                max_dofs=110000,max_factor_bytes=3500000000,retain_current=False)
                response,_=forward.evaluate(q if world else np.zeros_like(q),jac=False)
                data.append(response.copy())
                np.savez_compressed(out/f'WORLD{world}_COARSE.npz',data=response)
                write(out/'LIVE_COST.json',receipt())
            coarse0,coarse1=data
        changes=[float(np.linalg.norm(coarse0-fine0)/.01),float(np.linalg.norm(coarse1-fine1)/.01)]
        if reused and not all(np.isclose(changes[w],recorded_grid[f'world{w}_grid_distance'],rtol=1e-12,atol=1e-12) for w in (0,1)):
            raise RuntimeError('reused grid diagnostics disagree with saved raw responses')
        passed=bool(max(changes)<=.1)
        np.savez_compressed(out/'TASK_ARGMAX_GRID_CHECK.npz',coarse0=coarse0,coarse1=coarse1,
                            fine0=fine0,fine1=fine1,probe_index=index,original_q=q)
        result=dict(identity=IDENTITY,status='PASS_EMPIRICAL_ARGMAX_GRID_CHECK' if passed
                    else 'FAIL_EMPIRICAL_ARGMAX_GRID_CHECK',task_argmax=index,probe=probe,
                    fine_distance=float(np.linalg.norm(fine1-fine0)/.01),
                    coarse_distance=float(np.linalg.norm(coarse1-coarse0)/.01),
                    world0_grid_distance=changes[0],world1_grid_distance=changes[1],
                    max_world_grid_change=max(changes),task_argmax_numerical_pass=passed,
                    reused_existing_joint8_grid=reused,actual_cost=receipt(),
                    probabilities=None,original_validation_results_unchanged=True,
                    controlled_P1_prior_unchanged=True,min_feature_status='UNRESOLVED',
                    conditional_fixed_support_geometry=True,continuous_certification=False)
    except Exception as exc:
        failure=dict(error_type=type(exc).__name__,error=str(exc),no_recovery_solve=True)
        result=dict(identity=IDENTITY,status='INCOMPLETE_ARGMAX_GRID_CHECK',task_argmax=index,
                    task_argmax_numerical_pass=False,failure=failure,actual_cost=receipt(),probabilities=None)
    finally:
        write(out/'TASK_ARGMAX_GRID_CHECK.json',result)
        write(out/'RHS_CPU_WALL_RECEIPT.json',receipt())
        write(out/('COMPLETE.json' if failure is None else 'INCOMPLETE.json'),result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest')
    args=parser.parse_args()
    result=main(args.manifest)
    if not result['task_argmax_numerical_pass']:raise SystemExit(2)
