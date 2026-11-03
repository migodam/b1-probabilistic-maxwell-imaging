"""PUBLIC-only adapter to unchanged Maxwell/CIE implementations and R4 refine.

API for existing remote_runner.worker: solve(public, method, config, backend='cpu').
No files are opened here. No generator, private mapping, scoring or GT is imported.
"""
from collections import Counter
from datetime import datetime, timezone
import os
import sys
import time
import traceback
_IMPORT_WALL = time.perf_counter()
_IMPORT_CPU = time.process_time()
import numpy as np
from threadpoolctl import threadpool_limits
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src import cie_tsom_baseline as native
from research.research_max_struct_r2.src.cie_global import run_cie as global_cie
from research.research_max_struct_r4.src.cie_block_absorption import run_block_cie
from research.research_max_struct_r4.src.cie_onehop import run_onehop_cie
from research.research_max_struct_r4.experiments.cooperative_paths import refine

ALLOW = {'G','E0','S','Y','grid_shape','shape','sigma','noise_radius','host','amplitude',
         'material_bounds','smooth_weight','tv_weight','tv_epsilon','prior_width','coords',
         'sources','receivers','k','radius'}
METHODS = ('onehop', 'block4_cached', 'global', 'native_extended')
_IMPORT_WALL = time.perf_counter()-_IMPORT_WALL
_IMPORT_CPU = time.process_time()-_IMPORT_CPU
_IMPORT_CHARGED = False


def peak_memory():
    """Process high-water mark, including worker lifetime before this task."""
    try:
        if os.name == 'nt':
            import ctypes
            from ctypes import wintypes
            class Counters(ctypes.Structure):
                _fields_ = [('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[
                    (k,ctypes.c_size_t) for k in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage',
                     'QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
            value = Counters(); value.cb = ctypes.sizeof(value)
            kernel = ctypes.WinDLL('kernel32'); kernel.GetCurrentProcess.restype = wintypes.HANDLE
            psapi = ctypes.WinDLL('psapi')
            psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
            if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(value),value.cb):
                raise OSError('GetProcessMemoryInfo failed')
            return dict(peak_rss_bytes=int(value.PeakWorkingSetSize), scope='process lifetime high-water mark, not per-task allocation')
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return dict(peak_rss_bytes=int(rss if sys.platform == 'darwin' else rss*1024),
                    scope='process lifetime high-water mark, not per-task allocation')
    except Exception as exc:
        return dict(peak_rss_bytes=None, unavailable_reason=repr(exc))


def delta(after, before):
    return dict(Counter(after)-Counter(before))


def metric(model, Y, m):
    state = model.state(m)
    residual = state['Y']-Y
    return dict(m=np.asarray(m).copy(), objective=float(.5*np.vdot(residual,residual).real/np.vdot(Y,Y).real+model.reg.evaluate(m)[0]),
                relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)))


def history_arrays(history):
    """Columnar primitive NPZ arrays; full histories stay outside receipt JSON."""
    arrays = {}
    for key in sorted({key for row in history for key in row}):
        values = [row.get(key) for row in history]
        if all(value is None or isinstance(value, (int,float,bool,np.number)) for value in values):
            arrays['history_'+key] = np.array([np.nan if value is None else value for value in values])
        elif all(value is None or isinstance(value,str) for value in values):
            arrays['history_'+key] = np.array(['' if value is None else value for value in values], dtype=str)
    return arrays


def solve(public, method, config, backend='cpu'):
    global _IMPORT_CHARGED
    wall, cpu = time.perf_counter(), time.process_time()
    startup = dict(wall_seconds=0.,cpu_seconds=0.,charged_once_per_worker=False,
                   scope='Solver module/dependency imports only; earlier remote bootstrap imports not measured by unchanged runner')
    if not _IMPORT_CHARGED:
        startup.update(wall_seconds=_IMPORT_WALL,cpu_seconds=_IMPORT_CPU,charged_once_per_worker=True)
        _IMPORT_CHARGED = True
    model = None
    receipt = dict(status='RUNNING', method=method, backend=backend, truth_access='NONE',
                   started_utc=datetime.now(timezone.utc).isoformat(), configuration=config,
                   source='research.research_max_struct_r4.src.confirmation_solver',
                   failed_attempts=[], scientific_wallcap=None, solver_import_charge=startup, source_cache='Per-task mature receiver SVD/DPCS and Fourier caches; charged in auxiliary counts',
                   operator_cache='Existing public operators loaded by remote worker; model adjoint/identity copies charged in setup',
                   historical_operator_assembly='Not rerun; generation charged separately', trajectory_arrays={})
    try:
        if backend != 'cpu' or method not in METHODS:
            raise ValueError('Only frozen CPU confirmation methods enabled')
        if not set(public) <= ALLOW:
            raise ValueError('Unexpected public payload key; refuse data')
        if config.get('sourcewise_cg') is not (method == 'native_extended'):
            raise ValueError('Frozen sourcewise/global policy mismatch')
        forbidden = {'max_wall_seconds','budget_check','max_iterations','receipt_path'} & set(config)
        if forbidden:
            raise ValueError('No scientific wallcap/hidden schedule or receipt-path override allowed')
        with threadpool_limits(limits=1):
            setup, setup_cpu = time.perf_counter(), time.process_time()
            model = MaxwellModel(public)
            if tuple(model.shape) != (24,24) or (model.n,model.ns,model.nr) != (576,12,48):
                raise ValueError('Frozen grid/acquisition mismatch')
            expected = dict(host=.1, amplitude=12., smooth_weight=1e-4, tv_weight=1e-4, tv_epsilon=.01, prior_width=19.5)
            if any(float(public[key]) != value for key,value in expected.items()):
                raise ValueError('Frozen host/prior/amplitude mismatch')
            if not np.array_equal(public['material_bounds'], np.array([.02,20.])):
                raise ValueError('Frozen material box mismatch')
            Y = np.asarray(public['Y'])
            m0 = np.full(model.n, config['host'])
            receipt['setup_wall_seconds'] = time.perf_counter()-setup
            receipt['setup_cpu_seconds'] = time.process_time()-setup_cpu
            receipt['setup_physical_model_cost'] = dict(model.cost)
            tick,ct,before = time.perf_counter(),time.process_time(),dict(model.cost)
            initial = metric(model,Y,m0)
            receipt['initial_wall_seconds'] = time.perf_counter()-tick
            receipt['initial_cpu_seconds'] = time.process_time()-ct
            receipt['initial_physical_model_cost'] = delta(model.cost,before)
            numerical = dict(config,method='cie',regularizer=lambda x:model.reg.evaluate(x)[:2], physical_sigma=float(np.linalg.norm(Y)))
            # No algorithm edits: use mature implementations directly.
            start, ct, before = time.perf_counter(),time.process_time(),dict(model.cost)
            if method == 'onehop':
                if config.get('cache_material_value_gradient') is not True:
                    raise ValueError('Frozen onehop material cache policy required')
                aux = run_onehop_cie(model,Y,m0,numerical)
            elif method == 'block4_cached':
                if not (config.get('gin_mode') == 'block4' and config.get('block_solver') == 'cached_inverse'
                        and config.get('cache_material_value_gradient') is True):
                    raise ValueError('Frozen cached block mode required')
                aux = run_block_cie(model,Y,m0,numerical)
            elif method == 'global':
                aux = global_cie(model,Y,m0,numerical)
            else:
                aux = native.run_cie(model,Y,m0,numerical)
            receipt['auxiliary_wall_seconds'] = time.perf_counter()-start
            receipt['auxiliary_cpu_seconds'] = time.process_time()-ct
            receipt['auxiliary_physical_model_cost'] = delta(model.cost,before)
            receipt['auxiliary_operation_counts'] = aux['cost']
            receipt['auxiliary_execution_receipt'] = aux.get('execution_receipt')
            receipt['auxiliary_status'] = aux['status']
            receipt['history_rows'] = len(aux['history'])
            receipt['stage_history'] = [row for row in aux['history'] if row.get('event') == 'stage_end']
            receipt['trajectory_arrays'].update(history_arrays(aux['history']))
            start,before = time.perf_counter(),dict(model.cost)
            scored = [metric(model,Y,c['m']) for c in aux['candidates']]
            seed = min([initial]+scored,key=lambda row:row['objective'])
            receipt['ranking_wall_seconds'] = time.perf_counter()-start
            receipt['ranking_physical_model_cost'] = delta(model.cost,before)
            receipt['candidate_scores'] = [{k:v for k,v in row.items() if k != 'm'} for row in scored]
            receipt['trajectory_arrays'].update(candidate_m=np.array([c['m'] for c in scored]),
                candidate_original_Phi=np.array([c['objective'] for c in scored]), selected_before_corrector_m=seed['m'])
            start,ct,before = time.perf_counter(),time.process_time(),dict(model.cost)
            refined, corrector = refine(model,Y,seed['m'],maximum=config['corrector_maxiter'])
            receipt['corrector_wall_seconds'] = time.perf_counter()-start
            receipt['corrector_cpu_seconds'] = time.process_time()-ct
            receipt['corrector_physical_model_cost'] = delta(model.cost,before)
            receipt['trajectory_arrays']['corrector_history'] = np.array(corrector.pop('history'))
            receipt['corrector'] = corrector
            tick,ct,before = time.perf_counter(),time.process_time(),dict(model.cost)
            best = min([initial,seed,metric(model,Y,refined)],key=lambda row:row['objective'])
            state = model.state(best['m'])
            receipt['finalization_wall_seconds'] = time.perf_counter()-tick
            receipt['finalization_cpu_seconds'] = time.process_time()-ct
            receipt['finalization_physical_model_cost'] = delta(model.cost,before)
            receipt['trajectory_arrays'].update(m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'])
            receipt.update(status='SEALED_PREDICTION',physical_objective=best['objective'],
                physical_relative_residual=best['relative_residual'],
                physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
                noise_radius99=float(public['noise_radius']),linear_relative_residual=float(state['state_relative_residual']),
                material_min=float(best['m'].min()),material_max=float(best['m'].max()),
                selection='Original Phi including unchanged prior; host + true stage endpoints + selected-seed refine',
                complete_auxiliary_schedule=aux['status']=='SCHEDULE_COMPLETE')
    except Exception as exc:
        # Return charged failure receipt to remote worker; caller must use INNER status.
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),
                       auxiliary_cost_completeness='INCOMPLETE_UNLESS_FAILED_EXECUTION_RECEIPT_AVAILABLE')
        receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
        if hasattr(exc,'cie_execution_receipt'):
            receipt['failed_auxiliary_execution'] = exc.cie_execution_receipt
    finally:
        receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.perf_counter()-wall,
            cpu_seconds=time.process_time()-cpu,physical_model_cost=dict(model.cost) if model else {},
            physical_operation_wall=dict(model.timings) if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings) if model else {}, memory=peak_memory(),
            cost_scope='solve validation/model setup/auxiliary/ranking/corrector/final state; worker separately charges public loads/NPZ serialization; solver imports charged once separately, earlier remote bootstrap imports unavailable',
            rhs_semantics='physical source_RHS/adjoint_RHS/actual_complex_RHS separate from small4x4 inverse RHS and matrix-apply RHS; do not sum as equivalent full solves')
    return receipt
