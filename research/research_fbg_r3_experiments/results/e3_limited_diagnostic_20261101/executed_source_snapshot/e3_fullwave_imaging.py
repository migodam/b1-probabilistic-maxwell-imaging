"""Finite diagnostic E3 queue; does not run without explicit --execute.

Manifest/data/GT decisions belong to the parent. This runner reads only public
G/E0/S/Y and predefined public initializers; all 576 material cells remain free.
"""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import types

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'


def _module(package, name, file):
    if package not in sys.modules:
        pkg = types.ModuleType(package)
        pkg.__path__ = [str(file.parent)]
        sys.modules[package] = pkg
    spec = importlib.util.spec_from_file_location(package+'.'+name, file)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _resolve(root, value):
    path = Path(value)
    return path if path.is_absolute() else root/path


def idle_snapshot():
    """Read-only resource snapshot; never changes other processes or priority."""
    commands = [
        ['nvidia-smi', '--query-gpu=utilization.gpu,memory.used,memory.total', '--format=csv,noheader'],
        ['nvidia-smi', '--query-compute-apps=pid,process_name,used_memory', '--format=csv,noheader'],
    ]
    if os.name == 'nt':
        commands.append(['powershell', '-NoProfile', '-Command',
            'Get-CimInstance Win32_Process | Where-Object { $_.Name -match "python" } | Select-Object ProcessId,Name,CommandLine | ConvertTo-Json -Depth 3'])
    rows = []
    for command in commands:
        try:
            reply = subprocess.run(command, capture_output=True, text=True, timeout=15)
            rows.append(dict(command=command, returncode=reply.returncode, stdout=reply.stdout, stderr=reply.stderr))
        except (OSError, subprocess.TimeoutExpired) as exc:
            rows.append(dict(command=command, error=str(exc)))
    return rows


def run(manifest_path, execute=False):
    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    if not execute:
        print(json.dumps(dict(status='NOT_RUN', reason='Explicit --execute required after parent authorization',
                              manifest=str(manifest_path), tasks=len(manifest.get('tasks', []))), indent=2))
        return
    if not manifest.get('resource_idle_confirmed', False):
        raise RuntimeError('Parent must inspect current resource snapshot and set resource_idle_confirmed for this run; never interrupt other jobs')
    import numpy as np
    affinity = int(manifest.get('affinity_mask', 21845))
    if os.name == 'nt':
        import ctypes
        kernel = ctypes.windll.kernel32
        kernel.GetCurrentProcess.restype = ctypes.c_void_p
        kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        if not kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), affinity):
            raise OSError('Could not set this process affinity')
    pipeline = _resolve(root, manifest['pipeline_root']).resolve()
    r2root = _resolve(root, manifest['r2_root']).resolve()
    study = Path(__file__).resolve().parents[1]
    os.environ['CUDA_PIPELINE_ROOT'] = str(pipeline)
    if manifest.get('vendor_path'):
        sys.path.insert(0, str(_resolve(root, manifest['vendor_path'])))
    from threadpoolctl import threadpool_limits
    r2 = _module('fbg_r3_readonly_r2', 'adaptive_solver', r2root/'src/adaptive_solver.py')
    port = _module('fbg_r3_current_port', 'cie_tsom_cuda', study/'src/cie_tsom_cuda.py')
    torch = r2.torch
    original_lane = r2.FeedbackLane
    class CountedLane(original_lane):
        """Reserve final closure, charge it, never bypass the common caps."""
        _closing = False
        def _budget(self):
            if self.task_context is None:
                return
            counts = Counter()
            for other in self.task_context['lanes']:
                counts.update(other.counts)
            cfg = self.task_context['config']
            rhs = sum(counts[k] for k in ('source_RHS', 'adjoint_RHS', 'receiver_RHS', 'tangent_RHS'))
            s = self.shared
            # One full new evaluation plus one final physical closure can fit.
            need_lu = 1 if self._closing else 2
            need_rhs = s.ns if self._closing else 2*s.ns+max(s.ns, s.nr)
            if counts['LU']+need_lu > cfg.get('max_lu', float('inf')):
                raise r2.engine.BudgetStop('LU_CAP_WITH_PAID_CLOSURE')
            if rhs+need_rhs > cfg.get('max_rhs', float('inf')):
                raise r2.engine.BudgetStop('RHS_CAP_WITH_PAID_CLOSURE')
        def close(self, x):
            self._closing = True
            try:
                return r2.base.ObjectiveLane.close(self, x)
            finally:
                self._closing = False
    r2.FeedbackLane = CountedLane
    torch.set_num_threads(1)
    torch.cuda.set_device(int(manifest.get('device', 0)))
    lanes = int(manifest.get('lanes', 2))
    if lanes not in (2, 4):
        raise ValueError('Finite queue supports two or four existing stream lanes')
    tasks = list(manifest['tasks'])
    output = _resolve(root, manifest['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    if (output/'SUMMARY.json').exists() or any((output/(t['task_id']+'.json')).exists() for t in tasks):
        raise FileExistsError('Choose a fresh output directory; preserve old receipts')
    snapshot = idle_snapshot()
    started = time.perf_counter()
    shared_cache, observation_cache, setup = {}, {}, []

    def current_config(task, Y):
        cfg = {**manifest.get('config', {}), **task.get('config', {})}
        cfg.setdefault('dpcs_rtol', 1e-3)
        cfg.setdefault('coefficient_tol', 1e-8)
        cfg.setdefault('miniter', 50)
        cfg.setdefault('material_inner_maxiter', 60)
        cfg.setdefault('stages', [dict(mf=3, beta=1., maxiter=500),
                                  dict(mf=6, beta=.2, maxiter=500),
                                  dict(mf=10, beta=.2, maxiter=500)])
        cfg['stages'] = [dict(s) for s in cfg['stages']]
        method = task['method']
        cfg['method'] = {'cie_single':'cie', 'cie_dual':'cie', 'tsom':'tsom', 'csi':'csi'}.get(method, method)
        if method == 'cie_dual':
            for stage in cfg['stages'][1:]:
                stage.setdefault('beta_current', 1.)
        return cfg

    current_methods = {'cie_single', 'cie_dual', 'tsom', 'csi'}
    with threadpool_limits(limits=1):
        for task in tasks:
            operator_path = str(_resolve(root, task.get('operators_path', manifest['operators_path'])).resolve())
            if operator_path not in shared_cache:
                tick = time.perf_counter()
                with np.load(operator_path, allow_pickle=False) as source:
                    # Explicit public allowlist; no private/truth/label arrays.
                    operators = {key:source[key] for key in ('G', 'E0', 'S')}
                shared = r2.base.SharedOperators(operators)
                # Exactly the second paid R2 shared reference initialization.
                with torch.no_grad():
                    ref_lu = torch.linalg.lu_factor_ex(shared.eye-.1*shared.G, check_errors=False)
                    shared.feedback_reference_j = torch.linalg.lu_solve(ref_lu.LU, ref_lu.pivots, .1*shared.E0)
                    shared.feedback_reference_q = shared.E0+shared.G@shared.feedback_reference_j
                    shared.feedback_reference_lu, shared.feedback_reference_pivots = ref_lu.LU, ref_lu.pivots
                    shared.setup_counts.update(LU=2, source_RHS=2*shared.ns)
                    shared.ready.record()
                shared.cie_host_operators = operators
                shared.cie_objective_lane_class = r2.base.ObjectiveLane
                shared.ready.synchronize()
                shared_cache[operator_path] = shared
                setup.append(dict(kind='existing_shared_operators', operators_path=operator_path,
                                  wall_seconds=time.perf_counter()-tick, counts=shared.setup_counts))
            public_path = str(_resolve(root, task['public_path']).resolve())
            if public_path not in observation_cache:
                with np.load(public_path, allow_pickle=False) as source:
                    observation_cache[public_path] = {'Y': source['Y']}
            if task['method'] in current_methods-{'csi'}:
                cfg = current_config(task, observation_cache[public_path]['Y'])
                prepared, new = port.prepare_cie_shared(shared_cache[operator_path], cfg)
                if new:
                    setup.append(dict(kind='shared_CIE_public_basis', operators_path=operator_path, counts=prepared['cost']))
        setup_wall = time.perf_counter()-started
        hybrid = None
        if manifest.get('svd_backend') == 'hybrid':
            hybrid = r2._module('fbg_r3_existing_hybrid_svd', pipeline/'src/hybrid_svd.py').HybridSVD()
            hybrid.install()

        def worker(task):
            torch.cuda.set_device(int(manifest.get('device', 0)))
            began = time.perf_counter()
            context = dict(lanes=[])
            try:
                operator_path = str(_resolve(root, task.get('operators_path', manifest['operators_path'])).resolve())
                public_path = str(_resolve(root, task['public_path']).resolve())
                shared, public = shared_cache[operator_path], observation_cache[public_path]
                cfg = {**manifest.get('config', {}), **task.get('config', {})}
                context['config'] = cfg
                envelope = float(cfg.get('wall_seconds', 45.))
                initial = np.asarray(task.get('initial_pixels', np.full(shared.n, task.get('initial_material', .1))), float)
                current = None
                if task['method'] in current_methods:
                    ccfg = current_config(task, public['Y'])
                    ccfg['max_wall_seconds'] = envelope*float(cfg.get('proposal_fraction', .55))
                    ccfg['_context'] = context
                    current = port.run_cie_cuda(shared, public['Y'], initial, ccfg)
                    initial = current['m']
                    remaining = max(0., envelope-(time.perf_counter()-began))
                    if remaining > 0:
                        correction_cfg = {key:cfg[key] for key in ('smooth_weight', 'tv_weight', 'tv_epsilon', 'max_lu', 'max_rhs', 'lower', 'upper') if key in cfg}
                        correction_cfg.update(cfg.get('corrector', {}))
                        correction_cfg.update(wall_seconds=remaining, stages=[24], stage_fractions=[1.])
                        correction_task = dict(method=cfg.get('corrector_method', 'data_lbfgs'),
                                               initial_pixels=initial.tolist(), config=correction_cfg)
                        context['config'] = correction_cfg
                        correction = r2.solve(shared, public, correction_task, context)
                    else:
                        correction = None
                    chosen = min(current['candidates'], key=lambda row:row['objective'])
                    result = dict(m=chosen['m'], ranking_objective=chosen['objective'],
                                  relative_residual=chosen['physical_relative_residual'], selected_phase='current',
                                  current=current, corrector=correction)
                    if correction is not None and correction['ranking_objective'] < result['ranking_objective']:
                        result.update(m=correction['m'], ranking_objective=correction['ranking_objective'],
                                      relative_residual=correction['relative_residual'], selected_phase='corrector')
                else:
                    context['config'] = cfg
                    baseline_task = {**task, 'config':cfg}
                    result = r2.solve(shared, public, baseline_task, context)
                costs = Counter()
                for lane in context['lanes']:
                    costs.update(lane.counts)
                costs['actual_complex_RHS'] = sum(costs[key] for key in ('source_RHS', 'adjoint_RHS', 'receiver_RHS', 'tangent_RHS'))
                row = dict(task=task, status='OK', result=result, physical_counts=dict(costs),
                           current_costs=dict(context.get('current_costs', {})), truth_access='NONE_PUBLIC_ALLOWLIST',
                           all_material_cells_unknown=shared.n, wall_seconds=time.perf_counter()-began,
                           budget_envelope_seconds=envelope)
            except Exception as exc:
                partial = Counter()
                for lane in context['lanes']:
                    partial.update(lane.counts)
                row = dict(task=task, status='FAIL', error=str(exc), traceback=traceback.format_exc(),
                           partial_physical_counts=dict(partial), current_costs=dict(context.get('current_costs', {})),
                           wall_seconds=time.perf_counter()-began, truth_access='NONE_PUBLIC_ALLOWLIST')
            return row

        receipts = []
        try:
            # At most `lanes` futures exist. A completion immediately frees one
            # slot; no batch barrier, separate CUDA process, or permanent service.
            with ThreadPoolExecutor(max_workers=lanes, thread_name_prefix='fbg-r3-e3') as pool:
                iterator = iter(tasks)
                pending = set()
                for _ in range(lanes):
                    task = next(iterator, None)
                    if task is not None:
                        pending.add(pool.submit(worker, task))
                while pending:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        row = future.result()
                        receipts.append(row)
                        r2.engine.write_json(output/(row['task']['task_id']+'.json'), row)
                        print(row['task']['task_id'], row['status'], round(row['wall_seconds'], 3), flush=True)
                        task = next(iterator, None)
                        if task is not None:
                            pending.add(pool.submit(worker, task))
        finally:
            if hybrid is not None:
                hybrid.restore()
        r2.engine.write_json(output/'SUMMARY.json', dict(
            utc=datetime.now(timezone.utc).isoformat(), scope='LIMITED_E3_MATURE_BASELINE_DIAGNOSTIC',
            tasks=len(receipts), ok=sum(row['status']=='OK' for row in receipts),
            wall_seconds=time.perf_counter()-started, setup_wall_seconds=setup_wall,
            shared_setup=setup, resource_snapshot=snapshot, manifest=manifest,
            shader_or_driver_changes='NONE', shared_source_changes='NONE',
            scientific_expansion_gate='NOT_PASSED_BY_THIS_DIAGNOSTIC'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--inspect-resources', action='store_true')
    args = parser.parse_args()
    if args.inspect_resources:
        print(json.dumps(idle_snapshot(), indent=2))
        return
    run(args.manifest, args.execute)


if __name__ == '__main__':
    main()
