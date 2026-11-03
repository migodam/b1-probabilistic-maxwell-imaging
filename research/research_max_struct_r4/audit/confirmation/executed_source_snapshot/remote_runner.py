"""Isolated Windows P-core CPU/CUDA worker pools; no installation or global affinity changes."""
from __future__ import annotations
import argparse, ctypes, datetime, importlib, json, multiprocessing as mp
import os, pathlib, queue, subprocess, time, traceback
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent/"vendor"))

for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

def cpu_topology():
    if os.name != 'nt':
        raise RuntimeError('P-core detection requires Windows GetSystemCpuSetInformation')
    k = ctypes.WinDLL('kernel32', use_last_error=True)
    size = ctypes.c_ulong()
    k.GetSystemCpuSetInformation(None, 0, ctypes.byref(size), None, 0)
    data = ctypes.create_string_buffer(size.value)
    if not k.GetSystemCpuSetInformation(data, size, ctypes.byref(size), None, 0):
        raise ctypes.WinError(ctypes.get_last_error())
    raw = data.raw; offset = 0; groups = {}
    while offset < size.value:
        record_size = int.from_bytes(raw[offset:offset+4], 'little')
        if record_size < 8: raise RuntimeError('Invalid CPU set record')
        kind = int.from_bytes(raw[offset+4:offset+8], 'little')
        if kind == 0:
            group = int.from_bytes(raw[offset+12:offset+14], 'little')
            logical, core, efficiency = raw[offset+14], raw[offset+15], raw[offset+18]
            entry = groups.setdefault((group, core), {'group':group,'core':core,'logical':[],'efficiency_class':efficiency})
            entry['logical'].append(logical)
        offset += record_size
    rows = sorted(groups.values(), key=lambda x:(x['group'],x['core']))
    return {'cores':rows, 'p_cores':[r for r in rows if len(r['logical'])>1 and r['group']==0],
            'criterion':'group 0 physical cores with >1 logical siblings (observed i9 HT P cores)'}

def bind_core(core):
    k = ctypes.WinDLL('kernel32', use_last_error=True)
    k.GetCurrentProcess.restype = ctypes.c_void_p
    k.SetProcessAffinityMask.argtypes = [ctypes.c_void_p,ctypes.c_size_t]
    # One logical sibling per physical P core: one runnable thread per worker.
    mask = 1 << min(core['logical'])
    if not k.SetProcessAffinityMask(k.GetCurrentProcess(), mask):
        raise ctypes.WinError(ctypes.get_last_error())
    return mask

def probe():
    import numpy, scipy, torch
    result = {'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'python':__import__('sys').version, 'executable':__import__('sys').executable,
              'numpy':numpy.__version__, 'scipy':scipy.__version__, 'torch':torch.__version__,
              'torch_cuda_build':torch.version.cuda, 'cuda_available':torch.cuda.is_available(),
              'topology':cpu_topology()}
    if result['cuda_available']:
        result['cuda_devices']=[{'index':i,'name':torch.cuda.get_device_name(i),'memory_bytes':torch.cuda.get_device_properties(i).total_memory} for i in range(torch.cuda.device_count())]
    commands = {'cpu_model':['powershell','-NoProfile','-Command','Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json -Compress'],
                'gpu':['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu','--format=csv,noheader'],
                'compute_processes':['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],
                'python_processes':['powershell','-NoProfile','-Command',"Get-CimInstance Win32_Process | Where-Object {$_.Name -match 'python|pythonw|jupyter'} | Select-Object ProcessId,Name,CommandLine | ConvertTo-Json -Compress"]}
    result['activity']={}
    for name, cmd in commands.items():
        p = subprocess.run(cmd,capture_output=True,text=True,timeout=20)
        result['activity'][name]={'returncode':p.returncode,'stdout':p.stdout.strip(),'stderr':p.stderr.strip()}
    return result

def serial(value):
    if hasattr(value,'tolist'): return value.tolist()
    if isinstance(value,pathlib.Path): return str(value)
    raise TypeError(type(value).__name__)

def write_json(path, data):
    path=pathlib.Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,default=serial,ensure_ascii=False,indent=2),encoding='utf-8')
    temporary.replace(path)

def worker(slot, core, backend, module_name, incoming, outgoing):
    try:
        mask=bind_core(core)
        import numpy as np
        if backend=='cuda':
            import torch
            torch.set_num_threads(1); torch.set_num_interop_threads(1)
            torch.cuda.set_device(0)
        solve=importlib.import_module(module_name).solve
        outgoing.put({'event':'ready','slot':slot,'pid':os.getpid(),'mask':mask})
        while True:
            task=incoming.get()
            if task is None: return
            ident=task['task_id']; started=time.monotonic()
            try:
                with np.load(task['operators_path'],allow_pickle=False) as f:
                    public={k:f[k] for k in f.files}
                with np.load(task['public_path'],allow_pickle=False) as f:
                    public.update({k:f[k] for k in f.files})
                if backend=='cuda': torch.cuda.synchronize()
                outgoing.put({'event':'started','slot':slot,'pid':os.getpid(),'task_id':ident})
                begin=time.perf_counter()
                result=solve(public,task['method'],task['config'],backend=backend)
                if backend=='cuda': torch.cuda.synchronize()
                if 'trajectory_arrays' in result:
                    path=pathlib.Path(task['trajectory_path'])
                    path.parent.mkdir(parents=True,exist_ok=True)
                    np.savez_compressed(path,**result.pop('trajectory_arrays'))
                    result['trajectory_path']=str(path)
                outgoing.put({'event':'result','slot':slot,'pid':os.getpid(),'task_id':ident,'status':'OK',
                              'result':result,'solve_wall_s':time.perf_counter()-begin,'total_worker_s':time.monotonic()-started})
            except BaseException as exc:
                outgoing.put({'event':'result','slot':slot,'pid':os.getpid(),'task_id':ident,'status':'FAIL',
                              'error':str(exc),'traceback':traceback.format_exc(),'total_worker_s':time.monotonic()-started,
                              'costs_status':'INCOMPLETE_EXCEPTION','costs':None})
    except BaseException as exc:
        outgoing.put({'event':'worker_failed','slot':slot,'pid':os.getpid(),'error':str(exc),'traceback':traceback.format_exc()})

def run(manifest_path):
    manifest_path=pathlib.Path(manifest_path).resolve()
    manifest=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    if manifest.get('execution_authorized') is not True:
        raise RuntimeError('Frozen task manifest must explicitly set execution_authorized=true')
    topology=cpu_topology(); cores=topology['p_cores']
    if len(cores)<6: raise RuntimeError('Need six independent physical P cores')
    output=pathlib.Path(manifest['output_dir']); output.mkdir(parents=True,exist_ok=True)
    tasks=manifest['tasks']; ids=[t['task_id'] for t in tasks]
    if len(ids)!=len(set(ids)): raise ValueError('Duplicate task_id')
    pending={b:[dict(t) for t in tasks if t['backend']==b] for b in ('cpu','cuda')}
    for t in tasks:
        if t['backend'] not in pending: raise ValueError('Unknown backend')
        if pathlib.Path(t['task_id']).name != t['task_id']: raise ValueError('task_id must be a filename')
    ctx=mp.get_context('spawn'); messages=ctx.Queue(); slots={}; receipts=[]; affinity_records=[]
    module=manifest.get('solver_module','src.solver')
    def launch(i, backend):
        q=ctx.Queue(); p=ctx.Process(target=worker,args=(i,cores[i],backend,module,q,messages))
        p.start(); slots[i]={'process':p,'queue':q,'backend':backend,'ready':False,'task':None,'phase_time':time.monotonic(),'solve_started':None}
    for i in range(6):
        backend='cpu' if i<4 else 'cuda'
        if pending[backend]: launch(i,backend)
    def finish(i, receipt):
        state=slots[i]; task=state['task']; receipt.update({'task':task,'backend':state['backend'],'core':cores[i]})
        write_json(output/(task['task_id']+'.json'),receipt); receipts.append(receipt)
        state['task']=None; state['solve_started']=None
    try:
        while slots:
            now=time.monotonic()
            for i,state in list(slots.items()):
                if state['ready'] and state['task'] is None:
                    if pending[state['backend']]:
                        task=pending[state['backend']].pop(0); state['task']=task; state['phase_time']=now; state['queue'].put(task)
                    else:
                        state['queue'].put(None); state['process'].join(10)
                        if state['process'].is_alive(): state['process'].terminate(); state['process'].join()
                        del slots[i]; continue
                timeout=state['solve_started'] is not None and now-state['solve_started']>float(state['task'].get('hard_timeout_s',47))
                stalled=state['solve_started'] is None and now-state['phase_time']>120
                if timeout or stalled or not state['process'].is_alive():
                    state['process'].terminate(); state['process'].join()
                    if state['task'] is not None:
                        finish(i,{'status':'HARD_TIMEOUT' if timeout else 'WORKER_FAILED','wall_observed_s':now-state['phase_time'],
                                  'costs_status':'INCOMPLETE_WORKER_TERMINATED','costs':None})
                    elif not state['ready']:
                        raise RuntimeError(f'Worker {i} failed during initialization')
                    backend=state['backend']; del slots[i]
                    if pending[backend]: launch(i,backend)
            if not slots: break
            try: event=messages.get(timeout=.2)
            except queue.Empty: continue
            i=event['slot']
            if i not in slots: continue
            state=slots[i]
            if event['pid'] != state['process'].pid: continue
            if event['event']=='ready':
                state['ready']=True; affinity_records.append(event)
                write_json(output/'worker_affinity.json',affinity_records)
            elif event['event']=='started': state['solve_started']=time.monotonic()
            elif event['event']=='result': finish(i,event)
            elif event['event']=='worker_failed':
                raise RuntimeError(json.dumps(event))
    finally:
        for state in slots.values():
            if state['process'].is_alive(): state['process'].terminate()
            state['process'].join()
        write_json(output/'runner_summary.json',{'expected':len(tasks),'completed':len(receipts),
                   'status_counts':{s:sum(r['status']==s for r in receipts) for s in set(r['status'] for r in receipts)},
                   'cpu_workers':4,'cuda_workers':2,'cores':cores[:6],'topology':topology,
                   'worker_affinity':affinity_records})

if __name__=='__main__':
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='command',required=True)
    p=sub.add_parser('probe'); p.add_argument('--output',default='hardware/environment.json')
    p=sub.add_parser('run'); p.add_argument('manifest')
    args=ap.parse_args()
    if args.command=='probe':
        result=probe(); write_json(args.output,result); print(json.dumps(result,indent=2))
    else: run(args.manifest)
