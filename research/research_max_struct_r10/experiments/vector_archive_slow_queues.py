"""Root-authorized own slow-attempt stop/archive; no GT/NPZ reads.
Retained handle gives all consumed CPU, including any interrupted trial. Saved
stage/40/80 endpoints are preserved; an unsaved paid tail is explicit, not zero.
"""
import time,json,os,ctypes
from ctypes import wintypes as w
from pathlib import Path
import psutil
ROOT=Path(__file__).resolve().parents[1];ARCH=ROOT/'RESULTS/vector_engineering_attempts';ARCH.mkdir(parents=True,exist_ok=True)
T=time.perf_counter();C=time.process_time();summaries=[]
k=ctypes.WinDLL('kernel32',use_last_error=True)
k.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];k.OpenProcess.restype=w.HANDLE
k.TerminateProcess.argtypes=[w.HANDLE,w.UINT];k.TerminateProcess.restype=w.BOOL
k.WaitForSingleObject.argtypes=[w.HANDLE,w.DWORD];k.WaitForSingleObject.restype=w.DWORD
k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
for case,queue,processfolder,tag in [('v04','windows_v04_queue','case_v04_windows_process','v04_Windows_CPU_dense_attempt1'),('v08','windows_v08_gpu_queue','case_v08_gpu_process','v08_CUDA_FFT_GMRES_attempt1')]:
    destination=ARCH/tag;destination.mkdir(exist_ok=False)
    stages=[]
    for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
        folder=ROOT/'RESULTS/vector'/(case+'_'+arm)
        if not folder.exists():continue
        trace=folder/'TRACE.jsonl';rows=[]
        if trace.exists():
            for line in trace.read_text().splitlines():
                try:rows.append(json.loads(line))
                except json.JSONDecodeError:pass # An interrupted final append is retained raw.
        requests=[r for r in rows if r.get('kind')=='physical_evaluation'];accepted=[r for r in rows if r.get('kind')=='physical_accepted']
        saved=next((label for label in ('physical_accepted80','physical_accepted40','stage3_physical','stage2_physical','stage1_physical','physical_initial') if (folder/(label+'.json')).exists() and (folder/(label+'.npz')).exists()),None)
        stages.append(dict(task_id=folder.name,SEALED=(folder/'SEALED.json').exists(),requests=requests[-1]['request'] if requests else 0,
            accepted=accepted[-1]['accepted'] if accepted else 0,last_serialized_endpoint=saved,
            unsaved_paid_tail_is_not_claimed_serialized=True))
    marker=ROOT/'delegated/vector'/queue/'ACTIVE.json';active=json.loads(marker.read_text());pid=active['pid'];cpu=None;wall=None;terminated=False
    try:
        process=psutil.Process(pid);cmd=' '.join(process.cmdline())
        assert 'MAX_STRUCT_R10_VECTOR_20261107' in cmd and 'vector_run.py' in cmd and processfolder in cmd
        handle=k.OpenProcess(0x1000|0x100000|1,False,pid)
        if not handle:raise ctypes.WinError(ctypes.get_last_error())
        created=process.create_time()
        if not k.TerminateProcess(handle,78):raise ctypes.WinError(ctypes.get_last_error())
        k.WaitForSingleObject(handle,30000);a,b,c,d=(w.FILETIME() for _ in range(4))
        if not k.GetProcessTimes(handle,ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
        cpu=((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7;wall=time.time()-created;terminated=True
    except psutil.NoSuchProcess:pass
    # Existing parent is allowed to finish its failed-child receipt; no new arm
    # can launch because its child returned78 and its fixed queue stops on error.
    parent=ROOT/'delegated/vector'/queue/'CONTROLLER_COST.json'
    for _ in range(100):
        if parent.exists():break
        time.sleep(.1)
    old_process=ROOT/'RESULTS/vector'/processfolder
    record=dict(status='ROOT_COST_DRIVEN_WHOLE_COHORT_REPLACEMENT_NOT_GT',platform='XINAN',pid=pid,terminated_own_process=terminated,
        CPU_seconds=cpu,wall_seconds=wall,CPU_source='GetProcessTimes retained handle after stop; same child as parent queue final row, do not add twice',
        trigger='Public uniform CPU/GPU retained-dense vs sourcewiseFFT measurement, no GT',
        checkpoints=stages,scientific_parameters_changed=False,new_primary_restarts_from_uniform_m0=True,
        budget_policy='Full original 500x3/min20+80 corrector for replacement; old partial spend retained separately',
        interrupted_trial_cost_included=True,inflight_unrecorded_operations_unknown=True,other_processes_modified=False)
    (destination/'STOP_AND_MIGRATION_RECEIPT.json').write_text(json.dumps(record,indent=2)+'\n')
    mapping=[]
    for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
        source=ROOT/'RESULTS/vector'/(case+'_'+arm)
        if source.exists():
            target=destination/source.name;mapping.append(dict(old=str(source.relative_to(ROOT)),new=str(target.relative_to(ROOT))))
            source.rename(target)
    if old_process.exists():old_process.rename(destination/processfolder)
    (destination/'ARCHIVE_PATH_MAPPING.json').write_text(json.dumps(mapping,indent=2)+'\n')
    summaries.append(record);print(json.dumps(record),flush=True)
(ROOT/'delegated/vector/SLOW_QUEUE_MIGRATION.json').write_text(json.dumps(dict(summaries=summaries,helper_CPU_seconds=time.process_time()-C,helper_wall_seconds=time.perf_counter()-T,platform='XINAN',GT_read=False),indent=2)+'\n')
