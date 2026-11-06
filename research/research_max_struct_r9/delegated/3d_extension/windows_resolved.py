"""Portable public-only FIFO for the five resolved-hollow controls.
Uses unchanged StrictVectorInverse; no resource/threadpoolctl/GPU/truth imports.
Whole-process Windows CPU is psutil user+system (GetProcessTimes), peak_wset.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
ENTRY_WALL=time.perf_counter()
from pathlib import Path
import sys,json,argparse,subprocess,traceback,datetime,platform
ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
CONFIG=ROOT/'configs/vector_extension_resolved_hollow.json';PUBLIC=ROOT/'data_public/vector_extension';OUT=ROOT/'RESULTS/vector_extension';CASE='hollow_cube_resolved_control';LOGS=ROOT/'delegated/3d_extension/xinan_resolved_queue'
JOB_HANDLE=None

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def load(p):return json.loads(Path(p).read_text())
def dump(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def cpu_total(process):
    t=process.cpu_times();return t.user+t.system
def job_object():
    global JOB_HANDLE
    import ctypes
    from ctypes import wintypes as w
    class IO(ctypes.Structure):_fields_=[(x,ctypes.c_ulonglong) for x in ('ReadOperationCount','WriteOperationCount','OtherOperationCount','ReadTransferCount','WriteTransferCount','OtherTransferCount')]
    class BASIC(ctypes.Structure):
        _fields_=[('PerProcessUserTimeLimit',ctypes.c_longlong),('PerJobUserTimeLimit',ctypes.c_longlong),('LimitFlags',w.DWORD),('MinimumWorkingSetSize',ctypes.c_size_t),('MaximumWorkingSetSize',ctypes.c_size_t),('ActiveProcessLimit',w.DWORD),('Affinity',ctypes.c_size_t),('PriorityClass',w.DWORD),('SchedulingClass',w.DWORD)]
    class EXT(ctypes.Structure):
        _fields_=[('BasicLimitInformation',BASIC),('IoInfo',IO),('ProcessMemoryLimit',ctypes.c_size_t),('JobMemoryLimit',ctypes.c_size_t),('PeakProcessMemoryUsed',ctypes.c_size_t),('PeakJobMemoryUsed',ctypes.c_size_t)]
    k=ctypes.WinDLL('kernel32',use_last_error=True);k.CreateJobObjectW.argtypes=[w.LPVOID,w.LPCWSTR];k.CreateJobObjectW.restype=w.HANDLE
    k.SetInformationJobObject.argtypes=[w.HANDLE,ctypes.c_int,w.LPVOID,w.DWORD];k.SetInformationJobObject.restype=w.BOOL
    k.AssignProcessToJobObject.argtypes=[w.HANDLE,w.HANDLE];k.AssignProcessToJobObject.restype=w.BOOL;k.GetCurrentProcess.restype=w.HANDLE
    handle=k.CreateJobObjectW(None,None)
    if not handle:raise ctypes.WinError(ctypes.get_last_error())
    info=EXT();info.BasicLimitInformation.LimitFlags=0x2000
    if not k.SetInformationJobObject(handle,9,ctypes.byref(info),ctypes.sizeof(info)):raise ctypes.WinError(ctypes.get_last_error())
    if not k.AssignProcessToJobObject(handle,k.GetCurrentProcess()):raise ctypes.WinError(ctypes.get_last_error())
    JOB_HANDLE=handle;return dict(active=True,kill_on_job_close=True,parent_assigned=True,children_inherit=True,foreground_SSH=True)

def worker(arm,cpu):
    import psutil
    process=psutil.Process();process.cpu_affinity([cpu])
    if process.cpu_affinity()!=[cpu]:raise RuntimeError('Failed exact P-core affinity')
    # Affinity and environment are established before NumPy/SciPy import.
    import numpy as np
    from research.research_max_struct_r9.src.vector_extension import StrictVectorInverse
    cfg=load(CONFIG);folder=OUT/CASE/arm;folder.mkdir(parents=True,exist_ok=False)
    r=dict(status='RUNNING',case_id=CASE,arm=arm,started_utc=utc(),truth_access='NONE',scoring='NOT_RUN',
        config=stamp(CONFIG),input=[stamp(PUBLIC/CASE/'PUBLIC.json'),stamp(PUBLIC/CASE/'observations.npz')],
        implementation='Unchanged R9 StrictVectorInverse over inherited R1/R8 dyadicDDA; portableWindowswrapper only',
        platform='XINAN_Windows_i9',platform_details=platform.platform(),numpy=np.__version__,logical_cpu_affinity=[cpu],threads_per_process=1,
        prior_phase_files_modified=False,whole_process_CPU_source='psutil.cpu_times user+system, Windows GetProcessTimes',peak_RSS_source='psutil.memory_info.peak_wset')
    dump(folder/'RECEIPT.json',r);dump(folder/'FROZEN_CONFIG.json',cfg)
    try:
        meta=load(PUBLIC/CASE/'PUBLIC.json')
        with np.load(PUBLIC/CASE/'observations.npz',allow_pickle=False) as z:
            if z.files!=['Y']:raise ValueError('Unexpected public payload')
            Y=z['Y'].copy()
        model=StrictVectorInverse(8,meta['acquisition']['domain_width'],cfg,auxiliary=arm!='direct')
        prefix_cpu=cpu_total(process);prefix_wall=time.perf_counter()-ENTRY_WALL
        arrays,run=model.run(Y,meta,arm,lambda x:print(json.dumps(dict(case_id=CASE,arm=arm,**x)),flush=True))
        def save(path,state):
            values={k:state[k] for k in ('m','J','E','gradient') if k in state};values['Y_prediction']=state['Y'];np.savez_compressed(path,**values)
        for key,filename in [('initial','initial.npz'),('raw','raw.npz'),('final','prediction.npz')]:save(folder/filename,arrays[key])
        for name,state in arrays['snapshots'].items():save(folder/(name+'.npz'),state)
        for st in arrays['stages']:
            np.savez_compressed(folder/f'stage_{st["stage"]:02d}_joint.npz',**{k:st[k] for k in ('m','J','Q','gamma')})
            dump(folder/f'stage_{st["stage"]:02d}_joint.json',{k:v for k,v in st.items() if k not in ('m','J','Q','gamma')})
        dump(folder/'PAID_BUDGET.json',arrays['budget']);np.savez_compressed(folder/'PAID_BUDGET_MATERIAL.npz',m=arrays['budget_material'])
        dump(folder/'JOINT_HISTORY.json',run.pop('joint_history'));r.update(run,operator_setup=model.setup,prefix_CPU_seconds=prefix_cpu,prefix_wall_seconds=prefix_wall,
            cost_scope='Wholeprocess CPU includes interpreterstartup/imports/setup/online/export; online/prefix nested and not added twice; no borrowedMacwarmcost')
    except Exception as exc:r.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),counts=dict(model.counts) if 'model' in locals() else {})
    memory=process.memory_info()
    if not hasattr(memory,'peak_wset'):raise RuntimeError('Windows peak working set not exposed; do not write unknown zero')
    r.update(finished_utc=utc(),phase2_incremental_CPU_seconds=cpu_total(process),phase2_incremental_wall_seconds=time.perf_counter()-ENTRY_WALL,process_peak_RSS_bytes=memory.peak_wset)
    dump(folder/'RECEIPT.json',r);print(json.dumps(dict(event='WORKER_END',arm=arm,status=r['status'],error=r.get('error'))),flush=True)
    return 0 if r['status']=='SEALED_PUBLIC_PREDICTION' else 1

def queue():
    import psutil
    process=psutil.Process();process.cpu_affinity([8,10]);job=job_object();LOGS.mkdir(parents=True,exist_ok=False)
    cfg=load(CONFIG);pending=list(cfg['arms']);active={};free=[8,10];jobs=[];started=utc()
    while pending or active:
        while pending and free:
            arm=pending.pop(0);cpu=free.pop(0);path=LOGS/(arm+'.log');fh=path.open('x');command=[sys.executable,str(Path(__file__).resolve()),'worker','--arm',arm,'--cpu',str(cpu)]
            child=subprocess.Popen(command,cwd=REPO,stdout=fh,stderr=subprocess.STDOUT);row=dict(case_id=CASE,arm=arm,pid=child.pid,logical_cpu=cpu,command=command,started_utc=utc(),log=str(path.relative_to(REPO)))
            active[child]=(fh,time.perf_counter(),row);print(json.dumps(dict(event='START',**row)),flush=True)
        for child,(fh,t,row) in list(active.items()):
            code=child.poll()
            if code is None:continue
            fh.close();row.update(exit_code=code,ended_utc=utc(),wall_seconds=time.perf_counter()-t)
            path=OUT/CASE/row['arm']/'RECEIPT.json';r=load(path) if path.exists() else {};row.update(status=r.get('status','MISSING_RECEIPT'),error=r.get('error'))
            jobs.append(row);free.append(row['logical_cpu']);del active[child];print(json.dumps(dict(event='END',**row)),flush=True)
        dump(LOGS/'RECEIPT.json',dict(status='RUNNING' if active or pending else 'COMPLETE',started_utc=started,jobs=jobs,pending=len(pending),running=len(active),max_processes=2,
            platform='XINAN_Windows_i9',allocation_logical_Pcores=[8,10],threads_per_process=1,job_object=job,truth_access='NONE',scoring='NOT_RUN'))
        if active:time.sleep(.5)
    dump(LOGS/'PARENT_COST.json',dict(parent_CPU_seconds=cpu_total(process),parent_wall_seconds=time.perf_counter()-ENTRY_WALL,
        scope='Wholeparent process including imports/queue, excludeschildren separatelyreported; public inputs copied, no repeatgeneration',peak_RSS_bytes=process.memory_info().peak_wset))
    return 0 if all(j['exit_code']==0 for j in jobs) else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['queue','worker']);p.add_argument('--arm');p.add_argument('--cpu',type=int);args=p.parse_args()
    raise SystemExit(worker(args.arm,args.cpu) if args.command=='worker' else queue())
