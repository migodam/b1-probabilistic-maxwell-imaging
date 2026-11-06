"""One finite resident GPU process, three exposed-development physical starts.
No truth loads or score calls. Additive process receipt includes all paid work.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,datetime,traceback,importlib.util
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
OUT=ROOT/'RESULTS/original_vector3d/physical_refinement';CONFIG=ROOT/'configs/original_physical_refinement.json'
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,r):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)

def main():
    import psutil
    process=psutil.Process();process.cpu_affinity([12,14])
    spec=importlib.util.spec_from_file_location('joblib',ROOT/'delegated/3d_extension/windows_resolved.py');joblib=importlib.util.module_from_spec(spec);spec.loader.exec_module(joblib);job=joblib.job_object()
    import numpy as np
    import torch
    from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel
    from research.research_max_struct_r9.src.vector_original_cuda import CudaOriginalModel
    from research.research_max_struct_r9.src.vector_physical_refinement import PhysicalObjective,run_short
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    cfg=json.loads(CONFIG.read_text());OUT.mkdir(parents=True,exist_ok=False);dump(OUT/'FROZEN_CONFIG.json',cfg)
    model=None;objective=None;status='FAILED';failure=None;terminal=[]
    try:
        with np.load(ROOT/cfg['data'],allow_pickle=False) as z:
            if any(k.lower() in ('truth','mask','chi','m_true','y_clean') for k in z.files):raise ValueError('Unexpected private payload in publicY')
            Y=z['Y'].copy()
        inputs=[stamp(CONFIG),stamp(ROOT/cfg['data'])];starts=[]
        for start in cfg['starts']:
            if start['material']=='analytic_zero':m=np.zeros(30**3)
            else:
                path=ROOT/start['material'];inputs.append(stamp(path))
                with np.load(path,allow_pickle=False) as z:m=z[start['key']].ravel().astype(float)
            starts.append((start['name'],m))
        dump(OUT/'PUBLIC_INPUT_EVIDENCE.json',dict(inputs=inputs,classification=cfg['classification'],GT_access=False))
        print(json.dumps(dict(event='PUBLIC_INPUTS_READY',GT_access=False,affinity=process.cpu_affinity())),flush=True)
        model=CudaOriginalModel(OriginalVectorModel(30,gmres_rtol=2e-6,gmres_restart=60,gmres_cycles=20),fft_layout='packed',block_rhs=16)
        print(json.dumps(dict(event='SHARED_MODEL_READY',setup=model.setup)),flush=True)
        def progress(r):
            with (OUT/'PHYSICAL_RHS.jsonl').open('a') as file:file.write(json.dumps(r)+'\n')
            if r.get('source',0)%20==0:print(json.dumps(dict(event='RHS',**r)),flush=True)
        objective=PhysicalObjective(model,Y,progress)
        wiring=objective.wiring(cfg['air_wiring']);dump(OUT/'GRADIENT_WIRING.json',wiring);print(json.dumps(dict(event='AIR_GRADIENT_WIRING',**wiring)),flush=True)
        if wiring['status']!='CONSISTENT_FIXED_AIR_WIRING':raise RuntimeError('Gradient wiring requiresrootreview; no optimization launched')
        dump(OUT/'SHARED_SETUP.json',dict(setup=model.setup,counts=dict(model.counts),additive=False,included_in='RUN_PROCESS_COST.json',receiverSVD=0,prior=0))
        for name,m in starts:
            dest=OUT/name;dest.mkdir(exist_ok=False)
            def save(label,r):
                np.savez_compressed(dest/(label+'.npz'),material=r['m'],prediction=r['prediction'],gradient=r['gradient'])
                dump(dest/(label+'.json'),dict(state_id=r['id'],Phi=r['Phi'],PG_inf=r['projected_gradient_inf'],first_owner=r['first_owner'],uses=r['uses'].copy()))
            def emit(row):
                with (dest/'TRACE.jsonl').open('a') as file:file.write(json.dumps(row)+'\n')
                print(json.dumps(row),flush=True)
            receipt=run_short(objective,m,name,cfg['optimizer'],save,emit);dump(dest/'RECEIPT.json',receipt);terminal.append(dict(start=name,status=receipt['status']))
            dump(OUT/'PAID_EVALUATIONS.json',objective.events)
            dump(OUT/'CACHE_OWNERSHIP.json',dict(additive=False,included_in='RUN_PROCESS_COST.json',
                records=[dict(state_id=r['id'],first_owner=r['first_owner'],uses=r['uses'],Phi=r['Phi'],gradient_available=r.get('gradient') is not None) for r in objective.records.values()]))
        entries=[]
        for name,_ in starts:entries.append(dict(start=name,files=[stamp(p) for p in sorted((OUT/name).iterdir()) if p.is_file()]))
        dump(OUT/'SEALED_ALL_THREE.json',dict(status='ALL_THREE_TERMINAL_PREDICTIONS_SEALED_BEFORE_OFFLINE_GT_SCORE',created_utc=utc(),entries=entries,
            input_evidence=inputs,source_evidence=[stamp(ROOT/'src/vector_physical_refinement.py'),stamp(Path(__file__))],classification=cfg['classification'],GT_access=False))
        status='COMPLETE'
    except BaseException as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
    finally:
        if model:torch.cuda.synchronize()
        cp=process.cpu_times();r=dict(status=status,platform='XINAN',CPU_seconds=cp.user+cp.system,wall_seconds=time.perf_counter()-ENTRY,
            additive=True,included_in=None,scope='Wholeprocess imports/sharedCM G/E/S/airFDwiring/all3optimizationarms/rejectedtrials/cache/I/O; no GT',
            operation='SHORT_ORIGINAL_PHYSICAL_REFINEMENT_THREE_STARTS',counts=dict(model.counts) if model else{},
            CPU_operator_counts=dict(model.cpu.counts) if model else{},terminals=terminal,failure=failure,actual_source_solve_records=model.solve_records if model else[],
            peak_RSS_bytes=process.memory_info().peak_wset,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated() if model else None,
            peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved() if model else None,CPU_affinity=process.cpu_affinity(),job_object=job,
            GT_access=False,old_nine_modified=False)
        dump(OUT/'RUN_PROCESS_COST.json',r)
        if objective:dump(OUT/'PAID_EVALUATIONS.json',objective.events)
    return 0 if status=='COMPLETE' else 1
if __name__=='__main__':raise SystemExit(main())
