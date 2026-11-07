"""Root-authorized 3 exposed states ×2 existing CUDA LU libraries, no GT."""
import time,json,gc,traceback
from pathlib import Path
import numpy as np
import torch
from research.research_max_struct_r10.src.vector_physical import PhysicalObjective

def compare_public_linalg_libraries(model,state_path,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    with np.load(state_path,allow_pickle=False) as z:materials=z['material'];Y=z['Y']
    if materials.shape!=(3,model.N) or Y.shape!=(180,60) or np.any(materials<0) or np.any(materials>8):raise ValueError('Only the3 predeclared legal public states permitted')
    rows=[];totals=[];start,cpustart=time.perf_counter(),time.process_time();before=dict(model.counts)
    for library in ('default','cusolver'):
        setting_wall,setting_cpu=time.perf_counter(),time.process_time();setting_failure=None
        try:torch.backends.cuda.preferred_linalg_library(library)
        except BaseException as exc:setting_failure=dict(type=type(exc).__name__,message=str(exc))
        setup=dict(CPU_seconds=time.process_time()-setting_cpu,wall_seconds=time.perf_counter()-setting_wall)
        qualified=True;library_rows=[]
        for index,m in enumerate(materials):
            if setting_failure:
                row=dict(library=library,state=index,status='LIBRARY_UNAVAILABLE_NOT_RUN',failure=setting_failure)
                qualified=False;rows.append(row);library_rows.append(row);continue
            model.gate.wait('public_linalg_comparison_'+library+'_state'+str(index),torch.cuda.memory_reserved())
            model.dense_state=None;objective=PhysicalObjective(model,Y,upper=8.)
            w,cpu=time.perf_counter(),time.process_time();counts0=dict(model.counts);records0=len(model.solve_records);batch0=len(model.batch_records)
            try:
                result=objective.evaluate(m,library+'_public_state'+str(index),True)
                torch.cuda.synchronize()
                residual=max(r['relative_residual'] for r in model.solve_records[records0:])
                ok=bool(np.isfinite(result['Phi']) and np.isfinite(result['gradient']).all() and residual<=2e-6)
                row=dict(library=library,state=index,status='QUALIFIED' if ok else 'FAILED_QUALIFICATION',
                    Phi=result['Phi'],gradient_inf=float(np.max(abs(result['gradient']))),maximum_original_relative_residual=residual,
                    CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-w,
                    counts_delta={k:int(v-counts0.get(k,0)) for k,v in model.counts.items() if v-counts0.get(k,0)},
                    batch_records=model.batch_records[batch0:],GPU_peak_allocated_bytes_to_date=torch.cuda.max_memory_allocated(),GPU_peak_reserved_bytes_to_date=torch.cuda.max_memory_reserved())
                qualified=qualified and ok
            except BaseException as exc:
                row=dict(library=library,state=index,status='FAILED',failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()),
                    CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-w,
                    counts_delta={k:int(v-counts0.get(k,0)) for k,v in model.counts.items() if v-counts0.get(k,0)})
                qualified=False
            rows.append(row);library_rows.append(row);del objective;model.dense_state=None;gc.collect()
            (output/'PROGRESS.json').write_text(json.dumps(dict(rows=rows,truth_access=False),indent=2)+'\n')
        totals.append(dict(library=library,all3_qualified=qualified,setting_cost=setup,
            whole_evaluation_wall_seconds=sum(r.get('wall_seconds',0.) for r in library_rows)+setup['wall_seconds'],
            whole_evaluation_CPU_seconds=sum(r.get('CPU_seconds',0.) for r in library_rows)+setup['CPU_seconds']))
    eligible=[t for t in totals if t['all3_qualified']]
    if not eligible:raise RuntimeError('Neither declared library qualified on all3 public states; no confirmation started')
    selected=min(eligible,key=lambda r:(r['whole_evaluation_wall_seconds'],r['whole_evaluation_CPU_seconds']))['library']
    torch.backends.cuda.preferred_linalg_library(selected);model.dense_state=None
    receipt=dict(status='EXISTING_LIBRARY_SELECTED_BEFORE_ANY_CONFIRMATION_INVERSE',selected_library=selected,rows=rows,totals=totals,
        selection='Minimum summed3 whole value+gradient+original-residual wall among libraries qualifying all3, with CPU tie-break; no singleLU speed extrapolation',
        CPU_seconds=time.process_time()-cpustart,wall_seconds=time.perf_counter()-start,
        counts_delta={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)},
        shared_model_setup_outside_comparison=True,all60_primal_and_adjoint=True,precision='complex128/float64 unchanged',truth_access=False,
        public_states=str(state_path),new_confirmation_Y_access=False,new_confirmation_inverse_calls=0)
    (output/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(event='CONFIRMATION_LINALG_SELECTED',library=selected,totals=totals,counts=receipt['counts_delta'])),flush=True)
    return receipt
