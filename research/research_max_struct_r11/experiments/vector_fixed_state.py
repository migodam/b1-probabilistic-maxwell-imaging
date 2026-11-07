"""One frozen public3D state, all-source device readout cost; no physical solve."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import time,json,sys,argparse,traceback
from pathlib import Path
ENTRY=time.perf_counter()
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r11.src.vector_denominator import VectorLowRankDenominator
from research.research_max_struct_r10.src.vector_resources import CooperativeGate

def main():
    p=argparse.ArgumentParser();p.add_argument('--basis',required=True);p.add_argument('--state',required=True)
    p.add_argument('--output',required=True);p.add_argument('--device',default='cuda',choices=['cpu','cuda']);args=p.parse_args()
    out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False);torch.set_num_threads(1)
    gate=CooperativeGate(out/'GPU_BOUNDARIES.jsonl',memory_limit_bytes=2_000_000_000)
    if args.device=='cuda':gate.wait('before_fixed_state_context')
    cpu=time.process_time();failure=None;rows=[];status='FAILED'
    try:
        setup=time.perf_counter()
        with np.load(args.basis,allow_pickle=False) as z:
            U=torch.as_tensor(z['right_vectors'],device=args.device);GU=torch.as_tensor(z['images'],device=args.device)
        with np.load(args.state,allow_pickle=False) as z:
            J=torch.as_tensor(z['J'],device=args.device);Q=torch.as_tensor(z['Q'],device=args.device)
            m=torch.as_tensor(z['material'],device=args.device);omega=torch.as_tensor(z['weights'],device=args.device)
        if args.device=='cuda':torch.cuda.synchronize()
        setup_seconds=time.perf_counter()-setup;op=VectorLowRankDenominator(U,GU)
        for mode in ('full','lagged'):
            if args.device=='cuda':gate.wait('fixed_state_'+mode,torch.cuda.memory_reserved())
            start,cpustart=time.perf_counter(),time.process_time()
            result,info=op.readout(J,Q,m,1.,omega,mode,rho=1/60)
            if args.device=='cuda':torch.cuda.synchronize()
            cost_wall=time.perf_counter()-start;cost_cpu=time.process_time()-cpustart
            packed=info.pop('packed_metrics').cpu().tolist();trace=info.pop('device_step_trace').cpu().tolist()
            true_full=float(op.value_gradient(J,Q,result,1.,omega,1/60,gradient=False)[0].item())
            info.update(mode=mode,readout_wall_seconds=cost_wall,readout_CPU_seconds=cost_cpu,
                own_objective_before=packed[0],ordinary_seed_own_objective=packed[1],own_objective_after=packed[2],
                seed_selected=bool(packed[3]),max_small_LU_info=packed[4],
                true_full_objective_after=true_full,steps_trace=trace,
                trace_columns=['own_objective','alpha','projected_direction_slope','gradient_inf'],
                maximum_material=float(result.max().item()),minimum_material=float(result.min().item()))
            np.savez_compressed(out/(mode+'_material.npz'),material=result.cpu().numpy());rows.append(info)
        status='COMPLETE'
    except BaseException as exc:failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
    receipt=dict(status=status,device=args.device,rows=rows,failure=failure,GT_access=False,new_physical_RHS=0,
        shared_setup_wall_seconds=setup_seconds if 'setup_seconds' in locals() else None,
        CPU_seconds=time.process_time()-cpu,whole_entry_wall_seconds=time.perf_counter()-ENTRY,
        counts=dict(op.cost) if 'op' in locals() else {},
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_initialized() else 0,
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_initialized() else 0,
        gate_helper_CPU_seconds=gate.helper_CPU_seconds,gate_query_count=gate.query_count,
        wholeprocess_includes_imports=True,CPU_scope='Post-import Python process CPU; whole OS process separately measured by transport',
        basis=str(Path(args.basis)),state=str(Path(args.state)),source_shape=list(J.shape) if 'J' in locals() else None)
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n');print(json.dumps(receipt))
    return 0 if status=='COMPLETE' else 1

if __name__=='__main__':raise SystemExit(main())
