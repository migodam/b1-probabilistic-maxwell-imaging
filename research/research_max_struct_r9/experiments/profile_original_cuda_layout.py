"""One fixed-public-state full current/material replay: strided4 vs packed16.

Run only while the live inverse worker is safely paused. No full solves, GT,
parameter sweep, precision change, or scientific budget change.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import json,sys,time,platform
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
import psutil,torch
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel
from research.research_max_struct_r9.src.vector_original_cuda import CudaOriginalModel,CudaReceiverSubspace,CudaImplicitFourier,cuda_material_readout
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d'

def replay(model,receiver,Y):
    B=CudaImplicitFourier(model,receiver,6);J=receiver.Jd;Q=receiver.Qd;D=receiver.Dd
    m=torch.full((model.N,),.1,device='cuda',dtype=torch.float64);chi=torch.repeat_interleave(m,3)[:,None]
    wd=1/torch.sum(abs(Y)**2,dim=0);ws=1/torch.sum(abs(model.inc)**2,dim=0);p=6/(1+6*chi)
    v=p.conj()*p*(J-chi*Q)*ws
    gradient=B.adjoint(model.S.conj().T@(D*wd)+v-model.green(chi*v,True));direction=-gradient
    u=B.apply(direction);gu=model.green(u);su=model.S@u;hu=p*(u-chi*gu)
    numerator=-torch.real(torch.sum(gradient.conj()*direction,dim=0))
    denominator=torch.sum(abs(su)**2*wd,dim=0)+torch.sum(abs(hu)**2*ws,dim=0)
    alpha=torch.where(denominator>0,numerator/denominator,torch.zeros_like(numerator))
    gamma=direction*alpha;j=J+u*alpha;q=Q+gu*alpha;m=cuda_material_readout(j,q,ws,6)
    chi=torch.repeat_interleave(m,3)[:,None];p=6/(1+6*chi)
    loss=torch.mean(torch.sum(abs(D+su*alpha)**2,dim=0)*wd+torch.sum(abs(p*(j-chi*q))**2,dim=0)*ws)
    return dict(material=m,gamma=gamma,current=j,field=q,alpha=alpha,loss=loss)

def main():
    psutil.Process().cpu_affinity([12,14]);torch.set_num_threads(1);torch.set_num_interop_threads(1)
    wall,cput=time.perf_counter(),time.process_time();records=[];results=[]
    cpu=OriginalVectorModel(30);model=CudaOriginalModel(cpu);Y=model.array(np.load(OUT/'public'/'observations.npz')['Y'])
    receiver=CudaReceiverSubspace(model,Y,.1);torch.cuda.reset_peak_memory_stats()
    memory_before=torch.cuda.mem_get_info()
    for layout,batch in [('strided',4),('packed',16)]:
        warmwall,warmcpu=time.perf_counter(),time.process_time();model.configure_fft(layout,batch)
        B=CudaImplicitFourier(model,receiver,6);z=torch.zeros((B.K,batch),device='cuda',dtype=torch.complex128)
        model.green(receiver.Jd[:,:batch]);B.adjoint(B.apply(z));torch.cuda.synchronize()
        warm=dict(CPU_seconds=time.process_time()-warmcpu,wall_seconds=time.perf_counter()-warmwall)
        before=dict(model.counts);start,cc=time.perf_counter(),time.process_time()
        result=replay(model,receiver,Y);torch.cuda.synchronize()
        records.append(dict(layout=layout,FFT_block_rhs=batch,CPU_seconds=time.process_time()-cc,
            wall_seconds=time.perf_counter()-start,loss=float(result['loss'].item()),warmup=warm,
            counts={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}))
        results.append(result)
    differences={}
    for key in ['material','gamma','current','field','alpha']:
        a=results[0][key];b=results[1][key]
        differences[key]=dict(relative_L2=float((torch.linalg.vector_norm(a-b)/torch.clamp(torch.linalg.vector_norm(a),min=1e-300)).item()),
            max_absolute=float(torch.max(abs(a-b)).item()))
    result=dict(scope='ONE_BOUNDED_PUBLIC_CURRENT_STEP_REPLAY; no GT; no full solves; uniform0.1 material and same Jd/Qd/Dd',
        records=records,differences=differences,wall_speed_ratio=records[0]['wall_seconds']/records[1]['wall_seconds'],
        numerically_qualified=all(differences[k]['relative_L2']<1e-10 for k in differences),
        CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,counts=dict(model.counts),physics_solves=0,
        process_affinity=psutil.Process().cpu_affinity(),memory_free_total_before=list(memory_before),
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        torch_version=torch.__version__,platform=platform.platform())
    (OUT/'CUDA_LAYOUT_CURRENT_STEP_PROFILE.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
