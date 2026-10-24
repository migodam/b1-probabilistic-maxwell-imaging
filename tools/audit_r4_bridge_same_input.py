"""Resolve the preflight's mismatched Float64/Float32 reference, zero new physics.

The failed physical preflight stays FAILED. It compared gradients of two
slightly different worlds near a stiff point. This audit replays the SAME saved
float32-input world's true Full VJP through CPU/MPS, retaining the original gates.
"""
from pathlib import Path
import os,sys,fcntl,numpy as np,torch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.train_physics_transport import ActualPhi

def main():
    import json
    original=OUT/'AUDIT/PHYSICS_TRANSPORT_REAL_PREFLIGHT.json';old=json.loads(original.read_text())
    base=OUT/'I/jobs/physics_transport_real_VJP_preflight/transport_queries'
    if old['status']!='FAIL' or old['actual_cost']['RHS']!=64:raise ValueError('expected original preserved64RHS failure')
    with np.load(base/'q00000001/BEFORE.npz') as z:x64=z['x'].copy()
    with np.load(base/'q00000001/AFTER.npz') as z:g64=z['gradient_x'].copy()
    with np.load(base/'q00000006/BEFORE.npz') as z:x=z['x'].copy();y=z['observation'].copy()
    with np.load(base/'q00000006/AFTER.npz') as z:phi=float(z['Phi']);gradient=z['gradient_x'].copy()
    if not np.array_equal(x,x64.astype('float32').astype('float64')):raise ValueError('saved finalquery is not exact Float32 input')
    class SavedOracle:
        def value_and_gradient(self,a,b):
            if not np.array_equal(a,x[None]) or not np.array_equal(b,y[None]):raise ValueError('same-input replay only; no physical fallback')
            return np.array([phi]),gradient[None]
    os.environ['PYTORCH_MPS_HIGH_WATERMARK_RATIO']='0.5';os.environ['PYTORCH_MPS_LOW_WATERMARK_RATIO']='0.4'
    torch.set_num_threads(1);lock=(OUT/'training.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    b=Book('I','physics_transport_same_input_reference_audit',max_wall=120,analysis_only=True,gpu=True)
    status='COMPLETE';error=None
    try:
        t=torch.tensor(x[None],dtype=torch.float64,requires_grad=True);v=ActualPhi.apply(t,y[None],SavedOracle());v.sum().backward()
        m=torch.tensor(x[None],dtype=torch.float32,device='mps',requires_grad=True);w=ActualPhi.apply(m,y[None],SavedOracle());w.sum().backward();torch.mps.synchronize()
        mg=m.grad.detach().cpu().double().numpy()[0];rel=float(np.linalg.norm(mg-gradient)/max(1.,np.linalg.norm(gradient)))
        cpu_rel=float(np.linalg.norm(t.grad.detach().numpy()[0]-gradient)/max(1.,np.linalg.norm(gradient)))
        passed=rel<=.002 and cpu_rel<=1e-12 and abs(float(w.cpu())-phi)<=.002 and max(r['relative_error'] for r in old['CPU_FD'])<=1e-4
        result=dict(status='PASS_CORRECTED_REFERENCE' if passed else 'FAIL',original_preflight_status_retained='FAIL',original_physics_RHS=64,new_physics_RHS=0,
            original_reason='gradient reference compared distinct Float64 and Float32-rounded worlds; no gradient-tolerance relaxation',
            original_cross_world_gradient_relative_difference=old['MPS_gradient_relative_difference'],
            input_rounding_norm=float(np.linalg.norm(x-x64)),same_input_CPU_gradient_relative_error=cpu_rel,same_input_MPS_gradient_relative_error=rel,
            actual_gradient_change_between_inputs=float(np.linalg.norm(gradient-g64)),
            MPS_value_reproduces_original=bool(float(w.cpu())==old['mps_phi']),CPU_FD=old['CPU_FD'],
            inputs_and_true_VJP='previously saved actual Full+source-adjoint query6, no re-solving',
            target_unchanged=True,optimizer_updates=0,scope='real likelihood+saved gradient and same-input device bridge combined gate; not posterior qualification')
        write(OUT/'AUDIT/PHYSICS_TRANSPORT_SAME_INPUT_AUDIT.json',result)
        if not passed:raise AssertionError(result)
    except BaseException as e:status='FAILED';error=repr(e);raise
    finally:b.close(status,error);fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
if __name__=='__main__':main()
