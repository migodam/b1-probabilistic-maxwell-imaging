"""Paid64-RHS check of the actual CPU/MPS likelihood autograd bridge."""
from pathlib import Path
import sys,os,fcntl,json,numpy as np,torch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,Target,write
from p3r4.train_physics_transport import ActualPhi,TargetLikelihoodOracle

def main():
    release=OUT/'configs/PHYSICS_TRANSPORT_PREFLIGHT_RELEASE.json'
    if release.exists():raise FileExistsError('preflight already registered; no replay')
    write(release,dict(released=True,parent='R4DEV_00064',max_RHS=64,
        CPU_FD_directions=2,FD_h=1e-4,CPU_relative_tolerance=1e-4,MPS_relative_gradient_tolerance=.002,MPS_absolute_value_tolerance=.002,
        actual_new_Full_worlds=6,original_Y=True,GT=False,optimizer_updates=0,codec='stored'))
    os.environ['R4_ARRAY_CODEC']='stored'
    os.environ['PYTORCH_MPS_HIGH_WATERMARK_RATIO']='0.5';os.environ['PYTORCH_MPS_LOW_WATERMARK_RATIO']='0.4'
    torch.set_num_threads(1)
    if not torch.backends.mps.is_available():raise RuntimeError('released local MPS resource unavailable; do not silently replace')
    lock=(OUT/'training.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    book=None;status='COMPLETE';error=None
    try:
        book=Book('I','physics_transport_real_VJP_preflight',max_wall=120,gpu=True,kind='physics')
        costs=dict(RHS=0,states=0,LU=0);update=book._update
        def limited(event,**kw):
            if costs['RHS']+kw.get('rhs',0)>64:raise RuntimeError('preflight64RHS hardcap')
            result=update(event,**kw)
            for k,q in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(q,0)
            return result
        book._update=limited;target=Target(book,keep_fields=True);oracle=TargetLikelihoodOracle(target,costs,64)
        with np.load(OUT/'results/contexts/dev/R4DEV_00064.npz') as z:x0=z['x'].copy()
        y=target.observation('R4DEV_00064')
        tensor=torch.tensor(x0[None],dtype=torch.float64,requires_grad=True)
        value=ActualPhi.apply(tensor,y[None],oracle);value.sum().backward()
        cpu_phi=float(value.item());gradient=tensor.grad.detach().numpy()[0].copy()
        rng=np.random.default_rng(20261025141);fd=[]
        for _ in range(2):
            direction=rng.normal(size=66);direction/=np.linalg.norm(direction)
            plus=oracle.values((x0+1e-4*direction)[None],y[None])[0]
            minus=oracle.values((x0-1e-4*direction)[None],y[None])[0]
            estimate=(plus-minus)/2e-4;analytic=float(gradient@direction)
            fd.append(dict(analytic=analytic,FD=float(estimate),relative_error=float(abs(analytic-estimate)/max(1.,abs(analytic),abs(estimate)))))
        t=torch.tensor(x0[None],dtype=torch.float32,device='mps',requires_grad=True)
        v=ActualPhi.apply(t,y[None],oracle);v.sum().backward();torch.mps.synchronize()
        mps_phi=float(v.cpu().item());mps_gradient=t.grad.detach().cpu().double().numpy()[0]
        relative=float(np.linalg.norm(mps_gradient-gradient)/max(1.,np.linalg.norm(gradient)))
        passed=max(q['relative_error'] for q in fd)<=1e-4 and relative<=.002 and abs(mps_phi-cpu_phi)<=.002 and costs['RHS']==64
        result=dict(status='PASS' if passed else 'FAIL',cpu_phi=cpu_phi,mps_phi=mps_phi,CPU_FD=fd,MPS_gradient_relative_difference=relative,
            actual_cost=costs,no_GT=True,optimizer_updates=0,scope='first-derivative real Full likelihood bridge only; no posterior qualification')
        write(OUT/'AUDIT/PHYSICS_TRANSPORT_REAL_PREFLIGHT.json',result)
        if not passed:raise AssertionError(result)
    except BaseException as e:status='FAILED';error=repr(e);raise
    finally:
        if book is not None:book.close(status,error)
        fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
if __name__=='__main__':main()
