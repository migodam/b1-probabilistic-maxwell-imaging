"""Public full-physics curvature at frozen wrong endpoints; no GT access.

This is a finite-model diagnostic, not a global-convexity certificate. Use the
exact smooth-TV Hessian, not its positive IRLS majorizer.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.linalg import lu_solve,eigh
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r8.experiments.run_case import load_public

def exact_hessian(model,m,Y,free):
    value,g=model.data_value_gradient(m,Y);st=model.state(m)
    Q,lu=st['Q'],st['lu'];den=float(np.vdot(Y,Y).real)
    with model._physics.operation('CURVATURE_BASE_ADJOINT',rhs=model.ns,adjoint_RHS=model.ns):
        a=lu_solve(lu,model.S.conj().T@(st['Y']-Y),trans=2,check_finite=False)
    D=model.reg.D;d=D@m;q=np.sqrt(d*d+model.reg.epsilon**2)
    weights=(2*model.reg.smooth+model.reg.tv*model.reg.epsilon**2/q**3)/model.reg.edges
    Hess=np.empty((len(free),len(free)))
    for begin in range(0,len(free),32):
        idx=free[begin:begin+32];v=np.zeros((len(idx),model.n));v[np.arange(len(idx)),idx]=1
        rhs=(v[:,:,None]*Q[None]).transpose(1,0,2).reshape(model.n,-1)
        with model._physics.operation('CURVATURE_TANGENT',rhs=rhs.shape[1],tangent_RHS=rhs.shape[1]):
            dj=lu_solve(lu,rhs,check_finite=False)
        dq=(model.G@dj).reshape(model.n,len(idx),model.ns).transpose(1,0,2)
        va=(v[:,:,None]*a[None]).transpose(1,0,2).reshape(model.n,-1)
        arhs=model.S.conj().T@(model.S@dj)+model.G.conj().T@va
        with model._physics.operation('CURVATURE_ADJOINT_TANGENT',rhs=arhs.shape[1],adjoint_RHS=arhs.shape[1]):
            da=lu_solve(lu,arhs,trans=2,check_finite=False).reshape(model.n,len(idx),model.ns).transpose(1,0,2)
        hv=np.real(np.sum(a.conj()[None]*dq+da.conj()*Q[None],axis=2))/den
        hv+=(D.T@(weights[:,None]*(D@v.T))).T
        Hess[:,begin:begin+len(idx)]=hv[:,free].T
    return value,g,Hess

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser();ap.add_argument('--case-id',required=True);ap.add_argument('--arm',default='a_plus')
    args=ap.parse_args();case=args.case_id
    manifests=[ROOT/'configs/baseline_historical_fixed.json',ROOT/'configs/baseline_cost_controls.json']
    tasks=[t for p in manifests for t in json.loads(p.read_text())['tasks']]
    task=next(t for t in tasks if t['case_id']==case and t['arm']==args.arm)
    inp=Path(task['output_dir']);out=ROOT/'RESULTS/curvature'/(case+'__'+args.arm);out.mkdir(parents=True,exist_ok=False)
    p=load_public(task['public_files']);model=MaxwellModel(p);Y=p['Y']
    with np.load(inp/'prediction.npz') as z:m=z['m'].copy()
    free=np.flatnonzero((m>.02+1e-5)&(m<20-1e-5))
    value,g,H=exact_hessian(model,m,Y,free)
    sym=float(np.linalg.norm(H-H.T)/max(np.linalg.norm(H),np.finfo(float).tiny))
    vals,vecs=eigh((H+H.T)/2,check_finite=False)
    v=np.zeros(model.n);v[free]=vecs[:,0];ray=float(v[free]@H@v[free])
    residual=float(np.linalg.norm(H@vecs[:,0]-vals[0]*vecs[:,0]))
    # One new derivative wiring check along a deterministic free direction.
    w=np.sin(np.arange(len(free))+.7);w/=np.linalg.norm(w);d=np.zeros(model.n);d[free]=w
    h=min(1e-4,float(np.min(np.minimum(m[free]-.02,20-m[free])/(2*np.maximum(abs(w),1e-30)))))
    gp=model.data_value_gradient(m+h*d,Y)[1];gm=model.data_value_gradient(m-h*d,Y)[1]
    fd=(gp[free]-gm[free])/(2*h);hv=H@w
    check=float(np.linalg.norm(fd-hv)/max(np.linalg.norm(fd),np.linalg.norm(hv)))
    proposals=[];best=dict(Phi=value,m=m.copy(),step=0.)
    # This fixed finite probe is conditional on a genuinely negative computed
    # free-space Rayleigh quotient; no attempt to infer PSD from no improvement.
    if ray<0:
        v/=np.max(abs(v))
        for step in [-3.,-.8,-.2,-.05,.05,.2,.8,3.]:
            candidate=np.clip(m+step*v,.02,20.)
            phi=model.physical_score(candidate,Y)
            proposals.append(dict(step=step,Phi=phi))
            if phi<best['Phi']:best=dict(Phi=phi,m=candidate,step=step)
    np.savez_compressed(out/'curvature.npz',free=free,H=H,eigenvalues=vals,weakest_direction=v,anchor=m,proposal=best['m'])
    record=dict(status='COMPLETED_PUBLIC_DIAGNOSTIC',case_id=case,arm=args.arm,truth_access='NONE',
        Phi=value,projected_gradient_inf=float(np.linalg.norm(m-np.clip(m-g,.02,20.),np.inf)),
        free_coordinates=len(free),active_coordinates=model.n-len(free),eigenvalue_min=float(vals[0]),
        eigenvalue_max=float(vals[-1]),negative_eigenvalues=int(np.sum(vals<0)),rayleigh=ray,
        eigenpair_residual=residual,Hessian_symmetry_relative_defect=sym,HVP_directional_relative_error=check,
        proposals=proposals,selected_proposal_step=best['step'],selected_proposal_Phi=best['Phi'],
        scope='Exact finite-model smooth objective Hessian restricted to coordinates not within1e-5 ofbox; not continuous or global convexity, active-set opening not tested',
        cost=dict(physical=dict(model.cost),dense_eigendecompositions=1,dimension=len(free),
                  CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
    (out/'RECEIPT.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
