"""Paid P1 observation/physics-anchor stage; default allocation remains ZERO.

No call is authorized by this file's existence. It exits before importing a
physical backend or reading offline materials unless an explicit user
allocation is present. No network or neural model is involved.
"""
from pathlib import Path
import sys,json,argparse,traceback
import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.target import OUT,MaxwellTarget,PhysicsBook,BudgetDenied,pack,write_json


def generate_observations(target,protocol):
    """Dataset-only truth access; inversion below receives only y and sigma."""
    with np.load(OUT/'inputs/offline_materials/dev.npz') as z:latents=z['latent']
    with np.load(OUT/'inputs/noise_draws.npz') as z:noise=z['noise']
    for i in [0,1]:
        name=['P1_obs_A','P1_obs_B'][i];path=OUT/f'inputs/online/{name}.npz'
        if path.exists():raise RuntimeError('observation already generated: no repeat truth payment')
        path.parent.mkdir(parents=True,exist_ok=True)
        state=target.forward(latents[i],purpose='dataset_truth_'+name)
        with path.open('wb') as f:np.savez_compressed(f,observation=state.data+noise[i],sigma=target.sigma,target_id=np.array(target.target_id))
        write_json(OUT/f'inputs/offline_materials/{name}_provenance.json',dict(parent_id=protocol['preselected_physics_parent_ids'][i],development_index=i,noise_index=i,exposed_development=True,online_q_or_seed=False))


def reconstruct(target,observation_id):
    """No access to truth arrays, generator seeds, parent latent or truth error."""
    y=target.load_observation(observation_id);P=target.prior;q=np.r_[np.zeros(64),.5,.5]
    state=target.forward(q,purpose=observation_id+'_independent_prior_mean_init');best=(target.phi(state.data,y),q.copy(),state);history=[];used=1
    prior_hessian=np.diag(np.r_[np.ones(64),np.zeros(2)])
    for iteration in range(12):
        A=target.jacobian(q,state);r=pack(state.data-y)*np.sqrt(2)/target.sigma
        g=A.T@r+np.r_[q[:64],0.,0.];diag=np.diag(A.T@A)+np.diag(prior_hessian);damping=max(1e-6,.001*float(np.median(diag)))
        M=np.r_[A,np.diag(np.r_[np.ones(64),np.zeros(2)]),np.sqrt(damping)*np.eye(66)]
        rhs=np.r_[-r,-np.r_[q[:64],0.,0.],np.zeros(66)]
        lo=np.r_[np.full(64,-np.inf),-q[64:]];hi=np.r_[np.full(64,np.inf),1-q[64:]]
        fit=lsq_linear(M,rhs,bounds=(lo,hi),tol=1e-10,max_iter=1000);delta=fit.x
        delta*=min(1,.5/max(np.linalg.norm(delta[:64]),1e-30),.2/max(np.linalg.norm(delta[64:],np.inf),1e-30))
        pg=g.copy();pg[64:][((q[64:]<=1e-10)&(g[64:]>=0))|((q[64:]>=1-1e-10)&(g[64:]<=0))]=0
        phi=target.phi(state.data,y);objective=phi+.5*np.dot(q[:64],q[:64]);row=dict(iteration=iteration,Phi=phi,white_residual=float(np.sqrt(2*phi)),posterior_objective_up_to_constant=objective,projected_gradient_inf=float(max(abs(pg))),material_rank=int(np.linalg.matrix_rank(np.r_[P.derivative(q).real,P.derivative(q).imag]/24)),prior_logdensity=P.logprior(q),states_used=used)
        history.append(row);write_json(OUT/f'results/physical/{observation_id}_history.json',history)
        if used>=12:break
        accepted=False
        for alpha in [1.,.5,.25]:
            if used>=12:break
            qnew=q+alpha*delta
            if not P.legal(qnew):raise RuntimeError('bound-aware solver returned an illegal latent')
            new=target.forward(qnew,purpose=observation_id+f'_LM_{iteration}_{alpha}');used+=1;phinew=target.phi(new.data,y);objnew=phinew+.5*np.dot(qnew[:64],qnew[:64])
            row.setdefault('trials',[]).append(dict(alpha=alpha,Phi=phinew,objective=objnew))
            if phinew<best[0]:best=(phinew,qnew.copy(),new)
            if objnew<objective:q=qnew;state=new;accepted=True;break
        write_json(OUT/f'results/physical/{observation_id}_history.json',history)
        if not accepted:break
    phi,qbest,base=best;tau=float(np.sqrt(chi2.ppf(.95,512)))
    # A feasible anchor is not claimed to be a stationary MAP or posterior sample.
    report=dict(observation_id=observation_id,qualification='prior_legal + original fixed-noise feasibility + backward residual',qualified=bool(P.legal(qbest) and 2*phi<=tau*tau and base.residual<=1e-10),q=qbest,Phi=phi,white_residual=float(np.sqrt(2*phi)),tau=tau,prior_logdensity=P.logprior(qbest),states_used=used,GT_online=False,NN_used=False,N3_warm_start=False,stationary_MAP='NOT_CERTIFIED',posterior_samples=False)
    write_json(OUT/f'results/physical/{observation_id}_qualification.json',report)
    return report


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--authorization',type=Path,default=OUT/'configs/EXECUTION_AUTHORIZATION.json');args=ap.parse_args()
    book=PhysicsBook(args.authorization)
    book.require(states=1,RHS=8)  # fails before offline truth/backend access at zero allocation
    error=None
    try:
        target=MaxwellTarget(book)
        protocol=json.loads((OUT/'configs/PROTOCOL.json').read_text());generate_observations(target,protocol)
        results=[reconstruct(target,name) for name in ['P1_obs_A','P1_obs_B']]
        write_json(OUT/'results/physical/QUALIFICATION_SUMMARY.json',results)
    except BaseException as e:
        error=repr(e);write_json(OUT/'results/PHYSICS_FAILURE.json',dict(error=error,traceback=traceback.format_exc(),totals=book.totals()));raise
    finally:book.close(error)

if __name__=='__main__':main()
