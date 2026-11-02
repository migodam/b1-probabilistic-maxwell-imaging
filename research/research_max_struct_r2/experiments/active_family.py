"""Separate acquisition/fit entrypoints; restricted public tensor family, not imaging.
All solves reuse qualified R1 vector DDA, no new backend.
"""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
R1=ROOT.parent/'research_max_struct_r1'
sys.path.insert(0,str(R1/'src'))
from vector_backend import TensorDDA
from tensor_families import grid,bump,sphere,DIR,POL
OUT=ROOT/'results/active'; OUT.mkdir(parents=True,exist_ok=True)
def dump(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2))
def solve(t,h,tag):
    _,p=grid(h); eps=np.tile(np.eye(3),(len(p),1,1));eps[:,2,2]+=t*bump(p)
    model=TensorDDA(p,h**3)
    state=model.solve(eps,DIR,POL,OUT/(tag+'_solver.json'))
    return state.near(sphere(12,.6))['E_scattered'],eps,state.receipt
def main(stage):
    start=time.perf_counter();cpu=time.process_time();sigma=.001; epsnoise=6.
    if stage=='plan':
        # Predict using public endpoints only. Selection never reads actual t/Ynew.
        A,e0,r0=solve(0.,.06,'plan_t0');B,e1,r1=solve(3.,.06,'plan_t3')
        _,sv,vh=np.linalg.svd((B[:,2:]-A[:,2:])/sigma,full_matrices=False)
        c=vh[0].conj();np.savez(OUT/'selected_source.npz',c=c)
        dump('PLAN.json',dict(public_family='epsilon=I+t*b(x)zz^T, t in [0,3]',
             status='FROZEN_BEFORE_ACQUISITION',old_columns=[0,1],new_columns=[2,3,4,5],
             source_metric='identity ideal plane-wave coefficient power',sigma=sigma,
             noise='bounded complex ball norm<=sigma*sqrt(36); actual noise half radius',
             predicted_endpoint_gap=float(sv[0]),pair_disjoint_threshold=12.,
             geometry_prior='KNOWN smooth bump, radius .1 lambda; only scalar t unknown',
             grid_h=.06,N=64,source_coeff_real=c.real.tolist(),source_coeff_imag=c.imag.tolist(),
             LU=2,forward_RHS=12))
    elif stage=='acquire':
        # Dedicated generator is allowed truth access; fit is a different invocation.
        rng=np.random.default_rng(2026110202);t=float(rng.uniform(.6,2.4))
        c=np.load(OUT/'selected_source.npz')['c'];Y,eps,rec=solve(t,.06,'acquire')
        noise=rng.normal(size=36)+1j*rng.normal(size=36);noise*=.5*epsnoise*sigma/np.linalg.norm(noise)
        oldnoise=rng.normal(size=(36,2))+1j*rng.normal(size=(36,2));oldnoise*=.5*np.sqrt(72)*sigma/np.linalg.norm(oldnoise)
        np.savez(ROOT/'sealed_truth/active_truth.npz',t=t,epsilon=eps,noise=noise)
        np.savez(ROOT/'data_public/active_observation.npz',Yold=Y[:,:2]+oldnoise,
                 Ynew=Y[:,2:]@c+noise,c=c,sigma=sigma)
        dump('ACQUISITION.json',dict(stage='generation',truth_access='GENERATOR_ONLY',
             LU=1,forward_RHS=6,new_independent_source_combinations=1,
             actual_measurement='independent synthetic generation, same finite DDA; not laboratory',
             old_rhs=2,library_rhs=4,all_computational_rhs_charged=True))
    elif stage=='fit':
        d=np.load(ROOT/'data_public/active_observation.npz');c=d['c'];obs=d['Ynew'];records=[]
        def fun(t):
            Y,eps,rec=solve(t,.06,f'fit_{len(records):03d}')
            val=float(np.linalg.norm(Y[:,2:]@c-obs)/sigma)
            records.append(dict(t=float(t),residual=val,LU=1,RHS=6))
            return val*val
        opt=minimize_scalar(fun,bounds=(0,3),method='bounded',options=dict(xatol=1e-5,maxiter=40))
        np.savez(ROOT/'predictions/active_fit.npz',t=float(opt.x))
        dump('FIT.json',dict(t=float(opt.x),residual=float(np.sqrt(opt.fun)),
            noise_radius=epsnoise,old_family_residual=float(np.linalg.norm(d['Yold'])/sigma),
            old_family_radius=float(np.sqrt(72)),old_solution_set='entire [0,3], exactly same old forward',
            prior='uniform interval, no GT ranking',truth_access=False,status=str(opt.message),
            LU=len(records),forward_RHS=6*len(records),trials=records,
            limitation='bounded scalar search not certified global minimization; public 1D family only'))
    elif stage=='coverage':
        # One predeclared model-mismatch control: add a legal transverse component.
        # Old-source prediction of EVERY candidate in public t-family is zero.
        _,p=grid(.06);eps=np.tile(np.eye(3),(len(p),1,1))
        eps[:,2,2]+=1.5*bump(p);eps[:,0,0]+=bump(p)
        mod=TensorDDA(p,.06**3);st=mod.solve(eps,DIR[:2],POL[:2],OUT/'coverage_solver.json')
        Y=st.near(sphere(12,.6))['E_scattered'];rng=np.random.default_rng(2026110203)
        noise=rng.normal(size=Y.shape)+1j*rng.normal(size=Y.shape)
        noise*=.5*np.sqrt(Y.size)*sigma/np.linalg.norm(noise);obs=Y+noise
        resid=float(np.linalg.norm(obs)/sigma)
        np.savez(OUT/'coverage_observation.npz',Y=obs,epsilon=eps)
        dump('COVERAGE.json',dict(stage='OFFLINE_MODEL_MISMATCH_CONTROL',
             family_best_possible_old_residual=resid,noise_radius=float(np.sqrt(Y.size)),
             all_family_rejected=bool(resid>np.sqrt(Y.size)),LU=1,forward_RHS=2,
             conclusion='failure of restricted family coverage, not optimization failure',
             legal_real_eigen_min=float(np.linalg.eigvalsh(eps).min())))
    elif stage=='score':
        truth=np.load(ROOT/'sealed_truth/active_truth.npz');pred=np.load(ROOT/'predictions/active_fit.npz')
        fit=json.loads((OUT/'FIT.json').read_text())
        # Independent necessary discretization check: SAME fitted and true t on h=.03.
        c=np.load(OUT/'selected_source.npz')['c']
        Yt,et,rt=solve(float(truth['t']),.03,'validation_truth_fine')
        Yp,ep,rp=solve(float(pred['t']),.03,'validation_fit_fine')
        finegap=float(np.linalg.norm((Yt[:,2:]-Yp[:,2:])@c)/sigma)
        dump('SCORE.json',dict(truth_t=float(truth['t']),estimated_t=float(pred['t']),
            absolute_t_error=float(abs(truth['t']-pred['t'])),
            tensor_rms_fine=float(np.linalg.norm(et-ep)/np.sqrt(len(et))),
            fine_prediction_gap=finegap,LU=2,forward_RHS=12,
            interpretation='one realization compatible scalar estimate; no uniqueness certificate; finite DDA only, not CIE advantage',
            independent_solver=False,continuum_certificate=False,truth_access='POST_SEAL_SCORING'))
    dump(stage+'_PROCESS.json',dict(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,stage=stage))
if __name__=='__main__':
    with threadpool_limits(limits=1):main(sys.argv[1])
