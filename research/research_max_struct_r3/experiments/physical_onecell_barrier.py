"""One real Hankel disk, exact circular response and legal boundary trap.
A physical mechanism control, never an unknown-image benchmark.
"""
from pathlib import Path
import sys,json,time
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import build_operators,Model

def main():
    start,cpu=time.perf_counter(),time.process_time();out=ROOT/'results/onecell';out.mkdir(parents=True,exist_ok=True)
    if(out/'RECEIPT.json').exists():raise FileExistsError('preserve one physical control')
    ops=build_operators([[0.,0.]],[[-.7,-.6]],[[.2,.8]],.04,2*np.pi)
    model=Model(ops);truth=1.2;seed=20.;Y=model.forward(np.array([truth]))['y'].copy();scale=np.vdot(Y,Y).real
    def fun(x):
        st=model.forward(x);r=st['y']-Y
        return float(.5*np.vdot(r,r).real/scale),model.data_gradient(st,Y,scale)
    val,g=fun(np.array([seed]));native=minimize(fun,[seed],jac=True,method='L-BFGS-B',bounds=[(.02,20.)],options=dict(gtol=1e-10,ftol=1e-14,maxiter=100))
    # Same exact rank-one formula evaluated using the already retained anchor.
    st=model.forward(np.array([seed]));L=1-seed*ops['G'][0,0]
    c=ops['G'][0,0]/L;U=ops['S'][:,0,None]/L*st['Q'][0,None,:]
    r=st['y']-Y;A=np.vdot(r,U);B=np.vdot(U,U).real;u=A.real;v=B-2*np.real(c.conjugate()*A)
    roots=np.roots([-u*abs(c)**2-v*c.real,v,u]);candidates=[.02-seed,20-seed,0.]
    candidates += [float(z.real)for z in roots if abs(z.imag)<1e-10 and .02-seed<=z.real<=20-seed]
    scores=[float(.5*np.linalg.norm(r+(d/(1-c*d))*U)**2/scale)for d in candidates]
    d=candidates[int(np.argmin(scores))];recovered=seed+d
    pred=model.forward(np.array([recovered]));relative=float(np.linalg.norm(pred['y']-Y)/np.linalg.norm(Y))
    # Data-current readout: this one-cell receiver S has full column rank.
    jD=np.linalg.lstsq(ops['S'],Y,rcond=None)[0];Q=ops['E0']+ops['G']@jD
    a=np.vdot(jD,jD).real;b=np.vdot(Q,jD).real;q=np.vdot(Q,Q).real
    cie_readout=(b+a)/(q+b)
    # Saved analytic plotting samples cost no LU; explicitly separate from candidates.
    plot_m=np.linspace(.02,20,600);plot_y=np.array([(ops['S']@(mm*ops['E0']/(1-mm*ops['G'][0,0]))).item()for mm in plot_m])
    g0=ops['G'][0,0];maximum=(g0.real*truth-1)/(abs(g0)**2*truth-g0.real)
    np.savez_compressed(out/'arrays.npz',m=plot_m,y=plot_y,target=Y,seed=seed,true_m=truth,proposal=recovered,
        candidates=np.array(candidates)+seed,candidate_scores=scores,G=ops['G'],E0=ops['E0'],S=ops['S'])
    rec=dict(status='COMPLETED',classification='DISCLOSED_PHYSICAL_MECHANISM_CONTROL_NOT_BLIND_IMAGING',
        physical_model='one radius.04 lambda circular pulse, actual full Hankel self and external source/receiver',time_convention='exp(-i omega t)',
        material_box=[.02,20.],truth=truth,initial=seed,noise='clean',prior='none for onecell theorem/control',
        imaginary_G=float(g0.imag),other_data_stationary_material=float(maximum),
        initial_objective=val,initial_derivative=float(g[0]),local_solution=float(native.x[0]),local_stop=str(native.message),
        rational_candidate_solution=float(recovered),rational_final_relative_residual=relative,
        cie_full_current_beta1_readout=float(cie_readout),
        interpretation='legal passive radiation, no real pole, but upper box endpoint is non-global constrained minimum; exact circle search escapes; mature full-current CIE readout already also solves this control',
        physical_model_cost=dict(model.counts),analytic_candidate_count=len(candidates),closed_form_plot_evaluations=600,
        local_iterations=int(native.nit),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
        finished_utc=datetime.now(timezone.utc).isoformat(),truth_access='DISCLOSED_CONTROL_ONLY')
    (out/'RECEIPT.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
