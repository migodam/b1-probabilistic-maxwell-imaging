"""Ordinary Full-J box Gauss-Newton supplement; no primary B0 mutation.

All 110 material coordinates stay active. Shared receiver adjoints form the
full Jacobian; full nonlinear P3 FEM objectives accept every proposed step.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import minimize
from .fem_model import PublicWall,physical_lock
from skfem import asm
from .fwi import FWIOptions,regularization,dump,plain
from .risk import prior_hessian
from wslam.fem_backend import _contrast_mass,_forcing


@dataclass
class GNOptions:
    max_outer:int=20
    initial_damping:float=1e-2
    trial_alphas:tuple=(1.,.5,.25)
    kkt_tolerance:float=1e-4
    quadratic_maxiter:int=2000
    quadratic_gtol:float=1e-9
    stagnation_outer:int=5


def full_jacobian(model,z):
    """Three full LU states; 4 primal +16 shared receiver-adjoint RHS/freq."""
    z=model.validate(z);nw,nt=model.nwall,model.ntarget
    wall=PublicWall(model.wx,model.wy,z[:nw].reshape(model.wall_shape),z[nw:2*nw].reshape(model.wall_shape))
    pieces=[];data=[];diagnostics=[]
    with physical_lock(model.book):
        for ii,ev in enumerate(model.evaluators):
            model.book.check();p=ev.p;inside,cells=model.cell_data[ii];inc=model.incident[ii]
            cpu=time.process_time();started=time.perf_counter()
            ct=np.zeros(len(cells),complex)
            ct[inside]=z[2*nw:2*nw+nt][cells[inside]]+1j*z[2*nw+nt:][cells[inside]]/p.freq
            combined=wall.values(p.material_points,p.freq)-1.+ct
            mass=asm(_contrast_mass,p.material_basis,contrast=combined.reshape(p._material_shape))
            A=(p._background_matrix-p.k**2*mass[p._free_dofs][:,p._free_dofs]).tocsc()
            b=np.column_stack([asm(_forcing,p.material_basis,load=(p.k**2*combined*inc[:,src]).reshape(p._material_shape))[p._free_dofs] for src in range(4)])
            factor,_,storage=ev.factor(A,label='FullGN cached-combined full state LU')
            u=factor.solve(b,label='FullGN full state primal4');primal_defect=factor.last_defect
            d=ev.direct+np.asarray(ev.R@u).T
            adj=factor.solve(ev.R.conj().T.toarray(),trans='H',label='FullGN shared receiver adjoint16')
            cache=model.quadrature_cache[ii];P=cache['P'];total=cache['incident']+np.asarray(P@u);av=np.asarray(P@adj)
            product=p.k**2*np.einsum('qr,qs,q->qsr',av.conj(),total,cache['weights'],optimize=True)
            gw=np.zeros((nw,4,16),complex);gt=np.zeros((nt,4,16),complex)
            wm,tm=cache['wall'],cache['target']
            np.add.at(gw,cache['wall_cells'][wm],product[wm]);np.add.at(gt,cache['target_cells'][tm],product[tm])
            J=np.column_stack([gw.reshape(nw,-1).T,(1j/p.freq*gw).reshape(nw,-1).T,gt.reshape(nt,-1).T,(1j/p.freq*gt).reshape(nt,-1).T])
            pieces.append(J);data.append(d)
            diagnostics.append(dict(frequency=p.freq,ndofs=p.ndofs,factor_bytes=storage,primal_backward=primal_defect,adjoint_backward=factor.last_defect,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started))
    return np.vstack(pieces),np.asarray(data),diagnostics


def projected_gradient(q,g,lower,upper):
    result=g.copy();result[((q<=lower+1e-10)&(g>0))|((q>=upper-1e-10)&(g<0))]=0.
    return float(np.linalg.norm(result,np.inf))


def refine(model,y,baseline,out,prior_options,gn_options=None):
    baseline=Path(baseline);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'ONLINE_SEALED.json').exists():raise RuntimeError('supplement already sealed')
    if not (baseline/'ONLINE_SEALED.json').exists():raise RuntimeError('warm start must be sealed B0')
    b0seal=json.loads((baseline/'ONLINE_SEALED.json').read_text());b0receipt=json.loads((baseline/'RECEIPT.json').read_text())
    if not b0seal.get('full_objective_comparable'):raise ValueError('B0 endpoint is not full-band comparable')
    with np.load(baseline/'endpoint.npz',allow_pickle=False) as saved:z=model.validate(saved['z']);data=saved['data'].copy()
    y=np.asarray(y);o=GNOptions() if gn_options is None else gn_options
    if isinstance(o,dict):o=GNOptions(**o)
    prior=prior_options if isinstance(prior_options,FWIOptions) else FWIOptions(**prior_options)
    lo,hi,scale=model.bounds_and_scale();lower,upper=lo/scale,hi/scale
    def objective(zz,pred):
        pv,pg,parts=regularization(zz,model.nwall,model.target_shape,model.wall_shape,prior)
        dv=.5*float(np.sum(abs(pred-y)**2))/model.sigma**2
        return dv+pv,dv,pv,pg,parts
    value,dv,pv,_,_=objective(z,data);initial_value=value
    if abs(value-b0seal['objective'])>1e-7:raise ValueError('B0 public objective does not match frozen prior')
    dump(out/'ONLINE_INPUT_CONTRACT.json',dict(case_id=model.config.get('case_id'),method='FULL_GN_TV_SUPPLEMENT',cohort='DATA_EXPOSED_MEASUREMENT_ONLY_REFINEMENT',truth_access=False,baseline=str(baseline),initialization='already sealed B0 minimum measured objective endpoint',prior_options=asdict(prior),gn_options=asdict(o),parameter_coordinates='all110 independent material coordinates; no SVD coordinate cut',full_jacobian='4 primal +16 shared receiver adjoint RHS per frequency',trial_objective='full original cached-combined P3 FEM plus identical TV/spatial prior',config=model.config))
    damping=o.initial_damping;accepted=0;outer_count=0;trial_count=0;fullj_count=0;reject_streak=0;stop='MAX_OUTER_NONCONVERGED';fatal=None;last_J_z=None;last_kkt=None;last_outer=None;last_linearization=None
    def trace(row):
        with (out/'TRACE.jsonl').open('a') as handle:handle.write(json.dumps(dict(time_unix=time.time(),**row),default=plain,allow_nan=False)+'\n')
    def save_endpoint(zz,pred,phi,data_value):
        np.savez_compressed(out/'endpoint.npz',z=zz,data=pred,objective=phi,data_objective=data_value,frequency_indices=np.arange(len(model.evaluators)),**model.maps(zz))
        np.savez_compressed(out/'endpoint_last.npz',z=zz,data=pred,objective=phi,data_objective=data_value,frequency_indices=np.arange(len(model.evaluators)),**model.maps(zz))
    def state_jacobian(iteration,diagnostic_only=False):
        nonlocal fullj_count,last_J_z,last_kkt,last_outer,last_linearization,value,data,dv,pv
        J,data,diag=full_jacobian(model,z);fullj_count+=1
        value,dv,pv,pg,parts=objective(z,data)
        residual=(data-y).ravel();rr=np.r_[residual.real,residual.imag]/model.sigma;A=np.vstack([J.real,J.imag])/model.sigma
        gradient=(A.T@rr+pg)*scale;q=z/scale;kkt=projected_gradient(q,gradient,lower,upper)
        H=(A.T@A+prior_hessian(model,z,prior))*scale[:,None]*scale[None,:]
        parameter_names=np.asarray([f'wall_eps_{j}' for j in range(model.nwall)]+[f'wall_loss_{j}' for j in range(model.nwall)]+[f'target_contrast_{j}' for j in range(model.ntarget)]+[f'target_loss_{j}' for j in range(model.ntarget)])
        raw=dict(parameter_names=parameter_names,parameter_scale=scale,lower_bounds=lo,upper_bounds=hi,frequencies=np.asarray(model.config['frequencies']),sigma_real_imag=model.sigma,real_stack_order=np.asarray('Re(all frequency/source/receiver rows), then Im(all rows), divided by sigma'),z=z.copy(),data=data.copy(),J_complex=J,J_real_whitened=A,residual_real_whitened=rr,gradient_scaled=gradient,H_gn_plus_exact_prior_scaled=H,objective=value,scaled_projected_gradient_inf=kkt,iteration=iteration)
        if fullj_count==1:np.savez_compressed(out/'FULL_J_FIRST.npz',**raw)
        np.savez_compressed(out/'FULL_J_LAST.npz',**raw)
        last_J_z=z.copy();last_kkt=kkt;last_outer=iteration;last_linearization=(q.copy(),gradient.copy(),H.copy(),kkt)
        save_endpoint(z,data,value,dv)
        trace(dict(kind='FULL_J_STATE',outer=iteration,diagnostic_only=diagnostic_only,objective=value,data_objective=dv,prior=pv,prior_components=parts,scaled_projected_gradient_inf=kkt,damping=damping,diagnostics=diag))
        return q,gradient,H,kkt
    save_endpoint(z,data,value,dv)
    try:
        for iteration in range(o.max_outer):
            model.book.check();outer_count=iteration+1
            if last_J_z is not None and np.array_equal(last_J_z,z):
                q,g,H,kkt=last_linearization
                trace(dict(kind='REUSE_UNCHANGED_FULL_J',outer=iteration,paid_J_outer=last_outer,damping=damping))
            else:q,g,H,kkt=state_jacobian(iteration)
            if kkt<=o.kkt_tolerance:stop='LOCAL_PROJECTED_GRADIENT_TOLERANCE';break
            damped=H+damping*np.eye(model.nparam)
            def quadratic(direction):return .5*float(direction@damped@direction)+float(g@direction),damped@direction+g
            result=minimize(quadratic,np.zeros(model.nparam),jac=True,method='L-BFGS-B',bounds=list(zip(lower-q,upper-q)),options=dict(maxiter=o.quadratic_maxiter,gtol=o.quadratic_gtol,ftol=1e-14,maxls=30))
            direction=result.x;accepted_trial=False
            trace(dict(kind='BOX_GN_SUBPROBLEM',outer=iteration,optimizer_success=bool(result.success),optimizer_stop=str(result.message),iterations=int(result.nit),step_scaled_norm=float(np.linalg.norm(direction)),damping=damping))
            for alpha in o.trial_alphas:
                zz=np.clip((q+alpha*direction)*scale,lo,hi)
                if np.array_equal(zz,z):
                    trace(dict(kind='ZERO_STEP',outer=iteration,alpha=alpha));continue
                prediction=model.predict(zz);trial_count+=1
                vv,dd,pp,_,parts=objective(zz,prediction.data)
                actual=value-vv;predicted=-float(alpha*g@direction+.5*alpha*alpha*direction@H@direction)
                okay=actual>max(1e-10,1e-12*abs(value));ratio=actual/max(predicted,1e-15)
                trace(dict(kind='FULL_FEM_TRIAL',outer=iteration,alpha=alpha,objective=vv,data_objective=dd,prior=pp,actual_reduction=actual,predicted_reduction=predicted,ratio=ratio,accepted=okay,diagnostics=prediction.diagnostics))
                np.savez_compressed(out/'trial_last.npz',z=zz,data=prediction.data,objective=vv,alpha=alpha,outer=iteration)
                if okay:
                    z=zz;data=prediction.data;value=vv;dv=dd;pv=pp;accepted+=1;accepted_trial=True;reject_streak=0
                    damping=max(1e-10,damping*(.5 if ratio>.75 else 2. if ratio<.25 else 1.))
                    save_endpoint(z,data,value,dv);break
            if not accepted_trial:
                damping=min(1e12,damping*10.);reject_streak+=1
                if reject_streak>=o.stagnation_outer:stop='STAGNATION_NONCONVERGED';break
        # One endpoint diagnostic is explicit and counted separately from the
        # maximum twenty optimization outer iterations; no step follows it.
        if last_J_z is None or not np.array_equal(last_J_z,z):
            state_jacobian(outer_count,diagnostic_only=True)
        if last_kkt is not None and last_kkt<=o.kkt_tolerance:stop='LOCAL_PROJECTED_GRADIENT_TOLERANCE'
    except BaseException as error:
        fatal=dict(type=type(error).__name__,message=str(error));stop='INTERRUPTED_NONCONVERGED';model.book.failure(model.phase,'FullGN supplement',error)
        if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
    finally:
        save_endpoint(z,data,value,dv)
        cost=model.book.snapshot(pid=__import__('os').getpid())
        dump(out/'ONLINE_SEALED.json',dict(time_unix=time.time(),case_id=model.config.get('case_id'),method='FULL_GN_TV_SUPPLEMENT',cohort='DATA_EXPOSED_MEASUREMENT_ONLY_REFINEMENT',status='COMPLETE_BOUNDED_ATTEMPT' if fatal is None else 'PARTIAL_INTERRUPTED',truth_access=False,template_used=False,full_objective_comparable=True,initial_objective=initial_value,objective=value,selection='monotone full physical measured objective plus identical frozen TV/spatial prior; B0 warm start retained',stop=stop,outer_iterations=outer_count,accepted_steps=accepted,trial_objective_calls=trial_count,full_jacobian_states=fullj_count,final_J_matches_endpoint=last_J_z is not None and np.array_equal(last_J_z,z),selected_scaled_projected_gradient_inf=last_kkt if last_J_z is not None and np.array_equal(last_J_z,z) else None,local_kkt_tolerance=o.kkt_tolerance,global_convergence_claim=False,error=fatal,gn_cost=cost,b0_cost={key:b0receipt[key] for key in ['full_fe','rhs','wrapper_process_cpu_seconds','wall_seconds']},cumulative_full_fe=b0receipt['full_fe']+cost['full_fe'],cumulative_rhs=b0receipt['rhs']+cost['rhs'],not_cold_start=True,independent_blind_claim=False))
    return json.loads((out/'ONLINE_SEALED.json').read_text())
