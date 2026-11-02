"""Bounded OFFLINE_DIAGNOSTIC of alternating dynamics, not an imaging test.

Uses frozen genuine 4x4 Hankel--Green operators. No GT dataset is opened.
The deliberately chosen material root is development truth. All counters are
paid operations in this script; explicit dense source-metric products are used.
"""
from pathlib import Path
from collections import Counter
import json
import time
import sys
import importlib.util
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
_frozen_name = '_max_struct_frozen_diversity'
_frozen_path = ROOT/'research/delegated/audit/frozen_diversity_physics.py'
_spec = importlib.util.spec_from_file_location(_frozen_name, _frozen_path)
_frozen = importlib.util.module_from_spec(_spec)
sys.modules[_frozen_name] = _frozen
_spec.loader.exec_module(_frozen)
_adapter_spec = importlib.util.spec_from_file_location('_max_struct_frozen_r3', ROOT/'research/delegated/audit/frozen_r3_physics.py')
_adapter = importlib.util.module_from_spec(_adapter_spec)
_adapter_spec.loader.exec_module(_adapter)
MaxwellModel = _adapter.MaxwellModel

COUNTS = Counter()
STARTUP_FAILURES = [
    dict(attempt='initial_python',error='SyntaxError: unclosed offroot.append parenthesis; corrected before any import or numerical work',
         wall_seconds=None,cpu_seconds=None,reason='Startup process was not instrumented',physical_LU=0,physical_complex_RHS=0),
    dict(attempt='system_python',error='ModuleNotFoundError: numpy',wall_seconds=None,cpu_seconds=None,
         reason='Startup process was not instrumented; no numerical work',physical_LU=0,physical_complex_RHS=0),
    dict(attempt='existing_venv_first_import',error='Missing frozen research_z_alg_r4_diversity package in isolated worktree',
         wall_seconds=None,cpu_seconds=None,reason='Import failed before root preparation; no numerical work',physical_LU=0,physical_complex_RHS=0),
]


def metric_inner(X, Z, omega):
    COUNTS['source_metric_products'] += 1
    return float(np.vdot(X, Z @ omega).real)


def readout(J, Q, m, beta, omega, bounds=(.02, 20.), regularizer=None):
    """Exact shared-real r-LS with full fixed Hermitian source metric.

    With a material prior, solve the actual r residual + R(m), not a prior on r.
    """
    V, T = Q + beta*J, beta*J
    VO = V @ omega
    den = np.sum(VO.conj()*V, axis=1).real
    num = np.sum(VO.conj()*T, axis=1).real
    lo, hi = [np.broadcast_to(np.asarray(x), m.shape) for x in bounds]
    lower, upper = beta*lo/(1+beta*lo), beta*hi/(1+beta*hi)
    r = np.clip(num/den, lower, upper)
    seed = np.clip(r/(beta*(1-r)), lo, hi)
    COUNTS['material_readouts'] += 1
    COUNTS['source_metric_products'] += 1
    if regularizer is None:
        return seed
    def vg(z):
        rr = beta*z/(1+beta*z)
        residual = T-rr[:, None]*V
        dres = -beta*V/(1+beta*z[:, None])**2
        value, gradient = regularizer(z)[:2]
        COUNTS['material_inner_evaluations'] += 1
        return (.5*metric_inner(residual, residual, omega)+value,
                np.sum(dres.conj()*(residual@omega), axis=1).real+gradient)
    start = min((seed, m), key=lambda z: vg(z)[0])
    opt = minimize(vg, start, jac=True, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                   options=dict(maxiter=120, ftol=1e-14, gtol=1e-11))
    COUNTS['material_inner_solves'] += 1
    COUNTS['material_inner_iterations'] += int(opt.nit)
    if not opt.success:
        COUNTS['material_inner_not_success'] += 1
    return opt.x if vg(opt.x)[0] <= vg(start)[0] else start.copy()


def current_exact(G, E, S, Y, m, wd, ws):
    """Exact beta_current=1 quadratic block in fixed full current space."""
    p = 1/(1+m)
    H = p[:, None]*(np.eye(len(m))-m[:, None]*G)
    target = (p*m)[:, None]*E
    out = np.empty_like(E)
    for s in range(E.shape[1]):
        K = wd[s]*(S.conj().T@S)+ws[s]*(H.conj().T@H)
        rhs = wd[s]*(S.conj().T@Y[:, s])+ws[s]*(H.conj().T@target[:, s])
        factor = cho_factor(K, lower=True)
        out[:, s] = cho_solve(factor, rhs)
        COUNTS['current_block_Cholesky'] += 1
        COUNTS['current_block_complex_RHS'] += 1
    COUNTS['current_exact_blocks'] += 1
    return out


def current_step(G, E, S, Y, J, m, od, os, memory, global_pr):
    """Actual coupled objective/gradient and exact line/subspace minimization.

    Sourcewise PR uses one real coefficient per column and a coupled p*p
    line minimizer for full metrics. This is a covariance diagnostic extension,
    not a literal historical independent-source implementation when offdiagonal.
    """
    p = 1/(1+m)
    H = p[:, None]*(np.eye(len(m))-m[:, None]*G)
    R = H@J-(p*m)[:, None]*E
    D = S@J-Y
    grad = S.conj().T@(D@od)+H.conj().T@(R@os)
    old_grad, old_direction = memory if memory is not None else (None, None)
    axis = None if global_pr else 0
    pr = 0.
    if old_grad is not None:
        pr = np.maximum(0., np.sum((grad.conj()*(grad-old_grad)).real, axis=axis)/
                        np.maximum(np.sum(abs(old_grad)**2, axis=axis), 1e-300))
    direction = -grad if old_grad is None else -grad+old_direction*pr
    gd = np.sum((grad.conj()*direction).real, axis=axis)
    if global_pr:
        if gd >= 0:
            direction = -grad
    else:
        bad = gd >= 0
        direction[:, bad] = -grad[:, bad]
    DD, RR = S@direction, H@direction
    if global_pr:
        numerator = -float(np.vdot(grad, direction).real)
        denominator = metric_inner(DD, DD, od)+metric_inner(RR, RR, os)
        alpha = numerator/denominator if denominator > 1e-300 else 0.
    else:
        # Re trace((D diag(a))* (D diag(a)) Omega) has Hessian
        # Re[(D*D)_{st} Omega_{ts}], including every cross-source term.
        K = ((DD.conj().T@DD)*od.T+(RR.conj().T@RR)*os.T).real
        rhs = -np.sum((grad.conj()*direction).real, axis=0)
        alpha = np.linalg.lstsq(K, rhs, rcond=1e-13)[0]
        COUNTS['sourcewise_coupled_line_solves'] += 1
    COUNTS['PR_steps'] += 1
    return J+direction*alpha, (grad, direction), dict(pr=np.asarray(pr).tolist(), alpha=np.asarray(alpha).tolist())


def objective(G, E, S, Y, J, m, od, os):
    p = 1/(1+m)
    R = p[:, None]*(J-m[:, None]*(E+G@J))
    D = S@J-Y
    return .5*(metric_inner(D, D, od)+metric_inner(R, R, os))


def timed_receipt(name, operation, metadata):
    before = COUNTS.copy()
    wall, cpu = time.perf_counter(), time.process_time()
    result = operation()
    result['receipt'] = dict(run_id=name, status='COMPLETED', classification='OFFLINE_DIAGNOSTIC',
        code_hash=None, config_hash=None, hash_null_reason='User forbids hash checks; source path provenance only',
        git_provenance='Isolated max_struct_r1 worktree; frozen base and commit recorded by parent',
        code_path='research/research_max_struct_r1/src/dynamics_audit.py',
        config=metadata, wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu,
        operation_counts=dict(COUNTS-before), physical_LU=0, physical_complex_RHS=0,
        physical_zero_reason='Physical roots prepared once in shared_setup, no additional physical solves',
        FFT=0, FFT_zero_reason='Explicit dense matrices; no FFT called', GPU_seconds=None,
        GPU_null_reason='CPU-only diagnostic', truth_access='Chosen disclosed development root; no dataset truth opened',
        prediction_sealed_at=None, prediction_sealed_null_reason='Not a blind inversion or prediction task',
        online_selection=None, online_selection_null_reason='No imaging candidate or online ranking',
        final_physical_relative_residual=None,
        physical_residual_null_reason='Iteration-map diagnostic, not a solved imaging/prediction task',
        gradient_and_field_product_counts=None,
        product_counts_null_reason='Dense G/S/Gram/regularizer products paid in CPU/wall but not exhaustively counted',
        monetary_cost=None, monetary_cost_null_reason='Local CPU diagnostic; pricing/amortized infrastructure cost not evaluated',
        small_SVD=None, small_SVD_null_reason='Spectral norm and least-squares LAPACK internals not separately instrumented; line solves counted',
        elliptic_solves=0, elliptic_zero_reason='No elliptic subproblem or oracle in this diagnostic',
        failed_attempts=[])
    return result


def d1(model, mstar, Jstar, Y):
    G, E, S = model.G, model.E, model.S
    n, ns = E.shape
    wd, ws = 1/np.sum(abs(Y)**2, axis=0), 1/np.sum(abs(E)**2, axis=0)
    omega = np.diag(ws)
    Qstar = E+G@Jstar
    p = 1/(1+mstar)
    L = np.eye(n)-mstar[:, None]*G
    H = p[:, None]*L
    # Analytic derivative of exact-current LS, then the root material readout.
    M = np.empty((n, n))
    for k in range(n):
        h = np.eye(n)[:, k]
        V = np.empty_like(E)
        for s in range(ns):
            K = wd[s]*(S.conj().T@S)+ws[s]*(H.conj().T@H)
            V[:, s] = cho_solve(cho_factor(K), ws[s]*(H.conj().T@(p*h*Qstar[:, s])))
            COUNTS['root_JVP_Cholesky'] += 1
            COUNTS['root_JVP_complex_RHS'] += 1
        LV = L@V
        M[:, k] = np.sum(ws[None]*(Qstar.conj()*LV).real, axis=1)/np.sum(ws[None]*abs(Qstar)**2, axis=1)
    def amap(m, beta, prior=None, bounds=(.02, 20.)):
        J = current_exact(G, E, S, Y, m, wd, ws)
        return readout(J, E+G@J, m, beta, omega, bounds, prior)
    rng = np.random.default_rng(1951)
    directions = rng.normal(size=(3, n))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    rows = []
    roots = {}
    for beta in (.2, 1., 4.):
        roots[str(beta)] = float(np.linalg.norm(amap(mstar, beta)-mstar))
        for t in (1e-3, 1e-4, 1e-5):
            for k, h in enumerate(directions):
                numerical = (amap(mstar+t*h, beta)-amap(mstar-t*h, beta))/(2*t)
                rows.append(dict(beta_material=beta, t=t, direction=k,
                                 JVP_error=float(np.linalg.norm(numerical-M@h)),
                                 JVP_norm=float(np.linalg.norm(numerical))))
    offroot=[]
    for t in (.01, .05, .2):
        z=mstar+t*directions[0]
        values={str(b):amap(z,b) for b in (.2,1.,4.)}
        offroot.append(dict(t=t, beta_spread=float(max(np.linalg.norm(values[str(b)]-values['1.0']) for b in (.2,4.))),
                            beta_spread_over_t_squared=float(max(np.linalg.norm(values[str(b)]-values['1.0']) for b in (.2,4.))/t**2)))
    prior_rows=[]
    for beta in (.2,1.,4.):
        z=amap(mstar,beta,model.regularizer)
        prior_rows.append(dict(beta_material=beta, root_defect=float(np.linalg.norm(z-mstar)), material=z.tolist()))
    # Put one exact-root cell on the lower bound; projection is directionally
    # nonsmooth. Keep the same genuine physics, with new disclosed root data.
    active=mstar.copy(); active[0]=.02
    st=model.state(active); Ya, Ja=st['Y'].copy(),st['j'].copy()
    wda=1/np.sum(abs(Ya)**2,axis=0)
    def active_map(m,beta):
        J=current_exact(G,E,S,Ya,m,wda,ws)
        return readout(J,E+G@J,m,beta,omega)
    hh=np.zeros(n);hh[0]=1
    t=1e-5
    active_rows=[]
    for beta in (.2,1.,4.):
        # Both +/- input roots may be evaluated by an unconstrained current
        # block; the output is clipped to the box. Negative input is an external
        # diagnostic direction, not a legal optimizer step.
        pos=(active_map(active+t*hh,beta)-active)/t
        neg=(active_map(active-t*hh,beta)-active)/t
        active_rows.append(dict(beta_material=beta, opposite_direction_sum_norm=float(np.linalg.norm(pos+neg)),
                                inward_JVP=pos.tolist(), outside_input_JVP=neg.tolist()))
    # A single reset PR step replaces exact current minimization. At the root
    # alpha is a directional Rayleigh quotient, so the map need not be Frechet.
    def one_pr(h,t):
        m=mstar+t*h
        J,_,_=current_step(G,E,S,Y,Jstar.copy(),m,np.diag(wd),omega,None,True)
        return (readout(J,E+G@J,m,1.,omega)-mstar)/t
    a,b=directions[:2]
    pr_rows=[]
    for t in (1e-3,1e-4,1e-5):
        da,db,dab=one_pr(a,t),one_pr(b,t),one_pr(a+b,t)
        pr_rows.append(dict(t=t, directional_additivity_defect=float(np.linalg.norm(dab-da-db)),
                            exact_current_JVP_difference=float(np.linalg.norm(da-M@a))))
    # The extra active-root physical solve is in this run rather than setup.
    return dict(beta_current=1.,current_space='Fixed full C^(16x6)',material='16 shared real scalars',
        root_fixedpoint_defects=roots, exact_root_JVP=rows, root_material_M=M.tolist(),
        root_M_spectral_radius=float(max(abs(np.linalg.eigvals(M)))),
        root_M_spectral_norm=float(np.linalg.norm(M,2)),offroot=offroot,
        common_prior_control=prior_rows,active_box_control=active_rows,single_PR_control=pr_rows,
        active_root_physical_cost=dict(LU=1,complex_RHS=ns),
        interpretation='Root derivative beta_material invariance does not imply off-root equivalence, prior root preservation, or PR dynamics equality')


def d3(model, mstar, Y):
    G,E,S=model.G,model.E,model.S
    n,ns=E.shape
    rng=np.random.default_rng(829)
    U,_=np.linalg.qr(rng.normal(size=(ns,ns))+1j*rng.normal(size=(ns,ns)))
    od=np.diag(1/np.sum(abs(Y)**2,axis=0))
    os=np.diag(1/np.sum(abs(E)**2,axis=0))
    m0=np.full(n,.1)
    # Use nonzero public-data backprojection current; identical physical state
    # in both coordinates, reset history at the sole fixed stage boundary.
    J0=S.conj().T@Y
    J0/=max(np.linalg.norm(J0),1e-300)
    arms=[]
    arrays={'unitary':U}
    for name,metric,global_pr in [('native_renormalized','native',False),
                                  ('full_metric_sourcewise_PR','full',False),
                                  ('full_metric_global_PR','full',True)]:
        outputs=[]
        for mixed in (False,True):
            Em,Ym=(E@U,Y@U) if mixed else (E,Y)
            if mixed and metric=='full':
                odx,osx=U.conj().T@od@U,U.conj().T@os@U
            elif mixed:
                odx=np.diag(1/np.sum(abs(Ym)**2,axis=0))
                osx=np.diag(1/np.sum(abs(Em)**2,axis=0))
            else:
                odx,osx=od,os
            J=J0@U if mixed else J0.copy()
            m=m0.copy();memory=None
            history=[];mats=[];currents=[]
            initial_objective=objective(G,Em,S,Ym,J,m,odx,osx)
            for it in range(16):
                J,memory,info=current_step(G,Em,S,Ym,J,m,odx,osx,memory,global_pr)
                m=readout(J,Em+G@J,m,.2,osx)
                history.append(dict(iteration=it+1,objective_current_beta=objective(G,Em,S,Ym,J,m,odx,osx),**info))
                mats.append(m.copy());currents.append(J@U.conj().T if mixed else J.copy())
            outputs.append(dict(initial_objective=initial_objective,history=history,materials=np.array(mats),currents=np.array(currents)))
        orig,mixed=outputs
        md=np.linalg.norm(orig['materials']-mixed['materials'],axis=1)
        jd=np.linalg.norm(orig['currents']-mixed['currents'],axis=(1,2))
        arms.append(dict(name=name,global_PR=global_pr,initial_objective_original=orig['initial_objective'],
            initial_objective_mixed=mixed['initial_objective'],material_trajectory_difference=md.tolist(),
            current_trajectory_difference=jd.tolist(),original_history=orig['history'],mixed_history=mixed['history']))
        for side,out in zip(('original','mixed'),outputs):
            arrays[name+'_'+side+'_m']=out['materials'];arrays[name+'_'+side+'_j']=out['currents']
    np.savez(ROOT/'research/delegated/audit/D3_trajectories.npz',**arrays)
    return dict(arms=arms,unitary_columns=ns,iterations_each=16,stages=1,
                stopping='Fixed 16 iterations; no coefficient stop confound',
                sourcewise_line_search='Coupled real ns-by-ns subspace minimizer under full source metrics',
                beta_current=1.,beta_material=.2,history_reset='Once at initial stage, identically in original/mixed',
                trajectories_path='research/delegated/audit/D3_trajectories.npz')


def main():
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=1)
    started,cpu=time.perf_counter(),time.process_time()
    outdir=ROOT/'research/delegated/audit';outdir.mkdir(parents=True,exist_ok=True)
    path=ROOT/'data_public/dynamics_operators_4x4.npz'
    with np.load(path,allow_pickle=False) as data:
        public={k:data[k] for k in ('G','E0','S','grid_shape')}
    model=MaxwellModel(public)
    mstar=.55+.17*np.sin(np.arange(model.n)*.83)+.09*np.cos(np.arange(model.n)*.31)
    state=model.state(mstar);Jstar,Y=state['j'].copy(),state['Y'].copy()
    setup=dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
               physical_LU=1,physical_complex_RHS=model.ns,operators_reassembled=False,
               operator_load='Existing genuine Hankel-Green dev_calibration G/E0/S/grid_shape only',
               root_material=mstar.tolist(),truth_status='Chosen OFFLINE_DIAGNOSTIC root',
               assembly_wall_seconds=None,assembly_null_reason='Frozen operator assembly paid historically; not rerun')
    r1=timed_receipt('D1_exact_current_dynamics',lambda:d1(model,mstar,Jstar,Y),
        dict(beta_current=1.,beta_material=[.2,1.,4.],grid=[4,4],source_columns=model.ns))
    r1['receipt']['physical_LU']=1;r1['receipt']['physical_complex_RHS']=model.ns
    r1['receipt'].pop('physical_zero_reason')
    (outdir/'D1.json').write_text(json.dumps(r1,indent=2))
    r3=timed_receipt('D3_source_unitary_dynamics',lambda:d3(model,mstar,Y),
        dict(beta_current=1.,beta_material=.2,grid=[4,4],iterations=16,arms=3,source_columns=model.ns))
    (outdir/'D3.json').write_text(json.dumps(r3,indent=2))
    total=dict(shared_setup=setup,wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
               physical_LU=2,physical_complex_RHS=2*model.ns,operation_counts=dict(COUNTS),
               failed_attempts=STARTUP_FAILURES,status='OFFLINE_DIAGNOSTIC_COMPLETED',
               excluded_costs='Shell startup and external JSON inspection not timed; historical assembly excluded explicitly')
    (outdir/'dynamics_cost.json').write_text(json.dumps(total,indent=2))
    print(json.dumps(dict(status=total['status'],wall_seconds=total['wall_seconds'],
        D1_spectral_radius=r1['root_M_spectral_radius'],
        D3_final_material_difference={a['name']:a['material_trajectory_difference'][-1] for a in r3['arms']})))


if __name__=='__main__':
    main()
