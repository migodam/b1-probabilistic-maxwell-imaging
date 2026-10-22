"""Explicit Book-B conditional-spine build and full-target HMC diagnostic.

No action occurs on import. --build and --chain are separately charged jobs.
"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from p3l.core import OUT, Book, Target, write, save, from_probit, probit_derivative
from p3l.reference import Posterior
from p3l.spine import PolynomialSpine, SpinePosterior, pulled_metric
from p3l.kernels import ChainState
from p3l.hmc import PreconditionedHMC

PARENT = 'DEV_030'
BUILD_JOB = 'build_property_spine_DEV030'
BUILD = OUT / 'B/jobs' / BUILD_JOB
MODE = OUT / f'A/jobs/reference_{PARENT}_mala_mode_c0/PROBIT_MODE_INITIALIZATION.json'


def freeze(book, protocol):
    """Write-once inputs and controlled source snapshot before any query."""
    path = book.jobdir / 'PROTOCOL.json'
    if path.exists(): raise RuntimeError('existing protocol: no overwrite or automatic replay')
    protocol.update(stage='B', GT_used=False, full_target_unchanged=True,
                    method_label='FROZEN_NONLINEAR_TRIANGULAR_PROPERTY_SPINE_HMC',
                    probability_status='INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE')
    write(path, protocol)
    for rel in ['tools/run_p3l_spine.py', 'src/p3l/spine.py', 'src/p3l/hmc.py',
                'src/p3l/core.py', 'src/p3l/reference.py', 'src/p3l/kernels.py',
                'src/p1_material/prior.py']:
        dst = book.jobdir / 'SOURCE' / rel; dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, dst)
    shutil.copyfile(OUT / 'configs/TARGET.json', book.jobdir / 'TARGET_FROZEN.json')
    shutil.copyfile(ROOT / 'research/p1_material_target_20261018/inputs/geometry.json',
                    book.jobdir / 'GEOMETRY_FROZEN.json')


def center_source():
    doc = json.loads(MODE.read_text())
    rows = [a for a in doc['mode_history'] if a['endpoint'] == doc['mode_endpoint']]
    if len(rows) != 1: raise ValueError('mode endpoint must identify exactly one history x')
    x = np.asarray(rows[0]['x'], dtype=float)
    if x.shape != (66,) or not np.isfinite(x).all(): raise ValueError('finite 66-D center required')
    return x, doc['mode_endpoint']


class NodeLimit(RuntimeError): pass


def build():
    from scipy.optimize import least_squares
    xcenter, source_key = center_source()
    book = Book('B', BUILD_JOB, max_wall=7200); status = 'ERROR'; error = None
    try:
        offsets = np.array([(a, b) for a in [-1., 0., 1.] for b in [-1., 0., 1.]])
        freeze(book, dict(action='build', parent=PARENT, mode_source=str(MODE),
            mode_endpoint=source_key, center=xcenter, scale=.5, grid_offsets=offsets,
            conditional_field_dimension=64, max_fresh_states_per_node=24,
            residual='[whiten(F(z,eta)-y), z]', analytic_jacobian='full receiver J, 32 billed RHS',
            conditional_solver=dict(method='trf', tr_solver='exact', x_scale='jac',
                ftol=1e-8, xtol=1e-8, gtol=1e-6, max_nfev=24),
            beta_basis=['1','t1','t2','t1^2','t1*t2','t2^2'],
            beta_fit='unweighted least squares, globally defined polynomial',
            prior='N66(x), x=(v+beta(eta),eta); never N66(r)',
            approximate_endpoints_not_conditional_MAP_certificate=True,
            budget_scope='Book B only; no transfer from old budgets'))
        shutil.copyfile(MODE,book.jobdir/'MODE_SOURCE_FROZEN.json')
        target = Target(book); y = target.observation(PARENT)
        center = target.evaluate(from_probit(xcenter), 'SPINE_fresh_center_and_full_J')
        Jx = target.jacobian(center) * probit_derivative(xcenter)[None, :]
        Hx = Jx.T @ Jx + np.eye(66)
        linear = -np.linalg.solve(Hx[:64, :64], Hx[:64, 64:])
        target.persist_cache(center, book.jobdir / 'cache/CENTER.npz')
        save(book.jobdir / 'CENTER.npz', x=xcenter, endpoint_key=np.array(center.key),
             Hx=Hx, Blinear=linear, gradient_Phi_x=target.phi_gradient(center,y)*probit_derivative(xcenter))
        # Only small arrays survive the center; its full cache is on disk.
        del center, Jx, Hx
        target.backend()._transfers.clear()
        node_x=[]; node_obj=[]; node_gradient=[]; summaries=[]
        for index, offset in enumerate(offsets):
            eta = xcenter[64:] + .5 * offset
            zstart = xcenter[:64] + linear @ (eta - xcenter[64:])
            folder = book.jobdir / 'nodes' / f'{index:02d}'
            best = [np.inf, None, None]; current = [None, None]; fresh = [0]
            history=[]; result=None; node_status='ERROR'; node_error=None
            def point(z):
                z=np.asarray(z); key=z.tobytes()
                if current[0] == key: return current[1]
                if best[1] is not None and np.array_equal(z,best[1][:64]):
                    current[:] = [key,best[2]]; return best[2]
                if fresh[0] >= 24: raise NodeLimit('frozen 24 fresh state cap')
                fresh[0] += 1
                x=np.r_[z,eta]; p=target.evaluate(from_probit(x), f'SPINE_node_{index:02d}')
                objective=target.phi(p,y)+.5*float(z@z)
                history.append(dict(endpoint=p.key, objective=objective, x=x))
                current[:] = [key,p]
                if objective < best[0]: best[:] = [objective,x.copy(),p]
                return p
            def fun(z): return np.r_[target.white(point(z).state.data-y),z]
            def jac(z): return np.vstack((target.jacobian(point(z))[:,:64],np.eye(64)))
            fatal = None
            try:
                result=least_squares(fun,zstart,jac=jac,method='trf',tr_solver='exact',
                    x_scale='jac',max_nfev=24,ftol=1e-8,xtol=1e-8,gtol=1e-6)
                node_status='SOLVER_TERMINATED_APPROXIMATE_ENDPOINT'
            except NodeLimit as exc:
                node_status='FRESH_STATE_CAP_APPROXIMATE_ENDPOINT';node_error=repr(exc)
            except BaseException as exc:
                node_status='FAILED_PRESERVE_BEST';node_error=repr(exc);fatal=exc
            finally:
                summary=dict(index=index,eta=eta,start_z=zstart,fresh_states=fresh[0],
                    status=node_status,error=node_error,conditional_MAP_certified=False,
                    history=history,scipy_status=None if result is None else int(result.status),
                    scipy_message=None if result is None else str(result.message),
                    scipy_optimality=None if result is None else float(result.optimality))
                if best[2] is not None:
                    bx,bp=best[1:]
                    # Persist the paid best state even if a derivative's next
                    # reservation fails or the controlled wall limit expires.
                    target.persist_cache(bp,folder/'BEST_CACHE.npz')
                    gradient=np.full(64,np.nan)
                    try:
                        gradient=target.phi_gradient(bp,y)[:64]+bx[:64]
                        summary['field_gradient_status']='COMPUTED_TRUE_FULL_LIKELIHOOD'
                    except BaseException as exc:
                        summary.update(field_gradient_status='NOT_COMPUTED_PRESERVED_CACHE',
                            field_gradient_error=repr(exc))
                        if fatal is None: fatal=exc
                    save(folder/'BEST.npz',x=bx,q=bp.q,objective=np.array(best[0]),
                        field_gradient=gradient,endpoint_key=np.array(bp.key))
                    summary.update(best_objective=best[0],best_endpoint=bp.key,
                        field_gradient_inf=float(np.linalg.norm(gradient,np.inf)) if np.isfinite(gradient).all() else None)
                    node_x.append(bx);node_obj.append(best[0]);node_gradient.append(gradient)
                write(folder/'STATUS.json',summary);summaries.append(summary)
                current[:]=[None,None];best[:]=[np.inf,None,None]
                if target.model is not None: target.model._transfers.clear()
            if fatal is not None: raise fatal
        design=np.column_stack((np.ones(9),offsets[:,0],offsets[:,1],
            offsets[:,0]**2,offsets[:,0]*offsets[:,1],offsets[:,1]**2))
        endpoints=np.asarray(node_x)
        fitted,residuals,rank,singular=np.linalg.lstsq(design,endpoints[:,:64],rcond=None)
        errors=design@fitted-endpoints[:,:64]
        save(book.jobdir/'SPINE.npz',coefficients=fitted.T,center=xcenter[64:],scale=np.array(.5),
             offsets=offsets,design=design,node_x=endpoints,node_objective=np.asarray(node_obj),
             node_field_gradient=np.asarray(node_gradient),fit_errors=errors,rank=np.array(rank),
             singular_values=singular,least_squares_residuals=residuals,Blinear=linear)
        write(book.jobdir/'BUILD_SUMMARY.json',dict(status='COMPLETE_FROZEN_PROPOSAL_COORDINATE',
            rank=int(rank),fit_error_frobenius=float(np.linalg.norm(errors)),
            fit_error_max_abs=float(np.max(np.abs(errors))),nodes=summaries,
            conditional_MAP_certified=False,global_bijection=True,log_abs_det_DT=0.,
            likelihood_queries='ALL via Target and Book B',GT_used=False))
        status='COMPLETE'
    except BaseException as exc: error=repr(exc);raise
    finally: book.close(status,error)


def load_spine():
    closed=json.loads((BUILD/'CLOSED.json').read_text())
    if closed['status'] != 'COMPLETE': raise ValueError('spine build must be complete')
    with np.load(BUILD/'SPINE.npz',allow_pickle=False) as z:
        if int(z['rank']) != 6: raise ValueError('rank-six polynomial design required')
        return PolynomialSpine(z['coefficients'],z['center'],float(z['scale']))


def chain(c):
    spine=load_spine()
    source=OUT/f'A/jobs/reference_{PARENT}_mala_mode_c{c}/CHAIN.npz'
    with np.load(source,allow_pickle=False) as z: initial_x=z['x'][0].copy()
    book=Book('B',f'spine_hmc_DEV030_pilot_c{c}',max_wall=7200);status='ERROR';error=None
    try:
        warmup,production=512,1024;refresh=(128,256,384)
        freeze(book,dict(action='chain',parent=PARENT,chain=c,spine_source=str(BUILD/'SPINE.npz'),
            source_initialization=str(source),source_index=0,source_is_not_teacher=True,
            initial_x=initial_x,
            initialization_setup_must_be_disclosed=True,warmup=warmup,production=production,
            mass_refresh_warmup=refresh,trajectory_length_choices=[4,8,12],
            RNG_seed_sequence=[2026102290,c],initial_step_size=.05,dual_averaging_target=.8,
            production_parameters_fixed=True,coordinate='r=(v64,eta2), x=(v+beta(eta),eta)',
            log_density='-Phi(from_probit(x))-.5*x@x; dx=dr',
            covariance='inverse(DT.T @ (Jx.T@Jx+I66) @ DT), true full J',
            CHAIN_q='native (z,u)',CHAIN_x='physical probit (z,eta)',CHAIN_r='transport (v,eta)',
            full_support=True,retained_residual_dimension=64,budget_scope='Book B only'))
        shutil.copyfile(BUILD/'SPINE.npz',book.jobdir/'SPINE_FROZEN.npz')
        target=Target(book);y=target.observation(PARENT);base=Posterior(target,y,probit=True)
        post=SpinePosterior(base,spine)
        p=target.evaluate(from_probit(initial_x),'SPINE_HMC_fresh_paid_initial_state')
        r=spine.inverse(initial_x);x=spine.transform(r)
        # Inverse/forward roundoff is explicit: physics corresponds to exactly
        # the represented starting x used below, never a rounded cache alias.
        if not np.array_equal(x,initial_x):
            p=target.evaluate(from_probit(x),'SPINE_HMC_roundtrip_represented_initial_state')
        base.cache[x.tobytes()]=p;post.pin(r)
        C,Hx,Hr=pulled_metric(target,p,x,r,spine)
        save(book.jobdir/'MASS_INITIAL.npz',covariance=C,Hx=Hx,Hr=Hr,x=x,r=r,
            source_initial_x=initial_x,roundtrip_error=x-initial_x)
        current=ChainState(r,post(r),post.gradient(r))
        rng=np.random.default_rng(np.random.SeedSequence([2026102290,c]))
        h=.05;mu=np.log(10*h);da_t=0;hbar=0.;logbar=np.log(h)
        history=[];qs=[];xs=[];rs=[];phis=[];lps=[];keys=[]
        for i in range(warmup+production):
            if i in refresh:
                x=spine.transform(current.q)
                C,Hx,Hr=pulled_metric(target,post.point(current.q),x,current.q,spine)
                save(book.jobdir/f'MASS_WARMUP_{i}.npz',covariance=C,Hx=Hx,Hr=Hr,x=x,r=current.q)
                h=.05;mu=np.log(10*h);da_t=0;hbar=0.;logbar=np.log(h)
            length=int(rng.choice([4,8,12]))
            tr=PreconditionedHMC(post,post.gradient,h,C,length).step(current,rng)
            current=tr.state;post.pin(current.q);m=tr.metadata
            if i<warmup:
                ratio=m['log_mh_ratio'];alpha=0. if not np.isfinite(ratio) else np.exp(min(0.,ratio))
                da_t+=1;eta=1/(da_t+10);hbar=(1-eta)*hbar+eta*(.8-alpha)
                logh=mu-np.sqrt(da_t)/.05*hbar;weight=da_t**(-.75)
                logbar=weight*logh+(1-weight)*logbar
                h=float(np.clip(np.exp(np.clip(logbar if i==warmup-1 else logh,-20,0)),1e-6,1.))
            p=post.point(current.q);x=spine.transform(current.q)
            qs.append(from_probit(x));xs.append(x);rs.append(current.q.copy())
            phis.append(target.phi(p,y));lps.append(current.log_density);keys.append(p.key);history.append(m)
            if (i+1)%128==0:
                arrays=dict(q=np.asarray(qs),x=np.asarray(xs),r=np.asarray(rs),Phi=np.asarray(phis),
                    log_density=np.asarray(lps),endpoint_key=np.asarray(keys),warmup=np.array(warmup),
                    target_family=np.array('FULL'))
                for name in ['accepted','log_mh_ratio','log_hamiltonian_error','leapfrog_steps',
                             'target_evaluations','gradient_evaluations','divergences','large_energy_error']:
                    arrays[name]=np.asarray([v[name] for v in history])
                save(book.jobdir/'CHAIN.npz',**arrays);target.persist_cache(p,book.jobdir/'cache/CURRENT.npz')
                write(book.jobdir/'CHECKPOINT.json',dict(iteration=i+1,total=warmup+production,
                    RNG=rng.bit_generator.state,step_size=h,production_parameters_frozen=i>=warmup,
                    all_rejections_retained=True,Phi=phis[-1],
                    acceptance_last128=float(np.mean([v['accepted'] for v in history[-128:]]))))
                print(json.dumps(dict(chain=c,iteration=i+1,total=warmup+production,Phi=phis[-1],
                    step_size=h,acceptance_last128=float(np.mean([v['accepted'] for v in history[-128:]])))),flush=True)
        save(book.jobdir/'MASS.npz',covariance=C)
        write(book.jobdir/'TRANSITION_SUMMARY.json',dict(target_family='FULL',
            method='frozen_nonlinear_triangular_property_spine_full_target_HMC',
            production_acceptance=float(np.mean([v['accepted'] for v in history[warmup:]])),
            gradient_invocations=sum(v['gradient_evaluations'] for v in history),
            completed_leapfrog_steps=sum(v['leapfrog_steps'] for v in history),
            numerical_divergences=sum(v['divergences'] for v in history),
            large_energy_errors=sum(v['large_energy_error'] for v in history),step_size=h,
            probability_status='INVARIANT_KERNEL_NOT_A_MIXING_OR_GLOBAL_REFERENCE_CERTIFICATE',
            full_likelihood=True,original_probit_Gaussian66_prior=True,GT_used=False))
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:book.close(status,error)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--build',action='store_true');actions.add_argument('--chain',type=int,choices=range(4))
    args=parser.parse_args()
    if args.build:build()
    else:chain(args.chain)


if __name__=='__main__':main()
