"""Explicit parent-run B pilots. Importing this file never performs physics."""
from pathlib import Path
import argparse,json,os,sys
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,Target,write,save,plain,to_probit,from_probit,probit_derivative
from p3r2.non_gaussian import quadrature_moments,mixture_from_quadrature,independent_mh,NonlinearChartMixture
from p3l.diagnostics import material_diagnostics


def finite(value):
    if isinstance(value,dict): return {str(k):finite(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [finite(v) for v in value]
    if isinstance(value,np.ndarray): return finite(value.tolist())
    if isinstance(value,np.generic): return finite(value.item())
    if isinstance(value,float) and not np.isfinite(value): return None
    return plain(value)


def freeze(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f: json.dump(finite(value),f,indent=2,allow_nan=False);f.write('\n')


def summarize_chains(states,prior,*,slice_anchor=None,burn=256):
    production=states[:,burn:]
    if production.shape[1]<4: raise ValueError('four retained draws required')
    if slice_anchor is None: q=np.array([from_probit(x) for x in production.reshape(-1,66)]).reshape(production.shape)
    else:
        q=np.broadcast_to(slice_anchor,(*production.shape[:2],66)).copy();q[:,:,64:]=production
    chi=np.array([prior.decode(x)[0] for x in q.reshape(-1,66)]).reshape(*q.shape[:2],576)
    stats=material_diagnostics(chi)
    return chi,stats,dict(material_mean=chi.mean(axis=(0,1)),material_variance=np.mean(np.abs(chi-chi.mean(axis=(0,1)))**2,axis=(0,1)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('slice','full'),required=True);parser.add_argument('--job',required=True)
    parser.add_argument('--parent',default='DEV_030');parser.add_argument('--anchor',type=Path,required=True)
    parser.add_argument('--seed',type=int,default=2026102291);parser.add_argument('--rhs-cap',type=int,default=131072)
    parser.add_argument('--draws',type=int,default=1024);parser.add_argument('--chains',type=int,default=4)
    parser.add_argument('--burn',type=int,default=256);parser.add_argument('--max-wall',type=int,default=7200)
    parser.add_argument('--orders',default='9,17,33,65');parser.add_argument('--refinement-rms',type=float,default=1e-6)
    args=parser.parse_args();orders=tuple(int(v) for v in args.orders.split(','))
    if args.rhs_cap>131072 or args.rhs_cap<8 or args.draws<=args.burn+3 or args.chains<2: raise ValueError('bounded pilot and diagnostics require valid caps/draws/chains')
    if any(o<3 for o in orders) or any(b<=a for a,b in zip(orders,orders[1:])): raise ValueError('increasing quadrature orders >=3 required')
    maximum_states=args.chains*(args.draws+1)+(sum(o*o for o in orders) if args.mode=='slice' else 0)
    if maximum_states*8>args.rhs_cap: raise ValueError('frozen maximum pilot states exceed RHS cap')
    with np.load(args.anchor,allow_pickle=False) as z:
        q0=z['q'].copy();has_T='receiver_T' in z.files and bool(z['receiver_T'].size)
    if args.mode=='full' and not has_T: raise ValueError('FULL pilot requires already-paid receiver transfer; no unregistered J solves')
    book=Book('B',args.job,max_wall=args.max_wall);status='COMPLETE';error=None
    calls=0;dest=book.jobdir;rows=[];transcript=None
    try:
        target=Target(book);y=target.observation(args.parent)
        protocol=dict(mode=args.mode,job=args.job,parent=args.parent,target_id=target.cfg['target_id'],anchor=str(args.anchor.resolve()),
            cache_provenance='Old paid DEV_030 only; old ledger read-only; all new evaluations billed new B',GT_used=False,
            source_isolated_online_observation=True,seed=args.seed,chains=args.chains,draws=args.draws,burn=args.burn,
            maximum_full_states=maximum_states,maximum_new_RHS=maximum_states*8,pilot_RHS_cap=args.rhs_cap,
            original_P1_target_unchanged=True,reference_qualification=False,
            status='FROZEN_BEFORE_FIRST_NEW_LIKELIHOOD',orders=orders if args.mode=='slice' else None,
            quadrature_refinement_rms=args.refinement_rms,
            target_family='SLICE_U2_FULL_NATIVE_BOX_FIXED_Z64' if args.mode=='slice' else 'FULL66_ORIGINAL_PROBIT_PUSHFORWARD',
            original_FULL_posterior=args.mode=='full',
            reference_measure='du_real du_imag; Gaussian64 fixed' if args.mode=='slice' else 'dx66; original Gaussian64 x Uniform2 under exact probit',
            warning='Exact likelihood and MH correction do not establish finite-run mixing; SLICE is not a FULL marginal')
        rng=np.random.default_rng(args.seed)
        # count/evaluate enforces a second independent local cap before Book reserves.
        def evaluate(q,purpose):
            nonlocal calls
            if (calls+1)*8>args.rhs_cap: raise RuntimeError('local B pilot RHS cap')
            calls+=1;p=target.evaluate(q,purpose=purpose);phi=target.phi(p,y);key=p.key
            del p
            if target.model is not None: target.model._transfers.clear()
            return phi,key
        if args.mode=='full':
            anchor=target.load_cache(args.anchor,owner='P3L old paid DEV_030 anchor; immutable prior ledger')
            J=target.jacobian(anchor).copy();del anchor
            x0=to_probit(q0,initialization=True)
            Hzz=np.eye(64)+J[:,:64].T@J[:,:64]
            covariance=np.linalg.solve(Hzz,np.eye(64));L=np.linalg.cholesky((covariance+covariance.T)/2)
            B=-np.linalg.solve(Hzz,J[:,:64].T@J[:,64:])
            Jeta=J[:,64:]*probit_derivative(x0)[64:]
            cross=J[:,:64].T@Jeta
            schur=np.eye(2)+Jeta.T@Jeta-cross.T@np.linalg.solve(Hzz,cross)
            Ceta=np.linalg.inv((schur+schur.T)/2)
            scales=np.sqrt(np.maximum(np.diag(Ceta),1e-8))
            offsets=np.array([[0,0],[.6,0],[-.6,0],[0,.6],[0,-.6]])*scales
            centers=np.tile(x0,(5,1));centers[:,64:]+=offsets
            centers[:,:64]+=(ndtr(centers[:,64:])-ndtr(x0[64:]))@B.T
            proposal=NonlinearChartMixture(centers,L,scales,B,np.ones(5),df=7.)
            protocol.update(proposal='Five frozen nonlinear triangular Student-t charts, df7; analytical normalized density',
                nonlinear_map='z=c_z+L*t_z+B*(Phi(eta)-Phi(c_eta)); eta=c_eta+s*t_eta',
                coordinate_initialization_inward_only=bool(np.any(q0[64:]<=0)|np.any(q0[64:]>=1)),
                posterior_target='-Phi(FullMaxwell(from_probit(x)), online_y)-0.5*x@x; constants omitted',
                target_dimension=66,adaptation_after_freeze=False)
            freeze(dest/'PROTOCOL.json',protocol)
            save(dest/'FROZEN_PROPOSAL.npz',centers=centers,field_cholesky=L,property_scales=scales,compensation=B,
                 weights=proposal.weights,df=7.,J=J,x0=x0,native_anchor_q=q0)
            def logtarget(x):
                phi,key=evaluate(from_probit(x),'p3r2_full_nonlinear_chart_mh')
                return -phi-.5*float(x@x)
            starts=centers[np.arange(args.chains)%len(centers)].copy()
        else:
            protocol.update(fixed_z=q0[:64],bounds=[[0.,1.],[0.,1.]],
                prior='Original independent Uniform[0,1]^2, conditioned on field z fixed at anchor',
                quadrature='Adaptive global tensor Gauss-Legendre order doubling; no surrogate likelihood',
                stopping='Two successive material-mean RMS and absolute log-integral differences <= frozen tolerance',
                posterior_target='-Phi(FullMaxwell(q0_z,u), online_y); constant fixed-z Gaussian prior omitted')
            freeze(dest/'PROTOCOL.json',protocol)
            def slice_value(u):
                if np.any(u<=0) or np.any(u>=1): return -np.inf
                q=q0.copy();q[64:]=u;phi,key=evaluate(q,'p3r2_SLICE_U2_quadrature_or_mh')
                return -phi
            def quad_value(u):
                lp=slice_value(u);q=q0.copy();q[64:]=u;return lp,target.prior.decode(q)[0]
            previous=None;stable=0
            for order in orders:
                reference=quadrature_moments(quad_value,order)
                save(dest/f'QUADRATURE_{order}.npz',**{k:v for k,v in reference.items() if isinstance(v,(np.ndarray,int,float))})
                rms=None if previous is None else float(np.linalg.norm(reference['mean']-previous['mean'])/24)
                dz=None if previous is None else float(abs(reference['log_integral']-previous['log_integral']))
                stable=stable+1 if previous is not None and rms<=args.refinement_rms and dz<=args.refinement_rms else 0
                rows.append(dict(order=order,new_states=order*order,log_integral=reference['log_integral'],
                    material_mean_successive_RMS=rms,log_integral_successive_abs_difference=dz,stable_refinements=stable))
                write(dest/'QUADRATURE_PROGRESS.json',dict(rows=rows,status='NUMERICAL_REFINEMENT_ONLY',new_RHS=calls*8))
                print(json.dumps(finite(rows[-1])),flush=True);previous=reference
                if stable>=2: break
            proposal=mixture_from_quadrature(reference)
            save(dest/'FROZEN_PROPOSAL.npz',centers=proposal.centers,cholesky=proposal.cholesky,weights=proposal.weights,df=proposal.df)
            write(dest/'PROPOSAL_FREEZE.json',dict(source_order=reference['order'],freeze_before_MH=True,
                no_production_adaptation=True,quadrature_qualified=False,stopping_rule_met=stable>=2,
                proposal='Non-Gaussian nine-center normalized Student-t mixture; outside native box rejected before physics'))
            logtarget=slice_value;starts=proposal.centers[np.arange(args.chains)%len(proposal.centers)]
        transcript=(dest/'TRANSITIONS.jsonl').open('x')
        def callback(row):
            transcript.write(json.dumps(finite(row),allow_nan=False,separators=(',',':'))+'\n');transcript.flush()
        result=independent_mh(logtarget,proposal,starts,rng=rng,draws=args.draws,callback=callback)
        save(dest/'CHAINS.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        chi,diagnostics,moments=summarize_chains(result['states'],target.prior,
            slice_anchor=q0 if args.mode=='slice' else None,burn=args.burn)
        save(dest/'MATERIAL_SAMPLES.npz',chi=chi,**moments)
        save(dest/'MATERIAL_DIAGNOSTICS.npz',**{k:v for k,v in diagnostics.items() if isinstance(v,np.ndarray)})
        rhat=diagnostics['rhat'];ess=diagnostics['mean_ess'];mcse=diagnostics['mean_mcse']
        summary=dict(protocol=protocol,status='PILOT_COMPLETE_UNQUALIFIED',new_states=calls,new_RHS=calls*8,
            acceptance_per_chain=result['accepted'].mean(axis=1),quadrature_rows=rows,
            material_max_finite_rhat=float(np.max(rhat[np.isfinite(rhat)])) if np.isfinite(rhat).any() else None,
            material_min_finite_mean_ESS=float(np.min(ess[np.isfinite(ess)])) if np.isfinite(ess).any() else None,
            undefined_rhat_functions=int(np.sum(~np.isfinite(rhat))),undefined_ESS_functions=int(np.sum(~np.isfinite(ess))),
            material_mean_RMS_MCSE=float(np.sqrt(np.mean(np.sum(mcse**2,axis=-1)))),
            MCSE_scope='Complex material RMS combines real and imaginary component variances; asymptotic diagnostic only',
            posterior_qualification=False,reference_qualification=False,claim='Finite pilot diagnostics only; no FULL calibrated posterior claim')
        if args.mode=='slice':
            mean=moments['material_mean'];variance=moments['material_variance'];production=result['states'][:,args.burn:]
            summary.update(material_mean_RMS_error_vs_quadrature=float(np.linalg.norm(mean-reference['mean'])/24),
                material_variance_RMS_error_vs_quadrature=float(np.linalg.norm(variance-reference['variance'])/24),
                property_mean=production.mean(axis=(0,1)),quadrature_property_mean=reference['coordinate_mean'],
                property_mean_error=production.mean(axis=(0,1))-reference['coordinate_mean'],
                property_mean_ESS=material_diagnostics(production)['mean_ess'],
                property_rhat=material_diagnostics(production)['rhat'],
                property_mean_MCSE=material_diagnostics(production)['mean_mcse'],
                quadrature_stopping_rule_met=stable>=2,
                comparison_scope='Independent exact MH vs deterministic Full-likelihood quadrature on same SLICE only')
        write(dest/'RESULT.json',finite(summary));print(json.dumps(finite(summary)),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(dest/'ERROR.json',dict(error=error,new_states=calls,new_RHS=calls*8,partial_endpoints_preserved=True));raise
    finally:
        if transcript is not None: transcript.close()
        book.close(status,error)


if __name__=='__main__': main()
