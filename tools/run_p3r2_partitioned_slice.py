"""New registered full-box native-property SLICE with response-scale partitions.

Every patch uses FullMaxwell likelihood, with positive quadrature weights and
full [0,1]^2 coverage. Cached response derivatives only choose frozen cuts.
"""
from pathlib import Path
import argparse,json,os,sys
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import Book,Target,write,save
from p3r2.non_gaussian import tensor_rule_2d,mixture_from_quadrature,independent_mh
from run_p3r2_posterior import finite,freeze,summarize_chains
from p3l.diagnostics import material_diagnostics


def partition_cuts(anchor,sd):
    cuts=[]
    for center,scale in zip(anchor,sd):
        values=[0.,1.,float(center)]
        values.extend(float(center+sign*multiple*scale) for sign in (-1,1) for multiple in (2.,6.))
        cuts.append(np.array(sorted(set(v for v in values if 0<=v<=1))))
    return cuts


def summarize(nodes,qw,lp,chi,order):
    lz=float(logsumexp(lp+np.log(qw)));w=np.exp(lp+np.log(qw)-lz);mean=w@chi
    return dict(nodes=nodes,quadrature_weights=qw,logtarget=lp,weights=w,observables=chi,
        mean=mean,variance=np.einsum('i,ij->j',w,np.abs(chi-mean)**2),coordinate_mean=w@nodes,
        coordinate_covariance=(nodes-w@nodes).T@(w[:,None]*(nodes-w@nodes)),log_integral=lz,order=order)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job',required=True);parser.add_argument('--anchor',type=Path,required=True)
    parser.add_argument('--parent',default='DEV_030');parser.add_argument('--rhs-cap',type=int,default=262144)
    parser.add_argument('--orders',default='8,12,20');parser.add_argument('--draws',type=int,default=1024)
    parser.add_argument('--chains',type=int,default=4);parser.add_argument('--burn',type=int,default=256)
    parser.add_argument('--seed',type=int,default=2026102292);parser.add_argument('--max-wall',type=int,default=7200)
    parser.add_argument('--refinement-rms',type=float,default=1e-6)
    args=parser.parse_args();orders=tuple(int(x) for x in args.orders.split(','))
    if args.rhs_cap>262144 or args.chains<2 or args.draws<args.burn+4 or any(n<2 for n in orders) or any(b<=a for a,b in zip(orders,orders[1:])): raise ValueError('bounded independently registered pilot required')
    with np.load(args.anchor,allow_pickle=False) as z:
        if 'receiver_T' not in z.files or not z['receiver_T'].size: raise ValueError('already-paid anchor J transfer required')
    book=Book('B',args.job,max_wall=args.max_wall);status='COMPLETE';error=None;calls=0;transcript=None
    try:
        target=Target(book);y=target.observation(args.parent);anchor=target.load_cache(args.anchor,owner='P3L old paid DEV030 metric only')
        q0=anchor.q.copy();Ju=target.jacobian(anchor)[:,64:].copy();residual=target.white(anchor.state.data-y)
        H=Ju.T@Ju;H=(H+H.T)/2;C=np.linalg.inv(H);sd=np.sqrt(np.diag(C))
        localcenter=q0[64:]-np.linalg.solve(H,Ju.T@residual);cuts=partition_cuts(q0[64:],sd);del anchor
        patches=[((a,b),(c,d)) for a,b in zip(cuts[0][:-1],cuts[0][1:]) for c,d in zip(cuts[1][:-1],cuts[1][1:])]
        maximum_states=len(patches)*sum(n*n for n in orders)+args.chains*(args.draws+1)
        if maximum_states*8>args.rhs_cap: raise ValueError(f'frozen patch count needs {maximum_states*8} RHS > local cap')
        dest=book.jobdir
        protocol=dict(job=args.job,target_family='SLICE_U2_PARTITIONED_FULL_NATIVE_BOX_FIXED_Z64',original_FULL_posterior=False,
            target_id=target.cfg['target_id'],parent=args.parent,anchor=str(args.anchor.resolve()),q0=q0,
            cuts=cuts,patch_count=len(patches),orders=orders,fullbox_covered=True,
            property_response_H=H,local_property_response_SD=sd,local_linear_response_optimum=localcenter,
            cuts_rule='Boundaries 0,1; native anchor; anchor +/-2SD and +/-6SD clipped by discarding outside cuts',
            metric_source='Already-paid anchor receiver transfer; no new J calls; z fixed while computing SD',
            production_target='Original FullMaxwell likelihood times native Uniform2; original Gaussian64 fixed at anchor',
            adaptation='All cuts fixed before first new physics; only global quadrature order is refined',
            maximum_new_states=maximum_states,maximum_new_RHS=maximum_states*8,pilot_RHS_cap=args.rhs_cap,
            quadrature_stop='Two consecutive material-mean RMS and absolute log-integral differences <= frozen tolerance',
            tolerance=args.refinement_rms,GT_used=False,online_source_isolated=True,seed=args.seed,
            chains=args.chains,draws=args.draws,burn=args.burn,reference_qualification=False,
            caution='SLICE conditions on64 field coordinates; no FULL marginal or quadrature certification')
        freeze(dest/'PROTOCOL.json',protocol)
        save(dest/'PARTITION.npz',cut_real=cuts[0],cut_imag=cuts[1],Huu=H,conditional_SD=sd,q0=q0)
        def logtarget(u):
            nonlocal calls
            if np.any(u<=0) or np.any(u>=1): return -np.inf
            if (calls+1)*8>args.rhs_cap: raise RuntimeError('independent partitioned SLICE pilot cap')
            calls+=1;q=q0.copy();q[64:]=u;p=target.evaluate(q,purpose='p3r2_SLICE_partitioned_exactFull');phi=target.phi(p,y)
            del p
            if target.model is not None: target.model._transfers.clear()
            return -phi
        rows=[];previous=None;stable=0
        for order in orders:
            ns=[];ws=[];lps=[];materials=[]
            for patchindex,bounds in enumerate(patches):
                nodes,weights=tensor_rule_2d(order,bounds);ns.append(nodes);ws.append(weights)
                for u in nodes:
                    lp=logtarget(u);q=q0.copy();q[64:]=u;lps.append(lp);materials.append(target.prior.decode(q)[0])
                write(dest/'PATCH_PROGRESS.json',dict(order=order,patch_completed=patchindex+1,patch_count=len(patches),new_RHS=calls*8))
            reference=summarize(np.concatenate(ns),np.concatenate(ws),np.array(lps),np.array(materials),order)
            save(dest/f'QUADRATURE_{order}.npz',**reference)
            rms=None if previous is None else float(np.linalg.norm(reference['mean']-previous['mean'])/24)
            dz=None if previous is None else abs(reference['log_integral']-previous['log_integral'])
            stable=stable+1 if previous is not None and rms<=args.refinement_rms and dz<=args.refinement_rms else 0
            rows.append(dict(order=order,patches=len(patches),log_integral=reference['log_integral'],material_mean_successive_RMS=rms,log_integral_successive_abs_difference=dz,stable_refinements=stable))
            write(dest/'QUADRATURE_PROGRESS.json',dict(rows=rows,new_RHS=calls*8,stopping_rule_met=stable>=2,reference_qualification=False));print(json.dumps(finite(rows[-1])),flush=True)
            previous=reference
            if stable>=2: break
        proposal=mixture_from_quadrature(reference,components=9)
        save(dest/'FROZEN_PROPOSAL.npz',centers=proposal.centers,cholesky=proposal.cholesky,weights=proposal.weights,df=proposal.df)
        write(dest/'PROPOSAL_FREEZE.json',dict(source_order=reference['order'],frozen_before_MH=True,no_production_adaptation=True))
        transcript=(dest/'TRANSITIONS.jsonl').open('x')
        def callback(row): transcript.write(json.dumps(finite(row),allow_nan=False,separators=(',',':'))+'\n');transcript.flush()
        starts=proposal.centers[np.arange(args.chains)%len(proposal.centers)]
        result=independent_mh(logtarget,proposal,starts,rng=np.random.default_rng(args.seed),draws=args.draws,callback=callback)
        save(dest/'CHAINS.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        chi,diagnostics,moments=summarize_chains(result['states'],target.prior,slice_anchor=q0,burn=args.burn)
        save(dest/'MATERIAL_SAMPLES.npz',chi=chi,**moments);save(dest/'MATERIAL_DIAGNOSTICS.npz',**{k:v for k,v in diagnostics.items() if isinstance(v,np.ndarray)})
        production=result['states'][:,args.burn:];pd=material_diagnostics(production)
        summary=dict(protocol=protocol,status='PARTITIONED_SLICE_PILOT_COMPLETE_UNQUALIFIED',new_RHS=calls*8,new_states=calls,
            quadrature_rows=rows,quadrature_stopping_rule_met=stable>=2,acceptance_per_chain=result['accepted'].mean(axis=1),
            material_mean_RMS_error_vs_quadrature=float(np.linalg.norm(moments['material_mean']-reference['mean'])/24),
            material_variance_RMS_error_vs_quadrature=float(np.linalg.norm(moments['material_variance']-reference['variance'])/24),
            property_mean=production.mean(axis=(0,1)),quadrature_property_mean=reference['coordinate_mean'],
            property_mean_error=production.mean(axis=(0,1))-reference['coordinate_mean'],
            property_mean_ESS=pd['mean_ess'],property_rhat=pd['rhat'],property_mean_MCSE=pd['mean_mcse'],
            material_complex_RMS_mean_MCSE=float(np.sqrt(np.mean(np.sum(diagnostics['mean_mcse']**2,axis=-1)))),
            reference_qualification=False,posterior_qualification=False,
            interpretation='Independent exactMH versus Full-likelihood fullbox numerical quadrature on nativeU2 conditionalSLICE only')
        write(dest/'RESULT.json',finite(summary));print(json.dumps(finite(summary)),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(book.jobdir/'ERROR.json',dict(error=error,new_RHS=calls*8,paid_endpoints_preserved=True));raise
    finally:
        if transcript is not None: transcript.close()
        book.close(status,error)


if __name__=='__main__': main()
