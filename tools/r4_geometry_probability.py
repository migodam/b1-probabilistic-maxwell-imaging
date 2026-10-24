#!/usr/bin/env python3
"""Prepare frozen full-rank geometry pCN kernels, or explicitly released runs.

No GT, new J, geometry/teacher recomputation, context generation, or training.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,sys,time,resource,gc,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.frozen_geometry_pcn import *
from p3r4.evaluation import checkpoint_identity,write_new_json
STAGE=ROOT/'research/p3l_r4_conditional_transport';FIRST=STAGE/'results/III/dynamic_geometry_first'
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(8));STEPS=1024;CAP=262144;EXPERIMENT='frozen_geometry_probability_v1';SEED=20261025111
TEACHER_RHS=459712


def finite(value):
    if isinstance(value,dict):return {k:finite(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [finite(v) for v in value]
    if isinstance(value,np.ndarray):return finite(value.tolist())
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    if isinstance(value,np.generic):return value.item()
    return value


def atomic_npz(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix('.writing')
    with temp.open('xb') as f:np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    temp.replace(path)


def prepare(args):
    if args.out.exists():raise FileExistsError('new frozen preparation directory only')
    needed={parent:dict(frame=FIRST/parent/'d66_r4/BASE_FRAME.npz',bracket=FIRST/parent/'d66_r4/BRACKET_h0.001.npz',
                       cache=STAGE/'results/contexts/cache'/f'{parent}.npz',observation=STAGE/'inputs/online'/f'{parent}.npz') for parent in PARENTS}
    missing=[str(p) for row in needed.values() for p in row.values() if not p.exists()]
    if missing:print(json.dumps(dict(status='MISSING_SAVED_TEACHER_NO_RECOMPUTATION',paths=missing,new_RHS=0)));raise SystemExit(3)
    from p1_material.prior import MaterialPrior
    cfg=json.loads((STAGE/'configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];prior=MaterialPrior(cfg,points)
    first_receipt=STAGE/'III/jobs/dynamic_geometry_first/CLOSED.json'
    if not first_receipt.exists() or json.loads(first_receipt.read_text())['status']!='COMPLETE':raise ValueError('saved teacher owner must be closedCOMPLETE')
    args.out.mkdir(parents=True);rows=[]
    for index,parent in enumerate(PARENTS):
        sources=needed[parent]
        with np.load(sources['frame'],allow_pickle=False) as z:
            x0=z['x'];q0=z['q'];V=z['Vw'];Dchi=z['Dchi'];rank=int(z['rank']);savedJshape=z['J'].shape
        with np.load(sources['bracket'],allow_pickle=False) as z:b=z['bracket']
        with np.load(sources['cache'],allow_pickle=False) as z:cached_q=z['q']
        if rank!=4 or x0.shape!=(66,) or q0.shape!=(66,) or savedJshape!=(512,66) or not np.array_equal(q0,cached_q):raise ValueError('saved full66/rank4 online anchor teacher mismatch')
        from scipy.special import ndtr
        if not np.allclose(np.r_[x0[:64],ndtr(x0[64:])],q0,rtol=0,atol=1e-14):raise ValueError('teacher is not exact Gaussianized full66')
        source_D=prior.derivative(q0)*np.r_[np.ones(64),np.exp(-.5*x0[64:]**2)/np.sqrt(2*np.pi)][None,:]
        Derror=float(np.max(np.abs(source_D-Dchi)))
        if not np.allclose(source_D,Dchi,rtol=1e-10,atol=1e-11):raise ValueError('saved Dchi is not original Gaussianized decoder derivative; no recomputation of teacher')
        family,enrichment=freeze_kernel_family(V,b,Dchi,x0,np.random.SeedSequence([SEED,index,77]))
        folder=args.out/parent;folder.mkdir();atomic_npz(folder/'ANCHOR_TEACHER.npz',x0=x0,q0=q0,Dchi=Dchi,Vw=V,bracket=b,random_direction=enrichment['random_direction'])
        matrices={}
        for name,item in family.items():
            destination=folder/(name+'.npz');atomic_npz(destination,Q=item['Q'],rho=item['rho'],sqrtQ=item['sqrtQ']);matrices[name]=dict(file=str(destination.resolve()),identity=checkpoint_identity(destination),**item['metadata'])
        rows.append(dict(parent=parent,source_Dchi_validation_max_error=Derror,sources={k:checkpoint_identity(p) for k,p in sources.items()},kernels=matrices,
                         teacher_file=checkpoint_identity(folder/'ANCHOR_TEACHER.npz'),bracket_normalization='unitEuclidean,positive rank1relativeweight1',random_seed=[SEED,index,77]))
    protocol=dict(status='FROZEN_GEOMETRY_PROBABILITY_PREPARED',experiment=EXPERIMENT,parents=PARENTS,kernels=KERNELS,steps=STEPS,seed=SEED,
        epsilon_isotropic_weight=EPSILON,enrichment_relative_weight=1.,material_linear_RMS_target=RMS_TARGET,Q_eigenvalue_cap=Q_CAP,
        matching='anchor expected linear material RMS includes drift squared; ifcapbinds do notclaimmatchedRMS',
        covariance='frozen peronlineanchor; SPD full66; rho=sqrt(I-Q); prior-reversible exactGaussian proposal',
        likelihood_only_MH='logaccept=min(0,Phi_current-Phi_candidate)',source_teacher='already-paid first8DEV originalrank4/full66probit geometry pluscentralbracketh.001',
        earlier_teacher_setup_RHS=TEACHER_RHS,earlier_teacher_setup_owner='sealed StageIII geometry campaigns; disclosed sunkresearchcost, no rebilling',
        new_J_RHS=0,new_teacher_RHS=0,new_full_RHS_cap=CAP,new_full_LU_planned=8*4*STEPS,
        common_random_stream='each kernel resets SeedSequence([seed,parent_index]); every step normal66 thenuniform, same innovations/U and initialanchor',
        no_GT=True,posterior_qualification='UNQUALIFIED',common_reference='NONE',dynamic_bracket_nonzero_is_not_benefit=True,rows=rows)
    write_new_json(args.out/'FROZEN_PROTOCOL.json',finite(protocol));print(json.dumps(dict(status=protocol['status'],parents=8,kernels=4,newphysics=0,cap_bindings=sum(k['cap_binds'] for r in rows for k in r['kernels'].values()))))


def functions(chi,base,roi,projections,fraction):
    real=np.r_[chi.real,chi.imag]/24
    return np.r_[chi[roi].real.mean(),chi[roi].imag.mean(),np.mean(fraction[roi]>=.5),np.sqrt(np.mean(abs(chi-base)**2)),projections@real]


def run(args):
    if args.release is None:raise ValueError('parentrelease required; no automatic paidrun')
    release=json.loads(args.release.read_text());frozen=json.loads((args.prepared/'FROZEN_PROTOCOL.json').read_text())
    if release.get('released') is not True or release.get('experiment')!=EXPERIMENT or release.get('campaign_id')!='P3LR4-20261025':raise PermissionError('exact parentrelease required')
    if frozen['status']!='FROZEN_GEOMETRY_PROBABILITY_PREPARED' or tuple(frozen['parents'])!=PARENTS:raise ValueError('frozen first8DEV protocol')
    if not args.max_wall or args.max_wall<=15:raise ValueError('explicit bounded maxwall required')
    if args.out.exists():raise FileExistsError('fresh paidout only')
    for row in frozen['rows']:
        for identity in list(row['sources'].values())+[row['teacher_file']]+[k['identity'] for k in row['kernels'].values()]:
            if checkpoint_identity(identity['path'])!=identity:raise ValueError('frozen source/kernel changed; no recomputation')
    from p3r4.core import Book,Target,BudgetStop,save,write,from_probit
    book=Book('III',EXPERIMENT,max_wall=args.max_wall);original=book._update;costs=dict(RHS=0,states=0,LU=0);status='COMPLETE';error=None;started=time.monotonic();summaries=[]
    def limited(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>CAP:raise BudgetStop('262144new fullsourceRHS cap')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiB physicallimit')
        answer=original(event,**kw)
        for k,a in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(a,0)
        return answer
    book._update=limited
    try:
        args.out.mkdir(parents=True);write(args.out/'FROZEN_PROTOCOL.json',frozen);target=Target(book,keep_fields=True);roi=np.flatnonzero(target.prior.mask)
        pr=np.random.default_rng(SEED+99).normal(size=(4,1152));pr/=np.linalg.norm(pr,axis=1)[:,None]
        save(args.out/'MATERIAL_FUNCTIONS.npz',projections=pr,names=np.array(['roiReMean','roiImMean','fraction_ge_half','RMS_from_anchor']+[f'fixed_projection{i}' for i in range(4)]))
        for parent_index,row in enumerate(frozen['rows']):
            parent=row['parent'];base=target.load_cache(row['sources']['cache']['path'],owner='paid StageIonline anchor plus StageIIIframe/bracket owner459712RHS')
            with np.load(row['teacher_file']['path']) as z:x0=z['x0']
            y=target.observation(parent);base_phi=target.phi(base,y)
            for name in KERNELS:
                folder=args.out/parent/name;folder.mkdir(parents=True);initial=costs.copy();tick=time.monotonic()
                with np.load(row['kernels'][name]['file']) as z:kernel=FrozenMatrixPCN(z['Q'])
                rng=np.random.default_rng(np.random.SeedSequence([SEED,parent_index]));x=x0.copy();q=base.q.copy();chi=base.state.chi.copy();current=base;phi=base_phi
                xs=[x.copy()];material=[];phis=[];priors=[];accepts=[];moves=[];proposal_moves=[]
                for step in range(STEPS):
                    proposed=kernel.propose(x,rng.normal(size=66));logu=float(np.log(rng.uniform()));qnew=from_probit(proposed)
                    book.context=dict(parent=parent,kernel=name,step=step,purpose='fixed_fullrank_priorreversible_geometry_probability')
                    before=costs['RHS'];candidate=target.evaluate(qnew,purpose=f'{parent}_{name}_s{step:04d}');candidate_phi=target.phi(candidate,y)
                    if costs['RHS']-before!=8:raise AssertionError('all8complexsourceRHS perFull')
                    forward=kernel.log_proposal(proposed,x);reverse=kernel.log_proposal(x,proposed)
                    algebra=-candidate_phi+phi-.5*(proposed@proposed-x@x)+reverse-forward
                    simple=phi-candidate_phi
                    if abs(algebra-simple)>1e-7:raise AssertionError('fullGaussian proposal ratio failed prior cancellation')
                    accepted=logu<min(0.,simple);previous_chi=chi.copy();proposal_move=float(np.sqrt(np.mean(abs(candidate.state.chi-chi)**2)))
                    # Common immediate candidate/state persistence includes all
                    # rejections; no residual-based candidate filtering.
                    save(folder/'candidates'/f's{step:04d}.npz',x=proposed,q=qnew,chi=candidate.state.chi,current=candidate.state.current,field=candidate.state.field,prediction=candidate.state.data,Phi=np.array(candidate_phi),prior_energy=np.array(.5*(proposed@proposed)),logu=np.array(logu),accepted=np.array(accepted),proposal_material_RMS=np.array(proposal_move))
                    if accepted:x=proposed;q=qnew;chi=candidate.state.chi.copy();current=candidate;phi=candidate_phi
                    move=float(np.sqrt(np.mean(abs(chi-previous_chi)**2)));fraction=target.prior.decode(q)[1]['fraction'];fn=functions(chi,base.state.chi,roi,pr,fraction)
                    save(folder/'states'/f's{step:04d}.npz',x=x,q=q,chi=chi,current=current.state.current,field=current.state.field,prediction=current.state.data,Phi=np.array(phi),prior_energy=np.array(.5*(x@x)),accepted=np.array(accepted),material_functions=fn,material_move_RMS=np.array(move),endpoint=np.array(current.key))
                    xs.append(x.copy());material.append(fn);phis.append(phi);priors.append(.5*(x@x));accepts.append(accepted);moves.append(move);proposal_moves.append(proposal_move)
                save(folder/'CHAIN.npz',x_path=np.array(xs),material_functions=np.array(material),Phi=np.array(phis),prior_energy=np.array(priors),accept=np.array(accepts),material_move_RMS=np.array(moves),proposal_material_RMS=np.array(proposal_moves))
                summary=dict(parent=parent,kernel=name,accepts=int(sum(accepts)),steps=STEPS,costs={k:costs[k]-initial[k] for k in costs},wall_seconds=time.monotonic()-tick,
                    mean_material_move_RMS=float(np.mean(moves)),kernel_normalization=row['kernels'][name],material_diagnostics=finite_run_diagnostics(material),common_reference='NONE',posterior_qualification='UNQUALIFIED')
                write(folder/'SUMMARY.json',finite(summary));summaries.append(summary);write(args.out/'PROGRESS.json',finite(dict(completed=len(summaries),costs=costs,last=summary)))
                print(json.dumps(finite(summary)),flush=True);del candidate,current;target.backend()._transfers.clear();gc.collect()
            del base;gc.collect()
        write(args.out/'RESULTS.json',finite(dict(status='COMPLETE',costs=costs,earlier_teacher_setup_RHS=TEACHER_RHS,teacher_owner=frozen['earlier_teacher_setup_owner'],new_J_RHS=0,parents=summaries,
             wall_seconds=time.monotonic()-started,posterior_qualification='UNQUALIFIED',common_reference='NONE',bracket_nonzero_implies_benefit=False,no_GT=True)))
    except BaseException as exc:
        status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc);write(args.out/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:
        book.close(status,error);write(args.out/'FINAL_JOB_COST.json',dict(status=status,error=error,costs=costs,receipt=str(book.jobdir/'CLOSED.json')))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('prepare','run'));p.add_argument('--out',type=Path,required=True);p.add_argument('--prepared',type=Path);p.add_argument('--release',type=Path);p.add_argument('--max-wall',type=float);a=p.parse_args()
    if a.command=='prepare':prepare(a)
    else:run(a)

if __name__=='__main__':main()
