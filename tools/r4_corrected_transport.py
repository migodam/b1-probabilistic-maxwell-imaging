#!/usr/bin/env python3
"""Freeze/run explicitly corrected full66 q, pCN, and constant-mixture chains.

No GT, locked data, residual selection, projected/corrected proposals, CFM,
training, new Jacobian, or automatic follow-up pilot/production expansion.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,sys,time,resource,gc,traceback,re
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import ndtr
from p3r4.corrected_transport import *
from p3r4.evaluation import checkpoint_identity,write_new_json
from p3r4.probability_eval import source_major_pack
STAGE=ROOT/'research/p3l_r4_conditional_transport';EXPERIMENT='corrected_transport_full66_v1'


def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v


def get_model(path,loader):
    if loader=='raw':
        from p3r4.density import load_density
        model=load_density(path,'cpu')
    elif loader=='anchor':
        from p3r4.anchor_density import load_anchor_density
        model=load_anchor_density(path,'cpu')
    else:raise ValueError('raw/anchor explicitq loader only')
    if model.payload['kind'] not in ('spline','diag_gaussian') or model.model.dimension!=66:raise ValueError('explicitSpline/diagGaussian all66 only; CFM/numerical-density disallowed')
    if loader=='raw' and model.model.context_dim!=512:raise ValueError('raw loader must be rawy512; anchor checkpoints need explicit anchor loader')
    return model


def condition_from_export(file,loader):
    with np.load(file,allow_pickle=False) as z:
        observation=z['observation'];x0=z['context_x'];J=z['context_J'] if 'context_J' in z else None
    condition={} if loader=='raw' else dict(context_x=x0[None])
    if loader=='anchor' and J is not None:condition['context_J']=J[None]
    return source_major_pack(observation),condition,x0,observation


def sample_q(model,y,condition,count,seed):
    value=model.sample(y,n_samples=count,seed=seed,batch_size=16,**condition)[0]
    if value.shape!=(count,66) or not np.isfinite(value).all():raise ValueError('explicitq must return exactlyall66 finite independentdraws')
    return value


def density_q(model,x,y,condition):
    value=float(model.log_prob(np.asarray(x)[None],y,batch_size=1,**condition)[0])
    if not np.isfinite(value):raise ValueError('finite explicit rawdx66 logq required; no retry/repair')
    return value


def freeze(a):
    spec=json.loads(a.spec.read_text());parents=spec['parents']
    if len(parents) not in (4,8) or len(set(parents))!=len(parents) or any(not re.fullmatch(r'R4DEV_[0-9]{5}',p) for p in parents):raise ValueError('fixed4or8DEVparents,no locked')
    warm=int(spec.get('warmup',0));prod=int(spec.get('production',128));seed=int(spec.get('seed',20261025141));loader=spec.get('loader','raw')
    if warm<0 or prod<2 or warm>100000 or prod>100000:raise ValueError('bounded predeclared warmup/production, production>=2')
    if a.out.exists():raise FileExistsError('new frozen experiment only')
    model=get_model(spec['checkpoint'],loader);contract=model.payload.get('anchor_contract',{})
    rows=[];cfg=json.loads((STAGE/'configs/TARGET.json').read_text());a.out.mkdir(parents=True)
    for parent in parents:
        obs=STAGE/'inputs/online'/f'{parent}.npz';context=STAGE/'results/contexts/dev'/f'{parent}.npz';side=context.with_suffix('.json')
        online_metadata=json.loads(side.read_text())
        if online_metadata.get('GT_read') is not False or online_metadata.get('parent')!=parent:raise ValueError('actual noGT online context sidecar required')
        with np.load(obs,allow_pickle=False) as z:
            y=z['observation'];target_id=str(z['target_id']);sigma=float(z['sigma'])
        if target_id!=cfg['target_id'] or sigma!=cfg['noise']['sigma']:raise ValueError('original acquisition/noise target')
        source_major_pack(y)
        # Ambiguousx is accepted ONLY from the canonical StageI truth-free dev
        # context path just verified, then exported under explicitcontext_x.
        with np.load(context,allow_pickle=False) as z:x0=z['context_x'] if 'context_x' in z else z['x'];q0=z['q0'];J=z['context_J'] if 'context_J' in z else z['J']
        if x0.shape!=(66,) or q0.shape!=(66,) or not np.isfinite(x0).all() or np.max(abs(np.r_[x0[:64],ndtr(x0[64:])]-q0))>1e-12:raise ValueError('canonicalonlineexactprobitanchor')
        needJ=loader=='anchor' and (contract.get('j_summary','none')!='none' or contract.get('coordinate_preconditioner','none')!='none')
        exported=dict(parent=np.array(parent),observation=y,context_x=x0)
        if needJ:
            if J.shape!=(512,66) or not np.isfinite(J).all():raise ValueError('saved nativeonlineJ required, no newJ')
            exported['context_J']=J
        destination=a.out/'conditions'/f'{parent}.npz';destination.parent.mkdir(exist_ok=True)
        with destination.open('xb') as f:np.savez_compressed(f,**exported)
        rows.append(dict(parent=parent,export=checkpoint_identity(destination),observation=checkpoint_identity(obs),online_context=checkpoint_identity(context),online_sidecar=checkpoint_identity(side),q_condition='all rawy512 only' if loader=='raw' else 'rawy512+frozenobservation-derivedcontext_x/Jsummaries, never currentstate'))
    # Only validate the explicit q density interface; no source samples selected
    # by likelihood, no Full evaluation, and no q starts are generated here.
    protocol=dict(status='FROZEN_CORRECTED_TRANSPORT',experiment=EXPERIMENT,checkpoint=checkpoint_identity(spec['checkpoint']),loader=loader,kind=model.payload['kind'],
      parents=parents,rows=rows,chains=4,kernels=KERNELS,warmup=warm,production=prod,seed=seed,beta=BETA,mixture_q_probability=.5,
      initialization='same onlineanchor plus3 independentqdraws common acrosskernels; all4starts separatelyFullpaid ineachkernel',
      planned_Full_evaluations=len(parents)*3*4*(1+warm+prod),planned_RHS=8*len(parents)*3*4*(1+warm+prod),
      context_cost_owner='original paid online anchor/J and frozen q training checkpoint; not cold-free, no newJ',q_density='explicit algebraic full66 rawdx density, model transformsfloat32 plus savednormalizer/coordinate Jacobians; not machine-exact MH certification',
      MH='component-specific individually correctedqorpriorpCN; constantstateindependentselector, never wrongmixture-density ratio',
      no_GT=True,new_J_RHS=0,no_projection=True,no_corrector=True,no_residual_selection=True,no_autoexpand=True,posterior_qualification='NOT_ESTABLISHED',target_config=cfg)
    write_new_json(a.out/'FROZEN.json',protocol);print(json.dumps(dict(status=protocol['status'],parents=len(parents),planned_RHS=protocol['planned_RHS'],physics_RHS=0)))


def run(a):
    frozen=json.loads((a.prepared/'FROZEN.json').read_text());release=json.loads(a.release.read_text())
    if frozen['status']!='FROZEN_CORRECTED_TRANSPORT' or release.get('released') is not True or release.get('experiment')!=EXPERIMENT or release.get('campaign_id')!='P3LR4-20261025':raise ValueError('exact frozen scope+parentrelease required')
    if a.max_rhs is None or a.max_rhs<frozen['planned_RHS'] or a.max_wall is None or a.max_wall<=15 or not a.job or not a.stage:raise ValueError('explicitjob/stage/hardRHS/wall must cover entire frozenpilot')
    if a.out.exists():raise FileExistsError('no paidreplay')
    if checkpoint_identity(frozen['checkpoint']['path'])!=frozen['checkpoint']:raise ValueError('frozenqcheckpointchanged')
    for row in frozen['rows']:
        for k in ('export','observation','online_context','online_sidecar'):
            if checkpoint_identity(row[k]['path'])!=row[k]:raise ValueError('onlinecondition changed')
    from p3r4.core import Book,Target,BudgetStop,save,write,from_probit
    import torch
    torch.set_num_threads(1);book=Book(a.stage,a.job,max_wall=a.max_wall);original=book._update;costs=dict(RHS=0,states=0,LU=0);status='COMPLETE';error=None;start=time.monotonic();summaries=[]
    model_counts=dict(sample_calls=0,sample_draws=0,logq_calls=0,model_wall_seconds=0.)
    def limited(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>a.max_rhs:raise BudgetStop('correctedtransport hardRHS cap before everyattempt')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiB cap')
        value=original(event,**kw)
        for k,v in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(v,0)
        return value
    book._update=limited
    try:
        a.out.mkdir(parents=True);write(a.out/'FROZEN.json',frozen);target=Target(book,keep_fields=True);model=get_model(frozen['checkpoint']['path'],frozen['loader']);roi=np.flatnonzero(target.prior.mask)
        if target.cfg!=frozen['target_config']:raise ValueError('frozen statistical target changed')
        projection=np.random.default_rng(frozen['seed']+99).normal(size=(4,1152));projection/=np.linalg.norm(projection,axis=1)[:,None]
        save(a.out/'MATERIAL_FUNCTIONS.npz',projections=projection,names=np.array(['roiRe','roiIm','fraction_ge_half','RMS_from_anchor']+[f'projection{i}' for i in range(4)]+['Phi','prior_energy']))
        def draw(y,condition,count,seed):
            tick=time.perf_counter();v=sample_q(model,y,condition,count,seed);model_counts['sample_calls']+=1;model_counts['sample_draws']+=count;model_counts['model_wall_seconds']+=time.perf_counter()-tick;return v
        def logq(x,y,condition):
            tick=time.perf_counter();v=density_q(model,x,y,condition);model_counts['logq_calls']+=1;model_counts['model_wall_seconds']+=time.perf_counter()-tick;return v
        for parent_index,row in enumerate(frozen['rows']):
            parent=row['parent']
            if checkpoint_identity(row['export']['path'])!=row['export']:raise ValueError('frozencondition changed duringjob')
            y,condition,x0,original_observation=condition_from_export(row['export']['path'],frozen['loader'])
            starts=np.vstack((x0,draw(y,condition,3,seed_for(frozen['seed'],parent_index,0,0,90))))
            base_chi=target.prior.decode(from_probit(x0))[0];save(a.out/parent/'INITIAL_X_ALL4.npz',x=starts)
            for kernel in KERNELS:
                folder=a.out/parent/kernel;folder.mkdir(parents=True);before=costs.copy();tick=time.monotonic();production_functions=[];parent_chains=[]
                for chain in range(4):
                    chainfolder=folder/f'chain{chain}';chainfolder.mkdir();x=starts[chain].copy();q=from_probit(x)
                    if costs['RHS']+8>a.max_rhs:raise BudgetStop('initialFullpreflight8RHS')
                    book.context=dict(parent=parent,kernel=kernel,chain=chain,phase='independent_initialization')
                    current_logq=None if kernel=='pcn' else logq(x,y,condition)
                    initial_before=costs['RHS'];current=target.evaluate(q,purpose=f'{parent}_{kernel}_chain{chain}_initial');phi=target.phi(current,original_observation)
                    if costs['RHS']-initial_before!=8:raise AssertionError('allinitialstates includingonlineanchor charge8FullRHS')
                    save(chainfolder/'INITIAL_REFERENCE.npz',x=x,q=q,endpoint=np.array(current.key),Phi=np.array(phi))
                    xs=[x.copy()];flags=[];refs=[current.key];phis=[];prior=[];functions=[];component_history=[];draw_count=0
                    for step in range(frozen['warmup']+frozen['production']):
                        innovation,logu,selector,qseed=step_randomness(frozen['seed'],parent_index,chain,step);component=component_for(kernel,selector)
                        candidate=pcn_proposal(x,innovation) if component=='pcn' else draw(y,condition,1,qseed)[0]
                        qnew=from_probit(candidate);save(chainfolder/'attempts'/f's{step:06d}.npz',x=candidate,q=qnew,logu=np.array(logu),selector=np.array(selector),component=np.array(component),q_seed=np.array(qseed))
                        proposed_logq=logq(candidate,y,condition) if component=='q' else None
                        if costs['RHS']+8>a.max_rhs:raise BudgetStop('candidateFullpreflight8RHS')
                        book.context=dict(parent=parent,kernel=kernel,chain=chain,step=step,component=component)
                        previous=costs['RHS'];point=target.evaluate(qnew,purpose=f'{parent}_{kernel}_c{chain}_s{step:06d}')
                        if costs['RHS']-previous!=8:raise AssertionError('everyFull8complexsourceRHS')
                        newphi=target.phi(point,original_observation)
                        logalpha=log_acceptance(component,x,phi,candidate,newphi,current_logq,proposed_logq)
                        accepted=logu<logalpha
                        if accepted:
                            x=candidate;q=qnew;current=point;phi=newphi
                            if kernel!='pcn':current_logq=proposed_logq if component=='q' else logq(x,y,condition)
                        # Rejections keep exactsamex/ref/current cachedlogq.
                        chi=current.state.chi;fraction=target.prior.decode(q)[1]['fraction'];real=np.r_[chi.real,chi.imag]/24
                        fn=np.r_[chi[roi].real.mean(),chi[roi].imag.mean(),np.mean(fraction[roi]>=.5),np.sqrt(np.mean(abs(chi-base_chi)**2)),projection@real,phi,.5*(x@x)]
                        append=dict(step=step,production=step>=frozen['warmup'],component=component,selector=selector,logu=logu,logalpha=logalpha,
                            accepted=accepted,candidate_endpoint=point.key,current_endpoint=current.key,current_x=x,FullPhi=phi,candidatePhi=newphi,
                            current_logq=current_logq,proposed_logq=proposed_logq,precision='algebraicfloat32modeldensity, not machine-exact certificate')
                        with (chainfolder/'HISTORY.jsonl').open('a') as f:f.write(json.dumps(clean(append),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
                        xs.append(x.copy());flags.append(accepted);refs.append(current.key);phis.append(phi);prior.append(.5*(x@x));functions.append(fn);component_history.append(component)
                        del point
                    array=np.array(functions);production_functions.append(array[frozen['warmup']:])
                    save(chainfolder/'CHAIN.npz',x_path=np.array(xs),accept=np.array(flags),endpoint_reference=np.array(refs),Phi=np.array(phis),prior_energy=np.array(prior),material_functions=array,component=np.array(component_history),is_production=np.arange(len(flags))>=frozen['warmup'])
                    parent_chains.append(dict(chain=chain,accepts=int(sum(flags)),steps=len(flags),q_components=component_history.count('q'),pcn_components=component_history.count('pcn')))
                    del current;target.backend()._transfers.clear();gc.collect()
                diagnostics=multichain_diagnostics(np.array(production_functions))
                summary=dict(parent=parent,kernel=kernel,chains=parent_chains,costs={k:costs[k]-before[k] for k in costs},wall_seconds=time.monotonic()-tick,
                    diagnostics=diagnostics,full_reference='NONE',posterior_qualification='UNQUALIFIED')
                summaries.append(summary);write(folder/'SUMMARY.json',clean(summary));write(a.out/'PROGRESS.json',clean(dict(completed=len(summaries),costs=costs,last=summary)));print(json.dumps(clean(summary)),flush=True)
        if checkpoint_identity(frozen['checkpoint']['path'])!=frozen['checkpoint']:raise ValueError('frozenqcheckpoint changed duringjob')
        write(a.out/'RESULTS.json',clean(dict(status='COMPLETE',parents=summaries,costs=costs,model_costs=model_counts,wall_seconds=time.monotonic()-start,
            GT_read=False,new_J_RHS=0,no_autoexpand=True,kernel_correctness='formalMHratiosforfrozenfull66q; float32explicitdensityprecisionlimitation',posterior_qualification='UNQUALIFIED')))
    except BaseException as exc:
        status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc);write(a.out/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:book.close(status,error)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('freeze','run'))
    for key in ('spec','prepared','release','out'):p.add_argument('--'+key,type=Path,required=key=='out')
    p.add_argument('--stage',choices=('II','III','IV'));p.add_argument('--job');p.add_argument('--max-rhs',type=int);p.add_argument('--max-wall',type=float);a=p.parse_args()
    if a.command=='freeze':freeze(a)
    else:run(a)

if __name__=='__main__':main()
