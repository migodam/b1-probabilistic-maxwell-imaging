"""Bounded geometry-weighted anchor CFM; existing VI source is unchanged.

W_r is built from J_r=J_x T in the actual TRAIN-standardized network
coordinates. W_r=I+3 V diag(s²/(1+s²)) V.T has eigenvalues[1,4]. It weights
velocity ERROR, retaining every stochastic coordinate. No Maxwell action.
"""
from __future__ import annotations
import argparse,copy,json,os,resource,signal,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.linalg import qr,svd
import torch
from .anchor_density import FORMAT
from .conditional_coordinates import native_to_gaussianized_jacobian
from .density import make_model
from .train_anchor_density import load_anchor_data,prepare,cfm_error_loss
from .train_density import fit_normalization,write_json
from .train_physics_transport import identity,verify_identity

ROTATION_SEED=20261025191


def fixed_rotation(seed=ROTATION_SEED):
    r,u=qr(np.random.default_rng(seed).normal(size=(66,66)),mode='economic')
    return r*np.where(np.diag(u)>=0,1.,-1.)[None,:]


def network_metric(anchor_x,native_J,stats,rotation=None):
    """Only y-derived anchor/J and frozen global TRAIN statistics enter W."""
    jx=native_to_gaussianized_jacobian(np.asarray(anchor_x)[None],np.asarray(native_J)[None])[0]
    T=np.asarray(stats['frame_inverse'],np.float64).T*np.asarray(stats['x_scale'],np.float64)[None,:]
    jr=jx@T
    _,singular,vt=svd(jr,full_matrices=False,check_finite=False,lapack_driver='gesvd')
    shrink=(singular/np.hypot(singular,1.))**2
    W=np.eye(66)+3*(vt.T*shrink)@vt
    W=(W+W.T)/2
    if rotation is not None:W=rotation@W@rotation.T
    return (W+W.T)/2,shrink


def prepare_cache(train,validation,baseline,out):
    """Pure cached-array preprocessing; caller owns an analysis-only Book."""
    started,cpu=time.monotonic(),time.process_time()
    out=Path(out)
    if out.exists():raise FileExistsError('fresh frozen metric cache required')
    payload=torch.load(baseline,map_location='cpu',weights_only=False)
    ac=payload['anchor_contract'];stats=payload['normalization']
    if payload['format']!=FORMAT or payload['kind']!='cfm' or payload['model_config']['dimension']!=66 or ac['coordinate_preconditioner']!='none' or ac['j_summary']!='none' or not ac['residual_center']:
        raise ValueError('original anchor CFM_none full66/raw512+anchor66 baseline required')
    if not np.allclose(stats['frame'].T@stats['frame'],np.eye(66),atol=1e-10):raise ValueError('baseline frame must be orthogonal')
    rotation=fixed_rotation();out.mkdir(parents=True)
    np.save(out/'ROTATION.npy',rotation)
    np.savez(out/'NORMALIZATION.npz',**{k:v for k,v in stats.items() if isinstance(v,np.ndarray)})
    summaries={};allparents=[]
    for name,path,count in [('train',train,1024),('validation',validation,256)]:
        with np.load(path,allow_pickle=False) as z:
            # Never read x labels, material, decoder fields, or GT here.
            parents=z['parent'].astype(str);anchors=z['context_x'];J=z['context_J']
        if parents.shape!=(count,) or len(set(parents))!=count or anchors.shape!=(count,66) or J.shape!=(count,512,66):raise ValueError('fixed1024/256 full66 cached anchor/J shapes required')
        geometry=np.empty((count,66,66),np.float32);random=np.empty_like(geometry)
        shrinkages=np.empty((count,66),np.float64)
        eig_error=0.;lo,hi=4.,1.
        for i in range(count):
            W,shrink=network_metric(anchors[i],J[i],stats)
            control=rotation@W@rotation.T
            eig=np.linalg.eigvalsh(W)
            if eig.min()<1-1e-9 or eig.max()>4+1e-9:raise ValueError('unbounded metric')
            eig_error=max(eig_error,float(np.max(abs(np.linalg.eigvalsh(control)-eig))))
            lo=min(lo,float(eig.min()));hi=max(hi,float(eig.max()))
            geometry[i],random[i],shrinkages[i]=W,control,shrink
        np.savez(out/(name+'.npz'),parent=parents,geometry=geometry,random=random,shrinkage=shrinkages)
        summaries[name]=dict(count=count,minimum_eigenvalue=lo,maximum_eigenvalue=hi,
            random_spectrum_max_difference=eig_error,all66_retained=True,
            source=identity(path),cache=identity(out/(name+'.npz')))
        allparents.append(set(parents));del J,geometry,random
    if allparents[0]&allparents[1]:raise ValueError('TRAIN/VAL parent overlap')
    contract=dict(format='p3r4_geometry_CFM_metric_cache_v1',status='COMPLETE',baseline=identity(baseline),
        train=summaries['train'],validation=summaries['validation'],rotation=identity(out/'ROTATION.npy'),rotation_seed=ROTATION_SEED,
        normalization=identity(out/'NORMALIZATION.npz'),metric_coordinate='TRAIN-standardized residual velocity',
        map='rawx=center+T*r; T=frame_inverse.T@diag(TRAIN_x_scale); Jx=Jq*probitDerivative(anchor_x); Jr=Jx@T',
        metric='I+3*V*diag(singular(Jr)^2/(1+singular(Jr)^2))*V.T; stable SVD, no Gram eigenvalue clipping',
        random='R W_r R.T with fixed TRAIN-independent orthogonal R seed20261025191',
        bounds=[1,4],raw_metric='T^{-T}W_rT^{-1}; bounded relative to T^{-T}T^{-1}, not raw-Euclidean1..4',
        endpoint_GT_used=False,decoder_used=False,physical_RHS=0,optimizer_updates=0,
        preprocessing_wall_seconds=time.monotonic()-started,preprocessing_CPU_seconds=time.process_time()-cpu)
    write_json(out/'CACHE_COMPLETE.json',contract)
    return contract


def verify_cache(directory,train,validation):
    directory=Path(directory);c=json.loads((directory/'CACHE_COMPLETE.json').read_text())
    if c['status']!='COMPLETE' or c['bounds']!=[1,4] or c['endpoint_GT_used'] is not False:raise ValueError('frozen y-only bounded metric cache required')
    for key in ('baseline','rotation','normalization'):verify_identity(c[key])
    for name,path in [('train',train),('validation',validation)]:
        verify_identity(c[name]['source']);verify_identity(c[name]['cache'])
        if identity(path)!=c[name]['source']:raise ValueError('metric-cache dataset mismatch')
    return c


def run(args):
    wall,cpu=time.monotonic(),time.process_time()
    torch.set_num_threads(1)
    if not os.environ.get('R4_JOBDIR'):raise ValueError('use released training MPS Book wrapper, not unledgered execution')
    release=json.loads(args.release.read_text())
    if release.get('released') is not True or release.get('format')!='p3r4_geometry_weighted_CFM_release_v1':raise ValueError('inert until explicit parent release')
    if (args.seed,args.kind) not in [(r['seed'],r['kind']) for r in release['fits']]:raise ValueError('fit not released')
    if (args.updates,args.lr,args.batch_size)!=(3000,5e-5,128):raise ValueError('fixed3000/LR5e-5/batch128 required')
    paired=json.loads(Path(release['paired_VI_complete']).read_text())
    if paired.get('status')!='COMPLETE':raise ValueError('paired VI must be CLOSED COMPLETE before weighted CFM fits')
    for source in release['source_identities']:verify_identity(source)
    if args.device!='mps' or not torch.backends.mps.is_available():raise ValueError('released bounded MPS job required')
    cache=verify_cache(args.metrics,args.train,args.validation)
    if identity(args.metrics/'CACHE_COMPLETE.json')!=release['metric_cache']:raise ValueError('unreleased metric cache')
    if args.out.exists():raise FileExistsError('fresh weighted-CFM destination only')
    args.out.mkdir(parents=True)
    baseline_record=next(r for r in release['baselines'] if r['seed']==args.seed)
    verify_identity(baseline_record['checkpoint']);verify_identity(baseline_record['receipt'])
    if json.loads(Path(baseline_record['receipt']['path']).read_text())['status']!='COMPLETE':raise ValueError('baseline must be CLOSED COMPLETE')
    baseline=torch.load(baseline_record['checkpoint']['path'],map_location='cpu',weights_only=False)
    stats,ac=copy.deepcopy(baseline['normalization']),copy.deepcopy(baseline['anchor_contract'])
    if baseline['kind']!='cfm' or ac['coordinate_preconditioner']!='none':raise ValueError('matched anchorCFMnone baseline required')
    train=load_anchor_data(args.train);val=load_anchor_data(args.validation)
    if (len(train['x']),len(val['x']))!=(1024,256) or set(train['parent'])&set(val['parent']):raise ValueError('fixed independent1024/256 required')
    with np.load(args.metrics/'NORMALIZATION.npz',allow_pickle=False) as z:
        if any(not np.array_equal(stats[k],z[k]) for k in z.files):raise ValueError('baseline seed statistics differ from frozen TRAIN normalizers')
    matched=SimpleNamespace(coordinate_preconditioner='none',alignment_rotation=None,j_summary='none',
        anchor_scale=ac['anchor_scale'],residual_center=ac['residual_center'],preconditioner_base='raw',kind='cfm')
    # Recompute only to verify exact TRAIN provenance, never refit VAL stats.
    center=train['base_mean'] if train['base_mean'] is not None else train['context_x']
    from .anchor_density import context_features
    replay=fit_normalization(train['x']-center,context_features(train['y'],train['context_x']),stats['frame'])
    if any(not np.array_equal(stats[k],replay[k]) for k in ('x_mean','x_scale','y_mean','y_scale')):raise ValueError('TRAIN-statistics replay differs from original baseline')
    for data,name in [(train,'train'),(val,'validation')]:
        with np.load(args.metrics/(name+'.npz'),allow_pickle=False) as z:
            if not np.array_equal(z['parent'].astype(str),data['parent']):raise ValueError('metric row order mismatch')
            data['velocity_metric']=z[args.kind]
        eig=np.linalg.eigvalsh(data['velocity_metric'].astype(np.float64))
        if eig.min()<1-2e-6 or eig.max()>4+2e-6:raise ValueError('stored float32 metric eigenbound failure')
    tr,va=prepare(train,stats,matched),prepare(val,stats,matched)
    torch.manual_seed(args.seed);rng=np.random.default_rng(args.seed)
    generator=torch.Generator().manual_seed(args.seed+1103)
    config=copy.deepcopy(baseline['model_config']);model=make_model('cfm',**config).to(args.device)
    count=sum(p.numel() for p in model.parameters())
    if count!=baseline['contract']['parameter_count'] or config['dimension']!=66 or config['context_dim']!=578:raise ValueError('baseline architecture/capacity mismatch')
    ac.update(metric='bounded_geometry' if args.kind=='geometry' else 'bounded_random_alignment',
        metric_scope='TRAIN-standardized velocity ERROR; Jr=Jx*T; full66 SPD1..4; no hard projection',
        metric_contract=cache,metric_eigenvalue_checks={'bounds':[1,4]})
    criterion='validation_latent_energy_score_raw_x66'
    contract=dict(format=FORMAT,kind='cfm',weighted_kind=args.kind,seed=args.seed,model_config=config,parameter_count=count,
        train_size=1024,validation_size=256,full66_model=True,cfm_base='anchor',anchor_contract=ac,
        train=str(args.train),validation=str(args.validation),baseline=baseline_record,initialization='fresh same architecture/per-seed Torch initialization as original baseline',
        normalizer='bitwise original baseline TRAIN statistics; TRAIN replay checked',learning_rate=5e-5,updates=3000,
        batch_size=128,eval_every=100,validation_batch_size=128,validation_samples=4,validation_nfe=16,weight_decay=1e-4,
        criterion=criterion,update_zero_eligible=True,physical_RHS=0,posterior_fidelity='NOT_ESTABLISHED',
        CFM_density='numerical continuous flow only; weighting does not fix failed numeric-density gate',
        input='raw512+anchor66 only; cached y-only J enters training ERROR metric, not decoder or conditioning features',
        velocity_optimum='SPD W(y) preserves ideal conditional mean velocity; finite fit/calibration not guaranteed',device=args.device)
    write_json(args.out/'RUN_CONTRACT.json',contract)
    np.savez(args.out/'NORMALIZATION.npz',**{k:v for k,v in stats.items() if isinstance(v,np.ndarray)})
    optimizer=torch.optim.AdamW(model.parameters(),lr=5e-5,weight_decay=1e-4)
    history=[];best=float('inf');selected=0;state=None;update=0;trainloss=None;stop='UPDATE_LIMIT';requested=[False]
    old={s:signal.signal(s,lambda *_:requested.__setitem__(0,True)) for s in (signal.SIGTERM,signal.SIGINT)}
    def batch(data,ids):
        return [torch.from_numpy(data[k][ids]).to(args.device) for k in ('x','y','base_mean','base_chol','velocity_metric')]
    def payload(state,chosen,value):
        return dict(format=FORMAT,kind='cfm',seed=args.seed,model_config=config,normalization=stats,
            anchor_contract=ac,contract=contract,state_dict=state,selected_update=chosen,
            selection_metric=value,criterion=criterion,posterior_fidelity='NOT_ESTABLISHED')
    def evaluate():
        nonlocal best,selected,state
        model.eval();vg=torch.Generator().manual_seed(args.seed+7019);values=[];fm=[]
        for begin in range(0,256,128):
            x,y,mean,chol,metric=batch(va,slice(begin,begin+128))
            with torch.no_grad():
                fm.append((float(cfm_error_loss(model,x,y,vg,mean,chol,metric).cpu()),len(x)))
                samples=model.sample(y,4,vg,nfe=16,base_mean=mean,base_chol=chol)
                scale=torch.as_tensor(stats['x_scale'],dtype=x.dtype,device=args.device)
                value=torch.linalg.vector_norm((samples-x[:,None])*scale,dim=-1).mean(1)
                for i in range(4):
                    for j in range(i):value-=torch.linalg.vector_norm((samples[:,i]-samples[:,j])*scale,dim=-1)/12
            values.append(value.cpu().numpy())
        values=np.concatenate(values).astype(np.float64);score=float(values.mean())
        if not np.isfinite(values).all():raise FloatingPointError('nonfinite fixed validation score')
        row=dict(update=update,train_minibatch_loss=trainloss,validation_metric=score,
            validation_weighted_fm_mse=sum(v*n for v,n in fm)/256,learning_rate=5e-5,wall_seconds=time.monotonic()-wall)
        history.append(row)
        with (args.out/'loss_curve.jsonl').open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
        current={k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
        if score<best:
            best,selected,state=score,update,current
            torch.save(payload(state,selected,best),args.out/'BEST.pt')
            np.savez(args.out/'validation_metric.npz',parent=val['parent'],metric=values)
        torch.save(dict(**payload(current,update,score),optimizer=optimizer.state_dict(),generator_state=generator.get_state(),
            numpy_rng_state=rng.bit_generator.state,torch_rng_state=torch.get_rng_state()),args.out/'LAST.pt')
        print(json.dumps(dict(event='validation',**row,selected_update=selected)),flush=True)
    try:
        evaluate()
        for update in range(1,3001):
            if requested[0] or time.monotonic()-wall>=args.max_seconds:
                update-=1;stop='SIGNAL' if requested[0] else 'WALL_LIMIT';break
            model.train();ids=rng.integers(0,1024,128);x,y,mean,chol,metric=batch(tr,ids)
            optimizer.zero_grad(set_to_none=True)
            loss=cfm_error_loss(model,x,y,generator,mean,chol,metric)
            if not bool(torch.isfinite(loss)):raise FloatingPointError('nonfinite weighted CFM loss')
            loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.,error_if_nonfinite=True);optimizer.step()
            trainloss=float(loss.detach().cpu())
            if update%100==0:evaluate()
        if history[-1]['update']!=update:evaluate()
    finally:
        for s,h in old.items():signal.signal(s,h)
    result=dict(**contract,selected_update=selected,selected_initialization=selected==0,selected_validation_metric=best,
        completed_updates=update,stop_reason=stop,history=history,wall_seconds=time.monotonic()-wall,
        cpu_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    write_json(args.out/'RESULT.json',result)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('train','validation','metrics','out','release'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--kind',choices=['geometry','random'],required=True);p.add_argument('--seed',type=int,required=True)
    p.add_argument('--updates',type=int,default=3000);p.add_argument('--lr',type=float,default=5e-5)
    p.add_argument('--batch-size',type=int,default=128);p.add_argument('--device',choices=['cpu','mps'],default='mps')
    p.add_argument('--max-seconds',type=float,default=540.)
    run(p.parse_args())


if __name__=='__main__':main()
