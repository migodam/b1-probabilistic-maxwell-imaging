#!/usr/bin/env python3
"""Frozen raw-y provider moments: Spline/DiagGaussianK32 or learnedFieldK0.

Observation-only generation; no fitting/Maxwell or trueposterior/independence claim.
Legacy Spline CLI/checkpoints remain accepted; new artifacts declare provider/K.
"""
from pathlib import Path
import argparse,json,sys,os,time,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np,torch
from p3r4.density import load_density
from p3r4.risk_fusion import PROVIDERS,parent_seed
from p3r4.probability_eval import read_online_observation,decoder_from_geometry,decode_gaussianized
from p3r4.evaluation import checkpoint_identity,write_new_json


def load_provider(checkpoint,provider,device):
    if provider=='field_moments':
        from p3r4.field_moments import load_field_moments
        model=load_field_moments(checkpoint,device)
        if model.payload.get('architecture')!='mlp_small' or model.payload.get('second_moments') is not True:raise ValueError('supervised rawyMLPsmall jointsecond-moment checkpoint required')
    else:
        model=load_density(checkpoint,device)
        if model.payload['kind']!=provider or model.model.dimension!=66 or model.model.context_dim!=512:raise ValueError('matching rawy512 full66 Spline/DiagGaussian provider required')
    return model


def provider_image(model,provider,y,seed,decoder):
    if provider=='field_moments':
        result=model.predict_moments(y,batch_size=1);mean=result['material_mean'][0];cov=result['psd_pointwise_covariance'][0]
        image=np.stack((mean.real,mean.imag,cov[:,0,0],cov[:,1,1])).reshape(4,24,24).astype('float32')
        raw=result['raw_pointwise_covariance'][0]
        extras=dict(raw_variance_re=raw[:,0,0],raw_variance_im=raw[:,1,1],
                    variance_rule=np.array('deterministic PSD-projected marginal diagonals from learned joint ReIm moments; no independence assumption'),
                    projection_warning=np.array(json.dumps(result.get('psd_projection_warning',{}),sort_keys=True)),
                    moment_status=np.array(result['second_moment_status']))
        return image,extras
    x=model.sample(y,32,seed,batch_size=8)[0];chi,_,_=decode_gaussianized(x,decoder)
    image=np.stack((chi.real.mean(0),chi.imag.mean(0),chi.real.var(0),chi.imag.var(0))).reshape(4,24,24).astype('float32')
    return image,dict(x_samples=x,variance_rule=np.array('K32 empirical model-pushforward variance ddof0; no independence assumption'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('checkpoint','observations','parents','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--provider',choices=tuple(PROVIDERS),default='spline');p.add_argument('--device',choices=('cpu','mps'),default='cpu');p.add_argument('--sample-seed',type=int,default=20261025121);a=p.parse_args()
    if a.out.exists():raise FileExistsError('new momentdirectory required')
    parents=json.loads(a.parents.read_text())['parents']
    if len(parents)!=len(set(parents)):raise ValueError('unique original parents')
    start,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);identity=checkpoint_identity(a.checkpoint);model=load_provider(a.checkpoint,a.provider,a.device);count=PROVIDERS[a.provider]
    cfg=json.loads((ROOT/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];decoder=decoder_from_geometry(cfg,points)
    packed_y=None;packed_lookup=None
    if not a.observations.is_dir():
        with np.load(a.observations,allow_pickle=False) as z:packed_parents=z['parent'].astype(str);packed_y=z['y'];packed_lookup={parent:i for i,parent in enumerate(packed_parents)}
        if len(packed_lookup)!=len(packed_parents):raise ValueError('unique packed observationparents')
    a.out.mkdir(parents=True);records=[];artifacts=[]
    for parent in parents:
        try:
            if a.observations.is_dir():y,_=read_online_observation(a.observations/f'{parent}.npz',parent,target_id=cfg['target_id'])
            else:
                # Materialize only raw y andparent; neverTruex/material/s/q.
                if parent not in packed_lookup:raise ValueError('one matching observationparent')
                y=packed_y[packed_lookup[parent]:packed_lookup[parent]+1]
            seed=parent_seed(parent,a.sample_seed);image,extras=provider_image(model,a.provider,y,seed,decoder)
            if image.shape!=(4,24,24) or not np.isfinite(image).all() or np.any(image[2:]<-1e-12):raise ValueError('invalid finite mean/variance provider output')
            dest=a.out/f'{parent}.npz'
            arrays=dict(parent=np.array(parent),moment_image=image,provider=np.array(a.provider),K=np.array(count),
                checkpoint_identity=np.array(json.dumps(identity,sort_keys=True)),moment_scope=np.array('frozen provider features, not trueposterior/no independence claim'),**extras)
            if count:arrays['seed']=np.array(seed)
            else:arrays['deterministic']=np.array(True)
            with dest.open('xb') as f:np.savez_compressed(f,**arrays)
            artifacts.append(checkpoint_identity(dest));records.append(dict(parent=parent,status='COMPLETE',provider=a.provider,K=count,seed=seed if count else None))
        except Exception as e:records.append(dict(parent=parent,status='FAILED',error=repr(e)))
    if checkpoint_identity(a.checkpoint)!=identity:raise ValueError('provider changed duringfeaturegeneration')
    write_new_json(a.out/'MOMENTS_COMPLETE.json',dict(status='MOMENTS_COMPLETE',moment_contract_version=2,provider=a.provider,K=count,records=records,parents=parents,artifacts=artifacts,
        checkpoint_identity=identity,sample_seed=a.sample_seed if count else None,deterministic=count==0,
        truth_read=False,physics_calls=0,optimizer_updates=0,scope='coarse frozen provider features asriskCNN context, not posterior advantage/calibration/no independence assumption',
        wall_seconds=time.monotonic()-start,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
    print(json.dumps(dict(status='MOMENTS_COMPLETE',provider=a.provider,K=count,parents=len(parents),failures=sum(r['status']=='FAILED' for r in records))))

if __name__=='__main__':main()
