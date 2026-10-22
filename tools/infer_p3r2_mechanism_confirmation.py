"""Truth-free twenty-one-model inference on the separately frozen second cohort."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
import numpy as np
import torch
from p3r2.core import OUT,write,save
from p3r2.residual_risk import load_residual_predictor
from p3r2_mechanism.masked_risk import load_mechanism_predictor
from p3r2_mechanism.random_sensitivity import load_random_sensitivity_predictor
from p3r2_mechanism.contract import VERSION,MODEL_COUNT,PARENTS,validate_freeze,validate_online_inputs,read_online,file_record,checkpoint_metadata


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=OUT);p.add_argument('--manifest',type=Path);p.add_argument('--resume',action='store_true');a=p.parse_args()
    out=a.root.resolve();manifest=a.manifest or out/'MECHANISM_CONFIRMATION_FROZEN.json'
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    m,paths,freeze_record=validate_freeze(out,manifest);validate_online_inputs(out,m)
    dest=out/'L/mechanism_confirmation/estimates'
    if (dest/'COMPLETE.json').exists():raise FileExistsError('second-cohort mechanism inference already complete')
    dest.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    cfg=json.loads(paths['target'].read_text());points=np.array(json.loads(paths['geometry'].read_text())['points']);roi=np.flatnonzero(np.linalg.norm(points,axis=1)<=cfg['roi_radius'])
    if len(roi)!=256 or cfg['field_mean']!=0:raise ValueError('original symmetric P1 ROI/mean required')
    host=complex(*cfg['host_chi']);prior_mean=np.full(576,host,dtype=complex)
    prior_mean[roi]+=.5*complex(np.mean(cfg['inclusion_increment_real_interval']),np.mean(cfg['inclusion_increment_imag_interval']))
    models={};failures={};metadata={}
    for entry in m['models']:
        name=entry['name'];metadata[name]=checkpoint_metadata(paths['models'][name],entry)
        try: models[name]={'mechanism':load_mechanism_predictor,'random_sensitivity':load_random_sensitivity_predictor,'original':load_residual_predictor}[entry['loader']](paths['models'][name],device='cpu')
        except Exception as error:models[name]=None;failures[name]=repr(error)
    write(dest/'INFERENCE_PROTOCOL.json',dict(version=VERSION,manifest_record=freeze_record,source_revision=m['source_revision'],frozen_unix=m['frozen_unix'],
        parents=PARENTS,models=m['models'],checkpoint_metadata=metadata,GT_read=False,physics_calls=0,objective='SUPERVISED_RISK',
        estimator_action='raw_unprojected_no_blend',primary_metric='complex_fullgrid_MSE',
        fallback='Per-method anchor; if context fails, paid probit cache/native point, finally exact prior mean. No parent excluded.',load_failures=failures))
    statuses=[]
    for parent in PARENTS:
        path=dest/f'{parent}.npz'
        if path.exists():
            if not a.resume:raise FileExistsError(path)
            with np.load(path,allow_pickle=False) as z:
                if int(z['version'])!=VERSION or str(z['source_revision'])!=m['source_revision'] or float(z['frozen_unix'])!=float(m['frozen_unix']) or json.loads(str(z['manifest_record_json']))!=freeze_record:raise ValueError('existing estimate freeze mismatch')
                for entry in m['models']:
                    if z[entry['name']].shape!=(576,) or not np.isfinite(z[entry['name']]).all():raise ValueError('existing estimate invalid')
                statuses.append(json.loads(str(z['status_json'])))
            continue
        status=dict(parent=parent,failures={},method_seconds={},GT_read=False);anchor_source='online_probit_context'
        try:online=read_online(out/'L/contexts'/f'{parent}.npz');anchor=online['chi0'][0].copy()
        except Exception as error:
            online=None;anchor=None;status['failures']['context']=repr(error)
            for fallback,label in ((out/'L/contexts/cache'/f'{parent}.npz','paid_probit_cache'),(out/'L/contexts/details'/parent/'NATIVE.npz','native_point')):
                try:
                    with np.load(fallback,allow_pickle=False) as z:chi=z['chi'].copy()
                    if chi.shape!=(576,) or not np.isfinite(chi).all():raise ValueError('invalid fallback material')
                    anchor=chi;anchor_source=label;break
                except Exception:pass
            if anchor is None:anchor=prior_mean.copy();anchor_source='exact_prior_mean'
        status['anchor_source']=anchor_source;arrays=dict(anchor=anchor,estimator_action=np.array('raw_unprojected_no_blend'))
        for entry in m['models']:
            name=entry['name'];tick=time.monotonic()
            try:
                if online is None or models[name] is None:raise RuntimeError(failures.get(name,'online context unavailable'))
                chi=models[name].predict_full(online)[0]
                if chi.shape!=(576,) or not np.isfinite(chi).all():raise ValueError('invalid material point')
            except Exception as error:chi=anchor.copy();status['failures'][name]=repr(error)
            arrays[name]=chi;arrays[name+'__used_fallback']=np.array(name in status['failures']);status['method_seconds'][name]=time.monotonic()-tick
        save(path,**arrays,parent=np.array(parent),roi_indices=roi,version=np.array(VERSION),source_revision=np.array(m['source_revision']),
            frozen_unix=np.array(m['frozen_unix']),manifest_record_json=np.array(json.dumps(freeze_record)),status_json=np.array(json.dumps(status,allow_nan=False)))
        statuses.append(status);print(json.dumps(dict(parent=parent,complete=True,fallbacks=list(status['failures']))),flush=True)
    # Revalidate the full frozen source/checkpoint registry after all inference.
    _,_,after=validate_freeze(out,manifest)
    if after!=freeze_record or len(statuses)!=256:raise ValueError('freeze changed or incomplete inference')
    names=[e['name'] for e in m['models']]
    online_records=[file_record(path) for parent in PARENTS for path in (out/'L/contexts'/f'{parent}.npz',out/'L/contexts'/f'{parent}.json',out/'inputs/online'/f'{parent}.npz')]
    write(dest/'COMPLETE.json',dict(status='COMPLETE',version=VERSION,parents=PARENTS,count=256,model_count=MODEL_COUNT,model_names=names,
        source_revision=m['source_revision'],frozen_unix=m['frozen_unix'],manifest_record=freeze_record,checkpoint_metadata=metadata,
        estimate_files=[file_record(dest/f'{parent}.npz') for parent in PARENTS],online_files=online_records,
        fallback_counts={name:sum(name in s['failures'] for s in statuses) for name in names+['context']},
        GT_read=False,physics_calls=0,objective='SUPERVISED_RISK',estimator_action='raw_unprojected_no_blend',primary_metric='complex_fullgrid_MSE',wall_seconds=time.monotonic()-start,posterior_claim=False))


if __name__=='__main__':main()
