"""Fit real joint-supervision Gaussian controls; no physics, locked or GT IO."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write,save
from p3l.simple_density import SimpleConditionalGaussian

def main():
    folder=OUT/'D/simple_density'
    if folder.exists():raise FileExistsError('fit already present; no automatic rerun')
    folder.mkdir(parents=True)
    config=dict(train=512,validation=128,labels='weak8 relative to cheap16 context + exact probit2',
        contexts='scalar1154; train-only standardization; image not used',
        dimensions=[0,8,16,32],ridge=[1.,10.,100.,1000.,10000.],
        covariance_shrinkage=[0.,.1,.5],covariance_floor=1e-6,
        primary_selection='minimum completed validation raw-coordinate NLL; all fit attempts retained',
        reason='All neural families select update0; test low-variance normalized statistical controls',
        teacher='independent prior simulations, not posterior distillation',
        probability_status='NORMALIZED_PROPOSAL_NOT_FULL_POSTERIOR_CERTIFICATE',
        frozen_before_validation_fit=True,locked_read=False,new_Maxwell_RHS=0)
    write(folder/'FROZEN_CONFIG.json',config)
    b=Book('D','fit_simple_gaussian_controls',max_wall=300,analysis_only=True)
    status,error='COMPLETE',None
    try:
        with np.load(OUT/'D/training/TRAIN.npz',allow_pickle=False) as z:
            X=z['scalar'].copy();Y=z['labels'].copy()
        with np.load(OUT/'D/training/VALIDATION.npz',allow_pickle=False) as z:
            Xv=z['scalar'].copy();Yv=z['labels'].copy()
        sm,ss=X.mean(0),X.std(0);ss=np.where(ss>0,ss,1.)
        ym,ys=Y.mean(0),Y.std(0);ys=np.where(ys>0,ys,1.)
        X=(X-sm)/ss;Y=(Y-ym)/ys
        _,_,Vt=np.linalg.svd(X,full_matrices=False)
        rows=[];best=None;global_best=None
        for dim in config['dimensions']:
            P=Vt[:dim].T;C=X@P
            for ridge in ([1.] if dim==0 else config['ridge']):
                W=np.zeros((0,10)) if dim==0 else np.linalg.solve(C.T@C+ridge*np.eye(dim),C.T@Y)
                residual=Y-C@W
                empirical=residual.T@residual/len(Y)
                for shrink in config['covariance_shrinkage']:
                    cov=(1-shrink)*empirical+(shrink+config['covariance_floor'])*np.eye(10)
                    arrays=dict(scalar_mean=sm,scalar_scale=ss,label_mean=ym,label_scale=ys,
                        projection=P,weights=W,covariance=cov)
                    model=SimpleConditionalGaussian(arrays)
                    nll=float(-model.log_prob(Yv[:,None,:],Xv,None).mean())
                    ident=f'd{dim}_ridge{ridge:g}_shrink{shrink:g}'
                    save(folder/'fits'/f'{ident}.npz',**arrays)
                    row=dict(id=ident,dimension=dim,ridge=ridge,shrinkage=shrink,
                        validation_raw_NLL=nll,checkpoint=f'D/simple_density/fits/{ident}.npz')
                    rows.append(row)
                    if best is None or nll<best['validation_raw_NLL']:best=row
                    if dim==0 and (global_best is None or nll<global_best['validation_raw_NLL']):global_best=row
        write(folder/'RESULTS.json',dict(status='FITTED_NORMALIZED_JOINT_SUPERVISION_CONTROLS',rows=rows,
            selected=best,unconditional_selected=global_best,train_parents=512,validation_parents=128,
            posterior_qualified=False,GT_direct_read=False,locked_read=False))
        print(json.dumps(dict(selected=best,unconditional=global_best)))
    except BaseException as exc:status,error='ERROR',repr(exc);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
