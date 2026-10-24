"""Frozen cached-geometry preprocessing and zero-fit algebra/readiness audit."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,sys,time
from unittest.mock import patch
import numpy as np
import torch
from scipy.linalg import eigvalsh
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r4.core import Book,OUT,write
from p3r4.train_geometry_cfm import prepare_cache,network_metric,fixed_rotation
from p3r4.train_anchor_density import cfm_error_loss
from p3r4.density import make_model
from p3r4.conditional_coordinates import native_to_gaussianized_jacobian
from p3r4.train_physics_transport import verify_identity
import queue_r4_geometry_cfm as controller


def main():
    torch.set_num_threads(1)
    job='geometry_weighted_cfm_prepare_'+str(time.time_ns())
    book=Book('I',job,analysis_only=True,max_wall=600)
    folder=OUT/'AUDIT'/job;folder.mkdir(parents=True)
    status,error='COMPLETE',None
    try:
        rng=np.random.default_rng(319)
        frame=np.linalg.qr(rng.normal(size=(66,66)))[0]
        stats=dict(frame_inverse=frame.T,x_scale=np.exp(rng.normal(size=66)*.5))
        anchor=rng.normal(size=66);J=rng.normal(size=(512,66))*.1
        W,shrink=network_metric(anchor,J,stats);R=fixed_rotation();C=R@W@R.T
        eig=eigvalsh(W);control_error=float(np.max(abs(eigvalsh(C)-eig)))
        assert eig.min()>=1-1e-10 and eig.max()<=4+1e-10 and control_error<1e-10
        assert np.array_equal(fixed_rotation(),R)
        T=stats['frame_inverse'].T*stats['x_scale'][None,:]
        Ti=np.linalg.solve(T,np.eye(66));raw=Ti.T@W@Ti;baseline=Ti.T@Ti
        relative=eigvalsh(raw,baseline)
        assert relative.min()>=1-1e-10 and relative.max()<=4+1e-10
        # Baseline identity-metric loss and weighted ERROR differentiation.
        cfg=dict(dimension=66,context_dim=578,width=128,blocks=4,time_features=32,velocity_hidden=547)
        torch.manual_seed(2026102501);model=make_model('cfm',**cfg)
        first={k:v.detach().clone() for k,v in model.state_dict().items()}
        torch.manual_seed(2026102501);second=make_model('cfm',**cfg)
        assert all(torch.equal(first[k],v) for k,v in second.state_dict().items())
        n=4;x=torch.randn(n,66);y=torch.randn(n,578);mean=torch.zeros(n,66);chol=torch.eye(66).repeat(n,1,1)
        a=cfm_error_loss(model,x,y,torch.Generator().manual_seed(51),mean,chol,None)
        b=cfm_error_loss(model,x,y,torch.Generator().manual_seed(51),mean,chol,torch.eye(66).repeat(n,1,1))
        assert torch.allclose(a,b,atol=1e-6,rtol=1e-6)
        e=torch.randn(n,66,dtype=torch.double,requires_grad=True);w=torch.tensor(W).repeat(n,1,1)
        loss=torch.einsum('bi,bij,bj->b',e,w,e).mean()/66;loss.backward()
        expected=2*torch.einsum('bij,bj->bi',w,e.detach())/(n*66)
        assert torch.allclose(e.grad,expected,atol=1e-12,rtol=1e-12)
        base_loss=float(e.detach().square().mean())
        assert base_loss-1e-10<=float(loss.detach())<=4*base_loss+1e-10
        targets=rng.normal(size=(7,66));conditional_mean=targets.mean(0)
        optimum_residual=np.max(abs(((conditional_mean-targets)@W).mean(0)))
        assert optimum_residual<1e-12
        # Actual cache reads only parent/context_x/context_J and frozen TRAIN
        # normalization, never x/material labels or any physical backend.
        source=OUT/'inputs/strong_context_1024_256';cache=OUT/'inputs/geometry_weighted_cfm_v1'
        baseline_path=OUT/'models/anchor_density/anchor_cfm_none_s2026102501/BEST.pt'
        contract=prepare_cache(source/'density_train.npz',source/'density_validation.npz',baseline_path,cache)
        with np.load(cache/'NORMALIZATION.npz',allow_pickle=False) as z:
            for seed in controller.SEEDS:
                payload=torch.load(OUT/f'models/anchor_density/anchor_cfm_none_s{seed}/BEST.pt',map_location='cpu',weights_only=False)
                assert all(np.array_equal(payload['normalization'][k],z[k]) for k in z.files)
                assert payload['model_config']==cfg
        registry=controller.build_registry(cache)
        if controller.REGISTRY.exists():raise FileExistsError('inert registry already exists; do not overwrite')
        write(controller.REGISTRY,registry);controller.validate(registry)
        with patch.object(controller.subprocess,'run',side_effect=AssertionError('no fit launch')) as launch:
            try:controller.validate(registry,True)
            except ValueError as exc:assert 'inert' in str(exc)
            else:raise AssertionError('unreleased six-fit plan accepted')
            assert launch.call_count==0
        paired=json.loads((OUT/'configs/PHYSICS_TRANSPORT_PAIRED_V2_RUNTIME_RECOVERY.json').read_text())
        for record in paired['source_identities']:verify_identity(record)
        write(folder/'FROZEN_INERT_REGISTRY.json',registry)
        write(folder/'READINESS.json',dict(status='PASS_PREPARED_INERT',metric_cache=controller.identity(cache/'CACHE_COMPLETE.json'),
            preprocessing=contract,raw_to_network_metric_clarification='Jr=Jx*T, boundednetworkWr; raw generalized bounds relative to baseline metric only',
            mock_eigen_bounds=[float(eig.min()),float(eig.max())],random_eigenspectrum_error=control_error,
            generalized_raw_eigen_bounds=[float(relative.min()),float(relative.max())],
            weighted_gradient_max_error=float((e.grad-expected).abs().max()),identity_weight_loss_error=float(abs(a-b)),
            conditional_mean_optimum_gradient_error=float(optimum_residual),fresh_initialization_replay_bitwise=True,
            baseline_normalizers_match_all3seeds=True,model_parameter_count=sum(p.numel() for p in model.parameters()),
            seed_protocol=controller.SEEDS,fits=6,CFM_numeric_density_still_unqualified=True,
            all_VI_frozen_source_identities_unchanged=True,physical_RHS=0,optimizer_updates=0,child_launches=0,
            receipt=str(book.jobdir/'CLOSED.json')))
        print(folder)
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
