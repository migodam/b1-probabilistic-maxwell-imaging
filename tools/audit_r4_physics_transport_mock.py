"""Zero-physics numerical verification of the proposed transport trainer.

No real Target is constructed. Fake native-q likelihoods and fake counters
exercise the source-adjoint/probit bridge; these counters are not paid RHS.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['R4_ARRAY_CODEC']='stored'
from pathlib import Path
import sys,time,weakref,gc
from types import SimpleNamespace
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import Book,OUT,pack,write
from p3r4.density import make_model,DensityPredictor
from p3r4.anchor_density import AnchorDensityPredictor
from p3r4.conditional_coordinates import fixed_random_rotation
from p3r4.train_physics_transport import (ActualPhi,TargetLikelihoodOracle,DifferentiableDensity,
    unpack_observations,validation_updates,lambda_at,release_template,validate_release)
from p3r4.train_physics_transport import load_original_observations,persist_transport_batch,original_batch


class MockPoint:
    def __init__(self,q,key):self.q,self.key=q,key


class FakeNativeTarget:
    """Analytic native-q quadratic; never imports/constructs real physics."""
    def __init__(self,costs,jobdir):
        self.costs=costs
        self.book=SimpleNamespace(jobdir=jobdir,context={})
        self.model=SimpleNamespace(_transfers={})
        self.calls=0
        self.fail_gradient=False
        self.live=weakref.WeakSet()
        self.A=np.random.default_rng(117).normal(size=(66,66))/8
    def evaluate(self,q,purpose):
        self.calls+=1
        before=self.book.jobdir/f'transport_queries/q{self.calls:08d}/BEFORE.npz'
        assert before.exists(), 'rawx and nativeY were not saved before mock Full'
        with np.load(before,allow_pickle=False) as z:
            assert z['x'].shape==(66,) and z['observation'].dtype==np.complex128
            assert np.array_equal(z['q'],q)
        self.costs['RHS']+=8
        p=MockPoint(q,f'mock_{self.calls:08d}');self.live.add(p);return p
    def phi(self,p,y):
        residual=self.A@p.q-y.real.ravel()[:66]
        return .5*float(residual@residual)
    def phi_gradient(self,p,y):
        self.costs['RHS']+=8
        if self.fail_gradient:raise RuntimeError('intentional analytic mock gradient failure')
        return self.A.T@(self.A@p.q-y.real.ravel()[:66])


def payload(kind,anchor=False,mode='none'):
    cfg=dict(dimension=66,context_dim=578 if anchor else 512,width=8,layers=4,bins=5,output_rank=4,diag_hidden=12)
    model=make_model(kind,**cfg).double()
    # Exercise nonidentity transforms while keeping splines well conditioned.
    with torch.no_grad():
        if kind=='spline':
            model.flow.affine.bias.normal_(std=.02)
            for layer in model.flow.layers:
                layer.net[-1].weight.normal_(std=.01)
                layer.net[-1].bias.add_(torch.randn_like(layer.net[-1].bias)*.01)
        else:model.head[-1].bias.normal_(std=.02)
    rng=np.random.default_rng(442)
    frame=np.linalg.qr(rng.normal(size=(66,66)))[0]
    stats=dict(frame=frame,frame_inverse=frame.T,frame_logabsdet=float(torch.linalg.slogdet(torch.tensor(frame)).logabsdet),
        x_mean=rng.normal(size=66)*.2,x_scale=np.exp(rng.normal(size=66)*.1),
        y_mean=rng.normal(size=cfg['context_dim'])*.1,y_scale=np.exp(rng.normal(size=cfg['context_dim'])*.1))
    result=dict(format='p3r4_anchor_full66_density_v1' if anchor else 'p3r4_full66_density_v1',
        kind=kind,model_config=cfg,normalization=stats,state_dict=model.state_dict())
    if anchor:
        result['anchor_contract']=dict(j_summary='none',residual_center=True,anchor_scale=1.,
            coordinate_preconditioner=mode,alignment_rotation=fixed_random_rotation().tolist() if mode=='random' else None)
    return result


def main():
    torch.set_num_threads(1);torch.manual_seed(414)
    job='physics_transport_mock_'+str(time.time_ns())
    book=Book('I',job,analysis_only=True,max_wall=120)
    out=OUT/'AUDIT'/job;out.mkdir(parents=True)
    status,error='COMPLETE',None
    try:
        rng=np.random.default_rng(111)
        complex_y=rng.normal(size=(2,32,8))+1j*rng.normal(size=(2,32,8))
        packed=np.array([pack(y) for y in complex_y])
        assert np.array_equal(unpack_observations(packed),complex_y)
        costs=dict(RHS=0)
        target=FakeNativeTarget(costs,out/'native_mock')
        oracle=TargetLikelihoodOracle(target,costs,10**9)
        x=torch.tensor(rng.normal(size=(2,66)),dtype=torch.double,requires_grad=True)
        assert torch.autograd.gradcheck(lambda xx:ActualPhi.apply(xx,complex_y,oracle),(x,),atol=2e-6,rtol=2e-5)
        gc.collect();assert len(target.live)==0
        # Backward itself must never invoke an oracle or change fake RHS.
        value=ActualPhi.apply(x,complex_y,oracle).sum();before=costs['RHS'];value.backward();assert costs['RHS']==before
        rows=[]
        for kind in ('spline','diag_gaussian'):
            for anchor,mode in ((False,'none'),(True,'none'),(True,'gn'),(True,'random')):
                p=payload(kind,anchor,mode)
                d=DifferentiableDensity(p,'cpu',torch.double)
                data=dict(y=packed,context_x=rng.normal(size=(2,66)),context_J=rng.normal(size=(2,512,66))*.025,
                          base_mean=None,base_chol=None)
                # Any attempt to read a label during proposal preparation fails.
                class Poison:
                    def __getitem__(self,k):raise AssertionError('GT accessed by proposal preparation')
                data['x']=Poison()
                prep=d.prepare(data,np.arange(2))
                eps=torch.tensor(rng.normal(size=(2,66)),dtype=torch.double)
                sample,logq=d.reparameterize(eps,prep)
                direct=d.log_prob(sample,prep)
                agreement=float((direct-logq).abs().max().detach())
                assert agreement<2e-9
                # Compare to unchanged public inference path, which usesfloat32.
                pred=AnchorDensityPredictor(p) if anchor else DensityPredictor(p)
                kwargs=dict(context_x=data['context_x'],context_J=data['context_J']) if anchor else {}
                old=pred.log_prob(sample.detach().numpy(),packed,**kwargs)
                inference_error=float(np.max(abs(old-logq.detach().numpy())))
                assert inference_error<1e-4
                # VJP through all maps and density entropy. Central difference
                # perturbs one actual mean-head parameter; no optimizer step.
                theta=d.model.flow.affine.bias if kind=='spline' else d.model.head[-1].bias
                d.model.zero_grad(set_to_none=True)
                xx,lq=d.reparameterize(eps,prep)
                objective=(ActualPhi.apply(xx,complex_y,oracle)+.5*xx.square().sum(1)+lq).mean()
                objective.backward()
                derivative=float(theta.grad[0])
                initial=float(theta[0].detach());h=1e-5;values=[]
                for displacement in (h,-h):
                    with torch.no_grad():theta[0]=initial+displacement
                    xx,lq=d.reparameterize(eps,prep)
                    ph=oracle.values(xx.detach().numpy(),complex_y)
                    values.append(float((torch.tensor(ph)+.5*xx.square().sum(1)+lq).mean().detach()))
                with torch.no_grad():theta[0]=initial
                finite_difference=(values[0]-values[1])/(2*h)
                vjp_error=abs(derivative-finite_difference)
                assert vjp_error<2e-6
                assert all(p.grad is None or torch.isfinite(p.grad).all() for p in d.model.parameters())
                gc.collect();assert len(target.live)==0
                rows.append(dict(kind=kind,anchor=anchor,coordinate=mode,
                    sample_inverse_logq_max_error=agreement,unchanged_predictor_float32_max_error=inference_error,
                    mock_reparameterized_VJP=derivative,central_difference=finite_difference,VJP_absolute_error=vjp_error,
                    GT_read_by_proposal=False,live_Points_after_call=len(target.live)))
        assert validation_updates(3000,300)==list(range(0,3001,300))
        schedule=[dict(start=1,end=300,value=.001),dict(start=301,end=900,value=.01),dict(start=901,end=1800,value=.1),dict(start=1801,end=3000,value=1.)]
        assert [lambda_at(i,schedule) for i in (1,300,301,900,901,1800,1801,3000)]==[.001,.001,.01,.01,.1,.1,1.,1.]
        raw=OUT/'models/density/density_spline_32768_s2026102501/BEST.pt'
        template=release_template(raw,OUT/'inputs/joint/train32768.npz',OUT/'inputs/joint/validation.npz',[2026102501])
        try:validate_release(template,SimpleNamespace(seed=2026102501,mode='physics'))
        except ValueError:pass
        else:raise AssertionError('unreleased template accepted')
        template['released']=True
        physics_plan=validate_release(template,SimpleNamespace(seed=2026102501,mode='physics'))
        control_plan=validate_release(template,SimpleNamespace(seed=2026102501,mode='supervised'))
        assert physics_plan['planned_RHS']==203264 and control_plan['planned_RHS']==11264
        assert template['inclusive_NLL_batch']==128 and template['array_codec']=='stored'
        # Verify all1056 actual observation-only files and float32 NN packing;
        # only y and parent arrays are read, never labels/material/GT.
        joint=OUT/'inputs/joint'
        with np.load(joint/'train32768.npz',allow_pickle=False) as z:
            train=dict(y=z['y'][:1024],parent=z['parent'][:1024])
        with np.load(joint/'validation.npz',allow_pickle=False) as z:
            val=dict(y=z['y'][:32],parent=z['parent'][:32])
        original=load_original_observations(template,train,val,np.arange(1024),np.arange(32))
        assert len(original)==1056 and all(y.dtype==np.complex128 for y in original.values())
        assert any(np.any(y!=y.astype(np.complex64).astype(np.complex128)) for y in original.values())
        assert original_batch(original,train['parent'][:4]).dtype==np.complex128
        # A rounded physical Y is rejected rather than silently upcast/fallback.
        before_queries=oracle.query_count
        try:oracle.values(x.detach().numpy(),complex_y.astype(np.complex64))
        except ValueError:pass
        else:raise AssertionError('rounded physical Y accepted')
        assert oracle.query_count==before_queries
        # Exercise persisted failed-query metadata and release of its Point.
        target.fail_gradient=True
        try:oracle.value_and_gradient(x.detach().numpy()[:1],complex_y[:1])
        except RuntimeError:pass
        else:raise AssertionError('mock failure did not propagate')
        failed=target.book.jobdir/f'transport_queries/q{oracle.query_count:08d}/AFTER.json'
        import json
        assert json.loads(failed.read_text())['status']=='FAILED'
        target.fail_gradient=False
        gc.collect();assert len(target.live)==0
        batch=out/'batch_persistence_mock'
        persist_transport_batch(batch,'BEFORE',np.zeros((4,66)),np.ones(4),np.ones(4),np.array(['p']*4),1,.001,np.arange(1,5))
        persist_transport_batch(batch,'AFTER',np.zeros((4,66)),np.ones(4),np.ones(4),np.array(['p']*4),1,.001,np.arange(1,5),phi=np.arange(4))
        with np.load(batch/'AFTER.npz',allow_pickle=False) as z:
            assert {'x','logq_x','prior_energy','parent','update','lambda_value','query_index','Phi'}==set(z.files)
        # Independent inclusive/reverse RNGs retain paired128-label sequences.
        left,right=np.random.default_rng(123),np.random.default_rng(123)
        reverse_rng=np.random.default_rng(2360)
        for _ in range(8):
            assert np.array_equal(left.integers(0,32768,size=128),right.integers(0,32768,size=128))
            reverse_rng.integers(0,1024,size=4)
        write(out/'RESULT.json',dict(status='PASS',mock_cases=rows,custom_autograd_double_gradcheck='PASS',
            source_major_unpack_bitwise='PASS',backward_new_oracle_calls=0,no_retained_Points=True,
            release_gate_rejects_unreleased=True,physics_plan=physics_plan,control_plan=control_plan,
            original_observations_verified=1056,physical_Y_dtype='original complex128',rounded_Y_rejected=True,
            NN_sourcepack_float32_matches=True,original_observations_cached_bytes=sum(y.nbytes for y in original.values()),
            every_mock_query_saved_before_Full=True,failed_mock_query_metadata_preserved=True,
            batch_raw_x_logq_Phi_persistence=True,inclusive_NLL_batch=128,paired_inclusive_RNG_sequence='PASS',
            actual_physics_RHS=0,actual_Target_instances=0,optimizer_updates=0,checkpoint_modifications=0,
            fake_counter_scope='analytic mock8/16 counters only; not paid physical RHS',
            posterior_qualification='NOT_TESTED',receipt=str(book.jobdir/'CLOSED.json')))
        print(out)
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
