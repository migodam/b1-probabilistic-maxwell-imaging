"""Finite new wiring check only: no GT, sweeps, or historical solver edits."""
from pathlib import Path
from collections import Counter
import json,time,sys
import numpy as np
from scipy.linalg import svd
from scipy.fft import fftn
from threadpoolctl import threadpool_limits
REPO=Path(__file__).resolve().parents[3];sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import fourier_apcs,run_cie
from research.research_max_struct_r5.src.cie_factorial import run_factorial_cie
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS,DiskPulseToeplitzGreen
from research.research_max_struct_r6.src.cie_fft import run_cie_fft


def run():
    started,cpu=time.perf_counter(),time.process_time()
    path=REPO/'research/research_max_struct_r2/data_public/obs_b20d5a55677d4eca.npz'
    allowed={'G','E0','S','Y','grid_shape','coords','radius','k','host'}
    with np.load(path,allow_pickle=False) as d: public={k:d[k].copy() for k in d.files if k in allowed}
    out=dict(status='RUNNING',classification='COMPUTATIONAL_WIRING_ONLY',truth_access='NONE',
        input=str(path.relative_to(REPO)),keys_opened=sorted(public),readout_parameters_tuned=False,
        execution='CPU one BLAS thread; finite mf10 current update and mf12 null-action audit')
    model=MaxwellModel(public);shape=model.shape;ns=model.ns
    cost=Counter();green=DiskPulseToeplitzGreen(public['coords'],shape,float(public['radius']),float(public['k']),cost)
    _,sing,Vh=svd(model.S,full_matrices=False,check_finite=False)
    rank=int(np.count_nonzero(sing>1e-3*sing[0]));Vd=Vh[:rank].conj().T
    rng=np.random.default_rng(123)
    x=rng.normal(size=(model.n,ns))+1j*rng.normal(size=(model.n,ns))
    y=rng.normal(size=(model.n,ns))+1j*rng.normal(size=(model.n,ns))
    gx=green.apply(x)
    def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
    def adj_error(left,right):return float(abs(left-right)/max(abs(left),abs(right),1e-300))
    out['Green']=dict(shape=shape,padded_shape=green.pad,
        dense_action_relative_difference=rel(gx,model.G@x),
        adjoint_pairing_relative_difference=adj_error(np.vdot(y,gx),np.vdot(green.adjoint(y),x)),
        spectrum_bytes=green.spectrum.nbytes)
    out['Fourier']=[]
    for mf in (10,12):
        basis=NativeFourierAPCS(shape,mf,Vd,cost)
        dense,meta=fourier_apcs(shape,mf,Vd)
        a=rng.normal(size=(basis.k,ns))+1j*rng.normal(size=(basis.k,ns))
        ua=basis.apply(a)
        item=dict(mf=mf,coefficient_count=basis.k,historical_apcs_rank=meta['apcs_rank'],
            historical_roundoff_singular_actions_removed=basis.k-meta['apcs_rank'],
            native_vs_SVD_cleaned_action_relative_difference=rel(ua,dense@a),
            native_adjoint_pairing_relative_difference=adj_error(np.vdot(y,ua),np.vdot(basis.adjoint(y),a)),
            explicit_U_GU_H_SU_bytes=(3*model.n+model.nr)*basis.k*16,
            native_coefficients_bytes=basis.k*ns*16)
        if mf==12:
            full=fftn(Vd[:,:1].reshape(shape+(1,)),axes=(0,1),norm='ortho')
            null_a=full[tuple(basis.frequencies.T)]
            native_null=basis.apply(null_a);dense_null=dense@null_a
            item.update(unit_DPCS_null_coefficient_norm=float(np.linalg.norm(null_a)),
                native_null_current_absolute_norm=float(np.linalg.norm(native_null)),
                historical_null_current_absolute_norm=float(np.linalg.norm(dense_null)),
                null_action_absolute_difference=float(np.linalg.norm(native_null-dense_null)),
                note='Near-zero output: absolute norms are authoritative; relative output error is not meaningful')
        out['Fourier'].append(item)
    out['one_current_update']=[]
    geom={k:public[k] for k in ('coords','radius','k')}
    Y=public['Y'];m0=np.full(model.n,float(public['host']))
    common=dict(method='cie',dpcs_rtol=1e-3,
        stages=[dict(mf=10,beta_current=1.,beta_material=.2,maxiter=1)],
        material_box=(.02,20.),coefficient_tol=0.,miniter=50,
        material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,
        physical_sigma=float(np.linalg.norm(Y)),cache_material_value_gradient=True)
    for arm in ('cie_dual','a_material'):
        results=[]
        for mode in ('dense','fft'):
            mm=MaxwellModel(public);cfg=dict(common,regularizer=lambda z,mm=mm:mm.reg.evaluate(z)[:2])
            cfg['sourcewise_cg']=arm=='cie_dual'
            record=[];cfg['iteration_callback']=lambda row:record.append(row)
            wall,ccpu=time.perf_counter(),time.process_time()
            if mode=='fft':result=run_cie_fft(mm,Y,m0,dict(cfg,arm=arm,fft_geometry=geom))
            elif arm=='cie_dual':
                # mature observed copy only adds callbacks, not numerical changes
                from research.research_max_struct_r6.src.mature_cie_observed import run_cie as observed
                result=observed(mm,Y,m0,cfg)
            else:result=run_factorial_cie(mm,Y,m0,dict(cfg,current_transform='scalar',material_transform='onehop'))
            results.append((result,record,dict(mode=mode,wall_seconds=time.perf_counter()-wall,
                cpu_seconds=time.process_time()-ccpu,status=result['status'],
                original_Phi=min(c['objective'] for c in result['candidates']),
                final_candidate_Phi=result['candidates'][-1]['objective'],
                physical_relative_residual=result['candidates'][-1]['physical_relative_residual'],
                cost=result['cost'])))
        dense,fft=results
        out['one_current_update'].append(dict(arm=arm,executions=[dense[2],fft[2]],
            material_relative_difference=rel(fft[0]['last_joint_m'],dense[0]['last_joint_m']),
            current_relative_difference=rel(fft[1][0]['J'],dense[1][0]['J']),
            coefficients_relative_difference=rel(fft[1][0]['coefficients'],dense[1][0]['coefficients'])))
    out.update(status='PASS' if all(r['material_relative_difference']<1e-8 and r['current_relative_difference']<1e-10 for r in out['one_current_update']) else 'FAIL',
        wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
        scope='This verifies new action/adjoint and one actual update/readout. It is not a long nonlinear trajectory or end-to-end throughput claim.',
        operator_check_cost=dict(cost))
    target=REPO/'research/research_max_struct_r6/audit/FFT_WIRING_RECEIPT.json'
    target.write_text(json.dumps(out,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else list(x),allow_nan=False))
    print(json.dumps({k:out[k] for k in ('status','wall_seconds','one_current_update')},default=lambda x:x.item() if isinstance(x,np.generic) else list(x)))

if __name__=='__main__':
    with threadpool_limits(limits=1): run()
