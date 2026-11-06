"""Domain-augmented frozenCNN repair on public eps4 benchmark.

No GT access, no retraining, no candidate screening by initial physical loss.
Material-only initialization; legal TSOM state remains JD+APCS with zeroAPCS.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time
import numpy as np
import torch
from scipy.optimize import minimize
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,write
from research.research_max_struct_r8.src.paper_physics import acquisition,make_model
from research.research_max_struct_r8.src.paper_cie_fft import run_cie_fft
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
from research.research_max_struct_r8.src.warm_features import FeatureOperator
from research.research_max_struct_r8.src.warm_network import load_predictor

def run():
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    case='austria_eps4_400MHz';data=np.load(OUT/case/'PUBLIC_DATA.npz');Y=data['Y'];public=acquisition(32,400e6)
    c,t=time.process_time(),time.perf_counter();features=FeatureOperator(public);x,_,feature_receipt=features(Y)
    x[-2:]*=.48 # fixed public paperhalfwidth1m
    x=np.concatenate((x,np.broadcast_to(np.array([2*public['k']/(2*np.pi),1.],np.float32)[:,None,None],(2,32,32))),axis=0)
    weights=ROOT/'RESULTS/training/domain_augmented/weights.pt';predict,checkpoint=load_predictor(weights)
    m0=np.clip(predict(x)[0].reshape(-1),0,20);feature_receipt.update(setup=features.setup,total_cpu_seconds=time.process_time()-c,total_wall_seconds=time.perf_counter()-t,weights=str(weights),coordinates='x/1m,y/1m, fixed public paperDoI scaling',qualification='public development after domain augmentation on256 nonAustria training/64validation targets; noAustriaGT intraining; notblindconfirmation',material_only=True,source_conditioning='20/12 source count; aperture1')
    folder=OUT/case/'NN_AUG_INITIALIZATION';folder.mkdir(parents=True,exist_ok=True);np.savez_compressed(folder/'MATERIAL.npz',m=m0);write(folder/'RECEIPT.json',feature_receipt)
    for base in ('cie_dual','a_plus','a_plus_sourcewise'):
        arm='aug_nn_'+base;actual='a_plus' if base=='a_plus_sourcewise' else base;dest=OUT/case/arm;dest.mkdir(parents=True,exist_ok=True);model=make_model(public);wall,cpu=time.perf_counter(),time.process_time()
        cfg=dict(arm=actual,stages=[dict(mf=mf,beta_material=b,beta_current=8,maxiter=300,coefficient_tol=tol**2) for mf,b,tol in [(5,8,.01),(8,2,.001),(10,.5,.001)]],dpcs_rank=15,material_box=(0,20),miniter=2,sourcewise_cg=base in ('cie_dual','a_plus_sourcewise'),regularizer=lambda z:model.reg.evaluate(z)[:2],fft_geometry=dict(coords=public['centers'],radius=public['radius'],k=public['k']),physical_sigma=np.linalg.norm(Y),aplus_regularizer=model.reg,aplus_sweeps=2,aplus_arc_candidates=True,cache_material_value_gradient=True)
        if actual=='a_plus':cfg['material_readout']=mobius_material_readout
        result=run_cie_fft(model,Y,m0,cfg);m=result['m'];history=[]
        def fg(z):
            v,g=model.data_value_gradient(z,Y);history.append(dict(Phi=v,CPU_seconds=time.process_time()-cpu));return v,g
        opt=minimize(fg,m,jac=True,method='L-BFGS-B',bounds=[(0,20)]*model.n,options=dict(maxiter=80,maxfun=320,gtol=1e-8,ftol=1e-12,maxls=30))
        final=opt.x if model.physical_score(opt.x,Y)<model.physical_score(m,Y) else m;st=model.state(final)
        np.savez_compressed(dest/'PREDICTION.npz',m=final,m_auxiliary=m,Y_prediction=st['Y']);write(dest/'RECEIPT.json',dict(status='COMPLETED',arm=arm,wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,physical_cost=dict(model.cost),auxiliary=result['execution_receipt'],relative_data_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),initial_original_Phi=result['candidates'][0]['objective'],GT_loaded=False,feature_inference_cost_shared=feature_receipt))
        write(dest/'HISTORY.json',dict(auxiliary=result['history'],corrector=history));write(dest/'PREDICTION_SEAL.json',dict(status='SEALED_BEFORE_OFFLINE_GT',GT_loaded=False,time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())));print('NN_AUGMENTED',arm,time.perf_counter()-wall,flush=True)

if __name__=='__main__':run()
