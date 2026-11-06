"""Actually execute original150/64 Austria geometry with matrix-free physics.

Separate output from the32-grid nine-arm comparison. Uses the same noisy150²
observations, exact paper source/receiver/frequency/rank/mf/beta. Grid cell
quadrature, real material/prior/finite ceiling remain declared adaptations.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import time,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r8.src.paper_physics import acquisition,MatrixFreeModel,austria
from research.research_max_struct_r8.src.paper_cie_fft import run_cie_fft
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,write

def run():
    case='austria_eps3p5_400MHz';folder=OUT/'original64'/case/'cie_single';folder.mkdir(parents=True,exist_ok=True)
    Y=np.load(OUT/case/'PUBLIC_DATA.npz')['Y'];t,c=time.perf_counter(),time.process_time()
    public=acquisition(64,400e6,dense=False);model=MatrixFreeModel(public);m0=np.zeros(model.n)
    stages=[dict(mf=mf,beta_material=b,beta_current=b,maxiter=300,coefficient_tol=tol**2) for mf,b,tol in [(5,6,.01),(7,2,.001),(10,.5,.001)]]
    config=dict(arm='cie_single',stages=stages,dpcs_rank=15,material_box=(0,20),miniter=2,sourcewise_cg=True,regularizer=lambda z:model.reg.evaluate(z)[:2],fft_geometry=model.paper_geometry,physical_sigma=np.linalg.norm(Y))
    try:
        result=run_cie_fft(model,Y,m0,config);m=result['m'];aux=m.copy();history=[]
        def fg(z):
            v,g=model.data_value_gradient(z,Y);history.append(dict(Phi=v,cpu_seconds=time.process_time()-c,complex_RHS=model.cost['actual_complex_RHS']));return v,g
        opt=minimize(fg,m,jac=True,method='L-BFGS-B',bounds=[(0,20)]*model.n,options=dict(maxiter=80,maxfun=320,gtol=1e-8,ftol=1e-12,maxls=30))
        final=opt.x if model.physical_score(opt.x,Y)<model.physical_score(m,Y) else m;state=model.state(final)
        np.savez_compressed(folder/'PREDICTION.npz',m=final,m_auxiliary=aux,Y_prediction=state['Y'])
        write(folder/'PREDICTION_SEAL.json',dict(status='SEALED_BEFORE_OFFLINE_GT',time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        write(folder/'RECEIPT.json',dict(status='COMPLETED',wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,physical_cost=dict(model.cost),auxiliary=result['execution_receipt'],relative_data_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),corrector_iterations=int(opt.nit),corrector_evaluations=int(opt.nfev),stop=str(opt.message),forward_grid=[150,150],inverse_grid=[64,64],adaptations=['equal-area disk-pulse square-cell quadrature','real [0,20] all independent pixel contrasts','common smooth+TV material prior','finite300 auxiliary iteration ceilings and80 corrector iterations']))
        write(folder/'HISTORY.json',dict(auxiliary=result['history'],corrector=history))
        # Offline scoring begins only after the seal exists.
        truth=austria(public['centers'],3.5);pred=final>1.25;support=truth>0
        write(folder/'OFFLINE_SCORES.json',dict(NRMSE=float(np.linalg.norm(final-truth)/np.linalg.norm(truth)),auxiliary_NRMSE=float(np.linalg.norm(aux-truth)/np.linalg.norm(truth)),IoU=float(np.sum(pred&support)/np.sum(pred|support))))
        print('COMPLETED',json.loads((folder/'RECEIPT.json').read_text()),flush=True)
    except Exception as e:
        write(folder/'FAILURE.json',dict(status='FAILED',error=repr(e),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,physical_cost=dict(model.cost)));raise

if __name__=='__main__':run()
