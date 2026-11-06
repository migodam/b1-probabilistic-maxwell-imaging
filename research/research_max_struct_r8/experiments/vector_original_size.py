from pathlib import Path
import sys,os,json,time,resource,traceback
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r8.src.vector_fft3d import VectorFFT3D
OUT=ROOT/'RESULTS/vector3d/original_size';OUT.mkdir(exist_ok=True,parents=True)
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False))
def main():
    t=time.perf_counter();cpu=time.process_time();models=[];stage='setup';dump(OUT/'FROZEN_CONFIG.json',dict(status='FROZEN_BEFORE_RUN',domain_lambda=3,inverse_grid=12,synthetic_grid=16,sources=60,receivers=60,complex_observations=10800,original_four_cube_side_lambda=.8,epsilon=[2,3,4,5],noise=.1,seed=8034,method='direct full-vector physical L-BFGS-B',initial_chi=.25,box=[0,4],first_maxiter=20,common_long_corrector_maxiter=30,GMRES_rtol=2e-6,GMRES_restart=60,GMRES_maxiter=8,FFT='R1 full dyadic CM+RR; A18_C padded convolution form; no dense G',limitations='12³/16³ much coarser than paper30³/60³, insufficient internal wavelength at eps5; no continuum accuracy claim; direct baseline only original-size attempt'))
    try:
        fine=VectorFFT3D(16);models.append(fine);stage='fine16_forward';truth=fine.four_cubes();Y,J=fine.physical(truth);noise=np.random.default_rng(8034).normal(size=Y.shape)+1j*np.random.default_rng(8035).normal(size=Y.shape);noise*=.1*la.norm(Y)/la.norm(noise);np.savez_compressed(OUT/'public_observations.npz',Y=Y+noise,clean_Y=Y);dump(OUT/'FINE_FORWARD.json',dict(counts=fine.counts,setup_wall=fine.setup_wall,max_linear_residual=fine.max_residual,elapsed_s=time.perf_counter()-t));print('fine16 forward complete',time.perf_counter()-t,flush=True);del truth,J
        m=VectorFFT3D(12);models.append(m);stage='crossgrid12_forward';truth=m.four_cubes();clean12,J=m.physical(truth);np.savez_compressed(OUT/'truth_for_scoring.npz',chi=truth,points=m.points);dump(OUT/'CROSSGRID.json',dict(relative_data_difference=float(la.norm(clean12-Y)/la.norm(Y)),counts=m.counts.copy(),setup_wall=m.setup_wall,max_linear_residual=m.max_residual));del truth,J,clean12,fine;Y+=noise;z=np.full(m.N,.25);best=z.copy();bestphi=np.inf;trials=[];curve=[]
        def fun(x):return m.physical(x,Y)
        def callback(x):
            nonlocal best,bestphi
            phi,df,J=m.physical(x,Y,False);trials.append(x.copy());curve.append(dict(phi=phi,residual=np.sqrt(df),elapsed_s=time.perf_counter()-t,counts=m.counts.copy(),max_linear_residual=m.max_residual));
            if phi<bestphi:best=x.copy();bestphi=phi
            np.savez_compressed(OUT/'checkpoint_prediction.npz',prediction=best,trials=np.array(trials));dump(OUT/'PROGRESS.json',dict(stage=stage,iterations=len(trials),last=curve[-1]));print(json.dumps(dict(stage=stage,iteration=len(trials),phi=phi,elapsed_s=time.perf_counter()-t,fft_matvec=m.counts['fft_matvec'])),flush=True)
        logs=[]
        for stage,budget in [('direct_initial',20),('common_physical_corrector',30)]:
            res=minimize(fun,z,jac=True,method='L-BFGS-B',bounds=[(0,4)]*m.N,callback=callback,options=dict(maxiter=budget,ftol=1e-10,gtol=1e-7,maxls=12));callback(res.x);z=best.copy();logs.append(dict(stage=stage,nit=res.nit,nfev=res.nfev,success=bool(res.success),message=str(res.message)))
        np.savez_compressed(OUT/'prediction.npz',prediction=best,trials=np.array(trials));dump(OUT/'RECEIPT.json',dict(status='COMPLETED',optimizer=logs,counts=m.counts,max_linear_residual=m.max_residual,curve=curve));truth=np.load(OUT/'truth_for_scoring.npz')['chi'];dump(OUT/'SCORE.json',dict(material_NRMSE=float(la.norm(best-truth)/la.norm(truth)),physical_phi=bestphi))
    except Exception as ex:
        dump(OUT/'FAILED.json',dict(status='FAILED',stage=stage,error=repr(ex),traceback=traceback.format_exc(),model_costs=[dict(ns=m.ns,counts=m.counts,max_linear_residual=m.max_residual,failures=m.failures,setup_wall=m.setup_wall) for m in models]));raise
    finally:dump(OUT/'COSTS.json',dict(wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,model_costs=[dict(ns=m.ns,counts=m.counts,max_linear_residual=m.max_residual,setup_wall=m.setup_wall) for m in models]))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
