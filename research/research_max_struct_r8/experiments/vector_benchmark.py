from pathlib import Path
import sys,os,json,time,resource,traceback
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r8.src.vector_cie3d import VectorCIE3D,objects
OUT=ROOT/'RESULTS/vector3d';OUT.mkdir(exist_ok=True,parents=True)
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False))
METHODS=['LSIE_TSOM','CIE_TSOM','CIE_dualbeta_TSOM','A_TSOM','direct_physical']
def main():
    t=time.perf_counter();c=time.process_time();config=dict(status='FROZEN_BEFORE_RUN',inverse_grid=6,synthetic_grid=8,domain_width_lambda=.6,paper_width_lambda=3,adaptation_linear_scale=.2,source_count=60,receiver_count=60,complex_data=10800,noise_relative_Frobenius=.1,noise_seed=8033,material_real_box=[0,4],material_imag_box=[0,1.2],outer=16,inner=12,direct_maxiter=192,common_long_corrector=300,methods=METHODS,limits='scaled electrical size; cross-grid DDA not continuum; exact box volume fraction chi mixing; no original 30³/60³ reproduction; explicit dense propagation; exact reduced current solve differs paper CG',truth_policy='generate then delete truth variables; seal every prediction before scorer reads truth')
    dump(OUT/'FROZEN_CONFIG.json',config)
    fine=VectorCIE3D(8);truths=objects(fine);coarse=VectorCIE3D(6);ct=objects(coarse);rng=np.random.default_rng(8033);sim=[]
    for name,x in truths.items():
        fine.lossy=np.iscomplexobj(x);fine.reset();Y,J,E=fine.physical(fine.encode(x));noise=rng.normal(size=Y.shape)+1j*rng.normal(size=Y.shape);noise*=.1*la.norm(Y)/la.norm(noise);target=OUT/name;target.mkdir(exist_ok=True)
        np.savez_compressed(target/'public_observations.npz',Y=Y+noise,clean_Y=Y,noise=noise,receivers=fine.receivers,polarizations=fine.polarizations)
        np.savez_compressed(target/'truth_for_scoring.npz',chi=ct[name],fine_chi=x,points=coarse.points)
        # Cross-grid same analytic geometry: independent assembly and independent factorization.
        coarse.lossy=fine.lossy;coarse.reset();pred,*_=coarse.physical(coarse.encode(ct[name]));sim.append(dict(target=name,fine_counts=fine.counts.copy(),coarse_counts=coarse.counts.copy(),crossgrid_forward_discrepancy=float(la.norm(pred-Y)/la.norm(Y)),noise=float(la.norm(noise)/la.norm(Y))))
    dump(OUT/'SIMULATION.json',dict(fine_setup=fine.setup,coarse_setup=coarse.setup,rows=sim));del fine,truths,ct,x,coarse,pred,J,E
    # Necessary complex/nonroot derivative check; not a frozen formula sweep.
    m=VectorCIE3D(6,lossy=True);z=m.encode(np.linspace(.1,1,m.N)+.2j);Y=np.load(OUT/'lossy_coated_cube/public_observations.npz')['Y'];_,J,E=m.physical(z);J=J*(1.03+.02j);v=rng.normal(size=z.size);v/=la.norm(v);checks={}
    for kind in ['physical','LSIE_TSOM','CIE_TSOM','A_TSOM']:
        f=(lambda q:m.physical(q,Y)) if kind=='physical' else (lambda q,k=kind:m.material(q,J,k,1.))
        val,g=f(z);hh=1e-6;fd=(f(z+hh*v)[0]-f(z-hh*v)[0])/(2*hh);checks[kind]=dict(fd=fd,analytic=float(g@v),relative_error=float(abs(fd-g@v)/max(abs(fd),abs(g@v),1e-12)))
    dump(OUT/'DERIVATIVES.json',dict(checks=checks,counts=m.counts));assert max(x['relative_error'] for x in checks.values())<1e-5;del m
    for name in ['four_cubes','hollow_cube_rods','lossy_coated_cube']:
        target=OUT/name;Y=np.load(target/'public_observations.npz')['Y'];m=VectorCIE3D(6,lossy=name.startswith('lossy'))
        for kind in METHODS:
            if (target/(kind+'_receipt.json')).exists():continue
            try:
                arrays,receipt=m.optimize(Y,kind);np.savez_compressed(target/(kind+'_prediction.npz'),**arrays);dump(target/(kind+'_receipt.json'),receipt);print(json.dumps(dict(target=name,method=kind,residual=receipt['physical_data_residual'],cpu_s=receipt['cpu_s'],LU=receipt['counts']['physical_LU'])),flush=True)
            except Exception as ex:dump(target/(kind+'_FAILED.json'),dict(error=repr(ex),traceback=traceback.format_exc(),counts=m.counts));raise
        del m
    rows=[]
    for name in ['four_cubes','hollow_cube_rods','lossy_coated_cube']:
        target=OUT/name;truth=np.load(target/'truth_for_scoring.npz')['chi']
        for kind in METHODS:
            arr=np.load(target/(kind+'_prediction.npz'));r=json.loads((target/(kind+'_receipt.json')).read_text());rows.append(dict(target=name,method=kind,raw_NRMSE=float(la.norm(arr['raw']-truth)/la.norm(truth)),final_NRMSE=float(la.norm(arr['prediction']-truth)/la.norm(truth)),real_NRMSE=float(la.norm(arr['prediction'].real-truth.real)/la.norm(truth.real)),imag_NRMSE=float(la.norm(arr['prediction'].imag-truth.imag)/max(la.norm(truth.imag),1e-30)),residual=r['physical_data_residual'],cpu_s=r['cpu_s'],wall_s=r['wall_s'],counts=r['counts'],stop=r['optimizer'][-1]))
    dump(OUT/'SCORES.json',rows);dump(OUT/'COMPLETE.json',dict(status='COMPLETED',rows=len(rows),wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
