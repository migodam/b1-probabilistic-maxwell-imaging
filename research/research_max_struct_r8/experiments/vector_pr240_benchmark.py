from pathlib import Path
import sys,os,json,time,resource
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r8.src.vector_pr3d import VectorPR3D
DATA=ROOT/'RESULTS/vector3d';OUT=DATA/'pr240';OUT.mkdir(exist_ok=True);METHODS=['LSIE_sourcewise','CIE_sourcewise','Dual_sourcewise','CIE_global','A_global','Aplus_global']
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False))
def main():
    t=time.perf_counter();cpu=time.process_time();dump(OUT/'PR_FROZEN_CONFIG.json',dict(status='FROZEN_BEFORE_RUN',outer=240,inner=12,corrector=300,methods=METHODS,current_update='one sourcewise/global Polak-Ribiere exact-quadratic-line step per cycle; no exact current LS',Aplus='real scalar-isotropic only; state-stationary arc candidates + exact quadratic prior acceptance; complex NOT_AVAILABLE',uses_existing_public_data=True))
    for name in ['four_cubes']:
        target=OUT/name;target.mkdir(exist_ok=True);Y=np.load(DATA/name/'public_observations.npz')['Y'];m=VectorPR3D(6,lossy=name.startswith('lossy'))
        for kind in METHODS:
            if kind=='Aplus_global' and m.lossy:dump(target/(kind+'_NOT_AVAILABLE.json'),dict(status='NOT_AVAILABLE',reason='complex material is a two-real-dimensional domain; real Mobius arc invalid'));continue
            if (target/(kind+'_receipt.json')).exists():continue
            arrays,r=m.optimize_pr(Y,kind,outer=240);np.savez_compressed(target/(kind+'_prediction.npz'),**arrays);dump(target/(kind+'_receipt.json'),r);print(json.dumps(dict(target=name,method=kind,cpu_s=r['cpu_s'],residual=r['physical_data_residual'])),flush=True)
    rows=[]
    for name in ['four_cubes']:
        target=OUT/name;truth=np.load(DATA/name/'truth_for_scoring.npz')['chi']
        for kind in METHODS:
            p=target/(kind+'_prediction.npz')
            if not p.exists():continue
            a=np.load(p);r=json.loads((target/(kind+'_receipt.json')).read_text());rows.append(dict(target=name,method=kind,raw_NRMSE=float(la.norm(a['raw']-truth)/la.norm(truth)),final_NRMSE=float(la.norm(a['prediction']-truth)/la.norm(truth)),residual=r['physical_data_residual'],cpu_s=r['cpu_s'],counts=r['counts'],stop=r['optimizer'][-1]))
    dump(OUT/'PR_SCORES.json',rows);dump(OUT/'PR_COMPLETE.json',dict(status='COMPLETED',rows=len(rows),wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
