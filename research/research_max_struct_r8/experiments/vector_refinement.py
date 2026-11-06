from pathlib import Path
import sys,os,json,time,resource
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r8.src.vector_pr3d import VectorPR3D
from research.research_max_struct_r8.src.vector_cie3d import VectorCIE3D,objects
OUT=ROOT/'RESULTS/vector3d/refinement8';OUT.mkdir(exist_ok=True,parents=True)
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False))
METHODS=['LSIE_sourcewise','Dual_sourcewise','Aplus_global','direct_physical']
def main():
    t=time.perf_counter();cpu=time.process_time();dump(OUT/'FROZEN_CONFIG.json',dict(status='FROZEN_BEFORE_RUN',inverse_grid=8,synthetic_grid=10,target='four_cubes',domain_lambda=.6,sources=60,receivers=60,noise=.1,seed=8036,current_cycles=48,material_inner=12,common_long_corrector=300,direct_initial=192,methods=METHODS,physical_callback_cache='exact np.array_equal material and identical Y object; reuse paidvalue/J/E; no newLU; same300cap; engineeringattemptcostseparate',reason='necessary independent mesh refinement of15.16%coarsecrossgridmisfit; unchangedanalyticobject/acquisition/chi-box/smoothingformula; no parametersearch'))
    # Reuse alreadysealed independent10³data; no repeatedsimulation/GTonlineaccess.
    m=VectorPR3D(8);m.cache_physical_values=True;Y=np.load(OUT/'public_observations.npz')['Y'];sim=json.loads((OUT/'SIMULATION.json').read_text());sim['reused_previously_sealed_data_from']='../refinement8_uncached_attempt';sim['new_simulation_cost']=dict(physical_LU=0,source_RHS=0);sim['new_online_operator_setup']=m.setup;dump(OUT/'SIMULATION.json',sim)
    for kind in METHODS:
        if (OUT/(kind+'_receipt.json')).exists():continue
        arrays,r=m.optimize(Y,kind) if kind=='direct_physical' else m.optimize_pr(Y,kind);np.savez_compressed(OUT/(kind+'_prediction.npz'),**arrays);dump(OUT/(kind+'_receipt.json'),r);print(json.dumps(dict(method=kind,cpu_s=r['cpu_s'],residual=r['physical_data_residual'])),flush=True)
    truth=np.load(OUT/'truth_for_scoring.npz')['chi'];rows=[]
    for kind in METHODS:
        a=np.load(OUT/(kind+'_prediction.npz'));r=json.loads((OUT/(kind+'_receipt.json')).read_text());rows.append(dict(target='four_cubes_grid8',method=kind,raw_NRMSE=float(la.norm(a['raw']-truth)/la.norm(truth)),final_NRMSE=float(la.norm(a['prediction']-truth)/la.norm(truth)),residual=r['physical_data_residual'],cpu_s=r['cpu_s'],counts=r['counts'],stop=r['optimizer'][-1]))
    dump(OUT/'SCORES.json',rows);dump(OUT/'COMPLETE.json',dict(status='COMPLETED',rows=len(rows),wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
