"""Literal cubic-shell+rods correction; prior exposedcuboid results retained."""
from pathlib import Path
import sys,os,json,time,resource
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r8.src.vector_pr3d import VectorPR3D
from research.research_max_struct_r8.src.vector_cie3d import VectorCIE3D,box_fraction
OUT=ROOT/'RESULTS/vector3d/hollow_cubic_shell_rods';OUT.mkdir(exist_ok=True,parents=True)
METHODS=['LSIE_sourcewise','CIE_global','Aplus_global','direct_physical']
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False))
def cubic(m):
    x=box_fraction(m.points,m.h,[0,.05,0],[.15]*3)-box_fraction(m.points,m.h,[0,.05,0],[.09]*3)
    for xx in [-.2,.2]:x+=box_fraction(m.points,m.h,[xx,-.2,0],[.035,.035,.18])
    return 1.5*x

def main():
    t=time.perf_counter();cpu=time.process_time();dump(OUT/'FROZEN_CONFIG.json',dict(status='FROZEN_BEFORE_RUN',inverse_grid=6,synthetic_grid=8,domain_lambda=.6,outer_cube_side=.3,inner_cube_side=.18,shell_center=[0,.05,0],rods_centers=[[-.2,-.2,0],[.2,-.2,0]],rod_halfwidths=[.035,.035,.18],epsilon=2.5,sources=60,receivers=60,noise=.1,seed=8038,coefficient_coordinates='original projectedIFFTgamma, noSVDwhitening',PRcycles=48,inner=12,corrector=300,methods=METHODS,reason='Previous Austria adaptation wasa hollowcuboid; preserveitandrunliteralcubic-shell geometryinstead ofrelabelingpreviousarrays'))
    f=VectorCIE3D(8);Y,J,E=f.physical(f.encode(cubic(f)));fc=dict(counts=f.counts,setup=f.setup);del f,J,E;m=VectorPR3D(6);m.cache_physical_values=True;m.whiten_apcs=False;truth=cubic(m);pred,J,E=m.physical(m.encode(truth));noise=np.random.default_rng(8038).normal(size=Y.shape)+1j*np.random.default_rng(8039).normal(size=Y.shape);noise*=.1*la.norm(Y)/la.norm(noise);dump(OUT/'SIMULATION.json',dict(fine=fc,coarse_counts=m.counts.copy(),coarse_setup=m.setup,crossgrid_forward_discrepancy=float(la.norm(pred-Y)/la.norm(Y))));np.savez_compressed(OUT/'public_observations.npz',Y=Y+noise,clean_Y=Y);np.savez_compressed(OUT/'truth_for_scoring.npz',chi=truth,points=m.points);del truth,pred,J,E;Y+=noise
    for kind in METHODS:
        if (OUT/(kind+'_receipt.json')).exists():continue
        a,r=m.optimize(Y,kind) if kind=='direct_physical' else m.optimize_pr(Y,kind);np.savez_compressed(OUT/(kind+'_prediction.npz'),**a);dump(OUT/(kind+'_receipt.json'),r);print(json.dumps(dict(method=kind,cpu_s=r['cpu_s'],residual=r['physical_data_residual'])),flush=True)
    truth=np.load(OUT/'truth_for_scoring.npz')['chi'];rows=[]
    for kind in METHODS:
        a=np.load(OUT/(kind+'_prediction.npz'));r=json.loads((OUT/(kind+'_receipt.json')).read_text());rows.append(dict(target='hollow_cubic_shell_rods',method=kind,raw_NRMSE=float(la.norm(a['raw']-truth)/la.norm(truth)),final_NRMSE=float(la.norm(a['prediction']-truth)/la.norm(truth)),residual=r['physical_data_residual'],cpu_s=r['cpu_s'],counts=r['counts'],stop=r['optimizer'][-1]))
    dump(OUT/'SCORES.json',rows);dump(OUT/'COMPLETE.json',dict(status='COMPLETED',rows=len(rows),wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
