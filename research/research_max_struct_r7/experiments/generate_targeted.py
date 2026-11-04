"""Separate targeted data generator; online inverse does not import this file."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,shutil
from datetime import datetime,timezone
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.generate_benchmark import geometry
from research.research_z_alg_r4_diversity.src.physics import Model

def dump(p,o):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(o,f,indent=2)
def rel(p):return str(p.relative_to(REPO))
def main():
    start,cpu=time.perf_counter(),time.process_time()
    cfg=json.loads((ROOT/'configs/targeted_frozen.json').read_text());rng=np.random.default_rng(cfg['seed'])
    cases=[];mapping=[];costs=[];tasks=[]
    for cond in cfg['families']:
        opath=ROOT/'data_public/operators'/f"n24_s{cond['sources']}_{cond['aperture']}.npz"
        with np.load(opath,allow_pickle=False) as d:op={k:d[k] for k in d.files}
        model=Model(op)
        for _ in range(cfg['per_family']):
            ident=f'r7_target_{len(cases):02d}';p={'family':cond['family']}
            for key in cfg['draw_order']:p[key]=float(rng.uniform(*cfg['geometry_ranges'][key]))
            alpha=geometry(24,p);m=.1+cond['amplitude']*alpha.ravel()
            t,c=time.perf_counter(),time.process_time();clean=model.forward(m,label='R7_TARGETED_SEPARATE_DATA_GENERATION')['y']
            sigma=cond['noise']*np.linalg.norm(clean)/np.sqrt(clean.size)
            e=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));Y=clean+e
            obs=ROOT/'data_public/targeted'/(ident+'.npz');obs.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(obs,Y=Y,sigma=sigma,noise_radius=np.sqrt(chi2.ppf(.99,2*Y.size)/2),material_bounds=[.02,20.],grid_shape=[24,24])
            tr=ROOT/'sealed_truth/targeted'/(ident+'.npz');tr.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(tr,m=m,alpha=alpha,mask=alpha>=.5,Y_clean=clean,generation_material=m,noise=e)
            case=dict(case_id=ident,object_id=ident,cohort='R7_TARGETED_SECONDARY',family=cond['family'],amplitude=cond['amplitude'],noise_relative_complex_RMS=cond['noise'],sources=cond['sources'],aperture=cond['aperture'],inversion_grid=[24,24],generation_grid=[24,24],inverse_crime=True,public_files=[rel(opath),rel(obs)])
            cases.append(case);mapping.append(dict(case_id=ident,parameters=p,truth_path=rel(tr),generation_operators=rel(opath)))
            costs.append(dict(case_id=ident,dimension=576,LU=1,complex_RHS=cond['sources'],wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c))
            for arm in cfg['arms']:
                tasks.append(dict(task_id=ident+'__'+arm,case_id=ident,cohort=case['cohort'],arm=arm,config=rel(ROOT/cfg['arm_configs'][arm]),initial_material=.1,public_files=case['public_files'],output_dir=rel(ROOT/'RESULTS/targeted'/(ident+'__'+arm))))
    dump(ROOT/'data_public/targeted/cases.json',dict(status='TARGETED20_UNSCORED_PUBLIC',records=cases))
    dump(ROOT/'sealed_truth/targeted/mapping_PRIVATE.json',dict(seed=cfg['seed'],records=mapping,release='Only after all60 sealed and parent TARGETED_SEALED'))
    dump(ROOT/'configs/targeted_manifest.json',dict(config=rel(ROOT/'configs/primary.json'),tasks=tasks))
    dump(ROOT/'audit/targeted/GENERATION_RECEIPT.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),objects=20,independent_materials=20,GT_access='Separate generator only; inverse never opens sealed_truth',physical_counts=dict(LU=20,complex_RHS=sum(r['complex_RHS'] for r in costs)),per_case=costs,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,operator_reuse=True,source_modified=False))
    freeze=ROOT/'audit/targeted/source_freeze';freeze.mkdir(parents=True,exist_ok=False)
    for path in [ROOT/'TARGETED_PROTOCOL.md',ROOT/'configs/targeted_frozen.json',ROOT/'configs/targeted_manifest.json',ROOT/'configs/primary.json',ROOT/'configs/aplus_final.json',Path(__file__),ROOT/'experiments/run_case.py',ROOT/'experiments/run_queue.py',ROOT/'src/benchmark_solver.py',ROOT/'src/aplus_readout.py',ROOT/'experiments/generate_benchmark.py']:
        shutil.copy2(path,freeze/str(path.relative_to(ROOT)).replace('/','__'))
    print(json.dumps(dict(status='GENERATED20_TARGETED_UNSCORED',tasks=len(tasks),LU=20,complex_RHS=180)))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
