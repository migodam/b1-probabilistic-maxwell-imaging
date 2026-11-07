"""Generate the root-frozen final4 independent objects; no inverse or scoring."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
from datetime import datetime,timezone
import json,sys,time,shutil,argparse
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO));sys.path.insert(0,str(ROOT/'src'))


def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def main():
    wall,cpu=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
    cfg=json.loads((ROOT/'configs/twod_conditional_generator.json').read_text())
    freeze=ROOT/'delegated/twod/conditional_confirmation/FREEZE.json'
    if not args.execute or not freeze.is_file() or cfg['status']!='FROZEN_FOR_GENERATION':raise ValueError('Frozen rule and --execute required')
    if cfg['seed']!=2026110711 or len(cfg['cases'])!=4:raise ValueError('Frozen final4 seed/count mismatch')
    public=ROOT/'data_public/conditional_confirmation';private=ROOT/'sealed_truth/conditional_confirmation'
    receipt_path=ROOT/'delegated/twod/conditional_confirmation/GENERATION_RECEIPT.json'
    for p in [public,private,receipt_path]:
        if p.exists():raise FileExistsError(p)
    public.mkdir();private.mkdir();dump(private/'FROZEN_GENERATOR_CONFIG.json',cfg)
    receipt=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),objects_completed=0,
        inverse_runs=0,NN_calls=0,operator_builds=0,same_discretization=True,truth_access='GENERATOR_ONLY',costs=[])
    model=None;dump(receipt_path,receipt)
    try:
        import numpy as np
        from scipy.stats import chi2
        from threadpoolctl import threadpool_limits
        from twod_materials import draw_material
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        with threadpool_limits(limits=1):
            shutil.copy2(REPO/cfg['source_operator'],public/'operators.npz')
            with np.load(public/'operators.npz',allow_pickle=False) as z:ops={k:z[k] for k in z.files}
            model=MaxwellModel(ops);children=np.random.SeedSequence(cfg['seed']).spawn(4)
            records=[];meta=[];scores=[];tasks=[]
            for i,row in enumerate(cfg['cases']):
                t,c=time.perf_counter(),time.process_time();before=dict(model.cost)
                rng=np.random.default_rng(children[i]);m,parameters=draw_material(rng,row['family'],'high',row['high_chi'])
                clean=model.state(m)['Y'].copy();sigma=float(.01*np.linalg.norm(clean)/np.sqrt(clean.size))
                noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
                ident=row['case_id'];obs=public/(ident+'.npz')
                np.savez_compressed(obs,Y=clean+noise,sigma=sigma,noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2)),grid_shape=[24,24],material_bounds=[.02,20.])
                np.savez_compressed(private/(ident+'.npz'),m=m,Y_clean=clean)
                meta.append(dict(**row,parameters=parameters,spawn_key=list(children[i].spawn_key),realized_relative_noise=float(np.linalg.norm(noise)/np.linalg.norm(clean))))
                scores.append(dict(case_id=ident,amplitude=float(np.max(m)-.1),high_chi=row['high_chi'],family=row['family'],regime='high'))
                files=['research/research_max_struct_r10/data_public/conditional_confirmation/operators.npz',str(obs.relative_to(REPO))]
                records.append(dict(case_id=ident,public_files=files))
                task_id=ident+'__conditional_prefix_singleCIE'
                tasks.append(dict(task_id=task_id,case_id=ident,arm='cie_single',dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,input_base='.',public_files=files,
                    config='research/research_max_struct_r10/configs/twod_conditional_solver.json',experiment_scope='FRESH_CONDITIONAL_CONFIRMATION',
                    output_dir='research/research_max_struct_r10/RESULTS/conditional_confirmation/prefix/'+task_id))
                after=dict(model.cost);receipt['costs'].append(dict(case_id=ident,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,physical_cost={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)}));receipt['objects_completed']=i+1
            dump(private/'generator_metadata.json',dict(seed=cfg['seed'],records=meta));dump(private/'score_metadata.json',dict(records=scores,scope='OFFLINE_ONLY'))
            dump(public/'manifest.json',dict(path_base='REPOSITORY_ROOT',records=records,inverse_input_policy='Only public observations/operators; no generator labels/seed/contrast'))
            dump(ROOT/'configs/twod_conditional_prefix_tasks.json',dict(status='ROOT_FROZEN_FRESH_CONDITIONAL_CONFIRMATION',input_base='.',tasks=tasks))
            receipt['status']='GENERATED_FINAL4_SEPARATE_TRUTH'
    except BaseException as exc:
        receipt.update(status='FAILED_PRESERVED',error=repr(exc));raise
    finally:
        receipt.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost) if model else {},
            timing_scope='Generator entry/imports/operator copy/setup/new4 draws/forward/noise/serialization; nested object windows not additive',finished_utc=datetime.now(timezone.utc).isoformat())
        dump(receipt_path,receipt)
    print(json.dumps(dict(status=receipt['status'],objects=4,LU=receipt['physical_cost'].get('LU'),complex_RHS=receipt['physical_cost'].get('actual_complex_RHS'),CPU_seconds=receipt['CPU_seconds'])))


if __name__=='__main__':main()
