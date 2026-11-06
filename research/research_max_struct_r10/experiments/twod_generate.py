"""Generate frozen R10 observations with inherited MaxwellModel, no inverse/NN."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import argparse,json,sys,time,shutil
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO));sys.path.insert(0,str(ROOT/'src'))


def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    wall,cpu=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=ROOT/'configs/twod_generator.json')
    ap.add_argument('--execute',action='store_true')
    args=ap.parse_args();cfg=json.loads(args.config.read_text())
    if not args.execute or cfg['status']!='FROZEN_FOR_GENERATION':
        ap.error('Root must freeze generation, then explicitly pass --execute')
    if (cfg['generation_grid'],cfg['inverse_grid'],cfg['sources'],cfg['receivers'])!=(24,24,12,48):
        raise ValueError('Only the inherited 24-grid/full12x48 contract is supported')
    if len(cfg['cases'])!=24 or cfg['material_background']!=.1 or cfg['material_box']!=[.02,20.]:
        raise ValueError('Frozen material contract mismatch')
    public_dir=ROOT/'data_public/twod';truth_dir=ROOT/'sealed_truth/twod'
    receipt_path=ROOT/'delegated/twod/GENERATION_RECEIPT.json'
    for p in (public_dir,truth_dir,receipt_path):
        if p.exists():raise FileExistsError(p)
    public_dir.mkdir(parents=True);truth_dir.mkdir(parents=True)
    receipt_path.parent.mkdir(parents=True,exist_ok=True)
    receipt=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        objects_requested=24,objects_completed=0,inverse_runs=0,NN_calls=0,
        backend='Inherited MaxwellModel/R4 disk-pulse operators, complex128',
        same_discretization=True,operator_builds=0,truth_access='GENERATOR_ONLY',
        costs=[],shared_operator_setup=[],monetary_cost=None)
    dump(receipt_path,receipt);model=None;records=[]
    try:
        import numpy as np
        from scipy.stats import chi2
        from threadpoolctl import threadpool_limits
        from twod_materials import draw_material,FAMILIES
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        expected=[(f,r) for f in FAMILIES for r in ('low','high','mixed')]
        if [(r['family'],r['regime']) for r in cfg['cases']]!=expected:
            raise ValueError('Frozen family/regime order differs')
        dump(truth_dir/'FROZEN_GENERATOR_CONFIG.json',cfg)
        with threadpool_limits(limits=1):
            c,t=time.process_time(),time.perf_counter()
            source=REPO/cfg['operator_source'];destination=REPO/cfg['operator_destination']
            shutil.copy2(source,destination)
            with np.load(destination,allow_pickle=False) as z:
                public={key:z[key].copy() for key in z.files}
            if (public['G'].shape,public['E0'].shape,public['S'].shape)!=((576,576),(576,12),(48,576)):
                raise ValueError('Public operator dimensions differ')
            if {'truth','m_true','Y_clean','family','seed'}&set(public):raise ValueError('Private operator payload')
            model=MaxwellModel(public)
            receipt['shared_operator_setup'].append(dict(CPU_seconds=time.process_time()-c,
                wall_seconds=time.perf_counter()-t,operator_source=cfg['operator_source'],
                operator_destination=cfg['operator_destination'],cost=dict(model.cost),
                original_operator_generation='REUSED_EXTERNAL_NOT_RERUN'))
            children=np.random.SeedSequence(cfg['seed']).spawn(24)
            private_meta=[];score_meta=[]
            for index,row in enumerate(cfg['cases']):
                c,t=time.process_time(),time.perf_counter();before=dict(model.cost)
                ident=row['case_id'];rng=np.random.default_rng(children[index])
                m,parameters=draw_material(rng,row['family'],row['regime'],row['high_chi'])
                state=model.state(m);clean=state['Y'].copy()
                sigma=float(cfg['relative_complex_RMS_noise']*np.linalg.norm(clean)/np.sqrt(clean.size))
                noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
                obs=public_dir/(ident+'.npz')
                np.savez_compressed(obs,Y=clean+noise,sigma=sigma,
                    noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2)),
                    grid_shape=np.array([24,24]),material_bounds=np.array([.02,20.]))
                np.savez_compressed(truth_dir/(ident+'.npz'),m=m,Y_clean=clean)
                private_meta.append(dict(**row,parameters=parameters,
                    independent_child_spawn_key=list(children[index].spawn_key),
                    realized_relative_noise=float(np.linalg.norm(noise)/np.linalg.norm(clean))))
                score_meta.append(dict(case_id=ident,amplitude=float(np.max(m)-.1),
                    amplitude_semantics='max truth chi minus background .1; offline normalization only'))
                records.append(dict(case_id=ident,object_id=ident,
                    public_files=[cfg['operator_destination'],str(obs.relative_to(REPO))]))
                after=dict(model.cost)
                receipt['costs'].append(dict(case_id=ident,CPU_seconds=time.process_time()-c,
                    wall_seconds=time.perf_counter()-t,
                    physical_cost={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)}))
                receipt['objects_completed']=len(records)
                dump(receipt_path,receipt)
            dump(truth_dir/'generator_metadata.json',dict(seed=cfg['seed'],records=private_meta))
            dump(truth_dir/'score_metadata.json',dict(records=score_meta,scope='OFFLINE_ONLY'))
            dump(public_dir/'manifest.json',dict(status='GENERATED_PUBLIC_OBSERVATIONS',
                path_base='REPOSITORY_ROOT',records=records,
                inverse_input_policy='Only public_files; no geometry/regime/seed/amplitude/GT'))
            tasks=[]
            for row in records:
                for arm in ('cie_single','cie_dual','global_dual','a_plus','a_plus_sourcewise'):
                    task_id=row['case_id']+'__'+arm+'__noise_ball_99'
                    tasks.append(dict(task_id=task_id,case_id=row['case_id'],arm=arm,
                        dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,
                        public_files=row['public_files'],input_base='.',
                        config='research/research_max_struct_r10/configs/twod_solver.json',
                        output_dir='research/research_max_struct_r10/RESULTS/twod/'+task_id))
            dump(ROOT/'configs/twod_tasks.json',dict(status='PUBLIC_TASKS_NO_GT',input_base='.',tasks=tasks))
            receipt['status']='GENERATED_SEPARATE_TRUTH'
    except BaseException as exc:
        receipt.update(status='FAILED_PARTIAL_GENERATION_PRESERVED',error=repr(exc));raise
    finally:
        receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            physical_cost=dict(model.cost) if model else {},
            timing_scope='Whole generation process after entry: imports/copy/setup/draws/forward/noise/serialization; nested object/setup costs not additive')
        dump(receipt_path,receipt)
    print(json.dumps(dict(status=receipt['status'],objects=24,LU=receipt['physical_cost'].get('LU'),
        complex_RHS=receipt['physical_cost'].get('actual_complex_RHS'),CPU_seconds=receipt['CPU_seconds'])))


if __name__=='__main__':main()
