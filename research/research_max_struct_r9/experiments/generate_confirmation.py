"""Prepare six independent R9 observations with inherited R8 material/physics.

No work occurs on import. Parent must freeze the generator contract, then call
with --execute. This file does not run any inverse solver, feature or NN.
"""
from pathlib import Path
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse,json,sys,time
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]


def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def relative(path):
    return str(Path(path).resolve().relative_to(REPO))


def main():
    wall,cpu=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=ROOT/'configs/confirmation_generator.json')
    ap.add_argument('--execute',action='store_true')
    args=ap.parse_args()
    if not args.execute:
        ap.error('Preparation only. Root must freeze the contract and explicitly use --execute to generate.')
    config_path=args.config.resolve()
    cfg=json.loads(config_path.read_text())
    if cfg['status']!='FROZEN_FOR_GENERATION':
        raise ValueError('Generator contract is not frozen; no data were generated')
    expected=[('ellipse',2.5,.001,'full'),('disks',12.,.01,'half'),
              ('annulus',18.,.03,'full'),('comb',12.,.03,'half'),
              ('star',18.,.01,'full'),('perforated',12.,.03,'half')]
    cases=cfg['cases']
    if cfg['seed']!=2026110609 or len(cases)!=6:
        raise ValueError('Frozen six-object seed/count mismatch')
    if [(r['family'],r['amplitude'],r['relative_complex_RMS_noise'],r['aperture']) for r in cases]!=expected:
        raise ValueError('Frozen family/amplitude/noise/aperture schedule mismatch')
    if [r['case_id'] for r in cases]!=[f'r9_confirmation_{i:04d}' for i in range(6)]:
        raise ValueError('Opaque case IDs differ from the prepared public manifest')
    if (cfg['generation_grid'],cfg['inverse_grid'],cfg['sources'],cfg['receivers'])!=(24,24,12,48):
        raise ValueError('Only the frozen same-discretization 24-grid/12-source/48-receiver contract is supported')
    if cfg['material_box']!=[.02,20.] or cfg['material_background']!=.1 or not cfg['same_discretization']:
        raise ValueError('Material/discretization contract mismatch')
    public_dir=ROOT/'data_public/confirmation'
    truth_dir=ROOT/'sealed_truth/confirmation'
    receipt_path=ROOT/'audit/generation_confirmation.json'
    for p in (public_dir,truth_dir,receipt_path):
        if p.exists():raise FileExistsError(p)
    for ref in cfg['operators'].values():
        p=REPO/ref
        if not p.is_file():raise FileNotFoundError(p)
    public_dir.mkdir(parents=True);truth_dir.mkdir(parents=True)
    receipt_path.parent.mkdir(parents=True,exist_ok=True)
    dump(truth_dir/'FROZEN_GENERATOR_CONFIG.json',cfg)
    receipt=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        config=relative(config_path),seed=cfg['seed'],objects_requested=6,
        backend='Inherited MaxwellModel over copied R8 disk-pulse G/E/S; exp(-iwt); radius.45h; no operator rebuild',
        generation_grid=24,inverse_grid=24,same_discretization=True,
        feature_calls=0,NN_calls=0,inverse_runs=0,operator_builds=0,
        truth_access='GENERATOR_ONLY; separate sealed_truth; inverse manifest contains no truth or score metadata',
        costs=[],shared_operator_setup=[],monetary_cost=None)
    dump(receipt_path,receipt)
    models={};score_meta=[];offline_meta=[];records=[]
    try:
        import numpy as np
        from scipy.stats import chi2
        from threadpoolctl import threadpool_limits
        if str(REPO) not in sys.path:sys.path.insert(0,str(REPO))
        # Import only the existing parameter/material functions. R8 main and
        # its feature-generation path are never called; no FeatureOperator is
        # instantiated here.
        from research.research_max_struct_r8.experiments.generate_warm_data import parameters,material
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        with threadpool_limits(limits=1):
            for aperture,ref in cfg['operators'].items():
                t,c=time.perf_counter(),time.process_time()
                with np.load(REPO/ref,allow_pickle=False) as z:
                    forbidden={'truth','gt','m_true','m_truth','mask','family','seed','alpha','y_clean','yclean'}
                    if forbidden & {k.lower() for k in z.files}:raise ValueError('Private payload in public operator')
                    public={k:z[k].copy() for k in z.files}
                if public['G'].shape!=(576,576) or public['E0'].shape!=(576,12) or public['S'].shape!=(48,576):
                    raise ValueError('Copied operator dimensions differ from frozen acquisition')
                if tuple(public['grid_shape'])!=(24,24):raise ValueError('Operator grid mismatch')
                model=MaxwellModel(public);models[aperture]=model
                receipt['shared_operator_setup'].append(dict(aperture=aperture,operator=ref,
                    CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,
                    physical_cost=dict(model.cost),original_operator_generation='REUSED_R8_EXTERNAL_NOT_RERUN'))
            children=np.random.SeedSequence(cfg['seed']).spawn(6)
            for i,row in enumerate(cases):
                t,c=time.perf_counter(),time.process_time()
                rng=np.random.default_rng(children[i]);p=parameters(rng,row['family'])
                m=material(24,p,float(row['amplitude']))
                if not np.isfinite(m).all() or np.min(m)<.02 or np.max(m)>20.:
                    raise ValueError('Material outside frozen box; no reroll permitted')
                model=models[row['aperture']];before=dict(model.cost)
                try:
                    state=model.state(m);clean=state['Y'].copy()
                    if clean.shape!=(48,12):raise ValueError('Unexpected measurement dimensions')
                    sigma=float(row['relative_complex_RMS_noise']*np.linalg.norm(clean)/np.sqrt(clean.size))
                    noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
                    Y=clean+noise
                    ident=row['case_id'];obs=public_dir/(ident+'.npz')
                    # Exactly these five fields are allowed in inverse observations.
                    np.savez_compressed(obs,Y=Y,sigma=sigma,
                        noise_radius=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2)),
                        grid_shape=np.array([24,24]),material_bounds=np.array([.02,20.]))
                    np.savez_compressed(truth_dir/(ident+'.npz'),m=m,
                        mask=m>.1+.5*float(row['amplitude']),Y_clean=clean)
                    score_meta.append(dict(case_id=ident,amplitude=float(row['amplitude'])))
                    offline_meta.append(dict(case_id=ident,family=row['family'],parameters=p,
                        independent_child_spawn_key=list(children[i].spawn_key),
                        relative_complex_RMS_noise=row['relative_complex_RMS_noise']))
                    records.append(dict(case_id=ident,object_id=ident,
                        public_files=[cfg['operators'][row['aperture']],relative(obs)]))
                finally:
                    after=dict(model.cost)
                    delta={k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before)}
                    receipt['costs'].append(dict(case_id=row['case_id'],dimension=576,
                        physical_cost=delta,CPU_seconds=time.process_time()-c,
                        wall_seconds=time.perf_counter()-t))
            dump(truth_dir/'score_metadata.json',dict(records=score_meta,
                scope='OFFLINE_ONLY; amplitude denominator; not an inverse input'))
            dump(truth_dir/'generator_metadata.json',dict(seed=cfg['seed'],records=offline_meta,
                scope='GENERATOR_OFFLINE_ONLY; never pass to an inverse runner'))
            manifest=dict(status='GENERATED_PUBLIC_OBSERVATIONS',path_base='REPOSITORY_ROOT',
                inverse_input_policy='Only public_files; no GT/family/amplitude/seed/features',records=records)
            dump(public_dir/'manifest.json',manifest)
            receipt.update(status='GENERATED_SEPARATE_TRUTH',objects_completed=len(records),
                LU=sum(v['physical_cost'].get('LU',0) for v in receipt['costs']),
                complex_RHS=sum(v['physical_cost'].get('actual_complex_RHS',0) for v in receipt['costs']),
                model_cost_by_aperture={key:dict(model.cost) for key,model in models.items()},
                public_manifest=relative(public_dir/'manifest.json'),
                public_observation_keys=['Y','sigma','noise_radius','grid_shape','material_bounds'],
                scoring_metadata='sealed_truth/confirmation/score_metadata.json',
                confirmation_scope='Six diverse fixed independent realizations; not a population/statistical claim')
    except BaseException as exc:
        receipt.update(status='FAILED_PARTIAL_GENERATION_PRESERVED',error=repr(exc),objects_completed=len(records),
            model_cost_by_aperture={key:dict(model.cost) for key,model in models.items()})
        raise
    finally:
        receipt.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            LU=sum(model.cost.get('LU',0) for model in models.values()),
            complex_RHS=sum(model.cost.get('actual_complex_RHS',0) for model in models.values()),
            finished_utc=datetime.now(timezone.utc).isoformat(),
            timing_scope='Imports/public operator loading/model setup/material draws/forward/noise/array serialization; per-object and shared setup timers nested, not additive to total; final receipt write excluded')
        dump(receipt_path,receipt)


if __name__=='__main__':main()
