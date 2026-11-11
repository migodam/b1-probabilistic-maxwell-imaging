"""Generate only root-frozen public observations, with separate generator GT.

No inversions. No rejection/resampling on data or inverse performance.
The confirmation24 draft refuses execution until root freezes its JSON.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import argparse,json,sys,time,shutil
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent

def dump(path,obj):
    def convert(x):
        import numpy as np
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False,default=convert)+'\n')

def material(rng,family,regime,high,xy,low=2.5):
    """Four explicit distributions; geometry is generator-only context."""
    import numpy as np
    from scipy.ndimage import gaussian_filter
    x,y=xy[:,0],xy[:,1]
    theta=float(rng.uniform(-np.pi,np.pi));cx,cy=rng.uniform(-.055,.055,2)
    u=(x-cx)*np.cos(theta)+(y-cy)*np.sin(theta)
    v=-(x-cx)*np.sin(theta)+(y-cy)*np.cos(theta)
    parameters=dict(rotation=theta,center=[float(cx),float(cy)])
    if family=='levelsets':
        # Signed analytic field, fixed threshold0 and smooth enclosing ellipse.
        modes=np.array([[1,0],[0,1],[1,1],[2,-1],[1,2],[2,2]])
        amplitudes=rng.uniform(.25,1.,len(modes))/(1+np.sum(modes*modes,axis=1))**.5
        phases=rng.uniform(-np.pi,np.pi,len(modes))
        field=sum(a*np.cos(2*np.pi*(q[0]*u+q[1]*v)+p) for a,q,p in zip(amplitudes,modes,phases))
        outer=(u/.39)**2+(v/.34)**2<1
        support=(field>0)&outer
        phase=(u+.35*v)>0
        parameters.update(modes=modes,amplitudes=amplitudes,phases=phases,threshold=0.,envelope_axes=[.39,.34])
    elif family=='branches':
        # Curved three-branch junction; distance to sampled piecewise segments.
        points=[(-.30,-.24),(0.,.02),(-.28,.29),(.28,.27),(.30,-.20)]
        jitter=rng.uniform(-.035,.035,(5,2));points=np.asarray(points)+jitter
        pairs=((0,1),(1,2),(1,3),(1,4));width=float(rng.uniform(.028,.046))
        bend=rng.uniform(-.055,.055,len(pairs));distance=np.full(len(x),np.inf)
        for j,(a,b) in enumerate(pairs):
            t=np.linspace(0.,1.,17);delta=points[b]-points[a]
            normal=np.array([-delta[1],delta[0]])/np.linalg.norm(delta)
            curve=(1-t[:,None])*points[a]+t[:,None]*points[b]+bend[j]*np.sin(np.pi*t)[:,None]*normal
            for first,last in zip(curve[:-1],curve[1:]):
                d=last-first
                f=np.clip(((u-first[0])*d[0]+(v-first[1])*d[1])/np.dot(d,d),0.,1.)
                dist=np.hypot(u-first[0]-f*d[0],v-first[1]-f*d[1])
                distance=np.minimum(distance,dist)
        support=distance<width;phase=v>.015
        parameters.update(control_points=points,bends=bend,half_width=width,curve_segments_per_branch=16)
    elif family=='mixed_islands':
        # Five deliberately separated elliptic components; mixed phase alternates.
        angles=np.arange(5)*2*np.pi/5+rng.uniform(-.15,.15,5)
        radii=rng.uniform(.21,.27,5);axes=rng.uniform(.065,.10,(5,2))
        support=np.zeros(len(x),bool);phase=np.zeros(len(x),bool)
        for j in range(5):
            xx=u-radii[j]*np.cos(angles[j]);yy=v-radii[j]*np.sin(angles[j])
            inside=(xx/axes[j,0])**2+(yy/axes[j,1])**2<1
            support|=inside
            if j%2==0:phase|=inside
        parameters.update(angles=angles,radii=radii,ellipse_axes=axes,high_phase_component_indices=[0,2,4])
    elif family=='continuous_random':
        # No binary support model is supplied to inversions. Compact smooth
        # envelope times positive Gaussian-smoothed random field; exact peak high.
        raw=rng.normal(size=(24,24));scale=float(rng.uniform(1.4,2.2))
        field=gaussian_filter(raw,scale,mode='reflect').ravel()
        field=(field-field.min())/max(float(np.ptp(field)),1e-300)
        envelope=np.maximum(0.,1-(u/.40)**2-(v/.35)**2)**2
        field=(.25+.75*field)*envelope
        if regime=='mixed':
            phase_field=.5+.5*np.tanh((u-.06*np.sin(8*v))/.05)
            field*=low/high+(1-low/high)*phase_field
        field/=max(float(field.max()),1e-300)
        m=.1+(float(high)-.1)*field
        if not np.any(field>0):raise ValueError('Empty fixed continuous draw; no reroll')
        parameters.update(gaussian_sigma_pixels=scale,envelope_axes=[.40,.35],
            raw_gaussian=raw,continuous_mixed_rule='low/high+(1-low/high)*sigmoid partition; normalize final field to max1')
        return m,parameters
    else:raise ValueError('Unknown frozen family '+family)
    if not np.any(support):raise ValueError('Empty fixed structured draw; no reroll')
    m=np.full(len(x),.1)
    if regime=='high':m[support]=float(high)
    elif regime=='mixed':
        h=support&phase;l=support&~phase
        if not h.any() or not l.any():raise ValueError('Degenerate fixed mixed draw; no reroll')
        m[h]=float(high);m[l]=float(low)
    else:raise ValueError(regime)
    return m,parameters

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True);ap.add_argument('--execute',action='store_true')
    ap.add_argument('--legacy-repo',type=Path,default=Path(os.environ.get('R12_LEGACY_REPO',str(ROOT.parents[3]))))
    ap.add_argument('--vendor-path');ap.add_argument('--cpu',type=int);a=ap.parse_args()
    if a.vendor_path:sys.path.insert(0,a.vendor_path)
    legacy=a.legacy_repo.resolve();sys.path.insert(0,str(legacy))
    if a.cpu is not None:
        import psutil
        psutil.Process().cpu_affinity([a.cpu])
        if os.name=='nt':psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    cfg=json.loads(a.config.read_text());stage=cfg['stage']
    if not a.execute or cfg['status']!='FROZEN_FOR_GENERATION':
        raise ValueError('Frozen generator config plus explicit --execute required; confirm24 draft intentionally blocked')
    if stage not in ('dev6','confirm24'):raise ValueError('Unknown stage')
    if stage=='confirm24' and not cfg.get('root_core_freeze_receipt'):
        raise ValueError('Confirmation additionally requires explicit root_core_freeze_receipt in frozen generator config')
    expected=6 if stage=='dev6' else 24
    if len(cfg['cases'])!=expected:raise ValueError('Case count differs')
    if (cfg['generation_grid'],cfg['inverse_grid'],cfg['sources'],cfg['receivers'],cfg['material_background'],cfg['material_box'],cfg['relative_complex_RMS_noise'])!=(24,24,12,48,.1,[.02,20.],.01):raise ValueError('Physics contract differs')
    pub=ROOT/'data_public'/stage;private=ROOT/'sealed_truth'/stage
    receipt_path=ROOT/(stage.upper()+'_GENERATION_RECEIPT.json')
    for path in (pub,private,receipt_path):
        if path.exists():raise FileExistsError(path)
    pub.mkdir(parents=True);private.mkdir(parents=True)
    receipt=dict(status='RUNNING',stage=stage,requested=expected,completed=0,
        started_utc=datetime.now(timezone.utc).isoformat(),truth_access='GENERATOR_ONLY',inversions=0,
        scientific_tuning='NONE',seed_policy='Independent deterministic SeedSequence children; distinct geometry/noise substreams; no reroll',objects=[])
    model=None;wall=time.perf_counter()
    dump(receipt_path,receipt)
    try:
        import numpy as np
        from scipy.stats import chi2
        from threadpoolctl import threadpool_limits,threadpool_info
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        with threadpool_limits(limits=1):
            source=(legacy/cfg['operator_source']).resolve();dest=pub/'operators.npz'
            shutil.copy2(source,dest)
            shutil.copy2(Path(__file__),private/'GENERATOR_SOURCE_FROZEN.py')
            with np.load(dest,allow_pickle=False) as z:operators={k:z[k].copy() for k in z.files}
            forbidden={'truth','m_true','Y_clean','family','seed','mask','support'}
            if forbidden&set(operators):raise ValueError('Private payload in inherited operator')
            if (operators['G'].shape,operators['E0'].shape,operators['S'].shape)!=((576,576),(576,12),(48,576)):raise ValueError('Operator dimension differs')
            model=MaxwellModel(operators);xy=np.asarray(operators.get('coords',operators['centers']))
            children=np.random.SeedSequence(cfg['seed']).spawn(expected);records=[];metadata=[];scoremeta=[]
            dump(private/'FROZEN_GENERATOR_CONFIG.json',cfg)
            stat=source.stat()
            receipt.update(threads=threadpool_info(),operator_source=str(source),operator_source_bytes=stat.st_size,
                operator_source_mtime_ns=stat.st_mtime_ns,new_operator_builds=0)
            for index,row in enumerate(cfg['cases']):
                before=dict(model.cost);cpu,t=time.process_time(),time.perf_counter()
                geometry_seed,noise_seed=children[index].spawn(2)
                m,parameters=material(np.random.default_rng(geometry_seed),row['family'],row['regime'],row['high_chi'],xy,cfg['low_chi'])
                if m.shape!=(576,) or not np.isfinite(m).all() or m.min()<.02 or m.max()>20.:raise ValueError('Illegal material; no reroll')
                clean=model.state(m)['Y'].copy()
                sigma=float(.01*np.linalg.norm(clean)/np.sqrt(clean.size));rng=np.random.default_rng(noise_seed)
                noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
                obs=pub/(row['case_id']+'.npz')
                np.savez_compressed(obs,Y=clean+noise,sigma=sigma,
                    noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2)),noise_model=np.array('iid_circular_complex'),
                    grid_shape=np.array([24,24]),material_bounds=np.array([.02,20.]))
                np.savez_compressed(private/(row['case_id']+'.npz'),m=m,Y_clean=clean)
                metadata.append(dict(**row,parameters=parameters,geometry_spawn_key=list(geometry_seed.spawn_key),noise_spawn_key=list(noise_seed.spawn_key),realized_noise_relative=float(np.linalg.norm(noise)/np.linalg.norm(clean))))
                scoremeta.append(dict(case_id=row['case_id'],amplitude=float(m.max()-.1)))
                records.append(dict(case_id=row['case_id'],public_files=[str(dest),str(obs)]))
                after=dict(model.cost)
                receipt['objects'].append(dict(case_id=row['case_id'],CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-t,
                    physical_cost={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)}))
                receipt['completed']=len(records);dump(receipt_path,receipt)
            dump(private/'generator_metadata.json',dict(seed=cfg['seed'],records=metadata))
            dump(private/'score_metadata.json',dict(records=scoremeta,scope='OFFLINE_ONLY'))
            dump(pub/'manifest.json',dict(status='PUBLIC_ONLY',records=records))
            baseline=[];dr=[]
            for row in records:
                for method in ('cie_single','cie_dual','global_dual','a_plus','fullG19','DR'):
                    ident=row['case_id']+'__'+method+'__240CPU'
                    task=dict(task_id=ident,case_id=row['case_id'],arm='a_plus' if method=='fullG19' else method,
                        algorithm_label=method,dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,
                        public_files=row['public_files'],input_base=str(legacy),
                        config=str(ROOT/cfg.get('solver_config','DEV6_BASELINE_SOLVER.json')),output_dir=str(ROOT/'RESULTS'/stage/ident),
                        experiment_scope='R12_'+stage.upper())
                    if method=='fullG19':task['material_readout']=dict(module_file=str(legacy/'research/research_max_struct_r11/src/material_readouts.py'),
                        callable='full_material_readout',prepare='prepare_full',label='fullG19')
                    if method=='DR':
                        task.update(status='PENDING_ROOT_CORE_IMPLEMENTATION',runner=None);dr.append(task)
                    else:baseline.append(task)
            dump(pub/'baseline_tasks.json',dict(status='PUBLIC_TASKS_NO_GT',tasks=baseline))
            dump(pub/'DR_TASKS_PENDING_CORE.json',dict(status='DO_NOT_EXECUTE_PENDING_ROOT_CORE',tasks=dr))
            receipt['status']='GENERATED_PUBLIC_AND_SEPARATE_TRUTH'
    except BaseException as exc:
        receipt.update(status='FAILED_PARTIAL_GENERATION_PRESERVED',error=repr(exc));raise
    finally:
        receipt.update(whole_child_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,
            physical_model_cost=dict(model.cost) if model else {},
            cost_scope='Absolute child process_time including imports, setup/copy/generation/forward/noise/serialization; external lifetime CPU remains authoritative')
        dump(receipt_path,receipt)
    print(json.dumps(dict(status=receipt['status'],completed=receipt['completed'],whole_child_CPU_seconds=receipt['whole_child_CPU_seconds'])),flush=True)

if __name__=='__main__':main()
