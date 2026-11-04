"""R7 separate data-generation process; solvers never import this file/truth."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time
from datetime import datetime,timezone
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import make_grid
from research.research_z_alg_r4_diversity.src.physics import Model

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(obj,f,indent=2)
def rel(p):return str(p.relative_to(REPO))

def geometry(side,params):
    # Physical analytic shapes, area averages over subpixels, no font/support prior.
    fine=side*8;axis=(np.arange(fine)+.5)*.96/fine-.48
    y,x=np.meshgrid(axis,axis,indexing='ij')
    x=x-params['dx'];y=y-params['dy'];t=params['theta'];scale=params['scale']
    X=(x*np.cos(t)+y*np.sin(t))/scale;Y=(-x*np.sin(t)+y*np.cos(t))/scale
    r=np.hypot(X,Y);angle=np.arctan2(Y,X);f=params['family']
    if f=='curved_bar':inside=(r>.205)&(r<.325)&(angle>params['opening'])&(angle<2.8)
    elif f=='eccentric_cavity':inside=(r<.345)&(np.hypot(X-.035,Y+.015)>.233)
    elif f=='slit_shell':inside=(r>.210)&(r<.345)&~((X>0)&(abs(Y)<params['slit']))
    elif f=='thin_winding':inside=((abs(Y-.12*np.sin(11*X))<.050)&(abs(X)<.335))|((abs(X+.25)<.045)&(Y<.16)&(Y>-.25))
    elif f=='multiple_islands':
        inside=(np.hypot(X+.22,Y+.13)<.095)|(np.hypot(X-.21,Y+.10)<.12)|(np.hypot(X+.015,Y-.24)<.105)
    elif f=='coupled_capsules':
        inside=(np.hypot(np.maximum(abs(X+.153)-.04,0),Y+.025)<.093)|(np.hypot(np.maximum(abs(X-.153)-.025,0),Y-.025)<.095)
    else:raise ValueError(f)
    return inside.reshape(side,8,side,8).mean((1,3))

def main():
    start,cpu=time.perf_counter(),time.process_time();rng=np.random.default_rng(2026110407)
    families=['curved_bar','eccentric_cavity','slit_shell','thin_winding','multiple_islands','coupled_capsules']
    cases=[];mapping=[];costs=[];ops_cache={};total_lu=total_rhs=0
    # Five prespecified acquisition/contrast strata, six different objects each.
    conditions=[(2.5,.001,12,'full',24),(6.,.01,12,'full',24),(12.,.001,6,'full',24),(18.,.001,12,'half',24),(12.,.005,12,'full',36)]
    for variant,(amp,noise,ns,aperture,gen_side) in enumerate(conditions):
        for family in families:
            i=len(cases);ident=f'r7_obj_{i:02d}';p=dict(family=family,dx=float(rng.uniform(-.035,.035)),dy=float(rng.uniform(-.035,.035)),theta=float(rng.uniform(-.7,.7)),scale=float(rng.uniform(.88,1.04)),opening=float(rng.uniform(-.8,-.2)),slit=float(rng.uniform(.035,.055)))
            def operators(side):
                key=(side,ns,aperture)
                if key not in ops_cache:
                    op=make_grid(side,sources=12,receivers=48)
                    inds=np.arange(0,12,2) if ns==6 else np.arange(12)
                    if aperture=='half':inds=np.arange(12) # narrower external source arc, same legal source count
                    if aperture=='half':
                        from research.research_z_alg_r4_diversity.src.physics import build_operators
                        angles=np.linspace(-np.pi/2,np.pi/2,12)+.071
                        op['sources']=1.4*np.c_[np.cos(angles),np.sin(angles)]
                        op.update(build_operators(op['centers'],op['sources'],op['receivers'],op['radius'],op['k']))
                    else:
                        op['sources']=op['sources'][inds];op['E0']=op['E0'][:,inds]
                    op['E']=op['E0'];op['U0']=op['E0'];op['coords']=op['centers']
                    path=ROOT/'data_public/operators'/f'n{side}_s{ns}_{aperture}.npz';path.parent.mkdir(exist_ok=True,parents=True)
                    np.savez_compressed(path,**op);ops_cache[key]=(op,path)
                return ops_cache[key]
            inv,invpath=operators(24);gen,genpath=operators(gen_side)
            alpha=geometry(gen_side,p);m=.1+amp*alpha.ravel();model=Model(gen)
            t,c=time.perf_counter(),time.process_time()
            state=model.forward(m,label='R7_SEPARATE_DATA_GENERATION');clean=state['y']
            sigma=noise*np.linalg.norm(clean)/np.sqrt(clean.size)
            e=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
            Y=clean+e
            obs=ROOT/'data_public'/(ident+'.npz')
            np.savez_compressed(obs,Y=Y,sigma=sigma,noise_radius=np.sqrt(chi2.ppf(.99,2*Y.size)/2),material_bounds=[.02,20.],grid_shape=[24,24])
            tr=ROOT/'sealed_truth'/(ident+'.npz')
            coarse=geometry(24,p)
            np.savez_compressed(tr,m=.1+amp*coarse.ravel(),alpha=coarse,mask=coarse>=.5,Y_clean=clean,generation_material=m,noise=e)
            cases.append(dict(case_id=ident,object_id=ident,cohort='R7_NEW_CONFIRMATION',condition=variant,amplitude=amp,noise_relative_complex_RMS=noise,sources=ns,aperture=aperture,inversion_grid=[24,24],generation_grid=[gen_side,gen_side],inverse_crime=gen_side==24,public_files=[rel(invpath),rel(obs)]))
            mapping.append(dict(case_id=ident,parameters=p,truth_path=rel(tr),generation_operators=rel(genpath)))
            total_lu+=1;total_rhs+=ns;costs.append(dict(case_id=ident,dimension=gen_side**2,LU=1,complex_RHS=ns,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c))
    dump(ROOT/'data_public/cases.json',dict(status='NEW_30_UNSCORED_PUBLIC',records=cases))
    dump(ROOT/'sealed_truth/mapping_PRIVATE.json',dict(seed=2026110407,records=mapping,release='Only after all primary predictions sealed'))
    arms=['cie_single','cie_dual','global_dual','a_material','a_plus','physical_lbfgs']
    tasks=[dict(task_id=r['case_id']+'__'+arm,case_id=r['case_id'],cohort=r['cohort'],arm=arm,public_files=r['public_files'],output_dir='research/research_max_struct_r7/RESULTS/primary/'+r['case_id']+'__'+arm) for r in cases for arm in arms]
    manifest=dict(config='research/research_max_struct_r7/configs/primary.json',tasks=tasks)
    dump(ROOT/'configs/primary_manifest.json',manifest)
    dump(ROOT/'configs/first_four_manifest.json',dict(manifest,tasks=[t for t in tasks if t['arm'] in arms[:4]]))
    dump(ROOT/'audit/GENERATION_RECEIPT.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),objects=30,independent_materials=30,repeated_measurements_counted_as_materials=False,physical_counts=dict(LU=total_lu,complex_RHS=total_rhs),per_case=costs,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,GT_access='Separate generator only; inverse never opens sealed_truth',conditions=conditions))
    print(json.dumps(dict(status='GENERATED_30_PUBLIC_OBSERVATIONS_GT_HELD',LU=total_lu,complex_RHS=total_rhs)))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
