"""Independent physical generator for supervised training or frozen public test cases.

The inverse runner never imports this module. Test truth is written separately.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import argparse,json,time,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import make_grid
from research.research_z_alg_r4_diversity.src.physics import Model,build_operators
from research.research_max_struct_r8.src.warm_features import FeatureOperator,FEATURE_NAMES

FAMILIES=['ellipse','rectangles','disks','annulus','nested','cross','random_fourier','smooth_blobs']
OOD=['austria_variant','comb','star','perforated']

def parameters(rng,family):
    return dict(family=family,dx=float(rng.uniform(-.06,.06)),dy=float(rng.uniform(-.06,.06)),
                theta=float(rng.uniform(-np.pi,np.pi)),a=float(rng.uniform(.12,.29)),
                b=float(rng.uniform(.07,.20)),centers=rng.uniform(-.27,.27,(4,2)).tolist(),
                radii=rng.uniform(.055,.12,4).tolist(),weights=rng.uniform(.4,1,4).tolist(),
                phases=rng.uniform(0,2*np.pi,6).tolist(),scale=float(rng.uniform(.82,1.05)))

def material(side,p,amp):
    sub=4;z=(np.arange(side*sub)+.5)*.96/(side*sub)-.48;y,x=np.meshgrid(z,z,indexing='ij')
    x-=p['dx'];y-=p['dy'];co,si=np.cos(p['theta']),np.sin(p['theta'])
    X=(co*x+si*y)/p['scale'];Y=(-si*x+co*y)/p['scale'];r=np.hypot(X,Y)
    f=p['family'];a,b=p['a'],p['b'];v=np.zeros_like(X)
    if f=='ellipse':v=((X/a)**2+(Y/b)**2<1).astype(float)
    elif f=='rectangles':
        for (cx,cy),rad,w in zip(p['centers'][:3],p['radii'],p['weights']):
            v=np.maximum(v,w*((abs(X-cx)<rad*1.4)&(abs(Y-cy)<rad)))
    elif f=='disks':
        for (cx,cy),rad,w in zip(p['centers'],p['radii'],p['weights']):v=np.maximum(v,w*(np.hypot(X-cx,Y-cy)<rad))
    elif f=='annulus':v=((r<.29)&(r>.16+.2*b)).astype(float)
    elif f=='nested':v=.45*((abs(X)<.30)&(abs(Y)<.26))+.55*((X/a)**2+(Y/b)**2<1)
    elif f=='cross':v=(((abs(X)<.30)&(abs(Y)<.055))|((abs(Y)<.30)&(abs(X)<.055))).astype(float)
    elif f=='random_fourier':
        field=sum(np.cos(2*np.pi*(u*X+v*Y)+t)/(1+u*u+v*v)**.5 for (u,v),t in zip([(1,0),(0,1),(1,1),(2,1),(1,-2),(3,0)],p['phases']))
        v=((field>.25)&(r<.36)).astype(float)
    elif f=='smooth_blobs':
        for (cx,cy),rad,w in zip(p['centers'],p['radii'],p['weights']):v+=w*np.exp(-((X-cx)**2+(Y-cy)**2)/(2*rad**2))
        v=np.minimum(v,1.)*(r<.40)
    elif f=='austria_variant':v=(((r>.16)&(r<.28)&(Y<.12))|(np.hypot(X-.18,Y-.24)<.08)|(np.hypot(X+.18,Y-.24)<.08)).astype(float)
    elif f=='comb':v=(((abs(Y+.23)<.045)&(abs(X)<.30))|((abs(X+.24)<.035)&(abs(Y)<.25))|((abs(X)<.035)&(abs(Y)<.25))|((abs(X-.24)<.035)&(abs(Y)<.25))).astype(float)
    elif f=='star':v=(r<.23+.065*np.cos(5*np.arctan2(Y,X))).astype(float)
    elif f=='perforated':
        v=((abs(X)<.31)&(abs(Y)<.30)).astype(float)
        for cx,cy in [(-.15,-.12),(.15,-.12),(0,.15)]:v*=np.hypot(X-cx,Y-cy)>.075
    else:raise ValueError(f)
    return .1+amp*v.reshape(side,sub,side,sub).mean((1,3)).ravel()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--split',choices=['train','validation','test'],required=True);ap.add_argument('--count',type=int);args=ap.parse_args()
    seeds=dict(train=2026110601,validation=2026110602,test=2026110613);rng=np.random.default_rng(seeds[args.split])
    count=args.count or dict(train=1536,validation=192,test=42)[args.split]
    target=ROOT/('data_training' if args.split!='test' else 'data_public/warm_test');target.mkdir(parents=True,exist_ok=True)
    receipt_path=ROOT/'audit'/('generation_'+args.split+'.json')
    if receipt_path.exists():raise FileExistsError(receipt_path)
    begin,cpu=time.perf_counter(),time.process_time();cache={};features=[];labels=[];rows=[];costs=[];featurecost=[]
    def operator(side,ns,aperture):
        key=(side,ns,aperture)
        if key not in cache:
            p=make_grid(side,sources=ns,receivers=48)
            if aperture==.5:
                ang=np.linspace(-np.pi/2,np.pi/2,ns)+.071;p['sources']=1.4*np.c_[np.cos(ang),np.sin(ang)]
                p.update(build_operators(p['centers'],p['sources'],p['receivers'],p['radius'],p['k']))
            p.update(source_aperture_fraction=aperture,coords=p['centers']);p['E']=p['E0']
            out=ROOT/'data_public/operators'/f'warm_n{side}_s{ns}_a{aperture}.npz';out.parent.mkdir(exist_ok=True,parents=True)
            if not out.exists():np.savez_compressed(out,**p)
            cache[key]=(p,Model(p),FeatureOperator(p),str(out.relative_to(REPO)))
        return cache[key]
    truthdir=ROOT/'sealed_truth/warm_test'
    if args.split=='test':truthdir.mkdir(exist_ok=True,parents=True)
    for i in range(count):
        family=(FAMILIES[i%8] if i<24 else OOD[(i-24)%4] if i<36 else FAMILIES[(i-36)%8]) if args.split=='test' else FAMILIES[i%8]
        p=parameters(rng,family);amp=float(rng.choice([2.5,6,12,18]));ns=int(rng.choice([6,12]));aperture=float(rng.choice([.5,1.]));noise=float(rng.choice([.001,.01,.03]))
        gen_side=36 if args.split=='test' and i>=36 else 24
        inv,_,feat,opfile=operator(24,ns,aperture);gen,model,_,_=operator(gen_side,ns,aperture)
        m=material(gen_side,p,amp);t,c=time.perf_counter(),time.process_time();Yclean=model.forward(m,label='SUPERVISED_TRAIN' if args.split!='test' else 'SEALED_TEST_GENERATOR')['y']
        sigma=noise*np.linalg.norm(Yclean)/np.sqrt(Yclean.size);Y=Yclean+sigma/np.sqrt(2)*(rng.normal(size=Yclean.shape)+1j*rng.normal(size=Yclean.shape))
        X,analytic,fc=feat(Y);features.append(X);featurecost.append(fc)
        ident=f'r8_{args.split}_{i:04d}';coarse=material(24,p,amp)
        rows.append(dict(case_id=ident,object_id=ident,split=args.split,condition=dict(amplitude=amp,sources=ns,aperture_fraction=aperture,relative_complex_RMS_noise=noise,generation_grid=gen_side,inverse_grid=24),
                         cohort=('IN_DISTRIBUTION' if i<24 else 'HELDOUT_FAMILY' if i<36 else 'CROSS_GRID') if args.split=='test' else args.split))
        if args.split!='test':labels.append(coarse.reshape(24,24));rows[-1].update(family=family,parameters=p)
        else:
            obs=target/(ident+'.npz');np.savez_compressed(obs,Y=Y,sigma=sigma,noise_radius=np.sqrt(chi2.ppf(.99,2*Y.size)/2),grid_shape=[24,24],material_bounds=[.02,20.])
            np.savez_compressed(truthdir/(ident+'.npz'),m=coarse,mask=coarse>.1+.5*amp,Y_clean=Yclean,generation_material=m)
            np.savez_compressed(target/(ident+'_features.npz'),features=X,born=analytic['born'],dpcs=analytic['dpcs'])
            rows[-1].update(public_files=[opfile,str(obs.relative_to(REPO))],features_file=str((target/(ident+'_features.npz')).relative_to(REPO)))
        costs.append(dict(case_id=ident,dimension=gen_side**2,LU=1,complex_RHS=ns,cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t))
        if (i+1)%128==0:print(json.dumps(dict(split=args.split,completed=i+1,total=count,wall=time.perf_counter()-begin)),flush=True)
    if args.split!='test':np.savez_compressed(target/(args.split+'.npz'),X=np.array(features),m=np.array(labels,dtype=np.float32))
    (target/(args.split+'_manifest.json')).write_text(json.dumps(dict(seed=seeds[args.split],records=rows),indent=2))
    receipt=dict(status='GENERATED',split=args.split,seed=seeds[args.split],objects=count,LU=count,complex_RHS=sum(x['complex_RHS'] for x in costs),per_object=costs,
                 feature_names=FEATURE_NAMES,feature_cost=featurecost,operator_feature_setup=[dict(key=list(k),**v[2].setup) for k,v in cache.items()],
                 wall_seconds=time.perf_counter()-begin,cpu_seconds=time.process_time()-cpu,created_utc=datetime.now(timezone.utc).isoformat(),
                 truth_access='SUPERVISED_TRAINING_LABELS' if args.split!='test' else 'SEPARATE_GENERATOR_ONLY',backend='Inherited R7 disk-pulse scalar TM, exp(-iwt), radius .45h; finite discrete model, not continuum qualification')
    receipt_path.write_text(json.dumps(receipt,indent=2));print(json.dumps({k:receipt[k] for k in ('status','split','objects','LU','complex_RHS','wall_seconds','cpu_seconds')}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
