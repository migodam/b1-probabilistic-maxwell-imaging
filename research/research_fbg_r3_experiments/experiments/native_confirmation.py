"""Separate freeze/generate/infer actions for native-affine bounded confirmation."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,sys,time,zipfile,shutil,secrets
import numpy as np
from concurrent.futures import ProcessPoolExecutor,as_completed
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE/'experiments'))
from src.physics import load_model


def freeze():
    plan=dict(status='METHODS_FROZEN_BEFORE_NEW_GENERATION',utc=datetime.now(timezone.utc).isoformat(),
        side=6,sources=6,receivers=16,k=float(2*np.pi),material_box=[.02,20],
        new_objects=6,families=['ring','islands','curved'],amplitudes=[10.,16.],
        backgrounds=[.07,.16],layouts=['uniform','jittered'],observations=['clean','noisy_1e-3'],
        anonymous_solver_inputs='G/E0/S/Y/grid_shape/public shared starts only',
        starts=[0,1,4],wall_per_start=4.,max_lu=120,max_rhs=2400,
        methods=[dict(method='native',chart='csi'),dict(method='feedback',chart='csi'),
                 dict(method='cie_dual',chart='csi'),dict(method='material_seeds',chart='csi'),
                 dict(method='n1_fourier',chart='cie'),dict(method='n1_task',chart='cie'),
                 dict(method='n1_finite',chart='cie'),dict(method='n1_finite_fourier',chart='cie'),
                 dict(method='n1_finite_random',chart='cie')],
        n1=dict(affine_mode='native_tsom',native_mf=2,native_core_apcs=15,real_current_DOF=192,
                beta=.2,rho=1.,pool_groups=8,joint_steps=80,operation='replace',
                source_normalized_joint=True,original_physical_score=True,local_gain_gate=False,
                seed_rms_max=2.,seed_current_over_E_cap=3.,seed_defect_over_E_cap=3.,seed_scale_trials=8),
        gates=dict(nrmse=.15,iou=.7,data_clean=.001,data_noisy=.003),
        prior='unchanged R4 smooth/TV 1e-4 each; no shape/template support',
        selection='minimum original normalized full-data score plus shared prior; GT only after all predictions close',
        scope='new native-affine confirmation after development generic finite gain; not A17 attribution or 576 scaling; six varied objects across three families')
    path=HERE/'data/NATIVE_CONFIRMATION_FREEZE.json'
    path.write_text(json.dumps(plan,indent=2))
    with zipfile.ZipFile(HERE/'audit/NATIVE_CONFIRMATION_SOURCE.zip','w',zipfile.ZIP_DEFLATED) as z:
        for folder in ('src','experiments'):
            for f in (HERE/folder).glob('*.py'):z.write(f,f.relative_to(HERE))
    print(json.dumps(dict(status=plan['status'],utc=plan['utc'])))


def generate():
    plan=json.loads((HERE/'data/NATIVE_CONFIRMATION_FREEZE.json').read_text())
    private=HERE.parents[3]/'fbg_r3e_native_confirmation_private_20261101'
    private.mkdir(exist_ok=False)
    public=HERE/'data/native_confirmation_public';public.mkdir(exist_ok=False)
    rng=np.random.default_rng(secrets.randbits(128))
    side=plan['side'];old=HERE/'data/e1_public'
    operators={}
    for layout in plan['layouts']:
        shutil.copy2(old/f'operators_{side}_{layout}.npz',public/f'operators_{side}_{layout}.npz')
        operators[layout]=load_model(public/f'operators_{side}_{layout}.npz')[0]
    shutil.copy2(old/f'starts_{side}.npz',public/f'starts_{side}.npz')
    yy,xx=np.meshgrid((np.arange(side*16)+.5)/(side*16)*2-1,
                      (np.arange(side*16)+.5)/(side*16)*2-1,indexing='ij')
    records=[];cases=[];started=time.perf_counter()
    for family in plan['families']:
        for amplitude in plan['amplitudes']:
            dx,dy=rng.uniform(-.1,.1,2);angle=rng.uniform(-.32,.32)
            x=np.cos(angle)*(xx-dx)+np.sin(angle)*(yy-dy)
            y=-np.sin(angle)*(xx-dx)+np.cos(angle)*(yy-dy)
            if family=='ring':
                radius=np.hypot(x/rng.uniform(.94,1.08),y/rng.uniform(.92,1.1))
                outer=rng.uniform(.61,.75);width=rng.uniform(.21,.3)
                mask=(radius<outer)&(radius>outer-width)
            elif family=='islands':
                rr1,rr2=rng.uniform(.24,.34,2)
                mask=((x+.36)**2+(y+.17)**2<rr1**2)|((x-.34)**2+(y-.19)**2<rr2**2)
            else:
                radius=np.hypot(x,y);phi=np.arctan2(y,x)
                inner=rng.uniform(.3,.42);outer=inner+rng.uniform(.23,.3)
                mask=(radius>inner)&(radius<outer)&(phi>rng.uniform(-1.25,-.55))&(phi<rng.uniform(2.1,2.8))
            alpha=mask.reshape(side,16,side,16).mean((1,3)).ravel()
            bg=rng.uniform(*plan['backgrounds']);m=bg+amplitude*alpha
            oid=rng.bytes(8).hex();truth=private/f'{oid}.npz'
            np.savez_compressed(truth,m=m,alpha=alpha,mask=alpha>.35)
            for layout,model in operators.items():
                clean=model.forward(m).copy()
                noise=rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape)
                noise*=1e-3*np.linalg.norm(clean)/np.linalg.norm(noise)
                for name,Y in [('clean',clean),('noisy',clean+noise)]:
                    case='obs_'+rng.bytes(8).hex()
                    np.savez_compressed(public/(case+'.npz'),Y=Y,operators_path=np.array(f'operators_{side}_{layout}.npz'),grid_shape=np.array([side,side]))
                    cases.append(dict(case_id=case,side=side,layout=layout,noise=name))
                    records.append(dict(case_id=case,object_id=oid,side=side,shape=family,amplitude=amplitude,
                        layout=layout,noise=name,background=float(bg),truth_file=truth.name))
    (private/'mapping.json').write_text(json.dumps(records,indent=2))
    (public/'cases.json').write_text(json.dumps(cases,indent=2))
    receipt=dict(status='NEW_DATA_GENERATED_AFTER_METHOD_FREEZE',utc=datetime.now(timezone.utc).isoformat(),
        freeze_utc=plan['utc'],material_objects=6,observations=24,physical_forward_counts={k:dict(v.cost) for k,v in operators.items()},
        wall=time.perf_counter()-started,private_truth_read_by_inference=False)
    (HERE/'audit/NATIVE_CONFIRMATION_GENERATION.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt))


def run():
    # This branch has no private path, truth, amplitude, support, or RNG state.
    from e1_current_release import infer
    plan=json.loads((HERE/'data/NATIVE_CONFIRMATION_FREEZE.json').read_text())
    public=HERE/'data/native_confirmation_public'
    cases=json.loads((public/'cases.json').read_text())
    tasks=[dict(**c,**method,start=start,wall=plan['wall_per_start'],tag='e1_native_confirmation',
        public_dir=str(public),backbone='observable',joint_steps=80,max_lu=120,max_rhs=2400,affine_mode='native_tsom')
        for c in cases for start in plan['starts'] for method in plan['methods']]
    out=HERE/'results/e1_native_confirmation';out.mkdir(exist_ok=False)
    (out/'PLAN.json').write_text(json.dumps(dict(status='FROZEN_BEFORE_PREDICTIONS',tasks=tasks,truth_access=False),indent=2))
    started=time.perf_counter();rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(infer,t) for t in tasks]
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            with (out/'rows.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
            print(json.dumps({k:row[k] for k in ('run_id','status','relative_residual','wall_seconds')}),flush=True)
    summary=dict(status='PREDICTIONS_CLOSED_TRUTH_NOT_READ',utc=datetime.now(timezone.utc).isoformat(),
        rows=rows,count=len(rows),batch_wall=time.perf_counter()-started,workers=4,blas_threads_per_worker=1,
        classification='PROSPECTIVE_FROZEN_NATIVE_AFFINE_CONFIRMATION')
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','generate','run']);a=ap.parse_args()
    {'freeze':freeze,'generate':generate,'run':run}[a.action]()
