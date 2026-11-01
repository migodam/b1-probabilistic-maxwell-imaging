"""Generation is separate from inversion; private truth is not solver input."""
from pathlib import Path
import json, platform, sys, subprocess, time
import numpy as np
import scipy
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from src.physics import make_grid,MaxwellModel


def generate():
    private=HERE.parents[3]/'fbg_r3e_private_20261101'
    private.mkdir(exist_ok=True)
    public=HERE/'data/e1_public'
    public.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(202611010519)
    records=[]
    started=time.perf_counter()
    for side in (4,6):
        geom={layout:make_grid(side,layout) for layout in ('uniform','jittered')}
        for layout,op in geom.items():
            np.savez_compressed(public/f'operators_{side}_{layout}.npz',**op)
        yy,xx=np.meshgrid((np.arange(side*12)+.5)/(side*12)*2-1,
                          (np.arange(side*12)+.5)/(side*12)*2-1,indexing='ij')
        for label in ('ring','islands','curved'):
            if label=='ring':
                rad=np.hypot(xx+.04,yy-.03)
                mask=(rad>.34)&(rad<.73)
            elif label=='islands':
                mask=((xx+.4)**2+(yy+.15)**2<.29**2)|((xx-.31)**2+(yy-.2)**2<.33**2)
            else:
                rad=np.hypot(xx-.01,yy+.04)
                ang=np.arctan2(yy+.04,xx-.01)
                mask=(rad>.37)&(rad<.69)&(ang>-.8)&(ang<4.0)
            alpha=mask.reshape(side,12,side,12).mean((1,3)).ravel()
            for amplitude in (8.,18.):
                m=.1+amplitude*alpha
                oid=rng.bytes(8).hex()
                truthfile=private/f'{oid}.npz'
                np.savez_compressed(truthfile,m=m,mask=alpha>.35,alpha=alpha)
                for layout,op in geom.items():
                    case='obs_'+rng.bytes(8).hex()
                    model=MaxwellModel(op)
                    y=model.forward(m)
                    np.savez_compressed(public/f'{case}.npz',Y=y,
                        operators_path=np.array(f'operators_{side}_{layout}.npz'),grid_shape=np.array([side,side]))
                    records.append(dict(case_id=case,object_id=oid,shape=label,amplitude=amplitude,
                        side=side,layout=layout,truth_file=truthfile.name,forward_cost=dict(model.cost)))
        t=np.linspace(0,1,side)
        ya,xa=np.meshgrid(t,t,indexing='ij')
        starts=np.stack((np.full(side*side,.1),np.ones(side*side),np.full(side*side,4.),
            (1+.55*np.cos(2*np.pi*xa)*np.cos(np.pi*ya)).ravel(),
            (3.+1.2*np.sin(np.pi*xa+.3)*np.cos(2*np.pi*ya+.7)).ravel()))
        np.savez_compressed(public/f'starts_{side}.npz',starts=starts)
    op=make_grid(4,'uniform')
    yy,xx=np.meshgrid(np.linspace(-1,1,4),np.linspace(-1,1,4),indexing='ij')
    m=(.1+10*np.exp(-((xx-.1)**2+(yy+.1)**2)/.5)).ravel()
    np.savez_compressed(HERE/'data/dev_calibration.npz',**op,Y=MaxwellModel(op).forward(m),m_disclosed=m)
    (private/'mapping.json').write_text(json.dumps(records,indent=2))
    (HERE/'data/e1_case_ids.json').write_text(json.dumps([
        {k:r[k] for k in ('case_id','side','layout')} for r in records],indent=2))
    repo=HERE.parents[1]
    manifest=dict(started_utc='2026-11-01T05:19:08Z',
        current_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        historical_b1_commit='6659aa49dff2535996d40962791e9a5a6e974adc',
        historical_a17_commit='2ac81b989a2f13d62df9d5a7abd690a33c2189c7',
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
        physical_backend='frozen R4 build_operators and Model; thin composition wrapper',
        threads=1,material_box=[.02,20],physical_grids=[4,6],sources=6,receivers=16,
        k=float(2*np.pi),source_radius=1.4,receiver_radius=1.1,domain_width=.96,
        pulse_radius_spacing_ratio=.45,assembly='new assembly per grid, not matrix crop',
        generation_wall=time.perf_counter()-started,private_truth=str(private),online_truth_access=False)
    (HERE/'audit/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(dict(observations=len(records),objects=len({r['object_id'] for r in records}),
        wall=time.perf_counter()-started,public=str(public),private=str(private))))


if __name__=='__main__':generate()
