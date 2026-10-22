"""Independent field refinement of online-selected endpoints; no truth choice."""
from pathlib import Path
import json,gc
import numpy as np
from wq.search import Forward,reference,real
from wq.fem import FEM
from wq.accounting import SearchBook,encode
from wq.inverse_support import PixelGrid,pixel_mesh,pixel_material
root=Path(__file__).resolve().parents[1];out=root/'results/p0b/finer-endpoints';out.mkdir(exist_ok=False)
cfg=json.loads((root/'docs/P0_FROZEN_CONFIG.json').read_text());rows=json.loads((root/'results/p0b/search-corrected/RESULTS.json').read_text());book=SearchBook(root/'results/p0b');g=PixelGrid();base=pixel_mesh(g);mesh=base.refined(1)
saved=np.load(root/'results/p0b/search-corrected/ARRAYS.npz');task_saved=np.load(root/'results/p0b/task-endpoints/ARRAYS.npz');probe=cfg['probes'][8]
reports=[];arrays={}
contexts=sorted(set((r['seed'],r['band']) for r in rows if r['case']=='calibration_candidate'))
for seed,band in contexts:
    candidates=[r for r in rows if r['case']=='calibration_candidate' and r['seed']==seed and r['band']==band]
    chosen=min(candidates,key=lambda r:(r['d_cal'],r['method']));method=chosen['method'];key=f'{seed}_{band}_{method}';theta=saved[key];data={}
    for f in cfg['frequencies'][band]:
        p=FEM(f,np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)],np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)],book,phase='validation',h=.09,order=2,mesh=mesh,receiver_mode='volume')
        for label,pr in [('cal',None),('task',probe)]:
            basecoef=pixel_material(theta,f,g)
            if pr is None:coef=basecoef
            else:
                c,rad,val=pr
                def coef(pts,f=f,c=c,rad=rad,val=val,basecoef=basecoef):
                    a=basecoef(pts);m=np.linalg.norm(pts-np.asarray(c),axis=1)<=rad;a[m]=1+val+1j*.02/f;return a
            s=p.solve(coef);data.setdefault(label,[]).append(s.data.copy())
            if f==1.:
                arrays[key+'_'+label+'_points']=s.material_points;arrays[key+'_'+label+'_epsilon']=s.epsilon;arrays[key+'_'+label+'_field']=s.field
        del p,s;gc.collect()
    for label in ('cal','task'):
        current=np.stack(data[label]);old=saved[key+'_data'] if label=='cal' else task_saved[key+'_task_data']
        r=dict(seed=seed,band=band,method=method,state=label,finer_method='P2 uniform subdivision of every original P3 triangle; same pixel coefficient/point sources',mesh_white=float(np.linalg.norm(current-old)/.01),data_distance_to_reference=float(np.linalg.norm(real(current.ravel()-reference(seed,band,label)))),qualification='EMPIRICAL_ENDPOINT_CHECK_NOT_CONTINUUM_BOUND')
        reports.append(r);arrays[key+'_'+label+'_data']=current
        (out/'RESULTS.json').write_text(json.dumps(reports,indent=2,default=encode));np.savez_compressed(out/'ARRAYS.npz',**arrays);print(json.dumps(r),flush=True)
