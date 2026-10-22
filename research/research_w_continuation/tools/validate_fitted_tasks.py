"""Same selected endpoints/probe8, circle-fitted P3 polygon refinement only."""
from pathlib import Path
import gc, json
import numpy as np
from wq.accounting import SearchBook, encode
from wq.fem import FEM
from wq.inverse_support import PixelGrid, pixel_mesh, pixel_material
from wq.target_mesh import fit_probe_mesh
from wq.search import reference, real

root=Path(__file__).resolve().parents[1]
out=root/'results/p0b/fitted-task-endpoints';out.mkdir(exist_ok=False)
cfg=json.loads((root/'docs/P0_FROZEN_CONFIG.json').read_text())
rows=json.loads((root/'results/p0b/search-corrected/RESULTS.json').read_text())
saved=np.load(root/'results/p0b/search-corrected/ARRAYS.npz')
old=np.load(root/'results/p0b/task-endpoints/ARRAYS.npz')
book=SearchBook(root/'results/p0b');g=PixelGrid();base=pixel_mesh(g)
center,radius,contrast=cfg['probes'][8]
tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)]
rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
contexts=sorted(set((r['seed'],r['band']) for r in rows if r['case']=='calibration_candidate'))
chosen=[]
for seed,band in contexts:
    candidates=[r for r in rows if r['case']=='calibration_candidate' and r['seed']==seed and r['band']==band]
    chosen.append(min(candidates,key=lambda r:(r['d_cal'],r['method'])))
reports=[];data={};diagnostics=[]
for n in (256,512):
    mesh=fit_probe_mesh(base,center,radius,n=n,min_angle=None,max_area=None)
    diagnostics.append(mesh.probe_fit_diagnostics)
    (out/'MESH_DIAGNOSTICS.json').write_text(json.dumps(diagnostics,indent=2,default=encode))
    for ch in chosen:
        seed,band,method=ch['seed'],ch['band'],ch['method'];key=f'{seed}_{band}_{method}'
        theta=saved[key];current=[];residual=[];ndofs=[]
        for f in cfg['frequencies'][band]:
            p=FEM(f,tx,rx,book,phase='validation',h=.18,order=3,mesh=mesh,receiver_mode='volume')
            basecoef=pixel_material(theta,f,g)
            def coef(pts,basecoef=basecoef,f=f):
                a=basecoef(pts)
                inside=np.linalg.norm(pts-np.asarray(center),axis=1)<=radius
                a[inside]=1+contrast+1j*.02/f
                return a
            s=p.solve(coef);current.append(s.data.copy());residual.append(s.residual);ndofs.append(p.ndofs)
            if f==1.:
                np.savez_compressed(out/f'{key}_n{n}_FIELD.npz',points=s.material_points,epsilon=s.epsilon,field=s.field,data=s.data)
            del p,s;gc.collect()
        current=np.stack(current);data[f'{key}_n{n}']=current
        row=dict(seed=seed,band=band,method=method,state='fixed_probe8_task',polygon_n=n,field_order=3,
                 unfitted_difference_white=float(np.linalg.norm(current-old[key+'_task_data'])/.01),
                 data_distance_to_reference=float(np.linalg.norm(real(current.ravel()-reference(seed,band,'task')))),
                 maximum_rhs_relative_residual=max(residual),ndofs=max(ndofs),
                 geometry_scope='same pixel material and true circle; fitted inscribed polygon refined 256 to 512',
                 skinny_mesh=True,qualification='EMPIRICAL_DISCRETE_CHECK_NOT_CONTINUUM_BOUND')
        if n==512:
            row['fitted_refinement_white']=float(np.linalg.norm(current-data[f'{key}_n256'])/.01)
            row['passes_original_0_1_consistency']=bool(row['fitted_refinement_white']<=.1)
        reports.append(row)
        (out/'RESULTS.json').write_text(json.dumps(reports,indent=2,default=encode))
        np.savez_compressed(out/'DATA.npz',**data)
        print(json.dumps(row,default=encode),flush=True)
    del mesh;gc.collect()
