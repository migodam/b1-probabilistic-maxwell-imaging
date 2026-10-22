"""Matched geometry/order/output controls; all physics uses Q book."""
from pathlib import Path
import json,gc
import numpy as np
from p0.geometry import antennas
from p0.qualification import cylinder_reference
from wq.fem import FEM
from wq.meshes import interface_mesh,circle_loop
from wq.high_precision import line_mie
from wq.accounting import Book,encode
root=Path(__file__).resolve().parents[1];out=root/'results/circle-controls';out.mkdir(exist_ok=False)
book=Book(root/'results');tx,rx=antennas();center=np.array([0.,1.5]);ep=4+.03j;r=.3;angles=np.linspace(1.3,1.8,4)
ref=cylinder_reference(rx,2*np.pi,center,r,ep,angles);lr=line_mie(tx,rx,2*np.pi,center,r,ep)
def cyl(p):return 1+(ep-1)*(np.linalg.norm(p-center,axis=1)<r)
rows=[];arrays=dict(plane_reference=ref,line_reference=lr)
for h,degree,n in [(.16,2,256),(.16,3,256),(.13,3,512)]:
    mesh=interface_mesh(h=h,loops=[circle_loop(center,r,n)])
    p=FEM(1.,tx,rx,book,h=h,order=degree,mesh=mesh,receiver_mode='volume')
    s=p.solve(cyl,angles=angles,plane_wave_origin=center);l=p.solve(cyl)
    key=f'h{h}_p{degree}_n{n}'
    arrays[key+'_plane_volume']=s.scattered;arrays[key+'_plane_probe']=s.probe_scattered;arrays[key+'_line_volume']=l.scattered
    row=dict(h=h,degree=degree,n=n,ndofs=p.ndofs,plane_volume_relative=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),plane_probe_relative=float(np.linalg.norm(s.probe_scattered-ref)/np.linalg.norm(ref)),line_volume_relative=float(np.linalg.norm(l.scattered-lr)/np.linalg.norm(lr)),line_volume_white=float(np.linalg.norm(l.scattered-lr)/.01),factor_reused=l.metadata['factor_reused'],residual=max(s.residual,l.residual))
    rows.append(row);(out/'RESULTS.json').write_text(json.dumps(rows,indent=2,default=encode));np.savez_compressed(out/'ARRAYS.npz',**arrays);print(json.dumps(row,default=encode),flush=True)
    del p,s,l;gc.collect()
