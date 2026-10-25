"""Geometry-only DOF planning after explicit memory rejection; no PDE actions."""
from pathlib import Path
import json
from wslam.forward import public_mesh
root=Path(__file__).resolve().parents[1];rows=[]
for n in (256,384):
    for pose in ((.2,2.),(1.,2.6)):
        m=public_mesh(pose,.12,h=.17,n=n,grid_refinement=3)
        count=int(m.nvertices+2*m.nfacets+m.nelements)
        rows.append(dict(h=.17,n=n,grid_refinement=3,pose=pose,ndofs=count,before_basis=True,physical_solves=0,
                         within_original_60k_guard=count<=60000))
(root/'results/e0/MESH_RESOURCE_PLAN.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows))
