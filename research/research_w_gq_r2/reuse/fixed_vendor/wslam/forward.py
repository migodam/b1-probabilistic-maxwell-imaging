"""Paid full-wave forward and independently auditable material/shape J."""
import gc,time
from dataclasses import dataclass
import numpy as np
from .fem_backend import FEM,PMLConfig
from .physics import acquisition,point_receiver_adjoint,point_coefficient_jacobian,surface_shape_jacobian
from .world import circle_vertices,boundary_quadrature

def public_mesh(pose,radius,*,h=.22,n=128,pml=None,grid_refinement=1):
    import triangle
    from skfem import MeshTri
    pm=PMLConfig() if pml is None else pml
    vertices=[];segments=[];vertex_map={};edge_set=set()
    def vid(point):
        key=tuple(float(x) for x in point)
        if key not in vertex_map:vertex_map[key]=len(vertices);vertices.append(key)
        return vertex_map[key]
    def edge(a,b):
        a,b=vid(a),vid(b);key=tuple(sorted((a,b)))
        if a!=b and key not in edge_set:edge_set.add(key);segments.append((a,b))
    def loop(v):
        for a,b in zip(v,np.roll(v,-1,axis=0)):edge(a,b)
    def rect(box):
        a,b,c,d=box;loop(np.array([[a,c],[b,c],[b,d],[a,d]]))
    rect(pm.outer_box);rect(pm.physical_box);rect((-2,2,.25,3))
    # ALL public macro interfaces, not the unknown occupied subset. The same
    # mesh support is available to every occupancy model and optimizer.
    if grid_refinement not in (1,2,3):raise ValueError('Unsupported field-grid refinement')
    x=np.linspace(-1.5,1.5,24*grid_refinement+1);y=np.linspace(.5,1.25,6*grid_refinement+1)
    for yi in y:
        for a,b in zip(x[:-1],x[1:]):edge((a,yi),(b,yi))
    for xi in x:
        for a,b in zip(y[:-1],y[1:]):edge((xi,a),(xi,b))
    loop(circle_vertices(pose,radius,n))
    mesh=triangle.triangulate(dict(vertices=np.asarray(vertices),segments=np.asarray(segments,dtype=np.int32)),f'pq28a{h*h/2:.12g}Q')
    if 'triangles' not in mesh:raise RuntimeError('PSLG mesh failed')
    return MeshTri(mesh['vertices'].T,mesh['triangles'].T)

@dataclass
class Evaluation:
    data: object
    jacobian: object
    names: object
    backend: object
    state: object
    metadata: object

class Evaluator:
    def __init__(self,book,stage,*,h=.17,n=256,order=3,pml=None,tx=None,rx=None,grid_refinement=2):
        self.book,self.stage=book,stage;self.h,self.n,self.order=h,n,order
        self.pml=PMLConfig() if pml is None else pml
        self.grid_refinement=grid_refinement
        base_tx,base_rx=acquisition();self.tx=base_tx if tx is None else np.asarray(tx,float);self.rx=base_rx if rx is None else np.asarray(rx,float)
        self.backend=None;self.key=None;self.last=None;self.last_world_key=None
    def clear(self):self.backend=None;self.key=None;self.last=None;self.last_world_key=None;gc.collect()
    def get_backend(self,world,f):
        key=(*world.pose,world.radius,f)
        if self.backend is None or key!=self.key:
            self.clear();mesh=public_mesh(world.pose,world.radius,h=self.h,n=self.n,pml=self.pml,grid_refinement=self.grid_refinement)
            self.backend=FEM(f,self.tx,self.rx,self.book,h=self.h,order=self.order,mesh=mesh,pml=self.pml,phase=self.stage,receiver_mode='probe')
            self.key=key
        return self.backend
    def evaluate(self,world,f=1.,*,jacobian=False,kind='shared',include_target=True):
        started=time.perf_counter();p=self.get_backend(world,f)
        worldkey=(world.support.tobytes(),np.asarray(world.epsinf).tobytes(),np.asarray(world.ell).tobytes(),*world.pose,world.radius,world.contrast,world.target_loss,float(f),include_target)
        if self.last_world_key==worldkey and self.last is not None:
            if not jacobian or self.last.jacobian is not None and self.last.metadata.get('jacobian_kind')==kind:return self.last
            s=self.last.state;reused=True
        else:s=p.solve(world.coefficient(f,self.n,include_target));reused=False
        j=None;names=None
        if jacobian:
            q,residual=point_receiver_adjoint(p,s)
            d=point_coefficient_jacobian(p,s,lambda points,freq:world.material_directions(points,freq,self.n,kind,include_target),qnodes=q)
            if include_target:
                pts,normal,ds,radial=boundary_quadrature(world.pose,world.radius,self.n)
                shape=surface_shape_jacobian(p,s,q,pts,normal,ds,radial,contrast=world.contrast+1j*world.target_loss/f)
            else:shape=np.zeros((s.data.size,3),complex)
            j=np.column_stack((d[:,:2],shape[:,:2],d[:,2:],shape[:,2:]))
            names=(['epsinf','ell'] if kind=='shared' else ['eps_left','eps_right'])+['pose_x','pose_y','target_contrast','target_loss','target_radius']
            del q
        meta=dict(world=world.metadata(),freq=f,ndofs=p.ndofs,order=self.order,circle_polygon_n=self.n,residual=s.residual,wall_seconds=time.perf_counter()-started,
                  target_included=include_target,derivative='physical coefficient and shape, independently FD audited' if jacobian else None,jacobian_kind=kind,
                  forward_field_reused_for_J=reused)
        self.book.event(self.stage,'evaluation_complete',**meta)
        self.last=Evaluation(s.data.copy(),j,names,p,s,meta);self.last_world_key=worldkey
        return self.last

def volume_extraction(ev):
    p,s=ev.backend,ev.state
    value=p.k**2*np.einsum('rn,n,ns->sr',p._receiver_green,s.material_weights*(s.epsilon-1),s.field,optimize=True)
    return s.direct+value
